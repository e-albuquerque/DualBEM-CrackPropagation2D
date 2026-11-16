#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
preproc.py — Pré-processador Gmsh + meshio para o Solver Dual-BEM (Fracture 2D)
===============================================================================

Fluxo de trabalho (baseado nos projetos BEM2D_elastostatic e BEM3D_elastostatic):
  1. Lê as configurações do problema em `input_data.py`;
  2. Executa o Gmsh sobre `meshes/<file_name>.geo` para gerar a malha quadrática
     contínua `meshes/<file_name>.msh` (`line3`);
  3. Lê o arquivo `.msh` utilizando `meshio`;
  4. Converte TODOS os elementos contínuos `line3` do Gmsh em **elementos
     quadráticos descontínuos de 3 nós** (nós de colocação em eta = -lamb, 0, +lamb),
     e duplica automaticamente as curvas de trinca nas duas faces do método Dual-BEM:
       - Face 1 (bie = 2): Equação Integral de Contorno de Deslocamento;
       - Face 2 (bie = 3): Equação Integral de Contorno Hipersingular de Tração;
  5. Gera o arquivo de entrada `.dat` limpo e legível para o solver Fortran;
  6. Plota e salva a figura da malha e condições de contorno (`<id>_mesh_bc.png`).
"""

import argparse
import math
import pathlib
import subprocess
import sys

import meshio
import numpy as np

from input_data import input_data
from plot_mesh_bc import parse_input_data, comp_xyplot, detect_crack_elements, show_geometry


ROOT_DIR = pathlib.Path(__file__).resolve().parent
MESH_DIR = ROOT_DIR / "meshes"


def geom_shape_functions(eta: float) -> tuple[float, float, float]:
    """Funções de forma quadráticas contínuas M1, M2, M3 em eta in [-1, 1] (nós em -1, 0, +1)."""
    m1 = 0.5 * eta * (eta - 1.0)
    m2 = 1.0 - eta * eta
    m3 = 0.5 * eta * (eta + 1.0)
    return m1, m2, m3


def eval_quad_point(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray, eta: float) -> np.ndarray:
    """Avalia a posição cartesiana (x, y) em eta in [-1, 1] sobre o elemento quadrático (p1, p2, p3)."""
    m1, m2, m3 = geom_shape_functions(eta)
    return m1 * p1 + m2 * p2 + m3 * p3


def generate_msh_from_geo(geo_path: pathlib.Path, msh_path: pathlib.Path, gmsh_params: dict | None = None) -> None:
    """Gera o arquivo `.msh` de ordem 2 a partir do `.geo` usando a API Python do Gmsh (ou CLI)."""
    gmsh_params = gmsh_params or {}
    try:
        import gmsh

        argv = [""]
        for k, v in gmsh_params.items():
            argv.extend(["-setnumber", str(k), str(v)])
        gmsh.initialize(argv)
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(str(geo_path))
        gmsh.model.mesh.generate(2)
        gmsh.model.mesh.setOrder(2)
        gmsh.write(str(msh_path))
        gmsh.finalize()
    except ImportError:
        cmd = ["gmsh", str(geo_path), "-2", "-order", "2", "-o", str(msh_path)]
        for k, v in gmsh_params.items():
            cmd.extend(["-setnumber", str(k), str(v)])
        subprocess.run(cmd, check=True)


def order_chain_elements(elem_pts: list[tuple[np.ndarray, np.ndarray, np.ndarray]]) -> list[tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Garante que os elementos `line3` de uma trinca ou curva formem uma cadeia
    contínua ordenada do início ao fim (p3[k] == p1[k+1]).
    """
    if len(elem_pts) <= 1:
        return elem_pts

    tol = 1e-7
    n = len(elem_pts)
    # Verifica se já estão encadeados em sequência
    already_ordered = True
    for k in range(n - 1):
        if np.linalg.norm(elem_pts[k][2] - elem_pts[k + 1][0]) > tol:
            already_ordered = False
            break
    if already_ordered:
        return elem_pts

    # Encontra o primeiro elemento (cujo p1 não coincide com o p3 de nenhum outro)
    start_idx = 0
    for i in range(n):
        p1_i = elem_pts[i][0]
        is_start = True
        for j in range(n):
            if i != j and np.linalg.norm(p1_i - elem_pts[j][2]) <= tol:
                is_start = False
                break
        if is_start:
            start_idx = i
            break

    ordered = [elem_pts[start_idx]]
    used = {start_idx}
    while len(ordered) < n:
        last_p3 = ordered[-1][2]
        found_next = False
        for j in range(n):
            if j in used:
                continue
            if np.linalg.norm(elem_pts[j][0] - last_p3) <= tol:
                ordered.append(elem_pts[j])
                used.add(j)
                found_next = True
                break
        if not found_next:
            # Adiciona os restantes na ordem original caso sejam curvas disjuntas
            for j in range(n):
                if j not in used:
                    ordered.append(elem_pts[j])
                    used.add(j)
            break
    return ordered


def decode_bc_spec(bc_type: str, bc_val, s_param: float) -> tuple[int, float, float]:
    """
    Converte a especificação de condição de contorno do `input_data.py` para o
    formato numérico do Fortran: `(bck, val_a, val_b)` onde:
      - bck = 0 (deslocamento prescrito): val_a = u, val_b = 0.0
      - bck = 1 (tração prescrita):       val_a = t, val_b = 0.0
      - bck = 2 (mola elástica Robin):    val_a = k, val_b = t0
    """
    key = bc_type.strip().lower()
    if key in ("ux", "uy", "un", "ut", "disp"):
        bck = 0
    elif key in ("tx", "ty", "tn", "tt", "trac"):
        bck = 1
    elif key in ("spring_x", "spring_y", "spring_n", "spring_t", "spring", "robin"):
        bck = 2
    else:
        raise ValueError(f"Tipo de condição de contorno desconhecido: '{bc_type}'")

    if bck == 2:
        if isinstance(bc_val, (tuple, list)):
            return 2, float(bc_val[0]), float(bc_val[1]) if len(bc_val) > 1 else 0.0
        return 2, float(bc_val), 0.0

    # Para deslocamento (bck=0) ou tração (bck=1): pode ser escalar, tupla polinomial (b1, b2, b3, b4) ou função
    if callable(bc_val):
        val = float(bc_val(s_param))
    elif isinstance(bc_val, (tuple, list)):
        coeffs = list(bc_val) + [0.0, 0.0, 0.0, 0.0]
        b1, b2, b3, b4 = coeffs[:4]
        val = float(b1 + b2 * s_param + b3 * (s_param**2) + b4 * (s_param**3))
    else:
        val = float(bc_val)

    return bck, val, 0.0


def compute_curve_s_params(curve_elems: list[tuple[np.ndarray, np.ndarray, np.ndarray]], lamb: float) -> list[tuple[float, float, float]]:
    """
    Calcula a coordenada paramétrica normalizada s in [0, 1] nos 3 nós funcionais
    descontínuos (eta = -lamb, 0, +lamb) de cada elemento ao longo de uma curva,
    reproduzindo exatamente a parametrização por comprimento de corda de
    `get_parametric_coords` na dissertação.
    """
    n_el = len(curve_elems)
    if n_el == 0:
        return []

    # Constrói a lista de 2*n_el + 1 pontos geométricos da curva
    pts = [curve_elems[0][0]]
    for p1, p2, p3 in curve_elems:
        pts.append(p2)
        pts.append(p3)

    n_pts = len(pts)  # 2 * n_el + 1
    s_pts = np.zeros(n_pts, dtype=float)

    pa = eval_quad_point(pts[0], pts[1], pts[2], -lamb)
    r0 = float(np.linalg.norm(pa - pts[0]))
    s_pts[0] = r0

    for i in range(1, n_pts - 1):
        r0 += float(np.linalg.norm(pts[i] - pa))
        s_pts[i] = r0
        pa = pts[i]

    pb = eval_quad_point(pts[-3], pts[-2], pts[-1], +lamb)
    r0 += float(np.linalg.norm(pb - pts[-2]))
    s_pts[-1] = r0
    r0 += float(np.linalg.norm(pts[-1] - pb))

    if r0 > 0.0:
        s_pts /= r0

    # Para cada elemento i (0..n_el-1), associa s nos 3 nós
    s_per_elem = []
    for i in range(n_el):
        s1 = float(s_pts[i])
        s2 = float(s_pts[2 * i + 1])
        s3 = float(s_pts[2 * i + 2])
        s_per_elem.append((s1, s2, s3))
    return s_per_elem


def generate_internal_points(rects: list[dict]) -> list[tuple[float, float]]:
    """Gera a grade de pontos internos a partir de áreas retangulares (x1, y1) -> (x2, y2)."""
    int_pts = []
    for r in rects:
        x1, y1 = float(r["x1"]), float(r["y1"])
        x2, y2 = float(r["x2"]), float(r["y2"])
        nx, ny = int(r["nx"]), int(r["ny"])
        xs = np.linspace(x1, x2, nx) if nx > 1 else np.array([0.5 * (x1 + x2)])
        ys = np.linspace(y1, y2, ny) if ny > 1 else np.array([0.5 * (y1 + y2)])
        for yv in ys:
            for xv in xs:
                int_pts.append((float(xv), float(yv)))
    return int_pts


def build_discontinuous_model(data: dict, mesh: meshio.Mesh) -> dict:
    """
    Transforma a malha contínua `line3` do Gmsh em elementos quadráticos
    descontínuos de 3 nós (para o contorno externo e para as duas faces de
    cada trinca no método Dual-BEM).
    """
    points_2d = np.asarray(mesh.points[:, :2], dtype=float)
    if "line3" not in mesh.cells_dict:
        raise RuntimeError("O arquivo .msh não contém elementos quadráticos 1D ('line3'). Verifique Mesh.ElementOrder = 2.")

    raw_line3 = mesh.cells_dict["line3"]
    phys_tags = mesh.cell_data_dict["gmsh:physical"]["line3"]

    # No Gmsh, line3 = [n_ini, n_fim, n_meio] -> reordenamos para [n_ini (eta=-1), n_meio (eta=0), n_fim (eta=+1)]
    line3_ordered = raw_line3[:, [0, 2, 1]]

    # Mapeia nome do Physical Curve -> tag inteiro
    curve_name_to_tag = {}
    for name, (tag, dim) in mesh.field_data.items():
        if dim == 1:
            curve_name_to_tag[name] = int(tag)

    lamb = float(data.get("lamb", 2.0 / 3.0))
    elements = []   # Lista de dicts: bp_1, bp_2, lamb_1, lamb_2, incid(3)
    nodes = []      # Lista de dicts: point, bie, bc_dir, bck1, v1a, v1b, bck2, v2a, v2b
    cracks_out = [] # Lista de dicts com ext_el_num(4), length, tips_num, pairs

    def add_discontinuous_element(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray,
                                  bie: int, bc_dir: int,
                                  bc1_tuple: tuple[int, float, float],
                                  bc2_tuple: tuple[int, float, float],
                                  s_tuple: tuple[float, float, float] | None = None,
                                  bc_cfg: dict | None = None):
        """Adiciona 1 elemento quadrático descontínuo e seus 3 nós exclusivos de colocação."""
        n1_pt = eval_quad_point(p1, p2, p3, -lamb)
        n2_pt = p2.copy()
        n3_pt = eval_quad_point(p1, p2, p3, +lamb)

        inc = []
        for loc_idx, pt in enumerate((n1_pt, n2_pt, n3_pt)):
            if bc_cfg is not None and s_tuple is not None:
                s_val = s_tuple[loc_idx]
                b1 = decode_bc_spec(bc_cfg["type_x"], bc_cfg["value_x"], s_val)
                b2 = decode_bc_spec(bc_cfg["type_y"], bc_cfg["value_y"], s_val)
            else:
                b1 = bc1_tuple
                b2 = bc2_tuple

            node_id = len(nodes) + 1
            nodes.append({
                "id": node_id,
                "x": float(pt[0]),
                "y": float(pt[1]),
                "bie": int(bie),
                "bc_dir": int(bc_dir),
                "bck1": b1[0], "val1_a": b1[1], "val1_b": b1[2],
                "bck2": b2[0], "val2_a": b2[1], "val2_b": b2[2],
            })
            inc.append(node_id)

        elem_id = len(elements) + 1
        elements.append({
            "id": elem_id,
            "x1": float(p1[0]), "y1": float(p1[1]),
            "x2": float(p3[0]), "y2": float(p3[1]),
            "xm": float(p2[0]), "ym": float(p2[1]),
            "lamb_1": lamb,
            "lamb_2": lamb,
            "incid": inc,
        })

    # 1. Processa todos os elementos do contorno externo (bie = 1) na ordem do .msh
    crack_phys_names = {c["physical_group"] for c in data.get("cracks", [])}
    crack_phys_tags = {curve_name_to_tag[name] for name in crack_phys_names if name in curve_name_to_tag}

    # Pré-calcula os parâmetros s in [0, 1] para curvas com condição de contorno polinomial
    poly_s_map = {}
    for cname, bc_cfg in data["bound_cond"].items():
        if cname not in curve_name_to_tag:
            raise KeyError(f"Physical Curve '{cname}' não encontrada no arquivo .msh ({list(curve_name_to_tag.keys())}).")
        ctag = curve_name_to_tag[cname]
        is_poly = isinstance(bc_cfg["value_x"], (tuple, list)) or isinstance(bc_cfg["value_y"], (tuple, list)) or callable(bc_cfg["value_x"]) or callable(bc_cfg["value_y"])
        if is_poly and bc_cfg["type_x"] not in ("spring_x", "spring_n", "spring") and bc_cfg["type_y"] not in ("spring_y", "spring_t", "spring"):
            elem_indices = np.where(phys_tags == ctag)[0]
            curve_elems = [
                (points_2d[line3_ordered[idx, 0]], points_2d[line3_ordered[idx, 1]], points_2d[line3_ordered[idx, 2]])
                for idx in elem_indices
            ]
            s_list = compute_curve_s_params(curve_elems, lamb)
            for local_i, global_idx in enumerate(elem_indices):
                poly_s_map[int(global_idx)] = s_list[local_i]

    tag_to_bc = {curve_name_to_tag[cname]: bc_cfg for cname, bc_cfg in data["bound_cond"].items()}

    for idx in range(len(line3_ordered)):
        ctag = int(phys_tags[idx])
        if ctag in crack_phys_tags:
            continue
        if ctag not in tag_to_bc:
            continue

        bc_cfg = tag_to_bc[ctag]
        p1 = points_2d[line3_ordered[idx, 0]]
        p2 = points_2d[line3_ordered[idx, 1]]
        p3 = points_2d[line3_ordered[idx, 2]]
        bc_dir = int(bc_cfg.get("bc_dir", 0))

        if idx in poly_s_map:
            add_discontinuous_element(
                p1, p2, p3,
                bie=1, bc_dir=bc_dir,
                bc1_tuple=(1, 0.0, 0.0), bc2_tuple=(1, 0.0, 0.0),
                s_tuple=poly_s_map[idx], bc_cfg=bc_cfg,
            )
        else:
            b1 = decode_bc_spec(bc_cfg["type_x"], bc_cfg["value_x"], 0.0)
            b2 = decode_bc_spec(bc_cfg["type_y"], bc_cfg["value_y"], 0.0)
            add_discontinuous_element(
                p1, p2, p3,
                bie=1, bc_dir=bc_dir,
                bc1_tuple=b1, bc2_tuple=b2,
            )

    # 2. Processa cada trinca gerando automaticamente as duas faces do Dual-BEM (bie = 2 e bie = 3)
    free_bc = (1, 0.0, 0.0)
    for c_spec in data.get("cracks", []):
        cname = c_spec["physical_group"]
        tips_num = int(c_spec.get("tips_num", 1))
        if cname not in curve_name_to_tag:
            raise KeyError(f"Physical Curve da trinca '{cname}' não encontrada no arquivo .msh.")
        ctag = curve_name_to_tag[cname]
        elem_indices = np.where(phys_tags == ctag)[0]
        raw_crack_elems = [
            (points_2d[line3_ordered[idx, 0]], points_2d[line3_ordered[idx, 1]], points_2d[line3_ordered[idx, 2]])
            for idx in elem_indices
        ]
        crack_elems = order_chain_elements(raw_crack_elems)
        nbe_crack = len(crack_elems)

        ext_el_1 = len(elements) + 1
        ext_el_2 = len(elements) + nbe_crack
        ext_el_3 = len(elements) + nbe_crack + 1
        ext_el_4 = len(elements) + 2 * nbe_crack

        first_node_f1 = len(nodes) + 1
        cr_length = 0.0

        # Face 1 (bie = 2, EIC de Deslocamento): sentido original A -> B
        for p1, p2, p3 in crack_elems:
            cr_length += float(np.linalg.norm(p3 - p1))
            add_discontinuous_element(p1, p2, p3, bie=2, bc_dir=0, bc1_tuple=free_bc, bc2_tuple=free_bc)

        last_node_f1 = len(nodes)
        n_crack_nodes = last_node_f1 - first_node_f1 + 1

        # Face 2 (bie = 3, EIC Hipersingular de Tração): sentido reverso B -> A
        for p1, p2, p3 in reversed(crack_elems):
            add_discontinuous_element(p3, p2, p1, bie=3, bc_dir=0, bc1_tuple=free_bc, bc2_tuple=free_bc)

        last_node_f2 = len(nodes)

        # Pares de nós coincidentes (Face 1 pareada com Face 2 em ordem reversa)
        pairs = []
        for k in range(n_crack_nodes):
            nf1 = first_node_f1 + k
            nf2 = last_node_f2 - k
            pairs.append((nf1, nf2))

        cracks_out.append({
            "id": len(cracks_out) + 1,
            "tips_num": tips_num,
            "n_pairs": n_crack_nodes,
            "ext_el": (ext_el_1, ext_el_2, ext_el_3, ext_el_4),
            "length": cr_length,
            "pairs": pairs,
        })

    # 3. Pontos internos
    int_pts = generate_internal_points(data.get("internal_points_rects", []))

    return {
        "elements": elements,
        "nodes": nodes,
        "cracks": cracks_out,
        "internal_points": int_pts,
    }


def write_fortran_dat_file(data: dict, model: dict, dat_path: pathlib.Path) -> None:
    """Escreve o arquivo `.dat` formatado e documentado para o solver Fortran."""
    elems = model["elements"]
    nodes = model["nodes"]
    cracks = model["cracks"]
    int_pts = model["internal_points"]
    sif = data["sif"]
    fat = data["fatigue"]
    vamp = data["variable_amplitude"]

    lines = [
        "# =============================================================================",
        "# Arquivo de Entrada Dual-BEM (Fracture 2D - Elementos Quadraticos Descontinuos)",
        "# Gerado automaticamente por preproc.py a partir do Gmsh (.geo -> .msh) + meshio",
        "# =============================================================================",
        "# 1. Titulo da Analise e Identificador (Analysis ID)",
        f"{data['title']}",
        f"{data['analysis_id']}",
        "#",
        "# 2. Propriedades do Material e Parametros Gerais:",
        "#    E (Modulo de Young) | nu (Poisson) | plane_state (1=EPT, 2=EPD) | npow (Gao)",
        f"{data['E']:.8e}  {data['nu']:.8f}  {data['plane_state']}  {data['npow']}",
        "#",
        "# 3. Dimensoes da Malha Descontinua:",
        "#    be_num (Elementos) | bn_num (Nos de Colocacao) | cr_num (Trincas) | in_num (Pts Internos)",
        f"{len(elems)}  {len(nodes)}  {len(cracks)}  {len(int_pts)}",
        "#",
        "# =============================================================================",
        "# 4. Elementos de Contorno Quadraticos Descontinuos (3 nos por elemento)",
        "#    elem_id  x1(eta=-1)  y1(eta=-1)  x2(eta=+1)  y2(eta=+1)  lamb_1  lamb_2  inc1  inc2  inc3",
        "# =============================================================================",
    ]

    for e in elems:
        inc1, inc2, inc3 = e["incid"]
        lines.append(
            f"{e['id']:6d}  {e['x1']:20.12e}  {e['y1']:20.12e}  {e['x2']:20.12e}  {e['y2']:20.12e}  "
            f"{e['lamb_1']:14.10f}  {e['lamb_2']:14.10f}  {inc1:6d}  {inc2:6d}  {inc3:6d}"
        )

    lines.extend([
        "#",
        "# =============================================================================",
        "# 5. Nos Funcionais de Colocacao e Condicoes de Contorno",
        "#    bie:    1 = Contorno (Desloc.), 2 = Trinca Face 1 (Desloc.), 3 = Trinca Face 2 (Tracao)",
        "#    bc_dir: 0 = Cartesiano (x, y),  1 = Local (normal, tangencial)",
        "#    bck:    0 = Deslocamento prescrito, 1 = Tracao prescrita, 2 = Mola elastica (k, t0)",
        "#    node_id  x  y  bie  bc_dir  bck1  val1_a  val1_b  bck2  val2_a  val2_b",
        "# =============================================================================",
    ])

    for nd in nodes:
        lines.append(
            f"{nd['id']:6d}  {nd['x']:20.12e}  {nd['y']:20.12e}  {nd['bie']:2d}  {nd['bc_dir']:2d}  "
            f"{nd['bck1']:2d}  {nd['val1_a']:20.12e}  {nd['val1_b']:14.6e}  "
            f"{nd['bck2']:2d}  {nd['val2_a']:20.12e}  {nd['val2_b']:14.6e}"
        )

    lines.extend([
        "#",
        "# =============================================================================",
        "# 6. Topologia das Trincas (Dual-BEM: Pares de Nos Coincidentes Face 1 / Face 2)",
        "#    Cabecalho: crack_id  tips_num  n_pairs  ext_el_1  ext_el_2  ext_el_3  ext_el_4  length",
        "#    Linhas:    pair_idx  node_face1(bie=2)  node_face2(bie=3)",
        "# =============================================================================",
    ])

    for cr in cracks:
        e1, e2, e3, e4 = cr["ext_el"]
        lines.append(
            f"{cr['id']:4d}  {cr['tips_num']:2d}  {cr['n_pairs']:6d}  "
            f"{e1:6d}  {e2:6d}  {e3:6d}  {e4:6d}  {cr['length']:20.12e}"
        )
        for k, (nf1, nf2) in enumerate(cr["pairs"], start=1):
            lines.append(f"  {k:6d}  {nf1:6d}  {nf2:6d}")

    if int_pts:
        lines.extend([
            "#",
            "# =============================================================================",
            "# 7. Coordenadas dos Pontos Internos (pt_id  x  y)",
            "# =============================================================================",
        ])
        for idx, (px, py) in enumerate(int_pts, start=1):
            lines.append(f"{idx:6d}  {px:20.12e}  {py:20.12e}")

    lines.extend([
        "#",
        "# =============================================================================",
        "# 8. Parametros para Calculo dos Fatores de Intensidade de Tensao (KI e KII)",
        "#    Linha 1: K_key (1=Correl. Desloc., 2=Integral J)  method_key",
        "#    Linha 2: nbe_j  nbn_j  lb1_j  lb2_j  np_gauss  initial_node_rel",
        "# =============================================================================",
        f"{sif['k_key']}  {sif['method_key']}",
    ])

    if sif["k_key"] == 1:
        lines.append(f"{sif.get('nd', 1)}")
    elif sif["method_key"] in (1, 2):
        lines.append(f"{sif['np']}  {sif['initial_node']}")
    else:
        lines.append(
            f"{sif['nbe']}  {sif['nbn']}  {sif['lb1']:.10f}  {sif['lb2']:.10f}  {sif['np']}  {sif['initial_node']}"
        )

    lines.extend([
        "#",
        "# =============================================================================",
        "# 9. Parametros de Propagacao de Trinca e Fadiga (Lei de Paris-Walker)",
        "#    Linha 1: k_eff_key (1=Tanaka, 2=MTC/Portela)",
        "#    Linha 2: Da_0  r_da  Da_max  steps",
        "#    Linha 3: R  C0  m  gamma  K_th  K_IC",
        "# =============================================================================",
        f"{fat['k_eff_key']}",
        f"{fat['da0']:.8e}  {fat['r_da']:.8f}  {fat['da_max']:.8e}  {fat['steps']}",
        f"{fat['R']:.10f}  {fat['C0']:.8e}  {fat['m']:.8f}  {fat['gamma']:.8f}  {fat['K_th']:.8f}  {fat['K_IC']:.8f}",
        "#",
        "# =============================================================================",
        "# 10. Carregamento de Amplitude Variavel (0 = Constante, 1 = Variavel / Rainflow)",
        "# =============================================================================",
    ])

    if vamp.get("enabled", False):
        alphas = vamp["alphas"]
        lines.append("1")
        lines.append(f"{len(alphas)}")
        for a_val in alphas:
            lines.append(f"{float(a_val):.8f}")
    else:
        lines.append("0")

    dat_path.parent.mkdir(parents=True, exist_ok=True)
    dat_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_input_data_echo_for_plot(data: dict, model: dict, echo_path: pathlib.Path) -> None:
    """
    Gera também o arquivo `<id>_input_data.txt` compatível com `plot_mesh_bc.py`
    para visualização imediata da malha descontínua e das condições de contorno.
    """
    elems = model["elements"]
    nodes = model["nodes"]
    int_pts = model["internal_points"]
    psk_str = "Plane stress problem." if data["plane_state"] == 1 else "Plane strain problem."

    out = [
        "*" * 82,
        "Fracture 2D",
        "",
        "A Dual-BEM based program - Discontinuous Elements Edition",
        "",
        f"Analysis name: {data['title']}",
        "",
        f" {psk_str}",
        "",
        f" Elasticity modulus =   {data['E']:.12E}",
        f" Poisson ratio = {data['nu']:.8f}",
        "",
        f" Number of boundary elements = {len(elems)}",
        f" Number of boundary points = {len(nodes)}",
    ]
    if int_pts:
        out.append(f" Number of internal points = {len(int_pts)}")

    out.extend([
        "",
        "Boundary Elements:",
        "",
        "Elem.     X1            Y1            X2            Y2",
    ])
    for e in elems:
        out.append(f"{e['id']:5d}  {e['x1']:12.5E}  {e['y1']:12.5E}  {e['x2']:12.5E}  {e['y2']:12.5E}")

    out.extend([
        "",
        "Elem.     Lambda 1      Lambda 2      Incidences",
    ])
    for e in elems:
        i1, i2, i3 = e["incid"]
        out.append(f"{e['id']:5d}  {e['lamb_1']:12.5E}  {e['lamb_2']:12.5E}  {i1:7d}  {i2:7d}  {i3:7d}")

    out.extend([
        "",
        "Boundary Points:",
        "",
        "Point     X             Y          B.C. 1    Value        B.C. 2    Value",
    ])
    for nd in nodes:
        d1, d2 = ("_x", "_y") if nd["bc_dir"] == 0 else ("_n", "_t")
        map_bc = {0: "disp", 1: "trac", 2: "sprg"}
        bc1 = map_bc[nd["bck1"]] + d1
        bc2 = map_bc[nd["bck2"]] + d2
        out.append(
            f"{nd['id']:5d}  {nd['x']:12.5E}  {nd['y']:12.5E}  "
            f"{bc1:6s} {nd['val1_a']:14.7E}  {bc2:6s} {nd['val2_a']:14.7E}"
        )

    if int_pts:
        out.extend([
            "",
            "Internal Points:",
            "",
            "Point     X             Y",
        ])
        for i, (px, py) in enumerate(int_pts, start=1):
            out.append(f"{i:5d}  {px:12.5E}  {py:12.5E}")

    out.append("*" * 82)
    echo_path.parent.mkdir(parents=True, exist_ok=True)
    echo_path.write_text("\n".join(out) + "\n", encoding="utf-8")


def run_preproc(problem=1, aw: float = 0.4, s51_variant: str = "da0_01_r11",
                output_dir: pathlib.Path | None = None, plot_png: bool = True,
                dpi: int = 300) -> tuple[pathlib.Path, dict]:
    """
    Executa todo o pipeline de pré-processamento para um problema:
      .geo (Gmsh) -> .msh (meshio) -> modelo descontínuo (.dat) + gráfico (_mesh_bc.png).
    """
    data = input_data(problem=problem, aw=aw, s51_variant=s51_variant)
    file_name = data["file_name"]
    aid = data["analysis_id"]

    geo_path = MESH_DIR / f"{file_name}.geo"
    if not geo_path.exists():
        raise FileNotFoundError(f"Arquivo geométrico Gmsh não encontrado: {geo_path}")

    out_dir = pathlib.Path(output_dir) if output_dir else (ROOT_DIR / "inputs")
    out_dir.mkdir(parents=True, exist_ok=True)

    msh_path = MESH_DIR / f"{aid}.msh"
    generate_msh_from_geo(geo_path, msh_path, data.get("gmsh_params"))

    mesh = meshio.read(str(msh_path))
    model = build_discontinuous_model(data, mesh)

    dat_path = out_dir / f"{aid}.dat"
    write_fortran_dat_file(data, model, dat_path)

    if plot_png:
        echo_path = out_dir / f"{aid}_input_data.txt"
        write_input_data_echo_for_plot(data, model, echo_path)
        png_path = out_dir / f"{aid}_mesh_bc.png"
        parsed = parse_input_data(echo_path)
        xplot, yplot = comp_xyplot(parsed["elem_ends"], parsed["elem_incid"], parsed["b_nodes"])
        is_crack = detect_crack_elements(parsed["elem_ends"], parsed["elem_lambdas"], parsed["elem_incid"], parsed["b_nodes"])
        show_geometry(parsed, xplot, yplot, is_crack, png_path, dpi=dpi, show_nodes=True, show_plot=False)

    print(
        f"[preproc] {aid}: {len(model['elements'])} elementos quadráticos descontínuos, "
        f"{len(model['nodes'])} nós de colocação, {len(model['cracks'])} trinca(s) -> {dat_path}"
    )
    return dat_path, data


def main():
    parser = argparse.ArgumentParser(
        description="Pré-processador Gmsh (.geo -> .msh) + meshio para o Fracture 2D (Elementos Descontínuos)."
    )
    parser.add_argument(
        "problem",
        nargs="?",
        default="1",
        help="Número (1..8) ou nome do problema (ex.: Chambers_CMM_4, Dowling_11.6_variable_ampl, s521_single_edge, all).",
    )
    parser.add_argument("--aw", type=float, default=0.4, help="Razão a/w ou 2a/w para os problemas 7 e 8 (padrão: 0.4).")
    parser.add_argument("--s51-variant", default="da0_01_r11", choices=["da0_01_r11", "da0_0125_r105"], help="Esquema de avanço da Seção 5.1.")
    parser.add_argument("-o", "--output-dir", default=None, help="Diretório de saída para o arquivo .dat e a imagem .png (padrão: inputs/).")
    parser.add_argument("--no-plot", action="store_true", help="Não gera a figura <id>_mesh_bc.png.")
    parser.add_argument("--dpi", type=int, default=300, help="Resolução da imagem .png (padrão: 300).")
    args = parser.parse_args()

    out_dir = pathlib.Path(args.output_dir) if args.output_dir else None

    if str(args.problem).strip().lower() == "all":
        for p_id in range(1, 9):
            run_preproc(problem=p_id, aw=args.aw, s51_variant=args.s51_variant,
                        output_dir=out_dir, plot_png=not args.no_plot, dpi=args.dpi)
    else:
        run_preproc(problem=args.problem, aw=args.aw, s51_variant=args.s51_variant,
                    output_dir=out_dir, plot_png=not args.no_plot, dpi=args.dpi)


if __name__ == "__main__":
    main()
