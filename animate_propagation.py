#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
animate_propagation.py — Gerador de Animação MP4 da Propagação de Trincas (Dual-BEM)
================================================================================
Gera um vídeo `.mp4` (H.264 / yuv420p) animando passo a passo a propagação da(s)
trinca(s) calculada pelo solver Fortran `Fracture 2D` (Dual-BEM com elementos
quadráticos descontínuos).

Lê automaticamente todos os arquivos da pasta de resultados:
  - <id>-STEP-NNN.xml                 : Malha de contorno e deslocamentos por passo
  - <id>_geometry_crack_NNN.txt       : Coordenadas da(s) ponta(s) e comprimento a
  - <id>_results_crack_NNN.txt        : Fatores de Intensidade de Tensão (KI, KII, Keq)
  - <id>_fatigue_results_crack_NNN.txt: Ciclos acumulados N vs comprimento da trinca

Layouts disponíveis (`--layout`):
  - dashboard (padrão): Malha deformada global + Zoom na região da trinca +
                        Evolução sincronizada de KI/KII/Keq + Curva de fadiga N x Delta a
  - mesh              : Malha deformada global + Zoom na região da trinca
  - single            : Apenas a malha deformada global

Uso:
    python animate_propagation.py <diretorio_ou_prefixo> [opções]

Exemplos:
    python animate_propagation.py results/Chambers_CMM_4
    python animate_propagation.py results/Chambers_CMM_4 --scale 40 --fps 6
    python animate_propagation.py results/s521_single_edge --layout mesh --scale 80
    python animate_propagation.py results/Dowling_11.6_variable_ampl -o examples/dowling.mp4
================================================================================
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # Backend offscreen para renderização rápida e segura

import matplotlib.animation as animation
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter
import numpy as np

from plot_results import (
    find_simulation_files,
    parse_fatigue_results,
    parse_geometry_crack,
    parse_results_crack,
    parse_step_xml,
)


# ==============================================================================
# CONFIGURAÇÃO DO FFMPEG E ESTILO VISUAL
# ==============================================================================
def configure_ffmpeg() -> str:
    """
    Localiza o executável `ffmpeg` no PATH do sistema ou via `imageio_ffmpeg`
    e registra o caminho no Matplotlib (`animation.ffmpeg_path`).
    """
    ffmpeg_path = shutil.which("ffmpeg")
    if ffmpeg_path is None:
        try:
            import imageio_ffmpeg
            ffmpeg_path = imageio_ffmpeg.get_ffmpeg_exe()
        except Exception:
            ffmpeg_path = None

    if not ffmpeg_path or not os.path.exists(ffmpeg_path):
        raise RuntimeError(
            "Executável 'ffmpeg' não encontrado para gerar arquivos .mp4.\n"
            "Instale o ffmpeg no sistema ou instale o pacote Python: pip install imageio-ffmpeg"
        )

    plt.rcParams["animation.ffmpeg_path"] = ffmpeg_path
    return ffmpeg_path


def setup_animation_style() -> None:
    """Aplica estilo limpo e de alto contraste para vídeos didáticos e apresentações."""
    plt.rcParams.update({
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.titlesize": 12,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "legend.fontsize": 9,
        "figure.titlesize": 13.5,
        "grid.alpha": 0.35,
        "grid.linestyle": "--",
        "axes.facecolor": "#fbfcfe",
        "figure.facecolor": "white",
    })


# ==============================================================================
# RECONSTRUÇÃO DAS CURVAS DA MALHA XML E DETECÇÃO AUTOMÁTICA DAS FACES DA TRINCA
# ==============================================================================
def extract_step_number(xml_path: str) -> int:
    """Extrai o número inteiro NNN de `<prefix>-STEP-NNN.xml`."""
    m = re.search(r"-STEP-(\d+)\.xml$", os.path.basename(xml_path), re.IGNORECASE)
    return int(m.group(1)) if m else 0


def build_step_geometry(
    xml_data: dict,
    scale_factor: float,
    n_quad_pts: int = 11,
    model_span: float | None = None,
) -> dict:
    """
    Reconstrói os segmentos indeformados e deformados de um passo XML, separando
    automaticamente:
      - Contorno externo indeformado / deformado
      - Faces das trincas (indeformadas / deformadas)

    Nota sobre a convenção INSANE gravada em `src/io_mod.f90`:
      - `BL3`: incidência `[n1, n3, n2]` (início, fim, nó médio);
      - `BL2`: pseudo-elementos conectando os extremos geométricos aos nós
               descontínuos `n1` e `n3` do respectivo `BL3`.
    """
    nodes = xml_data["nodes"]
    elements = xml_data["elements"]

    # Determina a escala característica do modelo para tolerância de coincidência
    if model_span is None:
        xs = [nd["x"] for nd in nodes.values()]
        ys = [nd["y"] for nd in nodes.values()]
        model_span = max(max(xs) - min(xs), max(ys) - min(ys), 1e-12)
    tol = 1.0e-4 * model_span

    # 1. Identifica os elementos BL3 que possuem um par coincidente (faces opostas da trinca)
    bl3_indices = []
    bl3_mids = []
    for idx, el in enumerate(elements):
        inc = el["nodes"]
        if len(inc) == 3 and all(n in nodes for n in inc):
            nmid = inc[2]  # Convenção INSANE BL3: [n1, n3, n2]
            bl3_indices.append(idx)
            bl3_mids.append((nodes[nmid]["x"], nodes[nmid]["y"]))

    is_crack_elem = [False] * len(elements)
    crack_nodes = set()

    if bl3_mids:
        mids_arr = np.array(bl3_mids, dtype=float)
        # Matriz de distâncias entre pontos médios dos elementos BL3
        diff = mids_arr[:, None, :] - mids_arr[None, :, :]
        dist = np.hypot(diff[:, :, 0], diff[:, :, 1])
        np.fill_diagonal(dist, np.inf)
        has_twin = np.min(dist, axis=1) < tol

        for k, elem_idx in enumerate(bl3_indices):
            if has_twin[k]:
                is_crack_elem[elem_idx] = True
                crack_nodes.update(elements[elem_idx]["nodes"])

    # 2. Pseudo-elementos BL2 conectados a nós de um BL3 da trinca herdam is_crack = True
    for idx, el in enumerate(elements):
        inc = el["nodes"]
        if len(inc) == 2 and any(n in crack_nodes for n in inc):
            is_crack_elem[idx] = True

    # 3. Avalia polilinhas paramétricas de todos os elementos
    xi = np.linspace(-1.0, 1.0, n_quad_pts)
    N1 = 0.5 * xi * (xi - 1.0)
    Nmid = 1.0 - xi * xi
    N2 = 0.5 * xi * (xi + 1.0)

    ext_orig, ext_def = [], []
    crk_orig, crk_def = [], []

    for idx, el in enumerate(elements):
        inc = el["nodes"]
        if not all(n in nodes for n in inc):
            continue

        if len(inc) == 2:
            n1, n2 = inc
            p1, p2 = nodes[n1], nodes[n2]
            x0 = np.array([p1["x"], p2["x"]])
            y0 = np.array([p1["y"], p2["y"]])
            xd = np.array([p1["x"] + scale_factor * p1["ux"],
                           p2["x"] + scale_factor * p2["ux"]])
            yd = np.array([p1["y"] + scale_factor * p1["uy"],
                           p2["y"] + scale_factor * p2["uy"]])
        elif len(inc) == 3:
            # Convenção INSANE: [n1, n3, n2] -> n1 (início), n2 (fim), nmid (meio)
            n1, n2, nmid = inc
            p1, p2, pm = nodes[n1], nodes[n2], nodes[nmid]
            x0 = N1 * p1["x"] + Nmid * pm["x"] + N2 * p2["x"]
            y0 = N1 * p1["y"] + Nmid * pm["y"] + N2 * p2["y"]
            ux = N1 * p1["ux"] + Nmid * pm["ux"] + N2 * p2["ux"]
            uy = N1 * p1["uy"] + Nmid * pm["uy"] + N2 * p2["uy"]
            xd = x0 + scale_factor * ux
            yd = y0 + scale_factor * uy
        else:
            continue

        seg_orig = np.column_stack((x0, y0))
        seg_def = np.column_stack((xd, yd))

        if is_crack_elem[idx]:
            crk_orig.append(seg_orig)
            crk_def.append(seg_def)
        else:
            ext_orig.append(seg_orig)
            ext_def.append(seg_def)

    return {
        "ext_orig": ext_orig,
        "ext_def": ext_def,
        "crk_orig": crk_orig,
        "crk_def": crk_def,
    }


def fit_equal_limits(ax, xmin: float, xmax: float, ymin: float, ymax: float) -> None:
    """
    Ajusta `xlim` e `ylim` para preencher a caixa do subplot mantendo proporção
    geométrica exata 1:1 (`aspect='equal'`) sem cortar nenhuma parte do domínio.
    """
    fig_w, fig_h = ax.figure.get_size_inches()
    bbox = ax.get_position()
    box_aspect = (bbox.width * fig_w) / max(bbox.height * fig_h, 1e-9)

    dx = max(xmax - xmin, 1e-12)
    dy = max(ymax - ymin, 1e-12)
    cx = 0.5 * (xmin + xmax)
    cy = 0.5 * (ymin + ymax)

    if dx / dy < box_aspect:
        dx = dy * box_aspect
    else:
        dy = dx / box_aspect

    ax.set_xlim(cx - 0.5 * dx, cx + 0.5 * dx)
    ax.set_ylim(cy - 0.5 * dy, cy + 0.5 * dy)
    ax.set_aspect("equal", adjustable="box")


def compute_nice_scale_factor(xml_steps_data: list[dict], model_span: float) -> float:
    """
    Calcula automaticamente um fator de amplificação limpo (arredondado) para que o
    deslocamento máximo no passo final represente ~4% da dimensão característica da peça.
    """
    max_u = 0.0
    for xd in xml_steps_data:
        for nd in xd["nodes"].values():
            u_mag = float(np.hypot(nd["ux"], nd["uy"]))
            if u_mag > max_u:
                max_u = u_mag

    if max_u <= 1e-15:
        return 1.0

    raw = (0.045 * model_span) / max_u
    exp = np.floor(np.log10(raw))
    mant = raw / (10.0 ** exp)
    for nice in (1.0, 2.0, 2.5, 4.0, 5.0, 8.0, 10.0):
        if mant <= nice:
            return float(nice * (10.0 ** exp))
    return float(round(raw, 1))


# ==============================================================================
# CONSTRUÇÃO E RENDERIZAÇÃO DA ANIMAÇÃO MP4
# ==============================================================================
def create_propagation_animation(
    target_path: str,
    output_path: str | None = None,
    scale: float | None = None,
    fps: int = 5,
    dpi: int = 150,
    layout: str = "dashboard",
    hold_start: float = 0.8,
    hold_end: float = 2.0,
    follow_tip: bool = False,
    stride: int = 1,
) -> Path:
    """
    Carrega todos os passos da simulação em `target_path` e gera o vídeo `.mp4`
    (ou `.gif` caso a extensão de `output_path` seja `.gif`).
    """
    setup_animation_style()

    files_info = find_simulation_files(target_path)
    prefix = files_info["prefix"] or Path(target_path).name

    if not files_info["xml_steps"]:
        raise FileNotFoundError(
            f"Nenhum arquivo '<id>-STEP-NNN.xml' encontrado em '{files_info['dir']}'."
        )

    # Ordena os XMLs numericamente pelo índice do passo
    xml_entries = sorted(
        [(extract_step_number(p), p) for p in files_info["xml_steps"]],
        key=lambda item: item[0],
    )
    if stride > 1:
        selected = xml_entries[::stride]
        if selected[-1][0] != xml_entries[-1][0]:
            selected.append(xml_entries[-1])
        xml_entries = selected

    step_nums = [s for s, _ in xml_entries]
    n_steps = len(step_nums)

    print("=" * 72)
    print(" Animador de Propagação de Trincas — Fracture 2D (Dual-BEM)")
    print("=" * 72)
    print(f" Diretório : {files_info['dir']}")
    print(f" Prefixo   : {prefix}")
    print(f" Passos XML: {n_steps} quadros (Passo {step_nums[0]} até Passo {step_nums[-1]})")

    # Carrega arquivos auxiliares de trinca, FITs e fadiga
    geom_list = [parse_geometry_crack(f) for f in files_info["geometry_crack"]]
    res_list = [parse_results_crack(f) for f in files_info["results_crack"]]
    fat_list = [parse_fatigue_results(f) for f in files_info["fatigue_crack"]]

    # Carrega todos os arquivos XML da sequência
    print(" Lendo arquivos XML de todos os passos...")
    xml_steps_data = [parse_step_xml(p) for _, p in xml_entries]

    # Determina o bounding box global da geometria indeformada
    all_x0 = [nd["x"] for nd in xml_steps_data[0]["nodes"].values()]
    all_y0 = [nd["y"] for nd in xml_steps_data[0]["nodes"].values()]
    x_min0, x_max0 = min(all_x0), max(all_x0)
    y_min0, y_max0 = min(all_y0), max(all_y0)
    span_x = max(x_max0 - x_min0, 1e-12)
    span_y = max(y_max0 - y_min0, 1e-12)
    model_span = max(span_x, span_y)

    # Define fator de amplificação da deformada
    if scale is None:
        scale_factor = compute_nice_scale_factor(xml_steps_data, model_span)
    else:
        scale_factor = float(scale)

    show_deformed = scale_factor > 0.0
    print(f" Fator de amplificação da deformada: {scale_factor:g}x")

    # Pré-computa as polilinhas (contorno externo e faces da trinca) de cada passo
    step_geoms = [
        build_step_geometry(xd, scale_factor, n_quad_pts=11, model_span=model_span)
        for xd in xml_steps_data
    ]

    # Bounding box global considerando a malha deformada do passo final para evitar cortes
    all_x_glob = list(all_x0)
    all_y_glob = list(all_y0)
    if show_deformed:
        for seg in step_geoms[-1]["ext_def"] + step_geoms[-1]["crk_def"]:
            all_x_glob.extend(seg[:, 0].tolist())
            all_y_glob.extend(seg[:, 1].tolist())
    gx_min, gx_max = min(all_x_glob), max(all_x_glob)
    gy_min, gy_max = min(all_y_glob), max(all_y_glob)
    g_pad_x = 0.08 * max(gx_max - gx_min, model_span * 0.2)
    g_pad_bot = 0.08 * max(gy_max - gy_min, model_span * 0.2)
    g_pad_top = 0.22 * max(gy_max - gy_min, model_span * 0.2)

    # Coleta coordenadas das pontas das trincas (indeformadas e deformadas) ao longo de todos os passos
    # e calcula janela fixa de zoom sobre a região da(s) trinca(s)
    crack_tip_coords_per_step: list[list[tuple[float, float]]] = []
    crack_xs, crack_ys = [], []

    # Determina se a segunda extremidade de cada trinca também propaga (trinca interna)
    crack_two_tips = []
    for c_idx, g in enumerate(geom_list):
        has_t2 = False
        if c_idx < len(res_list) and res_list[c_idx]["has_tip_2"]:
            has_t2 = True
        elif len(g["x_init"]) > 1:
            moved = np.hypot(g["x_init"][-1] - g["x_init"][0], g["y_init"][-1] - g["y_init"][0])
            has_t2 = moved > 1e-6 * model_span
        crack_two_tips.append(has_t2)

    def deformed_tip_position(tx: float, ty: float, xml_data: dict) -> tuple[float, float]:
        """Calcula a posição deformada da ponta da trinca usando os nós XML coincidentes."""
        if not show_deformed:
            return tx, ty
        nodes = xml_data["nodes"]
        tol_tip = 1.0e-3 * model_span
        uxs, uys = [], []
        best_d = np.inf
        best_u = (0.0, 0.0)
        for nd in nodes.values():
            d = float(np.hypot(nd["x"] - tx, nd["y"] - ty))
            if d < best_d:
                best_d = d
                best_u = (nd["ux"], nd["uy"])
            if d <= tol_tip:
                uxs.append(nd["ux"])
                uys.append(nd["uy"])
        if uxs:
            ux_m = float(np.mean(uxs))
            uy_m = float(np.mean(uys))
        else:
            ux_m, uy_m = best_u
        return tx + scale_factor * ux_m, ty + scale_factor * uy_m

    for step_idx, s_num in enumerate(step_nums):
        tips_here = []
        xd_cur = xml_steps_data[step_idx]
        for c_idx, g in enumerate(geom_list):
            matches = np.where(g["step"] == s_num)[0]
            if len(matches) == 0:
                continue
            k = int(matches[0])
            tx1, ty1 = float(g["x_final"][k]), float(g["y_final"][k])
            tips_here.append(deformed_tip_position(tx1, ty1, xd_cur))
            crack_xs.extend([float(g["x_init"][k]), tx1])
            crack_ys.extend([float(g["y_init"][k]), ty1])
            if crack_two_tips[c_idx]:
                tx2, ty2 = float(g["x_init"][k]), float(g["y_init"][k])
                tips_here.append(deformed_tip_position(tx2, ty2, xd_cur))
        crack_tip_coords_per_step.append(tips_here)

    if crack_xs and crack_ys:
        # Inclui também as faces deformadas da trinca no passo final na janela de zoom
        zx_vals = list(crack_xs)
        zy_vals = list(crack_ys)
        for seg in step_geoms[-1]["crk_def"]:
            zx_vals.extend(seg[:, 0].tolist())
            zy_vals.extend(seg[:, 1].tolist())
        z_span = max(max(zx_vals) - min(zx_vals), max(zy_vals) - min(zy_vals), 0.12 * model_span)
        z_pad = 0.22 * z_span
        zx_min, zx_max = min(zx_vals) - z_pad, max(zx_vals) + z_pad
        zy_min, zy_max = min(zy_vals) - z_pad, max(zy_vals) + z_pad
        z_half = 0.5 * z_span + z_pad
    else:
        cx_mid = 0.5 * (x_min0 + x_max0)
        cy_mid = 0.5 * (y_min0 + y_max0)
        z_half = 0.25 * model_span
        zx_min, zx_max = cx_mid - z_half, cx_mid + z_half
        zy_min, zy_max = cy_mid - z_half, cy_mid + z_half

    # Prepara séries sincronizadas da trinca principal (Trinca 1) para os gráficos laterais e HUD
    has_sif = len(res_list) > 0 and len(res_list[0]["step"]) > 0
    has_fatigue = len(fat_list) > 0 and len(fat_list[0]["length"]) > 1

    # Mapeia passo -> (a, da, N, KI, KII, Keq)
    hud_info: dict[int, dict[str, float]] = {}
    if has_sif:
        r0 = res_list[0]
        a0_ref = float(r0["length"][0])
        for i, s_val in enumerate(r0["step"]):
            a_val = float(r0["length"][i])
            n_val = 0.0
            if has_fatigue:
                f0 = fat_list[0]
                n_val = float(np.interp(a_val, f0["length"], f0["n1_total"]))
            hud_info[int(s_val)] = {
                "a": a_val,
                "da": a_val - a0_ref,
                "N": n_val,
                "ki": float(r0["ki_1"][i]),
                "kii": float(r0["kii_1"][i]),
                "keq": float(r0["keq_1"][i]),
            }

    # Determina se layout dashboard tem dados suficientes; se não, simplifica
    effective_layout = layout
    if effective_layout == "dashboard" and not has_sif:
        effective_layout = "mesh"

    # ==========================================================================
    # MONTAGEM DA FIGURA E SUBPLOTS
    # ==========================================================================
    ax_zoom = None
    ax_sif = None
    ax_fat = None

    if effective_layout == "dashboard":
        fig = plt.figure(figsize=(14.4, 8.4), dpi=dpi)
        gs = GridSpec(
            3, 2, figure=fig,
            width_ratios=[1.22, 1.08],
            height_ratios=[1.15, 0.90, 0.90] if has_fatigue else [1.3, 1.0, 0.01],
            hspace=0.48, wspace=0.20,
            left=0.055, right=0.975, top=0.90, bottom=0.07,
        )
        ax_glob = fig.add_subplot(gs[:, 0])
        ax_zoom = fig.add_subplot(gs[0, 1])
        if has_fatigue:
            ax_sif = fig.add_subplot(gs[1, 1])
            ax_fat = fig.add_subplot(gs[2, 1])
        else:
            ax_sif = fig.add_subplot(gs[1:, 1])
    elif effective_layout == "mesh":
        fig, (ax_glob, ax_zoom) = plt.subplots(
            1, 2, figsize=(13.6, 6.8), dpi=dpi,
            gridspec_kw={"wspace": 0.18, "left": 0.06, "right": 0.97, "top": 0.87, "bottom": 0.09},
        )
    else:  # single
        fig, ax_glob = plt.subplots(
            1, 1, figsize=(9.6, 7.6), dpi=dpi,
            gridspec_kw={"left": 0.09, "right": 0.96, "top": 0.88, "bottom": 0.09},
        )

    fig.suptitle(
        f"Propagação de Trinca por Fadiga (Dual-BEM) — {prefix}",
        fontweight="bold",
        color="#1a2433",
    )

    # --- Configuração do painel global ---
    fit_equal_limits(ax_glob, gx_min - g_pad_x, gx_max + g_pad_x, gy_min - g_pad_bot, gy_max + g_pad_top)
    ax_glob.set_xlabel("X (unid. comp.)")
    ax_glob.set_ylabel("Y (unid. comp.)")
    ax_glob.set_title(
        f"Configuração Global (Deformada {scale_factor:g}x)" if show_deformed else "Configuração Global da Malha",
        fontweight="bold",
    )
    ax_glob.grid(True)

    # Coleções de linhas para atualização ultrarrápida quadro a quadro
    lc_glob_ext_orig = LineCollection(
        [], colors="#94a3b8", linewidths=0.95, linestyles="--",
        capstyle="round", joinstyle="round", alpha=0.65, zorder=2,
    )
    lc_glob_crk_orig = LineCollection(
        [], colors="#64748b", linewidths=1.0, linestyles=":",
        capstyle="round", joinstyle="round", alpha=0.65, zorder=3,
    )
    lc_glob_ext_def = LineCollection(
        [], colors="#0f2942", linewidths=1.5,
        capstyle="round", joinstyle="round", alpha=0.96, zorder=4,
    )
    lc_glob_crk_def = LineCollection(
        [], colors="#dc2626", linewidths=2.2,
        capstyle="round", joinstyle="round", alpha=1.0, zorder=5,
    )

    for lc in (lc_glob_ext_orig, lc_glob_crk_orig, lc_glob_ext_def, lc_glob_crk_def):
        ax_glob.add_collection(lc)

    # Marcador brilhante na(s) ponta(s) ativa(s) da trinca
    (tip_glob_marker,) = ax_glob.plot(
        [], [], "o", color="#f59e0b", markeredgecolor="#7f1d1d",
        markeredgewidth=1.4, markersize=7.5, zorder=7, label="Ponta da Trinca",
    )

    legend_handles = []
    if show_deformed:
        legend_handles.append(Line2D([0], [0], color="#94a3b8", lw=1.0, ls="--", label="Indeformada"))
        legend_handles.append(Line2D([0], [0], color="#0f2942", lw=1.5, label=f"Contorno ({scale_factor:g}x)"))
        legend_handles.append(Line2D([0], [0], color="#dc2626", lw=2.2, label=f"Faces da Trinca ({scale_factor:g}x)"))
    else:
        legend_handles.append(Line2D([0], [0], color="#0f2942", lw=1.5, label="Contorno Externo"))
        legend_handles.append(Line2D([0], [0], color="#dc2626", lw=2.2, label="Trinca"))
    legend_handles.append(tip_glob_marker)
    ax_glob.legend(handles=legend_handles, loc="lower left", framealpha=0.92)

    # Caixa de telemetria (HUD) no canto superior esquerdo do painel global
    hud_text = ax_glob.text(
        0.025, 0.975, "", transform=ax_glob.transAxes,
        va="top", ha="left", fontsize=9.0, family="monospace",
        bbox=dict(boxstyle="round,pad=0.45", facecolor="white", edgecolor="#cbd5e1", alpha=0.94),
        zorder=10,
    )

    # --- Configuração do painel de Zoom na trinca ---
    lc_zoom_ext_orig = lc_zoom_crk_orig = lc_zoom_ext_def = lc_zoom_crk_def = None
    tip_zoom_marker = None
    if ax_zoom is not None:
        fit_equal_limits(ax_zoom, zx_min, zx_max, zy_min, zy_max)
        ax_zoom.set_xlabel("X")
        ax_zoom.set_ylabel("Y")
        ax_zoom.set_title(
            f"Zoom na Região da Trinca ({scale_factor:g}x)" if show_deformed else "Zoom na Trajetória da Trinca",
            fontweight="bold",
        )
        ax_zoom.grid(True)

        lc_zoom_ext_orig = LineCollection(
            [], colors="#94a3b8", linewidths=1.0, linestyles="--",
            capstyle="round", joinstyle="round", alpha=0.65, zorder=2,
        )
        lc_zoom_crk_orig = LineCollection(
            [], colors="#64748b", linewidths=1.15, linestyles=":",
            capstyle="round", joinstyle="round", alpha=0.68, zorder=3,
        )
        lc_zoom_ext_def = LineCollection(
            [], colors="#0f2942", linewidths=1.65,
            capstyle="round", joinstyle="round", alpha=0.96, zorder=4,
        )
        lc_zoom_crk_def = LineCollection(
            [], colors="#dc2626", linewidths=2.6,
            capstyle="round", joinstyle="round", alpha=1.0, zorder=5,
        )
        for lc in (lc_zoom_ext_orig, lc_zoom_crk_orig, lc_zoom_ext_def, lc_zoom_crk_def):
            ax_zoom.add_collection(lc)

        (tip_zoom_marker,) = ax_zoom.plot(
            [], [], "o", color="#f59e0b", markeredgecolor="#7f1d1d",
            markeredgewidth=1.4, markersize=7.5, zorder=7,
        )

    # --- Configuração do painel de Fatores de Intensidade de Tensão (KI, KII, Keq) ---
    sif_lines = {}
    sif_dots = {}
    if ax_sif is not None and has_sif:
        r0 = res_list[0]
        s_arr = r0["step"]
        ki_arr = r0["ki_1"]
        kii_arr = r0["kii_1"]
        keq_arr = r0["keq_1"]

        # Curvas fantasma completas ao fundo como guia visual
        ax_sif.plot(s_arr, ki_arr, "--", color="#1d4ed8", lw=0.9, alpha=0.25)
        ax_sif.plot(s_arr, kii_arr, "--", color="#d97706", lw=0.9, alpha=0.25)
        ax_sif.plot(s_arr, keq_arr, "--", color="#dc2626", lw=1.0, alpha=0.25)

        (sif_lines["ki"],) = ax_sif.plot([], [], "-", color="#1d4ed8", lw=1.8, label=r"$K_I$")
        (sif_lines["kii"],) = ax_sif.plot([], [], "-", color="#d97706", lw=1.6, label=r"$K_{II}$")
        (sif_lines["keq"],) = ax_sif.plot([], [], "-", color="#dc2626", lw=2.0, label=r"$K_{\mathrm{eq}}$")

        (sif_dots["ki"],) = ax_sif.plot([], [], "o", color="#1d4ed8", ms=5.5, zorder=6)
        (sif_dots["kii"],) = ax_sif.plot([], [], "s", color="#d97706", ms=5.0, zorder=6)
        (sif_dots["keq"],) = ax_sif.plot([], [], "^", color="#dc2626", ms=6.0, zorder=6)

        y_all = np.concatenate([ki_arr, kii_arr, keq_arr])
        y_lo, y_hi = float(np.min(y_all)), float(np.max(y_all))
        y_pad = max(0.08 * (y_hi - y_lo), 1.0)
        ax_sif.set_xlim(float(s_arr[0]) - 0.5, float(s_arr[-1]) + 0.5)
        ax_sif.set_ylim(y_lo - y_pad, y_hi + y_pad)
        ax_sif.set_xlabel("Passo de Propagação")
        ax_sif.set_ylabel(r"FIT (MPa$\sqrt{\mathrm{m}}$ / unid.)")
        ax_sif.set_title("Evolução dos Fatores de Intensidade de Tensão", fontweight="bold")
        ax_sif.legend(loc="best", ncol=3, framealpha=0.9)
        ax_sif.grid(True)

    # --- Configuração do painel de Vida em Fadiga (Delta a vs N) ---
    fat_line = fat_dot = None
    if ax_fat is not None and has_fatigue:
        f0 = fat_list[0]
        N_full = f0["n1_total"]
        da_full = f0["length"] - f0["length"][0]

        ax_fat.plot(N_full, da_full, "--", color="#0f766e", lw=1.0, alpha=0.25)
        (fat_line,) = ax_fat.plot([], [], "-", color="#0f766e", lw=2.0, label=r"$\Delta a(N)$")
        (fat_dot,) = ax_fat.plot([], [], "o", color="#0f766e", markeredgecolor="black", ms=6.0, zorder=6)

        fmt = ScalarFormatter(useMathText=True)
        fmt.set_powerlimits((-2, 4))
        ax_fat.xaxis.set_major_formatter(fmt)

        n_max = max(float(np.max(N_full)), 1.0)
        da_max = max(float(np.max(da_full)), 1e-6)
        ax_fat.set_xlim(-0.03 * n_max, 1.05 * n_max)
        ax_fat.set_ylim(-0.05 * da_max, 1.08 * da_max)
        ax_fat.set_xlabel(r"Número de Ciclos $N$")
        ax_fat.set_ylabel(r"Avanço $\Delta a = a - a_0$")
        ax_fat.set_title("Crescimento da Trinca por Fadiga", fontweight="bold")
        ax_fat.grid(True)

    # ==========================================================================
    # SEQUÊNCIA DE QUADROS (com pausa inicial e final configuráveis)
    # ==========================================================================
    n_hold_start = max(0, int(round(hold_start * fps)))
    n_hold_end = max(0, int(round(hold_end * fps)))
    frame_indices = [0] * n_hold_start + list(range(n_steps)) + [n_steps - 1] * n_hold_end

    def update(frame_idx: int):
        s_num = step_nums[frame_idx]
        sg = step_geoms[frame_idx]

        # Atualiza coleções do painel global
        if show_deformed:
            lc_glob_ext_orig.set_segments(sg["ext_orig"])
            lc_glob_crk_orig.set_segments(sg["crk_orig"])
            lc_glob_ext_def.set_segments(sg["ext_def"])
            lc_glob_crk_def.set_segments(sg["crk_def"])
        else:
            lc_glob_ext_orig.set_segments([])
            lc_glob_crk_orig.set_segments([])
            lc_glob_ext_def.set_segments(sg["ext_orig"])
            lc_glob_crk_def.set_segments(sg["crk_orig"])

        # Pontas ativas da trinca
        tips = crack_tip_coords_per_step[frame_idx]
        if tips:
            tx = [p[0] for p in tips]
            ty = [p[1] for p in tips]
            tip_glob_marker.set_data(tx, ty)
        else:
            tip_glob_marker.set_data([], [])

        # Atualiza painel de zoom
        if ax_zoom is not None:
            if show_deformed:
                lc_zoom_ext_orig.set_segments(sg["ext_orig"])
                lc_zoom_crk_orig.set_segments(sg["crk_orig"])
                lc_zoom_ext_def.set_segments(sg["ext_def"])
                lc_zoom_crk_def.set_segments(sg["crk_def"])
            else:
                lc_zoom_ext_orig.set_segments([])
                lc_zoom_crk_orig.set_segments([])
                lc_zoom_ext_def.set_segments(sg["ext_orig"])
                lc_zoom_crk_def.set_segments(sg["crk_orig"])

            if tips and tip_zoom_marker is not None:
                tip_zoom_marker.set_data([p[0] for p in tips], [p[1] for p in tips])
                if follow_tip:
                    cx, cy = tips[0]
                    r_follow = max(0.18 * model_span, z_half * 0.55)
                    fit_equal_limits(ax_zoom, cx - r_follow, cx + r_follow, cy - r_follow, cy + r_follow)

        # Atualiza HUD
        lines_hud = [f"Passo : {s_num:3d} / {step_nums[-1]:d}"]
        if s_num in hud_info:
            info = hud_info[s_num]
            lines_hud.append(f"a     : {info['a']:.4g} (Δa={info['da']:.4g})")
            if has_fatigue:
                lines_hud.append(f"Ciclos: {info['N']:,.1f}")
            lines_hud.append(f"K_I   : {info['ki']:.2f}")
            lines_hud.append(f"K_II  : {info['kii']:.2f}")
            lines_hud.append(f"K_eq  : {info['keq']:.2f}")
        hud_text.set_text("\n".join(lines_hud))

        # Atualiza curvas de FITs
        if ax_sif is not None and has_sif:
            r0 = res_list[0]
            mask = r0["step"] <= s_num
            s_cur = r0["step"][mask]
            ki_cur = r0["ki_1"][mask]
            kii_cur = r0["kii_1"][mask]
            keq_cur = r0["keq_1"][mask]
            sif_lines["ki"].set_data(s_cur, ki_cur)
            sif_lines["kii"].set_data(s_cur, kii_cur)
            sif_lines["keq"].set_data(s_cur, keq_cur)
            if len(s_cur) > 0:
                sif_dots["ki"].set_data([s_cur[-1]], [ki_cur[-1]])
                sif_dots["kii"].set_data([s_cur[-1]], [kii_cur[-1]])
                sif_dots["keq"].set_data([s_cur[-1]], [keq_cur[-1]])

        # Atualiza curva de fadiga (Delta a vs N)
        if ax_fat is not None and has_fatigue and s_num in hud_info:
            f0 = fat_list[0]
            a_cur = hud_info[s_num]["a"]
            da_cur = hud_info[s_num]["da"]
            n_cur = hud_info[s_num]["N"]
            mask_f = f0["length"] <= a_cur + 1e-12
            n_trace = list(f0["n1_total"][mask_f])
            da_trace = list(f0["length"][mask_f] - f0["length"][0])
            if not n_trace or abs(n_trace[-1] - n_cur) > 1e-9:
                n_trace.append(n_cur)
                da_trace.append(da_cur)
            fat_line.set_data(n_trace, da_trace)
            fat_dot.set_data([n_cur], [da_cur])

        return ()

    # ==========================================================================
    # EXPORTAÇÃO DO VÍDEO (.MP4 ou .GIF)
    # ==========================================================================
    if output_path is None:
        out_dir = Path(__file__).resolve().parent / "examples"
        out_dir.mkdir(parents=True, exist_ok=True)
        out_file = out_dir / f"{prefix}_propagation.mp4"
    else:
        out_file = Path(output_path).resolve()
        out_file.parent.mkdir(parents=True, exist_ok=True)

    ext = out_file.suffix.lower()
    if ext == "":
        out_file = out_file.with_suffix(".mp4")
        ext = ".mp4"

    anim = animation.FuncAnimation(
        fig,
        update,
        frames=frame_indices,
        interval=1000.0 / max(fps, 1),
        blit=False,
    )

    print(f" Renderizando {len(frame_indices)} quadros a {fps} fps (DPI={dpi}, layout='{layout}')...")
    if ext == ".gif":
        writer = animation.PillowWriter(fps=fps)
        anim.save(str(out_file), writer=writer, dpi=dpi)
    else:
        ffmpeg_bin = configure_ffmpeg()
        print(f" Encoder FFmpeg: {ffmpeg_bin}")
        writer = animation.FFMpegWriter(
            fps=fps,
            codec="libx264",
            bitrate=2400,
            extra_args=[
                "-pix_fmt", "yuv420p",
                "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2",
                "-movflags", "+faststart",
            ],
        )
        anim.save(str(out_file), writer=writer, dpi=dpi)

    plt.close(fig)
    size_kb = out_file.stat().st_size / 1024.0
    print(f" [OK] Animação salva com sucesso: {out_file} ({size_kb:.1f} KB)")
    print("=" * 72)
    return out_file


# ==============================================================================
# INTERFACE DE LINHA DE COMANDO (CLI)
# ==============================================================================
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Gera uma animação em vídeo (.mp4) da propagação da trinca do Fracture 2D (Dual-BEM).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemplos de uso:
  python animate_propagation.py results/Chambers_CMM_4
  python animate_propagation.py results/Chambers_CMM_4 --scale 40 --fps 6
  python animate_propagation.py results/s521_single_edge --layout mesh --fps 5
  python animate_propagation.py results/Dowling_11.6_variable_ampl -o examples/dowling.mp4
        """,
    )
    parser.add_argument(
        "target",
        help="Diretório de resultados, prefixo ou arquivo da simulação (ex.: results/Chambers_CMM_4).",
    )
    parser.add_argument(
        "-o", "--output",
        default=None,
        help="Caminho do arquivo .mp4 (ou .gif) de saída (padrão: examples/<prefix>_propagation.mp4).",
    )
    parser.add_argument(
        "--scale",
        type=float,
        default=None,
        help="Fator de amplificação da malha deformada (padrão: automático; use 0 para malha indeformada).",
    )
    parser.add_argument(
        "--fps",
        type=int,
        default=5,
        help="Quadros por segundo do vídeo (padrão: 5).",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="Resolução em DPI de cada quadro (padrão: 150).",
    )
    parser.add_argument(
        "--layout",
        choices=["dashboard", "mesh", "single"],
        default="dashboard",
        help="Layout da animação: 'dashboard' (padrão), 'mesh' (global + zoom) ou 'single'.",
    )
    parser.add_argument(
        "--hold-start",
        type=float,
        default=0.8,
        help="Tempo (segundos) de pausa no passo inicial 0 (padrão: 0.8).",
    )
    parser.add_argument(
        "--hold-end",
        type=float,
        default=2.0,
        help="Tempo (segundos) de pausa no passo final (padrão: 2.0).",
    )
    parser.add_argument(
        "--follow-tip",
        action="store_true",
        help="Faz a janela de zoom acompanhar dinamicamente a ponta da trinca.",
    )
    parser.add_argument(
        "--stride",
        type=int,
        default=1,
        help="Salto entre passos XML (padrão: 1 = todos os passos).",
    )

    args = parser.parse_args()
    create_propagation_animation(
        target_path=args.target,
        output_path=args.output,
        scale=args.scale,
        fps=args.fps,
        dpi=args.dpi,
        layout=args.layout,
        hold_start=args.hold_start,
        hold_end=args.hold_end,
        follow_tip=args.follow_tip,
        stride=args.stride,
    )


if __name__ == "__main__":
    main()
