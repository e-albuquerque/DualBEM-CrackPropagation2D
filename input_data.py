#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
input_data.py — Definição dos Problemas de Fratura 2D e Propagação por Fadiga (Dual-BEM)
=========================================================================================

Segue o padrão modular dos projetos BEM2D_elastostatic e BEM3D_elastostatic:
cada problema referencia um arquivo `.geo` (Gmsh) na pasta `meshes/` e associa
as condições de contorno e trincas diretamente aos nomes dos `Physical Curve`
definidos no Gmsh.

Problemas disponíveis (por número 1..8 ou pelo identificador string):
  1 / "Chambers_CMM_4"           : Seção 5.3.2 — Corpo de prova CMM4 (alpha = 110°, theta = -4°)
  2 / "Dowling_11.6_variable_ampl": Seção 5.2.3 — Placa AISI 4340 sob amplitude variável (Rainflow)
  3 / "s521_single_edge"         : Seção 5.2.1 — Placa AISI 4340 com 1 trinca de borda (R = 1/3)
  4 / "s522_double_edge"         : Seção 5.2.2 — Placa AISI 4340 com 2 trincas simétricas de borda
  5 / "s51_internal_crack"       : Seção 5.1   — Placa retangular 3x1 com trinca central (Isida)
  6 / "Chambers_CMM_90"          : Seção 5.3.1 — Corpo de prova CMM1/2/3 (alpha = 90°, theta = -3°)
  7 / "s3131_edge_crack"         : Seção 3.1.3.1 — Placa quadrada 10x10 com trinca de borda (Modo I)
  8 / "s3132_inclined_crack"     : Seção 3.1.3.2 — Placa 10x20 com trinca central a 45° (Modo Misto)
"""

import math


def input_data(problem=1, aw: float = 0.4, s51_variant: str = "da0_01_r11") -> dict:
    """
    Retorna o dicionário completo de parâmetros físicos, geométricos, condições
    de contorno (associadas às Physical Curves do Gmsh), trincas e parâmetros
    de fadiga para o problema selecionado.
    """
    # Mapeamento de nomes amigáveis para o índice numérico do problema
    name_map = {
        "chambers_cmm_4": 1,
        "chambers_cmm4": 1,
        "cmm4": 1,
        "dowling_11.6_variable_ampl": 2,
        "dowling_11_6": 2,
        "dowling": 2,
        "s521_single_edge": 3,
        "s521": 3,
        "s522_double_edge": 4,
        "s522": 4,
        "s51_internal_crack": 5,
        "s51": 5,
        "chambers_cmm_90": 6,
        "s531_cmm_90": 6,
        "cmm90": 6,
        "s3131_edge_crack": 7,
        "s3131": 7,
        "s3132_inclined_crack": 8,
        "s3132": 8,
    }
    if isinstance(problem, str):
        key = problem.strip().lower()
        if key.isdigit():
            problem = int(key)
        elif key in name_map:
            problem = name_map[key]
        else:
            raise ValueError(f"Problema desconhecido: '{problem}'. Opções: {list(name_map.keys())} ou 1..8.")

    lamb_default = 2.0 / 3.0

    # =========================================================================
    # PROBLEMA 1: Seção 5.3.2 — Corpo de Prova Compact Mixed Mode (CMM4, 110°)
    # =========================================================================
    if problem == 1:
        return {
            "problem": 1,
            "title": "Compact Mixed Mode Specimen CMM4 (alpha=110 deg, theta=-4 deg)",
            "analysis_id": "Chambers_CMM_4",
            "file_name": "chambers_cmm4",
            "gmsh_params": {},
            "E": 200000.0,
            "nu": 0.25,
            "plane_state": 1,  # 1 = Estado Plano de Tensões (EPT)
            "npow": 9,         # Ordem da expansão polinomial radial de Gao
            "lamb": lamb_default,
            "bound_cond": {
                "free_boundary": {
                    "bc_dir": 0,
                    "type_x": "tx", "value_x": 0.0,
                    "type_y": "ty", "value_y": 0.0,
                },
                "right_symmetry": {
                    "bc_dir": 0,
                    "type_x": "ux", "value_x": 0.0,
                    "type_y": "ty", "value_y": 0.0,
                },
                "right_hole_supp": {
                    "bc_dir": 0,
                    "type_x": "tx", "value_x": 0.0,
                    "type_y": "uy", "value_y": 0.0,
                },
                "left_hole_q3": {
                    "bc_dir": 0,
                    "type_x": "tx", "value_x": (54.862812167336, -54.862812167336, 0.0, 0.0),
                    "type_y": "ty", "value_y": (0.0, 150.734337597536, 0.0, 0.0),
                },
                "left_hole_q4": {
                    "bc_dir": 0,
                    "type_x": "tx", "value_x": (0.0, 54.862812167336, 0.0, 0.0),
                    "type_y": "ty", "value_y": (150.734337597536, -150.734337597536, 0.0, 0.0),
                },
            },
            "cracks": [
                {"physical_group": "crack_1", "tips_num": 1},
            ],
            "internal_points_rects": [
                {"x1": 2.0, "y1": 2.0,   "x2": 37.105, "y2": 18.74, "nx": 8, "ny": 4},
                {"x1": 2.0, "y1": 32.94, "x2": 37.105, "y2": 56.74, "nx": 8, "ny": 4},
            ],
            "sif": {
                "k_key": 2, "method_key": 3,
                "nbe": 20, "nbn": 3, "lb1": lamb_default, "lb2": lamb_default,
                "np": 20, "initial_node": 2,
            },
            "fatigue": {
                "k_eff_key": 2,
                "da0": 0.1, "r_da": 1.1, "da_max": 1.0, "steps": 38,
                "R": 0.0625, "C0": 2.3459e-08, "m": 1.33, "gamma": 1.0,
                "K_th": 1.0, "K_IC": 1500.0,
            },
            "variable_amplitude": {
                "enabled": False,
                "alphas": [],
            },
        }

    # =========================================================================
    # PROBLEMA 2: Seção 5.2.3 — Placa com Trinca Central e Amplitude Variável
    # =========================================================================
    elif problem == 2:
        # Histórico de reversões de carga (Exemplo 11.6 de Dowling):
        # 101 ciclos entre 2.4 e 1.8, 25 ciclos entre 2.0 e 1.0, 40 ciclos entre 1.5 e 0.0
        alphas_dowling = [2.4, 1.8] * 101 + [2.0, 1.0] * 25 + [1.5, 0.0] * 40
        return {
            "problem": 2,
            "title": "Central Cracked Plate under Variable Amplitude Loading (Dowling 11.6)",
            "analysis_id": "Dowling_11.6_variable_ampl",
            "file_name": "dowling_11_6",
            "gmsh_params": {},
            "E": 207000.0,
            "nu": 0.30,
            "plane_state": 1,
            "npow": 3,
            "lamb": lamb_default,
            "bound_cond": {
                "bottom": {
                    "bc_dir": 0,
                    "type_x": "ux", "value_x": 0.0,
                    "type_y": "uy", "value_y": 0.0,
                },
                "right": {
                    "bc_dir": 0,
                    "type_x": "tx", "value_x": 0.0,
                    "type_y": "ty", "value_y": 0.0,
                },
                "top": {
                    "bc_dir": 0,
                    "type_x": "tx", "value_x": 0.0,
                    "type_y": "ty", "value_y": 219.298245614035,
                },
                "left": {
                    "bc_dir": 0,
                    "type_x": "tx", "value_x": 0.0,
                    "type_y": "ty", "value_y": 0.0,
                },
            },
            "cracks": [
                {"physical_group": "crack_1", "tips_num": 2},
            ],
            "internal_points_rects": [
                {"x1": 0.0038, "y1": 0.0038, "x2": 0.0722, "y2": 0.1482, "nx": 5, "ny": 10},
            ],
            "sif": {
                "k_key": 2, "method_key": 3,
                "nbe": 20, "nbn": 3, "lb1": lamb_default, "lb2": lamb_default,
                "np": 20, "initial_node": 2,
            },
            "fatigue": {
                "k_eff_key": 2,
                "da0": 0.000097, "r_da": 1.097, "da_max": 1.0, "steps": 46,
                "R": 0.3333333, "C0": 5.11e-13, "m": 3.24, "gamma": 0.420,
                "K_th": 1.0, "K_IC": 130.0,
            },
            "variable_amplitude": {
                "enabled": True,
                "alphas": alphas_dowling,
            },
        }

    # =========================================================================
    # PROBLEMA 3: Seção 5.2.1 — Placa AISI 4340 com Uma Trinca de Borda
    # =========================================================================
    elif problem == 3:
        return {
            "problem": 3,
            "title": "Single Edge Cracked Plate AISI 4340 (Sec. 5.2.1)",
            "analysis_id": "s521_single_edge",
            "file_name": "s521_single_edge",
            "gmsh_params": {},
            "E": 207000.0,
            "nu": 0.30,
            "plane_state": 1,
            "npow": 3,
            "lamb": lamb_default,
            "bound_cond": {
                "bottom": {
                    "bc_dir": 0,
                    "type_x": "tx", "value_x": 0.0,
                    "type_y": "ty", "value_y": -526.0,
                },
                "free_sides": {
                    "bc_dir": 0,
                    "type_x": "tx", "value_x": 0.0,
                    "type_y": "ty", "value_y": 0.0,
                },
                "right_support": {
                    "bc_dir": 0,
                    "type_x": "ux", "value_x": 0.0,
                    "type_y": "uy", "value_y": 0.0,
                },
                "top": {
                    "bc_dir": 0,
                    "type_x": "tx", "value_x": 0.0,
                    "type_y": "ty", "value_y": 526.0,
                },
            },
            "cracks": [
                {"physical_group": "crack_1", "tips_num": 1},
            ],
            "internal_points_rects": [],
            "sif": {
                "k_key": 2, "method_key": 3,
                "nbe": 20, "nbn": 3, "lb1": lamb_default, "lb2": lamb_default,
                "np": 20, "initial_node": 2,
            },
            "fatigue": {
                "k_eff_key": 2,
                "da0": 0.0001, "r_da": 1.1, "da_max": 1.0, "steps": 38,
                "R": 0.3333333, "C0": 5.11e-13, "m": 3.24, "gamma": 0.420,
                "K_th": 1.0, "K_IC": 130.0,
            },
            "variable_amplitude": {
                "enabled": False,
                "alphas": [],
            },
        }

    # =========================================================================
    # PROBLEMA 4: Seção 5.2.2 — Placa AISI 4340 com Duas Trincas de Borda
    # =========================================================================
    elif problem == 4:
        return {
            "problem": 4,
            "title": "Double Edge Cracked Plate AISI 4340 (Sec. 5.2.2)",
            "analysis_id": "s522_double_edge",
            "file_name": "s522_double_edge",
            "gmsh_params": {},
            "E": 207000.0,
            "nu": 0.30,
            "plane_state": 1,
            "npow": 3,
            "lamb": lamb_default,
            "bound_cond": {
                "bottom": {
                    "bc_dir": 0,
                    "type_x": "ux", "value_x": 0.0,
                    "type_y": "uy", "value_y": 0.0,
                },
                "free_sides": {
                    "bc_dir": 0,
                    "type_x": "tx", "value_x": 0.0,
                    "type_y": "ty", "value_y": 0.0,
                },
                "top": {
                    "bc_dir": 0,
                    "type_x": "tx", "value_x": 0.0,
                    "type_y": "ty", "value_y": 526.0,
                },
            },
            "cracks": [
                {"physical_group": "crack_1", "tips_num": 1},
                {"physical_group": "crack_2", "tips_num": 1},
            ],
            "internal_points_rects": [],
            "sif": {
                "k_key": 2, "method_key": 3,
                "nbe": 20, "nbn": 3, "lb1": lamb_default, "lb2": lamb_default,
                "np": 20, "initial_node": 2,
            },
            "fatigue": {
                "k_eff_key": 2,
                "da0": 0.0001, "r_da": 1.1, "da_max": 1.0, "steps": 44,
                "R": 0.3333333, "C0": 5.11e-13, "m": 3.24, "gamma": 0.420,
                "K_th": 1.0, "K_IC": 130.0,
            },
            "variable_amplitude": {
                "enabled": False,
                "alphas": [],
            },
        }

    # =========================================================================
    # PROBLEMA 5: Seção 5.1 — Placa Retangular 3x1 com Trinca Central (Isida)
    # =========================================================================
    elif problem == 5:
        if s51_variant == "da0_0125_r105":
            da0, r_da = 0.0125, 1.05
            aid = "s51_da0_0125_r105"
        else:
            da0, r_da = 0.01, 1.1
            aid = "s51_da0_01_r11"
        return {
            "problem": 5,
            "title": f"Rectangular Plate 3x1 Central Crack Isida ({s51_variant})",
            "analysis_id": aid,
            "file_name": "s51_internal_crack",
            "gmsh_params": {},
            "E": 30000.0,
            "nu": 0.0,
            "plane_state": 1,
            "npow": 3,
            "lamb": lamb_default,
            "bound_cond": {
                "free_horiz": {
                    "bc_dir": 0,
                    "type_x": "tx", "value_x": 0.0,
                    "type_y": "ty", "value_y": 0.0,
                },
                "right_load": {
                    "bc_dir": 0,
                    "type_x": "tx", "value_x": 10.0,
                    "type_y": "ty", "value_y": 0.0,
                },
                "left_fixed": {
                    "bc_dir": 0,
                    "type_x": "ux", "value_x": 0.0,
                    "type_y": "uy", "value_y": 0.0,
                },
            },
            "cracks": [
                {"physical_group": "crack_1", "tips_num": 2},
            ],
            "internal_points_rects": [],
            "sif": {
                "k_key": 2, "method_key": 3,
                "nbe": 20, "nbn": 3, "lb1": lamb_default, "lb2": lamb_default,
                "np": 20, "initial_node": 2,
            },
            "fatigue": {
                "k_eff_key": 2,
                "da0": da0, "r_da": r_da, "da_max": 1.0, "steps": 24,
                "R": 0.0, "C0": 1.5e-08, "m": 2.25, "gamma": 1.0,
                "K_th": 0.1, "K_IC": 1000.0,
            },
            "variable_amplitude": {
                "enabled": False,
                "alphas": [],
            },
        }

    # =========================================================================
    # PROBLEMA 6: Seção 5.3.1 — Corpo de Prova CMM1/2/3 (alpha=90°, theta=-3°)
    # =========================================================================
    elif problem == 6:
        fp = 4500.0 / (6.35 * 5.0 * (math.pi / 2.0))  # 90.229694923 MPa
        return {
            "problem": 6,
            "title": "Compact Mixed Mode Specimen CMM1-3 (alpha=90 deg, theta=-3 deg)",
            "analysis_id": "s531_cmm_90",
            "file_name": "chambers_cmm90",
            "gmsh_params": {},
            "E": 200000.0,
            "nu": 0.25,
            "plane_state": 1,
            "npow": 9,
            "lamb": lamb_default,
            "bound_cond": {
                "free_boundary": {
                    "bc_dir": 0,
                    "type_x": "tx", "value_x": 0.0,
                    "type_y": "ty", "value_y": 0.0,
                },
                "right_symmetry": {
                    "bc_dir": 0,
                    "type_x": "ux", "value_x": 0.0,
                    "type_y": "ty", "value_y": 0.0,
                },
                "right_hole_supp": {
                    "bc_dir": 0,
                    "type_x": "tx", "value_x": 0.0,
                    "type_y": "uy", "value_y": 0.0,
                },
                "left_hole_q3": {
                    "bc_dir": 0,
                    "type_x": "tx", "value_x": 0.0,
                    "type_y": "ty", "value_y": (0.0, fp, 0.0, 0.0),
                },
                "left_hole_q4": {
                    "bc_dir": 0,
                    "type_x": "tx", "value_x": 0.0,
                    "type_y": "ty", "value_y": (fp, -fp, 0.0, 0.0),
                },
            },
            "cracks": [
                {"physical_group": "crack_1", "tips_num": 1},
            ],
            "internal_points_rects": [
                {"x1": 2.0, "y1": 2.0,   "x2": 37.105, "y2": 18.74, "nx": 8, "ny": 4},
                {"x1": 2.0, "y1": 32.94, "x2": 37.105, "y2": 56.74, "nx": 8, "ny": 4},
            ],
            "sif": {
                "k_key": 2, "method_key": 3,
                "nbe": 20, "nbn": 3, "lb1": lamb_default, "lb2": lamb_default,
                "np": 20, "initial_node": 2,
            },
            "fatigue": {
                "k_eff_key": 2,
                "da0": 0.1, "r_da": 1.1, "da_max": 1.0, "steps": 38,
                "R": 1.0 / 9.0, "C0": 2.3459e-08, "m": 1.33, "gamma": 1.0,
                "K_th": 1.0, "K_IC": 1500.0,
            },
            "variable_amplitude": {
                "enabled": False,
                "alphas": [],
            },
        }

    # =========================================================================
    # PROBLEMA 7: Seção 3.1.3.1 — Placa Quadrada 10x10 com Trinca de Borda
    # =========================================================================
    elif problem == 7:
        tag = int(round(aw * 10))
        return {
            "problem": 7,
            "title": f"Square Plate 10x10 with Edge Crack (a/w = {aw:.1f})",
            "analysis_id": f"s3131_aw{tag}",
            "file_name": "s3131_edge_crack",
            "gmsh_params": {"aw": float(aw)},
            "E": 100.0,
            "nu": 0.25,
            "plane_state": 1,
            "npow": 3,
            "lamb": lamb_default,
            "bound_cond": {
                "bottom": {
                    "bc_dir": 0,
                    "type_x": "tx", "value_x": 0.0,
                    "type_y": "ty", "value_y": -1.0,
                },
                "free_sides": {
                    "bc_dir": 0,
                    "type_x": "tx", "value_x": 0.0,
                    "type_y": "ty", "value_y": 0.0,
                },
                "right_support": {
                    "bc_dir": 0,
                    "type_x": "ux", "value_x": 0.0,
                    "type_y": "uy", "value_y": 0.0,
                },
                "top": {
                    "bc_dir": 0,
                    "type_x": "tx", "value_x": 0.0,
                    "type_y": "ty", "value_y": 1.0,
                },
            },
            "cracks": [
                {"physical_group": "crack_1", "tips_num": 1},
            ],
            "internal_points_rects": [],
            "sif": {
                "k_key": 2, "method_key": 3,
                "nbe": 20, "nbn": 3, "lb1": lamb_default, "lb2": lamb_default,
                "np": 20, "initial_node": 2,
            },
            "fatigue": {
                "k_eff_key": 2,
                "da0": 0.1, "r_da": 1.1, "da_max": 1.0, "steps": 0,
                "R": 0.0, "C0": 1.0e-08, "m": 2.0, "gamma": 1.0,
                "K_th": 0.1, "K_IC": 1000.0,
            },
            "variable_amplitude": {
                "enabled": False,
                "alphas": [],
            },
        }

    # =========================================================================
    # PROBLEMA 8: Seção 3.1.3.2 — Placa 10x20 com Trinca Inclinada a 45°
    # =========================================================================
    elif problem == 8:
        tag = int(round(aw * 10))
        return {
            "problem": 8,
            "title": f"Rectangular Plate 10x20 with 45-deg Central Crack (2a/w = {aw:.1f})",
            "analysis_id": f"s3132_aw{tag}",
            "file_name": "s3132_inclined_crack",
            "gmsh_params": {"aw": float(aw)},
            "E": 100.0,
            "nu": 0.25,
            "plane_state": 1,
            "npow": 3,
            "lamb": lamb_default,
            "bound_cond": {
                "bottom": {
                    "bc_dir": 0,
                    "type_x": "ux", "value_x": 0.0,
                    "type_y": "uy", "value_y": 0.0,
                },
                "free_sides": {
                    "bc_dir": 0,
                    "type_x": "tx", "value_x": 0.0,
                    "type_y": "ty", "value_y": 0.0,
                },
                "top": {
                    "bc_dir": 0,
                    "type_x": "tx", "value_x": 0.0,
                    "type_y": "ty", "value_y": 1.0,
                },
            },
            "cracks": [
                {"physical_group": "crack_1", "tips_num": 2},
            ],
            "internal_points_rects": [],
            "sif": {
                "k_key": 2, "method_key": 3,
                "nbe": 20, "nbn": 3, "lb1": lamb_default, "lb2": lamb_default,
                "np": 20, "initial_node": 2,
            },
            "fatigue": {
                "k_eff_key": 2,
                "da0": 0.1, "r_da": 1.1, "da_max": 1.0, "steps": 0,
                "R": 0.0, "C0": 1.0e-08, "m": 2.0, "gamma": 1.0,
                "K_th": 0.1, "K_IC": 1000.0,
            },
            "variable_amplitude": {
                "enabled": False,
                "alphas": [],
            },
        }

    raise ValueError(f"Problema inválido: {problem}")
