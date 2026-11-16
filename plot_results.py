#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
plot_results.py - Visualizador e Pós-processador Gráfico para Fracture2D (Dual-BEM)
================================================================================
Este script realiza o pós-processamento completo dos resultados gerados pelo programa
Fortran Fracture2D (MECD - Método dos Elementos de Contorno Dual).

Lê os seguintes arquivos gerados pela simulação:
  - <id>_results_step_NNN.txt      : Resultados nos nós e pontos internos por passo
  - <id>_results_crack_NNN.txt     : Fatores de intensidade de tensão (K_I, K_II, K_eq)
  - <id>_geometry_crack_NNN.txt    : Coordenadas da ponta da trinca por passo
  - <id>_fatigue_results_crack_NNN.txt : Comprimento da trinca e ciclos N (passo e total)
  - <id>_cycles.txt                : Ciclos extraídos via Rainflow (amplitude variável)
  - <id>-STEP-NNN.xml              : Modelo de malha e solução em XML (formato INSANE)

Gera gráficos de alta qualidade (300 DPI) correspondentes às figuras da dissertação:
  - Figura 5.26: Trajetória de propagação da trinca (x vs y e coordenadas relativas)
  - Figura 5.27: Configuração final da malha com a trinca propagada
  - Figura 5.28: Configurações deformadas (ampliadas com fator de escala configurável)
  - Figura 5.29: Curva N (ciclos) x avanço da trinca (Delta a e a)
  - Figura 5.30: Taxa de propagação da/dN x Delta K em escala log-log (Lei de Paris)
  - Figuras 5.3/5.9/5.13/5.19: Variação dos Fatores de Intensidade de Tensão (K_I, K_II, K_eq) vs a
  - Figuras 4.4/D.1: Espectro de ciclos Rainflow (se arquivo _cycles.txt existir)
  - Painel Executivo Resumo: Dashboard consolidado com os principais resultados

Uso:
    python plot_results.py <diretorio_ou_prefixo> [opções]

Exemplos:
    python plot_results.py scratch/val_chambers
    python plot_results.py scratch/val_dowling --scale 50
    python plot_results.py scratch/val_chambers --step 38 --scale 40 --dpi 300
    python plot_results.py scratch/val_chambers/Chambers_CMM_4_geometry_crack_001.txt
================================================================================
"""

import os
import sys
import glob
import re
import argparse
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, LogFormatterSciNotation


# ==============================================================================
# CONFIGURAÇÃO DE ESTILO DOS GRÁFICOS (Dissertação / Artigo)
# ==============================================================================
def setup_plot_style():
    """Configura o estilo padrão do Matplotlib para qualidade de publicação."""
    plt.rcParams.update({
        'font.size': 11,
        'axes.labelsize': 12,
        'axes.titlesize': 13,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'figure.titlesize': 14,
        'lines.linewidth': 1.6,
        'lines.markersize': 5,
        'grid.alpha': 0.4,
        'grid.linestyle': '--',
        'axes.grid': True,
        'savefig.bbox': 'tight',
        'savefig.dpi': 300
    })


# ==============================================================================
# FUNÇÕES DE LOCALIZAÇÃO DE ARQUIVOS
# ==============================================================================
def find_simulation_files(target_path):
    """
    Localiza os arquivos de simulação a partir de um diretório, prefixo ou arquivo individual.
    Retorna um dicionário com os caminhos dos arquivos encontrados.
    """
    files_info = {
        'dir': '',
        'prefix': '',
        'geometry_crack': [],
        'results_crack': [],
        'fatigue_crack': [],
        'results_step': [],
        'xml_steps': [],
        'cycles': []
    }

    if os.path.isdir(target_path):
        target_dir = os.path.abspath(target_path)
        files_info['dir'] = target_dir
        # Encontra arquivos de geometria ou xml para inferir o prefixo
        geom_files = glob.glob(os.path.join(target_dir, "*_geometry_crack_*.txt"))
        xml_files = glob.glob(os.path.join(target_dir, "*-STEP-*.xml"))
        res_step_files = glob.glob(os.path.join(target_dir, "*_results_step_*.txt"))

        prefix = ""
        if geom_files:
            m = re.match(r"^(.*)_geometry_crack_\d+\.txt$", os.path.basename(geom_files[0]))
            if m:
                prefix = m.group(1)
        elif xml_files:
            m = re.match(r"^(.*)-STEP-\d+\.xml$", os.path.basename(xml_files[0]))
            if m:
                prefix = m.group(1)
        elif res_step_files:
            m = re.match(r"^(.*)_results_step_\d+\.txt$", os.path.basename(res_step_files[0]))
            if m:
                prefix = m.group(1)

        files_info['prefix'] = prefix
    else:
        # Caminho para arquivo ou prefixo de caminho
        target_path = os.path.abspath(target_path)
        if os.path.isfile(target_path):
            target_dir = os.path.dirname(target_path)
            fname = os.path.basename(target_path)
            files_info['dir'] = target_dir
            # Tentar extrair prefixo
            for pattern in [r"^(.*)_geometry_crack_\d+\.txt$",
                            r"^(.*)_results_crack_\d+\.txt$",
                            r"^(.*)_fatigue_results_crack_\d+\.txt$",
                            r"^(.*)_results_step_\d+\.txt$",
                            r"^(.*)-STEP-\d+\.xml$",
                            r"^(.*)_cycles\.txt$",
                            r"^(.*)_input_data\.txt$"]:
                m = re.match(pattern, fname)
                if m:
                    files_info['prefix'] = m.group(1)
                    break
            if not files_info['prefix']:
                files_info['prefix'] = os.path.splitext(fname)[0]
        else:
            # Prefixo fornecido diretamente
            target_dir = os.path.dirname(target_path)
            files_info['dir'] = target_dir if target_dir else os.getcwd()
            files_info['prefix'] = os.path.basename(target_path)

    # Coleta arquivos do prefixo (ou todos se prefixo vazio)
    prefix_filter = f"{files_info['prefix']}*" if files_info['prefix'] else "*"
    search_dir = files_info['dir']

    all_txt = glob.glob(os.path.join(search_dir, f"{prefix_filter}.txt"))
    all_xml = glob.glob(os.path.join(search_dir, f"{prefix_filter}.xml"))

    for f in all_txt:
        bn = os.path.basename(f)
        if "_geometry_crack_" in bn:
            files_info['geometry_crack'].append(f)
        elif "_fatigue_results_crack_" in bn:
            files_info['fatigue_crack'].append(f)
        elif "_results_crack_" in bn and "_fatigue_" not in bn:
            files_info['results_crack'].append(f)
        elif "_results_step_" in bn:
            files_info['results_step'].append(f)
        elif "cycles" in bn:
            files_info['cycles'].append(f)

    for f in all_xml:
        bn = os.path.basename(f)
        if "-STEP-" in bn:
            files_info['xml_steps'].append(f)

    # Ordenar listas
    files_info['geometry_crack'].sort()
    files_info['results_crack'].sort()
    files_info['fatigue_crack'].sort()
    files_info['results_step'].sort()
    files_info['xml_steps'].sort()
    files_info['cycles'].sort()

    return files_info


# ==============================================================================
# PARSERS DE DADOS
# ==============================================================================
def parse_geometry_crack(file_path):
    """
    Lê o arquivo <id>_geometry_crack_NNN.txt.
    Colunas: Step, Initial_x, Initial_y, Final_x, Final_y, Length.
    """
    steps = []
    x_init, y_init = [], []
    x_final, y_final = [], []
    length = []

    with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line_str = line.strip()
            if not line_str or line_str.startswith('Step'):
                continue
            parts = line_str.split()
            if len(parts) >= 6:
                try:
                    steps.append(int(parts[0]))
                    x_init.append(float(parts[1]))
                    y_init.append(float(parts[2]))
                    x_final.append(float(parts[3]))
                    y_final.append(float(parts[4]))
                    length.append(float(parts[5]))
                except ValueError:
                    continue

    return {
        'step': np.array(steps),
        'x_init': np.array(x_init),
        'y_init': np.array(y_init),
        'x_final': np.array(x_final),
        'y_final': np.array(y_final),
        'length': np.array(length),
        'file': file_path
    }


def parse_results_crack(file_path):
    """
    Lê o arquivo <id>_results_crack_NNN.txt.
    Colunas: Step, Length, K_I_1, K_II_1, Keq_1 (e opcionalmente K_I_2, K_II_2, Keq_2).
    """
    steps, length = [], []
    ki_1, kii_1, keq_1 = [], [], []
    ki_2, kii_2, keq_2 = [], [], []
    has_tip_2 = False

    with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line_str = line.strip()
            if not line_str or line_str.startswith('Step'):
                continue
            parts = line_str.split()
            if len(parts) >= 5:
                try:
                    steps.append(int(parts[0]))
                    length.append(float(parts[1]))
                    ki_1.append(float(parts[2]))
                    kii_1.append(float(parts[3]))
                    keq_1.append(float(parts[4]))
                    if len(parts) >= 8:
                        has_tip_2 = True
                        ki_2.append(float(parts[5]))
                        kii_2.append(float(parts[6]))
                        keq_2.append(float(parts[7]))
                except ValueError:
                    continue

    data = {
        'step': np.array(steps),
        'length': np.array(length),
        'ki_1': np.array(ki_1),
        'kii_1': np.array(kii_1),
        'keq_1': np.array(keq_1),
        'has_tip_2': has_tip_2,
        'file': file_path
    }
    if has_tip_2:
        data['ki_2'] = np.array(ki_2)
        data['kii_2'] = np.array(kii_2)
        data['keq_2'] = np.array(keq_2)

    return data


def parse_fatigue_results(file_path):
    """
    Lê o arquivo <id>_fatigue_results_crack_NNN.txt.
    Colunas: Length, N_1_step, N_1_total (e opcionalmente N_2_step, N_2_total).
    """
    length = []
    n1_step, n1_total = [], []
    n2_step, n2_total = [], []
    has_tip_2 = False

    with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line_str = line.strip()
            if not line_str or 'Length' in line_str:
                continue
            parts = line_str.split()
            if len(parts) >= 3:
                try:
                    length.append(float(parts[0]))
                    n1_step.append(float(parts[1]))
                    n1_total.append(float(parts[2]))
                    if len(parts) >= 5:
                        has_tip_2 = True
                        n2_step.append(float(parts[3]))
                        n2_total.append(float(parts[4]))
                except ValueError:
                    continue

    data = {
        'length': np.array(length),
        'n1_step': np.array(n1_step),
        'n1_total': np.array(n1_total),
        'has_tip_2': has_tip_2,
        'file': file_path
    }
    if has_tip_2:
        data['n2_step'] = np.array(n2_step)
        data['n2_total'] = np.array(n2_total)

    return data


def parse_cycles(file_path):
    """
    Lê o arquivo de contagem Rainflow (<id>_cycles.txt ou <id>_variable_ampl_cycles.txt).
    Colunas: Cycle, alpha_max, R.
    """
    cycles, alpha_max, r_ratio = [], [], []

    with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line_str = line.strip()
            if not line_str or 'Cycle' in line_str:
                continue
            parts = line_str.split()
            if len(parts) >= 3:
                try:
                    cycles.append(int(parts[0]))
                    alpha_max.append(float(parts[1]))
                    r_ratio.append(float(parts[2]))
                except ValueError:
                    continue

    return {
        'cycle': np.array(cycles),
        'alpha_max': np.array(alpha_max),
        'r': np.array(r_ratio),
        'file': file_path
    }


def parse_step_xml(file_path):
    """
    Lê o arquivo de malha e resultados em formato XML do INSANE (<id>-STEP-NNN.xml).
    Retorna nós, deslocamentos e elementos.
    """
    tree = ET.parse(file_path)
    root = tree.getroot()

    nodes = {}
    for elem in root.iter():
        if elem.tag.endswith('BemNode'):
            lbl = int(elem.attrib['label'])
            coord_elem = elem.find('{*}Coord')
            coord_txt = coord_elem.text.split() if coord_elem is not None else [0, 0, 0]
            coords = [float(c) for c in coord_txt]

            displ_elem = elem.find('.//{*}Displ')
            displ = [float(d) for d in displ_elem.text.split()] if displ_elem is not None else [0.0, 0.0]

            vm_elem = elem.find('.//{*}VonMisesStress')
            von_mises = float(vm_elem.text) if vm_elem is not None else 0.0

            nodes[lbl] = {
                'x': coords[0],
                'y': coords[1],
                'ux': displ[0],
                'uy': displ[1],
                'von_mises': von_mises
            }

    elements = []
    for elem in root.iter():
        if elem.tag.endswith('Element') and 'label' in elem.attrib:
            inc_elem = elem.find('{*}Incidence')
            if inc_elem is not None:
                inc = [int(x) for x in inc_elem.text.split()]
                elem_class = elem.attrib.get('class', '')
                elem_label = elem.attrib.get('label', '')
                elements.append({
                    'class': elem_class,
                    'label': elem_label,
                    'nodes': inc
                })

    return {
        'nodes': nodes,
        'elements': elements,
        'file': file_path
    }


# ==============================================================================
# FUNÇÕES DE PLOTAGEM (Figuras 5.26 a 5.30 e correlatas)
# ==============================================================================

def plot_crack_trajectory(geom_list, output_path, title_prefix=""):
    """
    FIGURA 5.26: Trajetória de propagação da trinca.
    Plota as coordenadas (x, y) da ponta da trinca ao longo de todos os passos.
    Gera dois subplots:
      1. Coordenadas Globais (X vs Y) no sistema da peça
      2. Coordenadas Relativas (Delta X vs Delta Y) a partir da ponta inicial
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 6))

    colors = ['#d62728', '#1f77b4', '#2ca02c', '#ff7f0e', '#9467bd']

    for idx, geom in enumerate(geom_list):
        c = colors[idx % len(colors)]
        crack_id = os.path.basename(geom['file']).split('_geometry_crack_')[-1].replace('.txt', '')

        # Coordenadas globais
        x = geom['x_final']
        y = geom['y_final']
        steps = geom['step']

        # Subplot 1: Global
        ax1.plot(x, y, '-o', color=c, lw=1.8, ms=4, label=f'Trinca {crack_id} (Passos 0-{steps[-1]})')
        ax1.plot(x[0], y[0], 's', color='green', ms=8, label=f'Início ({x[0]:.2f}, {y[0]:.2f})' if idx == 0 else "")
        ax1.plot(x[-1], y[-1], '^', color='black', ms=8, label=f'Ponta Final ({x[-1]:.2f}, {y[-1]:.2f})' if idx == 0 else "")

        # Anotar alguns passos
        n_pts = len(steps)
        annot_stride = max(1, n_pts // 8)
        for i in range(0, n_pts, annot_stride):
            ax1.annotate(f"P{steps[i]}", (x[i], y[i]), textcoords="offset points",
                         xytext=(5, 5), fontsize=8, alpha=0.85)

        # Subplot 2: Relativo à ponta inicial (Delta X, Delta Y) - Como na Fig. 5.26 da dissertação
        dx = x - x[0]
        dy = y - y[0]
        ax2.plot(dx, dy, '-o', color=c, lw=1.8, ms=4, label=f'Trinca {crack_id}')
        ax2.plot(0, 0, 's', color='green', ms=8, label='Origem (0, 0)' if idx == 0 else "")
        ax2.plot(dx[-1], dy[-1], '^', color='black', ms=8, label=f'Final ({dx[-1]:.2f}, {dy[-1]:.2f})' if idx == 0 else "")

        for i in range(0, n_pts, annot_stride):
            ax2.annotate(f"P{steps[i]}", (dx[i], dy[i]), textcoords="offset points",
                         xytext=(5, 5), fontsize=8, alpha=0.85)

    ax1.set_title(f"Trajetória Global da Trinca\n{title_prefix}")
    ax1.set_xlabel("X (mm / unid. comp.)")
    ax1.set_ylabel("Y (mm / unid. comp.)")
    ax1.set_aspect('equal', adjustable='datalim')
    ax1.legend(loc='best')
    ax1.grid(True)

    ax2.set_title(f"Trajetória Relativa da Ponta da Trinca (Fig. 5.26)\n{title_prefix}")
    ax2.set_xlabel(r"$\Delta X = X - X_0$ (mm / unid. comp.)")
    ax2.set_ylabel(r"$\Delta Y = Y - Y_0$ (mm / unid. comp.)")
    ax2.set_aspect('equal', adjustable='datalim')
    ax2.legend(loc='best')
    ax2.grid(True)

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"  [OK] Figura 5.26 salva: {output_path}")


def plot_final_mesh(xml_data, output_path, title_prefix="", step_num=None):
    """
    FIGURA 5.27: Configuração final da malha de contorno com a trinca propagada.
    Lê o modelo XML e desenha o contorno com elementos lineares e quadráticos.
    """
    fig, ax = plt.subplots(figsize=(9, 8))
    nodes = xml_data['nodes']
    elements = xml_data['elements']

    for elem in elements:
        inc = elem['nodes']
        if len(inc) == 2:
            n1, n2 = inc
            if n1 in nodes and n2 in nodes:
                ax.plot([nodes[n1]['x'], nodes[n2]['x']],
                        [nodes[n1]['y'], nodes[n2]['y']],
                        color='#2b2b2b', lw=1.0)
        elif len(inc) == 3:
            n1, nmid, n2 = inc
            if n1 in nodes and nmid in nodes and n2 in nodes:
                xi = np.linspace(-1, 1, 15)
                N1 = 0.5 * xi * (xi - 1)
                Nmid = 1 - xi**2
                N2 = 0.5 * xi * (xi + 1)
                x = N1 * nodes[n1]['x'] + Nmid * nodes[nmid]['x'] + N2 * nodes[n2]['x']
                y = N1 * nodes[n1]['y'] + Nmid * nodes[nmid]['y'] + N2 * nodes[n2]['y']
                ax.plot(x, y, color='#1f77b4', lw=1.2)

    step_txt = f" - Passo {step_num}" if step_num is not None else ""
    ax.set_title(f"Configuração Final da Malha de Contorno (Fig. 5.27)\n{title_prefix}{step_txt}")
    ax.set_xlabel("X (mm / unid. comp.)")
    ax.set_ylabel("Y (mm / unid. comp.)")
    ax.set_aspect('equal', adjustable='datalim')
    ax.grid(True, linestyle=':', alpha=0.5)

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"  [OK] Figura 5.27 salva: {output_path}")


def plot_deformed_mesh(xml_data, scale_factor, output_path, title_prefix="", step_num=None, crack_center=None):
    """
    FIGURA 5.28: Configuração deformada da malha (comparação indeformada x deformada ampliada).
    Linha preta: configuração original (indeformada).
    Linha vermelha: configuração deformada (u_x, u_y amplificados por scale_factor).
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 8))
    nodes = xml_data['nodes']
    elements = xml_data['elements']

    all_x, all_y = [], []
    all_x_def, all_y_def = [], []

    for elem in elements:
        inc = elem['nodes']
        if len(inc) == 2:
            n1, n2 = inc
            if n1 in nodes and n2 in nodes:
                x_orig = [nodes[n1]['x'], nodes[n2]['x']]
                y_orig = [nodes[n1]['y'], nodes[n2]['y']]
                x_def = [nodes[n1]['x'] + scale_factor * nodes[n1]['ux'],
                         nodes[n2]['x'] + scale_factor * nodes[n2]['ux']]
                y_def = [nodes[n1]['y'] + scale_factor * nodes[n1]['uy'],
                         nodes[n2]['y'] + scale_factor * nodes[n2]['uy']]
                for ax in (ax1, ax2):
                    ax.plot(x_orig, y_orig, color='black', lw=0.9, alpha=0.85)
                    ax.plot(x_def, y_def, color='crimson', lw=1.0, alpha=0.9)
                all_x.extend(x_orig)
                all_y.extend(y_orig)
                all_x_def.extend(x_def)
                all_y_def.extend(y_def)
        elif len(inc) == 3:
            n1, nmid, n2 = inc
            if n1 in nodes and nmid in nodes and n2 in nodes:
                xi = np.linspace(-1, 1, 20)
                N1 = 0.5 * xi * (xi - 1)
                Nmid = 1 - xi**2
                N2 = 0.5 * xi * (xi + 1)
                x0 = N1 * nodes[n1]['x'] + Nmid * nodes[nmid]['x'] + N2 * nodes[n2]['x']
                y0 = N1 * nodes[n1]['y'] + Nmid * nodes[nmid]['y'] + N2 * nodes[n2]['y']
                ux = N1 * nodes[n1]['ux'] + Nmid * nodes[nmid]['ux'] + N2 * nodes[n2]['ux']
                uy = N1 * nodes[n1]['uy'] + Nmid * nodes[nmid]['uy'] + N2 * nodes[n2]['uy']
                for ax in (ax1, ax2):
                    ax.plot(x0, y0, color='black', lw=0.9, alpha=0.85)
                    ax.plot(x0 + scale_factor * ux, y0 + scale_factor * uy, color='crimson', lw=1.0, alpha=0.9)
                all_x.extend(x0)
                all_y.extend(y0)
                all_x_def.extend(x_orig if len(inc) == 2 else [x0[0], x0[-1]])
                all_y_def.extend(y_orig if len(inc) == 2 else [y0[0], y0[-1]])

    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='black', lw=1.2, label='Indeformada'),
        Line2D([0], [0], color='crimson', lw=1.2, label=f'Deformada ({scale_factor:g}x)')
    ]

    step_txt = f" - Passo {step_num}" if step_num is not None else ""
    ax1.set_title(f"Configuração Deformada Global (Fig. 5.28)\n{title_prefix}{step_txt}")
    ax1.set_xlabel("X (mm / unid. comp.)")
    ax1.set_ylabel("Y (mm / unid. comp.)")
    ax1.set_aspect('equal', adjustable='datalim')
    ax1.legend(handles=legend_elements, loc='best')
    ax1.grid(True)

    # Subplot 2: Zoom na região da trinca
    if all_x:
        x_min, x_max = min(all_x), max(all_x)
        y_min, y_max = min(all_y), max(all_y)
        span_total = max(x_max - x_min, y_max - y_min)

        if crack_center is not None:
            cx, cy = crack_center
            zoom_radius = max(span_total * 0.22, 1.5 * abs(x_max - x_min) * 0.15)
        else:
            cx = 0.5 * (x_min + x_max)
            cy = 0.5 * (y_min + y_max)
            zoom_radius = span_total * 0.25

        ax2.set_xlim(cx - zoom_radius, cx + zoom_radius)
        ax2.set_ylim(cy - zoom_radius, cy + zoom_radius)
        ax2.set_aspect('equal', adjustable='box')

    ax2.set_title(f"Detalhe / Zoom na Trinca ({scale_factor:g}x)\n{title_prefix}")
    ax2.set_xlabel("X (mm / unid. comp.)")
    ax2.set_ylabel("Y (mm / unid. comp.)")
    ax2.legend(handles=legend_elements, loc='best')
    ax2.grid(True)

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"  [OK] Figura 5.28 salva: {output_path}")


def plot_fatigue_curve(fatigue_list, output_path, title_prefix=""):
    """
    FIGURA 5.29: Número de ciclos x Avanço da trinca (Delta a e a vs N).
    Plota N no eixo X e Delta a (e a) no eixo Y.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    colors = ['#1f77b4', '#d62728', '#2ca02c']

    for idx, fat in enumerate(fatigue_list):
        c = colors[idx % len(colors)]
        crack_id = os.path.basename(fat['file']).split('_fatigue_results_crack_')[-1].replace('.txt', '')

        L = fat['length']
        N1 = fat['n1_total']
        a0 = L[0]
        delta_a = L - a0

        # Subplot 1: Delta a vs N
        ax1.plot(N1, delta_a, '-o', color=c, lw=1.8, ms=4, label=f'Trinca {crack_id} - Ponta 1')
        if fat['has_tip_2']:
            N2 = fat['n2_total']
            ax1.plot(N2, delta_a, '--s', color='#ff7f0e', lw=1.6, ms=4, label=f'Trinca {crack_id} - Ponta 2')

        # Subplot 2: a (comprimento total) vs N
        ax2.plot(N1, L, '-o', color=c, lw=1.8, ms=4, label=f'Trinca {crack_id} - Ponta 1')
        if fat['has_tip_2']:
            ax2.plot(N2, L, '--s', color='#ff7f0e', lw=1.6, ms=4, label=f'Trinca {crack_id} - Ponta 2')

    # Notação científica no eixo N
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-2, 4))

    ax1.set_title(f"Avanço da Trinca vs Ciclos (Fig. 5.29)\n{title_prefix}")
    ax1.set_xlabel(r"Número de Ciclos $N$")
    ax1.set_ylabel(r"Avanço da Trinca $\Delta a = a - a_0$ (mm / unid. comp.)")
    ax1.xaxis.set_major_formatter(formatter)
    ax1.legend(loc='best')
    ax1.grid(True)

    ax2.set_title(f"Comprimento da Trinca vs Ciclos (Fig. 5.4 / 5.10)\n{title_prefix}")
    ax2.set_xlabel(r"Número de Ciclos $N$")
    ax2.set_ylabel(r"Comprimento da Trinca $a$ (mm / unid. comp.)")
    ax2.xaxis.set_major_formatter(formatter)
    ax2.legend(loc='best')
    ax2.grid(True)

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"  [OK] Figura 5.29 salva: {output_path}")


def plot_dadn_vs_dk(fatigue_list, results_list, output_path, title_prefix=""):
    """
    FIGURA 5.30: Taxa de propagação da/dN vs Delta K em escala log-log.
    Combina os dados de fadiga (comprimento e ciclos) com os fatores de intensidade
    de tensão (K_eq ou K_I) para traçar a curva de Paris e ajuste de potência.
    """
    fig, ax = plt.subplots(figsize=(8, 7))

    for fat, res in zip(fatigue_list, results_list):
        crack_id = os.path.basename(fat['file']).split('_fatigue_results_crack_')[-1].replace('.txt', '')

        L_fat = fat['length']
        N_fat = fat['n1_total']

        L_res = res['length']
        keq_res = res['keq_1']

        dadn_pts = []
        dk_pts = []

        for i in range(1, len(L_fat)):
            da = L_fat[i] - L_fat[i - 1]
            dN = N_fat[i] - N_fat[i - 1]
            if dN > 0 and da > 0:
                rate = da / dN
                # Pega o Keq correspondente ao comprimento
                idx = np.argmin(np.abs(L_res - L_fat[i]))
                sif = keq_res[idx]
                if sif > 0:
                    dadn_pts.append(rate)
                    dk_pts.append(sif)

        dadn_pts = np.array(dadn_pts)
        dk_pts = np.array(dk_pts)

        if len(dadn_pts) > 2:
            ax.loglog(dk_pts, dadn_pts, 'x', color='crimson', ms=7, mew=1.8,
                      label=f'Numérico - Trinca {crack_id} ({len(dadn_pts)} pts)')

            # Ajuste de Lei de Paris: da/dN = C * (Delta K)^m  =>  log10(da/dN) = log10(C) + m * log10(Delta K)
            log_x = np.log10(dk_pts)
            log_y = np.log10(dadn_pts)
            slope, intercept = np.polyfit(log_x, log_y, 1)
            C_fit = 10**intercept
            m_fit = slope

            fit_x = np.logspace(np.log10(dk_pts.min() * 0.9), np.log10(dk_pts.max() * 1.1), 50)
            fit_y = C_fit * (fit_x ** m_fit)
            ax.loglog(fit_x, fit_y, 'k-', lw=1.6,
                      label=f'Ajuste de Paris: $C={C_fit:.2e}, m={m_fit:.2f}$')

    ax.set_title(f"Taxa de Propagação $da/dN \\times \\Delta K$ (Fig. 5.30)\n{title_prefix}")
    ax.set_xlabel(r"$\Delta K$ ou $K_{eq}$ ($\mathrm{MPa\sqrt{m}}$ ou unid.)")
    ax.set_ylabel(r"$da/dN$ ($\mathrm{mm/ciclo}$ ou $\mathrm{m/ciclo}$)")
    ax.grid(True, which="both", ls="--", alpha=0.5)
    ax.legend(loc='best')

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"  [OK] Figura 5.30 salva: {output_path}")


def plot_sif_history(results_list, output_path, title_prefix=""):
    """
    FIGURAS 5.3 / 5.9 / 5.13 / 5.19: Variação dos Fatores de Intensidade de Tensão (FIT).
    Plota K_I, K_II e K_eq em função do avanço/comprimento da trinca.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    for idx, res in enumerate(results_list):
        crack_id = os.path.basename(res['file']).split('_results_crack_')[-1].replace('.txt', '')
        L = res['length']
        steps = res['step']

        # Subplot 1: SIF vs Comprimento da Trinca
        ax1.plot(L, res['keq_1'], '-o', color='black', lw=1.8, ms=4, label=f'$K_{{eq}}$ (Trinca {crack_id})')
        ax1.plot(L, res['ki_1'], '--', color='#1f77b4', lw=1.5, label=f'$K_I$ (Trinca {crack_id})')
        ax1.plot(L, np.abs(res['kii_1']), '-.', color='#d62728', lw=1.5, label=f'$|K_{{II}}|$ (Trinca {crack_id})')

        # Subplot 2: SIF vs Passo de Propagação
        ax2.plot(steps, res['keq_1'], '-o', color='black', lw=1.8, ms=4, label=f'$K_{{eq}}$')
        ax2.plot(steps, res['ki_1'], '--', color='#1f77b4', lw=1.5, label=f'$K_I$')
        ax2.plot(steps, np.abs(res['kii_1']), '-.', color='#d62728', lw=1.5, label=f'$|K_{{II}}|$')

    ax1.set_title(f"Fatores de Intensidade de Tensão vs $a$\n{title_prefix}")
    ax1.set_xlabel("Comprimento da Trinca $a$ (mm / unid. comp.)")
    ax1.set_ylabel(r"$K_I, K_{II}, K_{eq}$ ($\mathrm{MPa\sqrt{m}}$ ou unid.)")
    ax1.legend(loc='best')
    ax1.grid(True)

    ax2.set_title(f"Fatores de Intensidade de Tensão vs Passo\n{title_prefix}")
    ax2.set_xlabel("Passo de Propagação")
    ax2.set_ylabel(r"$K_I, K_{II}, K_{eq}$ ($\mathrm{MPa\sqrt{m}}$ ou unid.)")
    ax2.legend(loc='best')
    ax2.grid(True)

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"  [OK] Variação dos SIFs salva: {output_path}")


def plot_rainflow_cycles(cycles_data, output_path, title_prefix=""):
    """
    FIGURAS 4.4 / D.1: Espectro e distribuição de ciclos Rainflow.
    Plota o histograma de amplitudes de tensão alpha_max e a razão de tensões R.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    alpha = cycles_data['alpha_max']
    r_val = cycles_data['r']
    cyc = cycles_data['cycle']

    # Subplot 1: Histograma de alpha_max
    ax1.hist(alpha, bins=min(25, max(5, len(alpha) // 4)), color='#1f77b4', edgecolor='black', alpha=0.75)
    ax1.set_title(f"Distribuição de Amplitudes Rainflow\n{title_prefix}")
    ax1.set_xlabel(r"Amplitude de Tensão $\alpha_{max}$ (MPa ou unid.)")
    ax1.set_ylabel("Frequência de Ciclos")
    ax1.grid(True)

    # Subplot 2: Scatter alpha_max vs R
    scatter = ax2.scatter(alpha, r_val, c=cyc, cmap='viridis', s=35, edgecolors='k', alpha=0.85)
    cbar = plt.colorbar(scatter, ax=ax2)
    cbar.set_label("Número do Ciclo")
    ax2.set_title(f"Espectro Rainflow: $R$ vs $\\alpha_{{max}}$ (Total: {len(cyc)} ciclos)\n{title_prefix}")
    ax2.set_xlabel(r"$\alpha_{max}$")
    ax2.set_ylabel(r"Razão de Tensões $R$")
    ax2.grid(True)

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"  [OK] Espectro Rainflow salvo: {output_path}")


def plot_summary_dashboard(files_info, geom_list, res_list, fat_list, xml_data, scale_factor, output_path):
    """
    Painel Executivo Resumo: Dashboard consolidado 2x2 reunindo os 4 gráficos principais:
      (a) Trajetória da Trinca (Fig 5.26)
      (b) Deformada Ampliada (Fig 5.28)
      (c) Curva N x Delta a (Fig 5.29)
      (d) Curva da/dN x Delta K (Fig 5.30)
    """
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 14))
    prefix = files_info.get('prefix', 'Simulação Dual-BEM')

    # (a) Trajetória
    if geom_list:
        g = geom_list[0]
        dx = g['x_final'] - g['x_final'][0]
        dy = g['y_final'] - g['y_final'][0]
        ax1.plot(dx, dy, '-o', color='crimson', lw=1.8, ms=4, label='Ponta da Trinca')
        ax1.plot(0, 0, 's', color='green', ms=8, label='Início')
        ax1.plot(dx[-1], dy[-1], '^', color='black', ms=8, label='Final')
        ax1.set_title(f"(a) Trajetória da Trinca (Fig. 5.26)")
        ax1.set_xlabel(r"$\Delta X$ (mm)")
        ax1.set_ylabel(r"$\Delta Y$ (mm)")
        ax1.set_aspect('equal', adjustable='datalim')
        ax1.legend(loc='best')
        ax1.grid(True)

    # (b) Deformada
    if xml_data:
        nodes = xml_data['nodes']
        for elem in xml_data['elements']:
            inc = elem['nodes']
            if len(inc) == 2:
                n1, n2 = inc
                if n1 in nodes and n2 in nodes:
                    ax2.plot([nodes[n1]['x'], nodes[n2]['x']],
                             [nodes[n1]['y'], nodes[n2]['y']], 'k-', lw=0.8, alpha=0.7)
                    ax2.plot([nodes[n1]['x'] + scale_factor * nodes[n1]['ux'], nodes[n2]['x'] + scale_factor * nodes[n2]['ux']],
                             [nodes[n1]['y'] + scale_factor * nodes[n1]['uy'], nodes[n2]['y'] + scale_factor * nodes[n2]['uy']],
                             'r-', lw=0.9, alpha=0.85)
            elif len(inc) == 3:
                n1, nmid, n2 = inc
                if n1 in nodes and nmid in nodes and n2 in nodes:
                    xi = np.linspace(-1, 1, 15)
                    N1 = 0.5 * xi * (xi - 1)
                    Nmid = 1 - xi**2
                    N2 = 0.5 * xi * (xi + 1)
                    x0 = N1 * nodes[n1]['x'] + Nmid * nodes[nmid]['x'] + N2 * nodes[n2]['x']
                    y0 = N1 * nodes[n1]['y'] + Nmid * nodes[nmid]['y'] + N2 * nodes[n2]['y']
                    ux = N1 * nodes[n1]['ux'] + Nmid * nodes[nmid]['ux'] + N2 * nodes[n2]['ux']
                    uy = N1 * nodes[n1]['uy'] + Nmid * nodes[nmid]['uy'] + N2 * nodes[n2]['uy']
                    ax2.plot(x0, y0, 'k-', lw=0.8, alpha=0.7)
                    ax2.plot(x0 + scale_factor * ux, y0 + scale_factor * uy, 'r-', lw=0.9, alpha=0.85)
        ax2.set_title(f"(b) Deformada Ampliada {scale_factor:g}x (Fig. 5.28)")
        ax2.set_xlabel("X (mm)")
        ax2.set_ylabel("Y (mm)")
        ax2.set_aspect('equal', adjustable='datalim')
        ax2.grid(True)

    # (c) Curva N x Delta a
    if fat_list:
        fat = fat_list[0]
        L = fat['length']
        N = fat['n1_total']
        ax3.plot(N, L - L[0], '-o', color='#1f77b4', lw=1.8, ms=4, label='Numérico')
        ax3.set_title("(c) Avanço da Trinca $\\Delta a \\times N$ (Fig. 5.29)")
        ax3.set_xlabel(r"Número de Ciclos $N$")
        ax3.set_ylabel(r"$\Delta a = a - a_0$ (mm)")
        formatter = ScalarFormatter(useMathText=True)
        formatter.set_powerlimits((-2, 4))
        ax3.xaxis.set_major_formatter(formatter)
        ax3.legend(loc='best')
        ax3.grid(True)

    # (d) da/dN x Delta K
    if fat_list and res_list:
        fat, res = fat_list[0], res_list[0]
        L_fat, N_fat = fat['length'], fat['n1_total']
        L_res, keq_res = res['length'], res['keq_1']

        dadn_pts, dk_pts = [], []
        for i in range(1, len(L_fat)):
            da = L_fat[i] - L_fat[i - 1]
            dN = N_fat[i] - N_fat[i - 1]
            if dN > 0 and da > 0:
                idx = np.argmin(np.abs(L_res - L_fat[i]))
                sif = keq_res[idx]
                if sif > 0:
                    dadn_pts.append(da / dN)
                    dk_pts.append(sif)

        if len(dadn_pts) > 2:
            ax4.loglog(dk_pts, dadn_pts, 'x', color='crimson', ms=6, mew=1.5, label='Numérico')
            log_x = np.log10(dk_pts)
            log_y = np.log10(dadn_pts)
            slope, intercept = np.polyfit(log_x, log_y, 1)
            C_fit, m_fit = 10**intercept, slope
            fit_x = np.logspace(np.log10(min(dk_pts) * 0.95), np.log10(max(dk_pts) * 1.05), 50)
            ax4.loglog(fit_x, C_fit * (fit_x ** m_fit), 'k-', lw=1.5,
                       label=f'Paris: $m={m_fit:.2f}$')

        ax4.set_title(r"(d) Taxa $da/dN \times \Delta K$ Log-Log (Fig. 5.30)")
        ax4.set_xlabel(r"$\Delta K$ ou $K_{eq}$")
        ax4.set_ylabel(r"$da/dN$ (mm/ciclo)")
        ax4.grid(True, which="both", ls="--", alpha=0.5)
        ax4.legend(loc='best')

    fig.suptitle(f"Painel Resumo da Simulação: {prefix}", fontsize=16, fontweight='bold', y=0.995)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"  [OK] Painel Resumo salvo: {output_path}")


# ==============================================================================
# FUNÇÃO PRINCIPAL E CLI
# ==============================================================================
def main():
    parser = argparse.ArgumentParser(
        description="Visualizador de Resultados e Curvas de Fadiga Dual-BEM (Fracture2D)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemplos de uso:
  python plot_results.py scratch/val_chambers
  python plot_results.py scratch/val_dowling --scale 50
  python plot_results.py scratch/val_chambers --step 38 --scale 40
  python plot_results.py scratch/val_chambers/Chambers_CMM_4_geometry_crack_001.txt
        """
    )
    parser.add_argument("path", help="Caminho do diretório de saída, prefixo de análise ou arquivo de resultado")
    parser.add_argument("--scale", type=float, default=40.0,
                        help="Fator de amplificação para a configuração deformada (padrão: 40.0)")
    parser.add_argument("--step", type=int, default=None,
                        help="Passo específico para os gráficos de malha e deformada (padrão: último passo)")
    parser.add_argument("--output-dir", type=str, default=None,
                        help="Diretório de saída para salvar os gráficos PNG (padrão: mesmo diretório dos dados)")
    parser.add_argument("--dpi", type=int, default=300,
                        help="Resolução das imagens PNG (padrão: 300 DPI)")
    parser.add_argument("--trajectory", action="store_true", help="Gera apenas a trajetória da trinca (Fig 5.26)")
    parser.add_argument("--mesh", action="store_true", help="Gera apenas a malha final (Fig 5.27)")
    parser.add_argument("--deformed", action="store_true", help="Gera apenas a malha deformada (Fig 5.28)")
    parser.add_argument("--fatigue", action="store_true", help="Gera apenas a curva N vs a (Fig 5.29)")
    parser.add_argument("--dadn", action="store_true", help="Gera apenas da/dN vs dK em escala log-log (Fig 5.30)")
    parser.add_argument("--sif", action="store_true", help="Gera gráfico de variação de K_I, K_II, K_eq")
    parser.add_argument("--cycles", action="store_true", help="Gera gráfico do espectro de ciclos Rainflow")
    parser.add_argument("--dashboard", action="store_true", help="Gera painel executivo consolidado")
    parser.add_argument("--all", action="store_true", help="Gera todas as figuras disponíveis (padrão)")

    args = parser.parse_args()

    # Se nenhum flag individual foi escolhido, ativa --all
    selected_specific = (args.trajectory or args.mesh or args.deformed or
                         args.fatigue or args.dadn or args.sif or args.cycles or args.dashboard)
    if not selected_specific:
        args.all = True

    setup_plot_style()
    plt.rcParams['savefig.dpi'] = args.dpi

    print("=" * 70)
    print(" Pós-processador Fracture2D (Dual-BEM)")
    print("=" * 70)

    # 1. Localiza os arquivos
    files_info = find_simulation_files(args.path)
    prefix = files_info['prefix'] if files_info['prefix'] else "fracture2d"
    out_dir = os.path.abspath(args.output_dir) if args.output_dir else files_info['dir']
    os.makedirs(out_dir, exist_ok=True)

    print(f"Diretório dos dados: {files_info['dir']}")
    print(f"Prefixo identificado: {prefix}")
    print(f"Diretório de saída  : {out_dir}")
    print(f"Arquivos localizados:")
    print(f"  - Geometria da Trinca: {len(files_info['geometry_crack'])}")
    print(f"  - Resultados Trinca  : {len(files_info['results_crack'])}")
    print(f"  - Resultados Fadiga  : {len(files_info['fatigue_crack'])}")
    print(f"  - Passos XML         : {len(files_info['xml_steps'])}")
    print(f"  - Passos Resultados  : {len(files_info['results_step'])}")
    print(f"  - Arquivo de Ciclos  : {len(files_info['cycles'])}")
    print("-" * 70)

    # 2. Carrega os dados
    geom_list = [parse_geometry_crack(f) for f in files_info['geometry_crack']]
    res_list = [parse_results_crack(f) for f in files_info['results_crack']]
    fat_list = [parse_fatigue_results(f) for f in files_info['fatigue_crack']]
    cycles_list = [parse_cycles(f) for f in files_info['cycles']]

    # Seleciona o passo XML para a malha / deformada
    xml_data = None
    chosen_step = None
    if files_info['xml_steps']:
        if args.step is not None:
            pattern = f"-STEP-{args.step:03d}.xml"
            matching = [f for f in files_info['xml_steps'] if f.endswith(pattern)]
            chosen_xml = matching[0] if matching else files_info['xml_steps'][-1]
            chosen_step = args.step
        else:
            chosen_xml = files_info['xml_steps'][-1]
            m = re.search(r"-STEP-(\d+)\.xml$", chosen_xml)
            chosen_step = int(m.group(1)) if m else len(files_info['xml_steps']) - 1
        print(f"Carregando XML da malha (Passo {chosen_step}): {os.path.basename(chosen_xml)}")
        xml_data = parse_step_xml(chosen_xml)

    # 3. Gera os Gráficos
    print("\nGerando figuras...")

    # (a) Trajetória da Trinca (Figura 5.26)
    if (args.all or args.trajectory) and geom_list:
        out_file = os.path.join(out_dir, f"{prefix}_fig526_trajetoria.png")
        plot_crack_trajectory(geom_list, out_file, title_prefix=prefix)

    # (b) Malha Final (Figura 5.27)
    if (args.all or args.mesh) and xml_data:
        out_file = os.path.join(out_dir, f"{prefix}_fig527_malha_final.png")
        plot_final_mesh(xml_data, out_file, title_prefix=prefix, step_num=chosen_step)

    # (c) Malha Deformada (Figura 5.28)
    if (args.all or args.deformed) and xml_data:
        out_file = os.path.join(out_dir, f"{prefix}_fig528_deformada_{int(args.scale)}x.png")
        crack_center = None
        if geom_list and len(geom_list[0]['x_final']) > 0:
            crack_center = (geom_list[0]['x_final'][-1], geom_list[0]['y_final'][-1])
        plot_deformed_mesh(xml_data, args.scale, out_file, title_prefix=prefix, step_num=chosen_step, crack_center=crack_center)

    # (d) Ciclos vs Avanço da Trinca (Figura 5.29)
    if (args.all or args.fatigue) and fat_list:
        out_file = os.path.join(out_dir, f"{prefix}_fig529_ciclos_vs_trinca.png")
        plot_fatigue_curve(fat_list, out_file, title_prefix=prefix)

    # (e) Taxa da/dN vs Delta K em escala log-log (Figura 5.30)
    if (args.all or args.dadn) and fat_list and res_list:
        out_file = os.path.join(out_dir, f"{prefix}_fig530_dadn_vs_dK.png")
        plot_dadn_vs_dk(fat_list, res_list, out_file, title_prefix=prefix)

    # (f) Fatores de Intensidade de Tensão K_I, K_II, K_eq vs a (Figs 5.3 / 5.9 / 5.13 / 5.19)
    if (args.all or args.sif) and res_list:
        out_file = os.path.join(out_dir, f"{prefix}_sif_evolucao.png")
        plot_sif_history(res_list, out_file, title_prefix=prefix)

    # (g) Espectro Rainflow (Figuras 4.4 / D.1)
    if (args.all or args.cycles) and cycles_list:
        out_file = os.path.join(out_dir, f"{prefix}_rainflow_ciclos.png")
        plot_rainflow_cycles(cycles_list[0], out_file, title_prefix=prefix)

    # (h) Dashboard Resumo
    if (args.all or args.dashboard) and (geom_list or xml_data or fat_list):
        out_file = os.path.join(out_dir, f"{prefix}_dashboard_resumo.png")
        plot_summary_dashboard(files_info, geom_list, res_list, fat_list, xml_data, args.scale, out_file)

    print("\n" + "=" * 70)
    print(" Processamento concluído com sucesso!")
    print(f" Todas as figuras salvas em: {out_dir}")
    print("=" * 70)


if __name__ == "__main__":
    main()
