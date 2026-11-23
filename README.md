# Fracture 2D — Propagação Bidimensional de Trincas por Fadiga (Dual-BEM com Elementos Descontínuos)

Implementação pedagógica e modular em **Fortran 90 (`gfortran` + OpenMP + Netlib LAPACK/BLAS)** acoplada a pré-processamento de malhas em **Gmsh (`.geo` $\to$ `.msh`) + Python (`meshio`)**, referente à dissertação de mestrado de **Beatriz Borges Faria Fonseca** (UFMG, 2023), disponível neste repositório em [`DissertaçãoBeatrizFonseca.pdf`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/Disserta%C3%A7%C3%A3oBeatrizFonseca.pdf):

> *Modelagem computacional da propagação bidimensional de trincas por fadiga utilizando o Método dos Elementos de Contorno Dual.*

---

## 1. Visão Geral e Fluxo de Trabalho (Arquitetura Refatorada)

Para tornar o código mais didático, legível para os alunos e livre da entrada manual complexa de segmentos, arcos e polinômios no Fortran, o projeto foi completamente reestruturado em torno de **4 etapas integradas**:

1. **Geração da Geometria e Malha no Gmsh ([`meshes/*.geo`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/meshes)):**
   - A geometria do sólido e as trincas iniciais são modeladas em scripts paramétricos `.geo` do Gmsh com grupos físicos (`Physical Curve`) nomeados para as bordas externas e para cada trinca.
   - O Gmsh gera uma malha 1D curvilínea de segunda ordem (`line3`, elementos quadráticos contínuos de 3 nós) exportada em formato `.msh`.
2. **Pré-processamento e Conversão para Elementos Descontínuos ([`input_data.py`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/input_data.py) e [`preproc.py`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/preproc.py)):**
   - Como o Gmsh gera apenas elementos contínuos, o módulo [`preproc.py`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/preproc.py) lê o arquivo `.msh` através da biblioteca **`meshio`** e converte **todos os elementos** (tanto do contorno externo quanto das trincas) em **elementos quadráticos descontínuos de 3 nós** ($\lambda_1 = \lambda_2 = 2/3$, com nós funcionais em $\xi \in \{-2/3, 0, +2/3\}$).
   - Para cada curva física de trinca (`CRACK_*`), o pré-processador duplica automaticamente os elementos com orientações opostas para montar o par de faces coincidentes do **Dual-BEM**:
     - **Face 1 (`bie = 2`):** Equação Integral de Contorno (EIC) Singular de Deslocamento;
     - **Face 2 (`bie = 3`):** Equação Integral de Contorno (EIC) Hipersingular de Tração;
     - **Contorno externo (`bie = 1`):** EIC Singular de Deslocamento.
   - Gera automaticamente a figura da malha e condições de contorno (`<id>_mesh_bc.png`) e grava o arquivo `.dat` pronto para o solver Fortran.
3. **Solver Fortran 90 Exclusivamente com Elementos Quadráticos Descontínuos ([`src/*.f90`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/src)):**
   - Como **todos os elementos são descontínuos**, cada nó funcional pertence a um único elemento ($\text{nós} = 3 \times N_{\text{elem}}$), eliminando completamente:
     - ambiguidades de vetor normal em quinas ou bocas de trinca;
     - conflitos de condições de contorno em vértices (nunca são necessárias "condições de contorno em pontos particulares");
     - ramificações de código entre elementos contínuos, semidescontínuos e descontínuos nas funções de forma e nas integrais singulares.
   - Paralelizado com **OpenMP (`-fopenmp`)** nas rotinas de montagem das matrizes $[H]$ e $[G]$, avaliação de pontos internos e contorno circular da Integral $J$.
4. **Pós-processamento Gráfico ([`plot_results.py`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/plot_results.py) e [`plot_mesh_bc.py`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/plot_mesh_bc.py)):**
   - Geração automática de todas as figuras da dissertação (trajetória da trinca, malha final, configuração deformada ampliada, curvas de fadiga $N \times a$, lei de Paris $da/dN \times \Delta K$, evolução de $K_I, K_{II}, K_{\text{eq}}$ e espectro Rainflow).

---

## 2. Formulação Numérica Implementada

- **Método dos Elementos de Contorno Dual (Dual-BEM):**
  - EIC de deslocamentos nos nós do contorno externo (`bie = 1`) e da Face 1 da trinca (`bie = 2`);
  - EIC hipersingular de tração nos nós da Face 2 da trinca (`bie = 3`), utilizando as soluções fundamentais 2D de Kelvin ([`src/kelvin_sol_mod.f90`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/src/kelvin_sol_mod.f90)).
- **Integração Numérica Especializada ([`src/integration_mod.f90`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/src/integration_mod.f90)):**
  - *Regular / quase singular:* quadratura de Gauss-Legendre de 20 pontos ([`src/quadrature_mod.f90`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/src/quadrature_mod.f90));
  - *Fracamente singular ($\ln r$ na EIC de deslocamento):* transformação de coordenadas e quadratura de Gauss-Laguerre de 20 pontos;
  - *Forte ($1/r$) e hipersingular ($1/r^2$):* **Método da Expansão Radial de Gao** (3 termos de expansão polinomial).
- **Fatores de Intensidade de Tensão ($K_I, K_{II}$) ([`src/sif_mod.f90`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/src/sif_mod.f90)):**
  - *Correlação de Deslocamentos (`K_key = 1`):* 7 variantes de extrapolação a partir da abertura relativa (COD/CSD) nos nós descontínuos do elemento da ponta;
  - *Integral $J$ desacoplada (`K_key = 2`):* decomposição simétrica/antissimétrica dos campos elásticos sobre contorno circular ao redor da ponta da trinca (Regra de Simpson `method_key = 1`, pontos de Gauss `method_key = 2`, ou elementos quadráticos descontínuos no caminho circular `method_key = 3`).
- **Propagação Incremental e Vida em Fadiga ([`src/fatigue_prop_mod.f90`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/src/fatigue_prop_mod.f90)):**
  - Direção de propagação pelo critério da **Máxima Tensão Circunferencial (MTC)** (Erdogan & Sih, 1963) com esquema preditor-corretor iterativo de **Portela et al. (1993)**;
  - Montagem incremental das matrizes $[H]$ e $[G]$ durante a propagação (reaproveitando o bloco das etapas anteriores e integrando apenas as linhas/colunas dos novos elementos da trinca);
  - Lei de **Paris-Walker** para amplitude constante e contagem de ciclos **Rainflow** para amplitude variável, integrando o número de ciclos $N$ pela **Regra 1/3 de Simpson** a cada par de incrementos.

---

## 3. Estrutura do Repositório

```text
CrackPropagation2D/
├── DissertaçãoBeatrizFonseca.pdf             # Dissertação de referência (Beatriz Fonseca, UFMG, 2023)
├── Run_crack_propagation_2d.ipynb            # Jupyter Notebook didático passo a passo para os alunos
├── Makefile                                  # Compilação do solver Fortran 90 (gfortran + OpenMP + LAPACK)
├── README.md                                 # Documentação completa e validação numérica
├── input_data.py                             # Dicionário de casos e parâmetros físicos/numéricos dos exemplos
├── preproc.py                                # Pré-processador Gmsh (.geo -> .msh) + meshio -> .dat descontínuo
├── run_all_examples.py                       # Script automatizado que roda e valida todos os casos da dissertação
├── plot_mesh_bc.py                           # Visualização da malha descontínua e condições de contorno (.png)
├── plot_results.py                           # Pós-processador gráfico (Figuras 5.26 a 5.30, SIFs, Rainflow)
├── animate_propagation.py                    # Gerador de animações em vídeo (.mp4) da propagação da trinca
├── bin/
│   └── fracture_2d.exe                       # Executável compilado estaticamente (gfortran + OpenMP + LAPACK)
├── meshes/                                   # Arquivos geométricos .geo (Gmsh) e malhas .msh geradas
│   ├── s3131_edge_crack.geo                  # Seção 3.1.3.1: Placa quadrada com trinca de borda (Modo I)
│   ├── s3132_inclined_crack.geo              # Seção 3.1.3.2: Placa 10x20 com trinca central a 45° (Modo Misto)
│   ├── s51_internal_crack.geo                # Seção 5.1: Placa 3x1 com trinca central (Isida / Andrade)
│   ├── s521_single_edge.geo                  # Seção 5.2.1: Placa AISI 4340 com 1 trinca de borda
│   ├── s522_double_edge.geo                  # Seção 5.2.2: Placa AISI 4340 com 2 trincas de borda simétricas
│   ├── dowling_11_6.geo                      # Seção 5.2.3: Trinca central sob amplitude variável (Dowling 11.6)
│   ├── chambers_cmm90.geo                    # Seção 5.3.1: Corpo de prova CMM (α = 90°, CMM1, CMM2 e CMM3)
│   └── chambers_cmm4.geo                     # Seção 5.3.2: Corpo de prova CMM (α = 110°, CMM4)
├── src/                                      # Código-fonte Fortran 90 refatorado (10 módulos didáticos)
│   ├── constants_mod.f90                     # Precisão dupla (dp), constantes matemáticas e tratamento de erros
│   ├── quadrature_mod.f90                    # Quadraturas 1D de Gauss-Legendre e Gauss-Laguerre (1 a 20 pontos)
│   ├── model_mod.f90                         # Tipos de dados, malha descontínua, funções de forma 1D e DGESV
│   ├── kelvin_sol_mod.f90                    # Soluções fundamentais 2D de Kelvin (U, T, D, S e termos de Gao)
│   ├── integration_mod.f90                   # Integrais regulares, fracamente singulares e expansão radial de Gao
│   ├── assembler_mod.f90                     # Montagem Dual-BEM ([H], [G], [A]{x}={f}) e pontos internos
│   ├── sif_mod.f90                           # Cálculo de K_I e K_II (Correlação de Deslocamentos e Integral J)
│   ├── fatigue_prop_mod.f90                  # Rainflow, critério MTC, corretor de Portela e vida em fadiga (Simpson)
│   ├── io_mod.f90                            # Leitura do arquivo .dat e escrita dos relatórios .txt e .xml
│   └── main.f90                              # Programa principal (fluxo passo a passo da análise)
├── inputs/                                   # Arquivos .dat pré-processados e figuras <id>_mesh_bc.png
├── results/                                  # Resultados completos de cada simulação (<id>_results_*, XML)
└── examples/                                 # Figuras de validação e pós-processamento geradas (300 DPI)
```

---

## 4. Como Compilar e Executar

### 4.1. Compilação do Solver Fortran (`bin/fracture_2d.exe`)

Requisitos no Windows (**MSYS2 UCRT64**):
- `gfortran` (`mingw-w64-ucrt-x86_64-gcc-fortran`)
- `lapack` / `blas` (`mingw-w64-ucrt-x86_64-lapack`)
- `mingw32-make` (`mingw-w64-ucrt-x86_64-make`)

Para compilar:

```powershell
$env:PATH = "C:\msys64\ucrt64\bin;$env:PATH"
mingw32-make clean
mingw32-make
```

O executável `bin/fracture_2d.exe` é gerado com linkagem estática (`-static -fopenmp -llapack -lblas`), podendo ser executado diretamente em qualquer terminal Windows sem depender de DLLs externas.

### 4.2. Pré-processamento de um Exemplo (`.geo` $\to$ `.msh` $\to$ `.dat`)

Requisitos Python: `gmsh`, `meshio`, `numpy`, `matplotlib`.

```powershell
# Listar todos os casos configurados em input_data.py:
python preproc.py --list

# Gerar a malha no Gmsh, converter via meshio para elementos descontínuos (.dat)
# e salvar a figura da malha e CCs (<id>_mesh_bc.png) em inputs/:
python preproc.py Chambers_CMM_4
```

### 4.3. Execução do Solver Fortran e Pós-processamento

```powershell
# 1. Executar o solver Fortran passando o arquivo .dat gerado:
.\bin\fracture_2d.exe inputs\Chambers_CMM_4.dat

# 2. Gerar todas as figuras de pós-processamento (Figuras 5.26 a 5.30, SIFs, dashboard):
python plot_results.py results\Chambers_CMM_4 --output-dir examples --scale 40

# 3. Gerar a animação em vídeo (.mp4) da propagação da trinca:
python animate_propagation.py results\Chambers_CMM_4 --scale 40 --fps 6
```

### 4.4. Execução Automatizada de Todos os Exemplos da Dissertação

Para gerar todas as malhas no Gmsh, rodar todas as 16 simulações de validação dos Capítulos 3 e 5 e produzir todas as figuras de pós-processamento de uma só vez:

```powershell
python run_all_examples.py
```

---

## 5. Validação Numérica Completa (Novo Framework com Elementos Descontínuos)

Todos os exemplos da dissertação foram reexecutados utilizando o novo fluxo (`.geo` $\to$ `meshio` $\to$ elementos quadráticos descontínuos com $\lambda_1 = \lambda_2 = 2/3$ em todo o contorno e nas trincas).

### 5.1. Seção 3.1.3.1 — Placa Quadrada com Trinca de Borda ($10 \times 10$, Modo I Puro)

Placa sob tração uniforme ($\sigma = 1$), largura $w = 10$, razão $a/w \in \{0.2, 0.3, 0.4, 0.5, 0.6\}$, malha gerada a partir de [`meshes/s3131_edge_crack.geo`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/meshes/s3131_edge_crack.geo) (`64` elementos quadráticos descontínuos, `192` nós funcionais). Comparação de $K_I / (\sigma\sqrt{\pi a})$ via Integral $J$ com a solução de referência de Civelek & Erdogan (Tabela 3.4 da dissertação):

| $a/w$ | $K_I$ (Dual-BEM) | $K_I / (\sigma\sqrt{\pi a})$ (Novo Framework) | Ref. Civelek & Erdogan (Tab. 3.4) | Erro Relativo (%) |
|---:|---:|---:|---:|---:|
| 0,2 | 3,7274 | **1,487** | 1,488 | **-0,07%** |
| 0,3 | 5,6726 | **1,848** | 1,848 | **-0,01%** |
| 0,4 | 8,2349 | **2,323** | 2,324 | **-0,04%** |
| 0,5 | 11,9134 | **3,006** | 3,010 | **-0,14%** |
| 0,6 | 17,9737 | **4,140** | 4,152 | **-0,29%** |

> **Observação:** Com elementos descontínuos em todo o contorno, o erro máximo caiu de `0,78%` (código original) para apenas **`0,29%`**, pois não há singularidade ou descontinuidade de tração nos cantos nem na boca da trinca.

### 5.2. Seção 3.1.3.2 — Placa Retangular com Trinca Central Inclinada a $45^\circ$ ($10 \times 20$, Modo Misto)

Placa retangular ($w = 10, h = 20, \sigma = 1$) com trinca interna inclinada a $45^\circ$ ($2a/w \in \{0.2, 0.3, 0.4, 0.5, 0.6\}$), malha gerada a partir de [`meshes/s3132_inclined_crack.geo`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/meshes/s3132_inclined_crack.geo) (`36` elementos quadráticos descontínuos, `108` nós funcionais). Comparação com Murakami (Tabela 3.8 da dissertação):

| $2a/w$ | $K_{I,\text{norm}}$ (Dual-BEM) | Ref. $K_I$ | Erro $K_I$ (%) | $K_{II,\text{norm}}$ (Dual-BEM) | Ref. $K_{II}$ | Erro $K_{II}$ (%) |
|---:|---:|---:|---:|---:|---:|---:|
| 0,2 | **0,520** | 0,518 | +0,46% | **0,513** | 0,507 | +1,20% |
| 0,3 | **0,543** | 0,541 | +0,41% | **0,522** | 0,516 | +1,15% |
| 0,4 | **0,575** | 0,572 | +0,54% | **0,535** | 0,529 | +1,07% |
| 0,5 | **0,616** | 0,612 | +0,67% | **0,552** | 0,546 | +1,02% |
| 0,6 | **0,667** | 0,661 | +0,88% | **0,573** | 0,567 | +1,11% |

### 5.3. Seção 5.1 — Propagação de Trinca Central sob Amplitude Constante (Isida / Andrade)

Placa retangular ($2b = 1, 2h = 3, 2a_0 = 0,1$, $\sigma_{\max} = 10$, $R = 0$, $C_0 = 1,5 \times 10^{-8}$, $m = 2,25$), malha gerada por [`meshes/s51_internal_crack.geo`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/meshes/s51_internal_crack.geo) (`88` elementos descontínuos, `264` nós iniciais), propagada por 24 incrementos:

| Esquema de avanço | $2a_f$ (obtido / diss.) | $K_{I,0}$ (Dual-BEM / Isida) | $K_{I,f}$ (Dual-BEM / Isida) | $N_f$ (Dual-BEM) | $N_f$ (Eq. 5.4) | Erro $N_f$ |
|---|---:|---:|---:|---:|---:|---:|
| $\Delta a_0 = 0,01; \; r = 1,10$ | 0,9554 / 0,955 | 4,0169 / 3,9870 | 48,9424 / 47,5443 | **229 281,6** | 229 938,9 | **-0,29%** |
| $\Delta a_0 = 0,0125; \; r = 1,05$ | 0,8959 / 0,896 | 4,0169 / 3,9870 | 30,2906 / 29,9411 | **228 756,3** | 229 311,9 | **-0,24%** |

### 5.4. Seções 5.2.1 e 5.2.2 — Placas de Aço AISI 4340 (1 e 2 Trincas de Borda)

Placas de largura $76\text{ mm}$ em aço AISI 4340 ($E = 207\text{ GPa}, \nu = 0,3, R = 1/3, K_{IC} = 130\text{ MPa}\sqrt{\text{m}}$), malhas geradas por [`meshes/s521_single_edge.geo`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/meshes/s521_single_edge.geo) e [`meshes/s522_double_edge.geo`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/meshes/s522_double_edge.geo):

| Exemplo | $a_0 \to a_f$ (mm) | $\Delta K_0$ (Dual-BEM / Ref.) | $\Delta K_f$ (Dual-BEM / Ref.) | $K_{\max,f}$ ($\text{MPa}\sqrt{\text{m}}$) | $N_f$ (Dual-BEM) | $N_f$ (Analítico) | Erro $N_f$ |
|---|---:|---:|---:|---:|---:|---:|---:|
| **Sec. 5.2.1** (1 trinca de borda) | $2,000 \to 12,232$ | 31,526 / 31,616 | 88,150 / 89,493 | 132,22 ($> K_{IC}$) | **25 819,8** | 25 070,4 | **+2,99%** |
| **Sec. 5.2.2** (2 trincas de borda) | $1,000 \to 15,281$ | 21,334 / 22,051 | 86,697 / 86,969 | 130,05 ($> K_{IC}$) | **54 324,8** | 53 457,6 | **+1,62%** |

### 5.5. Seção 5.2.3 — Amplitude Variável com Contagem Rainflow (`Dowling_11.6_variable_ampl`)

Placa com trinca central ($2a_0 = 2,00\text{ mm} \to 2a_f = 31,637\text{ mm}$ em 46 passos), malha gerada por [`meshes/dowling_11_6.geo`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/meshes/dowling_11_6.geo) (`136` elementos descontínuos iniciais):
- **Ciclos Rainflow extraídos:** `166` ciclos, com fator equivalente de Walker $\alpha_q = 1,414222$ e pico $\alpha_{\max} = 2,4$;
- **Fator de intensidade de tensão máximo:** $K_{I,\max,0} = 29,6106\text{ MPa}\sqrt{\text{m}}$ (Feddersen: $29,5126$, erro `+0,33%`) e $K_{I,\max,f} = 131,3822\text{ MPa}\sqrt{\text{m}}$ (Feddersen: $131,6936$, erro `-0,24%`), interrompendo no passo 46 por atingir $K_{IC} = 130\text{ MPa}\sqrt{\text{m}}$;
- **Vida total prevista:** **`238 810,1` ciclos** (reproduzindo os `5` algarismos significativos da validação original).

### 5.6. Seções 5.3.1 e 5.3.2 — Corpos de Prova *Compact Mixed Mode* (CMM $90^\circ$ e $110^\circ$)

Ensaios de modo misto de Chambers et al. (1991) em aço Jethete M152, malhas curvilíneas geradas por [`meshes/chambers_cmm90.geo`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/meshes/chambers_cmm90.geo) e [`meshes/chambers_cmm4.geo`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/meshes/chambers_cmm4.geo) (`172` elementos quadráticos descontínuos iniciais, `516` nós funcionais, 38 passos de propagação):

| Configuração | $(K_I / K_{II})_0$ | $K_{II}$ no Passo 1 | Ensaio ($\Delta a$) | $N$ (Novo Framework) | $N_{\text{num}}$ (Diss.) | Dif. vs Diss. | $N_{\text{exp}}$ (Chambers) |
|---|---|---|---|---:|---:|---:|---:|
| **Sec. 5.3.1** ($\alpha=90^\circ, \theta=-3^\circ$) | $352,90 / 208,26 = 1,694$ | $0,215$ | **CMM1** ($8,52\text{ mm}$) | **68 958,4** | 66 550 | +3,62% | 51 166 |
|  |  |  | **CMM2** ($4,24\text{ mm}$) | **43 893,0** | 42 508 | +3,26% | 37 455 |
|  |  |  | **CMM3** ($3,64\text{ mm}$) | **38 979,4** | 38 059 | +2,42% | 30 007 |
| **Sec. 5.3.2** (`Chambers_CMM_4`, $\alpha=110^\circ, \theta=-4^\circ$) | $2,36 / 401,06 = 0,0059$ | $1,537$ | **CMM4** ($2,70\text{ mm}$) | **29 524,5** | 28 783 | +2,58% | 27 782 |

---

## 6. Figuras Geradas pelos Scripts Python

### 6.1. Malha Descontínua e Condições de Contorno (`inputs/*_mesh_bc.png`)
- [`inputs/Chambers_CMM_4_mesh_bc.png`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/inputs/Chambers_CMM_4_mesh_bc.png) — Corpo de prova CMM ($\alpha = 110^\circ$)
- [`inputs/s531_cmm_90_mesh_bc.png`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/inputs/s531_cmm_90_mesh_bc.png) — Corpo de prova CMM ($\alpha = 90^\circ$)
- [`inputs/Dowling_11.6_variable_ampl_mesh_bc.png`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/inputs/Dowling_11.6_variable_ampl_mesh_bc.png) — Placa com trinca central (Dowling 11.6)
- [`inputs/s521_single_edge_mesh_bc.png`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/inputs/s521_single_edge_mesh_bc.png) — Placa AISI 4340 com 1 trinca de borda
- [`inputs/s522_double_edge_mesh_bc.png`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/inputs/s522_double_edge_mesh_bc.png) — Placa AISI 4340 com 2 trincas de borda

### 6.2. Pós-processamento dos Resultados (`examples/*.png`)
- **Trajetória da Trinca (Fig. 5.26):** [`examples/Chambers_CMM_4_fig526_trajetoria.png`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/examples/Chambers_CMM_4_fig526_trajetoria.png)
- **Malha Final com Trinca Propagada (Fig. 5.27):** [`examples/Chambers_CMM_4_fig527_malha_final.png`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/examples/Chambers_CMM_4_fig527_malha_final.png)
- **Configuração Deformada Ampliada 40x (Fig. 5.28):** [`examples/Chambers_CMM_4_fig528_deformada_40x.png`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/examples/Chambers_CMM_4_fig528_deformada_40x.png)
- **Avanço da Trinca vs Ciclos (Fig. 5.29):** [`examples/Chambers_CMM_4_fig529_ciclos_vs_trinca.png`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/examples/Chambers_CMM_4_fig529_ciclos_vs_trinca.png)
- **Taxa de Propagação $da/dN \times \Delta K$ (Fig. 5.30):** [`examples/Chambers_CMM_4_fig530_dadn_vs_dK.png`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/examples/Chambers_CMM_4_fig530_dadn_vs_dK.png)
- **Painel Resumo CMM4:** [`examples/Chambers_CMM_4_dashboard_resumo.png`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/examples/Chambers_CMM_4_dashboard_resumo.png)
- **Painel Resumo Dowling 11.6:** [`examples/Dowling_11.6_variable_ampl_dashboard_resumo.png`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/examples/Dowling_11.6_variable_ampl_dashboard_resumo.png)
- **Espectro Rainflow Dowling 11.6:** [`examples/Dowling_11.6_variable_ampl_rainflow_ciclos.png`](file:///c:/Users/Adm/OneDrive%20-%20unb.br/Antigravity/CrackPropagation2D/examples/Dowling_11.6_variable_ampl_rainflow_ciclos.png)
