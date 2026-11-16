# Fracture 2D — Propagação Bidimensional de Trincas por Fadiga (Dual-BEM)

Implementação em Fortran 95 migrada para **GNU Fortran (`gfortran`)** e bibliotecas numéricas open-source (**Netlib LAPACK/BLAS**), referente à dissertação de mestrado de **Beatriz Borges Faria Fonseca** (UFMG, 2023), disponível neste repositório em [`DissertaçãoBeatrizFonseca.pdf`](DissertaçãoBeatrizFonseca.pdf):

> *Modelagem computacional da propagação bidimensional de trincas por fadiga utilizando o Método dos Elementos de Contorno Dual.*

---

## 1. O que o Código Faz

O programa **Fracture 2D** resolve problemas de elasticidade linear bidimensional (estado plano de tensões ou de deformações) para sólidos isotrópicos de região única contendo uma ou mais trincas (internas ou de borda), calculando fatores de intensidade de tensões ($K_I$ e $K_{II}$), trajetória incremental de propagação em modo I puro ou modo misto e vida em fadiga em número de ciclos ($N$).

Com base na dissertação, a formulação numérica reúne os seguintes métodos:

- **Método dos Elementos de Contorno Dual (Dual-BEM):** contorna a degeneração do sistema quando as duas faces da trinca coincidem geometricamente utilizando a **Equação Integral de Contorno (EIC) de deslocamentos** nos nós do contorno externo e em uma das faces da trinca, e a **EIC hipersingular de tensões/trações** na face oposta, com soluções fundamentais 2D de Kelvin (`kelvin_fund_sol_2d.f95`).
- **Discretização do Contorno e das Trincas:** elementos quadráticos contínuos ou semidescontínuos no contorno externo (usualmente com $\lambda_1 = \lambda_2 = 2/3$ nas quinas e em encontros com bocas de trinca) e elementos quadráticos descontínuos nas faces das trincas.
- **Integração Numérica Especializada:**
  - *Regular / quase singular:* quadratura de Gauss-Legendre (`fracture_2d_regular_integration.f95`, tipicamente 20 pontos de Gauss nos exemplos da dissertação).
  - *Fracamente singular:* quadratura de Gauss-Laguerre (`fracture_2d_weakly_sing_integration.f95`).
  - *Forte e hipersingular:* **Método da Expansão Radial de Gao** (`fracture_2d_rad_exp_intgr.f95`, usando 3 termos de expansão polinomial nos exemplos da dissertação).
- **Cálculo dos Fatores de Intensidade de Tensão ($K_I, K_{II}$):**
  - *Correlação de Deslocamentos* (`K_key = 1`, 7 variantes em `fracture_2d_displ_correlation.f95`).
  - *Integral $J$ desacoplada* (`K_key = 2`, 3 estratégias de contorno circular em `fracture_2d_j_integral.f95`, com contorno circular de 20 elementos como padrão na dissertação).
- **Critério de Direção e Propagação Incremental:** critério da **Máxima Tensão Circunferencial (MTC)** para obtenção do ângulo de deflexão $\theta$, refinado pelo **método iterativo corretivo de Portela** (`propagation_angle` em `fracture_2d_crack_propagation.f95`), e avaliação do fator equivalente $K_{\text{eq}}$ pelo critério MTC ou de Tanaka (`k_eff_key`).
- **Fadiga sob Amplitude Constante e Variável:**
  - Integração da equação de **Paris** com influência da razão de tensões $R$ pelo modelo de **Walker**, aplicando o fator equivalente $(1 - R)^\gamma$ em `get_aq`, e integração numérica da vida $N$ pela **Regra de Simpson** (`number_of_cycles`) com comprimento de avanço adaptativo $\Delta a$ multiplicado por $r_{\Delta a}$ a cada par de etapas até o limite $\Delta a_{\max}$.
  - Para carregamentos de amplitude variável (`vari_ampl_key = 1`), o módulo `fracture_2d_rainflow.f95` aplica o algoritmo **Rainflow Cycle Counting**, extraindo cada ciclo $(\alpha_{\max}, R)$ e calculando o fator de tensão equivalente $\alpha_q$ e o fator de pico $\alpha_{\max}$ para verificação à fratura estática ($K_{\text{eq}} \ge K_{IC}$).
- **Critérios de Parada e Pós-processamento:** o processo iterativo encerra ao atingir o número máximo de passos, quando $K_{\text{eq}} \ge K_{IC}$ ou quando $\Delta K \le \Delta K_{\text{th}}$. Além das tabelas `.txt`, o solver gera arquivos `-STEP-NNN.xml` compatíveis com visualização gráfica no software **INSANE**.

---

## 2. Estrutura do Repositório

```text
CrackPropagation2D/
├── DissertaçãoBeatrizFonseca.pdf
├── Makefile                                  # Makefile raiz (delega para Módulos e main)
├── README.md
├── plot_mesh_bc.py                           # Script Python para visualizar malha e CCs (.png)
├── bin/
│   └── fracture_2d.exe                       # Executável compilado (gfortran + LAPACK/BLAS estático)
├── examples/                                 # Exemplos de gráficos gerados por plot_mesh_bc.py
│   ├── Chambers_CMM_4_mesh_bc.png
│   ├── Dowling_11.6_variable_ampl_mesh_bc.png
│   └── s521_single_edge_mesh_bc.png
├── Arquivos de entrada-20241119T180208Z-001/
│   └── Arquivos de entrada/
│       ├── Chambers_CMM_4.txt                # Exemplo da Seção 5.3.2 (CMM α = 110°)
│       └── Dowling_11.6_variable_ampl.txt    # Exemplo da Seção 5.2.3 (amplitude variável)
└── Módulos e main-20241119T180145Z-001/
    └── Módulos e main/
        ├── Fracture_2D_v3-0.f95              # Programa principal
        ├── Makefile                          # Regras de compilação em ordem de dependência
        ├── build/                            # Objetos (.o) e módulos compilados (.mod)
        └── modules/                          # 22 módulos Fortran 95
```

### Arquitetura dos Módulos (`modules/`)

| Arquivo | Função principal no solver |
|---|---|
| `Fracture_2D_v3-0.f95` | Laço principal: leitura de entrada, Rainflow (quando ativo), montagem do sistema, chamada ao solver linear, campos de contorno/internos, avaliação de $K_I/K_{II}$, propagação e medição de tempos. |
| `error_messages.f95` | Rotinas padronizadas de mensagens de erro e encerramento seguro (`stop 1`). |
| `geometry.f95` | Tipos geométricos (`point_2d`, `point_3d`) e utilitários de distância/geometria. |
| `vector_algebra.f95` | Operações vetoriais básicas (como norma euclidiana `vector_len`). |
| `matrix_algebra.f95` | Álgebra matricial e interface explícita Fortran 95 para `dgesv` da LAPACK via `solve_linear_system(a,f)`. |
| `interpolation_functions.f95` | Funções de interpolação 1D para elementos contínuos, semidescontínuos ($\lambda_1, \lambda_2$) e descontínuos. |
| `numerical_integration_points.f95` | Pontos e pesos de quadratura numérica (Gauss-Legendre e afins). |
| `kelvin_fund_sol_2d.f95` | Núcleos das soluções fundamentais 2D de Kelvin em deslocamento e tração. |
| `fracture_2d_model.f95` | Estruturas globais da malha, nós, condições de contorno, trincas e pontos internos. |
| `fracture_2d_input.f95` | Leitura e validação do arquivo `.txt` de entrada, geração automática da discretização de segmentos/arcos/polinômios, trincas e grade de pontos internos. |
| `fracture_2d_rainflow.f95` | Algoritmo de contagem de ciclos *Rainflow* para séries temporais de amplitude variável. |
| `fracture_2d_regular_integration.f95` | Integração regular e quase singular via quadratura de Gauss-Legendre. |
| `fracture_2d_weakly_sing_integration.f95` | Integração fracamente singular (analítica em lineares/constantes e Gauss-Laguerre em quadráticos). |
| `fracture_2d_rad_exp_intgr.f95` | Integração forte e hipersingular pelo Método da Expansão Radial de Gao. |
| `fracture_2d_assembler.f95` | Montagem das matrizes globais $[A]$ e vetor $\{f\}$ do sistema Dual-BEM. |
| `fracture_2d_set_bound_results.f95` | Atribuição dos deslocamentos e forças de superfície incógnitas nos nós do contorno. |
| `fracture_2d_internal_points.f95` | Cálculo de deslocamentos, deformações e tensões nos pontos internos. |
| `fracture_2d_create_internal_points.f95` | Geração auxiliar de pontos internos ao redor da ponta da trinca. |
| `fracture_2d_displ_correlation.f95` | Avaliação de $K_I$ e $K_{II}$ por correlação de deslocamentos (7 métodos). |
| `fracture_2d_j_integral.f95` | Avaliação desacoplada de $K_I$ e $K_{II}$ pela Integral $J$ em contorno circular. |
| `fracture_2d_stress_intens_factor.f95` | Despacho entre Correlação de Deslocamentos (`K_key = 1`) e Integral $J$ (`K_key = 2`). |
| `fracture_2d_crack_propagation.f95` | Critério MTC, algoritmo corretivo de Portela, cálculo de $K_{\text{eq}}$, fator $(1-R)^\gamma$, incremento adaptativo $\Delta a$ e regra de Simpson para $N$. |
| `fracture_2d_output.f95` | Geração dos relatórios `.txt` de parâmetros, contorno, trinca, fadiga, ciclos Rainflow e malhas XML para o INSANE. |

---

## 3. Migração de Intel Fortran / MKL para `gfortran` + LAPACK/BLAS

Para eliminar qualquer dependência do compilador comercial Intel Fortran e da biblioteca Intel MKL, foram realizadas as seguintes adaptações e correções estruturais:

1. **Substituição da Intel MKL por Netlib LAPACK/BLAS (`matrix_algebra.f95`, `Fracture_2D_v3-0.f95`, `fracture_2d_crack_propagation.f95`):**
   - Adicionada interface explícita Fortran 95 para a rotina padrão LAPACK `dgesv(n,nrhs,a,lda,ipiv,b,ldb,info)`.
   - Implementada a subrotina `solve_linear_system(a,f)` em `matrix_algebra.f95`, validando dimensões, alocando `integer, allocatable :: ipiv(:)` e `real(kind=8), allocatable :: rhs(:,:)`, executando `call dgesv(...)`, checando `info` e atualizando `f`.
2. **Correção de bugs revelados na portabilidade para `gfortran`:**
   - **`ipiv` declarado incorretamente como `real(kind=8)`** na subrotina de ajuste iterativo de ângulo em `fracture_2d_crack_propagation.f95` (substituído por `call solve_linear_system(a,f)`).
   - **Variável `n` não inicializada no passo 0** ao gravar a geometria inicial em `Fracture_2D_v3-0.f95` (`call crack_geometry_output(n)` corrigido para `call crack_geometry_output(0)`).
   - **Erro de índice na coluna `N_2_total` (`fracture_2d_output.f95`, linha 665):** chamava `get_number_of_cycles_total(i,1,N_2_total)` em vez de `get_number_of_cycles_total(i,2,N_2_total)`, replicando os ciclos da ponta 1 na ponta 2.
   - **Inicialização explícita de acumuladores e temporizadores:** inicializados `t_sf`, `t_ss`, `t_ip`, `t_cr` em `Fracture_2D_v3-0.f95` e `N_step`, `N_total`, `delta_a_i`, `y` em `fracture_2d_crack_propagation.f95`.
   - **Eliminação de avisos `-Wmaybe-uninitialized` e `-Wdo-subscript`** em `fracture_2d_weakly_sing_integration.f95` (`eta0 = 0._8`) e `fracture_2d_j_integral.f95`.
   - **Codificação UTF-8:** conversão de comentários em Windows-1252 em `fracture_2d_crack_propagation.f95` e `fracture_2d_rainflow.f95`.
3. **Remoção de extensões não portáteis e automação de execução:**
   - Removidas todas as chamadas `system("pause")` em `error_messages.f95` e `Fracture_2D_v3-0.f95`; erros encerram imediatamente com `stop 1`.
   - Adicionado suporte opcional a **argumento de linha de comando** (`get_command_argument(1, input_file_name)` com buffer de 256 caracteres), mantendo o modo interativo via `stdin` quando nenhum argumento é informado.
4. **Sistema de build portável (`Makefile` raiz e `Makefile` do solver):**
   - Compilação em ordem topológica de módulos com `-O2 -Wall -Wextra` e linkagem com `-llapack -lblas` (padrão `STATIC=1` com `-static` no Windows).

---

## 4. Pré-requisitos e Como Compilar

### Windows (MSYS2 UCRT64 + MinGW-w64)

No ambiente **MSYS2 UCRT64**, os pacotes verificados para compilação são:
- `mingw-w64-ucrt-x86_64-gcc-fortran` (testado com GNU Fortran 16.2.0)
- `mingw-w64-ucrt-x86_64-lapack` (3.12.1 — fornece `liblapack` e `libblas`)
- `mingw-w64-ucrt-x86_64-make` (`mingw32-make` 4.4.1)

Caso precise instalá-los no terminal **MSYS2 UCRT64**:

```bash
pacman -S mingw-w64-ucrt-x86_64-gcc-fortran mingw-w64-ucrt-x86_64-lapack mingw-w64-ucrt-x86_64-make
```

> **Nota importante no Windows (`PATH`):** Ao compilar a partir do PowerShell ou CMD, certifique-se de que `C:\msys64\ucrt64\bin` apareça **no início** do `PATH` para evitar conflito de DLLs com outras ferramentas (como MATLAB, MiKTeX ou Emacs) que possam causar falha silenciosa no compilador `f951.exe`:

```powershell
$env:PATH = "C:\msys64\ucrt64\bin;$env:PATH"
mingw32-make clean
mingw32-make
```

Por padrão no Windows (`STATIC ?= 1`), a compilação usa `-static -llapack -lblas`. Conforme verificado com `objdump -p bin\fracture_2d.exe`, o executável resultante não depende de DLLs externas do `gfortran` ou da LAPACK/BLAS (apenas `KERNEL32.dll` e bibliotecas nativas do Windows CRT), podendo ser executado diretamente em qualquer terminal sem alterar o `PATH`.

Para compilação dinâmica no Windows:

```powershell
$env:PATH = "C:\msys64\ucrt64\bin;$env:PATH"
mingw32-make clean
mingw32-make STATIC=0
```

### Linux / macOS / POSIX

Com `gfortran`, `make` e pacotes `liblapack-dev` / `libblas-dev` instalados:

```bash
make STATIC=0 clean
make STATIC=0
```

O binário gerado ficará em `bin/fracture_2d` (ou `bin/fracture_2d.exe` no Windows).

---

## 5. Como Executar

O programa grava todos os arquivos de saída no **diretório de trabalho atual**, usando o identificador `Analysis Id` definido dentro do arquivo de entrada. Recomenda-se executar sempre a partir de uma subpasta dedicada aos resultados.

### Execução por argumento de linha de comando (recomendado)

```powershell
New-Item -ItemType Directory -Force -Path meus_resultados\chambers | Out-Null
Push-Location meus_resultados\chambers
..\..\bin\fracture_2d.exe "..\..\Arquivos de entrada-20241119T180208Z-001\Arquivos de entrada\Chambers_CMM_4.txt"
Pop-Location
```

### Execução interativa / redirecionada via entrada padrão (`stdin`)

```powershell
New-Item -ItemType Directory -Force -Path meus_resultados\dowling | Out-Null
Push-Location meus_resultados\dowling
"..\..\Arquivos de entrada-20241119T180208Z-001\Arquivos de entrada\Dowling_11.6_variable_ampl.txt" | ..\..\bin\fracture_2d.exe
Pop-Location
```

### Estrutura básica do arquivo de entrada e dica de modelagem

Cada seção do arquivo de entrada (ver exemplos completos em [`Chambers_CMM_4.txt`](Arquivos%20de%20entrada-20241119T180208Z-001/Arquivos%20de%20entrada/Chambers_CMM_4.txt) e [`Dowling_11.6_variable_ampl.txt`](Arquivos%20de%20entrada-20241119T180208Z-001/Arquivos%20de%20entrada/Dowling_11.6_variable_ampl.txt)) define:
1. Identificação (`Analysis name`, `Analysis Id`), constantes elásticas ($E$, $\nu$), estado plano (`1` = EPT, `2` = EPD) e termos da expansão de Gao (`npow`).
2. Contorno externo (`bk = 1` reta, `2` arco, `3` polinômio cúbico, `4` pontos dados), conectividade (`ck`), parâmetros de descontinuidade (`lb1`, `lb2`) e condições de contorno (`bck = 0` deslocamento, `1` tração, `2` apoio elástico).
3. Condições em pontos particulares e definição das trincas (`cn`, `ck`, elementos por face, `lb1`, `lb2`).
4. Grade de pontos internos (`ipk = 0..4`).
5. Método para cálculo dos FITs (`K_key = 1` Correlação de Deslocamentos ou `2` Integral $J$, mais `method_key` e parâmetros geométricos do caminho circular).
6. Parâmetros de fadiga e avanço (`Da_0`, `r_da`, `Da_max`, `steps`, razão de carga $R$, constantes $C_0$, $m$, expoente de Walker $\gamma$, limiar $\Delta K_{\text{th}}$, tenacidade $K_{IC}$, chave `k_eff_key` e opção de amplitude variável `vari_ampl_key` seguida do histórico de picos/vales quando `vari_ampl_key = 1`).

> **Regra prática para trincas de borda (*edge cracks*):** Nunca posicione um nó contínuo diretamente sobre a boca da trinca no contorno externo. Divida a aresta externa exatamente no ponto de abertura em dois segmentos com extremidades descontínuas/semidescontínuas (`lb = 2/3 ≈ 0.66666667`), exatamente como feito em `Chambers_CMM_4.txt`. Caso contrário, a sobreposição nodal na abertura produz matriz singular e valores `NaN`.

### Arquivos de saída gerados

| Padrão de arquivo | Conteúdo |
|---|---|
| `<id>_input_data.txt` | Eco dos dados de entrada, coordenadas dos nós, condições de contorno e pontos internos gerados. |
| `<id>_results_step_NNN.txt` | Deslocamentos, trações e tensões de contorno, além dos campos nos pontos internos no passo `NNN`. |
| `<id>_results_crack_NNN.txt` | Histórico de comprimento e valores de $K_I$, $K_{II}$ e $K_{\text{eq}}$ por ponta para cada etapa (indica `**** Keq greater than K_IC!****` ao atingir a fratura crítica). |
| `<id>_geometry_crack_NNN.txt` | Coordenadas sequenciais dos segmentos da trinca ao longo da propagação. |
| `<id>_fatigue_results_crack_NNN.txt` | Comprimento e contagem de ciclos por etapa (`N_1_step`, `N_1_total`, e `N_2_step`, `N_2_total` para trincas com duas pontas). |
| `<id>_cycles.txt` | Tabela de ciclos $(\alpha_{\max}, R)$ extraídos pelo *Rainflow* (quando `vari_ampl_key = 1`). |
| `<id>-STEP-NNN.xml` | Arquivo XML de malha e campos em cada etapa para pós-processamento gráfico no **INSANE**. |

---

## 6. Validação Numérica (Comparação com a Dissertação)

O executável `bin/fracture_2d.exe` compilado com `gfortran` 16.2.0 + Netlib LAPACK/BLAS foi validado em **8 problemas** abrangendo os capítulos de validação de FITs (Seção 3.1.3) e propagação por fadiga em modo I, modo misto e amplitude variável (Seções 5.1, 5.2 e 5.3).

### 6.1. Seção 3.1.3.1 — Placa quadrada com trinca de borda ($10 \times 10$, modo I puro)

Placa sob tração unitária ($\sigma = 1$), $w = 10$, variando $a/w$ de $0,2$ a $0,6$, com contorno dividido na boca da trinca ($\lambda = 2/3$). Comparação do FIT normalizado $K_I / (\sigma \sqrt{\pi a})$ frente à solução de referência de Civelek e Erdogan (Tabela 3.4 da dissertação):

| $a/w$ | Integral $J$ (20 elem.) | Ref. (Tab. 3.4) | Erro $J$ (%) | Correl. Desloc. (D–E) | Erro CD (%) |
|---:|---:|---:|---:|---:|---:|
| 0,2 | 1,488 | 1,488 | 0,00% | 1,556 | +4,56% |
| 0,3 | 1,848 | 1,848 | 0,00% | 1,940 | +4,97% |
| 0,4 | 2,322 | 2,324 | -0,11% | 2,449 | +5,38% |
| 0,5 | 3,000 | 3,010 | -0,35% | 3,183 | +5,73% |
| 0,6 | 4,120 | 4,152 | -0,78% | 4,403 | +6,04% |

A Integral $J$ reproduz os valores da Tabela 3.4 da dissertação em todas as razões $a/w$ testadas (erro máximo de `0,78%`), superando em precisão a correlação de deslocamentos (`4,56%` a `6,04%`), exatamente como discutido no Capítulo 3 da dissertação.

### 6.2. Seção 3.1.3.2 — Trinca interna inclinada a $45^\circ$ em modo misto ($10 \times 20$)

Placa retangular ($w = 10$, altura $20$, tração unitária $\sigma = 1$) com trinca centrada inclinada a $45^\circ$ ($2a/w = 0,2$ a $0,6$), normalizando $K_I$ e $K_{II}$ por $\sigma \sqrt{\pi (a_{\text{total}}/2)}$ via Integral $J$:

| $2a/w$ | $K_{I,\text{norm}}$ (código) | $K_{I,\text{ref}}$ | Erro $K_I$ (%) | $K_{II,\text{norm}}$ (código) | $K_{II,\text{ref}}$ | Erro $K_{II}$ (%) |
|---:|---:|---:|---:|---:|---:|---:|
| 0,2 | 0,520 | 0,518 | +0,44% | 0,513 | 0,507 | +1,19% |
| 0,3 | 0,543 | 0,541 | +0,38% | 0,522 | 0,516 | +1,13% |
| 0,4 | 0,575 | 0,572 | +0,50% | 0,534 | 0,529 | +1,02% |
| 0,5 | 0,616 | 0,612 | +0,60% | 0,551 | 0,546 | +0,94% |
| 0,6 | 0,666 | 0,661 | +0,76% | 0,572 | 0,567 | +0,94% |

### 6.3. Seção 5.1 — Placa retangular com trinca interna e amplitude constante ($h=3$, $b=0,5$, $2a_0=0,1$)

Comparação para 24 incrementos de propagação sob amplitude constante ($\sigma_{\max} = 10$, $R = 0$, $C_0 = 1,5\times 10^{-8}$, $m = 2,25$), avaliando $K_I$ contra a fórmula de Isida e o número de ciclos acumulado contra a integração direta da Equação (5.4):

| Esquema de avanço | $2a_f$ (código) | $2a_f$ (diss.) | $K_{I,0}$ (código / Isida) | $K_{I,f}$ (código / Isida) | $N_f$ (código) | $N_f$ (Eq. 5.4) | Erro $N_f$ |
|---|---:|---:|---:|---:|---:|---:|---:|
| $\Delta a_0 = 0,01; \; r = 1,1$ | 0,9554 | 0,955 | 4,0169 / 3,9870 | 48,2691 / 47,5443 | 229 439,6 | 229 938,9 | 0,22% |
| $\Delta a_0 = 0,0125; \; r = 1,05$ | 0,8959 | 0,896 | 4,0169 / 3,9870 | 29,6619 / 29,9411 | 228 880,6 | 229 311,9 | 0,19% |

Após a correção da rotina de saída em `fracture_2d_output.f95`, ambas as pontas apresentam contagem acumulada idêntica por simetria ($N_{1,\text{total}} = N_{2,\text{total}} = 229\,439,6$ e $228\,880,6$).

### 6.4. Seções 5.2.1 e 5.2.2 — Placas de aço AISI 4340 com amplitude constante ($R = 1/3$, $K_{IC} = 130\text{ MPa}\sqrt{\text{m}}$)

Ensaios numéricos para placa de largura $76\text{ mm}$ com uma trinca de borda ($a_0 = 2,000\text{ mm}$) e com duas trincas simétricas de borda ($a_0 = 1,000\text{ mm}$), comparados aos comprimentos finais da dissertação e às soluções analíticas das Equações (5.8), (5.9) e (5.10):

| Caso | $a_0$ (mm) | $a_f$ obtido (mm) | $a_f$ diss. (mm) | $\Delta K_0$ (código / analít.) | $\Delta K_f$ (código / analít.) | $K_{\max,f}$ ($\text{MPa}\sqrt{\text{m}}$) | $N_f$ (código) | $N_f$ (ref.) | Erro $N_f$ |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **Sec. 5.2.1** (1 trinca borda, 38 passos) | 2,000 | 12,232 | 12,232 | 31,576 / 31,616 | 88,258 / 89,493 | 132,39 | 25 601,1 | 25 070,4 | 2,12% |
| **Sec. 5.2.2** (2 trincas borda, 44 passos) | 1,000 | 15,281 | 15,281 | 21,276 / 22,051 | 86,695 / 86,969 | 130,04 | 54 349,7 / 54 349,2 | 53 457,6 | 1,67% |

Em ambos os casos, o solver interrompe exatamente no passo em que $K_{\text{eq}}$ ultrapassa $K_{IC} = 130\text{ MPa}\sqrt{\text{m}}$ (`**** Keq greater than K_IC!****`), reproduzindo os comprimentos finais críticos (`12,232 mm` e `15,281 mm`) até a terceira casa decimal e preservando simetria numérica superior a `99,998%` entre as duas trincas da Seção 5.2.2.

### 6.5. Seção 5.2.3 — Trinca interna e amplitude variável (`Dowling_11.6_variable_ampl.txt`)

- **Contagem Rainflow:** dos blocos históricos de carregamento são extraídos **166 ciclos**, resultando em $a_q = 1,414222$ ($\sigma_q = 219,298 \times 1,414222 = 310,14\text{ MPa}$, com diferença de apenas `0,37%` em relação ao valor arredondado `311,3 MPa` reportado na dissertação) e pico máximo $a_{\max} = 2,4$.
- **Fator de Intensidade de Tensões ($K_I$):** ao longo de todos os 47 estados (`Step 0` a `Step 46`, de $2a_0 = 2,000\text{ mm}$ a $2a_f = 31,637\text{ mm}$), o $K_I$ de pico calculado pelo Dual-BEM (`29,6106` até `131,3814 MPa·√m`) concorda com a solução analítica de Feddersen $\alpha_{\max}\sigma_{\text{ref}}\sqrt{\pi a \sec(\pi a/W)}$ com erro máximo absoluto de **`0,332%`**.
- **Comprimento crítico e vida:** a propagação cessa automaticamente no passo 46 por sobrecarga ($K_{\text{eq}} = 131,38 > K_{IC} = 130\text{ MPa}\sqrt{\text{m}}$) com $2a_f = 31,637\text{ mm}$ (dissertação: `31,64 mm`) e vida acumulada de **`238 810,1` ciclos**.

### 6.6. Seções 5.3.1 e 5.3.2 — Corpo de prova *Compact Mixed Mode* (CMM, aço Jethete M152 a $550^\circ\text{C}$)

Modelagem dos ensaios experimentais de modo misto de Chambers et al. (1991), avaliando a razão inicial $K_I/K_{II}$, a guinada para propagação em modo I puro e o número de ciclos comparado aos resultados numéricos e experimentais da Tabela 5.9 e da Seção 5.3.2:

| Configuração | $a_0 \to a_f$ ($\Delta a$) | $(K_I / K_{II})_0$ | $K_{II}$ no passo 1 | Ponto de comparação ($\Delta a$) | $N$ (código) | $N_{\text{num}}$ (diss.) | Dif. vs diss. | $N_{\text{exp}}$ (Chambers) |
|---|---|---|---|---:|---:|---:|---:|---:|
| **Sec. 5.3.1** ($\alpha=90^\circ, \theta=-3^\circ$) | $6,262 \to 16,494\text{ mm}$ ($10,232\text{ mm}$) | $373,72 / 204,40 = 1,828$ (diss. $\approx 1,8$) | $0,199$ | **CMM1** ($8,52\text{ mm}$) | 66 129,1 | 66 550 | 0,63% | 51 166 |
|  |  |  |  | **CMM2** ($4,24\text{ mm}$) | 42 182,5 | 42 508 | 0,77% | 37 455 |
|  |  |  |  | **CMM3** ($3,64\text{ mm}$) | 37 474,1 | 38 059 | 1,54% | 30 007 |
| **Sec. 5.3.2** (`Chambers_CMM_4`, $\alpha=110^\circ, \theta=-4^\circ$) | $6,262 \to 16,494\text{ mm}$ ($10,232\text{ mm}$) | $34,15 / 396,09 = 0,0862$ (diss. $\approx 0,086$) | $1,462$ | **$\Delta a = 2,70\text{ mm}$** | 28 392,8 | 28 783 | 1,36% | 27 782 |

Na configuração $\alpha = 110^\circ$ (fortemente dominada por modo II inicial), o critério MTC acoplado ao algoritmo corretivo de Portela captura a mudança brusca de trajetória logo no primeiro incremento ($K_{II}$ cai de `396,09` para `1,46`), totalizando `77 176,6` ciclos ao final dos 38 incrementos (`16,494 mm` vs `16,5 mm` na dissertação) e concordando dentro de `2,20%` com o valor experimental de Chambers et al. (1991) em $\Delta a = 2,7\text{ mm}$.

---

## 7. Visualização da Malha e Condições de Contorno (`plot_mesh_bc.py`)

O script [`plot_mesh_bc.py`](plot_mesh_bc.py) (baseado em `preproc.py`) lê o arquivo `<id>_input_data.txt` gerado pelo solver e exporta uma imagem `.png` contendo:

- **Elementos quadráticos de contorno** (linhas pretas, reconstruindo curvas retas, circulares e semidescontínuas via $\eta \in [-1, 1]$ e nó médio);
- **Faces da trinca** (destacadas em laranja pela detecção automática de pares coincidentes do Dual-BEM);
- **Pontos internos** (pontos verdes, quando presentes);
- **Condições de contorno** com marcadores e vetores proporcionais para valores não nulos:
  - **Deslocamentos prescritos (Dirichlet):** vermelho (`rx` para $u_x/u_n$, `rP` para $u_y/u_t$);
  - **Trações prescritas (Neumann):** azul (`bx` para $t_x/t_n$, `bP` para $t_y/t_t$);
  - **Apoios elásticos (molas):** verde (`gx` para $k_x/k_n$, `gP` para $k_y/k_t$).

### Como usar

Requer Python 3 com `numpy` e `matplotlib`:

```powershell
python plot_mesh_bc.py caminho\para\<id>_input_data.txt
```

Opções disponíveis:
- `-o`, `--output`: caminho do arquivo `.png` de saída (padrão: `<id>_mesh_bc.png` na mesma pasta do arquivo de entrada);
- `--dpi`: resolução em DPI (padrão: `300`);
- `--no-nodes`: oculta os marcadores dos nós e exibe apenas malha, trinca e setas não nulas;
- `--show`: abre a janela interativa do Matplotlib após salvar.

Exemplos gerados encontram-se em [`examples/Chambers_CMM_4_mesh_bc.png`](examples/Chambers_CMM_4_mesh_bc.png), [`examples/Dowling_11.6_variable_ampl_mesh_bc.png`](examples/Dowling_11.6_variable_ampl_mesh_bc.png) e [`examples/s521_single_edge_mesh_bc.png`](examples/s521_single_edge_mesh_bc.png).

---

## 8. Pós-processamento e Gráficos de Resultados (`plot_results.py`)

O script [`plot_results.py`](plot_results.py) realiza a leitura automática de todos os arquivos de saída gerados pelo programa Fortran durante a simulação de propagação e fadiga por elementos de contorno dual (MECD), gerando gráficos de alta resolução (300 DPI) correspondentes às figuras de resultados do Capítulo 5 da dissertação:

### Arquivos de Saída Lidos

| Arquivo | Conteúdo e Função |
| :--- | :--- |
| `<id>_geometry_crack_NNN.txt` | Histórico das coordenadas $(X, Y)$ da ponta da trinca e comprimento $a$ em cada passo $0 \le n \le n_p$. |
| `<id>_results_crack_NNN.txt` | Fatores de intensidade de tensão ($K_I, K_{II}, K_{\text{eq}}$) nas pontas em cada passo e verificação contra $K_{IC}$. |
| `<id>_fatigue_results_crack_NNN.txt` | Comprimento da trinca $a$, ciclos do incremento $\Delta N$ e ciclos totais acumulados $N_{\text{total}}$ (para pontas 1 e 2). |
| `<id>_cycles.txt` | Tabela de ciclos extraídos pelo algoritmo *Rainflow* para amplitude variável (ciclo, $\alpha_{\max}$, razão $R$). |
| `<id>_results_step_NNN.txt` | Campos nodais de contorno em cada passo (deslocamentos, trações, tensões principais e de von Mises). |
| `<id>-STEP-NNN.xml` | Modelo geométrico e solução completa no formato XML do INSANE (`BemSteadyState` / `BemModel`). |

---

### Gráficos Gerados (Figuras da Dissertação)

1. **Figura 5.26: Trajetória de Propagação da Trinca** (`<id>_fig526_trajetoria.png`)
   - Trajetória global no sistema de coordenadas da peça $(X, Y)$ com marcação dos passos $P_0, P_4, \dots, P_n$.
   - Trajetória relativa $(\Delta X, \Delta Y)$ a partir da ponta inicial da trinca, idêntica ao gráfico da p. 165 da dissertação.
2. **Figura 5.27: Configuração Final da Malha** (`<id>_fig527_malha_final.png`)
   - Reconstrução paramétrica precisa de todos os elementos lineares e quadráticos do modelo XML no passo final da análise.
3. **Figura 5.28: Configuração Deformada Ampliada** (`<id>_fig528_deformada_<scale>x.png`)
   - Sobreposição da malha indeformada (preto) e deformada (vermelho), amplificada por fator configurável (ex.: 40x).
   - Painel lateral com zoom e detalhe focado na abertura e deslizamento da fissura (COD e CSD).
4. **Figura 5.29: Avanço da Trinca vs Número de Ciclos** (`<id>_fig529_ciclos_vs_trinca.png`)
   - Curvas de vida útil à fadiga: $N \times \Delta a$ e $N \times a$, permitindo comparar o crescimento da trinca e ciclos totais.
5. **Figura 5.30: Taxa de Propagação $da/dN \times \Delta K$ em Escala Log-Log** (`<id>_fig530_dadn_vs_dK.png`)
   - Taxas numéricas calculadas por diferenciação dos passos e reta ajustada da Lei de Paris ($\log(da/dN) = \log(C) + m \log(\Delta K)$), extraindo automaticamente os expoentes $C$ e $m$ da simulação.
6. **Evolução dos FIT ($K_I, K_{II}, K_{\text{eq}}$) vs $a$** (`<id>_sif_evolucao.png`)
   - Correspondente às Figuras 5.3, 5.9, 5.13 e 5.19, ilustrando a transição de modo misto para modo I predominante.
7. **Espectro de Ciclos Rainflow** (`<id>_rainflow_ciclos.png`)
   - Histograma de amplitudes de pico $\alpha_{\max}$ e mapa de dispersão $R \times \alpha_{\max}$ (quando gerado via amplitude variável).
8. **Painel Resumo Integrado** (`<id>_dashboard_resumo.png`)
   - Painel 2x2 reunindo os 4 gráficos capitais da dissertação em uma única prancha executiva de apresentação.

---

### Como Usar

O script detecta automaticamente todos os arquivos a partir do diretório, de um arquivo ou de um prefixo:

```powershell
# Execução completa em um diretório de resultados (gera todas as figuras e o painel resumo):
python plot_results.py scratch/val_chambers

# Ajustar fator de escala de deformação (ex.: 50x) e resolução para 300 DPI:
python plot_results.py scratch/val_dowling --scale 50 --dpi 300

# Selecionar passo específico para visualização da deformada (ex.: passo 20):
python plot_results.py scratch/val_chambers --step 20 --scale 40

# Salvar figuras em uma pasta específica:
python plot_results.py scratch/val_chambers --output-dir examples --scale 40

# Gerar apenas figuras específicas:
python plot_results.py scratch/val_chambers --trajectory   # apenas Figura 5.26
python plot_results.py scratch/val_chambers --deformed    # apenas Figura 5.28
python plot_results.py scratch/val_chambers --fatigue     # apenas Figura 5.29
python plot_results.py scratch/val_chambers --dadn        # apenas Figura 5.30
python plot_results.py scratch/val_dowling --cycles       # apenas espectro Rainflow
```

### Exemplos Gerados na Pasta `examples/`

- **Trajetória da Trinca (Fig. 5.26):** [`examples/Chambers_CMM_4_fig526_trajetoria.png`](examples/Chambers_CMM_4_fig526_trajetoria.png)
- **Malha Final com Trinca (Fig. 5.27):** [`examples/Chambers_CMM_4_fig527_malha_final.png`](examples/Chambers_CMM_4_fig527_malha_final.png)
- **Malha Deformada 40x (Fig. 5.28):** [`examples/Chambers_CMM_4_fig528_deformada_40x.png`](examples/Chambers_CMM_4_fig528_deformada_40x.png)
- **Avanço $\Delta a \times N$ (Fig. 5.29):** [`examples/Chambers_CMM_4_fig529_ciclos_vs_trinca.png`](examples/Chambers_CMM_4_fig529_ciclos_vs_trinca.png)
- **Curva Log-Log $da/dN \times \Delta K$ (Fig. 5.30):** [`examples/Chambers_CMM_4_fig530_dadn_vs_dK.png`](examples/Chambers_CMM_4_fig530_dadn_vs_dK.png)
- **Histórico de $K_I, K_{II}, K_{\text{eq}}$:** [`examples/Chambers_CMM_4_sif_evolucao.png`](examples/Chambers_CMM_4_sif_evolucao.png)
- **Painel Resumo CMM:** [`examples/Chambers_CMM_4_dashboard_resumo.png`](examples/Chambers_CMM_4_dashboard_resumo.png)
- **Painel Resumo Amplitude Variável (Dowling):** [`examples/Dowling_11.6_variable_ampl_dashboard_resumo.png`](examples/Dowling_11.6_variable_ampl_dashboard_resumo.png)
- **Espectro Rainflow (Dowling):** [`examples/Dowling_11.6_variable_ampl_rainflow_ciclos.png`](examples/Dowling_11.6_variable_ampl_rainflow_ciclos.png)

