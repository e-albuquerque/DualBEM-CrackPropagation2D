!-----------------------------------------------------------------------
! PROGRAM: fracture_2d
!> @brief Programa principal para análise elástica bidimensional de fratura e
!>        propagação de trincas por fadiga utilizando o Método dos Elementos
!>        de Contorno Dual (Dual-BEM) com elementos quadráticos descontínuos.
!>
!> Fluxo de trabalho:
!>   1. Gmsh (.geo -> .msh) gera a malha contínua das curvas geométricas;
!>   2. `preproc.py` + `input_data.py` leem o `.msh` via `meshio` e geram o
!>      arquivo `.dat` com elementos quadráticos descontínuos;
!>   3. `fracture_2d` resolve o problema Dual-BEM passo a passo e exporta
!>      os resultados para pós-processamento com `plot_results.py`.
!-----------------------------------------------------------------------
program fracture_2d
  use constants_mod
  use model_mod
  use assembler_mod
  use sif_mod
  use fatigue_prop_mod
  use io_mod
  !$ use omp_lib, only: omp_get_max_threads
  implicit none

  character(len=512)    :: input_file_name
  real(dp), allocatable :: a(:,:), f(:)
  real(dp), allocatable :: rf_amax(:), rf_R(:)
  real                  :: t_sf, t_ss, t_ip, t_cr
  integer               :: c1, c2, c_rate
  integer               :: n, steps, rf_cycles, k_sif, m_sif, n_threads
  logical               :: bound_recovery, variable_ampl
  real(dp)              :: aq, amax

  t_sf = 0.0; t_ss = 0.0; t_ip = 0.0; t_cr = 0.0
  n_threads = 1
  !$ n_threads = omp_get_max_threads()
  call system_clock(count_rate=c_rate)

  write(*, '(78("*"))')
  write(*, '(A)') "Fracture 2D - Dual-BEM (Discontinuous Quadratic Elements Edition)"
  write(*, '(A)') "Propagacao Bidimensional de Trincas por Fadiga"
  write(*, '(A,I0,/)') "Threads OpenMP disponiveis: ", n_threads

  ! 1. Obtém o nome do arquivo de entrada (.dat) via argumento de linha de comando ou stdin
  if (command_argument_count() >= 1) then
    call get_command_argument(1, input_file_name)
  else
    write(*, '(A)') "Informe o caminho do arquivo de entrada (.dat):"
    read(*, '(A)') input_file_name
  end if
  input_file_name = adjustl(input_file_name)

  ! 2. Lê os dados pré-processados do modelo
  call input_from_text_file(trim(input_file_name), steps, variable_ampl)

  ! 3. Grava o eco dos dados de entrada (<id>_input_data.txt)
  call initial_data_output()

  ! 4. Contagem de ciclos Rainflow (amplitude variável) ou fator de Walker (amplitude constante)
  if (variable_ampl) then
    write(*, '(1X,A)') "Executando contagem de ciclos Rainflow..."
    call rainflow(aq, amax)
    call get_rainflow_cycles(rf_cycles, rf_amax, rf_R)
    call cycles_output(rf_cycles, rf_amax, rf_R)
    if (allocated(rf_amax)) deallocate(rf_amax, rf_R)
    call set_parameter_R(0.0_dp)
  else
    call get_aq(aq)
    amax = 1.0_dp / aq
    call set_parameter_R(0.0_dp)
  end if

  ! 5. Montagem do sistema de equações no Passo 0
  write(*, '(1X,A)') "Montando o sistema de equacoes Dual-BEM (Passo 0)..."
  call system_clock(c1)
  allocate(a(2*bn_num, 2*bn_num), f(2*bn_num))
  call assembler(aq, a, f)
  call system_clock(c2)
  t_sf = real(c2 - c1) / real(max(c_rate, 1))

  ! 6. Solução do sistema linear via LAPACK (DGESV)
  write(*, '(1X,A)') "Resolvendo o sistema linear..."
  call system_clock(c1)
  call solve_linear_system(a, f)
  call system_clock(c2)
  t_ss = real(c2 - c1) / real(max(c_rate, 1))
  deallocate(a)

  ! 7. Recuperação dos deslocamentos, trações e tensões no contorno
  call set_results_on_boundary_nodes(f)
  deallocate(f)
  call get_K_key(k_sif)
  call get_method_key(m_sif)
  bound_recovery = .not. (k_sif == 2 .and. m_sif == 1)
  call evaluate_stress_strain_on_boundary_nodes(bound_recovery)

  ! 8. Avaliação dos pontos internos no Passo 0
  if (cr_num > 0) call create_internal_points(20, 2)
  if (in_num > 0) then
    call system_clock(c1)
    call evaluate_internal_points()
    call system_clock(c2)
    t_ip = real(c2 - c1) / real(max(c_rate, 1))
  end if

  ! 9. Avaliação dos Fatores de Intensidade de Tensão e propagação incremental
  if (cr_num > 0) then
    write(*, '(1X,A)') "Avaliando Fatores de Intensidade de Tensao (Passo 0)..."
    call system_clock(c1)
    call stress_intensity_factor_evaluation()
    call propagation_angle()
    call equivalent_k()
    call system_clock(c2)
    t_cr = real(c2 - c1) / real(max(c_rate, 1))

    ! Arquivos de saída do Passo 0
    call results_output(t_sf, t_ss, t_ip, t_cr, 0)
    call crack_geometry_output(0)
    call results_per_crack_output(0, amax)
    call fatigue_results_per_crack_output(0)
    call output_file_for_insane(0)
    call original_internal_points()

    ! Laço incremental de propagação da trinca
    if (steps > 0) then
      write(*, '(1X,A,I0,A)') "Propagando trinca(s) por ", steps, " passos..."
      do n = 1, steps
        call crack_propagation(n)

        call create_internal_points(20, 2)
        if (in_num > 0) then
          call system_clock(c1)
          call evaluate_internal_points()
          call system_clock(c2)
          t_ip = real(c2 - c1) / real(max(c_rate, 1))
        end if

        if (mod(n, 2) == 0) then
          call number_of_cycles()
          call fatigue_results_per_crack_output(n)
        end if

        call results_output(t_sf, t_ss, t_ip, t_cr, n)
        call results_per_crack_output(n, amax)
        call crack_geometry_output(n)
        call output_file_for_insane(n)
        call original_internal_points()
      end do
      write(*, '(1X,A)') "Propagacao concluida com sucesso!"
    end if
  else
    ! Problema elastostático sem trinca: grava resultados de contorno/pontos internos e XML
    call results_output(t_sf, t_ss, t_ip, t_cr, 0)
    call output_file_for_insane(0)
  end if

  write(*, '(78("*"))')
end program fracture_2d
