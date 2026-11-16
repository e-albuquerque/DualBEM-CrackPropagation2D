!-----------------------------------------------------------------------
! MODULE: fatigue_prop_mod
!> @brief Propagação incremental de trincas e previsão de vida em fadiga:
!>        1. Critério da Máxima Tensão Circunferencial (MTC) para direção de
!>           propagação (theta_p) e algoritmo preditor-corretor iterativo de Portela;
!>        2. Adição automática de novos elementos quadráticos descontínuos nas
!>           duas faces da trinca (bie = 2 e bie = 3) a cada incremento;
!>        3. Integração da Lei de Paris-Walker pela Regra de Simpson 1/3;
!>        4. Contagem de ciclos Rainflow para carregamentos de amplitude variável.
!-----------------------------------------------------------------------
module fatigue_prop_mod
  use constants_mod
  use quadrature_mod
  use model_mod
  use assembler_mod
  use sif_mod
  implicit none

  type(boundary_elem), allocatable, private :: lbe(:)
  type(boundary_node), allocatable, private :: lbn(:)
  integer, private  :: tnbe = 0, tnbn = 0, k_eff_key = 2
  real(dp), private :: R_load = 0.0_dp, C0 = 0.0_dp, m_exp = 0.0_dp, gama = 1.0_dp
  real(dp), private :: K_th = 0.0_dp, K_IC = 0.0_dp
  real(dp), private :: Da_0 = 0.0_dp, Da_max = 0.0_dp, r_Da = 1.0_dp
  real(dp), private :: l_tip_elem = 0.0_dp, length_0 = 0.0_dp

  ! Acumuladores de vida em fadiga (até 1000 trincas, 2 pontas por trinca)
  real(dp), private :: N_step(1000, 2)    = 0.0_dp
  real(dp), private :: N_total(1000, 2)   = 0.0_dp
  real(dp), private :: delta_a_i(1000, 2) = 0.0_dp
  real(dp), private :: y_simp(1000, 2, 3) = 0.0_dp

  ! Espectro de amplitude variável para contagem Rainflow
  real(dp), allocatable, private :: alpha(:)
  integer, private               :: rf_cycles = 0
  real(dp), allocatable, private :: rf_alpha_max(:), rf_R(:)

contains

  !=====================================================================
  ! GETTERS E SETTERS DOS PARÂMETROS DE FADIGA E PROPAGAÇÃO
  !=====================================================================

  subroutine set_initial_delta_a(p)
    real(dp), intent(in) :: p
    Da_0 = p
  end subroutine set_initial_delta_a

  subroutine set_max_delta_a(p)
    real(dp), intent(in) :: p
    Da_max = p
  end subroutine set_max_delta_a

  subroutine set_r_Da(p)
    real(dp), intent(in) :: p
    r_Da = p
  end subroutine set_r_Da

  subroutine set_parameter_C(p)
    real(dp), intent(in) :: p
    C0 = p
  end subroutine set_parameter_C

  subroutine set_parameter_m(p)
    real(dp), intent(in) :: p
    m_exp = p
  end subroutine set_parameter_m

  subroutine get_parameter_m(p)
    real(dp), intent(out) :: p
    p = m_exp
  end subroutine get_parameter_m

  subroutine set_parameter_R(p)
    real(dp), intent(in) :: p
    R_load = p
  end subroutine set_parameter_R

  subroutine get_parameter_R(p)
    real(dp), intent(out) :: p
    p = R_load
  end subroutine get_parameter_R

  subroutine set_parameter_gamma(p)
    real(dp), intent(in) :: p
    gama = p
  end subroutine set_parameter_gamma

  subroutine get_parameter_gamma(p)
    real(dp), intent(out) :: p
    p = gama
  end subroutine get_parameter_gamma

  subroutine set_K_th(p)
    real(dp), intent(in) :: p
    K_th = p
  end subroutine set_K_th

  subroutine set_K_IC(p)
    real(dp), intent(in) :: p
    K_IC = p
  end subroutine set_K_IC

  subroutine get_K_IC(p)
    real(dp), intent(out) :: p
    p = K_IC
  end subroutine get_K_IC

  subroutine set_k_eff_key(p)
    integer, intent(in) :: p
    k_eff_key = p
  end subroutine set_k_eff_key

  subroutine get_aq(aq)
    real(dp), intent(out) :: aq
    aq = (1.0_dp - R_load)**gama
  end subroutine get_aq

  subroutine get_number_of_cycles_per_step(i, j, n_cyc)
    integer, intent(in)   :: i, j
    real(dp), intent(out) :: n_cyc
    n_cyc = N_step(i, j)
  end subroutine get_number_of_cycles_per_step

  subroutine get_number_of_cycles_total(i, j, n_cyc)
    integer, intent(in)   :: i, j
    real(dp), intent(out) :: n_cyc
    n_cyc = N_total(i, j)
  end subroutine get_number_of_cycles_total

  !=====================================================================
  ! CONTAGEM DE CICLOS RAINFLOW (CARREGAMENTO DE AMPLITUDE VARIÁVEL)
  !=====================================================================

  subroutine allocate_alpha(n)
    integer, intent(in) :: n
    if (allocated(alpha)) deallocate(alpha)
    allocate(alpha(n))
  end subroutine allocate_alpha

  subroutine set_alpha(i, a_val)
    integer, intent(in)  :: i
    real(dp), intent(in) :: a_val
    alpha(i) = a_val
  end subroutine set_alpha

  subroutine get_rainflow_cycles(cycles_out, alpha_max_out, R_out)
    integer, intent(out)               :: cycles_out
    real(dp), allocatable, intent(out) :: alpha_max_out(:), R_out(:)
    cycles_out = rf_cycles
    if (rf_cycles > 0) then
      allocate(alpha_max_out(rf_cycles), R_out(rf_cycles))
      alpha_max_out = rf_alpha_max
      R_out         = rf_R
    end if
  end subroutine get_rainflow_cycles

  !---------------------------------------------------------------------
  !> @brief Algoritmo Rainflow de contagem de ciclos (E548 / Dowling),
  !>        retornando o fator equivalente de Walker `aq` e o fator de pico `amax`.
  !---------------------------------------------------------------------
  subroutine rainflow(aq, amax)
    real(dp), intent(out) :: aq, amax
    integer :: i, j, m, cycles
    real(dp), allocatable :: aux_alpha_max(:), aux_R(:), aux_alpha(:)

    call rearrange_alpha()

    allocate(aux_alpha_max(100000), aux_R(100000))
    j      = size(alpha)
    i      = 1
    cycles = 0

    do while (j /= 1)
      m = size(alpha)
      if (abs(alpha(i) - alpha(i+1)) <= abs(alpha(i+1) - alpha(i+2))) then
        cycles                = cycles + 1
        aux_alpha_max(cycles) = max(alpha(i), alpha(i+1))
        aux_R(cycles)         = min(alpha(i), alpha(i+1)) / max(alpha(i), alpha(i+1))
        j                     = j - 2
        allocate(aux_alpha(m))
        aux_alpha = alpha
        deallocate(alpha)
        allocate(alpha(j))
        if (i == 1) then
          alpha(1:j) = aux_alpha(3:m)
        else
          alpha(1:i-1) = aux_alpha(1:i-1)
          alpha(i:j)   = aux_alpha(i+2:m)
        end if
        deallocate(aux_alpha)
        i = 1
      else
        i = i + 1
      end if
    end do

    rf_cycles = cycles
    if (allocated(rf_alpha_max)) deallocate(rf_alpha_max)
    if (allocated(rf_R))         deallocate(rf_R)
    allocate(rf_alpha_max(cycles), rf_R(cycles))
    rf_alpha_max = aux_alpha_max(1:cycles)
    rf_R         = aux_R(1:cycles)
    deallocate(aux_alpha_max, aux_R, alpha)

    call alpha_q_alpha_max(cycles, rf_alpha_max, rf_R, aq, amax)
  end subroutine rainflow

  subroutine rearrange_alpha()
    integer :: i, j, k, m
    real(dp), allocatable :: aux(:)
    real(dp) :: x

    x = alpha(1)
    j = 1
    do i = 1, size(alpha)
      if (abs(x) < abs(alpha(i))) then
        x = alpha(i)
        j = i
      end if
    end do

    m = size(alpha)
    allocate(aux(m + 1))
    k          = m - j + 1
    aux(1:k)   = alpha(j:m)
    aux(k+1:m) = alpha(1:j-1)
    aux(m+1)   = aux(1)
    deallocate(alpha)
    allocate(alpha(m + 1))
    alpha = aux
    deallocate(aux)
  end subroutine rearrange_alpha

  subroutine alpha_q_alpha_max(cycles, alpha_max, R_vec, aq, amax)
    integer, intent(in)   :: cycles
    real(dp), intent(in)  :: alpha_max(:), R_vec(:)
    real(dp), intent(out) :: aq, amax
    real(dp) :: alpha_sum, x_max
    integer  :: i

    alpha_sum = 0.0_dp
    do i = 1, cycles
      alpha_sum = alpha_sum + (alpha_max(i) * (1.0_dp - R_vec(i))**gama)**m_exp
    end do
    aq = (alpha_sum / real(cycles, dp))**(1.0_dp / m_exp)

    x_max = maxval(alpha_max(1:cycles))
    amax  = x_max / aq
  end subroutine alpha_q_alpha_max

  !=====================================================================
  ! PROPAGAÇÃO DE TRINCA E CÁLCULO DE VIDA EM FADIGA
  !=====================================================================

  !---------------------------------------------------------------------
  !> @brief Executa um passo de propagação incremental da trinca (etapa `n`),
  !>        incluindo a correção iterativa de direção de Portela.
  !---------------------------------------------------------------------
  subroutine crack_propagation(n)
    integer, intent(in) :: n
    real(dp) :: theta, aux_x1, aux_x2
    real(dp), allocatable :: beta_i(:), beta_i1(:), theta_p(:), coord(:,:)
    real(dp), allocatable :: a(:,:), f(:), h(:,:), g(:,:)
    real(dp), allocatable :: delta_k_eff_1(:,:), delta_k_eff_2(:,:), delta_a(:,:)
    logical :: bound_recovery, iteration
    logical, allocatable :: propagation_key(:)
    integer :: i, j, k, aux_be, aux_bn, b, c, d, tips_num, n_el, k_sif, m_sif
    type(boundary_elem) :: last_elem

    call get_K_key(k_sif)
    call get_method_key(m_sif)

    ! Salva a malha do passo anterior (antes de acrescentar o novo incremento)
    tnbe = be_num
    tnbn = bn_num
    allocate(lbe(tnbe), lbn(tnbn))
    lbe = be_list
    lbn = bn_list

    k      = 1
    aux_x1 = 0.0_dp
    aux_x2 = 0.0_dp

    call get_tips_num(tips_num)
    allocate(beta_i1(tips_num), beta_i(tips_num), theta_p(tips_num))
    allocate(coord(tips_num, 4), propagation_key(tips_num))
    allocate(delta_k_eff_1(cr_num, 2), delta_k_eff_2(cr_num, 2), delta_a(cr_num, 2))
    beta_i1 = 0.0_dp
    beta_i  = 0.0_dp

    call propagation_angle()

    b = 1
    do i = 1, cr_num
      do j = 1, cr_list(i)%tips_num
        if (cr_list(i)%keq(j) < K_th) then
          call error_message_2int_stop("Trinca ", i, " abaixo do limiar de propagacao na ponta ", j)
        end if

        last_elem  = be_list(cr_list(i)%ext_el_num(2*j))
        theta      = atan2(last_elem%bp_2%y - last_elem%bp_1%y, last_elem%bp_2%x - last_elem%bp_1%x)
        theta_p(b) = cr_list(i)%theta(j) + theta
        if (n == 1) then
          l_tip_elem = distance_only_2d(last_elem%bp_1, last_elem%bp_2)
        end if
        coord(b, 1) = last_elem%bp_2%x
        coord(b, 2) = last_elem%bp_2%y

        call get_delta_k_eff(i, j, delta_k_eff_1(i, j))
        if (mod(n, 2) == 1) then
          y_simp(i, j, 1) = 1.0_dp / (C0 * (delta_k_eff_1(i, j))**m_exp)
          if (j == 1) then
            length_0 = cr_list(i)%length / real(cr_list(i)%tips_num, dp)
          end if
        end if

        call get_delta_a_n_el(i, j, l_tip_elem, n, delta_a(i, j), n_el)
        cr_list(i)%length = cr_list(i)%length + delta_a(i, j)

        if (j == 2) then
          aux_x1 = coord(b, 1)
          aux_x2 = coord(b, 2)
        end if
        b = b + 1
      end do
    end do

    ! Laço iterativo corretor de direção (Portela et al., 1993)
    iteration       = .true.
    propagation_key = .true.

    do while (iteration)
      beta_i = beta_i1
      d      = 0
      b      = 1
      do i = 1, cr_num
        do j = 1, cr_list(i)%tips_num
          d = d + 1
          if (propagation_key(d)) then
            theta_p(d) = theta_p(d) + beta_i(d)
            if (j == 1) then
              coord(d, 3) = coord(d, 1) + delta_a(i, j) * cos(theta_p(d))
              coord(d, 4) = coord(d, 2) + delta_a(i, j) * sin(theta_p(d))
              if (k == 1) then
                cr_list(i)%ext_el_num(2) = tnbe + b * n_el
                cr_list(i)%ext_el_num(3) = tnbe + b * n_el + 1
                b = b + 1
              end if
            else if (j == 2) then
              coord(d, 3) = aux_x1
              coord(d, 4) = aux_x2
              coord(d, 1) = aux_x1 + delta_a(i, j) * cos(theta_p(d))
              coord(d, 2) = aux_x2 + delta_a(i, j) * sin(theta_p(d))
              if (k == 1) then
                cr_list(i)%ext_el_num(1) = tnbe + b * n_el + 1
                cr_list(i)%ext_el_num(4) = tnbe + 2 * b * n_el
                b = b + 1
              end if
            end if
          end if
        end do
        b = b + 1
      end do

      ! Atualiza a malha com os novos elementos quadráticos descontínuos da trinca
      aux_be = tnbe + 2 * d * n_el
      aux_bn = tnbn + 6 * d * n_el
      deallocate(be_list, bn_list)
      allocate(be_list(aux_be), bn_list(aux_bn))
      be_list(1:tnbe) = lbe
      bn_list(1:tnbn) = lbn
      be_num          = tnbe
      bn_num          = tnbn

      do c = 1, d
        call add_discontinuous_crack_elements(coord(c, 1), coord(c, 2), coord(c, 3), coord(c, 4), &
                                              n_el, last_elem%lamb_1, last_elem%lamb_2)
      end do

      call symetric_nodes(bn_num, tnbn, d)

      ! Resolve o problema elástico para a configuração estendida da trinca
      allocate(a(2*bn_num, 2*bn_num), f(2*bn_num))
      if (allocated(h)) deallocate(h)
      if (allocated(g)) deallocate(g)
      allocate(h(2*bn_num, 2*bn_num), g(2*bn_num, 2*bn_num))

      call assembler_for_propagation(tnbe, tnbn, a, f, h, g)
      call solve_linear_system(a, f)
      deallocate(a)
      call set_results_on_boundary_nodes(f)
      deallocate(f)

      if (k_sif == 2 .and. m_sif == 1) then
        bound_recovery = .false.
        call evaluate_stress_strain_on_boundary_nodes(bound_recovery)
      end if
      call stress_intensity_factor_evaluation()
      call propagation_angle()
      call equivalent_k()

      d = 0
      do i = 1, cr_num
        do j = 1, cr_list(i)%tips_num
          d          = d + 1
          beta_i1(d) = 0.5_dp * cr_list(i)%theta(j)
        end do
      end do

      if (k == 1) then
        do i = 1, tips_num
          propagation_key(i) = (abs(beta_i1(i)) > 1.0e-3_dp)
        end do
      else
        do i = 1, tips_num
          propagation_key(i) = (abs(beta_i1(i)) < abs(beta_i(i)) .and. &
                                abs(beta_i1(i)) > 1.0e-3_dp .and. propagation_key(i))
        end do
      end if

      iteration = any(propagation_key)
      k         = k + 1
    end do

    ! Recupera tensões/deformações no contorno (O(N)) após convergir a direção
    if (.not. (k_sif == 2 .and. m_sif == 1)) then
      bound_recovery = .true.
      call evaluate_stress_strain_on_boundary_nodes(bound_recovery)
    end if

    ! Avalia Delta K_eff após convergir a direção do incremento
    do i = 1, cr_num
      do j = 1, cr_list(i)%tips_num
        call get_delta_k_eff(i, j, delta_k_eff_2(i, j))
        if (mod(n, 2) == 1) then
          y_simp(i, j, 2) = 1.0_dp / (C0 * (delta_k_eff_2(i, j))**m_exp)
        else
          y_simp(i, j, 3) = 1.0_dp / (C0 * (delta_k_eff_2(i, j))**m_exp)
        end if
      end do
    end do

    call set_initial_h(h)
    call set_initial_g(g)
    deallocate(h, g, lbe, lbn, beta_i1, beta_i, theta_p, coord, propagation_key)
    deallocate(delta_k_eff_1, delta_k_eff_2, delta_a)
  end subroutine crack_propagation

  !---------------------------------------------------------------------
  !> @brief Integra o número de ciclos N pela Regra 1/3 de Simpson a cada
  !>        par de incrementos de propagação.
  !---------------------------------------------------------------------
  subroutine number_of_cycles()
    integer :: i, j
    do i = 1, cr_num
      do j = 1, cr_list(i)%tips_num
        N_step(i, j)  = delta_a_i(i, j) * (y_simp(i, j, 1) + 4.0_dp*y_simp(i, j, 2) + y_simp(i, j, 3)) / 3.0_dp
        N_total(i, j) = N_total(i, j) + N_step(i, j)
      end do
    end do
  end subroutine number_of_cycles

  subroutine get_tips_num(tips_num)
    integer, intent(out) :: tips_num
    integer :: i
    tips_num = 0
    do i = 1, cr_num
      tips_num = tips_num + cr_list(i)%tips_num
    end do
  end subroutine get_tips_num

  subroutine get_delta_k_eff(i, j, delta_k_eff)
    integer, intent(in)   :: i, j
    real(dp), intent(out) :: delta_k_eff
    real(dp) :: ang, dki, dkii

    dki  = cr_list(i)%ki(j)
    dkii = cr_list(i)%kii(j)
    if (k_eff_key == 1) then
      ! Critério de Tanaka
      delta_k_eff = (dki**4 + 8.0_dp * dkii**4)**0.25_dp
    else
      ! Critério da Máxima Tensão Circunferencial (MTC)
      ang         = cr_list(i)%theta(j)
      delta_k_eff = dki * (cos(0.5_dp*ang))**3 - 3.0_dp * dkii * sin(0.5_dp*ang) * (cos(0.5_dp*ang))**2
    end if
  end subroutine get_delta_k_eff

  subroutine get_delta_a_n_el(i, j, l_ref, n, delta_a, n_el)
    integer, intent(in)   :: i, j, n
    real(dp), intent(in)  :: l_ref
    real(dp), intent(out) :: delta_a
    integer, intent(out)  :: n_el
    real(dp) :: aux_n_el

    if (n == 1 .or. n == 2) then
      delta_a         = Da_0
      delta_a_i(i, j) = delta_a
    else
      if (mod(n, 2) == 1) then
        delta_a = r_Da * delta_a_i(i, j)
      else
        delta_a = delta_a_i(i, j)
      end if
      if (delta_a > Da_max) delta_a = Da_max
      delta_a_i(i, j) = delta_a
    end if

    aux_n_el = delta_a / l_ref
    n_el     = merge(1, nint(aux_n_el), aux_n_el < 1.0_dp)
  end subroutine get_delta_a_n_el

  !---------------------------------------------------------------------
  !> @brief Atualiza a matriz `cr_list(i)%nodes_num` pareando os novos nós
  !>        coincidentes da Face 1 (bie=2) e da Face 2 (bie=3).
  !---------------------------------------------------------------------
  subroutine symetric_nodes(bn_num_cur, tnbn_old, d)
    integer, intent(in) :: bn_num_cur, tnbn_old, d
    integer :: c, i, j
    integer :: aux1, aux2, aux3, aux4, aux5, aux6, aux7, aux8
    integer, allocatable :: nodes_num_temp(:,:)

    c = 0
    do i = 1, cr_num
      aux1 = (bn_num_cur - tnbn_old) / (2 * d)
      aux2 = size(cr_list(i)%nodes_num, 1)
      aux3 = aux1 + aux2
      aux4 = tnbn_old + 2*aux1 + c*aux1 + 1
      aux5 = tnbn_old + 4*aux1 + c*aux1
      aux6 = cr_list(i)%tips_num * aux1 + aux2
      aux7 = tnbn_old + c*aux1 + 1
      aux8 = tnbn_old + 2*aux1 + c*aux1

      allocate(nodes_num_temp(aux2, 2))
      nodes_num_temp = cr_list(i)%nodes_num
      deallocate(cr_list(i)%nodes_num)
      allocate(cr_list(i)%nodes_num(aux6, 2))

      if (cr_list(i)%tips_num == 2) then
        do j = 1, aux1
          cr_list(i)%nodes_num(j, 1) = aux4
          cr_list(i)%nodes_num(j, 2) = aux5
          aux4 = aux4 + 1
          aux5 = aux5 - 1
        end do
        cr_list(i)%nodes_num(aux1+1:aux3, :) = nodes_num_temp
        deallocate(nodes_num_temp)
        do j = aux3 + 1, aux6
          cr_list(i)%nodes_num(j, 1) = aux7
          cr_list(i)%nodes_num(j, 2) = aux8
          aux7 = aux7 + 1
          aux8 = aux8 - 1
        end do
      else
        cr_list(i)%nodes_num(1:aux2, :) = nodes_num_temp
        deallocate(nodes_num_temp)
        do j = aux2 + 1, aux6
          cr_list(i)%nodes_num(j, 1) = aux7
          cr_list(i)%nodes_num(j, 2) = aux8
          aux7 = aux7 + 1
          aux8 = aux8 - 1
        end do
      end if
      c = c + 2 * cr_list(i)%tips_num
    end do
  end subroutine symetric_nodes

  !---------------------------------------------------------------------
  !> @brief Calcula o ângulo local de propagação theta_p pelo critério da
  !>        Máxima Tensão Circunferencial (Erdogan & Sih, 1963).
  !---------------------------------------------------------------------
  subroutine propagation_angle()
    integer  :: i, j
    real(dp) :: theta_p, k_i, k_ii

    do i = 1, cr_num
      do j = 1, cr_list(i)%tips_num
        k_i  = cr_list(i)%ki(j)
        k_ii = cr_list(i)%kii(j)
        if (abs(k_ii) < 1.0e-15_dp) then
          theta_p = 0.0_dp
        else if (k_ii > 0.0_dp) then
          theta_p = 2.0_dp * atan(0.25_dp * ((k_i / k_ii) - sqrt((k_i / k_ii)**2 + 8.0_dp)))
        else
          theta_p = 2.0_dp * atan(0.25_dp * ((k_i / k_ii) + sqrt((k_i / k_ii)**2 + 8.0_dp)))
        end if
        cr_list(i)%theta(j) = theta_p
      end do
    end do
  end subroutine propagation_angle

  !---------------------------------------------------------------------
  !> @brief Calcula o fator de intensidade de tensão equivalente K_eq.
  !---------------------------------------------------------------------
  subroutine equivalent_k()
    real(dp) :: k_i, k_ii, theta_p
    integer  :: i, j

    do i = 1, cr_num
      do j = 1, cr_list(i)%tips_num
        k_i     = cr_list(i)%ki(j)
        k_ii    = cr_list(i)%kii(j)
        theta_p = cr_list(i)%theta(j)
        cr_list(i)%keq(j) = k_i * (cos(0.5_dp*theta_p))**3 - &
                            3.0_dp * k_ii * sin(0.5_dp*theta_p) * (cos(0.5_dp*theta_p))**2
      end do
    end do
  end subroutine equivalent_k

  !---------------------------------------------------------------------
  !> @brief Gera `nbe` elementos quadráticos descontínuos livres de tração
  !>        na Face 1 (bie = 2, de p1 para p2) e `nbe` elementos na Face 2
  !>        (bie = 3, de p2 para p1) durante o avanço da trinca.
  !---------------------------------------------------------------------
  subroutine add_discontinuous_crack_elements(x1, y1, x2, y2, nbe, lb1, lb2)
    real(dp), intent(in) :: x1, y1, x2, y2, lb1, lb2
    integer, intent(in)  :: nbe
    type(point_2d) :: p1, p2
    type(point_2d), allocatable :: pts(:), pts2(:)
    integer :: j, n_pts

    p1%x = x1; p1%y = y1
    p2%x = x2; p2%y = y2

    call divide_straight_line_quad(p1, p2, nbe, pts)

    ! Face 1 da trinca: EIC de Deslocamento (bie = 2)
    call append_discont_elements_from_points(pts, nbe, lb1, lb2, 2)

    ! Face 2 da trinca (orientação oposta): EIC Hipersingular de Tração (bie = 3)
    n_pts = size(pts)
    allocate(pts2(n_pts))
    do j = 1, n_pts
      pts2(n_pts - j + 1) = pts(j)
    end do
    call append_discont_elements_from_points(pts2, nbe, lb2, lb1, 3)

    deallocate(pts, pts2)
  end subroutine add_discontinuous_crack_elements

  subroutine append_discont_elements_from_points(pts, nbe, lb1, lb2, bie)
    type(point_2d), intent(in) :: pts(:)
    integer, intent(in)        :: nbe, bie
    real(dp), intent(in)       :: lb1, lb2
    type(point_2d) :: el_pts(3)
    integer :: i, k

    do i = 1, nbe
      el_pts(1) = pts(2*i - 1)
      el_pts(2) = pts(2*i)
      el_pts(3) = pts(2*i + 1)

      be_num = be_num + 1
      be_list(be_num)%bp_1      = el_pts(1)
      be_list(be_num)%bp_2      = el_pts(3)
      be_list(be_num)%nodes_num = 3
      be_list(be_num)%lamb_1    = lb1
      be_list(be_num)%lamb_2    = lb2

      ! Nó 1 (eta = -lb1)
      bn_num = bn_num + 1
      bn_list(bn_num)%point    = cartesian_point_from_intrp_func_1d(el_pts, -lb1)
      be_list(be_num)%incid(1) = bn_num

      ! Nó 2 (eta = 0)
      bn_num = bn_num + 1
      bn_list(bn_num)%point    = el_pts(2)
      be_list(be_num)%incid(2) = bn_num

      ! Nó 3 (eta = +lb2)
      bn_num = bn_num + 1
      bn_list(bn_num)%point    = cartesian_point_from_intrp_func_1d(el_pts, +lb2)
      be_list(be_num)%incid(3) = bn_num

      ! Condição de contorno livre de tração nas faces da trinca
      do k = bn_num - 2, bn_num
        bn_list(k)%bie    = bie
        bn_list(k)%bc_dir = 0
        bn_list(k)%bc     = [1, 1]
        bn_list(k)%disp   = 0.0_dp
        bn_list(k)%trac   = 0.0_dp
        bn_list(k)%h1     = 0.0_dp
        bn_list(k)%h2     = 0.0_dp
        bn_list(k)%strain = 0.0_dp
        bn_list(k)%stress = 0.0_dp
      end do
    end do
  end subroutine append_discont_elements_from_points

end module fatigue_prop_mod
