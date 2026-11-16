!-----------------------------------------------------------------------
! MODULE: sif_mod
!> @brief Cálculo dos Fatores de Intensidade de Tensão (K_I e K_II) nas
!>        pontas das trincas por dois métodos:
!>        1. Correlação de Deslocamentos (K_key = 1, 7 variantes de extrapolação);
!>        2. Integral J desacoplada em Modos I e II sobre contorno circular
!>           ao redor da ponta da trinca (K_key = 2, 3 variantes de integração).
!-----------------------------------------------------------------------
module sif_mod
  use constants_mod
  use quadrature_mod
  use model_mod
  use kelvin_sol_mod
  use integration_mod
  implicit none

  ! Tipo auxiliar para pontos do contorno circular da Integral J
  type :: j_integral_node
    type(point_2d) :: point
    real(dp)       :: stress(3) = 0.0_dp
    real(dp)       :: strain(3) = 0.0_dp
    real(dp)       :: grad(2)   = 0.0_dp
  end type j_integral_node

  integer, private  :: K_key            = 2
  integer, private  :: method_key       = 3
  integer, private  :: nd_correl        = 1
  integer, private  :: np_j             = 20
  integer, private  :: nbn_j            = 3
  integer, private  :: nbe_j            = 20
  real(dp), private :: lb1_j            = 2.0_dp / 3.0_dp
  real(dp), private :: lb2_j            = 2.0_dp / 3.0_dp
  integer, private  :: initial_node_num = 2
  integer, allocatable, private :: sym_node(:)

contains

  subroutine set_K_key(k)
    integer, intent(in) :: k
    K_key = k
  end subroutine set_K_key

  subroutine get_K_key(k)
    integer, intent(out) :: k
    k = K_key
  end subroutine get_K_key

  subroutine set_method_key(k)
    integer, intent(in) :: k
    method_key = k
  end subroutine set_method_key

  subroutine get_method_key(k)
    integer, intent(out) :: k
    k = method_key
  end subroutine get_method_key

  subroutine set_node(n)
    integer, intent(in) :: n
    nd_correl = n
  end subroutine set_node

  subroutine get_node(n)
    integer, intent(out) :: n
    n = nd_correl
  end subroutine get_node

  subroutine set_number_of_points(k)
    integer, intent(in) :: k
    np_j = k
  end subroutine set_number_of_points

  subroutine get_number_of_points(k)
    integer, intent(out) :: k
    k = np_j
  end subroutine get_number_of_points

  subroutine set_element_information(el, nd, lbd1, lbd2)
    integer, intent(in)  :: el, nd
    real(dp), intent(in) :: lbd1, lbd2
    nbe_j = el
    nbn_j = nd
    lb1_j = lbd1
    lb2_j = lbd2
  end subroutine set_element_information

  subroutine get_element_information(el, nd, lbd1, lbd2)
    integer, intent(out)  :: el, nd
    real(dp), intent(out) :: lbd1, lbd2
    el   = nbe_j
    nd   = nbn_j
    lbd1 = lb1_j
    lbd2 = lb2_j
  end subroutine get_element_information

  subroutine set_initial_node_num(k)
    integer, intent(in) :: k
    initial_node_num = k
  end subroutine set_initial_node_num

  subroutine get_initial_node_num(k)
    integer, intent(out) :: k
    k = initial_node_num
  end subroutine get_initial_node_num

  !---------------------------------------------------------------------
  !> @brief Despacha o cálculo de K_I e K_II conforme K_key.
  !---------------------------------------------------------------------
  subroutine stress_intensity_factor_evaluation()
    select case (K_key)
    case (1)
      call K_using_displ_correlation()
    case (2)
      call K_using_J_integral()
    case default
      call error_message_int_stop("Chave K_key invalida para calculo de FIT: ", K_key)
    end select
  end subroutine stress_intensity_factor_evaluation

  !=====================================================================
  ! 1. MÉTODO DA CORRELAÇÃO DE DESLOCAMENTOS (K_key = 1)
  !=====================================================================

  subroutine K_using_displ_correlation()
    integer  :: i
    real(dp) :: r, ki(2), kii(2)

    do i = 1, cr_num
      select case (method_key)
      case (1)
        call K_using_first_node(i, r, ki, kii)
      case (2)
        call K_using_mid_node(i, r, ki, kii)
      case (3)
        call K_using_last_node(i, r, ki, kii)
      case (4)
        call K_using_linear_extrap_last_and_mid_node(i, ki, kii)
      case (5)
        call K_using_linear_extrap_mid_and_first_node(i, ki, kii)
      case (6)
        call K_using_extrap_with_all_three_nodes(i, ki, kii)
      case (7)
        call K_using_any_node(i, ki, kii)
      case default
        call error_message_int_stop("Chave method_key invalida para correlacao de deslocamentos: ", method_key)
      end select
      cr_list(i)%ki  = ki
      cr_list(i)%kii = kii
    end do
  end subroutine K_using_displ_correlation

  subroutine K_using_first_node(i, r, ki, kii)
    integer, intent(in)   :: i
    real(dp), intent(out) :: r, ki(2), kii(2)
    integer :: j, k, nd1, nd2

    ki = 0.0_dp; kii = 0.0_dp
    do j = 1, cr_list(i)%tips_num
      if (j == 1) then
        k   = size(cr_list(i)%nodes_num, 1)
        nd1 = cr_list(i)%nodes_num(k, 1)
        nd2 = cr_list(i)%nodes_num(k, 2)
      else
        k   = 1
        nd1 = cr_list(i)%nodes_num(k, 2)
        nd2 = cr_list(i)%nodes_num(k, 1)
      end if
      call ki_and_kii(i, j, nd1, nd2, r, ki, kii)
    end do
  end subroutine K_using_first_node

  subroutine K_using_mid_node(i, r, ki, kii)
    integer, intent(in)   :: i
    real(dp), intent(out) :: r, ki(2), kii(2)
    integer :: j, k, nd1, nd2

    ki = 0.0_dp; kii = 0.0_dp
    do j = 1, cr_list(i)%tips_num
      if (j == 1) then
        k   = size(cr_list(i)%nodes_num, 1) - 1
        nd1 = cr_list(i)%nodes_num(k, 1)
        nd2 = cr_list(i)%nodes_num(k, 2)
      else
        k   = 2
        nd1 = cr_list(i)%nodes_num(k, 2)
        nd2 = cr_list(i)%nodes_num(k, 1)
      end if
      call ki_and_kii(i, j, nd1, nd2, r, ki, kii)
    end do
  end subroutine K_using_mid_node

  subroutine K_using_last_node(i, r, ki, kii)
    integer, intent(in)   :: i
    real(dp), intent(out) :: r, ki(2), kii(2)
    integer :: j, k, nd1, nd2

    ki = 0.0_dp; kii = 0.0_dp
    do j = 1, cr_list(i)%tips_num
      if (j == 1) then
        k   = size(cr_list(i)%nodes_num, 1) - 2
        nd1 = cr_list(i)%nodes_num(k, 1)
        nd2 = cr_list(i)%nodes_num(k, 2)
      else
        k   = 3
        nd1 = cr_list(i)%nodes_num(k, 2)
        nd2 = cr_list(i)%nodes_num(k, 1)
      end if
      call ki_and_kii(i, j, nd1, nd2, r, ki, kii)
    end do
  end subroutine K_using_last_node

  subroutine K_using_linear_extrap_last_and_mid_node(i, ki, kii)
    integer, intent(in)   :: i
    real(dp), intent(out) :: ki(2), kii(2)
    real(dp) :: r_2, ki_2(2), kii_2(2), r_3, ki_3(2), kii_3(2)

    call K_using_mid_node(i, r_2, ki_2, kii_2)
    call K_using_last_node(i, r_3, ki_3, kii_3)
    ki  = (-r_2 / (r_2 - r_3)) * (ki_2  - ki_3)  + ki_2
    kii = (-r_2 / (r_2 - r_3)) * (kii_2 - kii_3) + kii_2
  end subroutine K_using_linear_extrap_last_and_mid_node

  subroutine K_using_linear_extrap_mid_and_first_node(i, ki, kii)
    integer, intent(in)   :: i
    real(dp), intent(out) :: ki(2), kii(2)
    real(dp) :: r_2, ki_2(2), kii_2(2), r_1, ki_1(2), kii_1(2)

    call K_using_mid_node(i, r_2, ki_2, kii_2)
    call K_using_first_node(i, r_1, ki_1, kii_1)
    ki  = (-r_1 / (r_1 - r_2)) * (ki_1  - ki_2)  + ki_1
    kii = (-r_1 / (r_1 - r_2)) * (kii_1 - kii_2) + kii_1
  end subroutine K_using_linear_extrap_mid_and_first_node

  subroutine K_using_extrap_with_all_three_nodes(i, ki, kii)
    integer, intent(in)   :: i
    real(dp), intent(out) :: ki(2), kii(2)
    real(dp) :: r_3, ki_3(2), kii_3(2), r_2, ki_2(2), kii_2(2), r_1, ki_1(2), kii_1(2), ni(3)
    integer  :: j, k
    type(boundary_elem) :: last_elem

    ki = 0.0_dp; kii = 0.0_dp
    call K_using_last_node(i, r_3, ki_3, kii_3)
    call K_using_mid_node(i, r_2, ki_2, kii_2)
    call K_using_first_node(i, r_1, ki_1, kii_1)

    do j = 1, cr_list(i)%tips_num
      k = 2*j
      last_elem = be_list(cr_list(i)%ext_el_num(k))
      call intrp_func_1d(ni, 1.0_dp, last_elem%lamb_1, last_elem%lamb_2)
      ki(j)  = ni(1)*ki_3(j)  + ni(2)*ki_2(j)  + ni(3)*ki_1(j)
      kii(j) = ni(1)*kii_3(j) + ni(2)*kii_2(j) + ni(3)*kii_1(j)
    end do
  end subroutine K_using_extrap_with_all_three_nodes

  subroutine K_using_any_node(i, ki, kii)
    integer, intent(in)   :: i
    real(dp), intent(out) :: ki(2), kii(2)
    integer  :: j, k, nd1, nd2
    real(dp) :: r

    ki = 0.0_dp; kii = 0.0_dp
    do j = 1, cr_list(i)%tips_num
      if (j == 1) then
        k   = size(cr_list(i)%nodes_num, 1) - nd_correl + 1
        nd1 = cr_list(i)%nodes_num(k, 1)
        nd2 = cr_list(i)%nodes_num(k, 2)
      else
        k   = nd_correl
        nd1 = cr_list(i)%nodes_num(k, 2)
        nd2 = cr_list(i)%nodes_num(k, 1)
      end if
      call ki_and_kii(i, j, nd1, nd2, r, ki, kii)
    end do
  end subroutine K_using_any_node

  subroutine ki_and_kii(i, j, nd1, nd2, r, ki, kii)
    integer, intent(in)   :: i, j, nd1, nd2
    real(dp), intent(out) :: r
    real(dp), intent(inout) :: ki(2), kii(2)
    type(boundary_elem) :: last_elem
    real(dp) :: und1_1, und1_2, und2_1, und2_2, u1_1, u1_2, u2_1, u2_2, mu, kappa, normal(2)

    call get_shear_modulus(mu)
    call get_kappa(kappa)

    und1_1 = bn_list(nd1)%disp(1)
    und1_2 = bn_list(nd1)%disp(2)
    und2_1 = bn_list(nd2)%disp(1)
    und2_2 = bn_list(nd2)%disp(2)
    normal = get_normal_to_a_boundary_point(nd2)

    u1_1 = normal(2)*und1_1 - normal(1)*und1_2
    u1_2 = normal(1)*und1_1 + normal(2)*und1_2
    u2_1 = normal(2)*und2_1 - normal(1)*und2_2
    u2_2 = normal(1)*und2_1 + normal(2)*und2_2

    last_elem = be_list(cr_list(i)%ext_el_num(2*j))
    r         = distance_only_2d(bn_list(nd1)%point, last_elem%bp_2)
    ki(j)     = (u1_2 - u2_2) * mu * sqrt(two_pi / r) / (kappa + 1.0_dp)
    kii(j)    = (u1_1 - u2_1) * mu * sqrt(two_pi / r) / (kappa + 1.0_dp)
  end subroutine ki_and_kii

  !=====================================================================
  ! 2. MÉTODO DA INTEGRAL J EM CONTORNO CIRCULAR (K_key = 2)
  !=====================================================================

  subroutine K_using_J_integral()
    select case (method_key)
    case (1)
      call K_using_Simpson_rule()
    case (2)
      call K_using_Gauss_integration_points()
    case (3)
      call K_using_Gauss_integration_points_on_boundary()
    case default
      call error_message_int_stop("Chave method_key invalida para Integral J: ", method_key)
    end select
  end subroutine K_using_J_integral

  subroutine K_using_Simpson_rule()
    type(j_integral_node) :: cp, cp_sym, cp2, cp_sym2
    type(j_integral_node), allocatable :: j_int_list(:)
    type(boundary_elem) :: last_elem
    type(boundary_node) :: initial_node, final_node
    real(dp) :: J_i, J_ii, ang, rad, normal(2)
    real(dp) :: grad_i(2), grad_ii(2), W_i_normal, W_ii_normal, t_i(2), t_ii(2)
    real(dp) :: grad2_i(2), grad2_ii(2), W2_i_normal, W2_ii_normal, t2_i(2), t2_ii(2)
    real(dp) :: intgr1_i, intgr1_ii, intgr2_i, intgr2_ii
    integer  :: i, j, k, k_ii_sign

    do i = 1, cr_num
      do j = 1, cr_list(i)%tips_num
        J_i  = 0.0_dp
        J_ii = 0.0_dp
        call get_initial_node_and_last_element(i, j, last_elem, initial_node, final_node, normal)
        allocate(j_int_list(np_j))
        call Simpson_circular_line(last_elem%bp_2, initial_node%point, np_j, j_int_list, ang, rad)

        do k = 1, np_j
          if (k == 1) then
            call evaluate_first_and_last_node(j_int_list(k), initial_node, normal)
          else if (k == np_j) then
            call evaluate_first_and_last_node(j_int_list(k), final_node, normal)
          else
            call evaluate_j_integral_points(j_int_list(k), normal)
          end if
        end do

        do k = 1, np_j - 1
          cp      = j_int_list(k)
          cp_sym  = j_int_list(sym_node(k))
          call get_variables(cp, cp_sym, last_elem, normal, grad_i, grad_ii, W_i_normal, W_ii_normal, t_i, t_ii)
          cp2     = j_int_list(k+1)
          cp_sym2 = j_int_list(sym_node(k+1))
          call get_variables(cp2, cp_sym2, last_elem, normal, grad2_i, grad2_ii, W2_i_normal, W2_ii_normal, t2_i, t2_ii)

          intgr1_i  = W_i_normal  - t_i(1)*grad_i(1)   - t_i(2)*grad_i(2)
          intgr2_i  = W2_i_normal - t2_i(1)*grad2_i(1) - t2_i(2)*grad2_i(2)
          intgr1_ii = W_ii_normal  - t_ii(1)*grad_ii(1)   - t_ii(2)*grad_ii(2)
          intgr2_ii = W2_ii_normal - t2_ii(1)*grad2_ii(1) - t2_ii(2)*grad2_ii(2)

          J_i  = J_i  + (intgr1_i  + intgr2_i)  * ang * rad / 2.0_dp
          J_ii = J_ii + (intgr1_ii + intgr2_ii) * ang * rad / 2.0_dp
        end do

        call evaluate_sign_of_k_ii(initial_node, final_node, normal, k_ii_sign)
        call k_i_and_k_ii(i, j, k_ii_sign, J_i, J_ii)
        deallocate(j_int_list)
      end do
    end do
  end subroutine K_using_Simpson_rule

  subroutine K_using_Gauss_integration_points()
    type(j_integral_node) :: cp, cp_sym
    type(j_integral_node), allocatable :: j_int_list(:)
    type(boundary_elem) :: last_elem
    type(boundary_node) :: initial_node, final_node
    real(dp) :: J_i, J_ii, rad, normal(2)
    real(dp), allocatable :: ome(:)
    real(dp) :: grad_i(2), grad_ii(2), W_i_normal, W_ii_normal, t_i(2), t_ii(2)
    integer  :: i, j, k, k_ii_sign

    do i = 1, cr_num
      do j = 1, cr_list(i)%tips_num
        J_i  = 0.0_dp
        J_ii = 0.0_dp
        call get_initial_node_and_last_element(i, j, last_elem, initial_node, final_node, normal)
        allocate(j_int_list(np_j))
        call Gauss_circular_line(last_elem%bp_2, initial_node%point, np_j, j_int_list, rad, ome)

        !$omp parallel do default(shared) private(k) schedule(dynamic)
        do k = 1, np_j
          call evaluate_j_integral_points(j_int_list(k), normal)
        end do
        !$omp end parallel do

        do k = 1, np_j
          cp     = j_int_list(k)
          cp_sym = j_int_list(sym_node(k))
          call get_variables(cp, cp_sym, last_elem, normal, grad_i, grad_ii, W_i_normal, W_ii_normal, t_i, t_ii)
          J_i  = J_i  + (W_i_normal  - t_i(1)*grad_i(1)   - t_i(2)*grad_i(2))   * pi * rad * ome(k)
          J_ii = J_ii + (W_ii_normal - t_ii(1)*grad_ii(1) - t_ii(2)*grad_ii(2)) * pi * rad * ome(k)
        end do

        call evaluate_sign_of_k_ii(initial_node, final_node, normal, k_ii_sign)
        call k_i_and_k_ii(i, j, k_ii_sign, J_i, J_ii)
        deallocate(j_int_list, ome)
      end do
    end do
  end subroutine K_using_Gauss_integration_points

  !---------------------------------------------------------------------
  !> @brief Calcula a Integral J discretizando o contorno circular ao redor
  !>        da ponta da trinca em elementos quadráticos (method_key = 3).
  !---------------------------------------------------------------------
  subroutine K_using_Gauss_integration_points_on_boundary()
    type(j_integral_node) :: cp, cp_sym
    type(j_integral_node), allocatable :: j_int_list(:)
    type(boundary_elem) :: last_elem, be
    type(boundary_elem), allocatable :: el_list(:)
    type(boundary_node) :: initial_node, final_node
    type(point_2d) :: gi_pts(3)
    real(dp) :: J_i, J_ii, normal(2), v_normal(2), v_normal_x, v_normal_y, jac
    real(dp) :: grad_i(2), grad_ii(2), grad_i_t(2), grad_ii_t(2)
    real(dp) :: W_i, W_ii, W_i_t, W_ii_t
    real(dp) :: t_i(2), t_ii(2), sig_i_t(3), sig_ii_t(3), sig_i(3), sig_ii(3)
    real(dp) :: eta(20), ome(20), ni(3)
    integer  :: k_ii_sign, i, j, k, l, m, n, n_jpts

    call gauss_legendre_points_1d(np_j, eta(1:np_j), ome(1:np_j))

    do i = 1, cr_num
      do j = 1, cr_list(i)%tips_num
        J_i  = 0.0_dp
        J_ii = 0.0_dp
        call get_initial_node_and_last_element(i, j, last_elem, initial_node, final_node, normal)
        call circular_line_for_j_integral(last_elem%bp_2, initial_node%point, 360.0_dp, 1, &
                                          nbe_j, lb1_j, lb2_j, j_int_list, el_list)

        n_jpts = size(j_int_list)
        !$omp parallel do default(shared) private(k) schedule(dynamic)
        do k = 1, n_jpts
          call evaluate_j_integral_points(j_int_list(k), normal)
        end do
        !$omp end parallel do

        do n = 1, size(el_list)
          be = el_list(n)
          gi_pts(1) = be%bp_1
          gi_pts(2) = j_int_list(be%incid(2))%point
          gi_pts(3) = be%bp_2

          do l = 1, np_j
            call intrp_func_1d(ni, eta(l), be%lamb_1, be%lamb_2)
            call jac_normal_1d(jac, v_normal, gi_pts, eta(l))

            sig_i   = 0.0_dp; sig_ii  = 0.0_dp
            W_i     = 0.0_dp; W_ii    = 0.0_dp
            grad_i  = 0.0_dp; grad_ii = 0.0_dp

            do k = 1, 3
              cp     = j_int_list(be%incid(k))
              cp_sym = j_int_list(sym_node(be%incid(k)))
              call get_variables_key_3(cp, cp_sym, grad_i_t, grad_ii_t, W_i_t, W_ii_t, sig_i_t, sig_ii_t)
              do m = 1, 2
                grad_i(m)  = grad_i(m)  + grad_i_t(m)  * ni(k)
                grad_ii(m) = grad_ii(m) + grad_ii_t(m) * ni(k)
              end do
              do m = 1, 3
                sig_i(m)  = sig_i(m)  + sig_i_t(m)  * ni(k)
                sig_ii(m) = sig_ii(m) + sig_ii_t(m) * ni(k)
              end do
              W_i  = W_i  + W_i_t  * ni(k)
              W_ii = W_ii + W_ii_t * ni(k)
            end do

            v_normal_x = normal(2)*v_normal(1) - normal(1)*v_normal(2)
            v_normal_y = normal(1)*v_normal(1) + normal(2)*v_normal(2)
            t_i(1)  = sig_i(1)*v_normal_x  + sig_i(3)*v_normal_y
            t_i(2)  = sig_i(2)*v_normal_y  + sig_i(3)*v_normal_x
            t_ii(1) = sig_ii(1)*v_normal_x + sig_ii(3)*v_normal_y
            t_ii(2) = sig_ii(2)*v_normal_y + sig_ii(3)*v_normal_x

            J_i  = J_i  + (W_i*v_normal_x  - t_i(1)*grad_i(1)   - t_i(2)*grad_i(2))   * jac * ome(l)
            J_ii = J_ii + (W_ii*v_normal_x - t_ii(1)*grad_ii(1) - t_ii(2)*grad_ii(2)) * jac * ome(l)
          end do
        end do

        call evaluate_sign_of_k_ii(initial_node, final_node, normal, k_ii_sign)
        call k_i_and_k_ii(i, j, k_ii_sign, J_i, J_ii)
        deallocate(j_int_list, el_list)
      end do
    end do
  end subroutine K_using_Gauss_integration_points_on_boundary

  subroutine get_initial_node_and_last_element(i, j, last_elem, initial_node, final_node, normal)
    integer, intent(in)              :: i, j
    type(boundary_elem), intent(out) :: last_elem
    type(boundary_node), intent(out) :: initial_node, final_node
    real(dp), intent(out)            :: normal(2)
    integer :: aux1, initial_node_pos, final_node_pos

    if (j == 1) then
      aux1             = size(cr_list(i)%nodes_num, 1) - initial_node_num + 1
      initial_node_pos = cr_list(i)%nodes_num(aux1, 2)
      final_node_pos   = cr_list(i)%nodes_num(aux1, 1)
      last_elem        = be_list(cr_list(i)%ext_el_num(2))
    else
      initial_node_pos = cr_list(i)%nodes_num(initial_node_num, 1)
      final_node_pos   = cr_list(i)%nodes_num(initial_node_num, 2)
      last_elem        = be_list(cr_list(i)%ext_el_num(4))
    end if

    normal       = get_normal_to_a_boundary_point(initial_node_pos)
    initial_node = bn_list(initial_node_pos)
    final_node   = bn_list(final_node_pos)
  end subroutine get_initial_node_and_last_element

  pure subroutine evaluate_sign_of_k_ii(initial_node, final_node, normal, k_ii_sign)
    type(boundary_node), intent(in) :: initial_node, final_node
    real(dp), intent(in)            :: normal(2)
    integer, intent(out)            :: k_ii_sign
    real(dp) :: u1_sup, u1_inf

    u1_sup = normal(2)*final_node%disp(1)   - normal(1)*final_node%disp(2)
    u1_inf = normal(2)*initial_node%disp(1) - normal(1)*initial_node%disp(2)
    k_ii_sign = merge(0, 1, u1_sup > u1_inf)
  end subroutine evaluate_sign_of_k_ii

  subroutine k_i_and_k_ii(i, j, k_ii_sign, J_i, J_ii)
    integer, intent(in)  :: i, j, k_ii_sign
    real(dp), intent(in) :: J_i, J_ii
    real(dp) :: nu, E, Eb
    integer  :: psk

    call get_plane_state_key(psk)
    call get_poisson_ratio(nu)
    call get_young_modulus(E)
    Eb = merge(E, E / (1.0_dp - nu*nu), psk == 1)

    cr_list(i)%ki(j)  = sqrt(abs(J_i * Eb))
    cr_list(i)%kii(j) = merge(sqrt(abs(J_ii * Eb)), -sqrt(abs(J_ii * Eb)), k_ii_sign == 0)
  end subroutine k_i_and_k_ii

  subroutine evaluate_j_integral_points(cp, normal)
    type(j_integral_node), intent(inout) :: cp
    real(dp), intent(in)                 :: normal(2)
    real(dp) :: sig(3), grad(4), eps(3), hg(4, 6), gg(4, 6), hs(3, 6), gs(3, 6)
    real(dp) :: u_temp(6), t_temp(6), mu, nu, el
    real(dp) :: L(2, 2), sig_rot(3), grad_rot(4), eps_rot(3)
    integer  :: psk, i, j, k
    type(boundary_elem) :: be

    call get_plane_state_key(psk)
    call get_young_modulus(el)
    call get_shear_modulus(mu)
    call get_poisson_ratio(nu)

    sig  = 0.0_dp
    grad = 0.0_dp
    eps  = 0.0_dp

    do j = 1, be_num
      be = be_list(j)
      call reg_integr_j_integr_point_to_be(psk, mu, nu, cp%point, be, hg, gg, hs, gs)
      do k = 1, 3
        u_temp(2*k-1:2*k) = bn_list(be%incid(k))%disp
        t_temp(2*k-1:2*k) = bn_list(be%incid(k))%trac
      end do
      grad = grad - matmul(hg, u_temp) + matmul(gg, t_temp)
      sig  = sig  - matmul(hs, u_temp) + matmul(gs, t_temp)
    end do

    call get_strain_from_stress_elast(psk, nu, el, sig, eps)

    L(1, 1) =  normal(2); L(1, 2) = normal(1)
    L(2, 1) = -normal(1); L(2, 2) = normal(2)
    do i = 1, 2
      sig_rot(i)  = L(1,i)*L(1,i)*sig(1)  + 2.0_dp*L(1,i)*L(2,i)*sig(3) + L(2,i)*L(2,i)*sig(2)
      eps_rot(i)  = L(1,i)*L(1,i)*eps(1)  + 2.0_dp*L(1,i)*L(2,i)*eps(3) + L(2,i)*L(2,i)*eps(2)
      grad_rot(i) = L(1,i)*L(1,i)*grad(1) + L(1,i)*L(2,i)*grad(3) + L(2,i)*L(1,i)*grad(4) + L(2,i)*L(2,i)*grad(2)
    end do
    sig_rot(3)  = L(1,1)*L(1,2)*sig(1)  + L(1,1)*L(2,2)*sig(3)  + L(2,1)*L(1,2)*sig(3)  + L(2,1)*L(2,2)*sig(2)
    eps_rot(3)  = L(1,1)*L(1,2)*eps(1)  + L(1,1)*L(2,2)*eps(3)  + L(2,1)*L(1,2)*eps(3)  + L(2,1)*L(2,2)*eps(2)
    grad_rot(3) = L(1,1)*L(1,2)*grad(1) + L(1,1)*L(2,2)*grad(3) + L(2,1)*L(1,2)*grad(4) + L(2,1)*L(2,2)*grad(2)
    grad_rot(4) = L(1,2)*L(1,1)*grad(1) + L(1,2)*L(2,1)*grad(3) + L(2,2)*L(1,1)*grad(4) + L(2,2)*L(2,1)*grad(2)

    cp%grad(1) = grad_rot(1)
    cp%grad(2) = grad_rot(4)
    cp%stress  = sig_rot
    cp%strain  = eps_rot
  end subroutine evaluate_j_integral_points

  subroutine evaluate_first_and_last_node(cp, node, normal)
    type(j_integral_node), intent(inout) :: cp
    type(boundary_node), intent(in)      :: node
    real(dp), intent(in)                 :: normal(2)
    real(dp) :: sig(3), grad(4), eps(3), hi(4, 6), gi(4, 6), u_temp(6), t_temp(6)
    real(dp) :: mu, nu, el, L(2, 2), sig_rot(3), grad_rot(4), eps_rot(3)
    integer  :: psk, i, j, k, lcn
    logical  :: sing
    type(boundary_elem) :: be

    call get_plane_state_key(psk)
    call get_young_modulus(el)
    call get_shear_modulus(mu)
    call get_poisson_ratio(nu)

    grad = 0.0_dp
    do j = 1, be_num
      be = be_list(j)
      call be_contains_node(be, node, sing, lcn)
      if (sing) then
        call displ_grad_sing_intgr(node%point, be, lcn, hi, gi)
      else
        call reg_integr_grad_displ_to_be_point(psk, mu, nu, node%point, be, hi, gi)
      end if
      do k = 1, 3
        u_temp(2*k-1:2*k) = bn_list(be%incid(k))%disp
        t_temp(2*k-1:2*k) = bn_list(be%incid(k))%trac
      end do
      grad = grad - matmul(hi, u_temp) + matmul(gi, t_temp)
    end do

    sig  = node%stress
    eps  = node%strain
    grad = 2.0_dp * grad

    L(1, 1) =  normal(2); L(1, 2) = normal(1)
    L(2, 1) = -normal(1); L(2, 2) = normal(2)
    do i = 1, 2
      sig_rot(i)  = L(1,i)*L(1,i)*sig(1)  + 2.0_dp*L(1,i)*L(2,i)*sig(3) + L(2,i)*L(2,i)*sig(2)
      eps_rot(i)  = L(1,i)*L(1,i)*eps(1)  + 2.0_dp*L(1,i)*L(2,i)*eps(3) + L(2,i)*L(2,i)*eps(2)
      grad_rot(i) = L(1,i)*L(1,i)*grad(1) + L(1,i)*L(2,i)*grad(3) + L(2,i)*L(1,i)*grad(4) + L(2,i)*L(2,i)*grad(2)
    end do
    sig_rot(3)  = L(1,1)*L(1,2)*sig(1)  + L(1,1)*L(2,2)*sig(3)  + L(2,1)*L(1,2)*sig(3)  + L(2,1)*L(2,2)*sig(2)
    eps_rot(3)  = L(1,1)*L(1,2)*eps(1)  + L(1,1)*L(2,2)*eps(3)  + L(2,1)*L(1,2)*eps(3)  + L(2,1)*L(2,2)*eps(2)
    grad_rot(3) = L(1,1)*L(1,2)*grad(1) + L(1,1)*L(2,2)*grad(3) + L(2,1)*L(1,2)*grad(4) + L(2,1)*L(2,2)*grad(2)
    grad_rot(4) = L(1,2)*L(1,1)*grad(1) + L(1,2)*L(2,1)*grad(3) + L(2,2)*L(1,1)*grad(4) + L(2,2)*L(2,1)*grad(2)

    cp%grad(1) = grad_rot(1)
    cp%grad(2) = grad_rot(4)
    cp%stress  = sig_rot
    cp%strain  = eps_rot
  end subroutine evaluate_first_and_last_node

  pure subroutine get_variables(cp, cp_sym, last_elem, ni, grad_i, grad_ii, W_i_normal, W_ii_normal, t_i, t_ii)
    type(j_integral_node), intent(in) :: cp, cp_sym
    type(boundary_elem), intent(in)   :: last_elem
    real(dp), intent(in)              :: ni(2)
    real(dp), intent(out)             :: grad_i(2), grad_ii(2), W_i_normal, W_ii_normal, t_i(2), t_ii(2)
    real(dp) :: sig_i(3), sig_ii(3), eps_i(3), eps_ii(3), W_i, W_ii, normal(2), r(3), r_x, r_y

    sig_i(1)  = 0.5_dp * (cp%stress(1) + cp_sym%stress(1))
    sig_i(2)  = 0.5_dp * (cp%stress(2) + cp_sym%stress(2))
    sig_i(3)  = 0.5_dp * (cp%stress(3) - cp_sym%stress(3))

    sig_ii(1) = 0.5_dp * (cp%stress(1) - cp_sym%stress(1))
    sig_ii(2) = 0.5_dp * (cp%stress(2) - cp_sym%stress(2))
    sig_ii(3) = 0.5_dp * (cp%stress(3) + cp_sym%stress(3))

    eps_i(1)  = 0.5_dp * (cp%strain(1) + cp_sym%strain(1))
    eps_i(2)  = 0.5_dp * (cp%strain(2) + cp_sym%strain(2))
    eps_i(3)  = 0.5_dp * (cp%strain(3) - cp_sym%strain(3))

    eps_ii(1) = 0.5_dp * (cp%strain(1) - cp_sym%strain(1))
    eps_ii(2) = 0.5_dp * (cp%strain(2) - cp_sym%strain(2))
    eps_ii(3) = 0.5_dp * (cp%strain(3) + cp_sym%strain(3))

    grad_i(1)  = 0.5_dp * (cp%grad(1) + cp_sym%grad(1))
    grad_i(2)  = 0.5_dp * (cp%grad(2) - cp_sym%grad(2))
    grad_ii(1) = 0.5_dp * (cp%grad(1) - cp_sym%grad(1))
    grad_ii(2) = 0.5_dp * (cp%grad(2) + cp_sym%grad(2))

    W_i  = 0.5_dp * (sig_i(1)*eps_i(1)   + sig_i(2)*eps_i(2)   + 2.0_dp*sig_i(3)*eps_i(3))
    W_ii = 0.5_dp * (sig_ii(1)*eps_ii(1) + sig_ii(2)*eps_ii(2) + 2.0_dp*sig_ii(3)*eps_ii(3))

    r          = distance_2d(last_elem%bp_2, cp%point)
    r_x        = ni(2)*r(2) - ni(1)*r(3)
    r_y        = ni(1)*r(2) + ni(2)*r(3)
    normal(1)  = r_x / r(1)
    normal(2)  = r_y / r(1)
    W_i_normal  = W_i  * normal(1)
    W_ii_normal = W_ii * normal(1)

    t_i(1)  = sig_i(1)*normal(1)  + sig_i(3)*normal(2)
    t_i(2)  = sig_i(2)*normal(2)  + sig_i(3)*normal(1)
    t_ii(1) = sig_ii(1)*normal(1) + sig_ii(3)*normal(2)
    t_ii(2) = sig_ii(2)*normal(2) + sig_ii(3)*normal(1)
  end subroutine get_variables

  subroutine get_variables_key_3(cp, cp_sym, grad_i, grad_ii, W_i, W_ii, sig_i, sig_ii)
    type(j_integral_node), intent(in) :: cp, cp_sym
    real(dp), intent(out)             :: grad_i(2), grad_ii(2), W_i, W_ii, sig_i(3), sig_ii(3)
    real(dp) :: eps_i(3), eps_ii(3), el, nu
    integer  :: psk

    sig_i(1)  = 0.5_dp * (cp%stress(1) + cp_sym%stress(1))
    sig_i(2)  = 0.5_dp * (cp%stress(2) + cp_sym%stress(2))
    sig_i(3)  = 0.5_dp * (cp%stress(3) - cp_sym%stress(3))

    sig_ii(1) = 0.5_dp * (cp%stress(1) - cp_sym%stress(1))
    sig_ii(2) = 0.5_dp * (cp%stress(2) - cp_sym%stress(2))
    sig_ii(3) = 0.5_dp * (cp%stress(3) + cp_sym%stress(3))

    grad_i(1)  = 0.5_dp * (cp%grad(1) + cp_sym%grad(1))
    grad_i(2)  = 0.5_dp * (cp%grad(2) - cp_sym%grad(2))
    grad_ii(1) = 0.5_dp * (cp%grad(1) - cp_sym%grad(1))
    grad_ii(2) = 0.5_dp * (cp%grad(2) + cp_sym%grad(2))

    call get_plane_state_key(psk)
    call get_young_modulus(el)
    call get_poisson_ratio(nu)
    call get_strain_from_stress_elast(psk, nu, el, sig_i,  eps_i)
    call get_strain_from_stress_elast(psk, nu, el, sig_ii, eps_ii)

    W_i  = 0.5_dp * (sig_i(1)*eps_i(1)   + sig_i(2)*eps_i(2)   + 2.0_dp*sig_i(3)*eps_i(3))
    W_ii = 0.5_dp * (sig_ii(1)*eps_ii(1) + sig_ii(2)*eps_ii(2) + 2.0_dp*sig_ii(3)*eps_ii(3))
  end subroutine get_variables_key_3

  subroutine Simpson_circular_line(p1, p2, nip, j_int_list, d_theta, r)
    type(point_2d), intent(in)              :: p1, p2
    integer, intent(in)                     :: nip
    type(j_integral_node), intent(inout)    :: j_int_list(:)
    real(dp), intent(out)                   :: d_theta, r
    real(dp) :: rad(3), theta
    integer  :: i, n_seg

    n_seg   = nip - 1
    rad     = distance_2d(p1, p2)
    r       = rad(1)
    theta   = atan2(rad(3), rad(2))
    d_theta = two_pi / real(n_seg, dp)

    if (allocated(sym_node)) deallocate(sym_node)
    allocate(sym_node(nip))

    j_int_list(1)%point  = p2
    j_int_list(1)%stress = 0.0_dp
    j_int_list(1)%grad   = 0.0_dp
    sym_node(1)          = nip

    do i = 2, nip
      theta = theta + d_theta
      j_int_list(i)%point%x = p1%x + r * cos(theta)
      j_int_list(i)%point%y = p1%y + r * sin(theta)
      j_int_list(i)%stress  = 0.0_dp
      j_int_list(i)%grad    = 0.0_dp
      sym_node(i)           = nip + 1 - i
    end do
  end subroutine Simpson_circular_line

  subroutine Gauss_circular_line(p1, p2, nip, j_int_list, r, ome)
    type(point_2d), intent(in)              :: p1, p2
    integer, intent(in)                     :: nip
    type(j_integral_node), intent(inout)    :: j_int_list(:)
    real(dp), intent(out)                   :: r
    real(dp), allocatable, intent(out)      :: ome(:)
    real(dp), allocatable :: eta(:)
    real(dp) :: rad(3), theta, theta_initial
    integer  :: i

    rad           = distance_2d(p1, p2)
    r             = rad(1)
    theta_initial = atan2(rad(3), rad(2))

    allocate(eta(nip), ome(nip))
    call gauss_legendre_points_1d(nip, eta, ome)

    if (allocated(sym_node)) deallocate(sym_node)
    allocate(sym_node(nip))

    do i = 1, nip
      theta = theta_initial + pi * (eta(i) + 1.0_dp)
      j_int_list(i)%point%x = p1%x + r * cos(theta)
      j_int_list(i)%point%y = p1%y + r * sin(theta)
      j_int_list(i)%stress  = 0.0_dp
      j_int_list(i)%grad    = 0.0_dp

      if (mod(nip, 2) == 0) then
        sym_node(i) = merge(i + 1, i - 1, mod(i, 2) == 1)
      else
        if (i == 1) then
          sym_node(i) = 1
        else
          sym_node(i) = merge(i + 1, i - 1, mod(i, 2) == 0)
        end if
      end if
    end do
    deallocate(eta)
  end subroutine Gauss_circular_line

  !---------------------------------------------------------------------
  !> @brief Gera os elementos quadráticos do contorno circular da Integral J.
  !---------------------------------------------------------------------
  subroutine circular_line_for_j_integral(p1, p2, ang, dk, nbe, lb1, lb2, j_int_list, el_list)
    type(point_2d), intent(in)                         :: p1, p2
    real(dp), intent(in)                               :: ang, lb1, lb2
    integer, intent(in)                                :: dk, nbe
    type(j_integral_node), allocatable, intent(out)    :: j_int_list(:)
    type(boundary_elem), allocatable, intent(out)      :: el_list(:)

    type(point_2d), allocatable :: pts(:)
    type(point_2d) :: el_pts(3)
    real(dp)       :: ang_rad
    logical        :: ccw
    integer        :: i, tnjn

    ccw     = (dk >= 0)
    ang_rad = pi * ang / 180.0_dp
    call divide_circular_arc_quad(p1, p2, ang_rad, nbe, ccw, pts)

    allocate(el_list(nbe))
    allocate(j_int_list(2*nbe + 1))

    tnjn = 0
    do i = 1, nbe
      el_pts(1) = pts(2*i - 1)
      el_pts(2) = pts(2*i)
      el_pts(3) = pts(2*i + 1)

      el_list(i)%bp_1      = el_pts(1)
      el_list(i)%bp_2      = el_pts(3)
      el_list(i)%nodes_num = 3
      el_list(i)%lamb_1    = merge(lb1, 1.0_dp, i == 1)
      el_list(i)%lamb_2    = merge(lb2, 1.0_dp, i == nbe)

      if (i == 1) then
        tnjn = tnjn + 1
        j_int_list(tnjn)%point = cartesian_point_from_intrp_func_1d(el_pts, -lb1)
      end if
      el_list(i)%incid(1) = tnjn

      tnjn = tnjn + 1
      el_list(i)%incid(2)    = tnjn
      j_int_list(tnjn)%point = el_pts(2)

      tnjn = tnjn + 1
      el_list(i)%incid(3) = tnjn
      if (i == nbe) then
        j_int_list(tnjn)%point = cartesian_point_from_intrp_func_1d(el_pts, +lb2)
      else
        j_int_list(tnjn)%point = el_pts(3)
      end if
    end do

    if (allocated(sym_node)) deallocate(sym_node)
    allocate(sym_node(tnjn))
    do i = 1, tnjn
      j_int_list(i)%stress = 0.0_dp
      j_int_list(i)%grad   = 0.0_dp
      sym_node(i)          = tnjn - i + 1
    end do
    deallocate(pts)
  end subroutine circular_line_for_j_integral

end module sif_mod
