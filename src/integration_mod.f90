!-----------------------------------------------------------------------
! MODULE: integration_mod
!> @brief Integração numérica de elementos de contorno quadráticos descontínuos:
!>        1. Integração regular e quase-singular (Quadratura de Gauss-Legendre
!>           com subdivisão adaptativa de subelementos);
!>        2. Integração fracamente singular ln(1/r) (Quadratura logarítmica de
!>           Gauss-Laguerre em subelementos à esquerda e à direita do polo);
!>        3. Integração fortemente singular (1/r) e hipersingular (1/r^2) pelo
!>           Método de Expansão Radial de Gao (REIM).
!-----------------------------------------------------------------------
module integration_mod
  use constants_mod
  use quadrature_mod
  use model_mod
  use kelvin_sol_mod
  implicit none

  real(dp), parameter, private :: tol_pw = 1.0e-4_dp
  integer, private             :: npow   = 3  ! Ordem da expansão polinomial radial de Gao

contains

  subroutine set_npow(n)
    integer, intent(in) :: n
    npow = n
  end subroutine set_npow

  subroutine get_npow(n)
    integer, intent(out) :: n
    n = npow
  end subroutine get_npow

  !=====================================================================
  ! 1. INTEGRAÇÃO REGULAR E QUASE-SINGULAR (GAUSS-LEGENDRE COM SUBDIVISÃO)
  !=====================================================================

  !---------------------------------------------------------------------
  !> @brief Define o número de subdivisões `nsd` e pontos de Gauss `nip`
  !>        conforme a proximidade relativa entre o ponto fonte e o elemento.
  !---------------------------------------------------------------------
  pure subroutine near_sing_parameters(cp, be, nsd, nip)
    type(point_2d), intent(in)      :: cp
    type(boundary_elem), intent(in) :: be
    integer, intent(out)            :: nsd, nip
    real(dp) :: d0, dd, le
    integer  :: i

    nip = 20
    dd  = distance_only_2d(cp, bn_list(be%incid(1))%point)
    do i = 2, 3
      d0 = distance_only_2d(cp, bn_list(be%incid(i))%point)
      if (d0 < dd) dd = d0
    end do
    d0 = distance_only_2d(cp, be%bp_1)
    if (d0 < dd) dd = d0
    d0 = distance_only_2d(cp, be%bp_2)
    if (d0 < dd) dd = d0

    le = be_size(be)
    if (dd / le > 2.0_dp) then
      nsd = 1
    else if (dd / le < 0.01_dp) then
      nsd = 200
    else
      nsd = int(2.0_dp * le / dd)
    end if
  end subroutine near_sing_parameters

  !---------------------------------------------------------------------
  !> @brief Integração regular da EIC de Deslocamento (matrizes locais h e g 2x6).
  !---------------------------------------------------------------------
  subroutine displ_bie_reg_integr(psk, mu, nu, cp, be, h, g)
    integer, intent(in)             :: psk
    real(dp), intent(in)            :: mu, nu
    type(point_2d), intent(in)      :: cp
    type(boundary_elem), intent(in) :: be
    real(dp), intent(out)           :: h(2, 6), g(2, 6)

    real(dp)       :: eta_b(20), ome(20), ni(3), seb(2)
    type(point_2d) :: gi_pts(3), ipt
    real(dp)       :: eta, jac, v_normal(2), r_ri(3), w_eff
    real(dp)       :: u(2, 2), t(2, 2)
    integer        :: nsd, nip, i, j, k, kk

    call near_sing_parameters(cp, be, nsd, nip)
    call get_bound_elem_geom_points(be, gi_pts)
    call gauss_legendre_points_1d(nip, eta_b(1:nip), ome(1:nip))

    h = 0.0_dp
    g = 0.0_dp
    seb(2) = -1.0_dp

    do i = 1, nsd
      seb(1) = seb(2)
      if (i == nsd) then
        seb(2) = 1.0_dp
      else
        seb(2) = seb(1) + (2.0_dp / real(nsd, dp))
      end if

      do j = 1, nip
        eta = 0.5_dp * (seb(1) + seb(2)) + (eta_b(j) / real(nsd, dp))
        call intrp_func_1d(ni, eta, be%lamb_1, be%lamb_2)
        call jac_normal_1d(jac, v_normal, gi_pts, eta)
        ipt   = cartesian_point_from_intrp_func_1d(gi_pts, eta)
        r_ri  = distance_2d(cp, ipt)
        u     = u_ij(r_ri, mu, nu, psk)
        t     = t_ij(r_ri, nu, v_normal, psk)
        w_eff = jac * ome(j) / real(nsd, dp)

        do k = 1, 3
          kk = 2*k - 1
          h(:, kk:kk+1) = h(:, kk:kk+1) + (t * w_eff) * ni(k)
          g(:, kk:kk+1) = g(:, kk:kk+1) + (u * w_eff) * ni(k)
        end do
      end do
    end do
  end subroutine displ_bie_reg_integr

  !---------------------------------------------------------------------
  !> @brief Integração regular da EIC de Tração (matrizes locais h e g 2x6).
  !---------------------------------------------------------------------
  subroutine traction_bie_reg_integr(psk, mu, nu, cp, be, normal_cp, h, g)
    integer, intent(in)             :: psk
    real(dp), intent(in)            :: mu, nu, normal_cp(2)
    type(point_2d), intent(in)      :: cp
    type(boundary_elem), intent(in) :: be
    real(dp), intent(out)           :: h(2, 6), g(2, 6)

    real(dp)       :: eta_b(20), ome(20), ni(3), seb(2)
    type(point_2d) :: gi_pts(3), ipt
    real(dp)       :: eta, jac, v_normal(2), r_ri(3), w_eff
    real(dp)       :: u(2, 2), t(2, 2)
    integer        :: nsd, nip, i, j, k, kk

    call near_sing_parameters(cp, be, nsd, nip)
    call get_bound_elem_geom_points(be, gi_pts)
    call gauss_legendre_points_1d(nip, eta_b(1:nip), ome(1:nip))

    h = 0.0_dp
    g = 0.0_dp
    seb(2) = -1.0_dp

    do i = 1, nsd
      seb(1) = seb(2)
      if (i == nsd) then
        seb(2) = 1.0_dp
      else
        seb(2) = seb(1) + (2.0_dp / real(nsd, dp))
      end if

      do j = 1, nip
        eta = 0.5_dp * (seb(1) + seb(2)) + (eta_b(j) / real(nsd, dp))
        call intrp_func_1d(ni, eta, be%lamb_1, be%lamb_2)
        call jac_normal_1d(jac, v_normal, gi_pts, eta)
        ipt   = cartesian_point_from_intrp_func_1d(gi_pts, eta)
        r_ri  = distance_2d(cp, ipt)
        u     = normal_u_ijk(r_ri, normal_cp, nu, psk)
        t     = normal_t_ijk(r_ri, mu, nu, v_normal, normal_cp, psk)
        w_eff = jac * ome(j) / real(nsd, dp)

        do k = 1, 3
          kk = 2*k - 1
          h(:, kk:kk+1) = h(:, kk:kk+1) + (t * w_eff) * ni(k)
          g(:, kk:kk+1) = g(:, kk:kk+1) + (u * w_eff) * ni(k)
        end do
      end do
    end do
  end subroutine traction_bie_reg_integr

  !---------------------------------------------------------------------
  !> @brief Integração regular para tensões no contorno (hs e gs 3x6).
  !---------------------------------------------------------------------
  subroutine stress_reg_integr(psk, mu, nu, cp, be, h, g)
    integer, intent(in)             :: psk
    real(dp), intent(in)            :: mu, nu
    type(point_2d), intent(in)      :: cp
    type(boundary_elem), intent(in) :: be
    real(dp), intent(out)           :: h(3, 6), g(3, 6)

    real(dp)       :: eta_b(20), ome(20), ni(3), seb(2)
    type(point_2d) :: gi_pts(3), ipt
    real(dp)       :: eta, jac, v_normal(2), r_ri(3), w_eff
    real(dp)       :: u(3, 2), t(3, 2)
    integer        :: nsd, nip, i, j, k, kk

    call near_sing_parameters(cp, be, nsd, nip)
    call get_bound_elem_geom_points(be, gi_pts)
    call gauss_legendre_points_1d(nip, eta_b(1:nip), ome(1:nip))

    h = 0.0_dp
    g = 0.0_dp
    seb(2) = -1.0_dp

    do i = 1, nsd
      seb(1) = seb(2)
      if (i == nsd) then
        seb(2) = 1.0_dp
      else
        seb(2) = seb(1) + (2.0_dp / real(nsd, dp))
      end if

      do j = 1, nip
        eta = 0.5_dp * (seb(1) + seb(2)) + (eta_b(j) / real(nsd, dp))
        call intrp_func_1d(ni, eta, be%lamb_1, be%lamb_2)
        call jac_normal_1d(jac, v_normal, gi_pts, eta)
        ipt   = cartesian_point_from_intrp_func_1d(gi_pts, eta)
        r_ri  = distance_2d(cp, ipt)
        u     = u_ijk(r_ri, nu, psk)
        t     = t_ijk(r_ri, mu, nu, v_normal, psk)
        w_eff = jac * ome(j) / real(nsd, dp)

        do k = 1, 3
          kk = 2*k - 1
          h(:, kk:kk+1) = h(:, kk:kk+1) + (t * w_eff) * ni(k)
          g(:, kk:kk+1) = g(:, kk:kk+1) + (u * w_eff) * ni(k)
        end do
      end do
    end do
  end subroutine stress_reg_integr

  !---------------------------------------------------------------------
  !> @brief Integração regular para pontos internos (deslocamentos h,g 2x6
  !>        e tensões hs,gs 3x6).
  !---------------------------------------------------------------------
  subroutine reg_integr_int_point_to_be(psk, mu, nu, cp, be, h, g, hs, gs)
    integer, intent(in)             :: psk
    real(dp), intent(in)            :: mu, nu
    type(point_2d), intent(in)      :: cp
    type(boundary_elem), intent(in) :: be
    real(dp), intent(out)           :: h(2, 6), g(2, 6), hs(3, 6), gs(3, 6)

    real(dp)       :: eta_b(20), ome(20), ni(3), seb(2)
    type(point_2d) :: gi_pts(3), ipt
    real(dp)       :: eta, jac, v_normal(2), r_ri(3), w_eff
    real(dp)       :: u(2, 2), t(2, 2), us(3, 2), ts(3, 2)
    integer        :: nsd, nip, i, j, k, kk

    call near_sing_parameters(cp, be, nsd, nip)
    call get_bound_elem_geom_points(be, gi_pts)
    call gauss_legendre_points_1d(nip, eta_b(1:nip), ome(1:nip))

    h  = 0.0_dp; g  = 0.0_dp
    hs = 0.0_dp; gs = 0.0_dp
    seb(2) = -1.0_dp

    do i = 1, nsd
      seb(1) = seb(2)
      if (i == nsd) then
        seb(2) = 1.0_dp
      else
        seb(2) = seb(1) + (2.0_dp / real(nsd, dp))
      end if

      do j = 1, nip
        eta = 0.5_dp * (seb(1) + seb(2)) + (eta_b(j) / real(nsd, dp))
        call intrp_func_1d(ni, eta, be%lamb_1, be%lamb_2)
        call jac_normal_1d(jac, v_normal, gi_pts, eta)
        ipt   = cartesian_point_from_intrp_func_1d(gi_pts, eta)
        r_ri  = distance_2d(cp, ipt)
        u     = u_ij(r_ri, mu, nu, psk)
        t     = t_ij(r_ri, nu, v_normal, psk)
        us    = u_ijk(r_ri, nu, psk)
        ts    = t_ijk(r_ri, mu, nu, v_normal, psk)
        w_eff = jac * ome(j) / real(nsd, dp)

        do k = 1, 3
          kk = 2*k - 1
          h(:,  kk:kk+1) = h(:,  kk:kk+1) + (t  * w_eff) * ni(k)
          g(:,  kk:kk+1) = g(:,  kk:kk+1) + (u  * w_eff) * ni(k)
          hs(:, kk:kk+1) = hs(:, kk:kk+1) + (ts * w_eff) * ni(k)
          gs(:, kk:kk+1) = gs(:, kk:kk+1) + (us * w_eff) * ni(k)
        end do
      end do
    end do
  end subroutine reg_integr_int_point_to_be

  !---------------------------------------------------------------------
  !> @brief Integração regular para pontos do caminho da Integral J
  !>        (gradiente de deslocamento hg,gg 4x6 e tensões hs,gs 3x6).
  !---------------------------------------------------------------------
  subroutine reg_integr_j_integr_point_to_be(psk, mu, nu, cp, be, hg, gg, hs, gs)
    integer, intent(in)             :: psk
    real(dp), intent(in)            :: mu, nu
    type(point_2d), intent(in)      :: cp
    type(boundary_elem), intent(in) :: be
    real(dp), intent(out)           :: hg(4, 6), gg(4, 6), hs(3, 6), gs(3, 6)

    real(dp)       :: eta_b(20), ome(20), ni(3), seb(2)
    type(point_2d) :: gi_pts(3), ipt
    real(dp)       :: eta, jac, v_normal(2), r_ri(3), w_eff
    real(dp)       :: ug(4, 2), tg(4, 2), us(3, 2), ts(3, 2)
    integer        :: nsd, nip, i, j, k, kk

    call near_sing_parameters(cp, be, nsd, nip)
    call get_bound_elem_geom_points(be, gi_pts)
    call gauss_legendre_points_1d(nip, eta_b(1:nip), ome(1:nip))

    hg = 0.0_dp; gg = 0.0_dp
    hs = 0.0_dp; gs = 0.0_dp
    seb(2) = -1.0_dp

    do i = 1, nsd
      seb(1) = seb(2)
      if (i == nsd) then
        seb(2) = 1.0_dp
      else
        seb(2) = seb(1) + (2.0_dp / real(nsd, dp))
      end if

      do j = 1, nip
        eta = 0.5_dp * (seb(1) + seb(2)) + (eta_b(j) / real(nsd, dp))
        call intrp_func_1d(ni, eta, be%lamb_1, be%lamb_2)
        call jac_normal_1d(jac, v_normal, gi_pts, eta)
        ipt   = cartesian_point_from_intrp_func_1d(gi_pts, eta)
        r_ri  = distance_2d(cp, ipt)
        ug    = u_ij_k(r_ri, mu, nu, psk)
        tg    = t_ij_k(r_ri, nu, v_normal, psk)
        us    = u_ijk(r_ri, nu, psk)
        ts    = t_ijk(r_ri, mu, nu, v_normal, psk)
        w_eff = jac * ome(j) / real(nsd, dp)

        do k = 1, 3
          kk = 2*k - 1
          hg(:, kk:kk+1) = hg(:, kk:kk+1) + (tg * w_eff) * ni(k)
          gg(:, kk:kk+1) = gg(:, kk:kk+1) + (ug * w_eff) * ni(k)
          hs(:, kk:kk+1) = hs(:, kk:kk+1) + (ts * w_eff) * ni(k)
          gs(:, kk:kk+1) = gs(:, kk:kk+1) + (us * w_eff) * ni(k)
        end do
      end do
    end do
  end subroutine reg_integr_j_integr_point_to_be

  !---------------------------------------------------------------------
  !> @brief Integração regular apenas do gradiente de deslocamento (4x6).
  !---------------------------------------------------------------------
  subroutine reg_integr_grad_displ_to_be_point(psk, mu, nu, cp, be, hg, gg)
    integer, intent(in)             :: psk
    real(dp), intent(in)            :: mu, nu
    type(point_2d), intent(in)      :: cp
    type(boundary_elem), intent(in) :: be
    real(dp), intent(out)           :: hg(4, 6), gg(4, 6)

    real(dp)       :: eta_b(20), ome(20), ni(3), seb(2)
    type(point_2d) :: gi_pts(3), ipt
    real(dp)       :: eta, jac, v_normal(2), r_ri(3), w_eff
    real(dp)       :: ug(4, 2), tg(4, 2)
    integer        :: nsd, nip, i, j, k, kk

    call near_sing_parameters(cp, be, nsd, nip)
    call get_bound_elem_geom_points(be, gi_pts)
    call gauss_legendre_points_1d(nip, eta_b(1:nip), ome(1:nip))

    hg = 0.0_dp; gg = 0.0_dp
    seb(2) = -1.0_dp

    do i = 1, nsd
      seb(1) = seb(2)
      if (i == nsd) then
        seb(2) = 1.0_dp
      else
        seb(2) = seb(1) + (2.0_dp / real(nsd, dp))
      end if

      do j = 1, nip
        eta = 0.5_dp * (seb(1) + seb(2)) + (eta_b(j) / real(nsd, dp))
        call intrp_func_1d(ni, eta, be%lamb_1, be%lamb_2)
        call jac_normal_1d(jac, v_normal, gi_pts, eta)
        ipt   = cartesian_point_from_intrp_func_1d(gi_pts, eta)
        r_ri  = distance_2d(cp, ipt)
        ug    = u_ij_k(r_ri, mu, nu, psk)
        tg    = t_ij_k(r_ri, nu, v_normal, psk)
        w_eff = jac * ome(j) / real(nsd, dp)

        do k = 1, 3
          kk = 2*k - 1
          hg(:, kk:kk+1) = hg(:, kk:kk+1) + (tg * w_eff) * ni(k)
          gg(:, kk:kk+1) = gg(:, kk:kk+1) + (ug * w_eff) * ni(k)
        end do
      end do
    end do
  end subroutine reg_integr_grad_displ_to_be_point

  !=====================================================================
  ! 2. INTEGRAÇÃO FRACAMENTE SINGULAR (QUADRATURA LOGARÍTMICA)
  !=====================================================================

  !---------------------------------------------------------------------
  !> @brief Integra o núcleo fracamente singular u*_ij sobre o elemento
  !>        quadrático descontínuo `be` quando o ponto fonte é o nó local `lcp`.
  !---------------------------------------------------------------------
  subroutine displ_bie_weakly_sing_integr(psk, mu, nu, lcp, be, g)
    integer, intent(in)             :: psk, lcp
    real(dp), intent(in)            :: mu, nu
    type(boundary_elem), intent(in) :: be
    real(dp), intent(out)           :: g(2, 6)

    integer, parameter :: nip  = 20  ! Pontos de Gauss-Legendre (parte regular)
    integer, parameter :: nipl = 12  ! Pontos de Gauss-Laguerre (parte logarítmica)
    real(dp)       :: eta_g(nip), ome_g(nip), eta_l(nipl), ome_l(nipl)
    real(dp)       :: ni(3), intgr(3), r(3)
    real(dp)       :: nu_b, jac, eta0, eta_b, cte, f1, f2
    type(point_2d) :: gi_pts(3), se_pts(3), cp_pt, it_pt
    integer        :: i, j

    nu_b = merge(nu / (1.0_dp + nu), nu, psk == 1)
    g    = 0.0_dp

    select case (lcp)
    case (1)
      eta0 = -be%lamb_1
    case (2)
      eta0 = 0.0_dp
    case (3)
      eta0 = +be%lamb_2
    case default
      call error_message_int_stop("Indice local invalido em displ_bie_weakly_sing_integr: ", lcp)
    end select

    call get_bound_elem_geom_points(be, gi_pts)
    cp_pt = cartesian_point_from_intrp_func_1d(gi_pts, eta0)

    call gauss_legendre_points_1d(nip, eta_g, ome_g)
    call gauss_laguerre_points_1d(nipl, eta_l, ome_l)

    ! --- Subelemento 1: de eta = -1 até eta0 ---
    if (eta0 > -1.0_dp) then
      eta_b     = 0.5_dp * (eta0 - 1.0_dp)
      se_pts(1) = be%bp_1
      se_pts(2) = cartesian_point_from_intrp_func_1d(gi_pts, eta_b)
      se_pts(3) = cp_pt

      do i = 1, nip
        eta_b = 0.5_dp * (eta0 + 1.0_dp) * eta_g(i) + 0.5_dp * (eta0 - 1.0_dp)
        it_pt = cartesian_point_from_intrp_func_1d(gi_pts, eta_b)
        r     = distance_2d(cp_pt, it_pt)
        call intrp_func_1d(ni, eta_b, be%lamb_1, be%lamb_2)
        call jac_only_1d(jac, gi_pts, eta_b)
        call evaluate_f1_f2(3, se_pts, eta_g(i), f1, f2)

        ni = ni * 0.5_dp * (1.0_dp + eta0) * jac * ome_g(i)
        do j = 1, 2
          cte = 0.5_dp * (r(j+1)*r(j+1)/(r(1)*r(1)) - 0.5_dp*(3.0_dp - 4.0_dp*nu_b)*log(f1*f1 + f2*f2))
          g(j, j)   = g(j, j)   + cte * ni(1)
          g(j, j+2) = g(j, j+2) + cte * ni(2)
          g(j, j+4) = g(j, j+4) + cte * ni(3)
        end do
        cte     = 0.5_dp * r(2) * r(3) / (r(1) * r(1))
        g(1, 2) = g(1, 2) + cte * ni(1)
        g(1, 4) = g(1, 4) + cte * ni(2)
        g(1, 6) = g(1, 6) + cte * ni(3)
      end do

      intgr = 0.0_dp
      do i = 1, nipl
        eta_b = 0.5_dp * (eta0 + 1.0_dp) * (1.0_dp - 2.0_dp * eta_l(i)) + 0.5_dp * (eta0 - 1.0_dp)
        call intrp_func_1d(ni, eta_b, be%lamb_1, be%lamb_2)
        call jac_only_1d(jac, gi_pts, eta_b)
        intgr = intgr + ni * jac * ome_l(i)
      end do
      intgr = intgr * 0.5_dp * (1.0_dp + eta0)
      cte   = 3.0_dp - 4.0_dp * nu_b
      do i = 1, 2
        g(i, i)   = g(i, i)   + cte * intgr(1)
        g(i, i+2) = g(i, i+2) + cte * intgr(2)
        g(i, i+4) = g(i, i+4) + cte * intgr(3)
      end do
    end if

    ! --- Subelemento 2: de eta0 até eta = +1 ---
    if (eta0 < 1.0_dp) then
      eta_b     = 0.5_dp * (eta0 + 1.0_dp)
      se_pts(1) = cp_pt
      se_pts(2) = cartesian_point_from_intrp_func_1d(gi_pts, eta_b)
      se_pts(3) = be%bp_2

      do i = 1, nip
        eta_b = 0.5_dp * (1.0_dp - eta0) * eta_g(i) + 0.5_dp * (1.0_dp + eta0)
        it_pt = cartesian_point_from_intrp_func_1d(gi_pts, eta_b)
        r     = distance_2d(cp_pt, it_pt)
        call intrp_func_1d(ni, eta_b, be%lamb_1, be%lamb_2)
        call jac_only_1d(jac, gi_pts, eta_b)
        call evaluate_f1_f2(1, se_pts, eta_g(i), f1, f2)

        ni = ni * 0.5_dp * (1.0_dp - eta0) * jac * ome_g(i)
        do j = 1, 2
          cte = 0.5_dp * (r(j+1)*r(j+1)/(r(1)*r(1)) - 0.5_dp*(3.0_dp - 4.0_dp*nu_b)*log(f1*f1 + f2*f2))
          g(j, j)   = g(j, j)   + cte * ni(1)
          g(j, j+2) = g(j, j+2) + cte * ni(2)
          g(j, j+4) = g(j, j+4) + cte * ni(3)
        end do
        cte     = 0.5_dp * r(2) * r(3) / (r(1) * r(1))
        g(1, 2) = g(1, 2) + cte * ni(1)
        g(1, 4) = g(1, 4) + cte * ni(2)
        g(1, 6) = g(1, 6) + cte * ni(3)
      end do

      intgr = 0.0_dp
      do i = 1, nipl
        eta_b = 0.5_dp * (1.0_dp - eta0) * (2.0_dp * eta_l(i) - 1.0_dp) + 0.5_dp * (1.0_dp + eta0)
        call intrp_func_1d(ni, eta_b, be%lamb_1, be%lamb_2)
        call jac_only_1d(jac, gi_pts, eta_b)
        intgr = intgr + ni * jac * ome_l(i)
      end do
      intgr = intgr * 0.5_dp * (1.0_dp - eta0)
      cte   = 3.0_dp - 4.0_dp * nu_b
      do i = 1, 2
        g(i, i)   = g(i, i)   + cte * intgr(1)
        g(i, i+2) = g(i, i+2) + cte * intgr(2)
        g(i, i+4) = g(i, i+4) + cte * intgr(3)
      end do
    end if

    g(2, 1) = g(1, 2)
    g(2, 3) = g(1, 4)
    g(2, 5) = g(1, 6)
    g = g / (four_pi * mu * (1.0_dp - nu_b))
  end subroutine displ_bie_weakly_sing_integr

  ! Alias com o nome usado pelo montador
  subroutine evaluate_f1_f2(lcp, pts, eta, f1, f2)
    integer, intent(in)        :: lcp
    type(point_2d), intent(in) :: pts(3)
    real(dp), intent(in)       :: eta
    real(dp), intent(out)      :: f1, f2

    select case (lcp)
    case (1)
      f1 = -(2.0_dp - eta)*pts(1)%x + eta*pts(3)%x + 2.0_dp*(1.0_dp - eta)*pts(2)%x
      f2 = -(2.0_dp - eta)*pts(1)%y + eta*pts(3)%y + 2.0_dp*(1.0_dp - eta)*pts(2)%y
    case (2)
      f1 = 0.5_dp*((eta - 1.0_dp)*pts(1)%x + (eta + 1.0_dp)*pts(3)%x) - eta*pts(2)%x
      f2 = 0.5_dp*((eta - 1.0_dp)*pts(1)%y + (eta + 1.0_dp)*pts(3)%y) - eta*pts(2)%y
    case (3)
      f1 = -(2.0_dp + eta)*pts(3)%x - eta*pts(1)%x + 2.0_dp*(1.0_dp + eta)*pts(2)%x
      f2 = -(2.0_dp + eta)*pts(3)%y - eta*pts(1)%y + 2.0_dp*(1.0_dp + eta)*pts(2)%y
    case default
      call error_message_int_stop("Indice local invalido em evaluate_f1_f2: ", lcp)
    end select
  end subroutine evaluate_f1_f2

  !=====================================================================
  ! 3. INTEGRAÇÃO FORTE E HIPERSINGULAR (MÉTODO DE EXPANSÃO RADIAL DE GAO)
  !=====================================================================

  !---------------------------------------------------------------------
  !> @brief Coordenada intrínseca do ponto fonte local `lcp` no elemento `be`.
  !---------------------------------------------------------------------
  pure function get_source_eta0(be, lcp) result(eta0)
    type(boundary_elem), intent(in) :: be
    integer, intent(in)             :: lcp
    real(dp) :: eta0

    if (lcp == 1) then
      eta0 = -be%lamb_1
    else if (lcp == 3) then
      eta0 = +be%lamb_2
    else
      eta0 = 0.0_dp
    end if
  end function get_source_eta0

  !---------------------------------------------------------------------
  !> @brief Integração fortemente singular de t*_ij na EIC de Deslocamento.
  !---------------------------------------------------------------------
  subroutine displ_bie_sing_integr(sp, be, lcp, hi)
    type(point_2d), intent(in)      :: sp
    type(boundary_elem), intent(in) :: be
    integer, intent(in)             :: lcp
    real(dp), intent(inout)         :: hi(2, 6)
    real(dp) :: eta0

    hi   = 0.0_dp
    eta0 = get_source_eta0(be, lcp)
    if (eta0 > -1.0_dp) call displ_subreg_sing_intgr(sp, be, eta0, -1.0_dp, hi)
    if (eta0 <  1.0_dp) call displ_subreg_sing_intgr(sp, be, eta0,  1.0_dp, hi)
  end subroutine displ_bie_sing_integr

  subroutine displ_subreg_sing_intgr(sp, be, eta0, eta_final, hi)
    type(point_2d), intent(in)      :: sp
    type(boundary_elem), intent(in) :: be
    real(dp), intent(in)            :: eta0, eta_final
    real(dp), intent(inout)         :: hi(2, 6)

    type(point_2d) :: fp, begp(3)
    real(dp)       :: d_eta, normal(2), normal_eta0(2), ti(2), ri(3), drdx(2), r_dot_t
    real(dp)       :: c0(2, 6), ci(npow, 2, 6), func(2, 6), ni(3), amat(npow, npow)
    real(dp)       :: eta, pw, ts(2, 2), tk(2, 2), nu
    integer        :: i, j, k, psk

    d_eta = (eta_final - eta0) / real(npow, dp)
    call get_bound_elem_geom_points(be, begp)
    call get_plane_state_key(psk)
    call get_poisson_ratio(nu)

    normal_eta0 = normal_1d(begp, eta0)
    drdx(1) = -normal_eta0(2)
    drdx(2) =  normal_eta0(1)
    if (eta0 > eta_final) drdx = -drdx
    tk = t_ij_rp(drdx, normal_eta0, nu, psk)
    call intrp_func_1d(ni, eta0, be%lamb_1, be%lamb_2)
    do k = 1, 3
      c0(1:2, 2*k-1:2*k) = ni(k) * tk
    end do

    eta = eta0
    do i = 1, npow
      eta = eta + d_eta
      if (i == npow) eta = eta_final
      normal = normal_1d(begp, eta)
      ti(1) = -normal(2)
      ti(2) =  normal(1)
      if (eta0 > eta_final) ti = -ti
      fp   = cartesian_point_from_intrp_func_1d(begp, eta)
      ri   = distance_2d(sp, fp)
      drdx = ri(2:3) / ri(1)
      ts   = t_ij_rp(drdx, normal, nu, psk)
      call intrp_func_1d(ni, eta, be%lamb_1, be%lamb_2)
      do k = 1, 3
        func(1:2, 2*k-1:2*k) = ni(k) * ts
      end do
      r_dot_t     = dot_product(drdx, ti)
      ci(i, :, :) = ((func / abs(r_dot_t)) - c0) / ri(1)
      amat(i, 1)  = 1.0_dp
      do j = 2, npow
        amat(i, j) = amat(i, j-1) * ri(1)
      end do
    end do

    call mat_inverse(amat)
    do i = 1, 2
      do j = 1, 6
        ci(:, i, j) = matmul(amat, ci(:, i, j))
      end do
    end do

    hi = hi + c0 * log(ri(1))
    do i = 1, npow
      pw = real(i, dp)
      hi = hi + ci(i, :, :) * (ri(1)**pw) / pw
    end do
  end subroutine displ_subreg_sing_intgr

  !---------------------------------------------------------------------
  !> @brief Integração hipersingular (hi) e fortemente singular (gi) na
  !>        EIC de Tração (face da trinca com bie = 3).
  !---------------------------------------------------------------------
  subroutine traction_bie_sing_intgr(sp, be, lcp, normal_cp, hi, gi)
    type(point_2d), intent(in)      :: sp
    type(boundary_elem), intent(in) :: be
    integer, intent(in)             :: lcp
    real(dp), intent(in)            :: normal_cp(2)
    real(dp), intent(inout)         :: hi(2, 6), gi(2, 6)
    real(dp) :: eta0

    hi   = 0.0_dp
    gi   = 0.0_dp
    eta0 = get_source_eta0(be, lcp)
    if (eta0 > -1.0_dp) then
      call t_subreg_sing_intgr(sp, be, eta0, -1.0_dp, normal_cp, hi)
      call u_subreg_sing_intgr(sp, be, eta0, -1.0_dp, normal_cp, gi)
    end if
    if (eta0 < 1.0_dp) then
      call t_subreg_sing_intgr(sp, be, eta0,  1.0_dp, normal_cp, hi)
      call u_subreg_sing_intgr(sp, be, eta0,  1.0_dp, normal_cp, gi)
    end if
  end subroutine traction_bie_sing_intgr

  subroutine t_subreg_sing_intgr(sp, be, eta0, eta_final, normal_cp, hi)
    type(point_2d), intent(in)      :: sp
    type(boundary_elem), intent(in) :: be
    real(dp), intent(in)            :: eta0, eta_final, normal_cp(2)
    real(dp), intent(inout)         :: hi(2, 6)

    type(point_2d) :: fp, begp(3)
    real(dp)       :: d_eta, normal(2), normal_eta0(2), ti(2), ri(3), drdx(2), r_dot_t
    real(dp)       :: c0(2, 6), ci(npow, 2, 6), func(2, 6), ni(3), amat(npow, npow)
    real(dp)       :: eta, pw, ts(2, 2), tk(2, 2), mu, nu
    integer        :: i, j, k, psk

    d_eta = (eta_final - eta0) / real(npow, dp)
    call get_bound_elem_geom_points(be, begp)
    call get_plane_state_key(psk)
    call get_poisson_ratio(nu)
    call get_shear_modulus(mu)

    normal_eta0 = normal_1d(begp, eta0)
    drdx(1) = -normal_eta0(2)
    drdx(2) =  normal_eta0(1)
    if (eta0 > eta_final) drdx = -drdx
    tk = normal_t_ijk_rp(drdx, normal_eta0, normal_cp, mu, nu, psk)
    call intrp_func_1d(ni, eta0, be%lamb_1, be%lamb_2)
    do k = 1, 3
      c0(1:2, 2*k-1:2*k) = ni(k) * tk
    end do

    eta = eta0
    do i = 1, npow
      eta = eta + d_eta
      if (i == npow) eta = eta_final
      normal = normal_1d(begp, eta)
      ti(1) = -normal(2)
      ti(2) =  normal(1)
      if (eta0 > eta_final) ti = -ti
      fp   = cartesian_point_from_intrp_func_1d(begp, eta)
      ri   = distance_2d(sp, fp)
      drdx = ri(2:3) / ri(1)
      ts   = normal_t_ijk_rp(drdx, normal, normal_cp, mu, nu, psk)
      call intrp_func_1d(ni, eta, be%lamb_1, be%lamb_2)
      do k = 1, 3
        func(1:2, 2*k-1:2*k) = ni(k) * ts
      end do
      r_dot_t     = dot_product(drdx, ti)
      ci(i, :, :) = ((func / abs(r_dot_t)) - c0) / ri(1)
      amat(i, 1)  = 1.0_dp
      do j = 2, npow
        amat(i, j) = amat(i, j-1) * ri(1)
      end do
    end do

    call mat_inverse(amat)
    do i = 1, 2
      do j = 1, 6
        ci(:, i, j) = matmul(amat, ci(:, i, j))
      end do
    end do

    pw = -1.0_dp
    hi = hi + c0 * (ri(1)**pw) / pw
    do i = 1, npow
      pw = -1.0_dp + real(i, dp)
      if (abs(pw) < tol_pw) then
        hi = hi + ci(i, :, :) * log(ri(1))
      else
        hi = hi + ci(i, :, :) * (ri(1)**pw) / pw
      end if
    end do
  end subroutine t_subreg_sing_intgr

  subroutine u_subreg_sing_intgr(sp, be, eta0, eta_final, normal_cp, gi)
    type(point_2d), intent(in)      :: sp
    type(boundary_elem), intent(in) :: be
    real(dp), intent(in)            :: eta0, eta_final, normal_cp(2)
    real(dp), intent(inout)         :: gi(2, 6)

    type(point_2d) :: fp, begp(3)
    real(dp)       :: d_eta, normal(2), normal_eta0(2), ti(2), ri(3), drdx(2), r_dot_t
    real(dp)       :: c0(2, 6), ci(npow, 2, 6), func(2, 6), ni(3), amat(npow, npow)
    real(dp)       :: eta, pw, ui(2, 2), us(2, 2), nu
    integer        :: i, j, k, psk

    d_eta = (eta_final - eta0) / real(npow, dp)
    call get_bound_elem_geom_points(be, begp)
    call get_plane_state_key(psk)
    call get_poisson_ratio(nu)

    normal_eta0 = normal_1d(begp, eta0)
    drdx(1) = -normal_eta0(2)
    drdx(2) =  normal_eta0(1)
    if (eta0 > eta_final) drdx = -drdx
    us = normal_u_ijk_rp(drdx, normal_cp, nu, psk)
    call intrp_func_1d(ni, eta0, be%lamb_1, be%lamb_2)
    do k = 1, 3
      c0(1:2, 2*k-1:2*k) = ni(k) * us
    end do

    eta = eta0
    do i = 1, npow
      eta = eta + d_eta
      if (i == npow) eta = eta_final
      normal = normal_1d(begp, eta)
      ti(1) = -normal(2)
      ti(2) =  normal(1)
      if (eta0 > eta_final) ti = -ti
      fp   = cartesian_point_from_intrp_func_1d(begp, eta)
      ri   = distance_2d(sp, fp)
      drdx = ri(2:3) / ri(1)
      ui   = normal_u_ijk_rp(drdx, normal_cp, nu, psk)
      call intrp_func_1d(ni, eta, be%lamb_1, be%lamb_2)
      do k = 1, 3
        func(1:2, 2*k-1:2*k) = ni(k) * ui
      end do
      r_dot_t     = dot_product(drdx, ti)
      ci(i, :, :) = ((func / abs(r_dot_t)) - c0) / ri(1)
      amat(i, 1)  = 1.0_dp
      do j = 2, npow
        amat(i, j) = amat(i, j-1) * ri(1)
      end do
    end do

    call mat_inverse(amat)
    do i = 1, 2
      do j = 1, 6
        ci(:, i, j) = matmul(amat, ci(:, i, j))
      end do
    end do

    gi = gi + c0 * log(ri(1))
    do i = 1, npow
      pw = real(i, dp)
      gi = gi + ci(i, :, :) * (ri(1)**pw) / pw
    end do
  end subroutine u_subreg_sing_intgr

  !---------------------------------------------------------------------
  !> @brief Integração singular para o gradiente de deslocamento (4x6).
  !---------------------------------------------------------------------
  subroutine displ_grad_sing_intgr(sp, be, lcp, hi, gi)
    type(point_2d), intent(in)      :: sp
    type(boundary_elem), intent(in) :: be
    integer, intent(in)             :: lcp
    real(dp), intent(inout)         :: hi(4, 6), gi(4, 6)
    real(dp) :: eta0

    hi   = 0.0_dp
    gi   = 0.0_dp
    eta0 = get_source_eta0(be, lcp)
    if (eta0 > -1.0_dp) then
      call grad_t_subreg_sing_intgr(sp, be, eta0, -1.0_dp, hi)
      call grad_u_subreg_sing_intgr(sp, be, eta0, -1.0_dp, gi)
    end if
    if (eta0 < 1.0_dp) then
      call grad_t_subreg_sing_intgr(sp, be, eta0,  1.0_dp, hi)
      call grad_u_subreg_sing_intgr(sp, be, eta0,  1.0_dp, gi)
    end if
  end subroutine displ_grad_sing_intgr

  subroutine grad_t_subreg_sing_intgr(sp, be, eta0, eta_final, hi)
    type(point_2d), intent(in)      :: sp
    type(boundary_elem), intent(in) :: be
    real(dp), intent(in)            :: eta0, eta_final
    real(dp), intent(inout)         :: hi(4, 6)

    type(point_2d) :: fp, begp(3)
    real(dp)       :: d_eta, normal(2), normal_eta0(2), ti(2), ri(3), drdx(2), r_dot_t
    real(dp)       :: c0(4, 6), ci(npow, 4, 6), func(4, 6), ni(3), amat(npow, npow)
    real(dp)       :: eta, pw, ts(4, 2), tk(4, 2), nu
    integer        :: i, j, k, psk

    d_eta = (eta_final - eta0) / real(npow, dp)
    call get_bound_elem_geom_points(be, begp)
    call get_plane_state_key(psk)
    call get_poisson_ratio(nu)

    normal_eta0 = normal_1d(begp, eta0)
    drdx(1) = -normal_eta0(2)
    drdx(2) =  normal_eta0(1)
    if (eta0 > eta_final) drdx = -drdx
    tk = t_ij_k_rp(drdx, nu, normal_eta0, psk)
    call intrp_func_1d(ni, eta0, be%lamb_1, be%lamb_2)
    do k = 1, 3
      c0(1:4, 2*k-1:2*k) = ni(k) * tk
    end do

    eta = eta0
    do i = 1, npow
      eta = eta + d_eta
      if (i == npow) eta = eta_final
      normal = normal_1d(begp, eta)
      ti(1) = -normal(2)
      ti(2) =  normal(1)
      if (eta0 > eta_final) ti = -ti
      fp   = cartesian_point_from_intrp_func_1d(begp, eta)
      ri   = distance_2d(sp, fp)
      drdx = ri(2:3) / ri(1)
      ts   = t_ij_k_rp(drdx, nu, normal, psk)
      call intrp_func_1d(ni, eta, be%lamb_1, be%lamb_2)
      do k = 1, 3
        func(1:4, 2*k-1:2*k) = ni(k) * ts
      end do
      r_dot_t     = dot_product(drdx, ti)
      ci(i, :, :) = ((func / abs(r_dot_t)) - c0) / ri(1)
      amat(i, 1)  = 1.0_dp
      do j = 2, npow
        amat(i, j) = amat(i, j-1) * ri(1)
      end do
    end do

    call mat_inverse(amat)
    do i = 1, 4
      do j = 1, 6
        ci(:, i, j) = matmul(amat, ci(:, i, j))
      end do
    end do

    pw = -1.0_dp
    hi = hi + c0 * (ri(1)**pw) / pw
    do i = 1, npow
      pw = -1.0_dp + real(i, dp)
      if (abs(pw) < tol_pw) then
        hi = hi + ci(i, :, :) * log(ri(1))
      else
        hi = hi + ci(i, :, :) * (ri(1)**pw) / pw
      end if
    end do
  end subroutine grad_t_subreg_sing_intgr

  subroutine grad_u_subreg_sing_intgr(sp, be, eta0, eta_final, gi)
    type(point_2d), intent(in)      :: sp
    type(boundary_elem), intent(in) :: be
    real(dp), intent(in)            :: eta0, eta_final
    real(dp), intent(inout)         :: gi(4, 6)

    type(point_2d) :: fp, begp(3)
    real(dp)       :: d_eta, normal(2), normal_eta0(2), ti(2), ri(3), drdx(2), r_dot_t
    real(dp)       :: c0(4, 6), ci(npow, 4, 6), func(4, 6), ni(3), amat(npow, npow)
    real(dp)       :: eta, pw, ui(4, 2), us(4, 2), mu, nu
    integer        :: i, j, k, psk

    d_eta = (eta_final - eta0) / real(npow, dp)
    call get_bound_elem_geom_points(be, begp)
    call get_plane_state_key(psk)
    call get_poisson_ratio(nu)
    call get_shear_modulus(mu)

    normal_eta0 = normal_1d(begp, eta0)
    drdx(1) = -normal_eta0(2)
    drdx(2) =  normal_eta0(1)
    if (eta0 > eta_final) drdx = -drdx
    us = u_ij_k_rp(drdx, mu, nu, psk)
    call intrp_func_1d(ni, eta0, be%lamb_1, be%lamb_2)
    do k = 1, 3
      c0(1:4, 2*k-1:2*k) = ni(k) * us
    end do

    eta = eta0
    do i = 1, npow
      eta = eta + d_eta
      if (i == npow) eta = eta_final
      normal = normal_1d(begp, eta)
      ti(1) = -normal(2)
      ti(2) =  normal(1)
      if (eta0 > eta_final) ti = -ti
      fp   = cartesian_point_from_intrp_func_1d(begp, eta)
      ri   = distance_2d(sp, fp)
      drdx = ri(2:3) / ri(1)
      ui   = u_ij_k_rp(drdx, mu, nu, psk)
      call intrp_func_1d(ni, eta, be%lamb_1, be%lamb_2)
      do k = 1, 3
        func(1:4, 2*k-1:2*k) = ni(k) * ui
      end do
      r_dot_t     = dot_product(drdx, ti)
      ci(i, :, :) = ((func / abs(r_dot_t)) - c0) / ri(1)
      amat(i, 1)  = 1.0_dp
      do j = 2, npow
        amat(i, j) = amat(i, j-1) * ri(1)
      end do
    end do

    call mat_inverse(amat)
    do i = 1, 4
      do j = 1, 6
        ci(:, i, j) = matmul(amat, ci(:, i, j))
      end do
    end do

    gi = gi + c0 * log(ri(1))
    do i = 1, npow
      pw = real(i, dp)
      gi = gi + ci(i, :, :) * (ri(1)**pw) / pw
    end do
  end subroutine grad_u_subreg_sing_intgr

  !---------------------------------------------------------------------
  !> @brief Integração singular para tensões no contorno (3x6).
  !---------------------------------------------------------------------
  subroutine stress_sing_intgr(sp, be, lcp, hi, gi)
    type(point_2d), intent(in)      :: sp
    type(boundary_elem), intent(in) :: be
    integer, intent(in)             :: lcp
    real(dp), intent(inout)         :: hi(3, 6), gi(3, 6)
    real(dp) :: eta0

    hi   = 0.0_dp
    gi   = 0.0_dp
    eta0 = get_source_eta0(be, lcp)
    if (eta0 > -1.0_dp) then
      call stress_t_subreg_sing_intgr(sp, be, eta0, -1.0_dp, hi)
      call stress_u_subreg_sing_intgr(sp, be, eta0, -1.0_dp, gi)
    end if
    if (eta0 < 1.0_dp) then
      call stress_t_subreg_sing_intgr(sp, be, eta0,  1.0_dp, hi)
      call stress_u_subreg_sing_intgr(sp, be, eta0,  1.0_dp, gi)
    end if
  end subroutine stress_sing_intgr

  subroutine stress_t_subreg_sing_intgr(sp, be, eta0, eta_final, hi)
    type(point_2d), intent(in)      :: sp
    type(boundary_elem), intent(in) :: be
    real(dp), intent(in)            :: eta0, eta_final
    real(dp), intent(inout)         :: hi(3, 6)

    type(point_2d) :: fp, begp(3)
    real(dp)       :: d_eta, normal(2), normal_eta0(2), ti(2), ri(3), drdx(2), r_dot_t
    real(dp)       :: c0(3, 6), ci(npow, 3, 6), func(3, 6), ni(3), amat(npow, npow)
    real(dp)       :: eta, pw, ts(3, 2), tk(3, 2), mu, nu
    integer        :: i, j, k, psk

    d_eta = (eta_final - eta0) / real(npow, dp)
    call get_bound_elem_geom_points(be, begp)
    call get_plane_state_key(psk)
    call get_poisson_ratio(nu)
    call get_shear_modulus(mu)

    normal_eta0 = normal_1d(begp, eta0)
    drdx(1) = -normal_eta0(2)
    drdx(2) =  normal_eta0(1)
    if (eta0 > eta_final) drdx = -drdx
    tk = t_ijk_rp(drdx, normal_eta0, mu, nu, psk)
    call intrp_func_1d(ni, eta0, be%lamb_1, be%lamb_2)
    do k = 1, 3
      c0(1:3, 2*k-1:2*k) = ni(k) * tk
    end do

    eta = eta0
    do i = 1, npow
      eta = eta + d_eta
      if (i == npow) eta = eta_final
      normal = normal_1d(begp, eta)
      ti(1) = -normal(2)
      ti(2) =  normal(1)
      if (eta0 > eta_final) ti = -ti
      fp   = cartesian_point_from_intrp_func_1d(begp, eta)
      ri   = distance_2d(sp, fp)
      drdx = ri(2:3) / ri(1)
      ts   = t_ijk_rp(drdx, normal, mu, nu, psk)
      call intrp_func_1d(ni, eta, be%lamb_1, be%lamb_2)
      do k = 1, 3
        func(1:3, 2*k-1:2*k) = ni(k) * ts
      end do
      r_dot_t     = dot_product(drdx, ti)
      ci(i, :, :) = ((func / abs(r_dot_t)) - c0) / ri(1)
      amat(i, 1)  = 1.0_dp
      do j = 2, npow
        amat(i, j) = amat(i, j-1) * ri(1)
      end do
    end do

    call mat_inverse(amat)
    do i = 1, 3
      do j = 1, 6
        ci(:, i, j) = matmul(amat, ci(:, i, j))
      end do
    end do

    pw = -1.0_dp
    hi = hi + c0 * (ri(1)**pw) / pw
    do i = 1, npow
      pw = -1.0_dp + real(i, dp)
      if (abs(pw) < tol_pw) then
        hi = hi + ci(i, :, :) * log(ri(1))
      else
        hi = hi + ci(i, :, :) * (ri(1)**pw) / pw
      end if
    end do
  end subroutine stress_t_subreg_sing_intgr

  subroutine stress_u_subreg_sing_intgr(sp, be, eta0, eta_final, gi)
    type(point_2d), intent(in)      :: sp
    type(boundary_elem), intent(in) :: be
    real(dp), intent(in)            :: eta0, eta_final
    real(dp), intent(inout)         :: gi(3, 6)

    type(point_2d) :: fp, begp(3)
    real(dp)       :: d_eta, normal(2), normal_eta0(2), ti(2), ri(3), drdx(2), r_dot_t
    real(dp)       :: c0(3, 6), ci(npow, 3, 6), func(3, 6), ni(3), amat(npow, npow)
    real(dp)       :: eta, pw, ui(3, 2), us(3, 2), nu
    integer        :: i, j, k, psk

    d_eta = (eta_final - eta0) / real(npow, dp)
    call get_bound_elem_geom_points(be, begp)
    call get_plane_state_key(psk)
    call get_poisson_ratio(nu)

    normal_eta0 = normal_1d(begp, eta0)
    drdx(1) = -normal_eta0(2)
    drdx(2) =  normal_eta0(1)
    if (eta0 > eta_final) drdx = -drdx
    us = u_ijk_rp(drdx, nu, psk)
    call intrp_func_1d(ni, eta0, be%lamb_1, be%lamb_2)
    do k = 1, 3
      c0(1:3, 2*k-1:2*k) = ni(k) * us
    end do

    eta = eta0
    do i = 1, npow
      eta = eta + d_eta
      if (i == npow) eta = eta_final
      normal = normal_1d(begp, eta)
      ti(1) = -normal(2)
      ti(2) =  normal(1)
      if (eta0 > eta_final) ti = -ti
      fp   = cartesian_point_from_intrp_func_1d(begp, eta)
      ri   = distance_2d(sp, fp)
      drdx = ri(2:3) / ri(1)
      ui   = u_ijk_rp(drdx, nu, psk)
      call intrp_func_1d(ni, eta, be%lamb_1, be%lamb_2)
      do k = 1, 3
        func(1:3, 2*k-1:2*k) = ni(k) * ui
      end do
      r_dot_t     = dot_product(drdx, ti)
      ci(i, :, :) = ((func / abs(r_dot_t)) - c0) / ri(1)
      amat(i, 1)  = 1.0_dp
      do j = 2, npow
        amat(i, j) = amat(i, j-1) * ri(1)
      end do
    end do

    call mat_inverse(amat)
    do i = 1, 3
      do j = 1, 6
        ci(:, i, j) = matmul(amat, ci(:, i, j))
      end do
    end do

    gi = gi + c0 * log(ri(1))
    do i = 1, npow
      pw = real(i, dp)
      gi = gi + ci(i, :, :) * (ri(1)**pw) / pw
    end do
  end subroutine stress_u_subreg_sing_intgr

end module integration_mod
