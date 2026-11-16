!-----------------------------------------------------------------------
! MODULE: kelvin_sol_mod
!> @brief Soluções fundamentais 2D de Kelvin para elasticidade linear isotrópica
!>        (deslocamentos, trações, deformações, tensões, gradientes de
!>        deslocamento e suas partes regulares para integração singular de Gao).
!-----------------------------------------------------------------------
module kelvin_sol_mod
  use constants_mod
  implicit none

contains

  !---------------------------------------------------------------------
  !> @brief Tensor fundamental de deslocamento u*_ij(x_fonte, x_campo).
  !---------------------------------------------------------------------
  pure function u_ij(r, mu, nu, psk) result(u)
    real(dp), intent(in) :: r(3), mu, nu
    integer, intent(in)  :: psk
    real(dp) :: u(2, 2)
    real(dp) :: nu_b
    integer  :: i, j

    nu_b = merge(nu / (1.0_dp + nu), nu, psk == 1)

    do i = 1, 2
      do j = 1, 2
        u(i, j) = r(i+1) * r(j+1) / (r(1) * r(1))
        if (i == j) u(i, j) = u(i, j) + (4.0_dp * nu_b - 3.0_dp) * log(r(1))
      end do
    end do
    u = u / (eight_pi * mu * (1.0_dp - nu_b))
  end function u_ij

  !---------------------------------------------------------------------
  !> @brief Tensor fundamental de força de superfície (tração) t*_ij.
  !---------------------------------------------------------------------
  pure function t_ij(r, nu, ni, psk) result(t)
    real(dp), intent(in) :: r(3), ni(2), nu
    integer, intent(in)  :: psk
    real(dp) :: t(2, 2)
    real(dp) :: nu_b, drdn
    integer  :: i, j

    nu_b = merge(nu / (1.0_dp + nu), nu, psk == 1)
    drdn = dot_product(ni, r(2:3)) / r(1)

    do i = 1, 2
      do j = 1, 2
        t(i, j) = -2.0_dp * r(i+1) * r(j+1) / (r(1) * r(1))
        if (i == j) t(i, j) = t(i, j) + (2.0_dp * nu_b - 1.0_dp)
        t(i, j) = t(i, j) * drdn + (1.0_dp - 2.0_dp * nu_b) * (ni(j)*r(i+1) - ni(i)*r(j+1)) / r(1)
      end do
    end do
    t = t / (four_pi * (1.0_dp - nu_b) * r(1))
  end function t_ij

  !---------------------------------------------------------------------
  !> @brief Tensor u*_ijk para a Equação Integral de Tensões / Trações:
  !>        Sequência de linhas ij = (11, 22, 12).
  !---------------------------------------------------------------------
  pure function u_ijk(r, nu, psk) result(u)
    real(dp), intent(in) :: r(3), nu
    integer, intent(in)  :: psk
    real(dp) :: u(3, 2)
    real(dp) :: nu_b
    integer  :: i, j, k, ij

    nu_b = merge(nu / (1.0_dp + nu), nu, psk == 1)
    u = 0.0_dp

    do ij = 1, 3
      if (ij == 3) then
        i = 1; j = 2
      else
        i = ij; j = ij
      end if
      do k = 1, 2
        u(ij, k) = u(ij, k) + 2.0_dp * r(i+1) * r(j+1) * r(k+1) / (r(1)**3)
        if (i == j) u(ij, k) = u(ij, k) - (1.0_dp - 2.0_dp * nu_b) * r(k+1) / r(1)
        if (i == k) u(ij, k) = u(ij, k) + (1.0_dp - 2.0_dp * nu_b) * r(j+1) / r(1)
        if (j == k) u(ij, k) = u(ij, k) + (1.0_dp - 2.0_dp * nu_b) * r(i+1) / r(1)
      end do
    end do
    u = u / (four_pi * (1.0_dp - nu_b) * r(1))
  end function u_ijk

  !---------------------------------------------------------------------
  !> @brief Projeção normal n_j(fonte) * u*_ijk para a EIC de Tração.
  !---------------------------------------------------------------------
  pure function normal_u_ijk(r, nj, nu, psk) result(nu_ijk)
    real(dp), intent(in) :: r(3), nj(2), nu
    integer, intent(in)  :: psk
    real(dp) :: nu_ijk(2, 2), u(3, 2)

    u = u_ijk(r, nu, psk)
    nu_ijk(1, 1) = nj(1)*u(1, 1) + nj(2)*u(3, 1)
    nu_ijk(1, 2) = nj(1)*u(1, 2) + nj(2)*u(3, 2)
    nu_ijk(2, 1) = nj(2)*u(2, 1) + nj(1)*u(3, 1)
    nu_ijk(2, 2) = nj(2)*u(2, 2) + nj(1)*u(3, 2)
  end function normal_u_ijk

  !---------------------------------------------------------------------
  !> @brief Tensor t*_ijk para a Equação Integral de Tensões / Trações:
  !>        Sequência de linhas ij = (11, 22, 12).
  !---------------------------------------------------------------------
  pure function t_ijk(r, mu, nu, ni, psk) result(t)
    real(dp), intent(in) :: r(3), ni(2), mu, nu
    integer, intent(in)  :: psk
    real(dp) :: t(3, 2)
    real(dp) :: nu_b, drdn
    integer  :: i, j, k, ij

    nu_b = merge(nu / (1.0_dp + nu), nu, psk == 1)
    drdn = dot_product(ni, r(2:3)) / r(1)
    t = 0.0_dp

    do ij = 1, 3
      if (ij == 3) then
        i = 1; j = 2
      else
        i = ij; j = ij
      end if
      do k = 1, 2
        t(ij, k) = t(ij, k) + (1.0_dp - 2.0_dp*nu_b) * 2.0_dp * r(i+1)*r(j+1)*ni(k) / (r(1)*r(1))
        t(ij, k) = t(ij, k) - 8.0_dp * drdn * r(i+1)*r(j+1)*r(k+1) / (r(1)**3)
        t(ij, k) = t(ij, k) + 2.0_dp * nu_b * ni(j)*r(i+1)*r(k+1) / (r(1)*r(1))
        t(ij, k) = t(ij, k) + 2.0_dp * nu_b * ni(i)*r(j+1)*r(k+1) / (r(1)*r(1))
        if (i == j) t(ij, k) = t(ij, k) - (1.0_dp - 4.0_dp*nu_b)*ni(k) + 2.0_dp*(1.0_dp - 2.0_dp*nu_b)*drdn*r(k+1)/r(1)
        if (i == k) t(ij, k) = t(ij, k) + (1.0_dp - 2.0_dp*nu_b)*ni(j) + 2.0_dp*drdn*nu_b*r(j+1)/r(1)
        if (j == k) t(ij, k) = t(ij, k) + (1.0_dp - 2.0_dp*nu_b)*ni(i) + 2.0_dp*drdn*nu_b*r(i+1)/r(1)
      end do
    end do
    t = t * mu / (two_pi * (1.0_dp - nu_b) * r(1) * r(1))
  end function t_ijk

  !---------------------------------------------------------------------
  !> @brief Projeção normal n_j(fonte) * t*_ijk para a EIC Hipersingular.
  !---------------------------------------------------------------------
  pure function normal_t_ijk(r, mu, nu, ni, nj, psk) result(nt_ijk)
    real(dp), intent(in) :: r(3), ni(2), nj(2), mu, nu
    integer, intent(in)  :: psk
    real(dp) :: nt_ijk(2, 2), t(3, 2)

    t = t_ijk(r, mu, nu, ni, psk)
    nt_ijk(1, 1) = nj(1)*t(1, 1) + nj(2)*t(3, 1)
    nt_ijk(1, 2) = nj(1)*t(1, 2) + nj(2)*t(3, 2)
    nt_ijk(2, 1) = nj(2)*t(2, 1) + nj(1)*t(3, 1)
    nt_ijk(2, 2) = nj(2)*t(2, 2) + nj(1)*t(3, 2)
  end function normal_t_ijk

  !---------------------------------------------------------------------
  !> @brief Parte regular (independente de r) de t*_ij.
  !---------------------------------------------------------------------
  pure function t_ij_rp(drdx, ni, nu, psk) result(t_rp)
    real(dp), intent(in) :: drdx(2), ni(2), nu
    integer, intent(in)  :: psk
    real(dp) :: t_rp(2, 2)
    real(dp) :: nu_b, drdn, cte
    integer  :: i, j

    nu_b = merge(nu / (1.0_dp + nu), nu, psk == 1)
    drdn = dot_product(ni, drdx)
    cte  = 1.0_dp - 2.0_dp * nu_b

    do i = 1, 2
      do j = 1, 2
        t_rp(i, j) = -2.0_dp * drdx(i) * drdx(j)
        if (i == j) t_rp(i, j) = t_rp(i, j) - cte
        t_rp(i, j) = t_rp(i, j) * drdn + cte * (ni(j)*drdx(i) - ni(i)*drdx(j))
      end do
    end do
    t_rp = t_rp / (four_pi * (1.0_dp - nu_b))
  end function t_ij_rp

  !---------------------------------------------------------------------
  !> @brief Parte regular de u*_ijk.
  !---------------------------------------------------------------------
  pure function u_ijk_rp(drdx, nu, psk) result(u_rp)
    real(dp), intent(in) :: drdx(2), nu
    integer, intent(in)  :: psk
    real(dp) :: u_rp(3, 2)
    real(dp) :: nu_b, cte
    integer  :: i, j, k, ij

    nu_b = merge(nu / (1.0_dp + nu), nu, psk == 1)
    cte  = 1.0_dp - 2.0_dp * nu_b
    u_rp = 0.0_dp

    do ij = 1, 3
      if (ij == 3) then
        i = 1; j = 2
      else
        i = ij; j = ij
      end if
      do k = 1, 2
        u_rp(ij, k) = u_rp(ij, k) + 2.0_dp * drdx(i) * drdx(j) * drdx(k)
        if (i == j) u_rp(ij, k) = u_rp(ij, k) - cte * drdx(k)
        if (i == k) u_rp(ij, k) = u_rp(ij, k) + cte * drdx(j)
        if (j == k) u_rp(ij, k) = u_rp(ij, k) + cte * drdx(i)
      end do
    end do
    u_rp = u_rp / (four_pi * (1.0_dp - nu_b))
  end function u_ijk_rp

  !---------------------------------------------------------------------
  !> @brief Parte regular de t*_ijk.
  !---------------------------------------------------------------------
  pure function t_ijk_rp(drdx, ni, mu, nu, psk) result(t_rp)
    real(dp), intent(in) :: drdx(2), ni(2), mu, nu
    integer, intent(in)  :: psk
    real(dp) :: t_rp(3, 2)
    real(dp) :: nu_b, drdn, cte
    integer  :: i, j, k, ij

    nu_b = merge(nu / (1.0_dp + nu), nu, psk == 1)
    drdn = dot_product(ni, drdx)
    cte  = 1.0_dp - 2.0_dp * nu_b
    t_rp = 0.0_dp

    do ij = 1, 3
      if (ij == 3) then
        i = 1; j = 2
      else
        i = ij; j = ij
      end if
      do k = 1, 2
        t_rp(ij, k) = t_rp(ij, k) + 2.0_dp * cte * drdx(i) * drdx(j) * ni(k)
        t_rp(ij, k) = t_rp(ij, k) - 8.0_dp * drdn * drdx(i) * drdx(j) * drdx(k)
        t_rp(ij, k) = t_rp(ij, k) + 2.0_dp * nu_b * ni(j) * drdx(i) * drdx(k)
        t_rp(ij, k) = t_rp(ij, k) + 2.0_dp * nu_b * ni(i) * drdx(j) * drdx(k)
        if (i == j) t_rp(ij, k) = t_rp(ij, k) - (1.0_dp - 4.0_dp*nu_b)*ni(k) + 2.0_dp*cte*drdn*drdx(k)
        if (i == k) t_rp(ij, k) = t_rp(ij, k) + cte*ni(j) + 2.0_dp*drdn*nu_b*drdx(j)
        if (j == k) t_rp(ij, k) = t_rp(ij, k) + cte*ni(i) + 2.0_dp*drdn*nu_b*drdx(i)
      end do
    end do
    t_rp = t_rp * mu / (two_pi * (1.0_dp - nu_b))
  end function t_ijk_rp

  !---------------------------------------------------------------------
  !> @brief Parte regular de n_j(fonte) * u*_ijk.
  !---------------------------------------------------------------------
  pure function normal_u_ijk_rp(drdx, nj, nu, psk) result(nu_rp)
    real(dp), intent(in) :: drdx(2), nj(2), nu
    integer, intent(in)  :: psk
    real(dp) :: nu_rp(2, 2), u_rp(3, 2)

    u_rp = u_ijk_rp(drdx, nu, psk)
    nu_rp(1, 1) = nj(1)*u_rp(1, 1) + nj(2)*u_rp(3, 1)
    nu_rp(1, 2) = nj(1)*u_rp(1, 2) + nj(2)*u_rp(3, 2)
    nu_rp(2, 1) = nj(2)*u_rp(2, 1) + nj(1)*u_rp(3, 1)
    nu_rp(2, 2) = nj(2)*u_rp(2, 2) + nj(1)*u_rp(3, 2)
  end function normal_u_ijk_rp

  !---------------------------------------------------------------------
  !> @brief Parte regular de n_j(fonte) * t*_ijk.
  !---------------------------------------------------------------------
  pure function normal_t_ijk_rp(drdx, ni, nj, mu, nu, psk) result(nt_rp)
    real(dp), intent(in) :: drdx(2), ni(2), nj(2), mu, nu
    integer, intent(in)  :: psk
    real(dp) :: nt_rp(2, 2), t_rp(3, 2)

    t_rp = t_ijk_rp(drdx, ni, mu, nu, psk)
    nt_rp(1, 1) = nj(1)*t_rp(1, 1) + nj(2)*t_rp(3, 1)
    nt_rp(1, 2) = nj(1)*t_rp(1, 2) + nj(2)*t_rp(3, 2)
    nt_rp(2, 1) = nj(2)*t_rp(2, 1) + nj(1)*t_rp(3, 1)
    nt_rp(2, 2) = nj(2)*t_rp(2, 2) + nj(1)*t_rp(3, 2)
  end function normal_t_ijk_rp

  !---------------------------------------------------------------------
  !> @brief Tensor u*_ij,k para a Equação Integral do Gradiente de Deslocamento:
  !>        Sequência de linhas ik = (11, 22, 12, 21).
  !---------------------------------------------------------------------
  pure function u_ij_k(r, mu, nu, psk) result(u)
    real(dp), intent(in) :: r(3), mu, nu
    integer, intent(in)  :: psk
    real(dp) :: u(4, 2)
    real(dp) :: nu_b
    integer  :: i, j, k, ik

    nu_b = merge(nu / (1.0_dp + nu), nu, psk == 1)
    u = 0.0_dp

    do ik = 1, 4
      if (ik == 3) then
        i = 1; k = 2
      else if (ik == 4) then
        i = 2; k = 1
      else
        i = ik; k = ik
      end if
      do j = 1, 2
        u(ik, j) = u(ik, j) + 2.0_dp * r(i+1) * r(j+1) * r(k+1) / (r(1)**3)
        if (i == j) u(ik, j) = u(ik, j) + (3.0_dp - 4.0_dp * nu_b) * r(k+1) / r(1)
        if (i == k) u(ik, j) = u(ik, j) - r(j+1) / r(1)
        if (j == k) u(ik, j) = u(ik, j) - r(i+1) / r(1)
      end do
    end do
    u = u / (eight_pi * mu * (1.0_dp - nu_b) * r(1))
  end function u_ij_k

  !---------------------------------------------------------------------
  !> @brief Tensor t*_ij,k para a Equação Integral do Gradiente de Deslocamento:
  !>        Sequência de linhas ik = (11, 22, 12, 21).
  !---------------------------------------------------------------------
  pure function t_ij_k(r, nu, ni, psk) result(t)
    real(dp), intent(in) :: r(3), ni(2), nu
    integer, intent(in)  :: psk
    real(dp) :: t(4, 2)
    real(dp) :: nu_b, drdn
    integer  :: i, j, k, ik

    nu_b = merge(nu / (1.0_dp + nu), nu, psk == 1)
    drdn = dot_product(ni, r(2:3)) / r(1)
    t = 0.0_dp

    do ik = 1, 4
      if (ik == 3) then
        i = 1; k = 2
      else if (ik == 4) then
        i = 2; k = 1
      else
        i = ik; k = ik
      end if
      do j = 1, 2
        t(ik, j) = t(ik, j) + 8.0_dp * drdn * r(i+1) * r(j+1) * r(k+1) / (r(1)**3)
        t(ik, j) = t(ik, j) + (2.0_dp - 4.0_dp*nu_b) * ni(i) * r(j+1) * r(k+1) / (r(1)*r(1))
        t(ik, j) = t(ik, j) - (2.0_dp - 4.0_dp*nu_b) * ni(j) * r(i+1) * r(k+1) / (r(1)*r(1))
        t(ik, j) = t(ik, j) - 2.0_dp * ni(k) * r(i+1) * r(j+1) / (r(1)*r(1))
        if (i == j) t(ik, j) = t(ik, j) + (1.0_dp - 2.0_dp*nu_b) * (2.0_dp*drdn*(r(k+1)/r(1)) - ni(k))
        if (i == k) t(ik, j) = t(ik, j) + (1.0_dp - 2.0_dp*nu_b)*ni(j) - 2.0_dp*drdn*r(j+1)/r(1)
        if (j == k) t(ik, j) = t(ik, j) - (1.0_dp - 2.0_dp*nu_b)*ni(i) - 2.0_dp*drdn*r(i+1)/r(1)
      end do
    end do
    t = t / (four_pi * (nu_b - 1.0_dp) * r(1) * r(1))
  end function t_ij_k

  !---------------------------------------------------------------------
  !> @brief Parte regular de u*_ij,k.
  !---------------------------------------------------------------------
  pure function u_ij_k_rp(drdx, mu, nu, psk) result(u_rp)
    real(dp), intent(in) :: drdx(2), mu, nu
    integer, intent(in)  :: psk
    real(dp) :: u_rp(4, 2)
    real(dp) :: nu_b
    integer  :: i, j, k, ik

    nu_b = merge(nu / (1.0_dp + nu), nu, psk == 1)
    u_rp = 0.0_dp

    do ik = 1, 4
      if (ik == 3) then
        i = 1; k = 2
      else if (ik == 4) then
        i = 2; k = 1
      else
        i = ik; k = ik
      end if
      do j = 1, 2
        u_rp(ik, j) = u_rp(ik, j) + 2.0_dp * drdx(i) * drdx(j) * drdx(k)
        if (i == j) u_rp(ik, j) = u_rp(ik, j) + (3.0_dp - 4.0_dp * nu_b) * drdx(k)
        if (i == k) u_rp(ik, j) = u_rp(ik, j) - drdx(j)
        if (j == k) u_rp(ik, j) = u_rp(ik, j) - drdx(i)
      end do
    end do
    u_rp = u_rp / (eight_pi * mu * (1.0_dp - nu_b))
  end function u_ij_k_rp

  !---------------------------------------------------------------------
  !> @brief Parte regular de t*_ij,k.
  !---------------------------------------------------------------------
  pure function t_ij_k_rp(drdx, nu, ni, psk) result(t_rp)
    real(dp), intent(in) :: drdx(2), ni(2), nu
    integer, intent(in)  :: psk
    real(dp) :: t_rp(4, 2)
    real(dp) :: nu_b, drdn
    integer  :: i, j, k, ik

    nu_b = merge(nu / (1.0_dp + nu), nu, psk == 1)
    drdn = dot_product(ni, drdx)
    t_rp = 0.0_dp

    do ik = 1, 4
      if (ik == 3) then
        i = 1; k = 2
      else if (ik == 4) then
        i = 2; k = 1
      else
        i = ik; k = ik
      end if
      do j = 1, 2
        t_rp(ik, j) = t_rp(ik, j) + 8.0_dp * drdn * drdx(i) * drdx(j) * drdx(k)
        t_rp(ik, j) = t_rp(ik, j) + (2.0_dp - 4.0_dp*nu_b) * ni(i) * drdx(j) * drdx(k)
        t_rp(ik, j) = t_rp(ik, j) - (2.0_dp - 4.0_dp*nu_b) * ni(j) * drdx(i) * drdx(k)
        t_rp(ik, j) = t_rp(ik, j) - 2.0_dp * ni(k) * drdx(i) * drdx(j)
        if (i == j) t_rp(ik, j) = t_rp(ik, j) + (1.0_dp - 2.0_dp*nu_b) * (2.0_dp*drdn*drdx(k) - ni(k))
        if (i == k) t_rp(ik, j) = t_rp(ik, j) + (1.0_dp - 2.0_dp*nu_b)*ni(j) - 2.0_dp*drdn*drdx(j)
        if (j == k) t_rp(ik, j) = t_rp(ik, j) - (1.0_dp - 2.0_dp*nu_b)*ni(i) - 2.0_dp*drdn*drdx(i)
      end do
    end do
    t_rp = t_rp / (four_pi * (nu_b - 1.0_dp))
  end function t_ij_k_rp

end module kelvin_sol_mod
