module kelvin_fund_sol_2d
!----------------------------------------------------
! A module with Kelvin fundamental solution terms
! for two-dimensions.
!
! Date: June 2022
! Last update: July 2022
!----------------------------------------------------
!
implicit none
!
real(kind=8),private :: two_pi = 6.283185307179586_8
real(kind=8),private :: four_pi = 12.566370614359173_8
real(kind=8),private :: eight_pi = 25.132741228718346_8
!
contains
  function u_ij(r,mu,nu,psk)
    !----------------------------------------------------
    ! Compute the j-displacement components at a field
    ! point of an infinity isotropic elastic domain given
    ! by unitary concentrated loads acting at a source
    ! point in i-direction.
    !
    !  u_ij = [u_11 u_12]
    !         [u_21 u_22]
    !
    ! r(1) = distance from the point load
    ! r(2) = r_x, r(3) = r_y
    !
    ! mu = E/(2*(1+nu)), Shear Modulus
    !
    ! nu = Poisson ratio
    !
    ! plane_state_key (psk) = 1 for Plane Stress
    !                         2 for Plane Strain
    !
    ! Date: 24 June 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: r(3)
    real(kind=8),intent(in) :: mu,nu
    integer,intent(in) :: psk
    integer :: i,j
    real(kind=8) :: nu_b
    real(kind=8) :: u_ij(2,2)
    !
    if (psk.eq.1) then
      nu_b = nu/(1._8+nu)
    else
      nu_b = nu
    end if
    !
    do i = 1,2
      do j = 1,2
        u_ij(i,j) = r(i+1)*r(j+1)/(r(1)*r(1))
        if (i.eq.j) u_ij(i,j) = u_ij(i,j) + (4._8*nu_b-3._8)*log(r(1))
      end do
    end do
    !
    u_ij = u_ij / (eight_pi*mu*(1._8-nu_b))
    !
    return
  end function u_ij
  !****************************************************
  !****************************************************
  function t_ij(r,nu,ni,psk)
    !----------------------------------------------------
    ! Compute the j-traction components at a field
    ! point, on a surface oriented by the normal vector
    ! ni at this point, given by unitary concentrated
    ! loads acting in i-direction and at a source point
    ! of an infinity isotropic elastic domain.
    !
    !  t_ij = [t_11 t_12]
    !         [t_21 t_22]
    !
    ! r(1) = distance from the point load
    ! r(2) = r_x, r(3) = r_y
    !
    ! nu = Poisson ratio
    !
    ! plane_state_key (psk) = 1 for Plane Stress
    !                         2 for Plane Strain
    !
    ! Date: 24 June 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: r(3),ni(2)
    real(kind=8),intent(in) :: nu
    integer,intent(in) :: psk
    integer :: i,j
    real(kind=8) :: nu_b
    real(kind=8) :: t_ij(2,2)
    !
    if (psk.eq.1) then
      nu_b = nu/(1._8+nu)
    else
      nu_b = nu
    end if
    !
    do i = 1,2
      do j = 1,2
        t_ij(i,j) = -2._8*r(i+1)*r(j+1) / (r(1)*r(1))
        if (i.eq.j) t_ij(i,j) = t_ij(i,j) + (2._8*nu_b - 1._8)
        t_ij(i,j) = t_ij(i,j) * dot_product(ni,r(2:3)) / r(1)
        t_ij(i,j) = t_ij(i,j) + (1._8-(2._8*nu_b))*(ni(j)*r(i+1)-ni(i)*r(j+1)) / r(1)
      end do
    end do
    !
    t_ij = t_ij / (four_pi*(1._8-nu_b)*r(1))
    !
    return
  end function t_ij
  !****************************************************
  !****************************************************
  function ub_ijk(r,mu,nu,psk)
    !----------------------------------------------------
    ! A tensor for the strain integral equation:
    !
    ! ub_ijk = 0.5 * (u_ik,j + u_jk,i),
    !
    ! with the derivatives taken with respect to the
    ! source point.
    !
    !          [ub_111 ub_112]
    ! ub_ijk = [ub_221 ub_222]
    !          [ub_121 ub_122]
    !
    ! r(1) = distance from the point load
    ! r(2) = r_x, r(3) = r_y
    !
    ! mu = E/(2*(1+nu)), Shear Modulus
    !
    ! nu = Poisson ratio
    !
    ! plane_state_key (psk) = 1 for Plane Stress
    !                         2 for Plane Strain
    !
    ! Date: 24 June 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: r(3)
    real(kind=8),intent(in) :: mu,nu
    integer,intent(in) :: psk
    integer :: i,j,k,ij
    real(kind=8) :: nu_b
    real(kind=8) :: ub_ijk(3,2)
    !
    if (psk.eq.1) then
      nu_b = nu/(1._8+nu)
    else
      nu_b = nu
    end if
    !
    ub_ijk = 0.0
    !
    do ij = 1,3
      !The sequence is ij = 11, 22, 12
      if (ij.eq.3) then
        i = 1; j = 2
      else
        i = ij; j = ij
      end if
      do k = 1,2
        ub_ijk(ij,k) = ub_ijk(ij,k) + 2._8*r(i+1)*r(j+1)*r(k+1)/(r(1)*r(1)*r(1))
        if (i.eq.j) ub_ijk(ij,k) = ub_ijk(ij,k) - r(k+1)/r(1)
        if (i.eq.k) ub_ijk(ij,k) = ub_ijk(ij,k) + (1._8-2._8*nu_b)*r(j+1)/r(1)
        if (j.eq.k) ub_ijk(ij,k) = ub_ijk(ij,k) + (1._8-2._8*nu_b)*r(i+1)/r(1)
      end do
    end do
    !
    ub_ijk = ub_ijk / (eight_pi*mu*(1._8-nu_b)*r(1))
    !
    return
  end function ub_ijk
  !****************************************************
  !****************************************************
  function tb_ijk(r,nu,ni,psk)
    !----------------------------------------------------
    ! A tensor for the strain integral equation:
    !
    ! tb_ijk = 0.5 * (t_ik,j + t_jk,i),
    !
    ! with the derivatives taken with respect to the
    ! source point.
    !
    !          [tb_111 tb_112]
    ! tb_ijk = [tb_221 tb_222]
    !          [tb_121 tb_122]
    !
    ! r(1) = distance from the point load
    ! r(2) = r_x, r(3) = r_y
    !
    ! nu = Poisson ratio
    !
    ! plane_state_key (psk) = 1 for Plane Stress
    !                         2 for Plane Strain
    !
    ! Date: 24 June 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: r(3),ni(2)
    real(kind=8),intent(in) :: nu
    integer,intent(in) :: psk
    integer :: i,j,k,ij
    real(kind=8) :: nu_b,drdn
    real(kind=8) :: tb_ijk(3,2)
    !
    if (psk.eq.1) then
      nu_b = nu/(1._8+nu)
    else
      nu_b = nu
    end if
    !
    drdn = dot_product(ni,r(2:3)) / r(1)
    !
    tb_ijk = 0.0
    !
    do ij = 1,3
      !The sequence is ij = 11, 22, 12
      if (ij.eq.3) then
        i = 1; j = 2
      else
        i = ij; j = ij
      end if
      do k = 1,2
        tb_ijk(ij,k) = tb_ijk(ij,k) + (1._8-2._8*nu_b)*2._8*r(i+1)*r(j+1)*ni(k)/(r(1)*r(1))
        tb_ijk(ij,k) = tb_ijk(ij,k) - 8._8*drdn*r(i+1)*r(j+1)*r(k+1)/(r(1)*r(1)*r(1))
        tb_ijk(ij,k) = tb_ijk(ij,k) + 2._8*nu_b*ni(j)*r(i+1)*r(k+1)/(r(1)*r(1))
        tb_ijk(ij,k) = tb_ijk(ij,k) + 2._8*nu_b*ni(i)*r(j+1)*r(k+1)/(r(1)*r(1))
        if (i.eq.j) tb_ijk(ij,k) = tb_ijk(ij,k) - (1._8-2._8*nu_b)*ni(k) + 2._8*drdn*r(k+1)/r(1)
        if (i.eq.k) tb_ijk(ij,k) = tb_ijk(ij,k) + (1._8-2._8*nu_b)*ni(j) + 2._8*drdn*nu_b*r(j+1)/r(1)
        if (j.eq.k) tb_ijk(ij,k) = tb_ijk(ij,k) + (1._8-2._8*nu_b)*ni(i) + 2._8*drdn*nu_b*r(i+1)/r(1)
      end do
    end do
    !
    tb_ijk = tb_ijk / (four_pi*(1._8-nu_b)*r(1)*r(1))
    !
    return
  end function tb_ijk
  !****************************************************
  !****************************************************
  function u_ijk(r,nu,psk)
    !----------------------------------------------------
    ! A tensor for the stress integral equation:
    !
    ! u_ijk = C_ijrs : ub_rsk,
    !
    ! where C_ijrs is the fourth order isotropic elastic
    ! constitutive tensor.
    !
    !         [u_111 u_112]
    ! u_ijk = [u_221 u_222]
    !         [u_121 u_122]
    !
    ! r(1) = distance from the point load
    ! r(2) = r_x, r(3) = r_y
    !
    ! nu = Poisson ratio
    !
    ! plane_state_key (psk) = 1 for Plane Stress
    !                         2 for Plane Strain
    !
    ! Date: 24 June 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: r(3)
    real(kind=8),intent(in) :: nu
    integer,intent(in) :: psk
    integer :: i,j,k,ij
    real(kind=8) :: nu_b
    real(kind=8) :: u_ijk(3,2)
    !
    if (psk.eq.1) then
      nu_b = nu/(1._8+nu)
    else
      nu_b = nu
    end if
    !
    u_ijk = 0.0
    !
    do ij = 1,3
      !The sequence is ij = 11, 22, 12
      if (ij.eq.3) then
        i = 1; j = 2
      else
        i = ij; j = ij
      end if
      do k = 1,2
        u_ijk(ij,k) = u_ijk(ij,k) + 2._8*r(i+1)*r(j+1)*r(k+1)/(r(1)*r(1)*r(1))
        if (i.eq.j) u_ijk(ij,k) = u_ijk(ij,k) - (1._8-2._8*nu_b)*r(k+1)/r(1)
        if (i.eq.k) u_ijk(ij,k) = u_ijk(ij,k) + (1._8-2._8*nu_b)*r(j+1)/r(1)
        if (j.eq.k) u_ijk(ij,k) = u_ijk(ij,k) + (1._8-2._8*nu_b)*r(i+1)/r(1)
      end do
    end do
    !
    u_ijk = u_ijk / (four_pi*(1._8-nu_b)*r(1))
    !
    return
  end function u_ijk
  !****************************************************
  !****************************************************
  function normal_u_ijk(r,nj,nu,psk)
    !----------------------------------------------------
    ! A tensor for traction boundary integral equation:
    !
    ! n_j : u_ijk,
    !
    ! where n_j is the normal vector at the source point.
    !
    ! normal_u_ijk = [n_1*u_111+n_2*u_121 n_1*u_112+n_2*u_122]
    !                [n_2*u_221+n_1*u_121 n_2*u_222+n_1*u_122]
    !
    ! r(1) = distance from the point load
    ! r(2) = r_x, r(3) = r_y
    !
    ! nj = normal vector at the source point
    !
    ! nu = Poisson ratio
    !
    ! plane_state_key (psk) = 1 for Plane Stress
    !                         2 for Plane Strain
    !
    ! Date: 08 July 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: r(3)
    real(kind=8),intent(in) :: nj(2)
    real(kind=8),intent(in) :: nu
    integer,intent(in) :: psk
    real(kind=8) :: u(3,2)
    real(kind=8) :: normal_u_ijk(2,2)
    !
    u = u_ijk(r,nu,psk)
    !
    normal_u_ijk (1,1) = nj(1)*u(1,1) + nj(2)*u(3,1)
    normal_u_ijk (1,2) = nj(1)*u(1,2) + nj(2)*u(3,2)
    normal_u_ijk (2,1) = nj(2)*u(2,1) + nj(1)*u(3,1)
    normal_u_ijk (2,2) = nj(2)*u(2,2) + nj(1)*u(3,2)
    !
    return
  end function normal_u_ijk
  !****************************************************
  !****************************************************
  function t_ijk(r,mu,nu,ni,psk)
    !----------------------------------------------------
    ! A tensor for the stress integral equation:
    !
    ! t_ijk = C_ijrs : tb_rsk,
    !
    ! where C_ijrs is the fourth order isotropic elastic
    ! constitutive tensor.
    !
    !         [t_111 t_112]
    ! t_ijk = [t_221 t_222]
    !         [t_121 t_122]
    !
    ! r(1) = distance from the point load
    ! r(2) = r_x, r(3) = r_y
    !
    ! mu = E/(2*(1+nu)), Shear Modulus
    !
    ! nu = Poisson ratio
    !
    ! plane_state_key (psk) = 1 for Plane Stress
    !                         2 for Plane Strain
    !
    ! Date: 24 June 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: r(3),ni(2)
    real(kind=8),intent(in) :: mu,nu
    integer,intent(in) :: psk
    integer :: i,j,k,ij
    real(kind=8) :: nu_b,drdn
    real(kind=8) :: t_ijk(3,2)
    !
    if (psk.eq.1) then
      nu_b = nu/(1._8+nu)
    else
      nu_b = nu
    end if
    !
    drdn = dot_product(ni,r(2:3)) / r(1)
    !
    t_ijk = 0.0
    !
    do ij = 1,3
      !The sequence is ij = 11, 22, 12
      if (ij.eq.3) then
        i = 1; j = 2
      else
        i = ij; j = ij
      end if
      do k = 1,2
        t_ijk(ij,k) = t_ijk(ij,k) + (1._8-2._8*nu_b)*2._8*r(i+1)*r(j+1)*ni(k)/(r(1)*r(1))
        t_ijk(ij,k) = t_ijk(ij,k) - 8._8*drdn*r(i+1)*r(j+1)*r(k+1)/(r(1)*r(1)*r(1))
        t_ijk(ij,k) = t_ijk(ij,k) + 2._8*nu_b*ni(j)*r(i+1)*r(k+1)/(r(1)*r(1))
        t_ijk(ij,k) = t_ijk(ij,k) + 2._8*nu_b*ni(i)*r(j+1)*r(k+1)/(r(1)*r(1))
        if (i.eq.j) t_ijk(ij,k) = t_ijk(ij,k) - (1._8-4._8*nu_b)*ni(k)
        if (i.eq.j) t_ijk(ij,k) = t_ijk(ij,k) + 2._8*(1._8-2._8*nu_b)*drdn*r(k+1)/r(1)
        if (i.eq.k) t_ijk(ij,k) = t_ijk(ij,k) + (1._8-2._8*nu_b)*ni(j) + 2._8*drdn*nu_b*r(j+1)/r(1)
        if (j.eq.k) t_ijk(ij,k) = t_ijk(ij,k) + (1._8-2._8*nu_b)*ni(i) + 2._8*drdn*nu_b*r(i+1)/r(1)
      end do
    end do
    !
    t_ijk = t_ijk * mu/(two_pi*(1._8-nu_b)*r(1)*r(1))
    !
    return
  end function t_ijk
  !****************************************************
  !****************************************************
  function normal_t_ijk(r,mu,nu,ni,nj,psk)
    !----------------------------------------------------
    ! A tensor for traction boundary integral equation:
    !
    ! n_j : t_ijk,
    !
    ! where n_j is the normal vector at the source point.
    !
    ! normal_t_ijk = [n_1*t_111+n_2*t_121 n_1*t_112+n_2*t_122]
    !                [n_2*t_221+n_1*t_121 n_2*t_222+n_1*t_122]
    !
    ! r(1) = distance from the point load
    ! r(2) = r_x, r(3) = r_y
    !
    ! ni = normal vector at the field point
    !
    ! nj = normal vector at the source point
    !
    ! mu = E/(2*(1+nu)), Shear Modulus
    !
    ! nu = Poisson ratio
    !
    ! plane_state_key (psk) = 1 for Plane Stress
    !                         2 for Plane Strain
    !
    ! Date: 24 June 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: r(3),ni(2),nj(2)
    real(kind=8),intent(in) :: mu,nu
    integer,intent(in) :: psk
    real(kind=8) :: t(3,2)
    real(kind=8) :: normal_t_ijk(2,2)
    !
    t = t_ijk(r,mu,nu,ni,psk)
    !
    normal_t_ijk (1,1) = nj(1)*t(1,1) + nj(2)*t(3,1)
    normal_t_ijk (1,2) = nj(1)*t(1,2) + nj(2)*t(3,2)
    normal_t_ijk (2,1) = nj(2)*t(2,1) + nj(1)*t(3,1)
    normal_t_ijk (2,2) = nj(2)*t(2,2) + nj(1)*t(3,2)
    !
    return
  end function normal_t_ijk
  !****************************************************
  !****************************************************
  function t_ij_rp(drdx,ni,nu,psk)
    !----------------------------------------------------
    ! Compute the regular part of fundamental solution
    ! tensor t_ij
    !
    !  t_ij = [t_11 t_12]
    !         [t_21 t_22]
    !
    ! drdx(1) = r,_1, drdx(2) = r,_2
    !
    ! nu = Poisson ratio
    !
    ! plane_state_key (psk) = 1 for Plane Stress
    !                         2 for Plane Strain
    !
    ! Date: 25 January 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: drdx(2),ni(2)
    real(kind=8),intent(in) :: nu
    integer,intent(in) :: psk
    integer :: i,j
    real(kind=8) :: nu_b,drdn,cte
    real(kind=8) :: t_ij_rp(2,2)
    !
    if (psk.eq.1) then
      nu_b = nu/(1._8+nu)
    else
      nu_b = nu
    end if
    !
    drdn = dot_product(ni,drdx)
    cte = 1._8 - (2._8*nu_b)
    !
    do i = 1,2
      do j = 1,2
        t_ij_rp(i,j) = -2._8*drdx(i)*drdx(j)
        if (i.eq.j) t_ij_rp(i,j) = t_ij_rp(i,j) - cte
        t_ij_rp(i,j) = t_ij_rp(i,j) * drdn
        t_ij_rp(i,j) = t_ij_rp(i,j) + cte*(ni(j)*drdx(i) - ni(i)*drdx(j))
      end do
    end do
    !
    t_ij_rp = t_ij_rp / (four_pi*(1._8-nu_b))
    !
    return
  end function t_ij_rp
  !****************************************************
  !****************************************************
  function u_ijk_rp(drdx,nu,psk)
    !----------------------------------------------------
    ! Compute the regular part of fundamental solution
    ! tensor u_ijk
    !
    ! drdx(1) = r,_1, drdx(2) = r,_2
    !
    ! nu = Poisson ratio
    !
    ! plane_state_key (psk) = 1 for Plane Stress
    !                         2 for Plane Strain
    !
    ! Date: 19 July 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: drdx(2)
    real(kind=8),intent(in) :: nu
    integer,intent(in) :: psk
    integer :: i,j,k,ij
    real(kind=8) :: nu_b,cte
    real(kind=8) :: u_ijk_rp(3,2)
    !
    if (psk.eq.1) then
      nu_b = nu/(1._8+nu)
    else
      nu_b = nu
    end if
    !
    u_ijk_rp = 0.0
    !
    cte = 1._8 - (2._8*nu_b)
    !
    do ij = 1,3
      !The sequence is ij = 11, 22, 12
      if (ij.eq.3) then
        i = 1; j = 2
      else
        i = ij; j = ij
      end if
      do k = 1,2
        u_ijk_rp(ij,k) = u_ijk_rp(ij,k) + 2._8*drdx(i)*drdx(j)*drdx(k)
        if (i.eq.j) u_ijk_rp(ij,k) = u_ijk_rp(ij,k) - cte*drdx(k)
        if (i.eq.k) u_ijk_rp(ij,k) = u_ijk_rp(ij,k) + cte*drdx(j)
        if (j.eq.k) u_ijk_rp(ij,k) = u_ijk_rp(ij,k) + cte*drdx(i)
      end do
    end do
    !
    u_ijk_rp = u_ijk_rp / (four_pi*(1._8-nu_b))
    !
    return
  end function u_ijk_rp
  !****************************************************
  !****************************************************
  function t_ijk_rp(drdx,ni,mu,nu,psk)
    !----------------------------------------------------
    ! Compute the regular part of fundamental solution
    ! tensor t_ijk
    !
    ! drdx(1) = r,_1, drdx(2) = r,_2
    !
    ! ni = the normal vector at the field point
    !
    ! mu = E/(2*(1+nu)), Shear Modulus
    !
    ! nu = Poisson ratio
    !
    ! plane_state_key (psk) = 1 for Plane Stress
    !                         2 for Plane Strain
    !
    ! Date: 19 July 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: drdx(2),ni(2)
    real(kind=8),intent(in) :: mu,nu
    integer,intent(in) :: psk
    integer :: i,j,k,ij
    real(kind=8) :: nu_b,drdn,cte
    real(kind=8) :: t_ijk_rp(3,2)
    !
    if (psk.eq.1) then
      nu_b = nu/(1._8+nu)
    else
      nu_b = nu
    end if
    !
    drdn = dot_product(ni,drdx)
    cte = 1._8 - (2._8*nu_b)
    !
    t_ijk_rp = 0.0
    !
    do ij = 1,3
      !The sequence is ij = 11, 22, 12
      if (ij.eq.3) then
        i = 1; j = 2
      else
        i = ij; j = ij
      end if
      do k = 1,2
        t_ijk_rp(ij,k) = t_ijk_rp(ij,k) + 2._8*cte*drdx(i)*drdx(j)*ni(k)
        t_ijk_rp(ij,k) = t_ijk_rp(ij,k) - 8._8*drdn*drdx(i)*drdx(j)*drdx(k)
        t_ijk_rp(ij,k) = t_ijk_rp(ij,k) + 2._8*nu_b*ni(j)*drdx(i)*drdx(k)
        t_ijk_rp(ij,k) = t_ijk_rp(ij,k) + 2._8*nu_b*ni(i)*drdx(j)*drdx(k)
        if (i.eq.j) t_ijk_rp(ij,k) = t_ijk_rp(ij,k) - (1._8-4._8*nu_b)*ni(k)
        if (i.eq.j) t_ijk_rp(ij,k) = t_ijk_rp(ij,k) + 2._8*cte*drdn*drdx(k)
        if (i.eq.k) t_ijk_rp(ij,k) = t_ijk_rp(ij,k) + cte*ni(j) + 2._8*drdn*nu_b*drdx(j)
        if (j.eq.k) t_ijk_rp(ij,k) = t_ijk_rp(ij,k) + cte*ni(i) + 2._8*drdn*nu_b*drdx(i)
      end do
    end do
    !
    t_ijk_rp = t_ijk_rp * mu/(two_pi*(1._8-nu_b))
    !
    return
  end function t_ijk_rp
  !****************************************************
  !****************************************************
  function normal_u_ijk_rp(drdx,nj,nu,psk)
    !----------------------------------------------------
    ! Compute the regular part of fundamental solution
    ! tensor n_j : u_ijk
    !
    ! drdx(1) = r,_1, drdx(2) = r,_2
    !
    ! nj = normal vector at the source point
    !
    ! nu = Poisson ratio
    !
    ! plane_state_key (psk) = 1 for Plane Stress
    !                         2 for Plane Strain
    !
    ! Date: 08 July 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: drdx(2)
    real(kind=8),intent(in) :: nj(2)
    real(kind=8),intent(in) :: nu
    integer,intent(in) :: psk
    real(kind=8) :: normal_u_ijk_rp(2,2)
    real(kind=8) :: u_rp(3,2)
    !
    u_rp = u_ijk_rp(drdx,nu,psk)
    !
    normal_u_ijk_rp (1,1) = nj(1)*u_rp(1,1) + nj(2)*u_rp(3,1)
    normal_u_ijk_rp (1,2) = nj(1)*u_rp(1,2) + nj(2)*u_rp(3,2)
    normal_u_ijk_rp (2,1) = nj(2)*u_rp(2,1) + nj(1)*u_rp(3,1)
    normal_u_ijk_rp (2,2) = nj(2)*u_rp(2,2) + nj(1)*u_rp(3,2)
    !
    return
  end function normal_u_ijk_rp
  !****************************************************
  !****************************************************
  function normal_t_ijk_rp(drdx,ni,nj,mu,nu,psk)
    !----------------------------------------------------
    ! Compute the regular part of fundamental solution
    ! tensor t_ijk
    !
    !         [t_111 t_112]
    ! t_ijk = [t_221 t_222]
    !         [t_121 t_122]
    !
    ! drdx(1) = r,_1, drdx(2) = r,_2
    !
    ! ni = normal vector at the field point
    !
    ! nj = normal vector at the source point
    !
    ! mu = E/(2*(1+nu)), Shear Modulus
    !
    ! nu = Poisson ratio
    !
    ! plane_state_key (psk) = 1 for Plane Stress
    !                         2 for Plane Strain
    !
    ! Date: 19 July 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: drdx(2),ni(2),nj(2)
    real(kind=8),intent(in) :: mu,nu
    integer,intent(in) :: psk
    real(kind=8) :: t_rp(3,2)
    real(kind=8) :: normal_t_ijk_rp(2,2)
    !
    t_rp = t_ijk_rp(drdx,ni,mu,nu,psk)
    !
    normal_t_ijk_rp (1,1) = nj(1)*t_rp(1,1) + nj(2)*t_rp(3,1)
    normal_t_ijk_rp (1,2) = nj(1)*t_rp(1,2) + nj(2)*t_rp(3,2)
    normal_t_ijk_rp (2,1) = nj(2)*t_rp(2,1) + nj(1)*t_rp(3,1)
    normal_t_ijk_rp (2,2) = nj(2)*t_rp(2,2) + nj(1)*t_rp(3,2)
    !
    return
  end function normal_t_ijk_rp
  !****************************************************
  !****************************************************
  function u_ij_k(r,mu,nu,psk)
    !----------------------------------------------------
    ! A tensor for the displacement gradient equation:
    !
    !          [u_11_1 u_12_1]
    ! u_ij_k = [u_21_2 u_22_2]
    !          [u_11_2 u_12_2]
    !          [u_21_1 u_22_1]
    !
    ! r(1) = distance from the point load
    ! r(2) = r_x, r(3) = r_y
    !
    ! mu = E/(2*(1+nu)), Shear Elasticity Modulus
    !
    ! nu = Poisson ratio
    !
    ! plane_state_key (psk) = 1 for Plane Stress
    !                         2 for Plane Strain
    !
    ! Date: 16 September 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: r(3)
    real(kind=8),intent(in) :: mu,nu
    integer,intent(in) :: psk
    integer :: i,j,k,ik
    real(kind=8) :: nu_b
    real(kind=8) :: u_ij_k(4,2)
    !
    if (psk.eq.1) then
      nu_b = nu/(1._8+nu)
    else
      nu_b = nu
    end if
    !
    u_ij_k = 0.0
    !
    do ik = 1,4
      !The sequence is ik = 11, 22, 12, 21
      if (ik.eq.3) then
        i = 1; k = 2
      else if (ik.eq.4) then
        i = 2; k = 1
      else
        i = ik; k = ik
      end if
      do j = 1,2
        u_ij_k(ik,j) = u_ij_k(ik,j) + 2._8*r(i+1)*r(j+1)*r(k+1)/(r(1)*r(1)*r(1))
        if (i.eq.j) u_ij_k(ik,j) = u_ij_k(ik,j) + (3._8-4._8*nu_b)*r(k+1)/r(1)
        if (i.eq.k) u_ij_k(ik,j) = u_ij_k(ik,j) - r(j+1)/r(1)
        if (j.eq.k) u_ij_k(ik,j) = u_ij_k(ik,j) - r(i+1)/r(1)
      end do
    end do
    !
    u_ij_k = u_ij_k / (eight_pi*mu*(1._8-nu_b)*r(1))
    !
    return
  end function u_ij_k
  !****************************************************
  !****************************************************
  function t_ij_k(r,nu,ni,psk)
    !----------------------------------------------------
    ! A tensor for the displacement gradient equation:
    !
    !          [t_11_1 t_12_1]
    ! t_ij_k = [t_21_2 t_22_2]
    !          [t_11_2 t_12_2]
    !          [t_21_1 t_22_1]
    !
    ! r(1) = distance from the point load
    ! r(2) = r_x, r(3) = r_y
    !
    ! nu = Poisson ratio
    !
    ! plane_state_key (psk) = 1 for Plane Stress
    !                         2 for Plane Strain
    !
    ! Date: 16 September 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: r(3),ni(2)
    real(kind=8),intent(in) :: nu
    integer,intent(in) :: psk
    integer :: i,j,k,ik
    real(kind=8) :: nu_b,drdn
    real(kind=8) :: t_ij_k(4,2)
    !
    if (psk.eq.1) then
      nu_b = nu/(1._8+nu)
    else
      nu_b = nu
    end if
    !
    drdn = dot_product(ni,r(2:3)) / r(1)
    !
    t_ij_k = 0.0
    !
    do ik = 1,4
      !The sequence is ik = 11, 22, 12, 21
      if (ik.eq.3) then
        i = 1; k = 2
      else if (ik.eq.4) then
        i = 2; k = 1
      else
        i = ik; k = ik
      end if
      do j = 1,2
        t_ij_k(ik,j) = t_ij_k(ik,j) + 8._8*drdn*r(i+1)*r(j+1)*r(k+1)/(r(1)*r(1)*r(1))
        t_ij_k(ik,j) = t_ij_k(ik,j) + (2._8-4._8*nu_b)*ni(i)*r(j+1)*r(k+1)/(r(1)*r(1))
        t_ij_k(ik,j) = t_ij_k(ik,j) - (2._8-4._8*nu_b)*ni(j)*r(i+1)*r(k+1)/(r(1)*r(1))
        t_ij_k(ik,j) = t_ij_k(ik,j) - 2._8*ni(k)*r(i+1)*r(j+1)/(r(1)*r(1))
        if (i.eq.j) t_ij_k(ik,j) = t_ij_k(ik,j) + (1._8-2._8*nu_b)*(2._8*drdn*(r(k+1)/r(1))-ni(k))
        if (i.eq.k) t_ij_k(ik,j) = t_ij_k(ik,j) + (1._8-2._8*nu_b)*ni(j) - 2._8*drdn*r(j+1)/r(1)
        if (j.eq.k) t_ij_k(ik,j) = t_ij_k(ik,j) - (1._8-2._8*nu_b)*ni(i) - 2._8*drdn*r(i+1)/r(1)
      end do
    end do
    !
    t_ij_k = t_ij_k / (four_pi*(nu_b-1._8)*r(1)*r(1))
    !
    return
  end function t_ij_k
  !****************************************************
  !****************************************************
    function u_ij_k_rp(drdx,mu,nu,psk)
    !----------------------------------------------------
    ! A tensor for the displacement gradient equation:
    !
    !          [u_11_1 u_12_1]
    ! u_ij_k = [u_21_2 u_22_2]
    !          [u_11_2 u_12_2]
    !          [u_21_1 u_22_1]
    !
    ! r(1) = distance from the point load
    ! r(2) = r_x, r(3) = r_y
    !
    ! mu = E/(2*(1+nu)), Shear Elasticity Modulus
    !
    ! nu = Poisson ratio
    !
    ! plane_state_key (psk) = 1 for Plane Stress
    !                         2 for Plane Strain
    !
    ! Date: 22 September 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: drdx(2)
    real(kind=8),intent(in) :: mu,nu
    integer,intent(in) :: psk
    integer :: i,j,k,ik
    real(kind=8) :: nu_b
    real(kind=8) :: u_ij_k_rp(4,2)
    !
    if (psk.eq.1) then
      nu_b = nu/(1._8+nu)
    else
      nu_b = nu
    end if
    !
    u_ij_k_rp = 0.0
    !
    do ik = 1,4
      !The sequence is ik = 11, 22, 12, 21
      if (ik.eq.3) then
        i = 1; k = 2
      else if (ik.eq.4) then
        i = 2; k = 1
      else
        i = ik; k = ik
      end if
      do j = 1,2
        u_ij_k_rp(ik,j) = u_ij_k_rp(ik,j) + 2._8*drdx(i)*drdx(j)*drdx(k)
        if (i.eq.j) u_ij_k_rp(ik,j) = u_ij_k_rp(ik,j) + (3._8-4._8*nu_b)*drdx(k)
        if (i.eq.k) u_ij_k_rp(ik,j) = u_ij_k_rp(ik,j) - drdx(j)
        if (j.eq.k) u_ij_k_rp(ik,j) = u_ij_k_rp(ik,j) - drdx(i)
      end do
    end do
    !
    u_ij_k_rp = u_ij_k_rp / (eight_pi*mu*(1._8-nu_b))
    !
    return
  end function u_ij_k_rp
  !****************************************************
  !****************************************************
  function t_ij_k_rp(drdx,nu,ni,psk)
    !----------------------------------------------------
    ! A tensor for the displacement gradient equation:
    !
    !          [t_11_1 t_12_1]
    ! t_ij_k = [t_21_2 t_22_2]
    !          [t_11_2 t_12_2]
    !          [t_21_1 t_22_1]
    !
    ! r(1) = distance from the point load
    ! r(2) = r_x, r(3) = r_y
    !
    ! nu = Poisson ratio
    !
    ! plane_state_key (psk) = 1 for Plane Stress
    !                         2 for Plane Strain
    !
    ! Date: 16 September 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: drdx(2),ni(2)
    real(kind=8),intent(in) :: nu
    integer,intent(in) :: psk
    integer :: i,j,k,ik
    real(kind=8) :: nu_b,drdn
    real(kind=8) :: t_ij_k_rp(4,2)
    !
    if (psk.eq.1) then
      nu_b = nu/(1._8+nu)
    else
      nu_b = nu
    end if
    !
    drdn = dot_product(ni,drdx)
    !
    t_ij_k_rp = 0.0
    !
    do ik = 1,4
      !The sequence is ik = 11, 22, 12, 21
      if (ik.eq.3) then
        i = 1; k = 2
      else if (ik.eq.4) then
        i = 2; k = 1
      else
        i = ik; k = ik
      end if
      do j = 1,2
        t_ij_k_rp(ik,j) = t_ij_k_rp(ik,j) + 8._8*drdn*drdx(i)*drdx(j)*drdx(k)
        t_ij_k_rp(ik,j) = t_ij_k_rp(ik,j) + (2._8-4._8*nu_b)*ni(i)*drdx(j)*drdx(k)
        t_ij_k_rp(ik,j) = t_ij_k_rp(ik,j) - (2._8-4._8*nu_b)*ni(j)*drdx(i)*drdx(k)
        t_ij_k_rp(ik,j) = t_ij_k_rp(ik,j) - 2._8*ni(k)*drdx(i)*drdx(j)
        if (i.eq.j) t_ij_k_rp(ik,j) = t_ij_k_rp(ik,j) + (1._8-2._8*nu_b)*(2._8*drdn*(drdx(k))-ni(k))
        if (i.eq.k) t_ij_k_rp(ik,j) = t_ij_k_rp(ik,j) + (1._8-2._8*nu_b)*ni(j) - 2._8*drdn*drdx(j)
        if (j.eq.k) t_ij_k_rp(ik,j) = t_ij_k_rp(ik,j) - (1._8-2._8*nu_b)*ni(i) - 2._8*drdn*drdx(i)
      end do
    end do
    !
    t_ij_k_rp = t_ij_k_rp / (four_pi*(nu_b-1._8))
    !
    return
  end function t_ij_k_rp
  !****************************************************
  !****************************************************
end module kelvin_fund_sol_2d
