module fracture_2d_rad_exp_intgr
!----------------------------------------------------
! A module that performs numerical integrations by
! the Radial Expansion Integration Method.
!
! Part of REIM_2D program.
!
! Date: July 2022
! Last update: ---
!----------------------------------------------------
use error_messages
use geometry
use matrix_algebra
use fracture_2d_model
use numerical_integration_points
use kelvin_fund_sol_2d
!
implicit none
!
real(kind=8),private :: tol = 1.0e-4
!
integer,private :: npow !The order for radial polynomial expansion
!
contains
  subroutine displ_bie_sing_integr(sp,be,lcp,hi)
    !----------------------------------------------------
    ! Perform singular integrations for a given source
    ! point (sp) and a given boundary element (be).
    !
    ! lcp = location of the source point inside be.
    !
    ! Date: 23 July 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(inout) :: hi(:,:)
    type(boundary_elem),intent(in) :: be
    type(point_2d),intent(in) :: sp
    integer,intent(in) :: lcp
    real(kind=8) :: eta0,eta_final
    !
    hi = 0.0
    !
    !-----------------------------------------------------
    ! Storing the natural coordinate of the source point:
    !-----------------------------------------------------
    !
    if (lcp.eq.1) then
      eta0 = -be%lamb_1
    else if (lcp.eq.be%nodes_num) then
      eta0 = be%lamb_2
    else
      eta0 = 0.0_8 !Mid-node for quadratic elements
    end if
    !
    !-------------------------------------------
    ! Integrate region before the source point:
    !-------------------------------------------
    !
    if (eta0.gt.-1.0_8) then
      eta_final = -1.0_8
      call displ_subreg_sing_intgr(sp,be,eta0,eta_final,hi)
    end if
    !
    !------------------------------------------
    ! Integrate region after the source point:
    !------------------------------------------
    !
    if (eta0.lt.1.0_8) then
      eta_final = 1.0_8
      call displ_subreg_sing_intgr(sp,be,eta0,eta_final,hi)
    end if
    !
    return
  end subroutine displ_bie_sing_integr
  !****************************************************
  !****************************************************
  subroutine displ_subreg_sing_intgr(sp,be,eta0,eta_final,hi)
    !----------------------------------------------------
    ! Perform singular integrations for a subregion,
    ! delimited by [eta0,eta_final], of the boundary
    ! element be. variable eta0 is the parametric
    ! coordinate of the source point, sp.
    !
    ! Date: 23 July 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(inout) :: hi(:,:)
    type(boundary_elem),intent(in) :: be
    type(point_2d),intent(in) :: sp
    real(kind=8),intent(in) :: eta0,eta_final
    type(point_2d) :: fp !The current field point
    type(point_2d),allocatable :: begp(:) !Boundary element geometric points
    real(kind=8) :: d_eta !Interval size between two radial interpolation points
    real(kind=8) :: normal(2),normal_eta0(2) !Normal vector
    real(kind=8) :: ti(2) !Tangential vector
    real(kind=8) :: ri(3) ! r(1) = r, r(2) = r_1, r(3) = r_2
    real(kind=8) :: drdx(2) !r_i/r
    real(kind=8) :: r_dot_t !Internal product of drdx and ti
    real(kind=8) :: c0(2,2*be%nodes_num) !Coefficient c0 for radial expansion
    real(kind=8) :: ci(npow,2,2*be%nodes_num) !Coefficients c_i, i = 1,npow
    real(kind=8) :: func(2,2*be%nodes_num) !Current regular part of kernel
    real(kind=8),allocatable :: ni(:) !Interpolation function
    real(kind=8) :: eta !Current parametric coordinate
    real(kind=8),allocatable :: amat(:,:) !Matrix to evaluate expansion coefficients
    real(kind=8) :: pw
    real(kind=8) :: ts(2,2),tk(2,2)
    real(kind=8) :: nu
    integer :: i,j,k,psk
    !
    d_eta = (eta_final-eta0)/real(npow,8) !npow+1 radial interpolation points
    allocate(amat(npow,npow))
    !
    !-----------------------------------------------------------
    ! Store boundary element's geometry parametrization points:
    !-----------------------------------------------------------
    !
    call get_bound_elem_geom_points(be,begp)
    call get_plane_state_key(psk)
    call get_poisson_ratio(nu)
    !
    !---------------------------------------------------------
    ! Evaluating C0 component of the radial expansion series:
    !---------------------------------------------------------
    !
    allocate(ni(be%nodes_num))
    !
    tk = 0.0_8
    normal_eta0 = normal_1d(begp,eta0,be%nodes_num-1)
    drdx(1) = -normal_eta0(2)
    drdx(2) = normal_eta0(1)
    if (eta0.gt.eta_final) drdx = -1.0*drdx
    tk = t_ij_rp(drdx,normal_eta0,nu,psk)
    call intrp_func_1d(ni,eta0,be%lamb_1,be%lamb_2,be%nodes_num-1)
    do k = 1,size(ni)
        c0(1:2,2*k-1:2*k) = ni(k)*tk
    end do
    !
    !-----------------------------------------
    ! Evaluating the other coefficients, c_i:
    !-----------------------------------------
    !
    ts = 0.0_8
    eta = eta0
    do i = 1,npow
      eta = eta + d_eta
      if (i.eq.npow) eta = eta_final
      normal = normal_1d(begp,eta,be%nodes_num-1)
      ti(1) = -normal(2)
      ti(2) = normal(1)
      if (eta0.gt.eta_final) ti = -1.0*ti
      fp = cartesian_point_from_intrp_func_1d(begp,eta)
      ri = distance_2d(sp,fp)
      drdx = ri(2:3)/ri(1)
      ts = t_ij_rp(drdx,normal,nu,psk)
      call intrp_func_1d(ni,eta,be%lamb_1,be%lamb_2,be%nodes_num-1)
      do k = 1,size(ni)
        func(1:2,2*k-1:2*k) = ni(k)*ts
      end do
      r_dot_t = dot_product(drdx,ti)
      ci(i,:,:) = ((func/abs(r_dot_t)) - c0) / ri(1)
      amat(i,1) = 1.0_8
      do j = 2,npow
        amat(i,j) = amat(i,j-1)*ri(1)
      end do
    end do
    !
    call mat_inverse(amat)
    !
    do i = 1,size(ci(1,:,1))
      do j = 1,size(ci(1,1,:))
        ci(:,i,j) = matmul(amat,ci(:,i,j))
      end do
    end do
    !
    !-------------------------------
    ! Evaluating integrals results:
    !-------------------------------
    !
    pw = 0.0_8 ! 1.0 - 1.0 (beta)
    if (abs(pw).lt.tol) hi = hi + c0*log(ri(1))
    if (abs(pw).ge.tol) hi = hi + c0*ri(1)**pw / pw
    !
    !ci contribution:
    !
    do i = 1,npow
        pw = i ! 1.0 - 1.0 (beta) + i
        if (abs(pw).lt.tol) hi = hi + ci(i,:,:)*log(ri(1))
        if (abs(pw).ge.tol) hi = hi + ci(i,:,:)*ri(1)**pw / pw
    end do
    !
    deallocate(ni)
    !
    return
  end subroutine displ_subreg_sing_intgr
  !****************************************************
  !****************************************************
  subroutine traction_bie_sing_intgr(sp,be,lcp,normal_cp,hi,gi)
    !----------------------------------------------------
    ! Perform singular integrations for a given source
    ! point (sp) and a given boundary element (be).
    !
    ! lcp = location of the source point inside be.
    !
    ! Date: 06 July 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(inout) :: hi(:,:),gi(:,:)
    real(kind=8),intent(in) :: normal_cp(2)
    type(boundary_elem),intent(in) :: be
    type(point_2d),intent(in) :: sp
    integer,intent(in) :: lcp
    real(kind=8) :: eta0,eta_final
    !
    hi = 0.0
    gi = 0.0
    !
    !-----------------------------------------------------
    ! Storing the natural coordinate of the source point:
    !-----------------------------------------------------
    !
    if (lcp.eq.1) then
      eta0 = -be%lamb_1
    else if (lcp.eq.be%nodes_num) then
      eta0 = be%lamb_2
    else
      eta0 = 0.0_8 !Mid-node for quadratic elements
    end if
    !
    !-------------------------------------------
    ! Integrate region before the source point:
    !-------------------------------------------
    !
    if (eta0.gt.-1.0_8) then
      eta_final = -1.0_8
      call t_subreg_sing_intgr(sp,be,eta0,eta_final,normal_cp,hi)
      call u_subreg_sing_intgr(sp,be,eta0,eta_final,normal_cp,gi)
    end if
    !
    !------------------------------------------
    ! Integrate region after the source point:
    !------------------------------------------
    !
    if (eta0.lt.1.0_8) then
      eta_final = 1.0_8
      call t_subreg_sing_intgr(sp,be,eta0,eta_final,normal_cp,hi)
      call u_subreg_sing_intgr(sp,be,eta0,eta_final,normal_cp,gi)
    end if
    !
    return
  end subroutine traction_bie_sing_intgr
  !****************************************************
  !****************************************************
  subroutine t_subreg_sing_intgr(sp,be,eta0,eta_final,normal_cp,hi)
    !----------------------------------------------------
    ! Perform singular integrations for a subregion,
    ! delimited by [eta0,eta_final], of the boundary
    ! element be. variable eta0 is the parametric
    ! coordinate of the source point, sp.
    !
    ! Date: 06 July 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(inout) :: hi(:,:)
    type(boundary_elem),intent(in) :: be
    type(point_2d),intent(in) :: sp
    real(kind=8),intent(in) :: eta0,eta_final,normal_cp(2)
    type(point_2d) :: fp !The current field point
    type(point_2d),allocatable :: begp(:) !Boundary element geometric points
    real(kind=8) :: d_eta !Interval size between two radial interpolation points
    real(kind=8) :: normal(2),normal_eta0(2) !Normal vector
    real(kind=8) :: ti(2) !Tangential vector
    real(kind=8) :: ri(3) ! r(1) = r, r(2) = r_1, r(3) = r_2
    real(kind=8) :: drdx(2) !r_i/r
    real(kind=8) :: r_dot_t !Internal product of drdx and ti
    real(kind=8) :: c0(2,2*be%nodes_num) !Coefficient c0 for radial expansion
    real(kind=8) :: ci(npow,2,2*be%nodes_num) !Coefficients c_i, i = 1,npow
    real(kind=8) :: func(2,2*be%nodes_num) !Current regular part of kernel
    real(kind=8),allocatable :: ni(:) !Interpolation function
    real(kind=8) :: eta !Current parametric coordinate
    real(kind=8),allocatable :: amat(:,:) !Matrix to evaluate expansion coefficients
    real(kind=8) :: pw
    real(kind=8) :: ts(2,2),tk(2,2)
    real(kind=8) :: mu,nu
    integer :: i,j,k,psk
    !
    d_eta = (eta_final-eta0)/real(npow,8) !npow+1 radial interpolation points
    allocate(amat(npow,npow))
    !
    !-----------------------------------------------------------
    ! Store boundary element's geometry parametrization points:
    !-----------------------------------------------------------
    !
    call get_bound_elem_geom_points(be,begp)
    call get_plane_state_key(psk)
    call get_poisson_ratio(nu)
    call get_shear_modulus(mu)
    !
    !---------------------------------------------------------
    ! Evaluating C0 component of the radial expansion series:
    !---------------------------------------------------------
    !
    allocate(ni(be%nodes_num))
    !
    tk = 0.0_8
    normal_eta0 = normal_1d(begp,eta0,be%nodes_num-1)
    drdx(1) = -normal_eta0(2)
    drdx(2) = normal_eta0(1)
    if (eta0.gt.eta_final) drdx = -1.0*drdx
    tk = normal_t_ijk_rp(drdx,normal_eta0,normal_cp,mu,nu,psk)
    call intrp_func_1d(ni,eta0,be%lamb_1,be%lamb_2,be%nodes_num-1)
    do k = 1,size(ni)
        c0(1:2,2*k-1:2*k) = ni(k)*tk
    end do
    !
    !-----------------------------------------
    ! Evaluating the other coefficients, c_i:
    !-----------------------------------------
    !
    ts = 0.0_8
    eta = eta0
    do i = 1,npow
      eta = eta + d_eta
      if (i.eq.npow) eta = eta_final
      normal = normal_1d(begp,eta,be%nodes_num-1)
      ti(1) = -normal(2)
      ti(2) = normal(1)
      if (eta0.gt.eta_final) ti = -1.0*ti
      fp = cartesian_point_from_intrp_func_1d(begp,eta)
      ri = distance_2d(sp,fp)
      drdx = ri(2:3)/ri(1)
      ts = normal_t_ijk_rp(drdx,normal,normal_cp,mu,nu,psk)
      call intrp_func_1d(ni,eta,be%lamb_1,be%lamb_2,be%nodes_num-1)
      do k = 1,size(ni)
        func(1:2,2*k-1:2*k) = ni(k)*ts
      end do
      r_dot_t = dot_product(drdx,ti)
      ci(i,:,:) = ((func/abs(r_dot_t)) - c0) / ri(1)
      amat(i,1) = 1.0_8
      do j = 2,npow
        amat(i,j) = amat(i,j-1)*ri(1)
      end do
    end do
    !
    call mat_inverse(amat)
    !
    do i = 1,size(ci(1,:,1))
      do j = 1,size(ci(1,1,:))
        ci(:,i,j) = matmul(amat,ci(:,i,j))
      end do
    end do
    !
    !-------------------------------
    ! Evaluating integrals results:
    !-------------------------------
    !
    pw = - 1.0_8 ! 1.0 - 2.0 (beta)
    if (abs(pw).lt.tol) hi = hi + c0*log(ri(1))
    if (abs(pw).ge.tol) hi = hi + c0*ri(1)**pw / pw
    !
    !ci contribution:
    !
    do i = 1,npow
        pw = -1.0_8 + i ! 1.0 - 2.0 (beta) + i
        if (abs(pw).lt.tol) hi = hi + ci(i,:,:)*log(ri(1))
        if (abs(pw).ge.tol) hi = hi + ci(i,:,:)*ri(1)**pw / pw
    end do
    !
    deallocate(ni)
    !
    return
  end subroutine t_subreg_sing_intgr
  !****************************************************
  !****************************************************
  subroutine u_subreg_sing_intgr(sp,be,eta0,eta_final,normal_cp,gi)
    !----------------------------------------------------
    ! Perform singular integrations for a subregion,
    ! delimited by [eta0,eta_final], of the boundary
    ! element be. variable eta0 is the parametric
    ! coordinate of the source point, sp.
    !
    ! Date: 06 July 2022
    ! Last update: -
    !----------------------------------------------------
    real(kind=8),intent(inout) :: gi(:,:)
    type(boundary_elem),intent(in) :: be
    type(point_2d),intent(in) :: sp
    real(kind=8),intent(in) :: eta0,eta_final,normal_cp(2)
    type(point_2d) :: fp !The current field point
    type(point_2d),allocatable :: begp(:) !Boundary element geometric points
    real(kind=8) :: d_eta !Interval size between two radial interpolation points
    real(kind=8) :: normal(2), normal_eta0(2) !Normal vector
    real(kind=8) :: ti(2) !Tangential vector
    real(kind=8) :: ri(3) ! r(1) = r, r(2) = r_1, r(3) = r_2
    real(kind=8) :: drdx(2) !r_i/r
    real(kind=8) :: r_dot_t !Internal product of drdx and ti
    real(kind=8),allocatable :: c0(:,:) !Coefficient c0 for radial expansion
    real(kind=8),allocatable :: ci(:,:,:) !Coefficients c_i, i = 1,npow
    real(kind=8),allocatable :: func(:,:) !Current regular part of kernel
    real(kind=8),allocatable :: ni(:) !Interpolation function
    real(kind=8) :: eta !Current parametric coordinate
    real(kind=8) ,allocatable:: amat(:,:) !Matrix to evaluate expansion coefficients
    real(kind=8) :: pw
    real(kind=8) :: ui(2,2),us(2,2)
    real(kind=8) :: mu,nu
    integer :: i,j,k,psk
    !
    d_eta = (eta_final-eta0)/real(npow,8) !npow+1 radial interpolation points
    allocate(amat(npow,npow))
    !
    !-----------------------------------------------------------
    ! Store boundary element's geometry parametrization points:
    !-----------------------------------------------------------
    !
    call get_bound_elem_geom_points(be,begp)
    call get_plane_state_key(psk)
    call get_poisson_ratio(nu)
    call get_shear_modulus(mu)
    !
    !---------------------------------------------------------
    ! Evaluating C0 component of the radial expansion series:
    !---------------------------------------------------------
    !
    allocate(c0(size(gi(:,1)),size(gi(1,:))))
    allocate(ni(be%nodes_num))
    !
    us = 0.0_8
    normal_eta0 = normal_1d(begp,eta0,be%nodes_num-1)
    drdx(1) = -normal_eta0(2)
    drdx(2) = normal_eta0(1)
    if (eta0.gt.eta_final) drdx = -1.0*drdx
    us = normal_u_ijk_rp(drdx,normal_cp,nu,psk)
    call intrp_func_1d(ni,eta0,be%lamb_1,be%lamb_2,be%nodes_num-1)
    do k = 1,size(ni)
        c0(1:2,2*k-1:2*k) = ni(k)*us
    end do
    !
    !-----------------------------------------
    ! Evaluating the other coefficients, c_i:
    !-----------------------------------------
    !
    allocate(ci(npow,size(gi(:,1)),size(gi(1,:))))
    allocate(func(size(gi(:,1)),size(gi(1,:))))
    !
    ui = 0.0_8
    eta = eta0
    do i = 1,npow
      eta = eta + d_eta
      if (i.eq.npow) eta = eta_final
      normal = normal_1d(begp,eta,be%nodes_num-1)
      ti(1) = -normal(2)
      ti(2) = normal(1)
      if (eta0.gt.eta_final) ti = -1.0*ti
      fp = cartesian_point_from_intrp_func_1d(begp,eta)
      ri = distance_2d(sp,fp)
      drdx = ri(2:3)/ri(1)
      ui = normal_u_ijk_rp(drdx,normal_cp,nu,psk)
      call intrp_func_1d(ni,eta,be%lamb_1,be%lamb_2,be%nodes_num-1)
      do k = 1,size(ni)
        func(1:2,2*k-1:2*k) = ni(k)*ui
      end do
      r_dot_t = dot_product(drdx,ti)
      ci(i,:,:) = ((func/abs(r_dot_t)) - c0) / ri(1)
      amat(i,1) = 1.0_8
      do j = 2,npow
        amat(i,j) = amat(i,j-1)*ri(1)
      end do
      end do
    !
    call mat_inverse(amat)
    !
    do i = 1,size(ci(1,:,1))
      do j = 1,size(ci(1,1,:))
        ci(:,i,j) = matmul(amat,ci(:,i,j))
      end do
    end do
    !
    !-------------------------------
    ! Evaluating integrals results:
    !-------------------------------
    !
    !c0 contribution:
    !
    pw = 0.0_8 ! 1.0 - 1.0 (beta)
    if (abs(pw).lt.tol) gi = gi + c0*log(ri(1))
    if (abs(pw).ge.tol) gi = gi + c0*ri(1)**pw / pw
    !
    !ci contribution:
    !
    do i = 1,npow
        pw = i ! 1.0 - 1.0(beta) + i
        if (abs(pw).lt.tol) gi = gi + ci(i,:,:)*log(ri(1))
        if (abs(pw).ge.tol) gi = gi + ci(i,:,:)*ri(1)**pw / pw
    end do
    !
    deallocate(c0)
    deallocate(ci)
    deallocate(func)
    deallocate(ni)
    !
    return
  end subroutine u_subreg_sing_intgr
  !****************************************************
  !****************************************************
  subroutine displ_grad_sing_intgr(sp,be,lcp,hi,gi)
    !----------------------------------------------------
    ! Perform singular integrations for a given source
    ! point (sp) and a given boundary element (be).
    !
    ! lcp = location of the source point inside be.
    !
    ! Date: 22 July 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(inout) :: hi(:,:),gi(:,:)
    type(boundary_elem),intent(in) :: be
    type(point_2d),intent(in) :: sp
    integer,intent(in) :: lcp
    real(kind=8) :: eta0,eta_final
    !
    hi = 0.0
    gi = 0.0
    !
    !-----------------------------------------------------
    ! Storing the natural coordinate of the source point:
    !-----------------------------------------------------
    !
    if (lcp.eq.1) then
      eta0 = -be%lamb_1
    else if (lcp.eq.be%nodes_num) then
      eta0 = be%lamb_2
    else
      eta0 = 0.0_8 !Mid-node for quadratic elements
    end if
    !
    !-------------------------------------------
    ! Integrate region before the source point:
    !-------------------------------------------
    !
    if (eta0.gt.-1.0_8) then
      eta_final = -1.0_8
      call grad_t_subreg_sing_intgr(sp,be,eta0,eta_final,hi)
      call grad_u_subreg_sing_intgr(sp,be,eta0,eta_final,gi)
    end if
    !
    !------------------------------------------
    ! Integrate region after the source point:
    !------------------------------------------
    !
    if (eta0.lt.1.0_8) then
      eta_final = 1.0_8
      call grad_t_subreg_sing_intgr(sp,be,eta0,eta_final,hi)
      call grad_u_subreg_sing_intgr(sp,be,eta0,eta_final,gi)
    end if
    !
    return
  end subroutine displ_grad_sing_intgr
  !****************************************************
  !****************************************************
  subroutine grad_t_subreg_sing_intgr(sp,be,eta0,eta_final,hi)
    !----------------------------------------------------
    ! Perform singular integrations for a subregion,
    ! delimited by [eta0,eta_final], of the boundary
    ! element be. variable eta0 is the parametric
    ! coordinate of the source point, sp.
    !
    ! Date: 22 July 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(inout) :: hi(:,:)
    type(boundary_elem),intent(in) :: be
    type(point_2d),intent(in) :: sp
    real(kind=8),intent(in) :: eta0,eta_final
    type(point_2d) :: fp !The current field point
    type(point_2d),allocatable :: begp(:) !Boundary element geometric points
    real(kind=8) :: d_eta !Interval size between two radial interpolation points
    real(kind=8) :: normal(2),normal_eta0(2) !Normal vector
    real(kind=8) :: ti(2) !Tangential vector
    real(kind=8) :: ri(3) ! r(1) = r, r(2) = r_1, r(3) = r_2
    real(kind=8) :: drdx(2) !r_i/r
    real(kind=8) :: r_dot_t !Internal product of drdx and ti
    real(kind=8) :: c0(4,2*be%nodes_num) !Coefficient c0 for radial expansion
    real(kind=8) :: ci(npow,4,2*be%nodes_num) !Coefficients c_i, i = 1,npow
    real(kind=8) :: func(4,2*be%nodes_num) !Current regular part of kernel
    real(kind=8),allocatable :: ni(:) !Interpolation function
    real(kind=8) :: eta !Current parametric coordinate
    real(kind=8),allocatable :: amat(:,:) !Matrix to evaluate expansion coefficients
    real(kind=8) :: pw
    real(kind=8) :: ts(4,2),tk(4,2)
    real(kind=8) :: mu,nu
    integer :: i,j,k,psk
    !
    d_eta = (eta_final-eta0)/real(npow,8) !npow+1 radial interpolation points
    allocate(amat(npow,npow))
    !
    !-----------------------------------------------------------
    ! Store boundary element's geometry parametrization points:
    !-----------------------------------------------------------
    !
    call get_bound_elem_geom_points(be,begp)
    call get_plane_state_key(psk)
    call get_poisson_ratio(nu)
    call get_shear_modulus(mu)
    !
    !---------------------------------------------------------
    ! Evaluating C0 component of the radial expansion series:
    !---------------------------------------------------------
    !
    allocate(ni(be%nodes_num))
    !
    tk = 0.0_8
    normal_eta0 = normal_1d(begp,eta0,be%nodes_num-1)
    drdx(1) = -normal_eta0(2)
    drdx(2) = normal_eta0(1)
    if (eta0.gt.eta_final) drdx = -1.0*drdx
    tk = t_ij_k_rp(drdx,nu,normal_eta0,psk)
    call intrp_func_1d(ni,eta0,be%lamb_1,be%lamb_2,be%nodes_num-1)
    do k = 1,size(ni)
        c0(1:4,2*k-1:2*k) = ni(k)*tk
    end do
    !
    !-----------------------------------------
    ! Evaluating the other coefficients, c_i:
    !-----------------------------------------
    !
    ts = 0.0_8
    eta = eta0
    do i = 1,npow
      eta = eta + d_eta
      if (i.eq.npow) eta = eta_final
      normal = normal_1d(begp,eta,be%nodes_num-1)
      ti(1) = -normal(2)
      ti(2) = normal(1)
      if (eta0.gt.eta_final) ti = -1.0*ti
      fp = cartesian_point_from_intrp_func_1d(begp,eta)
      ri = distance_2d(sp,fp)
      drdx = ri(2:3)/ri(1)
      ts = t_ij_k_rp(drdx,nu,normal,psk)
      call intrp_func_1d(ni,eta,be%lamb_1,be%lamb_2,be%nodes_num-1)
      do k = 1,size(ni)
        func(1:4,2*k-1:2*k) = ni(k)*ts
      end do
      r_dot_t = dot_product(drdx,ti)
      ci(i,:,:) = ((func/abs(r_dot_t)) - c0) / ri(1)
      amat(i,1) = 1.0_8
      do j = 2,npow
        amat(i,j) = amat(i,j-1)*ri(1)
      end do
    end do
    !
    call mat_inverse(amat)
    !
    do i = 1,size(ci(1,:,1))
      do j = 1,size(ci(1,1,:))
        ci(:,i,j) = matmul(amat,ci(:,i,j))
      end do
    end do
    !
    !-------------------------------
    ! Evaluating integrals results:
    !-------------------------------
    !
    pw = - 1.0_8 ! 1.0 - 2.0 (beta)
    if (abs(pw).lt.tol) hi = hi + c0*log(ri(1))
    if (abs(pw).ge.tol) hi = hi + c0*ri(1)**pw / pw
    !
    !ci contribution:
    !
    do i = 1,npow
        pw = -1.0_8 + i ! 1.0 - 2.0 (beta) + i
        if (abs(pw).lt.tol) hi = hi + ci(i,:,:)*log(ri(1))
        if (abs(pw).ge.tol) hi = hi + ci(i,:,:)*ri(1)**pw / pw
    end do
    !
    deallocate(ni)
    !
    return
  end subroutine grad_t_subreg_sing_intgr
  !****************************************************
  !****************************************************
  subroutine grad_u_subreg_sing_intgr(sp,be,eta0,eta_final,gi)
    !----------------------------------------------------
    ! Perform singular integrations for a subregion,
    ! delimited by [eta0,eta_final], of the boundary
    ! element be. variable eta0 is the parametric
    ! coordinate of the source point, sp.
    !
    ! Date: 22 July 2022
    ! Last update: -
    !----------------------------------------------------
    real(kind=8),intent(inout) :: gi(:,:)
    type(boundary_elem),intent(in) :: be
    type(point_2d),intent(in) :: sp
    real(kind=8),intent(in) :: eta0,eta_final
    type(point_2d) :: fp !The current field point
    type(point_2d),allocatable :: begp(:) !Boundary element geometric points
    real(kind=8) :: d_eta !Interval size between two radial interpolation points
    real(kind=8) :: normal(2), normal_eta0(2) !Normal vector
    real(kind=8) :: ti(2) !Tangential vector
    real(kind=8) :: ri(3) ! r(1) = r, r(2) = r_1, r(3) = r_2
    real(kind=8) :: drdx(2) !r_i/r
    real(kind=8) :: r_dot_t !Internal product of drdx and ti
    real(kind=8),allocatable :: c0(:,:) !Coefficient c0 for radial expansion
    real(kind=8),allocatable :: ci(:,:,:) !Coefficients c_i, i = 1,npow
    real(kind=8),allocatable :: func(:,:) !Current regular part of kernel
    real(kind=8),allocatable :: ni(:) !Interpolation function
    real(kind=8) :: eta !Current parametric coordinate
    real(kind=8) ,allocatable:: amat(:,:) !Matrix to evaluate expansion coefficients
    real(kind=8) :: pw
    real(kind=8) :: ui(4,2),us(4,2)
    real(kind=8) :: mu,nu
    integer :: i,j,k,psk
    !
    d_eta = (eta_final-eta0)/real(npow,8) !npow+1 radial interpolation points
    allocate(amat(npow,npow))
    !
    !-----------------------------------------------------------
    ! Store boundary element's geometry parametrization points:
    !-----------------------------------------------------------
    !
    call get_bound_elem_geom_points(be,begp)
    call get_plane_state_key(psk)
    call get_poisson_ratio(nu)
    call get_shear_modulus(mu)
    !
    !---------------------------------------------------------
    ! Evaluating C0 component of the radial expansion series:
    !---------------------------------------------------------
    !
    allocate(c0(size(gi(:,1)),size(gi(1,:))))
    allocate(ni(be%nodes_num))
    !
    us = 0.0_8
    normal_eta0 = normal_1d(begp,eta0,be%nodes_num-1)
    drdx(1) = -normal_eta0(2)
    drdx(2) = normal_eta0(1)
    if (eta0.gt.eta_final) drdx = -1.0*drdx
    us = u_ij_k_rp(drdx,mu,nu,psk)
    call intrp_func_1d(ni,eta0,be%lamb_1,be%lamb_2,be%nodes_num-1)
    do k = 1,size(ni)
        c0(1:4,2*k-1:2*k) = ni(k)*us
    end do
    !
    !-----------------------------------------
    ! Evaluating the other coefficients, c_i:
    !-----------------------------------------
    !
    allocate(ci(npow,size(gi(:,1)),size(gi(1,:))))
    allocate(func(size(gi(:,1)),size(gi(1,:))))
    !
    ui = 0.0_8
    eta = eta0
    do i = 1,npow
      eta = eta + d_eta
      if (i.eq.npow) eta = eta_final
      normal = normal_1d(begp,eta,be%nodes_num-1)
      ti(1) = -normal(2)
      ti(2) = normal(1)
      if (eta0.gt.eta_final) ti = -1.0*ti
      fp = cartesian_point_from_intrp_func_1d(begp,eta)
      ri = distance_2d(sp,fp)
      drdx = ri(2:3)/ri(1)
      ui = u_ij_k_rp(drdx,mu,nu,psk)
      call intrp_func_1d(ni,eta,be%lamb_1,be%lamb_2,be%nodes_num-1)
      do k = 1,size(ni)
        func(1:4,2*k-1:2*k) = ni(k)*ui
      end do
      r_dot_t = dot_product(drdx,ti)
      ci(i,:,:) = ((func/abs(r_dot_t)) - c0) / ri(1)
      amat(i,1) = 1.0_8
      do j = 2,npow
        amat(i,j) = amat(i,j-1)*ri(1)
      end do
      end do
    !
    call mat_inverse(amat)
    !
    do i = 1,size(ci(1,:,1))
      do j = 1,size(ci(1,1,:))
        ci(:,i,j) = matmul(amat,ci(:,i,j))
      end do
    end do
    !
    !-------------------------------
    ! Evaluating integrals results:
    !-------------------------------
    !
    !c0 contribution:
    !
    pw = 0.0_8 ! 1.0 - 1.0 (beta)
    if (abs(pw).lt.tol) gi = gi + c0*log(ri(1))
    if (abs(pw).ge.tol) gi = gi + c0*ri(1)**pw / pw
    !
    !ci contribution:
    !
    do i = 1,npow
        pw = i ! 1.0 - 1.0(beta) + i
        if (abs(pw).lt.tol) gi = gi + ci(i,:,:)*log(ri(1))
        if (abs(pw).ge.tol) gi = gi + ci(i,:,:)*ri(1)**pw / pw
    end do
    !
    deallocate(c0)
    deallocate(ci)
    deallocate(func)
    deallocate(ni)
    !
    return
  end subroutine grad_u_subreg_sing_intgr
  !****************************************************
  !****************************************************
  subroutine stress_sing_intgr(sp,be,lcp,hi,gi)
    !----------------------------------------------------
    ! Perform singular integrations for a given source
    ! point (sp) and a given boundary element (be).
    !
    ! lcp = location of the source point inside be.
    !
    ! Date: 20 July 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(inout) :: hi(:,:),gi(:,:)
    type(boundary_elem),intent(in) :: be
    type(point_2d),intent(in) :: sp
    integer,intent(in) :: lcp
    real(kind=8) :: eta0,eta_final
    !
    hi = 0.0
    gi = 0.0
    !
    !-----------------------------------------------------
    ! Storing the natural coordinate of the source point:
    !-----------------------------------------------------
    !
    if (lcp.eq.1) then
      eta0 = -be%lamb_1
    else if (lcp.eq.be%nodes_num) then
      eta0 = be%lamb_2
    else
      eta0 = 0.0_8 !Mid-node for quadratic elements
    end if
    !
    !-------------------------------------------
    ! Integrate region before the source point:
    !-------------------------------------------
    !
    if (eta0.gt.-1.0_8) then
      eta_final = -1.0_8
      call stress_t_subreg_sing_intgr(sp,be,eta0,eta_final,hi)
      call stress_u_subreg_sing_intgr(sp,be,eta0,eta_final,gi)
    end if
    !
    !------------------------------------------
    ! Integrate region after the source point:
    !------------------------------------------
    !
    if (eta0.lt.1.0_8) then
      eta_final = 1.0_8
      call stress_t_subreg_sing_intgr(sp,be,eta0,eta_final,hi)
      call stress_u_subreg_sing_intgr(sp,be,eta0,eta_final,gi)
    end if
    !
    return
  end subroutine stress_sing_intgr
  !****************************************************
  !****************************************************
  subroutine stress_t_subreg_sing_intgr(sp,be,eta0,eta_final,hi)
    !----------------------------------------------------
    ! Perform singular integrations for a subregion,
    ! delimited by [eta0,eta_final], of the boundary
    ! element be. variable eta0 is the parametric
    ! coordinate of the source point, sp.
    !
    ! Date: 20 July 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(inout) :: hi(:,:)
    type(boundary_elem),intent(in) :: be
    type(point_2d),intent(in) :: sp
    real(kind=8),intent(in) :: eta0,eta_final
    type(point_2d) :: fp !The current field point
    type(point_2d),allocatable :: begp(:) !Boundary element geometric points
    real(kind=8) :: d_eta !Interval size between two radial interpolation points
    real(kind=8) :: normal(2),normal_eta0(2) !Normal vector
    real(kind=8) :: ti(2) !Tangential vector
    real(kind=8) :: ri(3) ! r(1) = r, r(2) = r_1, r(3) = r_2
    real(kind=8) :: drdx(2) !r_i/r
    real(kind=8) :: r_dot_t !Internal product of drdx and ti
    real(kind=8) :: c0(3,2*be%nodes_num) !Coefficient c0 for radial expansion
    real(kind=8) :: ci(npow,3,2*be%nodes_num) !Coefficients c_i, i = 1,npow
    real(kind=8) :: func(3,2*be%nodes_num) !Current regular part of kernel
    real(kind=8),allocatable :: ni(:) !Interpolation function
    real(kind=8) :: eta !Current parametric coordinate
    real(kind=8),allocatable :: amat(:,:) !Matrix to evaluate expansion coefficients
    real(kind=8) :: pw
    real(kind=8) :: ts(3,2),tk(3,2)
    real(kind=8) :: mu,nu
    integer :: i,j,k,psk
    !
    d_eta = (eta_final-eta0)/real(npow,8) !npow+1 radial interpolation points
    allocate(amat(npow,npow))
    !
    !-----------------------------------------------------------
    ! Store boundary element's geometry parametrization points:
    !-----------------------------------------------------------
    !
    call get_bound_elem_geom_points(be,begp)
    call get_plane_state_key(psk)
    call get_poisson_ratio(nu)
    call get_shear_modulus(mu)
    !
    !---------------------------------------------------------
    ! Evaluating C0 component of the radial expansion series:
    !---------------------------------------------------------
    !
    allocate(ni(be%nodes_num))
    !
    tk = 0.0_8
    normal_eta0 = normal_1d(begp,eta0,be%nodes_num-1)
    drdx(1) = -normal_eta0(2)
    drdx(2) = normal_eta0(1)
    if (eta0.gt.eta_final) drdx = -1.0*drdx
    tk = t_ijk_rp(drdx,normal_eta0,mu,nu,psk)
    call intrp_func_1d(ni,eta0,be%lamb_1,be%lamb_2,be%nodes_num-1)
    do k = 1,size(ni)
        c0(1:3,2*k-1:2*k) = ni(k)*tk
    end do
    !
    !-----------------------------------------
    ! Evaluating the other coefficients, c_i:
    !-----------------------------------------
    !
    ts = 0.0_8
    eta = eta0
    do i = 1,npow
      eta = eta + d_eta
      if (i.eq.npow) eta = eta_final
      normal = normal_1d(begp,eta,be%nodes_num-1)
      ti(1) = -normal(2)
      ti(2) = normal(1)
      if (eta0.gt.eta_final) ti = -1.0*ti
      fp = cartesian_point_from_intrp_func_1d(begp,eta)
      ri = distance_2d(sp,fp)
      drdx = ri(2:3)/ri(1)
      ts = t_ijk_rp(drdx,normal,mu,nu,psk)
      call intrp_func_1d(ni,eta,be%lamb_1,be%lamb_2,be%nodes_num-1)
      do k = 1,size(ni)
        func(1:3,2*k-1:2*k) = ni(k)*ts
      end do
      r_dot_t = dot_product(drdx,ti)
      ci(i,:,:) = ((func/abs(r_dot_t)) - c0) / ri(1)
      amat(i,1) = 1.0_8
      do j = 2,npow
        amat(i,j) = amat(i,j-1)*ri(1)
      end do
    end do
    !
    call mat_inverse(amat)
    !
    do i = 1,size(ci(1,:,1))
      do j = 1,size(ci(1,1,:))
        ci(:,i,j) = matmul(amat,ci(:,i,j))
      end do
    end do
    !
    !-------------------------------
    ! Evaluating integrals results:
    !-------------------------------
    !
    pw = - 1.0_8 ! 1.0 - 2.0 (beta)
    if (abs(pw).lt.tol) hi = hi + c0*log(ri(1))
    if (abs(pw).ge.tol) hi = hi + c0*ri(1)**pw / pw
    !
    !ci contribution:
    !
    do i = 1,npow
        pw = -1.0_8 + i ! 1.0 - 2.0 (beta) + i
        if (abs(pw).lt.tol) hi = hi + ci(i,:,:)*log(ri(1))
        if (abs(pw).ge.tol) hi = hi + ci(i,:,:)*ri(1)**pw / pw
    end do
    !
    deallocate(ni)
    !
    return
  end subroutine stress_t_subreg_sing_intgr
   !****************************************************
  !****************************************************
  subroutine stress_u_subreg_sing_intgr(sp,be,eta0,eta_final,gi)
    !----------------------------------------------------
    ! Perform singular integrations for a subregion,
    ! delimited by [eta0,eta_final], of the boundary
    ! element be. variable eta0 is the parametric
    ! coordinate of the source point, sp.
    !
    ! Date: 20 July 2022
    ! Last update: -
    !----------------------------------------------------
    real(kind=8),intent(inout) :: gi(:,:)
    type(boundary_elem),intent(in) :: be
    type(point_2d),intent(in) :: sp
    real(kind=8),intent(in) :: eta0,eta_final
    type(point_2d) :: fp !The current field point
    type(point_2d),allocatable :: begp(:) !Boundary element geometric points
    real(kind=8) :: d_eta !Interval size between two radial interpolation points
    real(kind=8) :: normal(2), normal_eta0(2) !Normal vector
    real(kind=8) :: ti(2) !Tangential vector
    real(kind=8) :: ri(3) ! r(1) = r, r(2) = r_1, r(3) = r_2
    real(kind=8) :: drdx(2) !r_i/r
    real(kind=8) :: r_dot_t !Internal product of drdx and ti
    real(kind=8),allocatable :: c0(:,:) !Coefficient c0 for radial expansion
    real(kind=8),allocatable :: ci(:,:,:) !Coefficients c_i, i = 1,npow
    real(kind=8),allocatable :: func(:,:) !Current regular part of kernel
    real(kind=8),allocatable :: ni(:) !Interpolation function
    real(kind=8) :: eta !Current parametric coordinate
    real(kind=8) ,allocatable:: amat(:,:) !Matrix to evaluate expansion coefficients
    real(kind=8) :: pw
    real(kind=8) :: ui(3,2),us(3,2)
    real(kind=8) :: mu,nu
    integer :: i,j,k,psk
    !
    d_eta = (eta_final-eta0)/real(npow,8) !npow+1 radial interpolation points
    allocate(amat(npow,npow))
    !
    !-----------------------------------------------------------
    ! Store boundary element's geometry parametrization points:
    !-----------------------------------------------------------
    !
    call get_bound_elem_geom_points(be,begp)
    call get_plane_state_key(psk)
    call get_poisson_ratio(nu)
    call get_shear_modulus(mu)
    !
    !---------------------------------------------------------
    ! Evaluating C0 component of the radial expansion series:
    !---------------------------------------------------------
    !
    allocate(c0(size(gi(:,1)),size(gi(1,:))))
    allocate(ni(be%nodes_num))
    !
    us = 0.0_8
    normal_eta0 = normal_1d(begp,eta0,be%nodes_num-1)
    drdx(1) = -normal_eta0(2)
    drdx(2) = normal_eta0(1)
    if (eta0.gt.eta_final) drdx = -1.0*drdx
    us = u_ijk_rp(drdx,nu,psk)
    call intrp_func_1d(ni,eta0,be%lamb_1,be%lamb_2,be%nodes_num-1)
    do k = 1,size(ni)
        c0(1:3,2*k-1:2*k) = ni(k)*us
    end do
    !
    !-----------------------------------------
    ! Evaluating the other coefficients, c_i:
    !-----------------------------------------
    !
    allocate(ci(npow,size(gi(:,1)),size(gi(1,:))))
    allocate(func(size(gi(:,1)),size(gi(1,:))))
    !
    ui = 0.0_8
    eta = eta0
    do i = 1,npow
      eta = eta + d_eta
      if (i.eq.npow) eta = eta_final
      normal = normal_1d(begp,eta,be%nodes_num-1)
      ti(1) = -normal(2)
      ti(2) = normal(1)
      if (eta0.gt.eta_final) ti = -1.0*ti
      fp = cartesian_point_from_intrp_func_1d(begp,eta)
      ri = distance_2d(sp,fp)
      drdx = ri(2:3)/ri(1)
      ui = u_ijk_rp(drdx,nu,psk)
      call intrp_func_1d(ni,eta,be%lamb_1,be%lamb_2,be%nodes_num-1)
      do k = 1,size(ni)
        func(1:3,2*k-1:2*k) = ni(k)*ui
      end do
      r_dot_t = dot_product(drdx,ti)
      ci(i,:,:) = ((func/abs(r_dot_t)) - c0) / ri(1)
      amat(i,1) = 1.0_8
      do j = 2,npow
        amat(i,j) = amat(i,j-1)*ri(1)
      end do
    end do
    !
    call mat_inverse(amat)
    !
    do i = 1,size(ci(1,:,1))
      do j = 1,size(ci(1,1,:))
        ci(:,i,j) = matmul(amat,ci(:,i,j))
      end do
    end do
    !
    !-------------------------------
    ! Evaluating integrals results:
    !-------------------------------
    !
    !c0 contribution:
    !
    pw = 0.0_8 ! 1.0 - 1.0 (beta)
    if (abs(pw).lt.tol) gi = gi + c0*log(ri(1))
    if (abs(pw).ge.tol) gi = gi + c0*ri(1)**pw / pw
    !
    !ci contribution:
    !
    do i = 1,npow
        pw = i ! 1.0 - 1.0(beta) + i
        if (abs(pw).lt.tol) gi = gi + ci(i,:,:)*log(ri(1))
        if (abs(pw).ge.tol) gi = gi + ci(i,:,:)*ri(1)**pw / pw
    end do
    !
    deallocate(c0)
    deallocate(ci)
    deallocate(func)
    deallocate(ni)
    !
    return
  end subroutine stress_u_subreg_sing_intgr
  !****************************************************
  !****************************************************
    subroutine set_npow(n)
    !----------------------------------------------------
    ! Set the radial expansion polynomial order.
    !
    ! Date: 20 July 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: n
    !
    npow = n
    !
    return
  end subroutine set_npow
  !****************************************************
  !****************************************************
  subroutine get_npow(n)
    !----------------------------------------------------
    ! Get the radial expansion polynomial order.
    !
    ! Date: 20 July 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(out) :: n
    !
    n = npow
    !
    return
  end subroutine get_npow
  !****************************************************
  !****************************************************
end module fracture_2d_rad_exp_intgr
