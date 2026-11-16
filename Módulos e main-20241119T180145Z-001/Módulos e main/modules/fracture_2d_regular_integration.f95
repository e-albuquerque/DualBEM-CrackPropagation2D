module fracture_2d_regular_integration
!----------------------------------------------------
! A module for regular (or near singular)
! integrations for two-dimensional isotropic
! fracture problems.
!
! Date: July 2022
! Last update: -
!----------------------------------------------------
use fracture_2d_model
use error_messages
use geometry
use kelvin_fund_sol_2d
use interpolation_functions
use numerical_integration_points
!
implicit none
!
contains
  subroutine displ_bie_reg_integr(psk,mu,nu,cp,be,h,g)
    !----------------------------------------------------
    ! Perform numerical integration for displacement
    ! boundary integral regular terms by Gauss-Legendre
    ! quadrature and store results in matrices h and g.
    !
    ! cp = the source point
    ! be = the boundary element
    ! h = Integral of [N_1*[T] ... N_n*[T]], n = order of be
    ! g = Integral of [N_1*[U] ... N_n*[U]], n = order of be
    ! N_k = interpolation functions
    ! [T] = 2x2 matrix with t* Kelvin fundamental solution
    ! [U] = 2x2 matrix with u* Kelvin fundamental solution
    !
    ! Date: 27 June 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: psk
    real(kind=8),intent(in) :: mu,nu
    type(point_2d),intent(in) :: cp
    type(boundary_elem),intent(in) :: be
    real(kind=8),intent(out) :: h(:,:),g(:,:)
    real(kind=8),allocatable :: eta_b(:),ome(:)
    real(kind=8),allocatable :: ni(:) !Interpolation functions values
    type(point_2d),allocatable :: gi_pts(:) !Geometry interpolation points in the boundary element
    real(kind=8) :: seb(2) !Current sub-element boundaries
    real(kind=8) :: eta !Natural coordinate in the original element (current integration point coordinate)
    type(point_2d) :: ipt !The current integration point
    real(kind=8) :: jac,v_normal(2),r_ri(3) !r_ri = {r, rx, ry}
    real(kind=8) :: u(2,2),t(2,2) !The current fundamental solutions values
    integer :: nsd,nip
    integer :: gn_num,i,j,k,kk
    !
    allocate (ni(be%nodes_num))
    if (2*size(ni).ne.size(h(1,:)) .or. 2*size(ni).ne.size(g(1,:)) &
       & .or. size(h(:,1)).ne.2 .or. size(g(:,1)).ne.2) then
      call error_message_stop("Incorrect sub-matrix h and/or g sizes for regular integration.")
    end if
    !
    !----------------------------------------------
    ! Evaluate the number of element subdivisions:
    !----------------------------------------------
    !
    call near_sing_parameters(cp,be,nsd,nip)
    !
    !--------------------------------------------------------------
    ! Fill array gi_pts with the be geometry interpolation points:
    !--------------------------------------------------------------
    !
    gn_num = be%nodes_num
    if (gn_num.eq.1) gn_num = 2
    allocate (gi_pts(gn_num))
    gi_pts(1) = be%bp_1; gi_pts(gn_num) = be%bp_2
    if (be%nodes_num.gt.2) then
      do i = 2,(be%nodes_num-1)
        gi_pts(i) = bn_list(be%incid(i))%point
      end do
    end if
    !
    !------------------------------------------------------
    ! Fill arrays with integration points and its weights:
    !------------------------------------------------------
    !
    allocate (eta_b(nip))
    allocate (ome(nip))
    call gauss_legendre_points_1d(nip,eta_b,ome)
    !
    !------------------------------------------------------------
    ! Evaluate Jacobian and normal vector for straight elements:
    !------------------------------------------------------------
    !
    if (gn_num.eq.2) call jac_normal_1d(jac,v_normal,gi_pts,0._8,be%nodes_num-1)
    !
    !-------------------
    ! Integration loop:
    !-------------------
    !
    h = 0._8; g = 0._8
    !
    seb(2) = -1._8
    !
    ! Loop on subdivisions:
    !
    do i = 1,nsd
      seb(1) = seb(2)
      if (i.eq.nsd) then
        seb(2) = 1._8
      else
        seb(2) = seb(1) + (2._8/nsd)
      end if
      !
      ! Loop on integration points:
      !
      do j = 1,nip
        eta = 0.5_8*(seb(1)+seb(2))+(eta_b(j)/nsd)
        call intrp_func_1d(ni,eta,be%lamb_1,be%lamb_2,be%nodes_num-1)
        ! Evaluate Jacobian and normal vector for non-straight elements:
        if (gn_num.gt.2) call jac_normal_1d(jac,v_normal,gi_pts,eta,be%nodes_num-1)
        ipt = cartesian_point_from_intrp_func_1d(gi_pts,eta)
        r_ri = distance_2d(cp,ipt)
        u = u_ij(r_ri,mu,nu,psk)
        t = t_ij(r_ri,nu,v_normal,psk)
        !
        kk = 1
        !
        ! Loop on interpolation functions:
        !
        do k = 1,be%nodes_num
          h(:,kk:kk+1) = h(:,kk:kk+1) + (t*jac*ome(j)/nsd)*ni(k)
          g(:,kk:kk+1) = g(:,kk:kk+1) + (u*jac*ome(j)/nsd)*ni(k)
          kk = kk + 2
        end do
        !
      end do
      !
    end do
    !
    deallocate (ni)
    deallocate (gi_pts)
    deallocate (eta_b)
    deallocate (ome)
    !
    return
  end subroutine displ_bie_reg_integr
  !****************************************************
  !****************************************************
  subroutine reg_integr_int_point_to_be(psk,mu,nu,cp,be,h,g,hs,gs)
    !----------------------------------------------------
    ! Perform numerical integration by Gauss-Legendre
    ! quadrature from internal points to boundary elements
    ! and store results in matrices h and g (values for
    ! displacement integral equation) and in matrices hs
    ! and gs (values for stress integral equation).
    !
    ! cp = the source point
    ! be = the boundary element
    ! h = Integral of [N_1*[T] ... N_n*[T]], n = order of be
    ! g = Integral of [N_1*[U] ... N_n*[U]], n = order of be
    ! hs = Integral of [N_1*[Ts] ... N_n*[Ts]], n = order of be
    ! gs = Integral of [N_1*[Us] ... N_n*[Us]], n = order of be
    ! N_k = interpolation functions
    ! [T] = 2x2 matrix with t_ij Kelvin fundamental solution
    ! [U] = 2x2 matrix with u_ij Kelvin fundamental solution
    ! [Ts] = 3x2 matrix with t_ijk
    ! [Us] = 3x2 matrix with u_ijk
    !
    ! Date: 27 June 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: psk
    real(kind=8),intent(in) :: mu,nu
    type(point_2d),intent(in) :: cp
    type(boundary_elem),intent(in) :: be
    real(kind=8),intent(out) :: h(:,:),g(:,:)
    real(kind=8),intent(out) :: hs(:,:),gs(:,:)
    real(kind=8),allocatable :: eta_b(:),ome(:)
    real(kind=8),allocatable :: ni(:) !Interpolation functions values
    type(point_2d),allocatable :: gi_pts(:) !Geometry interpolation points in the boundary element
    real(kind=8) :: seb(2) !Current sub-element boundaries
    real(kind=8) :: eta !Natural coordinate in the original element (current integration point coordinate)
    type(point_2d) :: ipt !The current integration point
    real(kind=8) :: jac,v_normal(2),r_ri(3) !r_ri = {r, rx, ry}
    real(kind=8) :: u(2,2),t(2,2) !The current fundamental solutions values
    real(kind=8) :: us(3,2),ts(3,2) !The current values of tensors u_ijk and t_ijk
    integer :: nsd,nip
    integer :: gn_num,i,j,k,kk
    !
    allocate (ni(be%nodes_num))
    !
    if (2*size(ni).ne.size(h(1,:)) .or. 2*size(ni).ne.size(g(1,:)) &
       & .or. size(h(:,1)).ne.2 .or. size(g(:,1)).ne.2) then
      call error_message_stop("Incorrect sub-matrix h and/or g sizes &
                              &for internal point integration.")
    end if
    !
    if (2*size(ni).ne.size(hs(1,:)) .or. 2*size(ni).ne.size(gs(1,:)) &
       & .or. size(hs(:,1)).ne.3 .or. size(gs(:,1)).ne.3) then
      call error_message_stop("Incorrect sub-matrix hs and/or gs sizes &
                              &for internal point integration.")
    end if
    !
    !
    !----------------------------------------------
    ! Evaluate the number of element subdivisions:
    !----------------------------------------------
    !
    call near_sing_parameters(cp,be,nsd,nip)
    !
    !--------------------------------------------------------------
    ! Fill array gi_pts with the be geometry interpolation points:
    !--------------------------------------------------------------
    !
    gn_num = be%nodes_num
    if (gn_num.eq.1) gn_num = 2
    allocate (gi_pts(gn_num))
    gi_pts(1) = be%bp_1; gi_pts(gn_num) = be%bp_2
    if (be%nodes_num.gt.2) then
      do i = 2,(be%nodes_num-1)
        gi_pts(i) = bn_list(be%incid(i))%point
      end do
    end if
    !
    !------------------------------------------------------
    ! Fill arrays with integration points and its weights:
    !------------------------------------------------------
    !
    allocate (eta_b(nip))
    allocate (ome(nip))
    call gauss_legendre_points_1d(nip,eta_b,ome)
    !
    !------------------------------------------------------------
    ! Evaluate Jacobian and normal vector for straight elements:
    !------------------------------------------------------------
    !
    if (gn_num.eq.2) call jac_normal_1d(jac,v_normal,gi_pts,0._8,be%nodes_num-1)
    !
    !-------------------
    ! Integration loop:
    !-------------------
    !
    h = 0._8; g = 0._8
    hs = 0._8; gs = 0._8
    !
    seb(2) = -1._8
    !
    ! Loop on subdivisions:
    !
    do i = 1,nsd
      seb(1) = seb(2)
      if (i.eq.nsd) then
        seb(2) = 1._8
      else
        seb(2) = seb(1) + (2._8/nsd)
      end if
      !
      ! Loop on integration points:
      !
      do j = 1,nip
        eta = 0.5_8*(seb(1)+seb(2))+(eta_b(j)/nsd)
        call intrp_func_1d(ni,eta,be%lamb_1,be%lamb_2,be%nodes_num-1)
        ! Evaluate Jacobian and normal vector for non-straight elements:
        if (gn_num.gt.2) call jac_normal_1d(jac,v_normal,gi_pts,eta,be%nodes_num-1)
        ipt = cartesian_point_from_intrp_func_1d(gi_pts,eta)
        r_ri = distance_2d(cp,ipt)
        u = u_ij(r_ri,mu,nu,psk)
        t = t_ij(r_ri,nu,v_normal,psk)
        us = u_ijk(r_ri,nu,psk)
        ts = t_ijk(r_ri,mu,nu,v_normal,psk)
        !
        kk = 1
        !
        ! Loop on interpolation functions:
        !
        do k = 1,be%nodes_num
          h(:,kk:kk+1) = h(:,kk:kk+1) + (t*jac*ome(j)/nsd)*ni(k)
          g(:,kk:kk+1) = g(:,kk:kk+1) + (u*jac*ome(j)/nsd)*ni(k)
          hs(:,kk:kk+1) = hs(:,kk:kk+1) + (ts*jac*ome(j)/nsd)*ni(k)
          gs(:,kk:kk+1) = gs(:,kk:kk+1) + (us*jac*ome(j)/nsd)*ni(k)
          kk = kk + 2
        end do
        !
      end do
      !
    end do
    !
    deallocate (ni)
    deallocate (gi_pts)
    deallocate (eta_b)
    deallocate (ome)
    !
    return
  end subroutine reg_integr_int_point_to_be
  !****************************************************
  !****************************************************
  subroutine reg_integr_j_integr_point_to_be(psk,mu,nu,cp,be,hg,gg,hs,gs)
    !----------------------------------------------------
    ! Perform numerical integration by Gauss-Legendre
    ! quadrature from J-integral points to boundary elements
    ! and store results in matrices h_ and g_ (values for
    ! displacement gradient integral equation) and in matrices hs
    ! and gs (values for stress integral equation).
    !
    ! cp = the source point
    ! be = the boundary element
    ! hg = Integral of [N_1*[T] ... N_n*[T]], n = order of be
    ! gg = Integral of [N_1*[U] ... N_n*[U]], n = order of be
    ! hs = Integral of [N_1*[Ts] ... N_n*[Ts]], n = order of be
    ! gs = Integral of [N_1*[Us] ... N_n*[Us]], n = order of be
    ! N_k = interpolation functions
    ! [T_] = 2x2 matrix with t_ij,k Kelvin fundamental solution
    ! [U_] = 2x2 matrix with u_ij,k Kelvin fundamental solution
    ! [Ts] = 3x2 matrix with t_ijk
    ! [Us] = 3x2 matrix with u_ijk
    !
    ! Date: 06 September 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: psk
    real(kind=8),intent(in) :: mu,nu
    type(point_2d),intent(in) :: cp
    type(boundary_elem),intent(in) :: be
    real(kind=8),intent(out) :: hg(:,:),gg(:,:)
    real(kind=8),intent(out) :: hs(:,:),gs(:,:)
    real(kind=8),allocatable :: eta_b(:),ome(:)
    real(kind=8),allocatable :: ni(:) !Interpolation functions values
    type(point_2d),allocatable :: gi_pts(:) !Geometry interpolation points in the boundary element
    real(kind=8) :: seb(2) !Current sub-element boundaries
    real(kind=8) :: eta !Natural coordinate in the original element (current integration point coordinate)
    type(point_2d) :: ipt !The current integration point
    real(kind=8) :: jac,v_normal(2),r_ri(3) !r_ri = {r, rx, ry}
    real(kind=8) :: ug(4,2),tg(4,2) !The current fundamental solutions values
    real(kind=8) :: us(3,2),ts(3,2) !The current values of tensors u_ijk and t_ijk
    integer :: nsd,nip
    integer :: gn_num,i,j,k,kk
    !
    allocate (ni(be%nodes_num))
    !
    if (2*size(ni).ne.size(hg(1,:)) .or. 2*size(ni).ne.size(gg(1,:)) &
       & .or. size(hg(:,1)).ne.4 .or. size(gg(:,1)).ne.4) then
      call error_message_stop("Incorrect sub-matrix hg and/or gg sizes &
                              &for internal point integration.")
    end if
    !
    if (2*size(ni).ne.size(hs(1,:)) .or. 2*size(ni).ne.size(gs(1,:)) &
       & .or. size(hs(:,1)).ne.3 .or. size(gs(:,1)).ne.3) then
      call error_message_stop("Incorrect sub-matrix hs and/or gs sizes &
                              &for internal point integration.")
    end if
    !
    !
    !----------------------------------------------
    ! Evaluate the number of element subdivisions:
    !----------------------------------------------
    !
    call near_sing_parameters(cp,be,nsd,nip)
    !
    !--------------------------------------------------------------
    ! Fill array gi_pts with the be geometry interpolation points:
    !--------------------------------------------------------------
    !
    gn_num = be%nodes_num
    if (gn_num.eq.1) gn_num = 2
    allocate (gi_pts(gn_num))
    gi_pts(1) = be%bp_1; gi_pts(gn_num) = be%bp_2
    if (be%nodes_num.gt.2) then
      do i = 2,(be%nodes_num-1)
        gi_pts(i) = bn_list(be%incid(i))%point
      end do
    end if
    !
    !------------------------------------------------------
    ! Fill arrays with integration points and its weights:
    !------------------------------------------------------
    !
    allocate (eta_b(nip))
    allocate (ome(nip))
    call gauss_legendre_points_1d(nip,eta_b,ome)
    !
    !------------------------------------------------------------
    ! Evaluate Jacobian and normal vector for straight elements:
    !------------------------------------------------------------
    !
    if (gn_num.eq.2) call jac_normal_1d(jac,v_normal,gi_pts,0._8,be%nodes_num-1)
    !
    !-------------------
    ! Integration loop:
    !-------------------
    !
    hg = 0._8; gg = 0._8
    hs = 0._8; gs = 0._8
    !
    seb(2) = -1._8
    !
    ! Loop on subdivisions:
    !
    do i = 1,nsd
      seb(1) = seb(2)
      if (i.eq.nsd) then
        seb(2) = 1._8
      else
        seb(2) = seb(1) + (2._8/nsd)
      end if
      !
      ! Loop on integration points:
      !
      do j = 1,nip
        eta = 0.5_8*(seb(1)+seb(2))+(eta_b(j)/nsd)
        call intrp_func_1d(ni,eta,be%lamb_1,be%lamb_2,be%nodes_num-1)
        ! Evaluate Jacobian and normal vector for non-straight elements:
        if (gn_num.gt.2) call jac_normal_1d(jac,v_normal,gi_pts,eta,be%nodes_num-1)
        ipt = cartesian_point_from_intrp_func_1d(gi_pts,eta)
        r_ri = distance_2d(cp,ipt)
        ug = u_ij_k(r_ri,mu,nu,psk)
        tg = t_ij_k(r_ri,nu,v_normal,psk)
        us = u_ijk(r_ri,nu,psk)
        ts = t_ijk(r_ri,mu,nu,v_normal,psk)
        !
        kk = 1
        !
        ! Loop on interpolation functions:
        !
        do k = 1,be%nodes_num
          hg(:,kk:kk+1) = hg(:,kk:kk+1) + (tg*jac*ome(j)/nsd)*ni(k)
          gg(:,kk:kk+1) = gg(:,kk:kk+1) + (ug*jac*ome(j)/nsd)*ni(k)
          hs(:,kk:kk+1) = hs(:,kk:kk+1) + (ts*jac*ome(j)/nsd)*ni(k)
          gs(:,kk:kk+1) = gs(:,kk:kk+1) + (us*jac*ome(j)/nsd)*ni(k)
          kk = kk + 2
        end do
        !
      end do
      !
    end do
    !
    deallocate (ni)
    deallocate (gi_pts)
    deallocate (eta_b)
    deallocate (ome)
    !
    return
  end subroutine reg_integr_j_integr_point_to_be
  !****************************************************
  !****************************************************
    subroutine reg_integr_grad_displ_to_be_point(psk,mu,nu,cp,be,hg,gg)
    !----------------------------------------------------
    ! Perform numerical integration by Gauss-Legendre
    ! quadrature from J-integral points to boundary elements
    ! and store results in matrices h_ and g_ (values for
    ! displacement gradient integral equation) and in matrices hs
    ! and gs (values for stress integral equation).
    !
    ! cp = the source point
    ! be = the boundary element
    ! hg = Integral of [N_1*[T] ... N_n*[T]], n = order of be
    ! gg = Integral of [N_1*[U] ... N_n*[U]], n = order of be
    ! hs = Integral of [N_1*[Ts] ... N_n*[Ts]], n = order of be
    ! gs = Integral of [N_1*[Us] ... N_n*[Us]], n = order of be
    ! N_k = interpolation functions
    ! [T_] = 2x2 matrix with t_ij,k Kelvin fundamental solution
    ! [U_] = 2x2 matrix with u_ij,k Kelvin fundamental solution
    ! [Ts] = 3x2 matrix with t_ijk
    ! [Us] = 3x2 matrix with u_ijk
    !
    ! Date: 06 September 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: psk
    real(kind=8),intent(in) :: mu,nu
    type(point_2d),intent(in) :: cp
    type(boundary_elem),intent(in) :: be
    real(kind=8),intent(out) :: hg(:,:),gg(:,:)
    real(kind=8),allocatable :: eta_b(:),ome(:)
    real(kind=8),allocatable :: ni(:) !Interpolation functions values
    type(point_2d),allocatable :: gi_pts(:) !Geometry interpolation points in the boundary element
    real(kind=8) :: seb(2) !Current sub-element boundaries
    real(kind=8) :: eta !Natural coordinate in the original element (current integration point coordinate)
    type(point_2d) :: ipt !The current integration point
    real(kind=8) :: jac,v_normal(2),r_ri(3) !r_ri = {r, rx, ry}
    real(kind=8) :: ug(4,2),tg(4,2) !The current fundamental solutions values
    integer :: nsd,nip
    integer :: gn_num,i,j,k,kk
    !
    allocate (ni(be%nodes_num))
    !
    if (2*size(ni).ne.size(hg(1,:)) .or. 2*size(ni).ne.size(gg(1,:)) &
       & .or. size(hg(:,1)).ne.4 .or. size(gg(:,1)).ne.4) then
      call error_message_stop("Incorrect sub-matrix hg and/or gg sizes &
                              &for internal point integration.")
    end if
    !
    !
    !
    !----------------------------------------------
    ! Evaluate the number of element subdivisions:
    !----------------------------------------------
    !
    call near_sing_parameters(cp,be,nsd,nip)
    !
    !--------------------------------------------------------------
    ! Fill array gi_pts with the be geometry interpolation points:
    !--------------------------------------------------------------
    !
    gn_num = be%nodes_num
    if (gn_num.eq.1) gn_num = 2
    allocate (gi_pts(gn_num))
    gi_pts(1) = be%bp_1; gi_pts(gn_num) = be%bp_2
    if (be%nodes_num.gt.2) then
      do i = 2,(be%nodes_num-1)
        gi_pts(i) = bn_list(be%incid(i))%point
      end do
    end if
    !
    !------------------------------------------------------
    ! Fill arrays with integration points and its weights:
    !------------------------------------------------------
    !
    allocate (eta_b(nip))
    allocate (ome(nip))
    call gauss_legendre_points_1d(nip,eta_b,ome)
    !
    !------------------------------------------------------------
    ! Evaluate Jacobian and normal vector for straight elements:
    !------------------------------------------------------------
    !
    if (gn_num.eq.2) call jac_normal_1d(jac,v_normal,gi_pts,0._8,be%nodes_num-1)
    !
    !-------------------
    ! Integration loop:
    !-------------------
    !
    hg = 0._8; gg = 0._8
    !
    seb(2) = -1._8
    !
    ! Loop on subdivisions:
    !
    do i = 1,nsd
      seb(1) = seb(2)
      if (i.eq.nsd) then
        seb(2) = 1._8
      else
        seb(2) = seb(1) + (2._8/nsd)
      end if
      !
      ! Loop on integration points:
      !
      do j = 1,nip
        eta = 0.5_8*(seb(1)+seb(2))+(eta_b(j)/nsd)
        call intrp_func_1d(ni,eta,be%lamb_1,be%lamb_2,be%nodes_num-1)
        ! Evaluate Jacobian and normal vector for non-straight elements:
        if (gn_num.gt.2) call jac_normal_1d(jac,v_normal,gi_pts,eta,be%nodes_num-1)
        ipt = cartesian_point_from_intrp_func_1d(gi_pts,eta)
        r_ri = distance_2d(cp,ipt)
        ug = u_ij_k(r_ri,mu,nu,psk)
        tg = t_ij_k(r_ri,nu,v_normal,psk)
        !
        kk = 1
        !
        ! Loop on interpolation functions:
        !
        do k = 1,be%nodes_num
          hg(:,kk:kk+1) = hg(:,kk:kk+1) + (tg*jac*ome(j)/nsd)*ni(k)
          gg(:,kk:kk+1) = gg(:,kk:kk+1) + (ug*jac*ome(j)/nsd)*ni(k)
          kk = kk + 2
        end do
        !
      end do
      !
    end do
    !
    deallocate (ni)
    deallocate (gi_pts)
    deallocate (eta_b)
    deallocate (ome)
    !
    return
  end subroutine reg_integr_grad_displ_to_be_point
  !****************************************************
  !****************************************************
  subroutine traction_bie_reg_integr(psk,mu,nu,cp,be,normal_cp,h,g)
    !----------------------------------------------------
    ! Perform numerical integration for traction boundary
    ! integral regular terms by Gauss-Legendre quadrature
    ! and store results in matrices h and g.
    !
    ! cp = the source point
    ! be = the boundary element
    ! h = Integral of [N_1*[T] ... N_n*[T]], n = be order
    ! g = Integral of [N_1*[U] ... N_n*[U]], n = be order
    ! N_k = interpolation functions
    ! [T] = 2x2 matrix with n_j : t_ijk
    ! [U] = 2x2 matrix with n_j : u_ijk
    !
    ! Date: 20 July 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: psk
    real(kind=8),intent(in) :: mu,nu,normal_cp(2)
    type(point_2d),intent(in) :: cp
    type(boundary_elem),intent(in) :: be
    real(kind=8),intent(out) :: h(:,:),g(:,:)
    real(kind=8),allocatable :: eta_b(:),ome(:)
    real(kind=8),allocatable :: ni(:) !Interpolation functions values
    type(point_2d),allocatable :: gi_pts(:) !Geometry interpolation points in the boundary element
    real(kind=8) :: seb(2) !Current sub-element boundaries
    real(kind=8) :: eta !Natural coordinate in the original element (current integration point coordinate)
    type(point_2d) :: ipt !The current integration point
    real(kind=8) :: jac,v_normal(2),r_ri(3) !r_ri = {r, rx, ry}
    real(kind=8) :: u(2,2),t(2,2) !The current fundamental solutions values
    integer :: nsd,nip
    integer :: gn_num,i,j,k,kk
    !
    allocate (ni(be%nodes_num))
    if (2*size(ni).ne.size(h(1,:)) .or. 2*size(ni).ne.size(g(1,:)) &
       & .or. size(h(:,1)).ne.2 .or. size(g(:,1)).ne.2) then
      call error_message_stop("Incorrect sub-matrix h and/or g sizes for regular integration.")
    end if
    !
    !----------------------------------------------
    ! Evaluate the number of element subdivisions:
    !----------------------------------------------
    !
    call near_sing_parameters(cp,be,nsd,nip)
    !
    !--------------------------------------------------------------
    ! Fill array gi_pts with the be geometry interpolation points:
    !--------------------------------------------------------------
    !
    gn_num = be%nodes_num
    if (gn_num.eq.1) gn_num = 2
    allocate (gi_pts(gn_num))
    gi_pts(1) = be%bp_1; gi_pts(gn_num) = be%bp_2
    if (be%nodes_num.gt.2) then
      do i = 2,(be%nodes_num-1)
        gi_pts(i) = bn_list(be%incid(i))%point
      end do
    end if
    !
    !------------------------------------------------------
    ! Fill arrays with integration points and its weights:
    !------------------------------------------------------
    !
    allocate (eta_b(nip))
    allocate (ome(nip))
    call gauss_legendre_points_1d(nip,eta_b,ome)
    !
    !------------------------------------------------------------
    ! Evaluate Jacobian and normal vector for straight elements:
    !------------------------------------------------------------
    !
    if (gn_num.eq.2) call jac_normal_1d(jac,v_normal,gi_pts,0._8,be%nodes_num-1)
    !
    !-------------------
    ! Integration loop:
    !-------------------
    !
    h = 0._8; g = 0._8
    !
    seb(2) = -1._8
    !
    ! Loop on subdivisions:
    !
    do i = 1,nsd
      seb(1) = seb(2)
      if (i.eq.nsd) then
        seb(2) = 1._8
      else
        seb(2) = seb(1) + (2._8/nsd)
      end if
      !
      ! Loop on integration points:
      !
      do j = 1,nip
        eta = 0.5_8*(seb(1)+seb(2))+(eta_b(j)/nsd)
        call intrp_func_1d(ni,eta,be%lamb_1,be%lamb_2,be%nodes_num-1)
        ! Evaluate Jacobian and normal vector for non-straight elements:
        if (gn_num.gt.2) call jac_normal_1d(jac,v_normal,gi_pts,eta,be%nodes_num-1)
        ipt = cartesian_point_from_intrp_func_1d(gi_pts,eta)
        r_ri = distance_2d(cp,ipt)
        u = normal_u_ijk(r_ri,normal_cp,nu,psk)
        t = normal_t_ijk(r_ri,mu,nu,v_normal,normal_cp,psk)
        !
        kk = 1
        !
        ! Loop on interpolation functions:
        !
        do k = 1,be%nodes_num
          h(:,kk:kk+1) = h(:,kk:kk+1) + (t*jac*ome(j)/nsd)*ni(k)
          g(:,kk:kk+1) = g(:,kk:kk+1) + (u*jac*ome(j)/nsd)*ni(k)
          kk = kk + 2
        end do
        !
      end do
      !
    end do
    !
    deallocate (ni)
    deallocate (gi_pts)
    deallocate (eta_b)
    deallocate (ome)
    !
    return
  end subroutine traction_bie_reg_integr
  !****************************************************
  !****************************************************
  subroutine stress_reg_integr(psk,mu,nu,cp,be,h,g)
    !----------------------------------------------------
    ! Perform numerical integration for stress boundary
    ! integral regular terms by Gauss-Legendre quadrature
    ! and store results in matrices h and g.
    !
    ! cp = the source point
    ! be = the boundary element
    ! h = Integral of [N_1*[T] ... N_n*[T]], n = be order
    ! g = Integral of [N_1*[U] ... N_n*[U]], n = be order
    ! N_k = interpolation functions
    ! [T] = 2x2 matrix with n_j : t_ijk
    ! [U] = 2x2 matrix with n_j : u_ijk
    !
    ! Date: 20 July 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: psk
    real(kind=8),intent(in) :: mu,nu
    type(point_2d),intent(in) :: cp
    type(boundary_elem),intent(in) :: be
    real(kind=8),intent(out) :: h(:,:),g(:,:)
    real(kind=8),allocatable :: eta_b(:),ome(:)
    real(kind=8),allocatable :: ni(:) !Interpolation functions values
    type(point_2d),allocatable :: gi_pts(:) !Geometry interpolation points in the boundary element
    real(kind=8) :: seb(2) !Current sub-element boundaries
    real(kind=8) :: eta !Natural coordinate in the original element (current integration point coordinate)
    type(point_2d) :: ipt !The current integration point
    real(kind=8) :: jac,v_normal(2),r_ri(3) !r_ri = {r, rx, ry}
    real(kind=8) :: u(3,2),t(3,2) !The current fundamental solutions values
    integer :: nsd,nip
    integer :: gn_num,i,j,k,kk
    !
    allocate (ni(be%nodes_num))
    if (2*size(ni).ne.size(h(1,:)) .or. 2*size(ni).ne.size(g(1,:)) &
       & .or. size(h(:,1)).ne.3 .or. size(g(:,1)).ne.3) then
      call error_message_stop("Incorrect sub-matrix h and/or g sizes for regular integration of stress components.")
    end if
    !
    !----------------------------------------------
    ! Evaluate the number of element subdivisions:
    !----------------------------------------------
    !
    call near_sing_parameters(cp,be,nsd,nip)
    !
    !--------------------------------------------------------------
    ! Fill array gi_pts with the be geometry interpolation points:
    !--------------------------------------------------------------
    !
    gn_num = be%nodes_num
    if (gn_num.eq.1) gn_num = 2
    allocate (gi_pts(gn_num))
    gi_pts(1) = be%bp_1; gi_pts(gn_num) = be%bp_2
    if (be%nodes_num.gt.2) then
      do i = 2,(be%nodes_num-1)
        gi_pts(i) = bn_list(be%incid(i))%point
      end do
    end if
    !
    !------------------------------------------------------
    ! Fill arrays with integration points and its weights:
    !------------------------------------------------------
    !
    allocate (eta_b(nip))
    allocate (ome(nip))
    call gauss_legendre_points_1d(nip,eta_b,ome)
    !
    !------------------------------------------------------------
    ! Evaluate Jacobian and normal vector for straight elements:
    !------------------------------------------------------------
    !
    if (gn_num.eq.2) call jac_normal_1d(jac,v_normal,gi_pts,0._8,be%nodes_num-1)
    !
    !-------------------
    ! Integration loop:
    !-------------------
    !
    h = 0._8; g = 0._8
    !
    seb(2) = -1._8
    !
    ! Loop on subdivisions:
    !
    do i = 1,nsd
      seb(1) = seb(2)
      if (i.eq.nsd) then
        seb(2) = 1._8
      else
        seb(2) = seb(1) + (2._8/nsd)
      end if
      !
      ! Loop on integration points:
      !
      do j = 1,nip
        eta = 0.5_8*(seb(1)+seb(2))+(eta_b(j)/nsd)
        call intrp_func_1d(ni,eta,be%lamb_1,be%lamb_2,be%nodes_num-1)
        ! Evaluate Jacobian and normal vector for non-straight elements:
        if (gn_num.gt.2) call jac_normal_1d(jac,v_normal,gi_pts,eta,be%nodes_num-1)
        ipt = cartesian_point_from_intrp_func_1d(gi_pts,eta)
        r_ri = distance_2d(cp,ipt)
        u = u_ijk(r_ri,nu,psk)
        t = t_ijk(r_ri,mu,nu,v_normal,psk)
        !
        kk = 1
        !
        ! Loop on interpolation functions:
        !
        do k = 1,be%nodes_num
          h(:,kk:kk+1) = h(:,kk:kk+1) + (t*jac*ome(j)/nsd)*ni(k)
          g(:,kk:kk+1) = g(:,kk:kk+1) + (u*jac*ome(j)/nsd)*ni(k)
          kk = kk + 2
        end do
        !
      end do
      !
    end do
    !
    deallocate (ni)
    deallocate (gi_pts)
    deallocate (eta_b)
    deallocate (ome)
    !
    return
  end subroutine stress_reg_integr
  !****************************************************
  !****************************************************
  subroutine near_sing_parameters(cp,be,nsd,nip)
    !----------------------------------------------------
    ! Compute the number of sub-elements and the number
    ! of Gauss points for regular integration, given a
    ! source point and a boundary element.
    !
    ! nsd = number of element's subdivisions for integration
    ! nip = number of integration points for each subdivision
    ! cp = the source point
    ! be = the boundary element
    !
    ! Date: 27 June 2022
    ! Last update: ---
    !----------------------------------------------------
    type(point_2d),intent(in) :: cp
    type(boundary_elem),intent(in) :: be
    integer,intent(out) :: nsd,nip
    real(kind=8) :: d0,dd,le
    integer :: i
    !
    nip = 20
    !
    dd = distance_only_2d(cp,bn_list(be%incid(1))%point)
    if (size(be%incid).gt.1) then
      do i = 2,size(be%incid)
        d0 = distance_only_2d(cp,bn_list(be%incid(i))%point)
        if (d0.lt.dd) dd = d0
      end do
    end if
    d0 = distance_only_2d(cp,be%bp_1)
    if (d0.lt.dd) dd = d0
    d0 = distance_only_2d(cp,be%bp_2)
    if (d0.lt.dd) dd = d0
    !
    le = be_size(be)
    !
    if (dd/le.gt.2.0) then
      nsd = 1
    else if (dd/le.lt.0.01) then
      nsd = 200
    else
      nsd = int(2*le/dd)
    end if
    !
    return
  end subroutine near_sing_parameters
  !****************************************************
  !****************************************************
 end module fracture_2d_regular_integration
