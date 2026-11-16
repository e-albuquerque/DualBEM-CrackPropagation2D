module interpolation_functions
!----------------------------------------------------
! A module for interpolation functions for
! continuous and discontinuous boundary elements
! and internal cells.
!
! Date: May 2018
! Last update: July 2021
!----------------------------------------------------
use vector_algebra
use error_messages
use geometry
implicit none
!
contains
  subroutine intrp_func_1d(ni,eta,lamb_1,lamb_2,io)
    !----------------------------------------------------
    ! Computes interpolation functions for line elements.
    ! ni(:) = vector with functions values at eta
    ! eta = the intrinsic coordinate
    ! lamb_1 = distance of the first interpolation point
    !          from the element center
    ! lamb_2 = distance of the last interpolation point
    !          from the element center
    ! io = interpolation order
    !
    ! Date: 24 May 2018
    ! Last update: 29 August 2018
    !----------------------------------------------------
    real(kind=8),intent(out) :: ni(:)
    real(kind=8),intent(in) :: eta,lamb_1,lamb_2
    integer,intent(in) :: io
    real(kind=8) :: a !temporary variable
    !
    select case (io)
    case(0) !constant approximation function
      ni(1) = 1.0_8
    case(1) !linear approximation function
      a = 1.0_8/(lamb_1+lamb_2)
      ni(1) = a*(lamb_2-eta)
      ni(2) = a*(lamb_1+eta)
    case(2) !quadratic approximation function
      a = 1.0_8/(lamb_1+lamb_2)
      ni(1) = a*eta*(eta-lamb_2)/lamb_1
      ni(2) = 1.0_8-(eta*(eta+lamb_1-lamb_2)/(lamb_1*lamb_2))
      ni(3) = a*eta*(eta+lamb_1)/lamb_2
    case default !error message
      call error_message_stop("Undefined interpolation order in an 1D element")
    end select
    return
  end subroutine intrp_func_1d
  !****************************************************
  !****************************************************
  subroutine deriv_intrp_func_1d(dni,eta,lamb_1,lamb_2,io)
    !----------------------------------------------------
    ! Computes the derivatives of interpolation functions
    ! for line elements.
    ! dni(:) = vector with functions derivative
    !          values at eta
    ! eta = the intrinsic coordinate
    ! lamb_1 = distance of the first interpolation point
    !          from the element center
    ! lamb_2 = distance of the last interpolation point
    !          from the element center
    ! io = interpolation order
    !
    ! Date: 29 May 2018
    ! Last update: 29 August 2018
    !----------------------------------------------------
    real(kind=8),intent(out) :: dni(:)
    real(kind=8),intent(in) :: eta,lamb_1,lamb_2
    integer,intent(in) :: io
    real(kind=8) :: a !temporary variable
    !
    select case (io)
    case(1) !linear approximation function
      a = 1.0_8/(lamb_1+lamb_2)
      dni(1) = -a
      dni(2) = a
    case(2) !quadratic approximation function
      a = 1.0_8/(lamb_1+lamb_2)
      dni(1) = a*(2.0_8*eta-lamb_2)/lamb_1
      dni(2) = -(2.0_8*eta+lamb_1-lamb_2)/(lamb_1*lamb_2)
      dni(3) = a*(2.0_8*eta+lamb_1)/lamb_2
    case default !error message
      call error_message_stop("Undefined interpolation order to calculate derivatives of 1D element")
    end select
    return
  end subroutine deriv_intrp_func_1d
  !****************************************************
  !****************************************************
  subroutine second_deriv_intrp_func_1d(d2ni,lamb_1,lamb_2,io)
    !----------------------------------------------------
    ! Computes the second derivatives of interpolation
    ! functions for line elements.
    ! d2ni(:) = vector with functions second derivative
    !           values at eta
    ! eta = the intrinsic coordinate
    ! lamb_1 = distance of the first interpolation point
    !          from the element center
    ! lamb_2 = distance of the last interpolation point
    !          from the element center
    ! io = interpolation order
    !
    ! Date: 19 September 2018
    ! Last update: 01 July 2021
    !----------------------------------------------------
    real(kind=8),intent(out) :: d2ni(:)
    real(kind=8),intent(in) :: lamb_1,lamb_2
    integer,intent(in) :: io
    real(kind=8) :: a !temporary variable
    !
    select case (io)
    case(1) !linear approximation function
      d2ni(1) = 0._8
      d2ni(2) = 0._8
    case(2) !quadratic approximation function
      a = 1.0_8/(lamb_1+lamb_2)
      d2ni(1) = 2.0_8*a/lamb_1
      d2ni(2) = 2.0_8/(lamb_1*lamb_2)
      d2ni(3) = 2.0_8*a/lamb_2
    case default !error message
      call error_message_stop("Undefined interpolation order to calculate second derivatives of 1D element")
    end select
    return
  end subroutine second_deriv_intrp_func_1d
  !****************************************************
  !****************************************************
  subroutine jac_normal_1d(jac,v_normal,gi_pts,eta,io)
    !----------------------------------------------------
    ! Computes the jacobian and normal vector at eta
    ! for line elements.
    ! jac = the element jacobian at eta
    ! v_normal(2) = normal vector at eta
    ! gi_pts(:) = geometry interpolation points
    ! eta = the intrinsic coordinate
    ! io = interpolation order
    !
    ! Date: 31 May 2018
    ! Last update: 29 August 2018
    !----------------------------------------------------
    real(kind=8),intent(out) :: jac
    real(kind=8),intent(out) :: v_normal(2)
    real(kind=8),intent(in) :: eta
    type(point_2d),intent(in) :: gi_pts(:)
    integer,intent(in) :: io
    real(kind=8),allocatable :: dni(:)
    real(kind=8),allocatable :: coords(:,:)
    real(kind=8) :: v_tang(2)
    integer :: gio !geometric interpolation order
    integer :: i
    !
    gio = io
    if (io.eq.0) gio = 1 !for constant elements
    allocate (dni(gio+1))
    call deriv_intrp_func_1d(dni,eta,1.0_8,1.0_8,gio)
    allocate (coords(size(gi_pts),2))
    do i=1,size(gi_pts)
      coords(i,1) = gi_pts(i)%x
      coords(i,2) = gi_pts(i)%y
    end do
    do i=1,2
      v_tang(i) = dot_product(dni,coords(:,i))
    end do
    jac = vector_len(v_tang)
    v_normal(1) = v_tang(2)
    v_normal(2) = -v_tang(1)
    call vector_norm(v_normal)
    deallocate(dni)
    deallocate(coords)
    return
  end subroutine jac_normal_1d
  !****************************************************
  !****************************************************
  subroutine jac_only_1d(jac,gi_pts,eta,io)
    !----------------------------------------------------
    ! Computes the jacobian at eta for line elements.
    ! jac = the element jacobian at eta
    ! gi_pts(:) = geometry interpolation points
    ! eta = the intrinsic coordinate
    ! io = interpolation order
    !
    ! Date: 19 September 2018
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(out) :: jac
    real(kind=8),intent(in) :: eta
    type(point_2d),intent(in) :: gi_pts(:)
    integer,intent(in) :: io
    real(kind=8),allocatable :: dni(:)
    real(kind=8),allocatable :: coords(:,:)
    real(kind=8) :: v_tang(2)
    integer :: gio !geometric interpolation order
    integer :: i
    !
    gio = io
    if (io.eq.0) gio = 1 !for constant elements
    allocate (dni(gio+1))
    call deriv_intrp_func_1d(dni,eta,1.0_8,1.0_8,gio)
    allocate (coords(size(gi_pts),2))
    do i=1,size(gi_pts)
      coords(i,1) = gi_pts(i)%x
      coords(i,2) = gi_pts(i)%y
    end do
    do i=1,2
      v_tang(i) = dot_product(dni,coords(:,i))
    end do
    jac = vector_len(v_tang)
    deallocate(dni)
    deallocate(coords)
    return
  end subroutine jac_only_1d
  !****************************************************
  !****************************************************
  function normal_1d(gi_pts,eta,io)
    !----------------------------------------------------
    ! Computes the normal vector at eta for line elements.
    !
    ! gi_pts(:) = geometry interpolation points
    ! eta = the intrinsic coordinate
    ! io = interpolation order
    !
    ! Date: 11 February 2019
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8) :: normal_1d(2)
    real(kind=8),intent(in) :: eta
    type(point_2d),intent(in) :: gi_pts(:)
    integer,intent(in) :: io
    real(kind=8),allocatable :: dni(:)
    real(kind=8),allocatable :: coords(:,:)
    real(kind=8) :: v_tang(2)
    integer :: gio !geometric interpolation order
    integer :: i
    !
    gio = io
    if (io.eq.0) gio = 1 !for constant elements
    allocate (dni(gio+1))
    call deriv_intrp_func_1d(dni,eta,1.0_8,1.0_8,gio)
    allocate (coords(size(gi_pts),2))
    do i=1,size(gi_pts)
      coords(i,1) = gi_pts(i)%x
      coords(i,2) = gi_pts(i)%y
    end do
    do i=1,2
      v_tang(i) = dot_product(dni,coords(:,i))
    end do
    normal_1d(1) = v_tang(2)
    normal_1d(2) = -v_tang(1)
    call vector_norm(normal_1d)
    deallocate(dni)
    deallocate(coords)
    return
  end function normal_1d
  !****************************************************
  !****************************************************
  function cartesian_point_from_intrp_func_1d(pts,eta)
    !----------------------------------------------------
    ! Returns the cartesian point from a given natural
    ! coordinate, eta.
    ! pts = an array of cartesian 2D points (the geometry
    !       interpolation points).
    ! eta = the unidimensional natural coordinate
    !
    ! Date: 26 July 2018
    ! Last update: ---
    !----------------------------------------------------
    type(point_2d),intent(in) :: pts(:)
    real(kind=8),intent(in) :: eta
    real(kind=8),allocatable :: ni(:)
    real(kind=8) :: x,y
    integer :: i
    type(point_2d) :: cartesian_point_from_intrp_func_1d
    !
    if (eta.lt.-1.0 .or. eta.gt.1.0) &
    call error_message_stop("invalid natural coordinate value.")
    allocate(ni(size(pts)))
    call intrp_func_1d(ni,eta,1.0_8,1.0_8,size(pts)-1)
    x = 0.0_8; y = 0.0_8;
    do i = 1,size(ni)
      x = x + pts(i)%x * ni(i)
      y = y + pts(i)%y * ni(i)
    end do
    cartesian_point_from_intrp_func_1d%x = x
    cartesian_point_from_intrp_func_1d%y = y
    return
  end function cartesian_point_from_intrp_func_1d
  !****************************************************
  !****************************************************
  subroutine intrp_func_2d(ni,eta_1,eta_2,lamb_1,lamb_2,lamb_3,lamb_4,io)
    !----------------------------------------------------
    ! Computes interpolation functions for quadrilateral
    ! surface elements.
    ! ni(:) = vector with functions values at (eta_1,eta_2)
    ! eta_1 and eta_2 = the intrinsic coordinates
    ! lamb_1 = distance of the interpolation points at the
    !          left edge from eta_1 = 0 axis
    ! lamb_2 = distance of the interpolation points at the
    !          right edge from eta_1 = 0 axis
    ! lamb_3 = distance of the interpolation points at the
    !          bottom edge from eta_2 = 0 axis
    ! lamb_4 = distance of the interpolation points at the
    !          upper edge from eta_2 = 0 axis
    ! io = interpolation order
    !
    ! Date: 29 May 2018
    ! Last update: 29 August 2018
    !----------------------------------------------------
    real(kind=8),intent(out) :: ni(:)
    real(kind=8),intent(in) :: eta_1,eta_2,lamb_1,lamb_2,lamb_3,lamb_4
    integer,intent(in) :: io
    real(kind=8) :: a !temporary variable
    !
    select case (io)
    case(0) !constant approximation function
      ni(1) = 1.0_8
    case(1) !linear approximation function
      a = 1.0_8/((lamb_1+lamb_2)*(lamb_3+lamb_4))
      ni(1) = a*(lamb_2-eta_1)*(lamb_4-eta_2)
      ni(2) = a*(lamb_1+eta_1)*(lamb_4-eta_2)
      ni(3) = a*(lamb_1+eta_1)*(lamb_3+eta_2)
      ni(4) = a*(lamb_2-eta_1)*(lamb_3+eta_2)
    case(2) !quadratic approximation function
      a = 1.0_8/((lamb_1+lamb_2)*(lamb_3+lamb_4))
      ni(1) = a*(lamb_2-eta_1)*(lamb_4-eta_2)*(-1._8-(eta_1/lamb_1)-(eta_2/lamb_3))
      ni(3) = a*(lamb_1+eta_1)*(lamb_4-eta_2)*(-1._8+(eta_1/lamb_2)-(eta_2/lamb_3))
      ni(5) = a*(lamb_1+eta_1)*(lamb_3+eta_2)*(-1._8+(eta_1/lamb_2)+(eta_2/lamb_4))
      ni(7) = a*(lamb_2-eta_1)*(lamb_3+eta_2)*(-1._8-(eta_1/lamb_1)+(eta_2/lamb_4))
      a = 1.0_8/(lamb_1*lamb_2*(lamb_3+lamb_4))
      ni(2) = a*(lamb_1+eta_1)*(lamb_4-eta_2)*(lamb_2-eta_1)
      ni(6) = a*(lamb_1+eta_1)*(lamb_3+eta_2)*(lamb_2-eta_1)
      a = 1.0_8/(lamb_3*lamb_4*(lamb_1+lamb_2))
      ni(4) = a*(lamb_1+eta_1)*(lamb_4-eta_2)*(lamb_3+eta_2)
      ni(8) = a*(lamb_2-eta_1)*(lamb_4-eta_2)*(lamb_3+eta_2)
    case default !error message
      call error_message_stop("Undefined interpolation order in an 2D element")
    end select
    return
  end subroutine intrp_func_2d
  !****************************************************
  !****************************************************
  subroutine deriv_intrp_func_2d(dni,eta_1,eta_2,lamb_1,lamb_2,lamb_3,lamb_4,io)
    !----------------------------------------------------
    ! Computes the derivatives of interpolation functions
    ! for quadrilateral surface elements.
    ! dni(:,:) = matrix with functions derivatives
    !            values at (eta_1,eta_2)
    ! eta_1 and eta_2 = the intrinsic coordinates
    ! lamb_1 = distance of the interpolation points at the
    !          left edge from eta_1 = 0 axis
    ! lamb_2 = distance of the interpolation points at the
    !          right edge from eta_1 = 0 axis
    ! lamb_3 = distance of the interpolation points at the
    !          bottom edge from eta_2 = 0 axis
    ! lamb_4 = distance of the interpolation points at the
    !          upper edge from eta_2 = 0 axis
    ! io = interpolation order
    !
    ! Date: 29 May 2018
    ! Last update: 29 August 2018
    !----------------------------------------------------
    real(kind=8),intent(out) :: dni(:,:) !row j: d(ni)/d(eta_j)
    real(kind=8),intent(in) :: eta_1,eta_2,lamb_1,lamb_2,lamb_3,lamb_4
    integer,intent(in) :: io
    real(kind=8) :: a !temporary variable
    !
    select case (io)
    case(1) !linear approximation function
      a = 1.0_8/((lamb_1+lamb_2)*(lamb_3+lamb_4))
      dni(1,1) = a*(eta_2-lamb_4)
      dni(1,2) = a*(lamb_4-eta_2)
      dni(1,3) = a*(lamb_3+eta_2)
      dni(1,4) = -a*(lamb_3+eta_2)
      dni(2,1) = a*(eta_1-lamb_2)
      dni(2,2) = -a*(lamb_1+eta_1)
      dni(2,3) = a*(lamb_1+eta_1)
      dni(2,4) = a*(lamb_2-eta_1)
    case(2) !quadratic approximation function
      a = 1.0_8/((lamb_1+lamb_2)*(lamb_3+lamb_4))
      dni(1,1) = a*((lamb_4-eta_2)* &
        (lamb_1*lamb_3-lamb_2*lamb_3+lamb_1*eta_2+2._8*lamb_3*eta_1))/(lamb_1*lamb_3)
      dni(1,3) = a*((lamb_4-eta_2)* &
        (lamb_1*lamb_3-lamb_2*lamb_3-lamb_2*eta_2+2._8*lamb_3*eta_1))/(lamb_2*lamb_3)
      dni(1,5) = a*((lamb_3+eta_2)* &
        (lamb_1*lamb_4-lamb_2*lamb_4+lamb_2*eta_2+2._8*lamb_4*eta_1))/(lamb_2*lamb_4)
      dni(1,7) = a*((lamb_3+eta_2)* &
        (lamb_1*lamb_4-lamb_2*lamb_4-lamb_1*eta_2+2._8*lamb_4*eta_1))/(lamb_1*lamb_4)
      dni(2,1) = a*((lamb_2-eta_1)* &
        (lamb_1*lamb_3-lamb_1*lamb_4+lamb_3*eta_1+2._8*lamb_1*eta_2))/(lamb_1*lamb_3)
      dni(2,3) = a*((lamb_1+eta_1)* &
        (lamb_2*lamb_3-lamb_2*lamb_4-lamb_3*eta_1+2._8*lamb_2*eta_2))/(lamb_2*lamb_3)
      dni(2,5) = a*((lamb_1+eta_1)* &
        (lamb_2*lamb_3-lamb_2*lamb_4+lamb_4*eta_1+2._8*lamb_2*eta_2))/(lamb_2*lamb_4)
      dni(2,7) = a*((lamb_2-eta_1)* &
        (lamb_1*lamb_3-lamb_1*lamb_4-lamb_4*eta_1+2._8*lamb_1*eta_2))/(lamb_1*lamb_4)
      a = 1.0_8/(lamb_1*lamb_2*(lamb_3+lamb_4))
      dni(1,2) = a*(eta_2-lamb_4)*(lamb_1-lamb_2+2._8*eta_1)
      dni(1,6) = -a*(eta_2+lamb_3)*(lamb_1-lamb_2+2._8*eta_1)
      dni(2,2) = a*(eta_1-lamb_2)*(lamb_1+eta_1)
      dni(2,6) = a*(lamb_2-eta_1)*(lamb_1+eta_1)
      a = 1.0_8/(lamb_3*lamb_4*(lamb_1+lamb_2))
      dni(1,4) = a*(lamb_4-eta_2)*(lamb_3+eta_2)
      dni(1,8) = a*(eta_2-lamb_4)*(lamb_3+eta_2)
      dni(2,4) = -a*(lamb_1+eta_1)*(lamb_3-lamb_4+2._8*eta_2)
      dni(2,8) = a*(eta_1-lamb_2)*(lamb_3-lamb_4+2._8*eta_2)
    case default !error message
      call error_message_stop("Undefined interpolation order to calculate derivatives of 2D element")
    end select
    return
  end subroutine deriv_intrp_func_2d
  !****************************************************
  !****************************************************
  subroutine jac_normal_2d(jac,v_normal,gi_pts,eta_1,eta_2,io)
    !----------------------------------------------------
    ! Computes the jacobian and normal vector at
    ! (eta_1,eta_2) for surface elements.
    ! jac = the element jacobian at (eta_1,eta_2)
    ! v_normal(3) = normal vector at (eta_1,eta_2)
    ! gi_pts(:) = geometry interpolation points
    ! eta_1 and eta_2 = the intrinsic coordinates
    ! io = interpolation order
    !
    ! Date: 31 May 2018
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(out) :: jac
    real(kind=8),intent(out) :: v_normal(3)
    type(point_3d),intent(in) :: gi_pts(:)
    real(kind=8),intent(in) :: eta_1,eta_2
    integer,intent(in) :: io
    real(kind=8),allocatable :: dni(:,:)
    real(kind=8),allocatable :: coords(:,:)
    real(kind=8) :: v_tang1(3),v_tang2(3)
    integer :: gio !geometric interpolation order
    integer :: i
    !
    gio = io
    if (io.eq.0) gio = 1 !for constant elements
    allocate (dni(2,gio*4))
    call deriv_intrp_func_2d(dni,eta_1,eta_2,1.0_8,1.0_8,1.0_8,1.0_8,gio)
    allocate (coords(size(gi_pts),3))
    do i = 1,size(gi_pts)
      coords(i,1) = gi_pts(i)%x
      coords(i,2) = gi_pts(i)%y
      coords(i,3) = gi_pts(i)%z
    end do
    do i=1,3
      v_tang1(i) = dot_product(dni(1,:),coords(:,i))
      v_tang2(i) = dot_product(dni(2,:),coords(:,i))
    end do
    v_normal = vector_cross_prod(v_tang1,v_tang2)
    jac = vector_len(v_normal)
    call vector_norm(v_normal)
    deallocate(dni)
    return
  end subroutine jac_normal_2d
  !****************************************************
  !****************************************************
end module interpolation_functions
