module fracture_2d_model
!----------------------------------------------------
! A module that defines a dual boundary element
! method model for two-dimensional linear elastic
! fracture analysis of isotropic media.
!
! Date: August 2022
! Last update: ---
!----------------------------------------------------
use geometry
use vector_algebra
use interpolation_functions
use error_messages
implicit none
!
character(len=200),private :: model_name !The model name (many strings in general)
character(len=200),private :: model_id !The model Id (one string)
!
integer,private :: plane_state_key !Plane state key: =1 for plane stress, =2 for plane strain
!
real(kind=8),private :: young_modulus !The material modulus of elasticity (E)
real(kind=8),private :: poisson_ratio !The material Poisson ratio (nu)
!
type :: boundary_node
  integer :: bie !=1 displ. BIE on boundary, =2 displ. BIE on cracks, =3 tract.BIE on cracks
  type(point_2d) :: point !Cartesian coordinates of the node
  integer :: bc_dir !=0, Displacements and tractions in Cartesian coordinates; =1, Displacements and tractions normal/tangential to boundary
  integer :: bc(2) !=0, Displacement component prescribed; =1, Traction component prescribed; =2, Robin boundary condition (elastic supports)
  real(kind=8) :: h1(2),h2(2)!Robin b.c. parameters: ti = h1(i)*ui + h2(i)
  real(kind=8) :: disp(2) !Displacement components
  real(kind=8) :: trac(2) !Traction components
  real(kind=8) :: stress(3) !Stress components: {s11 s22 s12}
  real(kind=8) :: strain(3) !Strain components: {e11 e22 e12}
end type boundary_node
!
type :: boundary_elem
  integer :: nodes_num !Number of nodes in the element
  integer,allocatable :: incid(:) !Incidences numbers
  real(kind=8) :: lamb_1,lamb_2 !Discontinuous element parameters
  type(point_2d) :: bp_1,bp_2 !Element geometry boundary points
end type boundary_elem
!
type :: internal_node
  type(point_2d) :: point !Cartesian coordinates of the node.
  real(kind=8) :: disp(2) !Potential value on node.
  real(kind=8) :: stress(3) !Stress components: {s11 s22 s12}
  real(kind=8) :: strain(3) !Strain components: {e11 e22 e12}
end type internal_node
!
type :: crack
  integer :: tips_num !Number of tips (1 or 2), if 1, the tip is at the last element of the list
  integer :: ext_el_num(4) !Position in be_list of the initial and final elements of the crack
  real(kind=8) :: ki(2) !K_I stress intensity factor for each tip, ki(1) refer to the last element
  real(kind=8) :: kii(2) !K_II stress intensity factor for each tip, kii(1) refer to the last element
  real(kind=8) :: theta(2) !Angle of crack propagation direction for each tip, theta(1) refer to the last element
  real(kind=8) :: keq(2) !Equivalent stress intensity factor
  real(kind=8) :: length
  integer,allocatable :: nodes_num(:,:) !Matrix (:,2) containing the coincidental nodes position on bn_list
end type crack
!
integer :: bn_num !Number of boundary points in the model
integer :: be_num !Number of boundary elements in the model
integer :: in_num !Number of internal points in the model
integer :: cr_num !Number of cracks
!
type(boundary_node),allocatable :: bn_list(:)
type(boundary_elem),allocatable :: be_list(:)
type(internal_node),allocatable :: in_list(:)
type(crack),allocatable :: cr_list(:)
!
contains
  subroutine be_contains_node(be,cp,sing,lcn)
    !----------------------------------------------------
    ! Verify if a a boundary element contains a
    ! source point and returns the local reference
    ! if true.
    !
    ! be = the boundary element
    ! cp = the source point
    ! sing = a boolean which is .true. if be contains cp
    ! lcn = an integer with the local reference
    !
    ! Date: 24 June 2022
    ! Last update: ---
    !----------------------------------------------------
    type(boundary_elem),intent(in) :: be
    type(boundary_node),intent(in) :: cp
    logical,intent(out) :: sing
    integer,intent(out) :: lcn
    type(point_2d) :: pt
    real(kind=8) :: tol
    integer :: i
    !
    sing = .false.; lcn = 0
    tol = distance_only_2d(be%bp_1,be%bp_2)*1.e-8_8
    !
    ! Verifying if the node is an element incidence:
    !
    do i = 1,size(be%incid)
      pt = bn_list(be%incid(i))%point
      if (distance_only_2d(cp%point,pt).lt.tol) then
        sing = .true.
        lcn = i
        exit
      end if
    end do
    !
    ! Verifying if the node is not an element incidence
    ! but coincides with one of its extreme points:
    !
    if (.not.sing) then
      !
      if (distance_only_2d(cp%point,be%bp_1).lt.tol) then
        sing = .true.
        lcn = 1
      end if
      !
      if (distance_only_2d(cp%point,be%bp_2).lt.tol) then
        sing = .true.
        lcn = size(be%incid)
      end if
      !
    end if
    !
    return
  end subroutine be_contains_node
  !****************************************************
  !****************************************************
  function be_size(be)
    !----------------------------------------------------
    ! Evaluate a typical size for the boundary element.
    !
    ! be = the boundary element
    !
    ! Date: 24 June 2022
    ! Last update: ---
    !----------------------------------------------------
    type(boundary_elem),intent(in) :: be
    real(kind=8) :: be_size
    real(kind=8) :: le
    integer :: i,nn
    !
    nn = size(be%incid)
    !
    if (nn.lt.3) then
      le = distance_only_2d(be%bp_1,be%bp_2)
    else
      le = distance_only_2d(be%bp_1,bn_list(be%incid(2))%point)
      le = le + distance_only_2d(be%bp_2,bn_list(be%incid(nn-1))%point)
      if (size(be%incid).gt.3) then
        do i = 2,nn-2
          le = le + distance_only_2d(bn_list(i)%point,bn_list(i+1)%point)
        end do
      end if
    end if
    !
    be_size = le
    !
    return
  end function be_size
  !****************************************************
  !****************************************************
  function get_normal_to_a_boundary_point(inc_num)
    !----------------------------------------------------
    ! Get the outward normal to a boundary source
    ! point, given its position in the boundary nodes
    ! array (total_bn_list).
    !
    ! Date: 24 June 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8) :: get_normal_to_a_boundary_point(2)
    integer,intent(in) :: inc_num
    integer :: bes(2),linc(2),nbe
    type(point_2d),allocatable :: gi_pts(:)
    real(kind=8) :: eta
    real(kind=8) :: ni(2)
    integer :: i,j
    !
    if (inc_num.gt.bn_num) then
      call error_2_messages_stop("Boundary node number given for outward &
        &normal evaluation","doesn't exist in the model.")
    end if
    !
    get_normal_to_a_boundary_point = 0.0
    !
    nbe = 0
    !
    !--------------------------------------------------------------
    ! Storing references for boundary elements that contain the
    ! given incidence with its local position inside the elements:
    !--------------------------------------------------------------
    !
    do i = 1,be_num
      do j = 1,size(be_list(i)%incid)
        if (inc_num.eq.be_list(i)%incid(j)) then
          nbe = nbe + 1
          bes(nbe) = i
          linc(nbe) = j
        end if
      end do
      if (nbe.eq.2) exit !(Only) This line limits the number of elements that share a node by two
    end do
    !
    if (nbe.eq.0) then
      call error_2_messages_stop("Boundary node number given for outward &
                          &normal evaluation","doesn't exist in the model.")
    end if
    !
    !-------------------------------------------------------------
    ! Loop on boundary elements that contain the given incidence:
    !-------------------------------------------------------------
    !
    do i = 1,nbe
      !
      if (size(be_list(bes(i))%incid).eq.1) then
        !Constant element
        allocate(gi_pts(2))
        gi_pts(1) = be_list(bes(i))%bp_1
        gi_pts(2) = be_list(bes(i))%bp_2
      else
        !non-constant element
        allocate(gi_pts(size(be_list(bes(i))%incid)))
        do j = 1,size(be_list(bes(i))%incid)
          gi_pts(j) = bn_list(be_list(bes(i))%incid(j))%point
        end do
      end if
      !
      ! Calculating the natural coordinate of the
      ! given incidence for the current boundary element:
      !
      select case (size(be_list(bes(i))%incid))
      case(1) !Constant element
        eta = 0.0
      case(2) !Linear element
        if (linc(i).eq.1) eta = -be_list(i)%lamb_1
        if (linc(i).eq.2) eta = be_list(i)%lamb_2
      case(3) !Quadratic element
        if (linc(i).eq.1) eta = -be_list(i)%lamb_1
        if (linc(i).eq.2) eta = 0.0
        if (linc(i).eq.3) eta = be_list(i)%lamb_2
      case default
        call error_message_stop("Not able to compute the normal &
                                &vector for a boundary element.")
      end select
      !
      ni = normal_1d(gi_pts,eta,size(be_list(i)%incid)-1)
      deallocate(gi_pts)
      get_normal_to_a_boundary_point = get_normal_to_a_boundary_point + ni
      !
    end do
    !
    !-----------------------------------------------------------------------------
    ! Averaging result if the given point is associated to more than one element:
    !-----------------------------------------------------------------------------
    !
    if (nbe.gt.1) then
      get_normal_to_a_boundary_point = get_normal_to_a_boundary_point/nbe
      call vector_norm(get_normal_to_a_boundary_point)
    end if
    !
    return
  end function get_normal_to_a_boundary_point
  !****************************************************
  !****************************************************
  subroutine get_bound_elem_geom_points(be,be_pts)
    !----------------------------------------------------
    ! Return the list of the geometric points of a
    ! given boundary element.
    !
    ! Date: 20 July 2022
    ! Last update: ---
    !----------------------------------------------------
    type(boundary_elem),intent(in) :: be
    type(point_2d),allocatable,intent(out) :: be_pts(:)
    integer :: ngp,i
    !
    ngp = be%nodes_num
    if (ngp.eq.1) ngp = 2
    allocate(be_pts(ngp))
    be_pts(1) = be%bp_1; be_pts(ngp) = be%bp_2
    if (be%nodes_num.gt.2) then
      do i = 2,(be%nodes_num-1)
        be_pts(i) = bn_list(be%incid(i))%point
      end do
    end if
    !
    return
  end subroutine get_bound_elem_geom_points
  !****************************************************
  !****************************************************
  subroutine get_strain_from_stress_elast(psk,nu,el,sig,eps)
    !----------------------------------------------------
    ! Get strain (eps) from given stress (sig).
    !
    ! Date: 27 June 2022
    ! Last update: 22 July 2022
    !----------------------------------------------------
    integer,intent(in) :: psk
    real(kind=8),intent(in) :: nu,el
    real(kind=8),intent(in) :: sig(3)
    real(kind=8),intent(out) :: eps(3)
    real(kind=8) :: dd(3,3)
    real(kind=8) :: nub,eb
    !
    if (psk.eq.1) then
      !Plane stress:
      nub = nu
      eb = el
    else
      !Plane strain:
      nub = nu/(1.0_8-nu)
      eb = el/(1.0_8-nu*nu)
    end if
    !
    dd = 0.0_8
    !
    dd(1,1) = 1.0_8/eb
    dd(2,2) = 1.0_8/eb
    dd(1,2) = -nub/eb
    dd(2,1) = -nub/eb
    dd(3,3) = (1.0_8+nu)/el
    !
    eps = matmul(dd,sig)
    !
    return
  end subroutine get_strain_from_stress_elast
  !****************************************************
  !****************************************************
  subroutine get_kappa(kappa)
    !----------------------------------------------------
    ! Get the variable kappa, defined by the plane state
    ! key.
    !
    ! Date: 10 September 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8), intent(out) :: kappa
    real(kind=8) :: nu
    integer :: psk
    !
    call get_plane_state_key(psk)
    call get_poisson_ratio(nu)
    !
    kappa = 0._8
    !
    if (psk .eq. 1) then
        kappa = 3 - 4*(nu/(1+nu))
        else if (psk .eq. 2) then
            kappa = 3 - 4*nu
    end if
    !
    return
  end subroutine get_kappa
  !****************************************************
  !****************************************************
  subroutine get_shear_modulus(g)
    !----------------------------------------------------
    ! Get the material shear modulus, G:
    !
    ! G = E/[2*(1+nu)]
    !
    ! E = Young modulus
    ! nu = Poisson ratio
    !
    ! Date: 24 June 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(out) :: g
    !
    g = 0.5_8*young_modulus/(1.0_8+poisson_ratio)
    !
    return
  end subroutine get_shear_modulus
  !****************************************************
  !****************************************************
  subroutine set_plane_state_key(psk)
    !----------------------------------------------------
    ! Set the plane state key:
    !
    ! 1 for plane stress, 2 for plane strain.
    !
    ! Date: 24 June 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: psk
    !
    if (psk.lt.1 .or. psk.gt.2) then
      call error_message_stop("Invalid key for set the plane state.")
    end if
    !
    plane_state_key = psk
    !
    return
  end subroutine set_plane_state_key
  !****************************************************
  !****************************************************
  subroutine get_plane_state_key(psk)
    !----------------------------------------------------
    ! Get the plane state key:
    !
    ! 1 for plane stress, 2 for plane strain.
    !
    ! Date: 24 June 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(out) :: psk
    !
    psk = plane_state_key
    !
    return
  end subroutine get_plane_state_key
  !****************************************************
  !****************************************************
  subroutine set_young_modulus(e)
    !----------------------------------------------------
    ! Set the material modulus of elasticity.
    !
    ! Date: 24 June 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: e
    young_modulus = e
    return
  end subroutine set_young_modulus
  !****************************************************
  !****************************************************
  subroutine get_young_modulus(e)
    !----------------------------------------------------
    ! Get the material modulus of elasticity.
    !
    ! Date: 24 June 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(out) :: e
    e = young_modulus
    return
  end subroutine get_young_modulus
  !****************************************************
  !****************************************************
  subroutine set_poisson_ratio(nu)
    !----------------------------------------------------
    ! Set the material Poisson ratio.
    !
    ! Date: 24 June 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: nu
    poisson_ratio = nu
    return
  end subroutine set_poisson_ratio
  !****************************************************
  !****************************************************
  subroutine get_poisson_ratio(nu)
    !----------------------------------------------------
    ! Get the material Poisson ratio.
    !
    ! Date: 24 June 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(out) :: nu
    nu = poisson_ratio
    return
  end subroutine get_poisson_ratio
  !****************************************************
  !****************************************************
  subroutine set_model_name(mod_name)
    !----------------------------------------------------
    ! Set a name for the model, storing it in
    ! variable "model_name".
    !
    ! Date: 24 June 2022
    ! Last update: ---
    !----------------------------------------------------
    character(len=*),intent(in) :: mod_name
    model_name = mod_name
    return
  end subroutine set_model_name
  !****************************************************
  !****************************************************
  subroutine get_model_name(mod_name)
    !----------------------------------------------------
    ! Get the name for the model, which is stored in
    ! variable "model_name".
    !
    ! Date: 24 June 2022
    ! Last update: ---
    !----------------------------------------------------
    character(len=*),intent(out) :: mod_name
    mod_name = model_name
    return
  end subroutine get_model_name
  !****************************************************
  !****************************************************
  subroutine set_model_id(mod_id)
    !----------------------------------------------------
    ! Set the model id (one string), storing it in
    ! variable "model_id".
    !
    ! Date: 24 June 2022
    ! Last update: ---
    !----------------------------------------------------
    character(len=*),intent(in) :: mod_id
    model_id = mod_id
    return
  end subroutine set_model_id
  !****************************************************
  !****************************************************
  subroutine get_model_id(mod_id)
    !----------------------------------------------------
    ! Get the model id, which is stored in variable
    ! "model_id".
    !
    ! Date: 24 June 2022
    ! Last update: ---
    !----------------------------------------------------
    character(len=*),intent(out) :: mod_id
    mod_id = model_id
    return
  end subroutine get_model_id
  !****************************************************
  !****************************************************
end module fracture_2d_model
