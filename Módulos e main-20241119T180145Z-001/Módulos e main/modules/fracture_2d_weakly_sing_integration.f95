module fracture_2d_weakly_sing_integration
!----------------------------------------------------
! A module for weakly singular integrations for
! two-dimensional isotropic linear elasticity problems.
!
! Date: June 2022
! Last update: ---
!----------------------------------------------------
use interpolation_functions
use numerical_integration_points
use fracture_2d_model
!
implicit none
!
real(kind=8),private :: four_pi = 12.566370614359173_8
real(kind=8),private :: eight_pi = 25.132741228718346_8
real(kind=8),private :: sixteen_pi = 50.265482457436692_8
!
contains
  subroutine displ_bie_weakly_sing_integr(psk,mu,nu,lcp,be,g)
    !----------------------------------------------------
    ! Perform the integration of weakly singular kernels,
    ! containing displacements from Kelvin's fundamental
    ! solution, associated to boundary elements in
    ! two-dimensional isotropic elasticity problems.
    !
    ! Constant and linear elements are analytically
    ! integrated.
    !
    ! For quadratic elements, the Gauss-Laguerre
    ! quadrature is used.
    !
    ! psk = plane state key (plane stress = 1, plane strain = 2)
    ! mu = shear modulus
    ! nu = Poisson ratio
    ! lcp = the local reference in the boundary element of the
    !       source point (an integer)
    ! be = the boundary element
    ! g = Integral of [N_1*[U] ... N_n*[U]], n = order of be
    ! N_k = interpolation functions
    ! [U] = 2x2 matrix with u_ij Kelvin fundamental solution
    !
    ! Date: 27 June 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: psk,lcp
    real(kind=8),intent(in) :: mu,nu
    type(boundary_elem),intent(in) :: be
    real(kind=8),intent(out) :: g(:,:)
    real(kind=8) :: nu_b
    !
    if (size(g(1,:)).ne.2*be%nodes_num .or. size(g(:,1)).ne.2) then
      call error_message_stop("Incorrect matrix g size for weakly singular integration.")
    end if
    !
    if (psk.eq.1) then
      nu_b = nu / (1._8 + nu)
    else
      nu_b = nu
    end if
    !
    select case (be%nodes_num)
    case(1) !be is a constant element
      call displ_bie_analytical_cte_be(mu,nu_b,be,g)
    case(2) !be is a linear element
      call displ_bie_analytical_lin_be(mu,nu_b,lcp,be,g)
    case(3) !be is a quadratic element
      call displ_bie_logarithmic_quadrature_quad_be(mu,nu_b,lcp,be,g)
    case default
      call error_message_stop("Undefined boundary element order during weakly singular integration.")
    end select
    !
    return
  end subroutine displ_bie_weakly_sing_integr
  !****************************************************
  !****************************************************
  subroutine displ_bie_analytical_cte_be(mu,nu,be,g)
    !----------------------------------------------------
    ! Analytically integrate weakly singular kernels,
    ! containing u_ij Kelvin's fundamental solution, for
    ! constant boundary elements.
    !
    ! mu = E/(2*(1+nu)), the Shear Modulus
    ! nu = the Poisson ratio
    !      (already modified for plane stress case)
    ! be = the constant boundary element
    ! g = Integral of [U]
    ! [U] = 2x2 matrix with u_ij Kelvin fundamental solution
    !
    ! Date: 27 June 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: mu,nu
    type(boundary_elem),intent(in) :: be
    real(kind=8),intent(inout) :: g(:,:)
    real(kind=8) :: le(3) !Boundary element length and its projections
    integer :: i
    !
    le = distance_2d(be%bp_1,be%bp_2)
    !
    do i = 1,2
      g(i,i) = le(i+1)*le(i+1)/le(1) + (3._8-4._8*nu)*le(1)*(1._8-log(0.5_8*le(1)))
    end do
    g(1,2) = le(2)*le(3)/le(1)
    g(2,1) = g(1,2)
    g = g/(eight_pi*mu*(1._8-nu))
    !
    return
  end subroutine displ_bie_analytical_cte_be
  !****************************************************
  !****************************************************
  subroutine displ_bie_analytical_lin_be(mu,nu,lcp,be,g)
    !----------------------------------------------------
    ! Analytically integrate weakly singular kernels,
    ! containing u_ij Kelvin's fundamental solution, for
    ! linear continuous or discontinuous boundary elements.
    !
    ! mu = E/(2*(1+nu)), the Shear Modulus
    ! nu = the Poisson ratio
    !      (already modified for plane stress case)
    ! lcp = the local reference in the boundary element of the
    !       source point
    ! be = the constant boundary element
    ! g = Integral of [N_1*[U] N_2*[U]]
    ! N_k = interpolation functions
    ! [U] = 2x2 matrix with u_ij Kelvin fundamental solution
    !
    ! Date: 27 June 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: mu,nu
    integer,intent(in) :: lcp
    type(boundary_elem),intent(in) :: be
    real(kind=8),intent(inout) :: g(:,:)
    real(kind=8) :: le(3) !Boundary element length and its projections
    real(kind=8) :: eta0 !The parametric coordinate of the source point
    real(kind=8) :: lb1,lb2
    real(kind=8) :: ak(2),bk(2)
    integer :: i,j
    !
    if (lcp.ne.1 .and. lcp.ne.2) call error_message_stop("Wrong source point reference &
    &during weakly singular integration of a linear element.")
    !
    le = distance_2d(be%bp_1,be%bp_2)
    !
    lb1 = be%lamb_1
    lb2 = be%lamb_2
    !
    eta0 = 0.0_8
    !
    select case (lcp)
    case (1)
      eta0 = -lb1
    case (2)
      eta0 = lb2
    case default
      call error_message_stop("Invalid reference for source point for weakly singular integration.")
    end select
    !
    bk(1) = 2.0_8*lb2/(lb1+lb2)
    bk(2) = 2.0_8*lb1/(lb1+lb2)
    !
    ak = 0._8
    !
    if (eta0.lt.1.0_8) then
      ak(1) = ak(1) + (lb2-eta0)*(1._8-eta0)*(log(0.5_8*le(1)*(1._8-eta0))-1._8) &
                  & - (1._8-eta0)*(1._8-eta0)*(0.5_8*log(0.5_8*le(1)*(1._8-eta0))-0.25_8)
      ak(2) = ak(2) + (lb1+eta0)*(1._8-eta0)*(log(0.5_8*le(1)*(1._8-eta0))-1._8) &
                  & + (1._8-eta0)*(1._8-eta0)*(0.5_8*log(0.5_8*le(1)*(1._8-eta0))-0.25_8)
    end if
    !
    if (eta0.gt.-1.0_8) then
      ak(1) = ak(1) + (lb2-eta0)*(1._8+eta0)*(log(0.5_8*le(1)*(1._8+eta0))-1._8) &
                  & + (1._8+eta0)*(1._8+eta0)*(0.5_8*log(0.5_8*le(1)*(1._8+eta0))-0.25_8)
      ak(2) = ak(2) + (lb1+eta0)*(1._8+eta0)*(log(0.5_8*le(1)*(1._8+eta0))-1._8) &
                  & - (1._8+eta0)*(1._8+eta0)*(0.5_8*log(0.5_8*le(1)*(1._8+eta0))-0.25_8)
    end if
    !
    ak = ak/(lb1+lb2)
    !
    do i = 1,2
      do j = 1,2
        g(i,j) = bk(1)*le(i+1)*le(j+1)/(le(1)*le(1))
        g(i,j+2) = bk(2)*le(i+1)*le(j+1)/(le(1)*le(1))
        if (i.eq.j) then
          g(i,j) = g(i,j) + (4._8*nu-3._8)*ak(1)
          g(i,j+2) = g(i,j+2) + (4._8*nu-3._8)*ak(2)
        end if
      end do
    end do
    !
    g = g*le(1)/(sixteen_pi*mu*(1._8-nu))
    !
    return
  end subroutine displ_bie_analytical_lin_be
  !****************************************************
  !****************************************************
  subroutine displ_bie_logarithmic_quadrature_quad_be(mu,nu,lcp,be,g)
    !----------------------------------------------------
    ! Perform the integration of weakly singular kernels,
    ! containing displacements from Kelvin's fundamental
    ! solution, associated to quadratic boundary elements
    ! in two-dimensional isotropic elasticity problems
    ! using the logarithmic quadrature procedures.
    !
    ! mu = E/(2*(1+nu)), the Shear Modulus
    ! nu = the Poisson ratio
    !      (already modified for plane stress case)
    ! lcp = the local reference in the boundary element of the
    !       source point
    ! be = the constant boundary element
    ! g = Integral of [N_1*[U] N_2*[U] N_3*[U]]
    ! N_k = interpolation functions
    ! [U] = 2x2 matrix with u_ij Kelvin fundamental solution
    !
    ! Date: 27 June 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: mu,nu
    integer,intent(in) :: lcp
    type(boundary_elem),intent(in) :: be
    real(kind=8),intent(inout) :: g(:,:)
    real(kind=8),allocatable :: eta(:),ome(:) !Integration points and its weights
    real(kind=8) :: ni(3) !Field interpolation functions
    real(kind=8) :: jac,eta0,eta_b !eta0 is the natural coordinate of the source point
    type(point_2d) :: gi_pts(3) !Geometry interpolation points of the original element
    type(point_2d) :: se_pts(3) !Sub-element interpolation points
    type(point_2d) :: cp_pt !The source point
    type(point_2d) :: it_pt !The integration point
    real(kind=8) :: intgr(3),r(3)
    real(kind=8) :: cte,f1,f2
    integer :: i,j
    !
    integer :: nip = 20 !Number of integration points for regular part
    integer :: nipl = 12 !Number of integration points for logarithmic quadrature
    !
    eta0 = 0._8
    if (lcp.ne.1 .and. lcp.ne.2 .and. lcp.ne.3) call error_message_stop("Wrong source &
    &point reference during weakly singular integration of a quadratic element.")
    !
    g = 0._8
    !
    !-----------------------------------------
    ! Natural coordinate of the source point:
    !-----------------------------------------
    !
    select case (lcp)
    case(1)
      eta0 = -be%lamb_1
    case(2)
      eta0 = 0._8
    case(3)
      eta0 = be%lamb_2
    case default
      call error_message_int_stop("Undefined local incidence number during weakly &
        &singular integration of a quadratic element. CP = ",lcp)
    end select
    !
    !--------------------------------------------------------------
    ! Fill array gi_pts with the be geometry interpolation points:
    !--------------------------------------------------------------
    !
    gi_pts(1) = be%bp_1; gi_pts(3) = be%bp_2
    gi_pts(2) = bn_list(be%incid(2))%point
    !
    !--------------------------------------------
    ! Cartesian coordinates of the source point:
    !--------------------------------------------
    !
    cp_pt = cartesian_point_from_intrp_func_1d(gi_pts,eta0)
    !
    ! -------------
    ! Sub-element 1
    ! -------------
    !
    if (eta0.ne.-1._8) then
      eta_b = 0.5*(eta0-1.0)
      se_pts(2) = cartesian_point_from_intrp_func_1d(gi_pts,eta_b)
      se_pts(1) = be%bp_1; se_pts(3) = cp_pt
      !
      ! -----------------------------
      ! Integrating the regular term:
      ! -----------------------------
      !
      allocate (eta(nip))
      allocate (ome(nip))
      call gauss_legendre_points_1d(nip,eta,ome)
      do i = 1,nip
        eta_b = 0.5_8*(eta0+1._8)*eta(i)+0.5_8*(eta0-1._8)
        it_pt = cartesian_point_from_intrp_func_1d(gi_pts,eta_b)
        r = distance_2d(cp_pt,it_pt)
        call intrp_func_1d(ni,eta_b,be%lamb_1,be%lamb_2,2)
        call jac_only_1d(jac,gi_pts,eta_b,2)
        call evaluate_f1_f2(3,se_pts,eta(i),f1,f2)
        !
        ni = ni*0.5_8*(1._8+eta0)*jac*ome(i)
        !
        do j = 1,2
          cte = 0.5_8*(r(j+1)*r(j+1)/(r(1)*r(1)) - 0.5_8*(3._8-4._8*nu)*log(f1*f1+f2*f2))
          g(j,j) = g(j,j) + cte*ni(1)
          g(j,j+2) = g(j,j+2) + cte*ni(2)
          g(j,j+4) = g(j,j+4) + cte*ni(3)
        end do
        !
        cte = 0.5_8*r(2)*r(3)/(r(1)*r(1))
        g(1,2) = g(1,2) + cte*ni(1)
        g(1,4) = g(1,4) + cte*ni(2)
        g(1,6) = g(1,6) + cte*ni(3)
      end do
      deallocate (eta)
      deallocate (ome)
      !
      ! ------------------------------
      ! Integrating the singular term:
      ! ------------------------------
      !
      allocate (eta(nipl))
      allocate (ome(nipl))
      call gauss_laguerre_points_1d(nipl,eta,ome)
      intgr = 0._8
      do i = 1,nipl
        eta_b = 0.5_8*(eta0+1._8)*(1._8-2._8*eta(i))+0.5_8*(eta0-1._8)
        call intrp_func_1d(ni,eta_b,be%lamb_1,be%lamb_2,2)
        call jac_only_1d(jac,gi_pts,eta_b,2)
        intgr = intgr + ni*jac*ome(i)
      end do
      intgr = intgr*0.5_8*(1._8+eta0)
      deallocate (eta)
      deallocate (ome)
      !
      cte = (3._8-4._8*nu)
      !
      do i = 1,2
        g(i,i) = g(i,i) + cte*intgr(1)
        g(i,i+2) = g(i,i+2) + cte*intgr(2)
        g(i,i+4) = g(i,i+4) + cte*intgr(3)
      end do
    end if
    !
    ! -------------
    ! Sub-element 2
    ! -------------
    !
    if (eta0.ne.1._8) then
      eta_b = 0.5*(eta0+1.0)
      se_pts(2) = cartesian_point_from_intrp_func_1d(gi_pts,eta_b)
      se_pts(1) = cp_pt; se_pts(3) = be%bp_2
      !
      ! -----------------------------
      ! Integrating the regular term:
      ! -----------------------------
      !
      allocate (eta(nip))
      allocate (ome(nip))
      call gauss_legendre_points_1d(nip,eta,ome)
      do i = 1,nip
        eta_b = 0.5_8*(1._8-eta0)*eta(i)+0.5_8*(1._8+eta0)
        it_pt = cartesian_point_from_intrp_func_1d(gi_pts,eta_b)
        r = distance_2d(cp_pt,it_pt)
        call intrp_func_1d(ni,eta_b,be%lamb_1,be%lamb_2,2)
        call jac_only_1d(jac,gi_pts,eta_b,2)
        call evaluate_f1_f2(1,se_pts,eta(i),f1,f2)
        !
        ni = ni*0.5_8*(1._8-eta0)*jac*ome(i)
        !
        do j = 1,2
          cte = 0.5_8*(r(j+1)*r(j+1)/(r(1)*r(1)) - 0.5_8*(3._8-4._8*nu)*log(f1*f1+f2*f2))
          g(j,j) = g(j,j) + cte*ni(1)
          g(j,j+2) = g(j,j+2) + cte*ni(2)
          g(j,j+4) = g(j,j+4) + cte*ni(3)
        end do
        !
        cte = 0.5_8*r(2)*r(3)/(r(1)*r(1))
        g(1,2) = g(1,2) + cte*ni(1)
        g(1,4) = g(1,4) + cte*ni(2)
        g(1,6) = g(1,6) + cte*ni(3)
      end do
      deallocate (eta)
      deallocate (ome)
      !
      ! ------------------------------
      ! Integrating the singular term:
      ! ------------------------------
      !
      allocate (eta(nipl))
      allocate (ome(nipl))
      call gauss_laguerre_points_1d(nipl,eta,ome)
      intgr = 0._8
      do i = 1,nipl
        eta_b = 0.5_8*(1._8-eta0)*(2._8*eta(i)-1._8)+0.5_8*(1._8+eta0)
        call intrp_func_1d(ni,eta_b,be%lamb_1,be%lamb_2,2)
        call jac_only_1d(jac,gi_pts,eta_b,2)
        intgr = intgr + ni*jac*ome(i)
      end do
      intgr = intgr*0.5_8*(1._8-eta0)
      deallocate (eta)
      deallocate (ome)
      !
      cte = (3._8-4._8*nu)
      !
      do i = 1,2
        g(i,i) = g(i,i) + cte*intgr(1)
        g(i,i+2) = g(i,i+2) + cte*intgr(2)
        g(i,i+4) = g(i,i+4) + cte*intgr(3)
      end do
    end if
    !
    ! ------------------------------------------------
    ! Complement matrix with already calculated terms:
    ! ------------------------------------------------
    !
    g(2,1) = g(1,2)
    g(2,3) = g(1,4)
    g(2,5) = g(1,6)
    !
    ! -------------------------------
    ! A constant common to all terms:
    ! -------------------------------
    !
    g = g/(four_pi*mu*(1._8-nu))
    !
    return
  end subroutine displ_bie_logarithmic_quadrature_quad_be
  !****************************************************
  !****************************************************
  subroutine evaluate_f1_f2(lcp,pts,eta,f1,f2)
    !----------------------------------------------------
    ! Evaluate parameters f1 and f2 for logarithmic
    ! quadrature procedure of quadratic elements.
    !
    ! cp = the local reference in the boundary element of the
    !      source point
    ! pts = an array with the element geometry interpolation
    !       points
    ! eta = the natural coordinate in the element
    !
    ! Date: 27 June 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: lcp
    real(kind=8),intent(in) :: eta
    type(point_2d),intent(in) :: pts(3)
    real(kind=8),intent(out) :: f1,f2
    !
    select case (lcp)
    case(1)
      f1 = -(2._8-eta)*pts(1)%x + eta*pts(3)%x + 2._8*(1._8-eta)*pts(2)%x
      f2 = -(2._8-eta)*pts(1)%y + eta*pts(3)%y + 2._8*(1._8-eta)*pts(2)%y
    case(2)
      f1 = 0.5_8*((eta-1._8)*pts(1)%x+(eta+1._8)*pts(3)%x) - eta*pts(2)%x
      f2 = 0.5_8*((eta-1._8)*pts(1)%y+(eta+1._8)*pts(3)%y) - eta*pts(2)%y
    case(3)
      f1 = -(2._8+eta)*pts(3)%x - eta*pts(1)%x + 2._8*(1._8+eta)*pts(2)%x
      f2 = -(2._8+eta)*pts(3)%y - eta*pts(1)%y + 2._8*(1._8+eta)*pts(2)%y
    case default
      call error_message_int_stop("Undefined local incidence number during weakly &
        &singular integration of a quadratic element. CP = ",lcp)
    end select
    !
    return
  end subroutine evaluate_f1_f2
  !****************************************************
  !****************************************************
end module fracture_2d_weakly_sing_integration
