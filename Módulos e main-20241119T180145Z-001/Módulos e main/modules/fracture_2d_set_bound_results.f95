module fracture_2d_set_bound_results
!----------------------------------------------------
! A module to set results on boundary nodes after
! algebraic system evaluation.
!
! Date: August 2022
! Last update: ---
!----------------------------------------------------
use error_messages
use interpolation_functions
use fracture_2d_model
use fracture_2d_rad_exp_intgr
use fracture_2d_regular_integration
!
implicit none
!
contains
  subroutine set_results_on_boundary_nodes(u)
    !----------------------------------------------------
    ! Set calculated values (displacement or traction)
    ! on boundary nodes from a given solution vector {u}.
    !
    ! This subroutine assumes that a unique traction
    ! value is permitted for each node.
    !
    ! Date: 09 August 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: u(:)
    real(kind=8) :: nv(2)
    real(kind=8) :: aux
    integer :: i,j,k
    !
    if (size(u).ne.2*bn_num) then
      call error_message_stop("Results vector size is different from &
                              &boundary degrees of freedom.")
    end if
    !
    do i = 1,bn_num
      !
      if (bn_list(i)%bc_dir.ne.0) then
        !
        !-----------------------------------------------------------------
        ! Boundary conditions prescribed in normal/tangential directions:
        !-----------------------------------------------------------------
        !
        nv = get_normal_to_a_boundary_point(i)
        !
        aux = bn_list(i)%disp(1)*nv(1) - bn_list(i)%disp(2)*nv(2)
        bn_list(i)%disp(2) = bn_list(i)%disp(1)*nv(2) + bn_list(i)%disp(2)*nv(1)
        bn_list(i)%disp(1) = aux
        !
        aux = bn_list(i)%trac(1)*nv(1) - bn_list(i)%trac(2)*nv(2)
        bn_list(i)%trac(2) = bn_list(i)%trac(1)*nv(2) + bn_list(i)%trac(2)*nv(1)
        bn_list(i)%trac(1) = aux
        !
        do k = 1,2
          !
          if (k.eq.2) then
            aux = nv(1)
            nv(1) = -nv(2)
            nv(2) = aux
          end if
          !
          select case (bn_list(i)%bc(k))
          case(0)
            !Displacement component prescribed (Dirichlet b.c.):
            do j = 1,2
              bn_list(i)%trac(j) = bn_list(i)%trac(j) + u(2*i+k-2)*nv(j)
            end do
          case(1)
            !Traction component prescribed (Neumann b.c.):
            do j = 1,2
              bn_list(i)%disp(j) = bn_list(i)%disp(j) + u(2*i+k-2)*nv(j)
            end do
          case(2)
            !Elastic support (Robin b.c.):
            aux = u(2*i+k-2)*bn_list(i)%h1(k) + bn_list(i)%h2(k)
            do j = 1,2
              bn_list(i)%disp(j) = bn_list(i)%disp(j) + u(2*i+k-2)*nv(j)
              bn_list(i)%trac(j) = bn_list(i)%trac(j) + aux*nv(j)
            end do
          case default
            !Invalid boundary condition
            call error_message_int_stop("Undefined boundary condition at node: ",i)
          end select
          !
        end do
        !
      else
        !
        !---------------------------------------------------------
        ! Boundary conditions prescribed in Cartesian directions:
        !---------------------------------------------------------
        !
        do k = 1,2
          select case (bn_list(i)%bc(k))
          case(0)
            !Displacement component prescribed (Dirichlet b.c.):
            bn_list(i)%trac(k) = u(2*i+k-2)
          case(1)
            !Traction component prescribed (Neumann b.c.):
            bn_list(i)%disp(k) = u(2*i+k-2)
          case(2)
            !Elastic support (Robin b.c.):
            bn_list(i)%disp(k) = u(2*i+k-2)
            bn_list(i)%trac(k) = u(2*i+k-2)*bn_list(i)%h1(k) + bn_list(i)%h2(k)
          case default
            !Invalid boundary condition
            call error_message_int_stop("Undefined boundary condition at node: ",i)
          end select
        end do
        !
      end if
      !
    end do
    !
    return
  end subroutine set_results_on_boundary_nodes
  !****************************************************
  !****************************************************
  subroutine evaluate_stress_strain_on_boundary_nodes(bound_recovery)
    !----------------------------------------------------
    ! Verify if the stress and strains components will be solved by the boundary
    ! recovery method or by solving the strongly singular and hyper singular integrations
    ! and chose the appropriate code.
    !
    ! Date: 09 August 2022
    ! Last update: ---
    !---------------------------------------------------
    logical, intent(in) :: bound_recovery
    !
    if (bound_recovery) then
        call stress_strain_using_boundary_recovery()
        else
            call stress_strain_using_singular_integration()
    end if
    !
    return
  end subroutine evaluate_stress_strain_on_boundary_nodes
  !****************************************************
  !****************************************************
  subroutine stress_strain_using_boundary_recovery()
    !----------------------------------------------------
    ! Evaluate and set strain and stress components on
    ! boundary nodes using the boundary recovery method.
    !
    ! Date: 27 June 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8) :: grad_u(2,2) !Displacement gradient tensor: grad_u(i,j) = u_i,j
    real(kind=8) :: l1(2) !l1 = [n_1 n_2], n_i = normal vector components
    real(kind=8) :: l2(2) !l2 = [-n_2 n_1], n_i = normal vector components
    real(kind=8) :: aa,bb !Temporary parameters: A_kmji and B_kmj
    real(kind=8),allocatable :: u(:,:) !Displacements at the element's nodes
    real(kind=8) :: t(2) !Traction at the current node
    real(kind=8) :: mu,nu,nub,lb,c1,c2,c3 !Material constants
    real(kind=8) :: cc !A temporary constant
    real(kind=8) :: eta,jac !The current natural coordinate and Jacobian
    real(kind=8) :: d_eta !Increments in the natural coordinate
    real(kind=8),allocatable :: dni(:) !interpolation functions derivatives
    type(point_2d),allocatable :: gi_pts(:) !interpolation points for an element
    integer,allocatable :: nefn(:) !Stores the number of elements per node
    integer :: psk !Plane state key (=1 for plane stress)
    integer :: inum !The current incidence number
    integer :: nn !Number of nodes in an element
    integer :: i,j,k,l,m,n
    !
    do k = 1,bn_num
      bn_list(k)%strain = 0.0_8
      bn_list(k)%stress = 0.0_8
    end do
    !
    allocate(nefn(bn_num))
    nefn = 0
    !
    !---------------------------------
    ! Setting the material constants:
    !---------------------------------
    !
    call get_plane_state_key(psk)
    call get_shear_modulus(mu)
    call get_poisson_ratio(nu)
    !
    if (psk.eq.1) then
      nub = nu/(1.0_8+nu)
    else
      nub = nu
    end if
    !
    lb = 2.0_8*mu*nub/(1.0_8-2.0_8*nub)
    c1 = 1.0_8/mu
    c2 = lb/(2.0_8*mu+lb)
    c3 = 1.0_8/(2.0_8*mu+lb)
    !
    !----------------------------
    ! Loop on boundary elements:
    !----------------------------
    !
    do l = 1,be_num
      !
      nn = be_list(l)%nodes_num
      if (nn.eq.1) cycle !Constant element
      !
      allocate(gi_pts(nn))
      allocate(u(2,nn))
      !
      gi_pts(1) = be_list(l)%bp_1
      gi_pts(nn) = be_list(l)%bp_2
      !
      do k = 1,nn
        !
        if (k.ne.1 .and. k.ne.nn) then
          gi_pts(k) = bn_list(be_list(l)%incid(k))%point
        end if
        !
        u(:,k) = bn_list(be_list(l)%incid(k))%disp
        !
      end do
      !
      d_eta = 2.0_8/(nn-1)
      allocate(dni(nn))
      !
      !--------------------------------------
      ! Loop on boundary element incidences:
      !--------------------------------------
      !
      do n = 1,nn
        !
        inum = be_list(l)%incid(n)
        nefn(inum) = nefn(inum) + 1
        !
        if (n.eq.1) then
          eta = -be_list(l)%lamb_1
        else if (n.eq.nn) then
          eta = be_list(l)%lamb_2
        else
          eta = -1.0_8 + d_eta*(n-1)
        end if
        !
        call deriv_intrp_func_1d(dni,eta,be_list(l)%lamb_1,be_list(l)%lamb_2,nn-1)
        call jac_only_1d(jac,gi_pts,eta,nn-1)
        !
        grad_u = 0.0_8
        l1 = get_normal_to_a_boundary_point(inum)
        l2(1) = -l1(2); l2(2) = l1(1)
        !
        t = bn_list(inum)%trac
        !
        !---------------------------------------
        ! Evaluating the displacement gradient:
        !---------------------------------------
        !
        do k = 1,2
          do m = 1,2
            do j = 1,2
              !
              ! B_kmj:t_j
              !
              bb = c1*l2(k)*l1(m)*l2(j) + c3*l1(k)*l1(m)*l1(j)
              grad_u(k,m) = grad_u(k,m) + bb*t(j)
              !
              do i = 1,nn
                !
                ! A_kmji:u_ij
                !
                aa = (l1(k)*l2(m)*l1(j) - l2(k)*l1(m)*l1(j) + l2(k)*l2(m)*l2(j) &
                   & - c2*l1(k)*l1(m)*l2(j)) * (1.0_8/jac) * dni(i)
                grad_u(k,m) = grad_u(k,m) + aa*u(j,i)
                !
              end do
            end do
          end do
        end do
        !
        bn_list(inum)%strain(1) = bn_list(inum)%strain(1) + grad_u(1,1)
        bn_list(inum)%strain(2) = bn_list(inum)%strain(2) + grad_u(2,2)
        bn_list(inum)%strain(3) = bn_list(inum)%strain(3) + 0.5_8*(grad_u(1,2)+grad_u(2,1))
        !
        cc = lb*(grad_u(1,1)+grad_u(2,2))
        bn_list(inum)%stress(1) = bn_list(inum)%stress(1) + 2.0_8*mu*grad_u(1,1) + cc
        bn_list(inum)%stress(2) = bn_list(inum)%stress(2) + 2.0_8*mu*grad_u(2,2) + cc
        bn_list(inum)%stress(3) = bn_list(inum)%stress(3) + mu*(grad_u(1,2)+grad_u(2,1))
        !
      end do
      !
      deallocate(gi_pts)
      deallocate(u)
      deallocate(dni)
      !
    end do
    !
    !--------------------
    ! Averaging process:
    !--------------------
    !
    do i = 1,bn_num
      if (nefn(i).gt.1) then
        bn_list(i)%strain = bn_list(i)%strain/nefn(i)
        bn_list(i)%stress = bn_list(i)%stress/nefn(i)
      end if
    end do
    !
    deallocate(nefn)
    !
    return
  end subroutine stress_strain_using_boundary_recovery
   !****************************************************
  !****************************************************
  subroutine stress_strain_using_singular_integration()
    !----------------------------------------------------
    ! Calculate strains and stresses at boundary points
    ! for two-dimensional isotropic linear elasticity problems.
    !
    ! Stress components are evaluated by Integral Equations
    ! and strain components are evaluated by Hooke's Law.
    !
    ! Date: 09 August 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8) :: eps(3),sig(3)
    real(kind=8),allocatable :: hs(:,:),gs(:,:)
    !real(kind=8),allocatable :: u(:),t(:)
    real(kind=8) :: u(6),t(6)
    real(kind=8) :: mu,nu,el
    integer :: psk
    type(boundary_node) :: cp
    type(boundary_elem) :: be
    integer :: i,j,k,lcn
    logical :: sing
    !
    !----------------------------------------
    ! Analysis type and material properties:
    !----------------------------------------
    !
    call get_plane_state_key(psk)
    call get_young_modulus(el)
    call get_shear_modulus(mu)
    call get_poisson_ratio(nu)
    !
    !--------------------------
    ! Loop on internal points:
    !--------------------------
    !
    do i = 1,bn_num
      !
      eps = 0.0_8
      sig = 0.0_8
      !
      cp = bn_list(i)
      !
      !----------------------------
      ! Loop on boundary elements:
      !----------------------------
      !
      do j = 1,be_num
        !
        be = be_list(j)
        !
        allocate(hs(3,2*be%nodes_num))
        allocate(gs(3,2*be%nodes_num))
        !
        call be_contains_node(be,cp,sing,lcn)
        !
        if (sing) then
          !Singular integration:
          call stress_sing_intgr(cp%point,be,lcn,hs,gs)
        else
          !Regular integration:
          call stress_reg_integr(psk,mu,nu,cp%point,be,hs,gs)
        end if
        !
        do k = 1,be%nodes_num
          u(2*k-1:2*k) = bn_list(be%incid(k))%disp
          t(2*k-1:2*k) = bn_list(be%incid(k))%trac
        end do
        !
        sig = sig - matmul(hs,u) + matmul(gs,t)
        !
        deallocate(hs); deallocate(gs)
        !
      end do
      !
      sig = 2*sig
      bn_list(i)%stress = sig
      !
      call get_strain_from_stress_elast(psk,nu,el,sig,eps)
      !
      bn_list(i)%strain = eps
      !
    end do
    !
    return
  end subroutine stress_strain_using_singular_integration
  !****************************************************
  !****************************************************
 end module fracture_2d_set_bound_results
