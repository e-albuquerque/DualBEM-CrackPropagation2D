module fracture_2d_assembler
!----------------------------------------------------
! A module with subroutines that mounts the final
! system of discrete equations, [A]{x} = {f}, for
! two-dimensional fracture problems.
!
! Date: August 2022
! Last update: -
!----------------------------------------------------
use fracture_2d_model
use fracture_2d_regular_integration
use fracture_2d_weakly_sing_integration
use fracture_2d_rad_exp_intgr
!
implicit none
!
real(kind=8),allocatable,private :: initial_h(:,:),initial_g(:,:),initial_f_vec(:)
!
contains
  subroutine assembler(alpha,a,f)
    !----------------------------------------------------
    ! Firstly mounts matrices [H] and [G] and, then,
    ! reorganizes according the boundary conditions to
    ! generate [A] and {f}.
    !
    ! Corners and traction discontinuities are treated by
    ! discontinuous elements or double nodes.
    !
    ! Date: 09 August 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: alpha ! It can be alpha q or alpha max
    real(kind=8),intent(out) :: a(2*bn_num,2*bn_num)
    real(kind=8),intent(out) :: f(2*bn_num)
    real(kind=8) :: h(2*bn_num,2*bn_num),g(2*bn_num,2*bn_num)
    real(kind=8) :: f_vec(2*bn_num)
    real(kind=8) :: aux(2*bn_num)
    real(kind=8) :: nv(2)
    type(boundary_node) :: cp
    integer :: i,j,k
    !
    do i = 1,bn_num
        bn_list(i)%trac = bn_list(i)%trac*alpha
    end do
    !
    call get_matrices_H_and_G(h,g)
    !
    allocate(initial_h(2*bn_num,2*bn_num))
    allocate(initial_g(2*bn_num,2*bn_num))
    call set_initial_h(h)
    call set_initial_g(g)
    !
    !----------------------------------
    ! Applying the boundary conditions:
    !----------------------------------
    !
    do j = 1,bn_num
      !
      cp = bn_list(j)
      !
      ! Modifying matrices for normal/tangential b.c.:
      !
      if (cp%bc_dir.ne.0) then
        !
        nv = get_normal_to_a_boundary_point(j)
        k = 2*j-1
        !
        aux = h(:,k)*nv(1) + h(:,k+1)*nv(2)
        h(:,k+1) = h(:,k+1)*nv(1) - h(:,k)*nv(2)
        h(:,k) = aux
        !
        aux = g(:,k)*nv(1) + g(:,k+1)*nv(2)
        g(:,k+1) = g(:,k+1)*nv(1) - g(:,k)*nv(2)
        g(:,k) = aux
        !
      end if
      !
      ! Reorganizing matrices H and G according the b.c.:
      !
      do k = 1,2
        !
        select case (cp%bc(k))
        case(0)
          !Displacement component prescribed:
          f_vec(2*j+k-2) = cp%disp(k)
          aux = g(:,2*j+k-2)
          g(:,2*j+k-2) = -h(:,2*j+k-2)
          h(:,2*j+k-2) = -aux
        case(1)
          !Traction component prescribed:
          f_vec(2*j+k-2) = cp%trac(k)
        case(2)
          !Elastic support (Robin b.c.):
          f_vec(2*j+k-2) = cp%h2(k)
          h(:,2*j+k-2) = h(:,2*j+k-2) - g(:,2*j+k-2)*cp%h1(k)
        case default
          !Invalid boundary condition
          call error_message_int_stop("Undefined boundary condition at node: ",j)
        end select
        !
      end do
      !
    end do
    !
    !----------------------------------
    ! Mounting matrices [A] and {f}:
    !----------------------------------
    !
    f = matmul(g,f_vec)
    a = h
    !
    return
  end subroutine assembler
  !****************************************************
  !****************************************************
  subroutine get_matrices_H_and_G(h,g)
    !----------------------------------------------------
    ! Verify if the problem will be solved by displacement
    ! or traction BIE and chose the appropriate code.
    !
    ! Date: 12 August 2022
    ! Last update: ---
    !---------------------------------------------------
    real(kind=8),intent(inout) :: h(2*bn_num,2*bn_num)
    real(kind=8),intent(inout) :: g(2*bn_num,2*bn_num)
    real(kind=8),allocatable :: gi(:,:),hi(:,:)
    type(boundary_node) :: cp
    type(boundary_elem) :: be
    real(kind=8) :: normal_cp(2)
    integer :: b,c,d,i,j,k,l,ii,kk,ll,bie
    !
    h = 0.0_8
    g = 0.0_8
    !
    do i = 1,bn_num
      !
      cp = bn_list(i)
      bie = cp%bie
      ii = 2*i-1
      normal_cp = get_normal_to_a_boundary_point(i)
      !
      !----------------------------
      ! Loop on boundary elements:
      !----------------------------
      !
      do j = 1,be_num
        !
        be = be_list(j)
        !
        allocate(hi(2,2*be%nodes_num))
        allocate(gi(2,2*be%nodes_num))
        !
        hi = 0.0_8
        gi = 0.0_8
        if (bie.eq.1 .or. bie.eq.2) then
            call displ_bie_matrices_H_and_G(cp,be,hi,gi)
            else if (bie.eq.3) then
                call tract_bie_matrices_H_and_G(cp,be,normal_cp,hi,gi)
            else
                call error_message_stop("Incorrect Boundary Integral Equation for node.")
        end if
        !
        do k = 1,be%nodes_num
          kk = 2*be%incid(k)-1
          h(ii:ii+1,kk:kk+1) = h(ii:ii+1,kk:kk+1) + hi(:,2*k-1:2*k)
          g(ii:ii+1,kk:kk+1) = g(ii:ii+1,kk:kk+1) + gi(:,2*k-1:2*k)
        end do
        !
        deallocate(hi)
        deallocate(gi)
        !
    end do
    !
    ll = 0.0_8
    if (bie.eq.1) then
        !
        ! Diagonal terms of matrix H
        !
        h(ii,ii) = h(ii,ii) + 0.5_8
        h(ii+1,ii+1) = h(ii+1,ii+1) + 0.5_8
        else if (bie.eq.2) then
            !
            ! Diagonal terms of matrix H on the source point
            !
            h(ii,ii) = h(ii,ii) + 0.5_8
            h(ii+1,ii+1) = h(ii+1,ii+1) + 0.5_8
            !
            ! Diagonal terms of matrix H on the coincidental node
            !
            do c = 1,cr_num
                d = size(cr_list(c)%nodes_num(:,1))
                do b = 1,d
                    if (cr_list(c)%nodes_num(b,1) .eq. i) then
                        l = cr_list(c)%nodes_num(b,2)
                        ll = 2*l-1
                    end if
                end do
            end do
            h(ii,ll) = h(ii,ll) + 0.5_8
            h(ii+1,ll+1) = h(ii+1,ll+1) + 0.5_8
            !
            else if(bie.eq.3) then
                !
                !
                ! Diagonal terms of matrix G on the source point
                !
                g(ii,ii) = g(ii,ii) - 0.5_8
                g(ii+1,ii+1) = g(ii+1,ii+1) - 0.5_8
                !
                ! Diagonal terms of matrix H on the coincidental node
                !
                do c = 1,cr_num
                    d = size(cr_list(c)%nodes_num(:,1))
                    do b = 1,d
                        if (cr_list(c)%nodes_num(b,2) .eq. i) then
                            l = cr_list(c)%nodes_num(b,1)
                            ll = 2*l-1
                        end if
                    end do
                end do
                g(ii,ll) = g(ii,ll) + 0.5_8
                g(ii+1,ll+1) = g(ii+1,ll+1) + 0.5_8
                !
                else
                call error_message_stop("Incorrect Boundary Integral Equation for node.")
    end if
    end do
    !
    return
  end subroutine get_matrices_H_and_G
  !****************************************************
  !****************************************************
  subroutine displ_bie_matrices_H_and_G(cp,be,hi,gi)
    !----------------------------------------------------
    ! Mounts matrices [H] and [G] using displacement BIE.
    !
    ! Corners and traction discontinuities are treated by
    ! discontinuous elements or double nodes.
    !
    ! Diagonal sub-matrices in matrix [H] are indirectly
    ! calculated by application of a rigid body motion.
    !
    ! Date: 09 August 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),allocatable,intent(inout) :: gi(:,:),hi(:,:)
    type(boundary_node),intent(in) :: cp
    type(boundary_elem),intent(in) :: be
    real(kind=8) :: mu,nu
    integer :: psk
    integer :: lcn
    logical :: sing
    !
    !
    !----------------------------------------
    ! Analysis type and material properties:
    !----------------------------------------
    !
    call get_plane_state_key(psk)
    call get_shear_modulus(mu)
    call get_poisson_ratio(nu)
    !
    !--------------------------------------
    ! Loop on boundary collocation points:
    !--------------------------------------
    !
    call be_contains_node(be,cp,sing,lcn)
    !
    if (sing) then
        !Singular integration:
        call displ_bie_sing_integr(cp%point,be,lcn,hi)
        call displ_bie_weakly_sing_integr(psk,mu,nu,lcn,be,gi)
    else
        !Regular integration:
        call displ_bie_reg_integr(psk,mu,nu,cp%point,be,hi,gi)
    end if
    !
    return
  end subroutine displ_bie_matrices_H_and_G
  !****************************************************
  !****************************************************
  subroutine tract_bie_matrices_H_and_G(cp,be,normal_cp,hi,gi)
    !----------------------------------------------------
    ! Mounts matrices [H] and [G] using traction BIE.
    !
    ! Corners and traction discontinuities are treated
    ! by discontinuous elements.
    !
    ! Date: 09 August 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),allocatable,intent(inout) :: gi(:,:),hi(:,:)
    type(boundary_node),intent(in) :: cp
    type(boundary_elem),intent(in) :: be
    real(kind=8),intent(in) :: normal_cp(2)
    real(kind=8) :: mu,nu
    integer :: psk
    integer :: lcn
    logical :: sing
    !
    !----------------------------------------
    ! Analysis type and material properties:
    !----------------------------------------
    !
    call get_plane_state_key(psk)
    call get_shear_modulus(mu)
    call get_poisson_ratio(nu)
    !
    !--------------------------------------
    ! Loop on boundary collocation points:
    !--------------------------------------
    !
    call be_contains_node(be,cp,sing,lcn)
    !
    if (sing) then
        !Singular integration:
        call traction_bie_sing_intgr(cp%point,be,lcn,normal_cp,hi,gi)
    else
        !Regular integration:
        call traction_bie_reg_integr(psk,mu,nu,cp%point,be,normal_cp,hi,gi)
    end if
    !
    return
  end subroutine tract_bie_matrices_H_and_G
  !****************************************************
  !****************************************************
  subroutine set_initial_h(n)
    !----------------------------------------------------
    ! Set the initial matrix H, before propagating the
    ! crack.
    !
    ! Date: 28 October 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: n(2*bn_num,2*bn_num)
    !
    initial_h = n
    !
  end subroutine set_initial_h
  !****************************************************
  !****************************************************
  subroutine set_initial_g(n)
    !----------------------------------------------------
    ! Set the initial matrix G, before propagating the
    ! crack.
    !
    ! Date: 28 October 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: n(2*bn_num,2*bn_num)
    !
    initial_g = n
    !
  end subroutine set_initial_g
  !****************************************************
  !****************************************************
  subroutine set_initial_f_vec(n)
    !----------------------------------------------------
    ! Set the initial vector F, before propagating the
    ! crack.
    !
    ! Date: 28 October 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: n(2*bn_num)
    !
    initial_f_vec = n
    !
  end subroutine set_initial_f_vec
  !****************************************************
  !****************************************************
  subroutine get_initial_h(n)
    !----------------------------------------------------
    ! Get the initial matrix H, before propagating the
    ! crack.
    !
    ! Date: 28 October 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),allocatable,intent(inout) :: n(:,:)
    !
    n = initial_h
    !
  end subroutine get_initial_h
  !****************************************************
  !****************************************************
  subroutine get_initial_g(n)
    !----------------------------------------------------
    ! Get the initial matrix G, before propagating the
    ! crack.
    !
    ! Date: 28 October 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),allocatable,intent(inout) :: n(:,:)
    !
    n = initial_g
    !
  end subroutine get_initial_g
  !****************************************************
  !****************************************************
  subroutine get_initial_f_vec(n)
    !----------------------------------------------------
    ! Get the initial vector F, before propagating the
    ! crack.
    !
    ! Date: 28 October 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(out) :: n(2*bn_num)
    !
    n = initial_f_vec
    !
  end subroutine get_initial_f_vec
  !****************************************************
  !****************************************************
end module fracture_2d_assembler
