module fracture_2d_internal_points
!----------------------------------------------------
! A module with subroutines that evaluate internal
! fields for two-dimensional fracture problems.
!
! Date: July 2022
! Last update: ---
!----------------------------------------------------
use fracture_2d_model
use fracture_2d_regular_integration
use kelvin_fund_sol_2d
use numerical_integration_points
!
implicit none
!
contains
  subroutine evaluate_internal_points()
    !----------------------------------------------------
    ! Calculate displacements, strains and stresses
    ! at internal points for two-dimensional isotropic
    ! linear elasticity problems.
    !
    ! Displacement and Stress components are evaluated
    ! by Integral Equations and strain components are
    ! evaluated by Hooke's Law.
    !
    ! Date: 27 June 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8) :: u(2),eps(3),sig(3)
    real(kind=8),allocatable :: h(:,:),g(:,:)
    real(kind=8),allocatable :: hs(:,:),gs(:,:)
    real(kind=8),allocatable :: u_temp(:),t_temp(:)
    real(kind=8) :: mu,nu,el
    integer :: psk
    type(internal_node) :: cp
    type(boundary_elem) :: be
    integer :: i,j,k
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
    do i = 1,in_num
      !
      u = 0.0_8
      eps = 0.0_8
      sig = 0.0_8
      !
      cp = in_list(i)
      !
      !----------------------------
      ! Loop on boundary elements:
      !----------------------------
      !
      do j = 1,be_num
        !
        be = be_list(j)
        !
        allocate(h(2,2*be%nodes_num))
        allocate(hs(3,2*be%nodes_num))
        allocate(g(2,2*be%nodes_num))
        allocate(gs(3,2*be%nodes_num))
        !
        allocate(u_temp(2*be%nodes_num))
        allocate(t_temp(2*be%nodes_num))
        !
        call reg_integr_int_point_to_be(psk,mu,nu,cp%point,be,h,g,hs,gs)
        !
        do k = 1,be%nodes_num
          u_temp(2*k-1:2*k) = bn_list(be%incid(k))%disp
          t_temp(2*k-1:2*k) = bn_list(be%incid(k))%trac
        end do
        !
        u = u - matmul(h,u_temp) + matmul(g,t_temp)
        sig = sig - matmul(hs,u_temp) + matmul(gs,t_temp)
        !
        deallocate(h); deallocate(g)
        deallocate(hs); deallocate(gs)
        deallocate(u_temp); deallocate(t_temp)
        !
      end do
      !
      in_list(i)%disp = u
      in_list(i)%stress = sig
      !
      call get_strain_from_stress_elast(psk,nu,el,sig,eps)
      !
      in_list(i)%strain = eps
      !
    end do
    !
    return
  end subroutine evaluate_internal_points
  !****************************************************
  !****************************************************
end module fracture_2d_internal_points
