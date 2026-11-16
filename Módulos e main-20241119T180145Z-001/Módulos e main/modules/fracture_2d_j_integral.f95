module fracture_2d_j_integral
!----------------------------------------------------
! A module with subroutines that calculate the
! J-integral for two-dimensional fracture problems.
!
! Date: September 2022
! Last update: ---
!----------------------------------------------------
use fracture_2d_model
use fracture_2d_regular_integration
use fracture_2d_rad_exp_intgr
use kelvin_fund_sol_2d
use numerical_integration_points
!
implicit none
!
type :: j_integral_node
  type(point_2d) :: point !Cartesian coordinates of the node.
  real(kind=8) :: stress(3) !Stress components: {s11 s22 s12}
  real(kind=8) :: strain(3) !Strain components: {s11 s22 s12}
  real(kind=8) :: grad(2) !Displacement gradient components: {u1,1 u2,2 u1,2}
end type j_integral_node
!
integer,private :: np ! Number of points that the circular boundary will be divided for the J-integral evaluation
integer,private :: nbn ! Number of boundary nodes per element
integer,private :: nbe ! Number of boundary elements
real(kind=8),private :: lb1, lb2 ! Lambda 1 and lambda 2
integer,private :: initial_node_num ! Position of the node where the boundary begins in relation to the crack tip
integer,allocatable :: sym_node(:) !Matrix with the position of the symmetric node on a J-integral node list
!
real(kind=8),private :: pi = 3.1415926535897932_8
!
contains
  subroutine K_using_Simpson_rule()
    !----------------------------------------------------
    ! Calculate the J-integral by Simpson's rule direct
    ! integration.
    !
    !
    ! Date: 06 September 2022
    ! Last update: ---
    !----------------------------------------------------
    type(j_integral_node) :: cp,cp_sym
    real(kind=8) :: J_i,J_ii
    type(boundary_elem) :: last_elem
    type(boundary_node) :: initial_node,final_node
    real(kind=8) :: ang,rad,normal(2)
    type(j_integral_node),allocatable :: j_int_list(:)
    real(kind=8) :: grad_i(2),grad_ii(2),W_i_normal,W_ii_normal,t_i(2),t_ii(2)
    integer :: i,j,k,k_ii_sign
    !
    type(j_integral_node) :: cp2,cp_sym2
    real(kind=8) :: grad2_i(2),grad2_ii(2),W2_i_normal,W2_ii_normal,t2_i(2),t2_ii(2),intgr1_i,intgr1_ii,intgr2_i,intgr2_ii
    !
    !--------------------------
    ! Evaluating the crack tips:
    !--------------------------
    !
    do i = 1,cr_num
        do j = 1,cr_list(i)%tips_num
            J_i = 0.0_8
            J_ii = 0.0_8
            call get_initial_node_and_last_element(i,j,last_elem,initial_node,final_node,normal)
            allocate(j_int_list(np))
            call Simpson_circular_line(last_elem%bp_2,initial_node%point,np,j_int_list,ang,rad)
            do k = 1,np
                if (k.eq.1) then
                    call evaluate_first_and_last_node(j_int_list(k),initial_node,normal)
                    else if (k.eq.np) then
                        call evaluate_first_and_last_node(j_int_list(k),final_node,normal)
                        else
                            call evaluate_j_integral_points(j_int_list(k),normal)
                end if
            end do
            do k = 1,np-1
                cp = j_int_list(k)
                cp_sym = j_int_list(sym_node(k))
                call get_variables(cp,cp_sym,last_elem,normal,grad_i,grad_ii,W_i_normal,W_ii_normal,t_i,t_ii)
                cp2 = j_int_list(k+1)
                cp_sym2 = j_int_list(sym_node(k+1))
                call get_variables(cp2,cp_sym2,last_elem,normal,grad2_i,grad2_ii,W2_i_normal,W2_ii_normal,t2_i,t2_ii)
                intgr1_i = W_i_normal - t_i(1)*grad_i(1) - t_i(2)*grad_i(2)
                intgr2_i = W2_i_normal - t2_i(1)*grad2_i(1) - t2_i(2)*grad2_i(2)
                intgr1_ii = W_ii_normal - t_ii(1)*grad_ii(1) - t_ii(2)*grad_ii(2)
                intgr2_ii = W2_ii_normal - t2_ii(1)*grad2_ii(1) - t2_ii(2)*grad2_ii(2)
                J_i = J_i + (intgr1_i+intgr2_i)*ang*rad/2._8
                J_ii = J_ii + (intgr1_ii+intgr2_ii)*ang*rad/2._8
            end do
            call evaluate_sign_of_k_ii(initial_node,final_node,normal,k_ii_sign)
            call k_i_and_k_ii(i,j,k_ii_sign,J_i,J_ii)
            deallocate(j_int_list)
        end do
    end do
    !
    return
  end subroutine K_using_Simpson_rule
  !****************************************************
  !****************************************************
  subroutine K_using_Gauss_integration_points()
    !----------------------------------------------------
    ! Calculate the J-integral by using Gauss integration
    ! points.
    !
    !
    ! Date: 07 September 2022
    ! Last update: ---
    !----------------------------------------------------
    type(j_integral_node) :: cp,cp_sym
    real(kind=8) :: J_i,J_ii
    type(boundary_elem) :: last_elem
    type(boundary_node) :: initial_node,final_node
    real(kind=8) :: rad,normal(2)
    real(kind=8),allocatable :: ome(:)
    type(j_integral_node),allocatable :: j_int_list(:)
    real(kind=8) :: grad_i(2),grad_ii(2),W_i_normal,W_ii_normal,t_i(2),t_ii(2)
    integer :: i,j,k,k_ii_sign
    !
    !--------------------------
    ! Evaluating the crack tips:
    !--------------------------
    !
    do i = 1,cr_num
        do j = 1,cr_list(i)%tips_num
            J_i = 0.0_8
            J_ii = 0.0_8
            call get_initial_node_and_last_element(i,j,last_elem,initial_node,final_node,normal)
            allocate(j_int_list(np))
            call Gauss_circular_line(last_elem%bp_2,initial_node%point,np,j_int_list,rad,ome)
            do k = 1,np
                call evaluate_j_integral_points(j_int_list(k),normal)
                !print*,j_int_list(k)%point
            end do
            do k = 1,np
                cp = j_int_list(k)
                cp_sym = j_int_list(sym_node(k))
                call get_variables(cp,cp_sym,last_elem,normal,grad_i,grad_ii,W_i_normal,W_ii_normal,t_i,t_ii)
                J_i = J_i + (W_i_normal - t_i(1)*grad_i(1) - t_i(2)*grad_i(2))*pi*rad*ome(k)
                J_ii = J_ii + (W_ii_normal - t_ii(1)*grad_ii(1) - t_ii(2)*grad_ii(2))*pi*rad*ome(k)
            end do
            call evaluate_sign_of_k_ii(initial_node,final_node,normal,k_ii_sign)
            call k_i_and_k_ii(i,j,k_ii_sign,J_i,J_ii)
            deallocate(j_int_list)
        end do
    end do
    !
    return
  end subroutine K_using_Gauss_integration_points
  !****************************************************
  !****************************************************
  subroutine K_using_Gauss_integration_points_on_boundary()
    !----------------------------------------------------
    ! Calculate the J-integral by using Gauss integration
    ! points on a circular boundary, with elements.
    !
    !
    ! Date: 07 September 2022
    ! Last update: ---
    !----------------------------------------------------
    type(j_integral_node) :: cp,cp_sym
    real(kind=8) :: J_i,J_ii
    real(kind=8) :: grad_i(2),grad_ii(2),grad_i_t(2),grad_ii_t(2)
    real(kind=8) :: W_i,W_ii,W_i_t,W_ii_t
    real(kind=8) :: t_i(2),t_ii(2),sig_i_t(3),sig_ii_t(3),sig_i(3),sig_ii(3)
    type(boundary_elem) :: last_elem
    type(boundary_node) :: initial_node,final_node
    type(boundary_elem),allocatable :: el_list(:),be
    type(point_2d),allocatable :: gi_pts(:)
    real(kind=8) :: lb1,lb2,v_normal(2),v_normal_x,v_normal_y,jac,normal(2)
    real(kind=8),allocatable :: ome(:),eta(:),ni(:)
    type(j_integral_node),allocatable :: j_int_list(:)
    integer :: gn_num,k_ii_sign
    integer :: i,j,k,l,m,o,n,nbe,nbn
    !
    !--------------------------
    ! Evaluating the crack tips:
    !--------------------------
    !
    !
    do i = 1,cr_num
        do j = 1,cr_list(i)%tips_num
            J_i = 0.0_8
            J_ii = 0.0_8
            call get_initial_node_and_last_element(i,j,last_elem,initial_node,final_node,normal)
            call get_element_information(nbe,nbn,lb1,lb2)
            call circular_line_for_j_integral(last_elem%bp_2,initial_node%point,360._8,1,nbe,1._8,nbn,lb1,lb2,j_int_list,el_list)
            do k = 1,size(j_int_list(:))
                call evaluate_j_integral_points(j_int_list(k),normal)
            end do
            do n = 1,size(el_list(:)) !Elements loop
                be = el_list(n)
                !
                gn_num = be%nodes_num
                if (gn_num.eq.1) gn_num = 2
                allocate (gi_pts(gn_num))
                gi_pts(1) = be%bp_1; gi_pts(gn_num) = be%bp_2
                if (be%nodes_num.gt.2) then
                    do o = 2,(be%nodes_num-1)
                        gi_pts(o) = j_int_list(be%incid(o))%point
                    end do
                end if
                !
                allocate(eta(np))
                allocate(ome(np))
                call gauss_legendre_points_1d(np,eta,ome)
                !
                if (gn_num.eq.2) then
                    call jac_normal_1d(jac,v_normal,gi_pts,0._8,be%nodes_num-1)
                end if
                do l = 1,np
                    allocate(ni(be%nodes_num))
                    call intrp_func_1d(ni,eta(l),be%lamb_1,be%lamb_2,be%nodes_num-1)
                    ! Evaluate Jacobian and normal vector for non-straight elements:
                    if (gn_num.gt.2) then
                        call jac_normal_1d(jac,v_normal,gi_pts,eta(l),be%nodes_num-1)
                    end if
                    sig_i = 0.0_8
                    sig_ii = 0.0_8
                    W_i = 0.0_8
                    W_ii = 0.0_8
                    grad_i = 0.0_8
                    grad_ii = 0.0_8
                    do k = 1,be%nodes_num
                        cp = j_int_list(be%incid(k))
                        cp_sym = j_int_list(sym_node(be%incid(k)))
                        call get_variables_key_3(cp,cp_sym,grad_i_t,grad_ii_t,W_i_t,W_ii_t,sig_i_t,sig_ii_t)
                        do m = 1,2
                            grad_i(m) = grad_i(m) + grad_i_t(m)*ni(k)
                            grad_ii(m) = grad_ii(m) + grad_ii_t(m)*ni(k)
                        end do
                        do m = 1,3
                            sig_i(m) = sig_i(m) + sig_i_t(m)*ni(k)
                            sig_ii(m) = sig_ii(m) + sig_ii_t(m)*ni(k)
                        end do
                        W_i = W_i + W_i_t*ni(k)
                        W_ii = W_ii + W_ii_t*ni(k)
                    end do
                    v_normal_x = normal(2)*v_normal(1)-normal(1)*v_normal(2)
                    v_normal_y = normal(1)*v_normal(1)+normal(2)*v_normal(2)
                    t_i(1) = sig_i(1)*v_normal_x+sig_i(3)*v_normal_y
                    t_i(2) = sig_i(2)*v_normal_y+sig_i(3)*v_normal_x
                    t_ii(1) = sig_ii(1)*v_normal_x+sig_ii(3)*v_normal_y
                    t_ii(2) = sig_ii(2)*v_normal_y+sig_ii(3)*v_normal_x
                    J_i = J_i + (W_i*v_normal_x-t_i(1)*grad_i(1)-t_i(2)*grad_i(2))*jac*ome(l)
                    J_ii = J_ii + (W_ii*v_normal_x-t_ii(1)*grad_ii(1)-t_ii(2)*grad_ii(2))*jac*ome(l)
                    !
                    deallocate(ni)
            end do
            deallocate(ome)
            deallocate(eta)
            deallocate(gi_pts)
            end do
            call evaluate_sign_of_k_ii(initial_node,final_node,normal,k_ii_sign)
            call k_i_and_k_ii(i,j,k_ii_sign,J_i,J_ii)
            deallocate(j_int_list)
        end do
    end do
    !
    return
  end subroutine K_using_Gauss_integration_points_on_boundary
  !****************************************************
  !****************************************************
  subroutine get_initial_node_and_last_element(i,j,last_elem,initial_node,final_node,normal)
    !----------------------------------------------------
    ! Get the last element and the node where the circular
    ! boundary for J-integral begins.
    !
    !
    ! Date: 08 September 2022
    ! Last update: 22 September 2022
    !----------------------------------------------------
    integer,intent(in) :: i,j
    type(boundary_elem),intent(out) :: last_elem
    type(boundary_node),intent(out) :: initial_node,final_node
    real(kind=8),intent(out) :: normal(2)
    integer :: aux1,initial_node_pos,final_node_pos
    !
    initial_node_pos = 0
    final_node_pos = 0
    !
    if (j.eq.1) then
        aux1 = size(cr_list(i)%nodes_num(:,1))-initial_node_num+1
        initial_node_pos = cr_list(i)%nodes_num(aux1,2)
        final_node_pos = cr_list(i)%nodes_num(aux1,1)
        last_elem = be_list(cr_list(i)%ext_el_num(2))
        else if (j.eq.2) then
            initial_node_pos = cr_list(i)%nodes_num(initial_node_num,1)
            final_node_pos = cr_list(i)%nodes_num(initial_node_num,2)
            last_elem = be_list(cr_list(i)%ext_el_num(4))
    end if
    normal = get_normal_to_a_boundary_point(initial_node_pos)
    initial_node = bn_list(initial_node_pos)
    final_node = bn_list(final_node_pos)
    !
  end subroutine get_initial_node_and_last_element
  !****************************************************
  !****************************************************
  subroutine evaluate_sign_of_k_ii(initial_node,final_node,normal,k_ii_sign)
    !----------------------------------------------------
    ! Evaluate if K_ii sign is going to be positive or
    ! negative.
    !
    ! K_ii_sign = 0 - positive
    !             1 - negative
    !
    ! Date: 20 October 2022
    ! Last update: ---
    !----------------------------------------------------
    type(boundary_node),intent(in) :: initial_node,final_node
    real(kind=8) :: normal(2)
    integer,intent(out) :: k_ii_sign
    real(kind=8) :: u1_sup,u1_inf
    !
    u1_sup = normal(2)*final_node%disp(1)-normal(1)*final_node%disp(2)
    u1_inf = normal(2)*initial_node%disp(1)-normal(1)*initial_node%disp(2)
    !
    if (u1_sup .gt. u1_inf) then
        k_ii_sign = 0
        else
            k_ii_sign = 1
    end if
    !
  end subroutine evaluate_sign_of_k_ii
  !****************************************************
  !****************************************************
  subroutine k_i_and_k_ii(i,j,k_ii_sign,J_i,J_ii)
    !----------------------------------------------------
    ! Calculate k_i and k_ii for a specific crack and a
    ! specific crack tip
    !
    !
    ! Date: 08 September 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: i,j,k_ii_sign
    real(kind=8) :: J_i,J_ii
    real(kind=8) :: nu,E,Eb
    integer :: psk
    !
    call get_plane_state_key(psk)
    call get_poisson_ratio(nu)
    call get_young_modulus(E)
    Eb = 0.0
    if (psk.eq.1) then
        Eb = E
        else if (psk.eq.2) then
            Eb = E/(1-nu**2)
    end if
    cr_list(i)%ki(j) = sqrt(abs(J_i*Eb))
    if (k_ii_sign .eq. 0) then
        cr_list(i)%kii(j) = sqrt(abs(J_ii*Eb))
        else if (k_ii_sign .eq. 1) then
            cr_list(i)%kii(j) = -sqrt(abs(J_ii*Eb))
    end if
    !
  end subroutine k_i_and_k_ii
  !****************************************************
  !****************************************************
  subroutine evaluate_j_integral_points(cp,normal)
    !----------------------------------------------------
    ! Calculate traction, strains, stresses, displacement
    ! gradient and strain energy density at internal points
    ! used for J-integral evaluation.
    !
    ! Stress and displacement gradient components are evaluated
    ! by Integral Equations and strain components are
    ! evaluated by Hooke's Law.
    !
    ! Date: 06 September 2022
    ! Last update: 17 September 2022
    !----------------------------------------------------
    type(j_integral_node),intent(inout) :: cp
    real(kind=8),intent(in) :: normal(2)
    real(kind=8) :: sig(3),grad(4),eps(3)
    real(kind=8),allocatable :: hg(:,:),gg(:,:)
    real(kind=8),allocatable :: hs(:,:),gs(:,:)
    real(kind=8),allocatable :: u_temp(:),t_temp(:)
    real(kind=8) :: mu,nu,el
    real(kind=8) :: L(2,2),sig_rot(3),grad_rot(4),eps_rot(3)
    integer :: psk
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
    ! Evaluating the node:
    !--------------------------
    !
    sig = 0._8
    grad = 0._8
    eps = 0._8
    !----------------------------
    ! Loop on boundary elements:
    !----------------------------
    !
    do j = 1,be_num
        !
        be = be_list(j)
        !
        allocate(hg(4,2*be%nodes_num))
        allocate(hs(3,2*be%nodes_num))
        allocate(gg(4,2*be%nodes_num))
        allocate(gs(3,2*be%nodes_num))
        !
        allocate(u_temp(2*be%nodes_num))
        allocate(t_temp(2*be%nodes_num))
        !
        call reg_integr_j_integr_point_to_be(psk,mu,nu,cp%point,be,hg,gg,hs,gs)
        !
        do k = 1,be%nodes_num
          u_temp(2*k-1:2*k) = bn_list(be%incid(k))%disp
          t_temp(2*k-1:2*k) = bn_list(be%incid(k))%trac
        end do
        !
        grad = grad - matmul(hg,u_temp) + matmul(gg,t_temp)
        sig = sig - matmul(hs,u_temp) + matmul(gs,t_temp)
        !
        deallocate(hg); deallocate(gg)
        deallocate(hs); deallocate(gs)
        deallocate(u_temp); deallocate(t_temp)
        !
    end do
    !
    call get_strain_from_stress_elast(psk,nu,el,sig,eps)
    !
    L(1,1) = normal(2)
    L(1,2) = normal(1)
    L(2,1) = -normal(1)
    L(2,2) = normal(2)
    do i = 1,2
        sig_rot(i) = L(1,i)*L(1,i)*sig(1)+L(1,i)*L(2,i)*sig(3)+L(2,i)*L(1,i)*sig(3)+L(2,i)*L(2,i)*sig(2)
        eps_rot(i) = L(1,i)*L(1,i)*eps(1)+L(1,i)*L(2,i)*eps(3)+L(2,i)*L(1,i)*eps(3)+L(2,i)*L(2,i)*eps(2)
        grad_rot(i) = L(1,i)*L(1,i)*grad(1)+L(1,i)*L(2,i)*grad(3)+L(2,i)*L(1,i)*grad(4)+L(2,i)*L(2,i)*grad(2)
    end do
    sig_rot(3) = L(1,1)*L(1,2)*sig(1)+L(1,1)*L(2,2)*sig(3)+L(2,1)*L(1,2)*sig(3)+L(2,1)*L(2,2)*sig(2)
    eps_rot(3) = L(1,1)*L(1,2)*eps(1)+L(1,1)*L(2,2)*eps(3)+L(2,1)*L(1,2)*eps(3)+L(2,1)*L(2,2)*eps(2)
    grad_rot(3) = L(1,1)*L(1,2)*grad(1)+L(1,1)*L(2,2)*grad(3)+L(2,1)*L(1,2)*grad(4)+L(2,1)*L(2,2)*grad(2)
    grad_rot(4) = L(1,2)*L(1,1)*grad(1)+L(1,2)*L(2,1)*grad(3)+L(2,2)*L(1,1)*grad(4)+L(2,2)*L(2,1)*grad(2)
    !
    cp%grad(1) = grad_rot(1)
    cp%grad(2) = grad_rot(4)
    cp%stress = sig_rot
    cp%strain = eps_rot
    !
    return
  end subroutine evaluate_j_integral_points
  !****************************************************
  !****************************************************
  subroutine evaluate_first_and_last_node(cp,node,normal)
    !----------------------------------------------------
    ! Calculate stresses and displacement gradient at
    ! boundary nodes used for J-integral evaluation.
    !
    ! Stress and displacement gradient components are evaluated
    ! by Integral Equations.
    !
    ! Date: 22 September 2022
    ! Last update: ---
    !----------------------------------------------------
    type(j_integral_node),intent(inout) :: cp
    type(boundary_node),intent(in) :: node
    real(kind=8),intent(in) :: normal(2)
    real(kind=8) :: sig(3),grad(4),eps(3)
    real(kind=8),allocatable :: hg(:,:),gg(:,:),hi(:,:),gi(:,:)
    real(kind=8),allocatable :: u_temp(:),t_temp(:)
    real(kind=8) :: mu,nu,el
    real(kind=8) :: L(2,2),sig_rot(3),grad_rot(4),eps_rot(3)
    integer :: psk
    type(boundary_elem) :: be
    integer :: i,j,k,lcn
    logical :: sing
    !
    allocate(hg(4,2*bn_num))
    allocate(gg(4,2*bn_num))
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
    ! Evaluating the node:
    !--------------------------
    !
    grad = 0._8
    !----------------------------
    ! Loop on boundary elements:
    !----------------------------
    !
    !
    do j = 1,be_num
        !
        be = be_list(j)
        call be_contains_node(be,node,sing,lcn)
        !
        allocate(hi(4,2*be%nodes_num))
        allocate(gi(4,2*be%nodes_num))
        allocate(u_temp(2*be%nodes_num))
        allocate(t_temp(2*be%nodes_num))
        !
        !
        if (sing) then
            call displ_grad_sing_intgr(node%point,be,lcn,hi,gi)
            else
                call reg_integr_grad_displ_to_be_point(psk,mu,nu,node%point,be,hi,gi)
        end if
        !
        do k = 1,be%nodes_num
          u_temp(2*k-1:2*k) = bn_list(be%incid(k))%disp
          t_temp(2*k-1:2*k) = bn_list(be%incid(k))%trac
        end do
        !
        grad = grad - matmul(hi,u_temp) + matmul(gi,t_temp)
        !
        deallocate(hi)
        deallocate(gi)
        deallocate(u_temp)
        deallocate(t_temp)
        !
        !
    end do
    sig = node%stress
    eps = node%strain
    grad = 2*grad
    !
    L(1,1) = normal(2)
    L(1,2) = normal(1)
    L(2,1) = -normal(1)
    L(2,2) = normal(2)
    do i = 1,2
        sig_rot(i) = L(1,i)*L(1,i)*sig(1)+L(1,i)*L(2,i)*sig(3)+L(2,i)*L(1,i)*sig(3)+L(2,i)*L(2,i)*sig(2)
        eps_rot(i) = L(1,i)*L(1,i)*eps(1)+L(1,i)*L(2,i)*eps(3)+L(2,i)*L(1,i)*eps(3)+L(2,i)*L(2,i)*eps(2)
        grad_rot(i) = L(1,i)*L(1,i)*grad(1)+L(1,i)*L(2,i)*grad(3)+L(2,i)*L(1,i)*grad(4)+L(2,i)*L(2,i)*grad(2)
    end do
    sig_rot(3) = L(1,1)*L(1,2)*sig(1)+L(1,1)*L(2,2)*sig(3)+L(2,1)*L(1,2)*sig(3)+L(2,1)*L(2,2)*sig(2)
    eps_rot(3) = L(1,1)*L(1,2)*eps(1)+L(1,1)*L(2,2)*eps(3)+L(2,1)*L(1,2)*eps(3)+L(2,1)*L(2,2)*eps(2)
    grad_rot(3) = L(1,1)*L(1,2)*grad(1)+L(1,1)*L(2,2)*grad(3)+L(2,1)*L(1,2)*grad(4)+L(2,1)*L(2,2)*grad(2)
    grad_rot(4) = L(1,2)*L(1,1)*grad(1)+L(1,2)*L(2,1)*grad(3)+L(2,2)*L(1,1)*grad(4)+L(2,2)*L(2,1)*grad(2)
    cp%grad(1) = grad_rot(1)
    cp%grad(2) = grad_rot(4)
    cp%stress = sig_rot
    cp%strain = eps_rot
    !
    return
  end subroutine evaluate_first_and_last_node
  !****************************************************
  !****************************************************
  subroutine get_variables (cp,cp_sym,last_elem,ni,grad_i,grad_ii,W_i_normal,W_ii_normal,t_i,t_ii)
    !----------------------------------------------------
    ! Get the variables needed to calculate the J-integral.
    !
    !
    ! Date: 06 September 2022
    ! Last update: ---
    !----------------------------------------------------
    type(j_integral_node),intent(in) :: cp, cp_sym
    type(boundary_elem),intent(in) :: last_elem
    real(kind=8) :: ni(2)
    real(kind=8) :: sig_i(3),sig_ii(3),eps_i(3),eps_ii(3)
    real(kind=8),intent(out) :: grad_i(2),grad_ii(2),W_i_normal,W_ii_normal,t_i(2),t_ii(2)
    real(kind=8) :: W_i,W_ii,normal(2),r(3)
    real(kind=8) :: r_x,r_y
    !
    sig_i (1) = (cp%stress(1)+cp_sym%stress(1))/2._8
    sig_i (2) = (cp%stress(2)+cp_sym%stress(2))/2._8
    sig_i (3) = (cp%stress(3)-cp_sym%stress(3))/2._8
    !
    sig_ii (1) = (cp%stress(1)-cp_sym%stress(1))/2._8
    sig_ii (2) = (cp%stress(2)-cp_sym%stress(2))/2._8
    sig_ii (3) = (cp%stress(3)+cp_sym%stress(3))/2._8
    !
    eps_i (1) = (cp%strain(1)+cp_sym%strain(1))/2._8
    eps_i (2) = (cp%strain(2)+cp_sym%strain(2))/2._8
    eps_i (3) = (cp%strain(3)-cp_sym%strain(3))/2._8
    !
    eps_ii (1) = (cp%strain(1)-cp_sym%strain(1))/2._8
    eps_ii (2) = (cp%strain(2)-cp_sym%strain(2))/2._8
    eps_ii (3) = (cp%strain(3)+cp_sym%strain(3))/2._8
    !
    grad_i (1) = (cp%grad(1)+cp_sym%grad(1))/2._8
    grad_i (2) = (cp%grad(2)-cp_sym%grad(2))/2._8
    !
    grad_ii (1) = (cp%grad(1)-cp_sym%grad(1))/2._8
    grad_ii (2) = (cp%grad(2)+cp_sym%grad(2))/2._8
    !
    !
    W_i = (sig_i(1)*eps_i(1)+sig_i(2)*eps_i(2)+2*sig_i(3)*eps_i(3))/2._8
    W_ii = (sig_ii(1)*eps_ii(1)+sig_ii(2)*eps_ii(2)+2*sig_ii(3)*eps_ii(3))/2._8
    !
    r = distance_2d(last_elem%bp_2,cp%point)
    r_x = ni(2)*r(2)-ni(1)*r(3)
    r_y = ni(1)*r(2)+ni(2)*r(3)
    normal(1) = r_x/r(1)
    normal(2) = r_y/r(1)
    W_i_normal = W_i*normal(1)
    W_ii_normal = W_ii*normal(1)
    !
    t_i(1) = sig_i(1)*normal(1) + sig_i(3)*normal(2)
    t_i(2) = sig_i(2)*normal(2) + sig_i(3)*normal(1)
    !
    t_ii(1) = sig_ii(1)*normal(1) + sig_ii(3)*normal(2)
    t_ii(2) = sig_ii(2)*normal(2) + sig_ii(3)*normal(1)
    !
    !
    return
  end subroutine get_variables
  !****************************************************
  !****************************************************
  subroutine get_variables_key_3 (cp,cp_sym,grad_i,grad_ii,W_i,W_ii,sig_i,sig_ii)
    !----------------------------------------------------
    ! Get the variables needed to calculate the J-integral.
    !
    !
    ! Date: 15 September 2022
    ! Last update: ---
    !----------------------------------------------------
    type(j_integral_node),intent(in) :: cp, cp_sym
    real(kind=8),intent(out) :: grad_i(2),grad_ii(2),W_i,W_ii,sig_i(3),sig_ii(3)
    real(kind=8) :: eps_i(3),eps_ii(3)
    real(kind=8) :: el,nu
    integer :: psk
    !
    sig_i (1) = (cp%stress(1)+cp_sym%stress(1))/2._8
    sig_i (2) = (cp%stress(2)+cp_sym%stress(2))/2._8
    sig_i (3) = (cp%stress(3)-cp_sym%stress(3))/2._8
    !
    sig_ii (1) = (cp%stress(1)-cp_sym%stress(1))/2._8
    sig_ii (2) = (cp%stress(2)-cp_sym%stress(2))/2._8
    sig_ii (3) = (cp%stress(3)+cp_sym%stress(3))/2._8
    !
    grad_i (1) = (cp%grad(1)+cp_sym%grad(1))/2._8
    grad_i (2) = (cp%grad(2)-cp_sym%grad(2))/2._8
    !
    grad_ii (1) = (cp%grad(1)-cp_sym%grad(1))/2._8
    grad_ii (2) = (cp%grad(2)+cp_sym%grad(2))/2._8
    !
    call get_plane_state_key(psk)
    call get_young_modulus(el)
    call get_poisson_ratio(nu)
    call get_strain_from_stress_elast(psk,nu,el,sig_i,eps_i)
    call get_strain_from_stress_elast(psk,nu,el,sig_ii,eps_ii)
    !
    W_i = (sig_i(1)*eps_i(1)+sig_i(2)*eps_i(2)+2*sig_i(3)*eps_i(3))/2._8
    W_ii = (sig_ii(1)*eps_ii(1)+sig_ii(2)*eps_ii(2)+2*sig_ii(3)*eps_ii(3))/2._8
    !
    return
  end subroutine get_variables_key_3
  !****************************************************
  !****************************************************
  subroutine Simpson_circular_line(p1,p2,nip,j_int_list,d_theta,r)
    !----------------------------------------------------
    ! Generate points along a circular line that will be
    ! used to the J-integral evaluation.
    !
    ! p1 = the circumference central point
    ! p2 = the line initial point
    ! p3 = the line final point
    ! nip = the number of points to be generated
    ! j_int_list = list of points for the J-integral
    ! d_theta = angle that the boundary was divided in
    ! rad = radius
    !
    ! Date: 06 September 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: nip
    type(point_2d), intent(in) :: p1,p2
    type(j_integral_node),intent(inout) :: j_int_list(:)
    real(kind=8),intent(out) :: d_theta,r
    type(point_2d),allocatable :: pts(:)
    real(kind=8) :: ang_rad
    integer :: i,n_seg
    real(kind=8) :: rad(3)
    real(kind=8) :: theta
    !
    n_seg = nip-1
    !
    ang_rad = pi*360._8/180.0_8
    !
    if (allocated(pts)) deallocate(pts)
    allocate(pts(nip))
    !
    rad = distance_2d(p1,p2)
    r = rad(1)
    theta = atan2(rad(3),rad(2)) !Initial angle
    !
    pts(1) = p2
    d_theta = ang_rad/(n_seg)
    theta = theta + d_theta
    do i = 2,nip
      pts(i)%x = p1%x + rad(1)*cos(theta)
      pts(i)%y = p1%y + rad(1)*sin(theta)
      if (i.eq.nip) exit
      theta = theta + d_theta
    end do
    !
    if (allocated(sym_node)) deallocate(sym_node)
    allocate(sym_node(nip))
    !
    do i = 1,nip
      j_int_list(i)%point = pts(i)
      j_int_list(i)%stress = 0.0
      j_int_list(i)%grad = 0.0
      sym_node(i) = nip + 1 - i
    end do
    !
    return
  end subroutine Simpson_circular_line
  !****************************************************
  !****************************************************
  subroutine Gauss_circular_line(p1,p2,nip,j_int_list,r,ome)
    !----------------------------------------------------
    ! Generate points along a circular line that will be
    ! used to the J-integral evaluation using Gauss integration
    ! points.
    !
    ! p1 = the circumference central point
    ! p2 = the line initial point
    ! p3 = the line final point
    ! nip = the number of integration points that will be used
    ! j_int_list = list of points for the J-integral
    !
    ! Date: 07 September 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: nip
    type(point_2d), intent(in) :: p1,p2
    type(j_integral_node),intent(inout) :: j_int_list(:)
    real(kind=8),intent(out) :: r
    real(kind=8),allocatable,intent(out) :: ome(:)
    type(point_2d),allocatable :: pts(:)
    real(kind=8),allocatable :: eta(:)
    integer :: i
    real(kind=8) :: rad(3)
    real(kind=8) :: theta,theta_initial
    !
    !
    !
    if (allocated(pts)) deallocate(pts)
    allocate(pts(nip))
    !
    rad = distance_2d(p1,p2)
    r = rad(1)
    theta_initial = atan2(rad(3),rad(2)) !Initial angle
    !
    allocate (eta(nip))
    allocate (ome(nip))
    call gauss_legendre_points_1d(nip,eta,ome)
    !
    do i = 1,nip
        theta = theta_initial + pi*(eta(i)+1)
        pts(i)%x = p1%x + rad(1)*cos(theta)
        pts(i)%y = p1%y + rad(1)*sin(theta)
        if (i.eq.nip) exit
    end do
    !
    if (allocated(sym_node)) deallocate(sym_node)
    allocate(sym_node(nip))
    !
    do i = 1,nip
      j_int_list(i)%point = pts(i)
      j_int_list(i)%stress = 0.0
      j_int_list(i)%grad = 0.0
      if (mod(nip,2) == 0) then
        if (mod(i,2) == 1) then
            sym_node(i) = i+1
            else
                sym_node(i) = i-1
        end if
        else if (mod(nip,2) == 1) then
            if (i .eq. 1) then
                sym_node(i) = i
                else
                    if (mod(i,2) == 0) then
                        sym_node(i) = i+1
                        else
                            sym_node(i) = i-1
                    end if
            end if
      end if
    end do
    !
    return
  end subroutine Gauss_circular_line
  !****************************************************
  !****************************************************
  subroutine circular_line_for_j_integral(p1,p2,ang,dk,nbe,scl,nbn,lb1,lb2,j_int_list,el_list)
    !----------------------------------------------------
    ! Generate elements and boundary nodes for a circular
    ! line.
    !
    ! The extreme elements (the first and last element in
    ! the line) are semi-discontinuous elements with the
    ! discontinuous parameters lb1 and lb2, respectively.
    !
    ! Date: 27 August 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: dk,nbe,nbn
    real(kind=8),intent(in) :: lb1,lb2,ang,scl
    type(point_2d),intent(in) :: p1,p2
    type(boundary_elem),allocatable,intent(out) :: el_list(:)
    type(j_integral_node),allocatable,intent(out) :: j_int_list(:)
    real(kind=8) :: ang_rad
    logical :: cc
    type(point_2d),allocatable :: pts(:)
    integer :: i
    !
    !
    !-------------------------------------------------------------
    ! Evaluating a list of points to generate elements and nodes:
    !-------------------------------------------------------------
    !
    if (dk.lt.0) then
      cc = .false. !The line goes clockwise direction
    else
      cc = .true. !The line goes counter-clockwise direction
    end if
    !
    ang_rad = pi*ang/180.0_8
    !
    !
    if (nbn.lt.3) then
      call divide_circular_line_2d_scaling(p1,p2,ang_rad,scl,nbe,pts,cc)
    else
      call divide_circular_line_2d_semi_scaling(p1,p2,ang_rad,scl,nbe,nbn-1,pts,cc)
    end if
    !
    !
    !--------------------------------
    ! Generating elements and nodes:
    !--------------------------------
    !
    !
    select case (nbn)
    case(1)
      call constant_elements_for_j_integral(pts,j_int_list,el_list)
    case(2)
      call linear_elements_for_j_integral(pts,lb1,lb2,j_int_list,el_list)
    case(3)
      call quadratic_elements_for_j_integral(pts,lb1,lb2,j_int_list,el_list)
    case default
      call error_message_stop("Undefined element order.")
    end select
    !
    !
    !------------------------------------
    ! Checking for closed circumference:
    !------------------------------------
    !
    if (allocated(sym_node)) deallocate(sym_node)
    allocate(sym_node(size(j_int_list(:))))
    do i = 1,size(j_int_list(:))
      j_int_list(i)%stress = 0.0
      j_int_list(i)%grad = 0.0
      sym_node(i) = size(j_int_list(:))-i+1
    end do
    !
    return
  end subroutine circular_line_for_j_integral
  !****************************************************
  !****************************************************
  subroutine constant_elements_for_j_integral(pts,j_int_list,el_list)
    !----------------------------------------------------
    ! Generate constant elements for a given list of
    ! points.
    !
    ! Date: 07 September 2022
    ! Last update: ---
    !----------------------------------------------------
    type(point_2d),intent(in) :: pts(:)
    type(j_integral_node),allocatable,intent(out) :: j_int_list(:)
    type(boundary_elem),allocatable,intent(out) :: el_list(:)
    type(point_2d) :: p1,p2
    type(j_integral_node) :: ljn(10000)
    type(boundary_elem) :: lel(10000)
    integer :: i,tnel,tnjn
    !
    !--------------------------------------------------------------
    ! Evaluating the parametric coordinate for each boundary node:
    !--------------------------------------------------------------
    !
    tnel = 0
    tnjn = 0
    !
    !---------------------------------------------
    ! Generating the boundary elements and nodes:
    !---------------------------------------------
    !
    do i = 1,size(pts)-1
      !
      p1 = pts(i)
      p2 = pts(i+1)
      !
      tnel = tnel + 1 !Adding a boundary element to the model
      !
      lel(tnel)%bp_1 = p1
      lel(tnel)%bp_2 = p2
      allocate(lel(tnel)%incid(1))
      lel(tnel)%incid(1) = tnel
      lel(tnel)%lamb_1 = 1.0_8
      lel(tnel)%lamb_2 = 1.0_8
      lel(tnel)%nodes_num = 1
      !
      tnjn = tnjn + 1 !Adding a boundary node to the model
      !
      ljn(tnjn)%point = midpoint_2d(p1,p2)
      !
    end do
    !
    allocate(el_list(tnel))
    allocate(j_int_list(tnjn))
    el_list = lel(1:tnel)
    j_int_list = ljn(1:tnjn)
    !
    return
  end subroutine constant_elements_for_j_integral
  !****************************************************
  !****************************************************
  subroutine linear_elements_for_j_integral(pts,lb1,lb2,j_int_list,el_list)
    !----------------------------------------------------
    ! Generate linear elements for a given list of
    ! points.
    !
    ! Date: 07 September 2022
    ! Last update: ---
    !----------------------------------------------------
    type(point_2d),intent(in) :: pts(:)
    real(kind=8),intent(in) :: lb1,lb2
    type(j_integral_node),allocatable,intent(out) :: j_int_list(:)
    type(boundary_elem),allocatable,intent(out) :: el_list(:)
    type(point_2d) :: p1,p2
    type(point_2d) :: el_pts(2)
    type(j_integral_node) :: ljn(10000)
    type(boundary_elem) :: lel(10000)
    integer :: i,tnel,tnjn
    !
    !--------------------------------------------------------------
    ! Evaluating the parametric coordinate for each boundary node:
    !--------------------------------------------------------------
    !
    tnel = 0
    tnjn = 0
    !
    !---------------------------------------------
    ! Generating the boundary elements and nodes:
    !---------------------------------------------
    !
    do i = 1,size(pts)-1
      !
      p1 = pts(i)
      p2 = pts(i+1)
      !
      tnel = tnel + 1 !Adding a boundary element to the model
      !
      lel(tnel)%bp_1 = p1
      lel(tnel)%bp_2 = p2
      lel(tnel)%nodes_num = 2
      allocate(lel(tnel)%incid(2))
      if (i.eq.1) then
        lel(tnel)%lamb_1 = lb1
      else
        lel(tnel)%lamb_1 = 1.0_8
      end if
      if (i.eq.size(pts)-1) then
        lel(tnel)%lamb_2 = lb2
      else
        lel(tnel)%lamb_2 = 1.0_8
      end if
      !
      if (i.eq.1) then
        tnjn = tnjn + 1 !Adding a boundary node to the model
        el_pts(1) = p1
        el_pts(2) = p2
        ljn(tnjn)%point = cartesian_point_from_intrp_func_1d(el_pts,-lb1)
      end if
      !
      lel(tnel)%incid(1) = tnjn
      !
      tnjn = tnjn + 1 !Adding a boundary node to the model
      if (i.eq.size(pts)-1) then
        el_pts(1) = p1
        el_pts(2) = p2
        ljn(tnjn)%point = cartesian_point_from_intrp_func_1d(el_pts,+lb2)
      else
        ljn(tnjn)%point = p2
      end if
      lel(tnel)%incid(2) = tnjn
      !
    end do
    !
    allocate(el_list(tnel))
    allocate(j_int_list(tnjn))
    el_list = lel(1:tnel)
    j_int_list = ljn(1:tnjn)
    !
    return
  end subroutine linear_elements_for_j_integral
  !****************************************************
  !****************************************************
  subroutine quadratic_elements_for_j_integral(pts,lb1,lb2,j_int_list,el_list)
    !----------------------------------------------------
    ! Generate quadratic elements for a given list of
    ! points.
    !
    ! Date: 07 September 2022
    ! Last update: ---
    !----------------------------------------------------
    type(point_2d),intent(in) :: pts(:)
    real(kind=8),intent(in) :: lb1,lb2
    type(j_integral_node),allocatable,intent(out) :: j_int_list(:)
    type(boundary_elem),allocatable,intent(out) :: el_list(:)
    type(point_2d) :: p1,p2,p3
    type(point_2d) :: el_pts(3)
    type(j_integral_node) :: ljn(10000)
    type(boundary_elem) :: lel(10000)
    integer :: i,n,tnel,tnjn
    !
    !--------------------------------------------------------------
    ! Evaluating the parametric coordinate for each boundary node:
    !--------------------------------------------------------------
    !
    !
    tnel = 0
    tnjn = 0
    !---------------------------------------------
    ! Generating the boundary elements and nodes:
    !---------------------------------------------
    !
    n = int(0.5*(size(pts)-1))
    !
    do i = 1,n
      !
      p1 = pts(2*i-1)
      p2 = pts(2*i)
      p3 = pts(2*i+1)
      !
      tnel = tnel + 1 !Adding a boundary element to the model
      !
      lel(tnel)%bp_1 = p1
      lel(tnel)%bp_2 = p3
      lel(tnel)%nodes_num = 3
      allocate(lel(tnel)%incid(3))
      if (i.eq.1) then
        lel(tnel)%lamb_1 = lb1
      else
        lel(tnel)%lamb_1 = 1.0_8
      end if
      if (i.eq.n) then
        lel(tnel)%lamb_2 = lb2
      else
        lel(tnel)%lamb_2 = 1.0_8
      end if
      !
      if (i.eq.1) then
        tnjn = tnjn + 1 !Adding a boundary node to the model
        el_pts(1) = p1
        el_pts(2) = p2
        el_pts(3) = p3
        ljn(tnjn)%point = cartesian_point_from_intrp_func_1d(el_pts,-lb1)
      end if
      !
      lel(tnel)%incid(1) = tnjn
      !
      tnjn = tnjn + 1 !Adding a boundary node to the model
      lel(tnel)%incid(2) = tnjn
      ljn(tnjn)%point = p2
      !
      tnjn = tnjn + 1 !Adding a boundary node to the model
      if (i.eq.n) then
        el_pts(1) = p1
        el_pts(2) = p2
        el_pts(3) = p3
        ljn(tnjn)%point = cartesian_point_from_intrp_func_1d(el_pts,+lb2)
      else
        ljn(tnjn)%point = p3
      end if
      lel(tnel)%incid(3) = tnjn
      !
    end do
    !
    allocate(el_list(tnel))
    allocate(j_int_list(tnjn))
    el_list = lel(1:tnel)
    j_int_list = ljn(1:tnjn)
    !
    return
  end subroutine quadratic_elements_for_j_integral
  !****************************************************
  !****************************************************
  subroutine set_number_of_points(k)
    !----------------------------------------------------
    ! Set the number of points that the circular boundary
    ! will be divided for the J-integral evaluation
    !
    !
    ! Date: 06 September 2022
    ! Last update: 26 September 2022
    !----------------------------------------------------
    integer,intent(in) :: k
    !
    np = k
    !
    return
  end subroutine set_number_of_points
  !****************************************************
  !****************************************************
  subroutine get_number_of_points(k)
    !----------------------------------------------------
    ! Get the number of points that the circular boundary
    ! will be divided for the J-integral evaluation.
    !
    ! Date: 06 September 2022
    ! Last update: 26 September 2022
    !----------------------------------------------------
    integer,intent(out) :: k
    !
    k = np
    !
    return
  end subroutine get_number_of_points
  !****************************************************
  !****************************************************
  subroutine set_element_information(el,nd,lbd1,lbd2)
    !----------------------------------------------------
    ! Set the number of elements, number of nodes per element,
    ! lambda 1 e lambda 2.
    !
    ! Date: 07 September 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: el,nd
    real(kind=8),intent(in) :: lbd1,lbd2
    !
    nbe = el
    nbn = nd
    lb1 = lbd1
    lb2 = lbd2
    !
    return
  end subroutine set_element_information
  !****************************************************
  !****************************************************
  subroutine get_element_information(el,nd,lbd1,lbd2)
    !----------------------------------------------------
    ! Get the number of elements, number of nodes per element,
    ! lambda 1 e lambda 2.
    !
    ! Date: 07 September 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(out) :: el,nd
    real(kind=8),intent(out) :: lbd1,lbd2
    !
    el = nbe
    nd = nbn
    lbd1 = lb1
    lbd2 = lb2
    !
    return
  end subroutine get_element_information
  !****************************************************
  !****************************************************
  subroutine set_initial_node_num(k)
    !----------------------------------------------------
    ! Set the position of the node where the circular
    ! boundary begins in relation to the crack tip.
    !
    !
    ! Date: 06 September 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: k
    !
    initial_node_num = k
    !
    return
  end subroutine set_initial_node_num
  !****************************************************
  !****************************************************
  subroutine get_initial_node_num(k)
    !----------------------------------------------------
    ! Get the position of the node where the circular
    ! boundary begins in relation to the crack tip.
    !
    ! Date: 06 September 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(out) :: k
    !
    k = initial_node_num
    !
    return
  end subroutine get_initial_node_num
  !****************************************************
  !****************************************************
end module fracture_2d_j_integral
