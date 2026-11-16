module fracture_2d_displ_correlation
!----------------------------------------------------
! A module with subroutines that calculates the stress
! intensity factors(Ki e Kii) by the displacement
! correlation method.
!
! Date: September 2022
! Last update: -
!----------------------------------------------------
use fracture_2d_model
!
implicit none
!
real(kind=8),private :: two_pi = 6.283185307179586_8
integer,private :: nd ! Node for evaluation of the displacement correlation method
!
contains
  subroutine K_using_first_node(i,r,ki,kii)
    !----------------------------------------------------
    ! Evaluate the stress intensity factor by using the
    ! displacements obtained in the first pair of nodes
    ! before the crack tip being analyzed.
    !
    ! Date: 02 September 2022
    ! Last update: ---
    !----------------------------------------------------
    integer, intent(in) :: i
    real(kind=8),intent(out) :: r,ki(2),kii(2)
    integer :: j,k,tips_num,nd1,nd2
    !
    tips_num = cr_list(i)%tips_num
    do j = 1,tips_num
        if (j.eq.1) then
            k = size(cr_list(i)%nodes_num(:,1))
            nd1 = cr_list(i)%nodes_num(k,1)
            nd2 = cr_list(i)%nodes_num(k,2)
            else if (j.eq.2) then
                k = 1
                nd1 = cr_list(i)%nodes_num(k,2)
                nd2 = cr_list(i)%nodes_num(k,1)
        end if
        call ki_and_kii(i,j,nd1,nd2,r,ki,kii)
    end do
    !
    return
  end subroutine K_using_first_node
  !****************************************************
  !****************************************************
  subroutine K_using_mid_node(i,r,ki,kii)
    !----------------------------------------------------
    ! Evaluate the stress intensity factor by using the
    ! displacements obtained in the second pair of nodes
    ! before the crack tip being analyzed.
    !
    ! Date: 02 September 2022
    ! Last update: ---
    !----------------------------------------------------
    integer, intent(in) :: i
    real(kind=8),intent(out) :: r,ki(2),kii(2)
    integer :: j,k,nd1,nd2,tips_num
    !
    tips_num = cr_list(i)%tips_num
    do j = 1,tips_num
        if (j.eq.1) then
            k = size(cr_list(i)%nodes_num(:,1))-1
            nd1 = cr_list(i)%nodes_num(k,1)
            nd2 = cr_list(i)%nodes_num(k,2)
            else if (j.eq.2) then
                k = 2
                nd1 = cr_list(i)%nodes_num(k,2)
                nd2 = cr_list(i)%nodes_num(k,1)
        end if
        call ki_and_kii(i,j,nd1,nd2,r,ki,kii)
    end do
    !
    return
  end subroutine K_using_mid_node
  !****************************************************
  !****************************************************
  subroutine K_using_last_node(i,r,ki,kii)
    !----------------------------------------------------
    ! Evaluate the stress intensity factor by using the
    ! displacements obtained in the third pair of nodes
    ! before the crack tip being analyzed.
    !
    ! Date: 02 September 2022
    ! Last update: ---
    !----------------------------------------------------
    integer, intent(in) :: i
    real(kind=8),intent(out) :: r,ki(2),kii(2)
    integer :: j,k,nd1,nd2,tips_num
    !
    tips_num = cr_list(i)%tips_num
    do j = 1,tips_num
        if (j.eq.1) then
            k = size(cr_list(i)%nodes_num(:,1))-2
            nd1 = cr_list(i)%nodes_num(k,1)
            nd2 = cr_list(i)%nodes_num(k,2)
            else if (j.eq.2) then
                k = 3
                nd1 = cr_list(i)%nodes_num(k,2)
                nd2 = cr_list(i)%nodes_num(k,1)
        end if
        call ki_and_kii(i,j,nd1,nd2,r,ki,kii)
    end do
    !
    return
  end subroutine K_using_last_node
  !****************************************************
  !****************************************************
  subroutine K_using_linear_extrap_last_and_mid_node(i,ki,kii)
    !----------------------------------------------------
    ! Evaluate the stress intensity factor by doing a
    ! linear extrapolation between F-G and D-E.
    !
    ! Date: 02 September 2022
    ! Last update: ---
    !----------------------------------------------------
    integer, intent(in) :: i
    real(kind=8),intent(out) :: ki(2),kii(2)
    real(kind=8) :: r_2,ki_2(2),kii_2(2)
    real(kind=8) :: r_3,ki_3(2),kii_3(2)
    !
    !
    call K_using_mid_node(i,r_2,ki_2,kii_2)
    call K_using_last_node(i,r_3,ki_3,kii_3)
    !
    ki = (-r_2/(r_2-r_3))*(ki_2-ki_3) + ki_2
    kii = (-r_2/(r_2-r_3))*(kii_2-kii_3) + kii_2
    !
    return
  end subroutine K_using_linear_extrap_last_and_mid_node
  !****************************************************
  !****************************************************
  subroutine K_using_linear_extrap_mid_and_first_node(i,ki,kii)
    !----------------------------------------------------
    ! Evaluate the stress intensity factor by doing a
    ! linear extrapolation between D-E and B-C.
    !
    ! Date: 02 September 2022
    ! Last update: ---
    !----------------------------------------------------
    integer, intent(in) :: i
    real(kind=8),intent(out) :: ki(2),kii(2)
    real(kind=8) :: r_2,ki_2(2),kii_2(2)
    real(kind=8) :: r_1,ki_1(2),kii_1(2)
    !
    !
    call K_using_mid_node(i,r_2,ki_2,kii_2)
    call K_using_first_node(i,r_1,ki_1,kii_1)
    !
    ki = (-r_1/(r_1-r_2))*(ki_1-ki_2) + ki_1
    kii = (-r_1/(r_1-r_2))*(kii_1-kii_2) + kii_1
    !
    return
  end subroutine K_using_linear_extrap_mid_and_first_node
  !****************************************************
  !****************************************************
  subroutine K_using_extrap_with_all_three_nodes(i,ki,kii)
    !----------------------------------------------------
    ! Evaluate the stress intensity factor by doing a
    ! linear extrapolation between D-E and B-C.
    !
    ! Date: 02 September 2022
    ! Last update: ---
    !----------------------------------------------------
    integer, intent(in) :: i
    real(kind=8),intent(out) :: ki(2),kii(2)
    real(kind=8) :: r_3,ki_3(2),kii_3(2)
    real(kind=8) :: r_2,ki_2(2),kii_2(2)
    real(kind=8) :: r_1,ki_1(2),kii_1(2)
    real(kind=8),allocatable :: ni(:)
    integer :: j,k
    type(boundary_elem) :: last_elem
    !
    !
    call K_using_last_node(i,r_3,ki_3,kii_3)
    call K_using_mid_node(i,r_2,ki_2,kii_2)
    call K_using_first_node(i,r_1,ki_1,kii_1)
    !
    do j = 1,cr_list(i)%tips_num
        k = 2*j
        last_elem = be_list(cr_list(i)%ext_el_num(k))
        allocate(ni(last_elem%nodes_num))
        call intrp_func_1d(ni,1.0_8,last_elem%lamb_1,last_elem%lamb_2,last_elem%nodes_num-1)
        ki(j) = ni(1)*ki_3(j)+ni(2)*ki_2(j)+ni(3)*ki_1(j)
        kii(j) = ni(1)*kii_3(j)+ni(2)*kii_2(j)+ni(3)*kii_1(j)
        deallocate(ni)
    end do
    !
    return
  end subroutine K_using_extrap_with_all_three_nodes
  !****************************************************
  !****************************************************
  subroutine K_using_any_node(i,ki,kii)
    !----------------------------------------------------
    ! Evaluate the stress intensity factor by using the
    ! displacements obtained in the first pair of nodes
    ! before the crack tip being analyzed.
    !
    ! Date: 02 September 2022
    ! Last update: ---
    !----------------------------------------------------
    integer, intent(in) :: i
    real(kind=8),intent(out) :: ki(2),kii(2)
    integer :: j,k,tips_num,nd1,nd2
    real(kind=8) :: r
    !
    tips_num = cr_list(i)%tips_num
    do j = 1,tips_num
        if (j.eq.1) then
            k = size(cr_list(i)%nodes_num(:,1)) - nd + 1
            nd1 = cr_list(i)%nodes_num(k,1)
            nd2 = cr_list(i)%nodes_num(k,2)
            else if (j.eq.2) then
                k = nd
                nd1 = cr_list(i)%nodes_num(k,2)
                nd2 = cr_list(i)%nodes_num(k,1)
        end if
        call ki_and_kii(i,j,nd1,nd2,r,ki,kii)
    end do
    !
    return
  end subroutine K_using_any_node
  !****************************************************
  !****************************************************
  subroutine ki_and_kii(i,j,nd1,nd2,r,ki,kii)
    !----------------------------------------------------
    ! Get the ki and kii for a specific crack tip.
    !
    ! Date: 17 September 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8), intent(out) :: r,ki(2),kii(2)
    integer, intent(in) :: nd1,nd2
    integer,intent(in) :: i,j
    type(boundary_elem) :: last_elem
    real(kind=8) :: und1_1,und1_2,und2_1,und2_2,u1_1,u1_2,u2_1,u2_2,mu,kappa,normal(2)
    !
    call get_shear_modulus(mu)
    call get_kappa (kappa)
    !
    und1_1 = bn_list(nd1)%disp(1)
    und1_2 = bn_list(nd1)%disp(2)
    und2_1 = bn_list(nd2)%disp(1)
    und2_2 = bn_list(nd2)%disp(2)
    normal = get_normal_to_a_boundary_point(nd2)
    u1_1 = normal(2)*und1_1-normal(1)*und1_2
    u1_2 = normal(1)*und1_1+normal(2)*und1_2
    u2_1 = normal(2)*und2_1-normal(1)*und2_2
    u2_2 = normal(1)*und2_1+normal(2)*und2_2
    last_elem = be_list(cr_list(i)%ext_el_num(2*j))
    r = distance_only_2d(bn_list(nd1)%point,last_elem%bp_2)
    ki(j) = (u1_2-u2_2)*mu*sqrt(two_pi/r)/(kappa+1)
    kii(j) = (u1_1-u2_1)*mu*sqrt(two_pi/r)/(kappa+1)
    !
    return
  end subroutine ki_and_kii
  !****************************************************
  !****************************************************
  subroutine set_node(n)
    !----------------------------------------------------
    ! Set the number of the pair of nodes used.
    !
    !
    ! Date: 27 September 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: n
    !
    nd = n
    !
    return
  end subroutine set_node
  !****************************************************
  !****************************************************
  subroutine get_node(n)
    !----------------------------------------------------
    ! Get the number of the pair of nodes used.
    !
    !
    ! Date: 27 September 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(out) :: n
    !
    n = nd
    !
    return
  end subroutine get_node
  !****************************************************
  !****************************************************
end module fracture_2d_displ_correlation
