module fracture_2d_stress_intens_factor
!----------------------------------------------------
! A module with subroutines that calculates the stress
! intensity factors (K).
!
! Date: September 2022
! Last update: -
!----------------------------------------------------
use fracture_2d_displ_correlation
use fracture_2d_j_integral
!
implicit none
!
integer,private :: K_key ! 1 = using displ. correl., 2 = using J-integral
integer,private :: method_key ! 1 to 6 for displ. correl. and 1 to 3 for j-integral and displ. fitting
!
contains
  subroutine stress_intensity_factor_evaluation()
    !----------------------------------------------------
    ! Decides if the stress intensity factor will be evaluated
    ! by the displacement correlation method or by the J-integral.
    !
    ! Date: 02 September 2022
    ! Last update: ---
    !----------------------------------------------------
    !
    !
    select case (K_key)
    case (1)
        call K_using_displ_correlation()
    case (2)
        call K_using_J_integral()
    case default
        call error_message_stop("Incorrect key for stress intensity factor evaluation.")
    end select
    !
    return
  end subroutine stress_intensity_factor_evaluation
  !****************************************************
  !****************************************************
  subroutine K_using_displ_correlation()
    !----------------------------------------------------
    ! Verify which mode of displacement correlation will be
    ! used to calculate the stress intensity factor.
    !
    ! 1- Using displacements B-C (first nodes before the crack tip)
    ! 2- Using displacements D-E (second nodes before the crack tip)
    ! 3- Using displacements F-G (third nodes before the crack tip)
    ! 4- Linear extrapolation F-G and D-E
    ! 5- Linear extrapolation D-E and B-C
    ! 6- Extrapolation using the interpolation functions, using B-C, D-E and F-G
    !
    ! Date: 02 September 2022
    ! Last update: ---
    !---------------------------------------------------
    !
    integer :: i
    real(kind=8) :: r,ki(2),kii(2)
    !
    do i = 1,cr_num
        select case(method_key)
            case(1)
                call K_using_first_node(i,r,ki,kii)
                cr_list(i)%ki = ki
                cr_list(i)%kii = kii
            case(2)
                call K_using_mid_node(i,r,ki,kii)
                cr_list(i)%ki = ki
                cr_list(i)%kii = kii
            case(3)
                call K_using_last_node(i,r,ki,kii)
                cr_list(i)%ki = ki
                cr_list(i)%kii = kii
            case(4)
                call K_using_linear_extrap_last_and_mid_node(i,ki,kii)
                cr_list(i)%ki = ki
                cr_list(i)%kii = kii
            case(5)
                call K_using_linear_extrap_mid_and_first_node(i,ki,kii)
                cr_list(i)%ki = ki
                cr_list(i)%kii = kii
           case(6)
                call K_using_extrap_with_all_three_nodes(i,ki,kii)
                cr_list(i)%ki = ki
                cr_list(i)%kii = kii
            case(7)
                call K_using_any_node(i,ki,kii)
                cr_list(i)%ki = ki
                cr_list(i)%kii = kii
            case default
                call error_message_stop("Incorrect key for evaluation of the stress intensity factor &
                &using the displacement correlation method.")
        end select
    end do
    !
    return
  end subroutine K_using_displ_correlation
  !****************************************************
  !****************************************************
  subroutine K_using_J_integral()
    !----------------------------------------------------
    ! Verify which mode of J-integration will be
    ! used to calculate the stress intensity factor.
    !
    ! 1- Using Simpson's rule
    ! 2- Using Gauss integration points
    ! 3- Creating a circular boundary and integrating it
    !
    ! Date: 06 September 2022
    ! Last update: ---
    !---------------------------------------------------
    !
    select case(method_key)
        case(1)
            call K_using_Simpson_rule()
        case(2)
            call K_using_Gauss_integration_points()
        case(3)
            call K_using_Gauss_integration_points_on_boundary()
            case default
                call error_message_stop("Incorrect key for evaluation of the stress intensity factor &
                &using the J-integral method.")
        end select
    !
    return
  end subroutine K_using_J_integral
  !****************************************************
  !****************************************************
  subroutine set_K_key(k)
    !----------------------------------------------------
    ! Set the stress intensity factor evaluation key.
    !
    !
    ! Date: 02 September 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: k
    !
    K_key = k
    !
    return
  end subroutine set_K_key
  !****************************************************
  !****************************************************
  subroutine get_K_key(k)
    !----------------------------------------------------
    ! Get the stress intensity factor evaluation key.
    !
    ! Date: 02 September 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(out) :: k
    !
    k = K_key
    !
    return
  end subroutine get_K_key
  !****************************************************
  !****************************************************
  subroutine set_method_key(k)
    !----------------------------------------------------
    ! Set the chosen method key.
    !
    ! For displacement correlation - 1 to 6
    ! For J-integral - 1 to 3
    ! For displacement fitting - 1 to 3
    !
    !
    ! Date: 25 September 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: k
    !
    method_key = k
    !
    return
  end subroutine set_method_key
  !****************************************************
  !****************************************************
  subroutine get_method_key(k)
    !----------------------------------------------------
    ! Get the chosen method key.
    !
    ! For displacement correlation - 1 to 6
    ! For J-integral - 1 to 3
    ! For displacement fitting - 1 to 3
    !
    !
    ! Date: 25 September 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(out) :: k
    !
    k = method_key
    !
    return
  end subroutine get_method_key
  !****************************************************
  !****************************************************
end module fracture_2d_stress_intens_factor
