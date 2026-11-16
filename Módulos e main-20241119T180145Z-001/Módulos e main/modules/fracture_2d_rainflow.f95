module fracture_2d_rainflow
!----------------------------------------------------
! A module with subroutines that calculates the number
! of cycles by the rainflow cycle counting.
!
! Date: December 2022
! Last update: -
!----------------------------------------------------
use fracture_2d_model
use fracture_2d_crack_propagation
use fracture_2d_output
!
implicit none
!
real(kind=8),allocatable,private :: alpha(:) ! List with the different alphas given by the user
!
contains
  subroutine rainflow(aq,amax)
    !-----------------------------------------------------
    ! Calculate the number of cycles by the rainflow cycle
    ! counting.
    !
    ! Date: 19 December 2022
    ! Last update: ---
    !-----------------------------------------------------
    real(kind=8),intent(out) :: aq,amax
    integer :: i,j,m,cycles
    real(kind=8),allocatable :: alpha_max(:),R(:),aux_alpha_max(:),aux_R(:),aux_alpha(:)
    !
    call rearrange_alpha()
    !
    allocate(aux_alpha_max(10000))
    allocate(aux_R(10000))
    j = size(alpha(:))
    i = 1
    cycles = 0
    do while (j .ne. 1)
        m = size(alpha(:))
        if (abs(alpha(i)-alpha(i+1)) .le. abs(alpha(i+1)-alpha(i+2))) then
            cycles = cycles + 1
            aux_alpha_max(cycles) = max(alpha(i),alpha(i+1))
            aux_R(cycles) = min(alpha(i),alpha(i+1))/max(alpha(i),alpha(i+1))
            j = j - 2
            allocate(aux_alpha(size(alpha(:))))
            aux_alpha = alpha
            deallocate(alpha)
            allocate(alpha(j))
            if (i.eq.1) then
                alpha(1:j) = aux_alpha(3:m)
                else
                    alpha(1:i-1) = aux_alpha(1:i-1)
                    alpha(i:j) = aux_alpha(i+2:m)
            end if
            deallocate(aux_alpha)
            i = 1
            else
                i = i + 1
        end if
    end do
    !
    allocate(alpha_max(cycles))
    alpha_max = aux_alpha_max(1:cycles)
    deallocate(aux_alpha_max)
    !
    allocate(R(cycles))
    R = aux_R(1:cycles)
    deallocate(aux_R)
    !
    call cycles_output(cycles,alpha_max,R)
    !
    !call deallocate_alpha()
    deallocate(alpha) !Problema está aqui
    !
    call alpha_q_alpha_max(cycles,alpha_max,R,aq,amax)
    !
    deallocate(alpha_max)
    deallocate(R)
    !
  end subroutine rainflow
  !****************************************************
  !****************************************************
  subroutine rearrange_alpha()
    !----------------------------------------------------
    ! Rearrange alpha, so the maximum value is the first
    ! of the list.
    !
    ! Date: 19 December 2022
    ! Last update: ---
    !----------------------------------------------------
    integer :: i,j,k,m
    real(kind=8),allocatable :: aux(:)
    real(kind=8) :: x
    !
    x = alpha(1)
    j = 1
    do i = 1,size(alpha(:))
        if (abs(x) .ge. abs(alpha(i))) then
            x = x
            else if (abs(x) .lt. abs(alpha(i))) then
                x = alpha(i)
                j = i
        end if
    end do
    !
    m = size(alpha(:))
    allocate(aux(m+1))
    k = m - j + 1
    aux(1:k) = alpha(j:m)
    aux(k+1:m) = alpha(1:j-1)
    aux(m+1) = aux(1)
    deallocate(alpha)
    allocate(alpha(m+1))
    alpha = aux
    deallocate(aux)
    !
  end subroutine rearrange_alpha
  !****************************************************
  !****************************************************
  subroutine alpha_q_alpha_max(cycles,alpha_max,R,aq,amax)
    !----------------------------------------------------
    ! Calculate the average alpha and the maximum alpha.
    !
    ! Date: 19 December 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: cycles
    real(kind=8),intent(in) :: alpha_max(:),R(:)
    real(kind=8),intent(out) :: aq,amax
    real(kind=8) :: gama,m,alpha_,x
    integer :: i
    !
    call get_parameter_gamma(gama)
    call get_parameter_m(m)
    !
    alpha_ = 0.0
    do i = 1,cycles
        alpha_ = alpha_ + (alpha_max(i)*(1-R(i))**gama)**m
    end do
    aq = (alpha_/cycles)**(1/m)
    !
    x = alpha_max(1)
    do i = 1,cycles
        if (x .gt. alpha_max(i)) then
            x = x
            else if (x .le. alpha_max(i)) then
                x = alpha_max(i)
        end if
    end do
    amax = x/aq
    !
  end subroutine alpha_q_alpha_max
  !****************************************************
  !****************************************************
  subroutine set_alpha(i,a)
    !----------------------------------------------------
    ! Set the values on the list alpha.
    !
    !
    ! Date: 19 December 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: i
    real(kind=8),intent(in) :: a
    !
    alpha(i) = a
    !
    return
  end subroutine set_alpha
  !****************************************************
  !****************************************************
  subroutine allocate_alpha(n)
    !----------------------------------------------------
    ! Allocate the list alpha.
    !
    !
    ! Date: 19 December 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: n
    !
    allocate(alpha(n))
    !
    return
  end subroutine allocate_alpha
  !****************************************************
  !****************************************************
end module fracture_2d_rainflow
