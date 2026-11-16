module vector_algebra
!----------------------------------------------------
! A module for vector algebra.
!
! Date: May 2018
! Last update: ---
!----------------------------------------------------
use error_messages
implicit none
!
contains
  function vector_len(v)
    !----------------------------------------------------
    ! Computes the norm of a vector of any dimension.
    !
    ! Date: 30 May 2018
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: v(:)
    real(kind=8) :: vector_len
    vector_len=sqrt(dot_product(v,v))
    return
  end function vector_len
  !****************************************************
  !****************************************************
  subroutine vector_norm(v)
    !----------------------------------------------------
    ! Normalizes a vector of any dimension.
    !
    ! Date: 30 May 2018
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(inout) :: v(:)
    real(kind=8) :: vlen
    real(kind=8) :: tol = 1.0e-12
    !
    vlen = vector_len(v)
    if (abs(vlen).lt.tol) return
    v = v/vlen
    return
  end subroutine vector_norm
  !****************************************************
  !****************************************************
  function vector_cross_prod(v,u)
    !----------------------------------------------------
    ! Calculates v x u (the cross product).
    !
    ! Date: 30 May 2018
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: v(3),u(3)
    real(kind=8) :: vector_cross_prod(3)
    vector_cross_prod(1) = v(2)*u(3)-v(3)*u(2)
    vector_cross_prod(2) = v(3)*u(1)-v(1)*u(3)
    vector_cross_prod(3) = v(1)*u(2)-v(2)*u(1)
    return
  end function
  !****************************************************
  !****************************************************
end module vector_algebra
