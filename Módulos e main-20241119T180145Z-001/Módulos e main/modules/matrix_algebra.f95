module matrix_algebra
!----------------------------------------------------
! A module for matrix algebra.
!
! Date: September 2019
! Last update: July 2020
!----------------------------------------------------
use error_messages
implicit none
!
interface
  subroutine dgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
    integer,intent(in) :: n,nrhs,lda,ldb
    integer,intent(out) :: ipiv(*),info
    real(kind=8),intent(inout) :: a(lda,*),b(ldb,*)
  end subroutine dgesv
end interface
!
contains
  subroutine mat_inverse(mat)
    !----------------------------------------------------
    ! Substitutes the given matrix, mat, by its inverse.
    !
    ! Date: 10 September 2019
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(inout) :: mat(:,:)
    integer,allocatable :: alo(:),amo(:)
    real(kind=8) :: biga,hold
    integer :: i,j,k,n
    !
    ! Check if mat is a square matrix:
    !
    if (size(mat(1,:)).ne.size(mat(:,1))) then
      call error_message_stop("Trying to invert a non-square matrix.")
    end if
    !
    n = size(mat(1,:))
    allocate(alo(n))
    allocate(amo(n))
    !
    do k = 1,n
      !
      ! Search for the major element:
      !
      alo(k) = k
      amo(k) = k
      biga = mat(k,k)
      do j = k,n
        do i = k,n
          if (abs(biga).ge.abs(mat(i,j))) cycle
          biga = mat(i,j)
          alo(k) = i
          amo(k) = j
        end do
      end do
      !
      ! Lines change:
      !
      j = alo(k)
      if (j.gt.k) then
        do i = 1,n
          hold = -mat(k,i)
          mat(k,i) = mat(j,i)
          mat(j,i) = hold
        end do
      end if
      !
      ! Columns change:
      !
      i = amo(k)
      if (i.gt.k) then
        do j = 1,n
          hold = -mat(j,k)
          mat(j,k) = mat(j,i)
          mat(j,i) = hold
        end do
      end if
      !
      ! Check if mat is a singular matrix:
      !
      if (abs(biga).lt.1e-28_8) then
        call error_message_stop("Trying to invert a singular matrix.")
      end if
      !
      ! Divide column by negative pivot:
      !
      do i = 1,n
        if (i.eq.k) cycle
        mat(i,k) = mat(i,k)/(-biga)
      end do
      !
      ! Matrix reduction:
      !
      do i = 1,n
        hold = mat(i,k)
        do j = 1,n
          if (i.eq.k .or. j.eq.k) cycle
          mat(i,j) = hold*mat(k,j) + mat(i,j)
        end do
      end do
      !
      ! Divide line by the pivot:
      !
      do j = 1,n
        if (j.eq.k) cycle
        mat(k,j) = mat(k,j)/biga
      end do
      !
      ! substitute the pivot by its inverse:
      !
      mat(k,k)=1._8/biga
      !
    end do
    !
    ! final change of lines and columns:
    !
    k = n-1
    do while (k.gt.0)
      i = alo(k)
      if (i.gt.k) then
        do j = 1,n
          hold = mat(j,k)
          mat(j,k) = -mat(j,i)
          mat(j,i) = hold
        end do
      end if
      j = amo(k)
      if (j.gt.k) then
        do i = 1,n
          hold = mat(k,i)
          mat(k,i) = -mat(j,i)
          mat(j,i) = hold
        end do
      end if
      k = k-1
    end do
    !
    deallocate(alo)
    deallocate(amo)
    !
    return
  end subroutine mat_inverse
  !****************************************************
  !****************************************************
  subroutine jacobi_eigenvalues(mat,evt,tol)
    !----------------------------------------------------
    ! Calculates the eigenvalues and eigenvectors of a
    ! symmetric matrix. The given matrix is substituted by
    ! the diagonalized form.
    !
    ! mat = the given matrix.
    ! evt = eigenvectors (each column).
    ! tol = tolerance [sum of (off-diagonal elements)^2]
    !
    ! Date: 04 July 2020
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(inout) :: mat(:,:)
    real(kind=8),intent(out) :: evt(:,:)
    real(kind=8) :: tol,b2,bar
    real(kind=8) :: beta,coeff,c,s,cs,sc
    real(kind=8) :: aux_value
    real(kind=8),allocatable :: aux_vector(:)
    integer :: i,j,k,n
    !
    if (size(mat(1,:)).ne.size(mat(:,1))) then
      call error_message_stop("Non-square matrix for eigenvalues evaluation.")
    end if
    !
    n = size(mat(1,:))
    !
    ! Check in the given matrix is symmetric:
    !
    do i = 1,n
      do j = 1,n
        !
        if (abs(mat(i,j)).lt.tol) then
          if (abs(mat(j,i)).lt.tol) then
            cycle
          else
            call error_message_stop("Non-symmetric matrix for Jacobi &
                                    &eigenvalues evaluation")
          end if
        end if
        !
        if (abs((mat(i,j)-mat(j,i))/mat(i,j)).gt.tol) then
          call error_message_stop("Non-symmetric matrix for Jacobi &
                                    &eigenvalues evaluation")
        end if
      end do
    end do
    !
    if (size(evt(1,:)).lt.n .or. size(evt(:,1)).lt.n) then
      call error_message_stop("No memory allocated for eigenvectors evaluation.")
    end if
    !
    ! Initialize evt as an identity matrix:
    !
    evt = 0.0_8
    do i = 1,n
      evt(i,i) = 1.0
    end do
    !
    ! Check if the given  matrix is already diagonal:
    !
    b2 = 0.0
    do i = 1,n
      do j = 1,n
        if (i.ne.j) b2 = b2 + mat(i,j)**2
      end do
    end do
    if (b2.le.tol) go to 1
    !
    ! Evaluate eigenvalues and eigenvectors:
    !
    bar = 0.5*b2/float(n*n)
    !
    do while (b2.gt.tol)
      do i = 1,n-1
        do j = i+1,n
          if (mat(j,i)**2 .le. bar) cycle  ! do not touch small elements
          b2 = b2 - 2.0*mat(j,i)**2
          bar = 0.5*b2/float(n*n)
          !
          beta = (mat(j,j)-mat(i,i))/(2.0*mat(j,i))
          coeff = 0.5*beta/sqrt(1.0+beta**2)
          s = sqrt(max(0.5+coeff,0.0))
          c = sqrt(max(0.5-coeff,0.0))
          ! recalculate rows i and j:
          do k = 1,n
            cs = c*mat(i,k)+s*mat(j,k)
            sc = -s*mat(i,k)+c*mat(j,k)
            mat(i,k) = cs
            mat(j,k) = sc
          end do
          ! new matrix mat_{k+1} from mat_{k}, and eigenvectors:
          do k = 1,n
            cs = c*mat(k,i)+s*mat(k,j)
            sc = -s*mat(k,i)+c*mat(k,j)
            mat(k,i) = cs
            mat(k,j) = sc
            cs =  c*evt(k,i)+s*evt(k,j)
            sc = -s*evt(k,i)+c*evt(k,j)
            evt(k,i) = cs
            evt(k,j) = sc
          end do
        end do
      end do
    end do
    !
    ! Re-ordering eigenvalues and eigenvectors:
    !
    1 allocate(aux_vector(n))
    !
    do i = 1,(n-1)
      k = i
      do j = (k+1),n
        if (mat(j,j).gt.mat(k,k)) k = j
      end do
      if (i.ne.k) then
        !
        aux_value = mat(k,k)
        mat(k,k) = mat(i,i)
        mat(i,i) = aux_value
        !
        aux_vector = evt(:,k)
        evt(:,k) = evt(:,i)
        evt(:,i) = aux_vector
        !
      end if
    end do
    !
    deallocate(aux_vector)
    !
    ! Taking the major eigenvector component as positive:
    !
    do j = 1,n
      aux_value = 0.0
      do i = 1,n
        if (abs(evt(i,j)).gt.abs(aux_value)) aux_value = evt(i,j)
      end do
      if (aux_value.lt.0.0) evt(:,j) = -evt(:,j)
    end do
    !
    return
  end subroutine jacobi_eigenvalues
  !****************************************************
  !****************************************************
  subroutine solve_linear_system(a,f)
    !----------------------------------------------------
    ! Solves the square linear system a*x = f in-place
    ! using open-source LAPACK (dgesv).
    !----------------------------------------------------
    real(kind=8),intent(inout) :: a(:,:),f(:)
    real(kind=8),allocatable :: rhs(:,:)
    integer,allocatable :: ipiv(:)
    integer :: n,info
    !
    n = size(a,1)
    if (size(a,2).ne.n .or. size(f).ne.n) then
      call error_message_stop("Incompatible matrix/vector dimensions in linear solver.")
    end if
    !
    allocate(ipiv(n))
    allocate(rhs(n,1))
    rhs(:,1) = f
    !
    call dgesv(n,1,a,n,ipiv,rhs,n,info)
    if (info.gt.0) then
      call error_message_stop("Singular matrix.")
    else if (info.lt.0) then
      call error_message_int_stop("Illegal value at argument ",-info)
    end if
    !
    f = rhs(:,1)
    deallocate(rhs)
    deallocate(ipiv)
    !
    return
  end subroutine solve_linear_system
  !****************************************************
  !****************************************************
end module matrix_algebra
