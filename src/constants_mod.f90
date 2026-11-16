!-----------------------------------------------------------------------
! MODULE: constants_mod
!> @brief Constantes matemáticas, tipos geométricos 2D, álgebra vetorial/matricial
!>        e interface para o solver linear LAPACK (DGESV).
!>
!> Código reestruturado e simplificado para fins didáticos, utilizando
!> exclusivamente elementos de contorno quadráticos descontínuos.
!-----------------------------------------------------------------------
module constants_mod
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none

  ! Constantes matemáticas de dupla precisão
  real(dp), parameter :: pi          = 3.14159265358979323846_dp
  real(dp), parameter :: two_pi      = 6.28318530717958647692_dp
  real(dp), parameter :: four_pi     = 12.56637061435917295385_dp
  real(dp), parameter :: eight_pi    = 25.13274122871834590770_dp
  real(dp), parameter :: sixteen_pi  = 50.26548245743669181540_dp

  ! Tipo para ponto cartesiano 2D
  type :: point_2d
    real(dp) :: x = 0.0_dp
    real(dp) :: y = 0.0_dp
  end type point_2d

  ! Interface explícita para o solver linear denso LAPACK (DGESV)
  interface
    subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
      import :: dp
      integer, intent(in)     :: n, nrhs, lda, ldb
      integer, intent(out)    :: ipiv(*), info
      real(dp), intent(inout) :: a(lda, *), b(ldb, *)
    end subroutine dgesv
  end interface

contains

  !---------------------------------------------------------------------
  !> @brief Exibe mensagem de erro e encerra a execução com código não nulo.
  !---------------------------------------------------------------------
  subroutine error_message_stop(e_msg)
    character(len=*), intent(in) :: e_msg
    write(*, '(/,A,A,/)') ' [ERRO] ', trim(e_msg)
    stop 1
  end subroutine error_message_stop

  !---------------------------------------------------------------------
  !> @brief Exibe mensagem de erro com valor inteiro e encerra a execução.
  !---------------------------------------------------------------------
  subroutine error_message_int_stop(e_msg, val)
    character(len=*), intent(in) :: e_msg
    integer, intent(in)          :: val
    write(*, '(/,A,A,I0,/)') ' [ERRO] ', trim(e_msg), val
    stop 1
  end subroutine error_message_int_stop

  !---------------------------------------------------------------------
  !> @brief Exibe mensagem de erro com dois valores inteiros e encerra.
  !---------------------------------------------------------------------
  subroutine error_message_2int_stop(msg1, v1, msg2, v2)
    character(len=*), intent(in) :: msg1, msg2
    integer, intent(in)          :: v1, v2
    write(*, '(/,A,A,I0,1X,A,I0,/)') ' [ERRO] ', trim(msg1), v1, trim(msg2), v2
    stop 1
  end subroutine error_message_2int_stop

  !---------------------------------------------------------------------
  !> @brief Exibe um aviso na tela sem interromper a execução.
  !---------------------------------------------------------------------
  subroutine warning_message(w_msg)
    character(len=*), intent(in) :: w_msg
    write(*, '(/,A,A,/)') ' [AVISO] ', trim(w_msg)
  end subroutine warning_message

  !---------------------------------------------------------------------
  !> @brief Calcula o vetor distância r = [|r|, rx, ry] de p1 (fonte) até p2 (campo).
  !---------------------------------------------------------------------
  pure function distance_2d(p1, p2) result(d)
    type(point_2d), intent(in) :: p1, p2
    real(dp) :: d(3)
    d(2) = p2%x - p1%x
    d(3) = p2%y - p1%y
    d(1) = sqrt(d(2)*d(2) + d(3)*d(3))
  end function distance_2d

  !---------------------------------------------------------------------
  !> @brief Calcula apenas a distância escalar euclidiana entre p1 e p2.
  !---------------------------------------------------------------------
  pure function distance_only_2d(p1, p2) result(dist)
    type(point_2d), intent(in) :: p1, p2
    real(dp) :: dist
    dist = sqrt((p2%x - p1%x)**2 + (p2%y - p1%y)**2)
  end function distance_only_2d

  !---------------------------------------------------------------------
  !> @brief Calcula o ponto médio entre dois pontos 2D.
  !---------------------------------------------------------------------
  pure function midpoint_2d(p1, p2) result(mid)
    type(point_2d), intent(in) :: p1, p2
    type(point_2d) :: mid
    mid%x = 0.5_dp * (p1%x + p2%x)
    mid%y = 0.5_dp * (p1%y + p2%y)
  end function midpoint_2d

  !---------------------------------------------------------------------
  !> @brief Divide um segmento de reta (p1 -> p2) em n_el elementos quadráticos
  !>        (retornando 2*n_el + 1 pontos geométricos igualmente espaçados).
  !---------------------------------------------------------------------
  subroutine divide_straight_line_quad(p1, p2, n_el, pts)
    type(point_2d), intent(in)               :: p1, p2
    integer, intent(in)                      :: n_el
    type(point_2d), allocatable, intent(out) :: pts(:)
    integer :: i, n_pts
    real(dp) :: dx, dy

    n_pts = 2 * n_el + 1
    allocate(pts(n_pts))
    dx = (p2%x - p1%x) / real(2 * n_el, dp)
    dy = (p2%y - p1%y) / real(2 * n_el, dp)
    do i = 1, n_pts
      pts(i)%x = p1%x + real(i - 1, dp) * dx
      pts(i)%y = p1%y + real(i - 1, dp) * dy
    end do
  end subroutine divide_straight_line_quad

  !---------------------------------------------------------------------
  !> @brief Divide um arco circular centrado em p_center, iniciando em p_start,
  !>        com ângulo total ang_rad (radianos), em n_el elementos quadráticos
  !>        (retornando 2*n_el + 1 pontos sobre o arco circular).
  !---------------------------------------------------------------------
  subroutine divide_circular_arc_quad(p_center, p_start, ang_rad, n_el, ccw, pts)
    type(point_2d), intent(in)               :: p_center, p_start
    real(dp), intent(in)                     :: ang_rad
    integer, intent(in)                      :: n_el
    logical, intent(in)                      :: ccw
    type(point_2d), allocatable, intent(out) :: pts(:)
    integer :: i, n_pts
    real(dp) :: rad(3), theta0, dtheta, theta, sgn

    n_pts = 2 * n_el + 1
    allocate(pts(n_pts))
    rad = distance_2d(p_center, p_start)
    theta0 = atan2(rad(3), rad(2))
    sgn = merge(1.0_dp, -1.0_dp, ccw)
    dtheta = sgn * ang_rad / real(2 * n_el, dp)

    pts(1) = p_start
    do i = 2, n_pts
      theta = theta0 + real(i - 1, dp) * dtheta
      pts(i)%x = p_center%x + rad(1) * cos(theta)
      pts(i)%y = p_center%y + rad(1) * sin(theta)
    end do
  end subroutine divide_circular_arc_quad

  !---------------------------------------------------------------------
  !> @brief Norma euclidiana de um vetor.
  !---------------------------------------------------------------------
  pure function vector_len(v) result(vlen)
    real(dp), intent(in) :: v(:)
    real(dp) :: vlen
    vlen = sqrt(dot_product(v, v))
  end function vector_len

  !---------------------------------------------------------------------
  !> @brief Normaliza um vetor in-place.
  !---------------------------------------------------------------------
  pure subroutine vector_norm(v)
    real(dp), intent(inout) :: v(:)
    real(dp) :: vlen
    vlen = vector_len(v)
    if (vlen > 1.0e-12_dp) v = v / vlen
  end subroutine vector_norm

  !---------------------------------------------------------------------
  !> @brief Inverte uma matriz quadrada pequena in-place (usada no método
  !>        de expansão radial de Gao, tipicamente 3x3 a 9x9).
  !---------------------------------------------------------------------
  subroutine mat_inverse(mat)
    real(dp), intent(inout) :: mat(:,:)
    integer, allocatable :: alo(:), amo(:)
    real(dp) :: biga, hold
    integer :: i, j, k, n

    if (size(mat, 1) /= size(mat, 2)) then
      call error_message_stop("Tentativa de inverter matriz nao quadrada.")
    end if

    n = size(mat, 1)
    allocate(alo(n), amo(n))

    do k = 1, n
      alo(k) = k
      amo(k) = k
      biga = mat(k, k)
      do j = k, n
        do i = k, n
          if (abs(biga) < abs(mat(i, j))) then
            biga = mat(i, j)
            alo(k) = i
            amo(k) = j
          end if
        end do
      end do

      j = alo(k)
      if (j > k) then
        do i = 1, n
          hold = -mat(k, i)
          mat(k, i) = mat(j, i)
          mat(j, i) = hold
        end do
      end if

      i = amo(k)
      if (i > k) then
        do j = 1, n
          hold = -mat(j, k)
          mat(j, k) = mat(j, i)
          mat(j, i) = hold
        end do
      end if

      if (abs(biga) < 1.0e-28_dp) then
        call error_message_stop("Matriz singular durante inversao local.")
      end if

      do i = 1, n
        if (i /= k) mat(i, k) = mat(i, k) / (-biga)
      end do

      do i = 1, n
        hold = mat(i, k)
        do j = 1, n
          if (i /= k .and. j /= k) mat(i, j) = hold * mat(k, j) + mat(i, j)
        end do
      end do

      do j = 1, n
        if (j /= k) mat(k, j) = mat(k, j) / biga
      end do

      mat(k, k) = 1.0_dp / biga
    end do

    k = n - 1
    do while (k > 0)
      i = alo(k)
      if (i > k) then
        do j = 1, n
          hold = mat(j, k)
          mat(j, k) = -mat(j, i)
          mat(j, i) = hold
        end do
      end if
      j = amo(k)
      if (j > k) then
        do i = 1, n
          hold = mat(k, i)
          mat(k, i) = -mat(j, i)
          mat(j, i) = hold
        end do
      end if
      k = k - 1
    end do

    deallocate(alo, amo)
  end subroutine mat_inverse

  !---------------------------------------------------------------------
  !> @brief Resolve o sistema linear denso [A]{x} = {f} in-place usando
  !>        a rotina DGESV da biblioteca open-source LAPACK.
  !---------------------------------------------------------------------
  subroutine solve_linear_system(a, f)
    real(dp), intent(inout) :: a(:,:), f(:)
    real(dp), allocatable   :: rhs(:,:)
    integer, allocatable    :: ipiv(:)
    integer :: n, info

    n = size(a, 1)
    if (size(a, 2) /= n .or. size(f) /= n) then
      call error_message_stop("Dimensoes incompativeis no solver linear.")
    end if

    allocate(ipiv(n), rhs(n, 1))
    rhs(:, 1) = f

    call dgesv(n, 1, a, n, ipiv, rhs, n, info)
    if (info > 0) then
      call error_message_stop("Matriz global singular em DGESV.")
    else if (info < 0) then
      call error_message_int_stop("Argumento invalido em DGESV: ", -info)
    end if

    f = rhs(:, 1)
    deallocate(rhs, ipiv)
  end subroutine solve_linear_system

end module constants_mod
