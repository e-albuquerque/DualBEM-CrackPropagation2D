module geometry
!----------------------------------------------------
! A module for geometric definitions.
!
! Date: July 2018
! Last update: December 2021
!----------------------------------------------------
implicit none
!
type :: point_2d
  real(kind=8) :: x,y
end type point_2d
!
type :: point_3d
  real(kind=8) :: x,y,z
end type point_3d
!
contains
  function distance_2d(p1,p2)
    !----------------------------------------------------
    ! Computes the distance between points p1 (initial)
    ! and p2 (final).
    ! p1 and p2 are of point_2d type.
    !
    ! Returns a vector with:
    ! distance_2d(1) = r (the distance)
    ! distance_2d(2) = rx (the x projection)
    ! distance_2d(3) = ry (the y projection)
    !
    ! Date: 13 July 2018
    ! Last update: ---
    !----------------------------------------------------
    type(point_2d),intent(in) :: p1,p2
    real(kind=8) :: distance_2d(3)
    distance_2d(2)=(p2%x-p1%x)
    distance_2d(3)=(p2%y-p1%y)
    distance_2d(1)=sqrt(distance_2d(2)**2+distance_2d(3)**2)
    return
  end function distance_2d
  !****************************************************
  !****************************************************
  function distance_only_2d(p1,p2)
    !----------------------------------------------------
    ! Compute the distance between points p1 and p2.
    ! p1 and p2 are of point_2d type.
    !
    ! Returns the distance as a real(8)
    !
    ! Date: 30 July 2018
    ! Last update: ---
    !----------------------------------------------------
    type(point_2d),intent(in) :: p1,p2
    real(kind=8) :: distance_only_2d
    distance_only_2d=sqrt((p2%x-p1%x)**2+(p2%y-p1%y)**2)
    return
  end function distance_only_2d
  !****************************************************
  !****************************************************
  function midpoint_2d(p1,p2)
    !----------------------------------------------------
    ! Compute the midpoint in the straight segment
    ! from point p1 to point p2 in 2D.
    !
    ! Date: 27 July 2018
    ! Last update: ---
    !----------------------------------------------------
    type(point_2d),intent(in) :: p1,p2
    type(point_2d) :: midpoint_2d
    midpoint_2d%x = 0.5*(p1%x+p2%x)
    midpoint_2d%y = 0.5*(p1%y+p2%y)
    return
  end function midpoint_2d
  !****************************************************
  !****************************************************
  function circular_midpoint_2d(p1,p2,theta,counter_clockwise)
    !----------------------------------------------------
    ! Compute the midpoint in the arc-length, with p1
    ! as central point, beginning at point p2 and with
    ! an angle of theta.
    !
    ! counter_clockwise = .true.
    ! theta goes counter-clockwise
    !
    ! counter_clockwise = .false.
    ! theta goes clockwise
    !
    ! Date: 08 February 2018
    ! Last update: 01 July 2021
    !----------------------------------------------------
    type(point_2d),intent(in) :: p1,p2
    real(kind=8),intent(in) :: theta
    logical,intent(in) :: counter_clockwise
    real(kind=8) :: rad(3),alpha
    type(point_2d) :: circular_midpoint_2d
    !
    rad = distance_2d(p1,p2)
    alpha = atan2(rad(3),rad(2)) !The angular coordinate of p2
    if (counter_clockwise) then
      alpha = alpha + 0.5_8*theta
    else
      alpha = alpha - 0.5_8*theta
    end if
    circular_midpoint_2d%x = p1%x + rad(1)*cos(alpha)
    circular_midpoint_2d%y = p1%y + rad(1)*sin(alpha)
    !
    return
  end function circular_midpoint_2d
  !****************************************************
  !****************************************************
  function point_to_point_dist_less_than_2d(p1,p2,tol)
    !----------------------------------------------------
    ! Verify if 2D points, p1 and p2, are closer (to
    ! each other) than a given real parameter, tol.
    !
    ! Returns a logical variable.
    !
    ! Date: 30 July 2018
    ! Last update: ---
    !----------------------------------------------------
    type(point_2d),intent(in) :: p1,p2
    real(kind=8),intent(in) :: tol
    real(kind=8) :: dist
    logical :: point_to_point_dist_less_than_2d
    dist = distance_only_2d(p1,p2)
    if (dist.lt.abs(tol)) then
      point_to_point_dist_less_than_2d = .true.
    else
      point_to_point_dist_less_than_2d = .false.
    end if
    return
  end function point_to_point_dist_less_than_2d
  !****************************************************
  !****************************************************
  function point_to_point_dist_greater_than_2d(p1,p2,tol)
    !----------------------------------------------------
    ! Verify if the distance between two points, p1 and
    ! p2, is greater than a given real value, tol.
    !
    ! Returns a logical variable.
    !
    ! Date: 30 July 2018
    ! Last update: ---
    !----------------------------------------------------
    type(point_2d),intent(in) :: p1,p2
    real(kind=8),intent(in) :: tol
    real(kind=8) :: dist
    logical :: point_to_point_dist_greater_than_2d
    dist = distance_only_2d(p1,p2)
    if (dist.gt.abs(tol)) then
      point_to_point_dist_greater_than_2d = .true.
    else
      point_to_point_dist_greater_than_2d = .false.
    end if
    return
  end function point_to_point_dist_greater_than_2d
  !****************************************************
  !****************************************************
  subroutine divide_straight_line_2d(p1,p2,n,pts)
    !----------------------------------------------------
    ! Divide a straight line (between points p1 and p2)
    ! in "n" segments and return array "pts" with the
    ! segments boundary points.
    !
    ! Date: 26 July 2018
    ! Last update: ---
    !----------------------------------------------------
    type(point_2d),intent(in) :: p1,p2
    integer,intent(in) :: n
    type(point_2d),allocatable,intent(out) :: pts(:)
    real(kind=8) :: dx,dy
    integer :: i
    !
    if (allocated(pts)) deallocate(pts)
    allocate(pts(n+1))
    pts(1) = p1
    pts(n+1) = p2
    dx = (p2%x-p1%x)/n
    dy = (p2%y-p1%y)/n
    do i = 2,n
      pts(i)%x = pts(i-1)%x + dx
      pts(i)%y = pts(i-1)%y + dy
    end do
    return
  end subroutine divide_straight_line_2d
  !****************************************************
  !****************************************************
  subroutine divide_straight_line_2d_scaling(p1,p2,n,scl,pts)
    !----------------------------------------------------
    ! Divide a straight line (between points p1 and p2)
    ! in "n" segments, scaling its size, and return
    ! array "pts" with the segments boundary points.
    !
    ! The last segment has a size equal the first
    ! multiplied by parameter "scl".
    !
    ! Date: 17 December 2021
    ! Last update: ---
    !----------------------------------------------------
    type(point_2d),intent(in) :: p1,p2
    real(kind=8),intent(in) :: scl
    integer,intent(in) :: n
    type(point_2d),allocatable,intent(out) :: pts(:)
    real(kind=8) :: dx,dy,dx0,dy0
    integer :: i
    !
    if (allocated(pts)) deallocate(pts)
    allocate(pts(n+1))
    !
    pts(1) = p1
    pts(n+1) = p2
    dx0 = 2.0_8*(p2%x-p1%x)/(n*(scl+1.0_8))
    dy0 = 2.0_8*(p2%y-p1%y)/(n*(scl+1.0_8))
    dx = dx0
    dy = dy0
    do i = 2,n
      pts(i)%x = pts(i-1)%x + dx
      pts(i)%y = pts(i-1)%y + dy
      if (i.eq.n) exit
      dx = dx0*(1.0_8 + (scl-1.0_8)*(i-1)/(n-1))
      dy = dy0*(1.0_8 + (scl-1.0_8)*(i-1)/(n-1))
    end do
    return
  end subroutine divide_straight_line_2d_scaling
  !****************************************************
  !****************************************************
  subroutine divide_straight_line_2d_semi_scaling(p1,p2,n,scl,nsd,pts)
    !----------------------------------------------------
    ! Divide a straight line (between points p1 and p2)
    ! in "n" segments, scaling its size. The last segment
    ! has a size equal the first multiplied by parameter
    ! "scl". Each segment is also divided into "nsd"
    ! subdivisions, equally spaced. Return array "pts"
    ! with the points that delimit these segments.
    !
    ! Date: 20 December 2021
    ! Last update: ---
    !----------------------------------------------------
    type(point_2d),intent(in) :: p1,p2
    real(kind=8),intent(in) :: scl
    integer,intent(in) :: n,nsd
    type(point_2d),allocatable,intent(out) :: pts(:)
    type(point_2d),allocatable :: pt(:)
    real(kind=8) :: dx,dy,dx0,dy0
    integer :: i
    !
    if (allocated(pts)) deallocate(pts)
    allocate(pts(n*nsd+1))
    !
    pts(1) = p1
    pts(n*nsd+1) = p2
    dx0 = 2.0_8*(p2%x-p1%x)/(n*(scl+1.0_8))
    dy0 = 2.0_8*(p2%y-p1%y)/(n*(scl+1.0_8))
    dx = dx0
    dy = dy0
    do i = 2,n+1
      pts(nsd*(i-1)+1)%x = pts(nsd*(i-2)+1)%x + dx
      pts(nsd*(i-1)+1)%y = pts(nsd*(i-2)+1)%y + dy
      call divide_straight_line_2d(pts(nsd*(i-2)+1),pts(nsd*(i-1)+1),nsd,pt)
      pts(nsd*(i-2)+1:nsd*(i-1)+1) = pt
      if (i.eq.n+1) exit
      dx = dx0*(1.0_8 + (scl-1.0_8)*(i-1)/(n-1))
      dy = dy0*(1.0_8 + (scl-1.0_8)*(i-1)/(n-1))
    end do
    return
  end subroutine divide_straight_line_2d_semi_scaling
  !****************************************************
  !****************************************************
  subroutine divide_circular_line_2d(p1,p2,ang,n,pts,counter_clockwise)
    !----------------------------------------------------
    ! Divide a circular line (centered at point p1,
    ! beginning at point p2 and with a total angle of
    ! ang) in "n" segments with the same arc-length and
    ! returns array "pts" with the points that delimit
    ! these segments.
    !
    ! counter_clockwise = .true.
    ! theta goes counter-clockwise
    !
    ! counter_clockwise = .false.
    ! theta goes clockwise
    !
    ! Date: 07 February 2019
    ! Last update: 01 July 2021
    !----------------------------------------------------
    type(point_2d),intent(in) :: p1,p2
    real(kind=8),intent(in) :: ang
    integer,intent(in) :: n
    logical,intent(in) :: counter_clockwise
    type(point_2d),allocatable,intent(out) :: pts(:)
    real(kind=8) :: rad(3)
    real(kind=8) :: theta,d_theta
    integer :: i
    !
    if (allocated(pts)) deallocate(pts)
    allocate(pts(n+1))
    pts(1) = p2
    !
    rad = distance_2d(p1,p2)
    theta = atan2(rad(3),rad(2)) !Initial angle
    if (counter_clockwise) then
      d_theta = ang/n
    else
      d_theta = -ang/n
    end if
    do i = 2,n+1
      theta = theta + d_theta
      pts(i)%x = p1%x + rad(1)*cos(theta)
      pts(i)%y = p1%y + rad(1)*sin(theta)
    end do
    return
  end subroutine divide_circular_line_2d
  !****************************************************
  !****************************************************
  subroutine divide_circular_line_2d_scaling(p1,p2,ang,scl,n,pts,cc)
    !----------------------------------------------------
    ! Divide a circular line (centered at point p1,
    ! beginning at point p2 and with a total angle of
    ! ang) in "n" segments, scaling its arc-length, and
    ! returns array "pts" with the points that delimit
    ! these segments.
    !
    ! cc = counter_clockwise
    !
    ! counter_clockwise = .true.
    ! theta goes counter-clockwise
    !
    ! counter_clockwise = .false.
    ! theta goes clockwise
    !
    ! The last arc has a size equal the first
    ! multiplied by parameter "scl".
    !
    ! Date: 20 December 2021
    ! Last update: ---
    !----------------------------------------------------
    type(point_2d),intent(in) :: p1,p2
    real(kind=8),intent(in) :: ang,scl
    integer,intent(in) :: n
    logical,intent(in) :: cc
    type(point_2d),allocatable,intent(out) :: pts(:)
    real(kind=8) :: rad(3)
    real(kind=8) :: theta,d_theta0,d_theta,sig
    integer :: i
    !
    if (allocated(pts)) deallocate(pts)
    allocate(pts(n+1))
    pts(1) = p2
    !
    rad = distance_2d(p1,p2)
    theta = atan2(rad(3),rad(2)) !Initial angle
    !
    if (cc) then
      sig = +1._8
    else
      sig = -1._8
    end if
    !
    d_theta0 = sig*2.0_8*ang/(n*(scl+1.0_8))
    d_theta = d_theta0
    do i = 2,n+1
      theta = theta + d_theta
      pts(i)%x = p1%x + rad(1)*cos(theta)
      pts(i)%y = p1%y + rad(1)*sin(theta)
      if (i.eq.n+1) exit
      d_theta = d_theta0*(1.0_8 + (scl-1.0_8)*(i-1)/(n-1))
    end do
    !
    return
  end subroutine divide_circular_line_2d_scaling
  !****************************************************
  !****************************************************
  subroutine divide_circular_line_2d_semi_scaling(p1,p2,ang,scl,n,nsd,pts,cc)
    !----------------------------------------------------
    ! Divide a circular line (centered at point p1,
    ! beginning at point p2 and with a total angle of
    ! ang) in "n" segments, scaling its arc-length. Each
    ! segment is also divided into "nsd" subdivisions with
    ! the same arc-length. Returns array "pts" with the
    ! points that delimit these segments.
    !
    ! cc = counter_clockwise
    !
    ! counter_clockwise = .true.
    ! theta goes counter-clockwise
    !
    ! counter_clockwise = .false.
    ! theta goes clockwise
    !
    ! The last arc has a size equal the first
    ! multiplied by parameter "scl".
    !
    ! Date: 20 December 2021
    ! Last update: ---
    !----------------------------------------------------
    type(point_2d),intent(in) :: p1,p2
    real(kind=8),intent(in) :: ang,scl
    integer,intent(in) :: n,nsd
    logical,intent(in) :: cc
    type(point_2d),allocatable,intent(out) :: pts(:)
    type(point_2d),allocatable :: pt(:)
    real(kind=8) :: rad(3)
    real(kind=8) :: theta,d_theta0,d_theta,sig
    integer :: i
    !
    if (allocated(pts)) deallocate(pts)
    allocate(pts(n*nsd+1))
    pts(1) = p2
    !
    rad = distance_2d(p1,p2)
    theta = atan2(rad(3),rad(2)) !Initial angle
    !
    if (cc) then
      sig = +1._8
    else
      sig = -1._8
    end if
    !
    d_theta0 = sig*2.0_8*ang/(n*(scl+1.0_8))
    d_theta = d_theta0
    do i = 2,n+1
      theta = theta + d_theta
      pts(nsd*(i-1)+1)%x = p1%x + rad(1)*cos(theta)
      pts(nsd*(i-1)+1)%y = p1%y + rad(1)*sin(theta)
      call divide_circular_line_2d(p1,pts(nsd*(i-2)+1),abs(d_theta),nsd,pt,cc)
      pts(nsd*(i-2)+1:nsd*(i-1)+1) = pt
      if (i.eq.n+1) exit
      d_theta = d_theta0*(1.0_8 + (scl-1.0_8)*(i-1)/(n-1))
    end do
    !
    return
  end subroutine divide_circular_line_2d_semi_scaling
  !****************************************************
  !****************************************************
  function distance_3d(p1,p2)
    !----------------------------------------------------
    ! Compute the distance between points p1 and p2.
    ! p1 and p2 are of point_3d type.
    !
    ! Returns a vector with:
    ! distance_2d(1) = r (the distance)
    ! distance_2d(2) = rx (the x projection)
    ! distance_2d(3) = ry (the y projection)
    ! distance_2d(4) = rz (the z projection)
    !
    ! Date: 13 July 2018
    ! Last update: ---
    !----------------------------------------------------
    type(point_3d),intent(in) :: p1,p2
    real(kind=8) :: distance_3d(4)
    distance_3d(2)=(p2%x-p1%x)
    distance_3d(3)=(p2%y-p1%y)
    distance_3d(4)=(p2%z-p1%z)
    distance_3d(1)=sqrt(distance_3d(2)**2+distance_3d(3)**2+distance_3d(4)**2)
    return
  end function distance_3d
  !****************************************************
  !****************************************************
end module geometry
