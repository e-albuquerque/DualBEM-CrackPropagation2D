module fracture_2d_create_internal_points
!----------------------------------------------------
! A module with subroutines that creates internal
! points in front of the crack.
!
! Date: January 2023
! Last update: ---
!----------------------------------------------------
use fracture_2d_model
use numerical_integration_points
!
implicit none
!
integer :: b, n_ip
real(kind=8),private :: pi = 3.1415926535897932_8
type(internal_node), allocatable :: in_list_original(:)
!
contains
  subroutine create_internal_points(nip,initial_node_num)
    !----------------------------------------------------
    ! Generate points along a circular line that will be
    ! used as internal points.
    !
    ! Date: 17 January 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: nip,initial_node_num
    type(point_2d) :: p1,p2
    type(point_2d),allocatable :: pts(:)
    type(internal_node),allocatable :: in_list_temp(:)
    real(kind=8),allocatable :: eta(:),ome(:),r
    integer :: i,m,n,j,k,l,aux1,initial_node_pos
    real(kind=8) :: rad(3)
    real(kind=8) :: theta,theta_initial
    type(boundary_elem) :: last_elem
    !
    b = 0
    n_ip = nip
    allocate(in_list_original(size(in_list(:))))
    in_list_original = in_list
    do k = 1,cr_num
        do j = 1,cr_list(k)%tips_num
            b = b + 1
            if (allocated(pts)) deallocate(pts)
            allocate(pts(nip))
            !
            initial_node_pos = 0
            !
            if (j.eq.1) then
                aux1 = size(cr_list(k)%nodes_num(:,1))-initial_node_num+1
                initial_node_pos = cr_list(k)%nodes_num(aux1,2)
                last_elem = be_list(cr_list(k)%ext_el_num(2))
                else if (j.eq.2) then
                    initial_node_pos = cr_list(k)%nodes_num(initial_node_num,1)
                    last_elem = be_list(cr_list(k)%ext_el_num(4))
            end if
            p1 = last_elem%bp_2
            p2 = bn_list(initial_node_pos)%point
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
            m = size(in_list(:))
            allocate(in_list_temp(m))
            in_list_temp = in_list
            deallocate(in_list)
            n = m + nip
            allocate(in_list(n))
            in_list(1:m) = in_list_temp
            l = 1
            do i = m+1, n
                in_list(i)%point = pts(l)
                l = l + 1
            end do
            deallocate(in_list_temp)
            deallocate(eta)
            deallocate(ome)
        end do
    end do
    in_num = size(in_list(:))
    !
    return
  end subroutine create_internal_points
  !****************************************************
  !****************************************************
  subroutine original_internal_points()
    !----------------------------------------------------
    ! Return the original in_list.
    !
    ! Date: 17 January 2022
    ! Last update: ---
    !----------------------------------------------------
    integer :: i
    !
    i = size(in_list_original(:))
    deallocate(in_list)
    allocate(in_list(i))
    in_list = in_list_original
    deallocate(in_list_original)
    in_num = size(in_list(:))
    !
    return
  end subroutine original_internal_points
  !****************************************************
  !****************************************************
end module fracture_2d_create_internal_points
