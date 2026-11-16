module fracture_2d_input
!----------------------------------------------------
! A module with routines to read input data and
! fill the model for two-dimensional fracture
! analysis.
!
! Date: August 2022
! Last update: ---
!----------------------------------------------------
use error_messages
use geometry
use interpolation_functions
use fracture_2d_model
use fracture_2d_rad_exp_intgr
use fracture_2d_stress_intens_factor
use fracture_2d_crack_propagation
use fracture_2d_rainflow
!
implicit none
!
integer,private :: max_num_bn = 100000  !Maximum number of boundary nodes
integer,private :: max_num_be = 100000  !Maximum number of boundary elements
integer,private :: max_num_in = 100000  !Maximum number of internal nodes
!
integer,private :: tnbe !Total number of boundary elements
integer,private :: tnbn !Total number of boundary nodes
integer,private :: tnip !Total number of internal points
!
type(boundary_elem),allocatable,private :: lbe(:) !Local list of boundary elements
type(boundary_node),allocatable,private :: lbn(:) !Local list of boundary nodes
type(internal_node),allocatable,private :: lip(:) !Local list of internal points
!
!
contains
  subroutine input_from_text_file(file_name,steps,variable_ampl)
    !----------------------------------------------------
    ! Read the input data from a given text file and
    ! fill the model.
    !
    ! Date: 12 August 2022
    ! Last update: ---
    !----------------------------------------------------
    character(len=*),intent(in) :: file_name
    integer,intent(out) :: steps
    logical,intent(out) :: variable_ampl
    character(len=200) :: an_name
    character(len=200) :: an_id
    real(kind=8) :: el_mod,p_ratio
    integer :: psk,npow
    integer :: ppbc,ps
    integer :: bk,ck,dk,nbe,nbn,bc_dir,tips_num,number_of_cracks
    integer :: bck(2),bcdf(2)
    integer :: ipk,nip,niph,nipv
    real(kind=8) :: x1,y1,x2,y2,lb1,lb2
    real(kind=8) :: bi(2,4)
    real(kind=8) :: ai(4)
    real(kind=8) :: ang,scl,sclh,sclv
    logical :: cip !connected initial point
    logical :: flcb !First line of a closed boundary
    logical :: crack_key != true if it is a crack, = false if it is boundary
    type(point_2d),allocatable :: bep(:)
    type(point_2d) :: first_pt
    integer :: first_pt_pos
    integer :: i,j,K_key,k,n,nd,initial_node,k_eff_key,vari_ampl_key
    real(kind=8) :: R,C0,m,gama,K_th,K_IC,da_0,da_max,r_da,a
    !
    open(10,file=file_name,status="old") !Opening the input file
    !
    !----------------------------------------
    ! Reading and storing model name and id:
    !----------------------------------------
    !
    100 format(200(a))
    !
    read(10,*) !Jump a line in the input file
    read(10,*) !Jump a line in the input file
    read(10,100) an_name
    read(10,*) !Jump a line in the input file
    read(10,100) an_id
    !
    call set_model_name(an_name)
    call set_model_id(an_id)
    !
    print*
    print*,"Model name: ",an_name
    print*,"Constructing the model... "
    !
    !------------------------------------------
    ! Reading and storing material properties:
    !------------------------------------------
    !
    read(10,*) !Jump a line in the input file
    read(10,*) !Jump a line in the input file
    read(10,*) el_mod
    read(10,*) !Jump a line in the input file
    read(10,*) p_ratio
    !
    call set_young_modulus(el_mod)
    call set_poisson_ratio(p_ratio)
    !
    !--------------------------
    ! Reading plane state key:
    !--------------------------
    !
    read(10,*) !Jump a line in the input file
    read(10,*) !Jump a line in the input file
    read(10,*) psk
    !
    call set_plane_state_key(psk)
    !
    !------------------------
    ! Reading option of BIE:
    !------------------------
    !
    read(10,*) !Jump a line in the input file
    read(10,*) !Jump a line in the input file
    read(10,*) npow
    !
    call set_npow(npow)
    !
    !------------------------
    ! Reading boundary data:
    !------------------------
    !
    tnbe = 0
    tnbn = 0
    !
    allocate(lbe(max_num_be))
    allocate(lbn(max_num_bn))
    !
    read(10,200) !Jump 86 lines in the input file
    200 format(85(/))
    !
    !Loop on information to generate each boundary line:
    !
    flcb = .true.
    cip = .false.
    !
    do
      !
      read(10,*) bk
      if (bk.eq.0) exit !No more lines
      !
      read(10,*) ck
      !
      if (flcb .and. ck.ne.0) cip = .true.
      crack_key = .false.
      !
      select case (bk)
      !
      case(1)
        !
        !Straight line:
        !
        read(10,*) x1,y1,x2,y2,nbe,scl,nbn,lb1,lb2,bc_dir
        !
        call check_max_num_bound_elem_and_nodes(nbe,nbn)
        if (ck.ne.0 .and. nbn.gt.1) call check_connection_by_node_with_previous_line(x1,y1,lb1,flcb)
        !
        do i = 1,2
          read(10,*) bck(i),bcdf(i),bi(i,1),bi(i,2),bi(i,3),bi(i,4)
        end do
        !
        call straight_line(x1,y1,x2,y2,nbe,scl,nbn,lb1,lb2,bc_dir,bck,bcdf,bi,crack_key,0)
        !
        if (nbn.gt.1) then
          call check_connection_by_node_for_closed_bound(first_pt,first_pt_pos,flcb,cip,nbe,nbn)
        end if
        !
      case(2)
        !
        !Circular  line:
        !
        read(10,*) x1,y1,x2,y2,ang,dk,nbe,scl,nbn,lb1,lb2,bc_dir
        !
        call check_max_num_bound_elem_and_nodes(nbe,nbn)
        if (ck.ne.0 .and. nbn.gt.1) call check_connection_by_node_with_previous_line(x2,y2,lb1,flcb)
        !
        do i = 1,2
          read(10,*) bck(i),bcdf(i),bi(i,1),bi(i,2),bi(i,3),bi(i,4)
        end do
        !
        call circular_line(x1,y1,x2,y2,ang,dk,nbe,scl,nbn,lb1,lb2,bc_dir,bck,bcdf,bi,crack_key,0)
        !
        if (nbn.gt.1 .and. ang.lt.360._8) then
          call check_connection_by_node_for_closed_bound(first_pt,first_pt_pos,flcb,cip,nbe,nbn)
        end if
        !
      case(3)
        !
        !Cubic polynomial (y = a1 + a2*x + a3*x**2 + a4*x**3):
        !
        read(10,*) x1,x2,ai(1),ai(2),ai(3),ai(4),nbe,nbn,lb1,lb2,bc_dir
        !
        call check_max_num_bound_elem_and_nodes(nbe,nbn)
        if (ck.ne.0 .and. nbn.gt.1) then
          y1 = ai(1) + ai(2)*x1 + ai(3)*x1*x1 + ai(4)*x1*x1*x1
          call check_connection_by_node_with_previous_line(x1,y1,lb1,flcb)
        end if
        !
        do i = 1,2
          read(10,*) bck(i),bcdf(i),bi(i,1),bi(i,2),bi(i,3),bi(i,4)
        end do
        !
        call cubic_polynomial(x1,x2,nbe,nbn,lb1,lb2,bc_dir,bck,bcdf,ai,bi,crack_key,0)
        !
        if (nbn.gt.1) then
          call check_connection_by_node_for_closed_bound(first_pt,first_pt_pos,flcb,cip,nbe,nbn)
        end if
        !
      case(4)
        !
        !A given set of points which define the elements extreme points:
        !
        read(10,*) nbe,nbn,lb1,lb2,bc_dir
        !
        call check_max_num_bound_elem_and_nodes(nbe,nbn)
        !
        allocate(bep(nbe+1))
        do i = 1,nbe+1
          read(10,*) bep(i)%x,bep(i)%y
        end do
        !
        if (ck.ne.0 .and. nbn.gt.1) call check_connection_by_node_with_previous_line(bep(1)%x,bep(1)%y,lb1,flcb)
        !
        do i = 1,2
          read(10,*) bck(i),bcdf(i),bi(i,1),bi(i,2),bi(i,3),bi(i,4)
        end do
        !
        call given_points(bep,nbn,lb1,lb2,bc_dir,bck,bcdf,bi,crack_key,0)
        deallocate(bep)
        !
        if (nbn.gt.1) then
          call check_connection_by_node_for_closed_bound(first_pt,first_pt_pos,flcb,cip,nbe,nbn)
        end if
        !
      case default
        call error_message_stop("Undefined boundary key.")
      !
      end select
      !
    end do
    !
    !---------------------------------------------------
    ! Reading boundary conditions at particular points:
    !---------------------------------------------------
    !
    read(10,201) !Jump 22 lines in the input file
    201 format(21(/))
    !
    read(10,*) ppbc
    !
    if (ppbc.ne.0) then
      !
      do i = 1,ppbc
        read(10,*) x1,y1,ps
        read(10,*) bc_dir,bck(1),bck(2),bi(1,1),bi(1,2),bi(1,3),bi(1,4)
        call set_particular_bc(ps,x1,y1,bc_dir,bck,bi(1,:))
      end do
      !
    end if
    !
    !-------------------------------
    ! Reading crack data:
    !-------------------------------
    !
    cr_num = 0
    !
    read(10,202) !Jump 85 lines in the input file
    202 format(84(/))
    !
    read(10,*) number_of_cracks
    cr_num = number_of_cracks
    crack_key = .true.
    !
    !Loop on information to generate each crack line:
    !
    allocate(cr_list(cr_num))
    !
    do j = 1,cr_num
      !
      read(10,*) bk,tips_num
      cr_list(j)%tips_num = tips_num
      !
      select case (bk)
      !
      case(1)
        !
        !Straight line:
        !
        read(10,*) x1,y1,x2,y2,nbe,scl,nbn,lb1,lb2,bc_dir
        cr_list(j)%ext_el_num(1) = tnbe+1
        cr_list(j)%ext_el_num(2) = tnbe+nbe
        cr_list(j)%ext_el_num(3) = tnbe+nbe+1
        cr_list(j)%ext_el_num(4) = tnbe+2*nbe
        !        !
        do i = 1,2
          read(10,*) bck(i),bcdf(i),bi(i,1),bi(i,2),bi(i,3),bi(i,4)
        end do
        !
        call straight_line(x1,y1,x2,y2,nbe,scl,nbn,lb1,lb2,bc_dir,bck,bcdf,bi,crack_key,j)
        !
        cr_list(j)%length = 0.0_8
        do i = cr_list(j)%ext_el_num(1),cr_list(j)%ext_el_num(2)
            cr_list(j)%length = cr_list(j)%length + distance_only_2d(lbe(i)%bp_1,lbe(i)%bp_2)
        end do
        !
        !
      case(2)
        !
        !Circular  line:
        !
        read(10,*) x1,y1,x2,y2,ang,dk,nbe,scl,nbn,lb1,lb2,bc_dir
        cr_list(j)%ext_el_num(1) = tnbe+1
        cr_list(j)%ext_el_num(2) = tnbe+nbe
        cr_list(j)%ext_el_num(3) = tnbe+nbe+1
        cr_list(j)%ext_el_num(4) = tnbe+2*nbe
        !
        do i = 1,2
          read(10,*) bck(i),bcdf(i),bi(i,1),bi(i,2),bi(i,3),bi(i,4)
        end do
        !
        call circular_line(x1,y1,x2,y2,ang,dk,nbe,scl,nbn,lb1,lb2,bc_dir,bck,bcdf,bi,crack_key,j)
        !
        cr_list(j)%length = 0.0_8
        do i = cr_list(j)%ext_el_num(1),cr_list(j)%ext_el_num(2)
            cr_list(j)%length = cr_list(j)%length + distance_only_2d(lbe(i)%bp_1,lbe(i)%bp_2)
        end do
        !
      case(3)
        !
        !Cubic polynomial (y = a1 + a2*x + a3*x**2 + a4*x**3):
        !
        read(10,*) x1,x2,ai(1),ai(2),ai(3),ai(4),nbe,nbn,lb1,lb2,bc_dir
        cr_list(j)%ext_el_num(1) = tnbe+1
        cr_list(j)%ext_el_num(2) = tnbe+nbe
        cr_list(j)%ext_el_num(3) = tnbe+nbe+1
        cr_list(j)%ext_el_num(4) = tnbe+2*nbe
        !
        do i = 1,2
          read(10,*) bck(i),bcdf(i),bi(i,1),bi(i,2),bi(i,3),bi(i,4)
        end do
        !
        call cubic_polynomial(x1,x2,nbe,nbn,lb1,lb2,bc_dir,bck,bcdf,ai,bi,crack_key,j)
        !
        cr_list(j)%length = 0.0_8
        do i = cr_list(j)%ext_el_num(1),cr_list(j)%ext_el_num(2)
            cr_list(j)%length = cr_list(j)%length + distance_only_2d(lbe(i)%bp_1,lbe(i)%bp_2)
        end do
        !
      case(4)
        !
        !A given set of points which define the elements extreme points:
        !
        read(10,*) nbe,nbn,lb1,lb2,bc_dir
        cr_list(j)%ext_el_num(1) = tnbe+1
        cr_list(j)%ext_el_num(2) = tnbe+nbe
        cr_list(j)%ext_el_num(3) = tnbe+nbe+1
        cr_list(j)%ext_el_num(4) = tnbe+2*nbe
        !
        allocate(bep(nbe+1))
        do i = 1,nbe+1
          read(10,*) bep(i)%x,bep(i)%y
        end do
        !
        do i = 1,2
          read(10,*) bck(i),bcdf(i),bi(i,1),bi(i,2),bi(i,3),bi(i,4)
        end do
        !
        call given_points(bep,nbn,lb1,lb2,bc_dir,bck,bcdf,bi,crack_key,j)
        !
        deallocate(bep)
        !
        cr_list(j)%length = 0.0_8
        do i = cr_list(j)%ext_el_num(1),cr_list(j)%ext_el_num(2)
            cr_list(j)%length = cr_list(j)%length + distance_only_2d(lbe(i)%bp_1,lbe(i)%bp_2)
        end do
        !
      case default
        call error_message_stop("Undefined crack key.")
      !
      end select
      !
    end do
    !-------------------------------
    ! Reading internal points data:
    !-------------------------------
    !
    tnip = 0
    !
    allocate(lip(max_num_in))
    !
    read(10,203) !Jump 52 lines in the input file
    203 format(51(/))
    !
    !Loop on information to generate internal points:
    !
    do
      !
      read(10,*) ipk
      !
      select case (ipk)
      !
      case(0)
        !No more points:
        exit
        !
      case(1)
        !Individual points:
        read(10,*) nip
        call check_max_num_int_points(nip)
        do i = 1,nip
          read(10,*) x1,y1
          tnip = tnip + 1
          lip(tnip)%point%x = x1
          lip(tnip)%point%y = y1
        end do
        !
      case(2)
        !A straight line:
        read(10,*) x1,y1,x2,y2,nip,scl
        call check_max_num_int_points(nip)
        call int_pts_straight_line(x1,y1,x2,y2,nip,scl)
        !
      case(3)
        !A circular line:
        read(10,*) x1,y1,x2,y2,ang,dk,nip,scl
        call check_max_num_int_points(nip)
        call int_pts_circular_line(x1,y1,x2,y2,ang,dk,nip,scl)
        !
      case(4)
        !A rectangular area with distributed points:
        read(10,*) x1,y1,x2,y2,niph,nipv,sclh,sclv
        call check_max_num_int_points(niph*nipv)
        call int_pts_rectangular_area(x1,y1,x2,y2,niph,nipv,sclh,sclv)
        !
      case default
        call error_message_stop("Undefined internal points key.")
      !
      end select
      !
    end do
    !
    !--------------------------------------------------------
    ! Reading the stress intensity factor evaluation data:
    !--------------------------------------------------------
    read(10,204) !Jump 34 lines in the input file
    204 format(33(/))
    read(10,*) K_key
    call set_K_key(K_key)
    read(10,*)
    read(10,*) k
    call set_method_key(k)
    !
    if (K_key.eq.1) then
        if (k.eq.7) then
            read(10,*)
            read(10,*) nd
            call set_node(nd)
            else
                read(10,*)
                read(10,*)
        end if
    end if
    !
    if (K_key.eq.2) then
        if (k.eq.1) then
            read(10,*)
            read(10,*) n, initial_node
            call set_number_of_points(n)
            call set_initial_node_num(initial_node)
            else if (k.eq.2) then
                read(10,*)
                read(10,*) n, initial_node
                call set_number_of_points(n)
                call set_initial_node_num(initial_node)
            else if (k.eq.3) then
                read(10,*)
                read(10,*) nbe,nbn,lb1,lb2,n,initial_node
                call set_element_information(nbe,nbn,lb1,lb2)
                call set_number_of_points(n)
                call set_initial_node_num(initial_node)
            end if
    end if
    !
    !---------------------------------------------
    ! Reading the fatigue crack propagation data:
    !---------------------------------------------
    read(10,205) !Jump 10 lines in the input file
    205 format(9(/))
    read(10,*) k_eff_key
    call set_k_eff_key(k_eff_key)
    read(10,*)
    read(10,*) Da_0, r_da, Da_max, steps
    if (mod(steps,2) == 1) then
        steps = steps + 1
    end if
    call set_initial_delta_a(Da_0)
    call set_r_Da(r_da)
    call set_max_delta_a(Da_max)
    read(10,*)
    read(10,*) R,C0,m,gama,K_th,K_IC
    call set_parameter_R(R)
    call set_parameter_C(C0)
    call set_parameter_m(m)
    call set_parameter_gamma(gama)
    call set_K_th(K_th)
    call set_K_IC(K_IC)
    !
    !--------------------------------------
    ! Reading the variable amplitude data:
    !--------------------------------------
    read(10,206) !Jump 12 lines in the input file
    206 format(11(/))
    read(10,*) vari_ampl_key
    if (vari_ampl_key .eq. 0) then
        variable_ampl = .false.
        else if (vari_ampl_key .eq. 1) then
            variable_ampl = .true.
            else
                call error_message_stop("Incorrect variable amplitude key")
    end if
    if (variable_ampl) then
        read(10,*) n
        call allocate_alpha(n)
        do i = 1,n
            read(10,*) a
            call set_alpha(i,a)
        end do
    end if
    !
    !---------------------------------------
    ! Setting discrete parameters on model:
    !---------------------------------------
    !
    be_num = tnbe
    bn_num = tnbn
    in_num = tnip
    !
    allocate(be_list(be_num))
    allocate(bn_list(bn_num))
    allocate(in_list(in_num))
    !
    be_list = lbe(1:tnbe)
    bn_list = lbn(1:tnbn)
    in_list = lip(1:in_num)
    !
    print*,"                          ...done!"
    print*
    !
    close(10) !Closing the input file
    !
    print 120, be_num
    120 format(1x,"Number of boundary elements = ",i0)
    print 121, bn_num
    121 format(1x,"Number of boundary points = ",i0)
    if (in_num.gt.0) print 122, in_num
    122 format(1x,"Number of internal points = ",i0)
    print*
    !
    return
  end subroutine input_from_text_file
  !****************************************************
  !****************************************************
  subroutine straight_line(x1,y1,x2,y2,nbe,scl,nbn,lb1,lb2,bc_dir,bck,bcdf,bi,crack_key,i)
    !----------------------------------------------------
    ! Generate elements and boundary nodes for a straight
    ! line.
    !
    ! The extreme elements (the first and last element in
    ! the line) are semi-discontinuous elements with the
    ! discontinuous parameters lb1 and lb2, respectively.
    !
    ! Date: 27 August 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: nbe,nbn,bc_dir,i!,bie
    integer,intent(in) :: bck(2),bcdf(2)
    real(kind=8),intent(in) :: x1,y1,x2,y2,lb1,lb2,scl
    real(kind=8),intent(in) :: bi(2,4)
    logical,intent(in) :: crack_key
    type(point_2d) :: p1,p2
    type(point_2d),allocatable :: pts(:),pts2(:)
    integer :: j,k,l,n,bie,aux1,aux2,aux3,ck
    !
    !-------------------------------------------------------------
    ! Evaluating a list of points to generate elements and nodes:
    !-------------------------------------------------------------
    !
    p1%x = x1; p1%y = y1
    p2%x = x2; p2%y = y2
    !
    if(crack_key) then
        bie = 2
        ck = 2
        else
            bie = 1
            ck = 1
    end if
    !
    if (lb1.eq.1.0_8 .or. lb2.eq.1.0_8) then
        call error_message_stop("Invalid lambda 1 or lambda 2, generating a continuous boundary elements.")
    end if
    !
    if (nbn.lt.3) then
      call divide_straight_line_2d_scaling(p1,p2,nbe,scl,pts)
    else
      call divide_straight_line_2d_semi_scaling(p1,p2,nbe,scl,nbn-1,pts)
    end if
    aux1 = tnbn+1
    !
    !--------------------------------
    ! Generating elements and nodes:
    !--------------------------------
    !
    select case (nbn)
    case(1)
      call constant_elements(pts,bc_dir,bck,bcdf,bi,bie)
    case(2)
      call linear_elements(pts,lb1,lb2,bc_dir,bck,bcdf,bi,bie,ck)
    case(3)
      call quadratic_elements(pts,lb1,lb2,bc_dir,bck,bcdf,bi,bie,ck)
    case default
      call error_message_stop("Undefined element order.")
    end select
    !
    aux2 = tnbn
    !
    !-------------------------------------------
    ! Generating elements and nodes for cracks:
    !-------------------------------------------
    !
    if (crack_key) then
        aux3 = aux2-aux1+1
        allocate(cr_list(i)%nodes_num(aux3,2))
        n = size(pts)
        allocate(pts2(n))
        do j = 1,size(pts)
            pts2(n) = pts(j)
            n = n - 1
        end do
        bie = 3
        ck = 2
        select case (nbn)
            case(1)
                call constant_elements(pts2,bc_dir,bck,bcdf,bi,bie)
            case(2)
                call linear_elements(pts2,lb2,lb1,bc_dir,bck,bcdf,bi,bie,ck)
            case(3)
                call quadratic_elements(pts2,lb2,lb1,bc_dir,bck,bcdf,bi,bie,ck)
            case default
                call error_message_stop("Undefined element order.")
        end select
        j = 0
        l = aux1
        do k = 1, aux3
            cr_list(i)%nodes_num(k,1) = l
            cr_list(i)%nodes_num(k,2) = tnbn-j
            j = j+1
            l = l+1
        end do
    end if
    !
    return
  end subroutine straight_line
  !****************************************************
  !****************************************************
  subroutine circular_line(x1,y1,x2,y2,ang,dk,nbe,scl,nbn,lb1,lb2,bc_dir,bck,bcdf,bi,crack_key,i)
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
    integer,intent(in) :: dk,nbe,nbn,bc_dir,i
    integer,intent(in) :: bck(2),bcdf(2)
    real(kind=8),intent(in) :: x1,y1,x2,y2,lb1,lb2,ang,scl
    real(kind=8),intent(in) :: bi(2,4)
    logical,intent(in) :: crack_key
    real(kind=8) :: ang_rad
    logical :: cc
    type(point_2d) :: p1,p2
    type(point_2d),allocatable :: pts(:),pts2(:)
    integer :: first_node,bie,j,k,l,n,aux1,aux2,aux3,ck
    !
    real(kind=8) :: tol = 1.0e-8_8
    !
    real(kind=8) :: pi = 3.1415926535897932_8
    !
    !-------------------------
    ! Checking maximum angle:
    !-------------------------
    !
    if (ang.gt.360.) then
      call error_message_stop("Trying to create an arc grater than 360 degrees.")
    end if
    !
    !if (lb1.eq.1.0_8 .or. lb2.eq.1.0_8) then
    !    call error_message_stop("Invalid lambda 1 or lambda 2, generating a continuous boundary elements.")
    !end if
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
    p1%x = x1; p1%y = y1
    p2%x = x2; p2%y = y2
    !
    if(crack_key) then
        bie = 2
        ck = 2
        else
            bie = 1
            ck = 1
    end if
    !
    if (nbn.lt.3) then
      call divide_circular_line_2d_scaling(p1,p2,ang_rad,scl,nbe,pts,cc)
    else
      call divide_circular_line_2d_semi_scaling(p1,p2,ang_rad,scl,nbe,nbn-1,pts,cc)
    end if
    !
    aux1 = tnbn+1
    !
    !--------------------------------
    ! Generating elements and nodes:
    !--------------------------------
    !
    first_node = tnbn + 1
    !
    select case (nbn)
    case(1)
      call constant_elements(pts,bc_dir,bck,bcdf,bi,bie)
    case(2)
      call linear_elements(pts,lb1,lb2,bc_dir,bck,bcdf,bi,bie,ck)
    case(3)
      call quadratic_elements(pts,lb1,lb2,bc_dir,bck,bcdf,bi,bie,ck)
    case default
      call error_message_stop("Undefined element order.")
    end select
    !
    aux2 = tnbn
    !
    !-------------------------------------------
    ! Generating elements and nodes for cracks:
    !-------------------------------------------
    !
    if (crack_key) then
        aux3 = aux2-aux1+1
        allocate(cr_list(i)%nodes_num(aux3,2))
        n = size(pts)
        allocate(pts2(n))
        do j = 1,size(pts)
            pts2(n) = pts(j)
            n = n - 1
        end do
        bie = 3
        ck = 2
        select case (nbn)
            case(1)
                call constant_elements(pts2,bc_dir,bck,bcdf,bi,bie)
            case(2)
                call linear_elements(pts2,lb2,lb1,bc_dir,bck,bcdf,bi,bie,ck)
            case(3)
                call quadratic_elements(pts2,lb2,lb1,bc_dir,bck,bcdf,bi,bie,ck)
            case default
                call error_message_stop("Undefined element order.")
        end select
        j = 0
        l = aux1
        do k = 1, aux3
            cr_list(i)%nodes_num(k,1) = l
            cr_list(i)%nodes_num(k,2) = tnbn-j
            j = j+1
            l = l+1
        end do
    end if
    !
    !------------------------------------
    ! Checking for closed circumference:
    !------------------------------------
    !
    if (nbn.ne.1 .and. ang.gt.(360._8-tol) .and. &
      & lb1.gt.(1._8-tol) .and. lb2.gt.(1._8-tol)) then
      lbn(tnbn)%bc = 0
      lbn(tnbn)%bc_dir = 0
      lbn(tnbn)%disp = 0.0
      lbn(tnbn)%trac = 0.0
      lbn(tnbn)%h1 = 0.0
      lbn(tnbn)%h2 = 0.0
      lbn(tnbn)%strain = 0.0
      lbn(tnbn)%stress = 0.0
      lbn(tnbn)%point%x = 0.0
      lbn(tnbn)%point%y = 0.0
      tnbn = tnbn - 1
      lbe(tnbe)%incid(lbe(tnbe)%nodes_num) = first_node
      lbe(tnbe)%bp_2 = lbn(first_node)%point
    end if
    !
    return
  end subroutine circular_line
  !****************************************************
  !****************************************************
  subroutine cubic_polynomial(x1,x2,nbe,nbn,lb1,lb2,bc_dir,bck,bcdf,ai,bi,crack_key,m)
    !----------------------------------------------------
    ! Generate elements and boundary nodes for a line
    ! defined by the following cubic polynomial:
    !
    ! y = a1 + a2*x + a3*x**2 + a4*x**3
    !
    ! The extreme elements (the first and last element in
    ! the line) are semi-discontinuous elements with the
    ! discontinuous parameters lb1 and lb2, respectively.
    !
    ! Date: 27 August 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: nbe,nbn,bc_dir,m
    integer,intent(in) :: bck(2),bcdf(2)
    real(kind=8),intent(in) :: x1,x2,lb1,lb2
    real(kind=8),intent(in) :: ai(4)
    logical,intent(in) :: crack_key
    real(kind=8),intent(inout) :: bi(2,4)
    real(kind=8) :: dx,x,y
    integer :: n_seg
    integer :: i,j,k,l,n,bie,aux1,aux2,aux3,ck
    type(point_2d),allocatable :: pts(:),pts2(:)
    !
    if (lb1.eq.1.0_8 .or. lb2.eq.1.0_8) then
        call error_message_stop("Invalid lambda 1 or lambda 2, generating a continuous boundary elements.")
    end if
    !
    !-------------------------------------------------------------
    ! Evaluating a list of points to generate elements and nodes:
    !-------------------------------------------------------------
    !
    if (nbn.lt.3) then
      n_seg = nbe
    else
      n_seg = (nbn-1)*nbe
    end if
    !
    allocate(pts(n_seg+1))
    !
    dx = (x2-x1)/n_seg
    x = x1 - dx
    !
    do i = 1,n_seg+1
      x = x + dx
      y = ai(1) + ai(2)*x + ai(3)*x**2 + ai(4)*x**3
      pts(i)%x = x
      pts(i)%y = y
    end do
    !
    !--------------------------------
    ! Generating elements and nodes:
    !--------------------------------
    !
    if(crack_key) then
        bie = 2
        ck = 2
        else
            bie = 1
            ck = 1
    end if
    !
    aux1 = tnbn+1
    !
    select case (nbn)
    case(1)
      call constant_elements(pts,bc_dir,bck,bcdf,bi,bie)
    case(2)
      call linear_elements(pts,lb1,lb2,bc_dir,bck,bcdf,bi,bie,ck)
    case(3)
      call quadratic_elements(pts,lb1,lb2,bc_dir,bck,bcdf,bi,bie,ck)
    case default
      call error_message_stop("Undefined element order.")
    end select
    !
    aux2 = tnbn
    !
    !-------------------------------------------
    ! Generating elements and nodes for cracks:
    !-------------------------------------------
    !
    if (crack_key) then
        aux3 = aux2-aux1+1
        allocate(cr_list(m)%nodes_num(aux3,2))
        n = size(pts)
        allocate(pts2(n))
        do j = 1,size(pts)
            pts2(n) = pts(j)
            n = n - 1
        end do
        bie = 3
        ck = 2
        select case (nbn)
            case(1)
                call constant_elements(pts2,bc_dir,bck,bcdf,bi,bie)
            case(2)
                call linear_elements(pts2,lb2,lb1,bc_dir,bck,bcdf,bi,bie,ck)
            case(3)
                call quadratic_elements(pts2,lb2,lb1,bc_dir,bck,bcdf,bi,bie,ck)
            case default
                call error_message_stop("Undefined element order.")
        end select
        j = 0
        l = aux1
        do k = 1, aux3
            cr_list(m)%nodes_num(k,1) = l
            cr_list(m)%nodes_num(k,2) = tnbn-j
            j = j+1
            l = l+1
        end do
    end if
    !
    return
  end subroutine cubic_polynomial
  !****************************************************
  !****************************************************
  subroutine given_points(bep,nbn,lb1,lb2,bc_dir,bck,bcdf,bi,crack_key,b)
    !----------------------------------------------------
    ! Generate elements and boundary nodes from a given
    ! set of points which represent the extreme points
    ! of the boundary elements.
    !
    ! The extreme elements (the first and last element in
    ! the line) are semi-discontinuous elements with the
    ! discontinuous parameters lb1 and lb2, respectively.
    !
    ! Date: 27 August 2022
    ! Last update: ---
    !----------------------------------------------------
    type(point_2d),intent(in) :: bep(:)
    integer,intent(in) :: nbn,bc_dir,b
    integer,intent(in) :: bck(2),bcdf(2)
    real(kind=8),intent(in) :: lb1,lb2
    real(kind=8),intent(in) :: bi(2,4)
    logical,intent(in) :: crack_key
    type(point_2d),allocatable :: pts(:),bep2(:)
    integer :: i,j,k,l,n,bie,aux1,aux2,aux3,ck
    !
    if (lb1.eq.1.0_8 .or. lb2.eq.1.0_8) then
        call error_message_stop("Invalid lambda 1 or lambda 2, generating a continuous boundary elements.")
    end if
    !
    if(crack_key) then
        bie = 2
        ck = 2
        else
            bie = 1
            ck = 1
    end if
    !
    aux1 = tnbn+1
    !
    select case (nbn)
    case(1)
      call constant_elements(bep,bc_dir,bck,bcdf,bi,bie)
    case(2)
      call linear_elements(bep,lb1,lb2,bc_dir,bck,bcdf,bi,bie,ck)
    case(3)
      allocate(pts(2*size(bep)-1))
      pts(1) = bep(1)
      do i = 1,size(bep)-1
        pts(2*i+1) = bep(i+1)
        pts(2*i) = midpoint_2d(pts(2*i-1),pts(2*i+1))
      end do
      call quadratic_elements(pts,lb1,lb2,bc_dir,bck,bcdf,bi,bie,ck)
      deallocate(pts)
    case default
      call error_message_stop("Undefined element order.")
    end select
    !
    aux2 = tnbn
    !
    !-------------------------------------------
    ! Generating elements and nodes for cracks:
    !-------------------------------------------
    !
    if (crack_key) then
        aux3 = aux2-aux1+1
        allocate(cr_list(b)%nodes_num(aux3,2))
        n = size(bep)
        allocate (bep2(n))
        do j = 1,size(bep)
            bep2(n) = bep(j)
            n = n - 1
        end do
        bie = 3
        ck = 2
        select case (nbn)
            case(1)
                call constant_elements(bep2,bc_dir,bck,bcdf,bi,bie)
            case(2)
                call linear_elements(bep2,lb2,lb1,bc_dir,bck,bcdf,bi,bie,ck)
            case(3)
                allocate(pts(2*size(bep2)-1))
                pts(1) = bep2(1)
                do i = 1,size(bep2)-1
                    pts(2*i+1) = bep2(i+1)
                    pts(2*i) = midpoint_2d(pts(2*i-1),pts(2*i+1))
                end do
                call quadratic_elements(pts,lb2,lb1,bc_dir,bck,bcdf,bi,bie,ck)
                deallocate(pts)
            case default
                call error_message_stop("Undefined element order.")
        end select
        j = 0
        l = aux1
        do k = 1, aux3
            cr_list(b)%nodes_num(k,1) = l
            cr_list(b)%nodes_num(k,2) = tnbn-j
            j = j+1
            l = l+1
        end do
    end if
    !
    return
  end subroutine given_points
  !****************************************************
  !****************************************************
  subroutine constant_elements(pts,bc_dir,bck,bcdf,bi,bie)
    !----------------------------------------------------
    ! Generate constant elements for a given list of
    ! points.
    !
    ! Date: 12 August 2022
    ! Last update: ---
    !----------------------------------------------------
    type(point_2d),intent(in) :: pts(:)
    integer,intent(in) :: bc_dir,bie
    integer,intent(in) :: bck(2),bcdf(2)
    real(kind=8),intent(in) :: bi(2,4)
    real(kind=8),allocatable :: s(:) !Parametric coordinates for boundary nodes s = [0,1]
    type(point_2d) :: p1,p2
    integer :: i
    !
    !--------------------------------------------------------------
    ! Evaluating the parametric coordinate for each boundary node:
    !--------------------------------------------------------------
    !
    call get_parametric_coords(pts,1.0_8,1.0_8,1,s)
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
      tnbe = tnbe + 1 !Adding a boundary element to the model
      !
      lbe(tnbe)%bp_1 = p1
      lbe(tnbe)%bp_2 = p2
      allocate(lbe(tnbe)%incid(1))
      lbe(tnbe)%incid(1) = tnbe
      lbe(tnbe)%lamb_1 = 1.0_8
      lbe(tnbe)%lamb_2 = 1.0_8
      lbe(tnbe)%nodes_num = 1
      !
      tnbn = tnbn + 1 !Adding a boundary node to the model
      !
      lbn(tnbn)%point = midpoint_2d(p1,p2)
      lbn(tnbn)%bie = bie
      call set_boundary_condition_on_node(lbn(tnbn),bc_dir,bck,bcdf,s(i),bi)
      !
    end do
    !
    return
  end subroutine constant_elements
  !****************************************************
  !****************************************************
  subroutine linear_elements(pts,lb1,lb2,bc_dir,bck,bcdf,bi,bie,ck)
    !----------------------------------------------------
    ! Generate linear elements for a given list of
    ! points.
    !
    ! Date: 12 August 2022
    ! Last update: ---
    !----------------------------------------------------
    type(point_2d),intent(in) :: pts(:)
    integer,intent(in) :: bc_dir,bie,ck
    integer,intent(in) :: bck(2),bcdf(2)
    real(kind=8),intent(in) :: lb1,lb2
    real(kind=8),intent(in) :: bi(2,4)
    real(kind=8),allocatable :: s(:) !Parametric coordinates for boundary nodes s = [0,1]
    type(point_2d) :: p1,p2
    type(point_2d) :: el_pts(2)
    integer :: i
    !
    !--------------------------------------------------------------
    ! Evaluating the parametric coordinate for each boundary node:
    !--------------------------------------------------------------
    !
    call get_parametric_coords(pts,lb1,lb2,2,s)
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
      tnbe = tnbe + 1 !Adding a boundary element to the model
      !
      lbe(tnbe)%bp_1 = p1
      lbe(tnbe)%bp_2 = p2
      lbe(tnbe)%nodes_num = 2
      allocate(lbe(tnbe)%incid(2))
      if (ck.eq.1) then
        if (i.eq.1) then
            lbe(tnbe)%lamb_1 = lb1
        else
            lbe(tnbe)%lamb_1 = 1.0_8
        end if
        if (i.eq.size(pts)-1) then
            lbe(tnbe)%lamb_2 = lb2
        else
            lbe(tnbe)%lamb_2 = 1.0_8
        end if
      end if
      if (ck.eq.2) then
        lbe(tnbe)%lamb_1 = lb1
        lbe(tnbe)%lamb_2 = lb2
      end if
      !
      if (ck.eq.2) then
        if (i.eq.1) then
            tnbn = tnbn + 1 !Adding a boundary node to the model
            el_pts(1) = p1
            el_pts(2) = p2
            lbn(tnbn)%point = cartesian_point_from_intrp_func_1d(el_pts,-lb1)
            lbn(tnbn)%bie = bie
            call set_boundary_condition_on_node(lbn(tnbn),bc_dir,bck,bcdf,s(i),bi)
        end if
        else if (ck.eq.2) then
            tnbn = tnbn + 1 !Adding a boundary node to the model
            el_pts(1) = p1
            el_pts(2) = p2
            lbn(tnbn)%point = cartesian_point_from_intrp_func_1d(el_pts,-lb1)
            lbn(tnbn)%bie = bie
            call set_boundary_condition_on_node(lbn(tnbn),bc_dir,bck,bcdf,s(i),bi)
      end if
      !
      lbe(tnbe)%incid(1) = tnbn
      !
      tnbn = tnbn + 1 !Adding a boundary node to the model
      if (ck.eq.1) then
        if (i.eq.size(pts)-1 .and. ck.eq.1) then
            el_pts(1) = p1
            el_pts(2) = p2
            lbn(tnbn)%point = cartesian_point_from_intrp_func_1d(el_pts,+lb2)
        else
            lbn(tnbn)%point = p2
        end if
      else if (ck.eq.2) then
            el_pts(1) = p1
            el_pts(2) = p2
            lbn(tnbn)%point = cartesian_point_from_intrp_func_1d(el_pts,+lb2)
      end if
      lbe(tnbe)%incid(2) = tnbn
      lbn(tnbn)%bie = bie
      call set_boundary_condition_on_node(lbn(tnbn),bc_dir,bck,bcdf,s(i+1),bi)
      !
    end do
    !
    return
  end subroutine linear_elements
  !****************************************************
  !****************************************************
  subroutine quadratic_elements(pts,lb1,lb2,bc_dir,bck,bcdf,bi,bie,ck)
    !----------------------------------------------------
    ! Generate quadratic elements for a given list of
    ! points.
    !
    ! Date: 12 August 2022
    ! Last update: ---
    !----------------------------------------------------
    type(point_2d),intent(in) :: pts(:)
    integer,intent(in) :: bc_dir,bie,ck
    integer,intent(in) :: bck(2),bcdf(2)
    real(kind=8),intent(in) :: lb1,lb2
    real(kind=8),intent(in) :: bi(2,4)
    real(kind=8),allocatable :: s(:) !Parametric coordinates for boundary nodes s = [0,1]
    type(point_2d) :: p1,p2,p3
    type(point_2d) :: el_pts(3)
    integer :: i,n
    !
    !--------------------------------------------------------------
    ! Evaluating the parametric coordinate for each boundary node:
    !--------------------------------------------------------------
    !
    call get_parametric_coords(pts,lb1,lb2,3,s)
    !
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
      tnbe = tnbe + 1 !Adding a boundary element to the model
      !
      lbe(tnbe)%bp_1 = p1
      lbe(tnbe)%bp_2 = p3
      lbe(tnbe)%nodes_num = 3
      allocate(lbe(tnbe)%incid(3))
      if (ck.eq.1) then
        if (i.eq.1) then
            lbe(tnbe)%lamb_1 = lb1
        else
            lbe(tnbe)%lamb_1 = 1.0_8
        end if
        if (i.eq.n) then
            lbe(tnbe)%lamb_2 = lb2
        else
            lbe(tnbe)%lamb_2 = 1.0_8
        end if
      end if
      if (ck.eq.2) then
        lbe(tnbe)%lamb_1 = lb1
        lbe(tnbe)%lamb_2 = lb2
      end if
      !
      if (i.eq.1 .and. ck.eq.1) then
        tnbn = tnbn + 1 !Adding a boundary node to the model
        el_pts(1) = p1
        el_pts(2) = p2
        el_pts(3) = p3
        lbn(tnbn)%point = cartesian_point_from_intrp_func_1d(el_pts,-lb1)
        lbn(tnbn)%bie = bie
        call set_boundary_condition_on_node(lbn(tnbn),bc_dir,bck,bcdf,s(i),bi)
        else if (ck.eq.2) then
            tnbn = tnbn + 1 !Adding a boundary node to the model
            el_pts(1) = p1
            el_pts(2) = p2
            el_pts(3) = p3
            lbn(tnbn)%point = cartesian_point_from_intrp_func_1d(el_pts,-lb1)
            lbn(tnbn)%bie = bie
            call set_boundary_condition_on_node(lbn(tnbn),bc_dir,bck,bcdf,s(i),bi)
      end if
      !
      lbe(tnbe)%incid(1) = tnbn
      !
      tnbn = tnbn + 1 !Adding a boundary node to the model
      lbe(tnbe)%incid(2) = tnbn
      lbn(tnbn)%point = p2
      lbn(tnbn)%bie = bie
      call set_boundary_condition_on_node(lbn(tnbn),bc_dir,bck,bcdf,s(2*i),bi)
      !
      tnbn = tnbn + 1 !Adding a boundary node to the model
      if (ck.eq.1) then
        if (i.eq.n) then
            el_pts(1) = p1
            el_pts(2) = p2
            el_pts(3) = p3
            lbn(tnbn)%point = cartesian_point_from_intrp_func_1d(el_pts,+lb2)
        else
            lbn(tnbn)%point = p3
        end if
        lbe(tnbe)%incid(3) = tnbn
        lbn(tnbn)%bie = bie
        call set_boundary_condition_on_node(lbn(tnbn),bc_dir,bck,bcdf,s(2*i+1),bi)
      else if (ck.eq.2) then
        lbn(tnbn)%point = cartesian_point_from_intrp_func_1d(el_pts,+lb2)
        lbe(tnbe)%incid(3) = tnbn
        lbn(tnbn)%bie = bie
        call set_boundary_condition_on_node(lbn(tnbn),bc_dir,bck,bcdf,s(2*i+1),bi)
      end if
      !
    end do
    !
    return
  end subroutine quadratic_elements
  !****************************************************
  !****************************************************
  subroutine get_parametric_coords(pts,lb1,lb2,nbn,s)
    !----------------------------------------------------
    ! Return a vector with approximate parametric
    ! coordinates, s = [0,1], for a given set of points.
    !
    ! Date: 27 June 2022
    ! Last update: ---
    !----------------------------------------------------
    type(point_2d),intent(in) :: pts(:)
    real(kind=8),intent(in) :: lb1,lb2
    integer,intent(in) :: nbn
    real(kind=8),allocatable,intent(out) :: s(:)
    real(kind=8) :: r0
    type(point_2d),allocatable :: pt(:)
    type(point_2d) :: pa,pb
    integer :: i,n
    !
    n = size(pts)
    !
    if (allocated(s)) deallocate(s)
    if (allocated(pt)) deallocate(pt)
    !
    select case (nbn)
    case(1)
      !Constant elements:
      allocate(s(n-1))
      pa = pts(1)
      r0 = 0.0
      do i = 1,n-1
        pb = midpoint_2d(pts(i),pts(i+1))
        r0 = r0 + distance_only_2d(pa,pb)
        s(i) = r0
        pa = pb
      end do
      r0 = r0 + distance_only_2d(pa,pts(n))
      s = s/r0
    case default
      !Linear or quadratic elements:
      allocate(s(n))
      allocate(pt(nbn))
      !
      do i = 1,nbn
        pt(i) = pts(i)
      end do
      pa = cartesian_point_from_intrp_func_1d(pt,-lb1)
      !
      r0 = distance_only_2d(pts(1),pa)
      s(1) = r0
      !
      do i = 2,n-1
        r0 = r0 + distance_only_2d(pa,pts(i))
        s(i) = r0
        pa = pts(i)
      end do
      !
      do i = 1,nbn
        pt(i) = pts(n-nbn+i)
      end do
      pb = cartesian_point_from_intrp_func_1d(pt,+lb2)
      !
      r0 = r0 + distance_only_2d(pts(n-1),pb)
      s(n) = r0
      r0 = r0 + distance_only_2d(pb,pts(n))
      s = s/r0
    end select
    !
    return
  end subroutine get_parametric_coords
  !****************************************************
  !****************************************************
  subroutine set_boundary_condition_on_node(bn,bc_dir,bck,bcdf,s,bi)
    !----------------------------------------------------
    ! Set boundary condition parameters on a given node.
    !
    ! Date: 27 June 2022
    ! Last update: ---
    !----------------------------------------------------
    type(boundary_node),intent(inout) :: bn
    integer,intent(in) :: bc_dir
    integer,intent(in) :: bck(2),bcdf(2)
    real(kind=8),intent(in) :: bi(2,4)
    real(kind=8),intent(in) :: s
    integer :: i
    !
    bn%bc_dir = bc_dir
    bn%bc = bck
    bn%strain = 0.0
    bn%stress = 0.0
    !
    do i = 1,2
      !
      select case (bn%bc(i))
      case(0)
        !Displacement prescribed:
        bn%disp(i) = prescribed_value(bcdf(i),s,bi(i,:))
        bn%trac(i) = 0.0
        bn%h1(i) = 0.0
        bn%h2(i) = 0.0
      case(1)
        !Traction prescribed:
        bn%trac(i) = prescribed_value(bcdf(i),s,bi(i,:))
        bn%disp(i) = 0.0
        bn%h1(i) = 0.0
        bn%h2(i) = 0.0
      case(2)
        !Elastic support [ti = h1(i)*ui + h2(i)]:
        bn%disp(i) = 0.0
        bn%trac(i) = 0.0
        bn%h1(i) = bi(1,1)
        bn%h2(i) = bi(1,2)
      case default
        call error_message_stop("Undefined boundary condition key.")
      end select
      !
    end do
    !
    return
  end subroutine set_boundary_condition_on_node
  !****************************************************
  !****************************************************
  function prescribed_value(bcdf,s,bi)
    !----------------------------------------------------
    ! Evaluate the prescribed boundary value for a given
    ! parametric coordinate s = [0,1].
    !
    ! Date: 27 June 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8) :: prescribed_value
    integer,intent(in) :: bcdf
    real(kind=8),intent(in) :: s
    real(kind=8),intent(in) :: bi(4)
    !
    if (s.lt.0.0 .or. s.gt.1.0) then
      call error_message_stop("Undefined parametric coordinate to interpolate boundary condition.")
    end if
    !
    prescribed_value = 0.0_8
    !
    select case (bcdf)
    case(0)
      prescribed_value = bi(1) + bi(2)*s + bi(3)*s*s + bi(4)*s*s*s
    case default
      call error_message_stop("Undefined boundary condition distribution function key.")
    end select
    !
    return
  end function prescribed_value
  !****************************************************
  !****************************************************
  subroutine set_particular_bc(ps,x,y,bc_dir,bck,pv)
    !----------------------------------------------------
    ! Set a particular thermal condition on a given node.
    !
    ! ps = 1 or 2, respectively for the 1st or 2nd
    !      appearance of the point in the boundary
    !      nodes list.
    !
    ! x,y = coordinates of the boundary node.
    !
    ! bc_dir: direction for prescribed mechanical fields
    !         = 0, (x,y) directions
    !         = 1, normal/tangential directions
    !
    ! bck(1) = 0 for displacement prescribed [ux or un],
    !          1 for traction prescribed [tx or tn].
    ! bck(2) = 0 for displacement prescribed [uy or ut],
    !          1 for traction prescribed [ty or tt].
    !
    ! pv(1) = presc. value (or h1 for elastic support) for x or normal direction
    ! pv(2) = presc. value (or h1 for elastic support) for y or tangential direction
    ! pv(3) = h2 for x or normal direction [ti = h1(i)*ui + h2(i)]
    ! pv(4) = h2 for y or tangential direction [ti = h1(i)*ui + h2(i)]
    !
    ! Date: 27 June 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: ps,bc_dir,bck(2)
    real(kind=8),intent(in) :: x,y,pv(4)
    type(point_2d) :: pt
    real(kind=8) :: dist
    real(kind=8) :: tol = 1.0e-6_8
    integer :: i,j,cps
    !
    cps = 0
    pt%x = x
    pt%y = y
    !
    do i = 1,tnbn
      !
      dist = distance_only_2d(lbn(i)%point,pt)
      !
      if (dist.lt.tol) then
        !
        cps = cps + 1
        if (cps.ne.ps) cycle
        !
        lbn(i)%bc_dir = bc_dir
        lbn(i)%bc = bck
        !
        do j = 1,2
          !
          select case (bck(j))
          case(0)
            !Displacement prescribed:
            lbn(i)%disp(j) = pv(j)
            lbn(i)%trac(j) = 0.0
            lbn(i)%h1(j) = 0.0
            lbn(i)%h2(j) = 0.0
          case(1)
            !Traction prescribed:
            lbn(i)%trac(j) = pv(j)
            lbn(i)%disp(j) = 0.0
            lbn(i)%h1(j) = 0.0
            lbn(i)%h2(j) = 0.0
          case(2)
            !Elastic support [ti = h1(i)*ui + h2(i)]:
            lbn(i)%trac(j) = 0.0
            lbn(i)%disp(j) = 0.0
            lbn(i)%h1(j) = pv(j)
            lbn(i)%h2(j) = pv(j+2)
          case default
            call error_message_stop("Undefined particular boundary condition key.")
          end select
          !
        end do
        !
      end if
      !
    end do
    !
    return
  end subroutine set_particular_bc
  !****************************************************
  !****************************************************
  subroutine int_pts_straight_line(x1,y1,x2,y2,nip,scl)
    !----------------------------------------------------
    ! Generate internal points along a straight line.
    !
    ! (x1,y1) = the initial point
    ! (x2,y2) = the final point
    ! nip = the number of points to be generated
    !
    ! Date: 27 June 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: nip
    real(kind=8),intent(in) :: x1,y1,x2,y2,scl
    type(point_2d) :: p1,p2
    type(point_2d),allocatable :: pts(:)
    integer :: i
    !
    p1%x = x1; p1%y = y1
    p2%x = x2; p2%y = y2
    call divide_straight_line_2d_scaling(p1,p2,nip-1,scl,pts)
    !
    do i = 1,nip
      tnip = tnip + 1
      lip(tnip)%point = pts(i)
      lip(tnip)%disp = 0.0
      lip(tnip)%strain = 0.0
      lip(tnip)%stress = 0.0
    end do
    !
    return
  end subroutine int_pts_straight_line
  !****************************************************
  !****************************************************
  subroutine int_pts_circular_line(x1,y1,x2,y2,ang,dk,nip,scl)
    !----------------------------------------------------
    ! Generate internal points along a circular line.
    !
    ! (x1,y1) = the circumference central point
    ! (x2,y2) = the line initial point
    ! ang = total angle for the line
    ! dk = direction key (negative for clockwise direction)
    ! nip = the number of points to be generated
    !
    ! Date: 27 June 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: nip,dk
    real(kind=8),intent(in) :: x1,y1,x2,y2,ang,scl
    type(point_2d) :: p1,p2
    type(point_2d),allocatable :: pts(:)
    real(kind=8) :: ang_rad,tol
    logical :: cc
    integer :: i,n_seg
    !
    real(kind=8) :: pi = 3.1415926535897932_8
    !
    if (ang.gt.360.) then
      call error_message_stop("Arc grater than 360 degrees for internal points generation.")
    end if
    !
    tol = 1.0e-8_8
    if (ang.gt.(360.0_8-tol)) then
      n_seg = nip
    else
      n_seg = nip-1
    end if
    !
    if (dk.lt.0) then
      cc = .false. !The line goes clockwise direction
    else
      cc = .true. !The line goes counter-clockwise direction
    end if
    !
    ang_rad = pi*ang/180.0_8
    !
    p1%x = x1; p1%y = y1
    p2%x = x2; p2%y = y2
    call divide_circular_line_2d_scaling(p1,p2,ang_rad,scl,n_seg,pts,cc)
    !
    do i = 1,nip
      tnip = tnip + 1
      lip(tnip)%point = pts(i)
      lip(tnip)%disp = 0.0
      lip(tnip)%strain = 0.0
      lip(tnip)%stress = 0.0
    end do
    !
    return
  end subroutine int_pts_circular_line
  !****************************************************
  !****************************************************
  subroutine int_pts_rectangular_area(x1,y1,x2,y2,niph,nipv,sclh,sclv)
    !----------------------------------------------------
    ! Generate internal points distributed over a
    ! rectangular area (aligned to x,y axis).
    !
    ! (x1,y1) = left-bottom vertex point
    ! (x2,y2) = right-upper vertex point
    ! niph = number of points to be generated on each horizontal line
    ! nipv = number of points to be generated on each vertical line
    !
    ! Date: 27 June 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: niph,nipv
    real(kind=8),intent(in) :: x1,y1,x2,y2,sclh,sclv
    real(kind=8) :: dy,dy0,y
    integer :: i
    !
    if (nipv.lt.2 .and. niph.lt.2) then
      call error_message_stop("At least two vertical/horizontal lines necessary to &
                              &generate a rectangular area of internal points.")
    end if
    !
    dy0 = 2.0_8*(y2-y1)/((nipv-1)*(sclv+1.0_8))
    dy = dy0
    y = y1 - dy
    !
    do i = 1,nipv
      y = y + dy
      call int_pts_straight_line(x1,y,x2,y,niph,sclh)
      if (i.eq.nipv) exit
      dy = dy0*(1.0_8 + (sclv-1.0_8)*(i-1)/(nipv-1))
    end do
    !
    return
  end subroutine int_pts_rectangular_area
  !****************************************************
  !****************************************************
  subroutine check_max_num_bound_elem_and_nodes(nbe,nbn)
    !----------------------------------------------------
    ! Check if the maximum numbers of boundary nodes or
    ! elements have been exceeded.
    !
    ! Date: 27 June 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: nbe,nbn
    integer :: new_bn
    !
    if (nbe+tnbe.gt.max_num_be) then
      call error_message_stop("Maximum number of boundary elements has been exceeded.")
    end if
    !
    if (nbn.eq.1) then
      new_bn = nbe
    else
      new_bn = nbe*(nbn-1) + 1
    end if
    !
    if (new_bn+tnbn.gt.max_num_bn) then
      call error_message_stop("Maximum number of boundary nodes has been exceeded.")
    end if
    !
    return
  end subroutine check_max_num_bound_elem_and_nodes
  !****************************************************
  !****************************************************
  subroutine check_connection_by_node_with_previous_line(x1,y1,lb1,flcb)
    !----------------------------------------------------
    ! Check if the current boundary line (which begins
    ! on point x1,y1) has the last node of previous line
    ! as its first node. If true, pertinent adjustments
    ! are performed.
    !
    ! Date: 27 June 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: x1,y1,lb1
    logical,intent(in) :: flcb
    type(point_2d) :: p1
    real(kind=8) :: dist
    real(kind=8) :: tol = 1.0e-8_8
    !
    if (lb1.lt.1.0_8) return
    if (flcb) return
    !
    p1%x = x1; p1%y = y1
    dist = distance_only_2d(lbn(tnbn)%point,p1)
    !
    if (dist.lt.tol) then
      tnbn = tnbn - 1
    end if
    !
    return
  end subroutine check_connection_by_node_with_previous_line
  !****************************************************
  !****************************************************
  subroutine check_connection_by_node_for_closed_bound(first_pt,first_pt_pos,flcb,cip,nbe,nbn)
    !----------------------------------------------------
    ! Check if the last node of the current line is
    ! also the first node of the current closed boundary.
    ! If true, pertinent adjustments are performed.
    !
    ! If the current line is the first line of the
    ! boundary, the first point is stored.
    !
    ! Date: 27 June 2022
    ! Last update: ---
    !----------------------------------------------------
    type(point_2d),intent(inout) :: first_pt
    integer,intent(inout) :: first_pt_pos
    logical,intent(inout) :: flcb,cip
    integer,intent(in) :: nbe,nbn
    real(kind=8) :: dist
    real(kind=8) :: tol = 1.0e-8_8
    !
    if (flcb) then
      first_pt_pos = tnbn-(nbe*(nbn-1))
      first_pt = lbn(first_pt_pos)%point
      flcb = .false.
    else
      dist = distance_only_2d(lbn(tnbn)%point,first_pt)
      if (dist.lt.tol .and. cip) then
        tnbn = tnbn - 1
        lbe(tnbe)%incid(nbn) = first_pt_pos
        flcb = .true.
        cip = .false.
      end if
    end if
    !
    return
  end subroutine check_connection_by_node_for_closed_bound
  !****************************************************
  !****************************************************
  subroutine check_max_num_int_points(nip)
    !----------------------------------------------------
    ! Check if the maximum number of internal points has
    ! been exceeded.
    !
    ! Date: 27 June 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: nip
    !
    if (nip+tnip.gt.max_num_in) then
      call error_message_stop("Maximum number of internal points has been exceeded.")
    end if
    !
    return
  end subroutine check_max_num_int_points
  !****************************************************
  !****************************************************
end module fracture_2d_input
