module fracture_2d_crack_propagation
!----------------------------------------------------
! A module with subroutines that propagates the crack.
!
!
! Date: October 2022
! Last update: ---
!----------------------------------------------------
use fracture_2d_model
use matrix_algebra
use fracture_2d_j_integral
use fracture_2d_assembler
use fracture_2d_set_bound_results
use fracture_2d_stress_intens_factor
!
implicit none
!
type(boundary_elem),allocatable,private :: lbe(:)
type(boundary_node),allocatable,private :: lbn(:)
integer,private :: tnbe,tnbn,k_eff_key
real(kind=8),private :: R,C0,m,gama,K_th,K_IC,Da_0,Da_max,r_Da,l,length_0
real(kind=8),private :: N_step(1000,2) = 0.0_8,N_total(1000,2) = 0.0_8,delta_a_i(1000,2) = 0.0_8,y(1000,2,3) = 0.0_8 !Aloca com mais memória, pois não é possível acessar quantas trincas aqui
!
contains
  subroutine crack_propagation(n)
    !----------------------------------------------------
    ! Subroutine that propagates the crack.
    !
    ! Date: 28 October 2022
    ! Last update: ---
    !---------------------------------------------------
    integer,intent(in) :: n
    real(kind=8) :: theta
    real(kind=8),allocatable :: beta_i(:),beta_i1(:),theta_p(:),coord(:,:)
    real(kind=8) :: aux_x1,aux_x2
    integer :: i,j,k,aux_be,aux_bn,b,c,d,tips_num,n_el
    type(boundary_elem) :: last_elem
    real(kind=8),allocatable :: a(:,:),f(:)
    logical :: bound_recovery,iteration
    logical,allocatable :: propagation_key(:)
    real(kind=8),allocatable :: h(:,:),g(:,:)
    real(kind=8),allocatable :: delta_k_eff_1(:,:), delta_k_eff_2(:,:),delta_a(:,:)
    !
    !-------------------------------------------------------------
    ! Saving the information of the previous step:
    !-------------------------------------------------------------
    !
    tnbe = be_num
    tnbn = bn_num
    allocate(lbe(tnbe))
    allocate(lbn(tnbn))
    lbe = be_list
    lbn = bn_list
    !
    k = 1
    aux_x1 = 0.0
    aux_x2 = 0.0
    !
    call get_tips_num(tips_num)
    !
    allocate(beta_i1(tips_num))
    allocate(beta_i(tips_num))
    allocate(theta_p(tips_num))
    allocate(coord(tips_num,4))
    allocate(propagation_key(tips_num))
    allocate(delta_k_eff_1(cr_num,2))
    allocate(delta_k_eff_2(cr_num,2))
    allocate(delta_a(cr_num,2))
    beta_i1 = 0._8
    beta_i = 0._8
    !
    call propagation_angle()
    !
    !----------------------------------------------------------------------
    ! Saving the information of the initial crack (before adding delta a):
    !----------------------------------------------------------------------
    !
    b = 1
    do i = 1,cr_num
        do j = 1,cr_list(i)%tips_num
            if (cr_list(i)%keq(j) .lt. K_th) then
                call error_message_2int_stop("Crack ",i, "is not propagating on tip ",j) !Melhorar
            end if
            last_elem = be_list(cr_list(i)%ext_el_num(2*j))
            theta = atan2((last_elem%bp_2%y-last_elem%bp_1%y),(last_elem%bp_2%x-last_elem%bp_1%x))
            theta_p(b) = cr_list(i)%theta(j) + theta
            if (n .eq. 1) then
                l = distance_only_2d(last_elem%bp_1,last_elem%bp_2)
            end if
            coord(b,1) = last_elem%bp_2%x
            coord(b,2) = last_elem%bp_2%y
            !
            !-----------------------------------------
            ! Calculating delta_k_eff(a0) and y_1:
            !-----------------------------------------
            !
            call get_delta_k_eff(i,j,delta_k_eff_1(i,j))
            !
            if ((mod(n,2) == 1)) then
                y(i,j,1) = 1/(C0*(delta_k_eff_1(i,j))**m)
                if (j .eq. 1) then
                    length_0 = cr_list(i)%length/cr_list(i)%tips_num
                end if
            end if
            !
            !-----------------------------------------
            ! Getting delta_a and number of elements:
            !-----------------------------------------
            !
            call get_delta_a_n_el(i,j,l,n,delta_a(i,j),n_el)
            !
            cr_list(i)%length = cr_list(i)%length + delta_a(i,j)
            if (j .eq. 2) then
                aux_x1 = coord(b,1)
                aux_x2 = coord(b,2)
            end if
            b = b + 1
        end do
    end do
    !
    !--------------------------
    ! Loop on the active tips:
    !--------------------------
    !
    iteration = .true.
    propagation_key = .true.
    do while (iteration)
        beta_i = beta_i1
        d = 0
        b = 1
        do i = 1,cr_num
            do j = 1,cr_list(i)%tips_num
                d = d + 1
                if (propagation_key(d) .eqv. .true.) then
                    theta_p(d) = theta_p(d)+beta_i(d)
                    if (j.eq.1) then
                        coord(d,3) = coord(d,1)+delta_a(i,j)*cos(theta_p(d))
                        coord(d,4) = coord(d,2)+delta_a(i,j)*sin(theta_p(d))
                        if (k.eq.1) then
                            cr_list(i)%ext_el_num(2) = tnbe+b*n_el
                            cr_list(i)%ext_el_num(3) = tnbe+b*n_el+1
                            b = b + 1
                        end if
                        else if (j.eq.2) then
                            coord(d,3) = aux_x1
                            coord(d,4) = aux_x2
                            coord(d,1) = aux_x1+delta_a(i,j)*cos(theta_p(d))
                            coord(d,2) = aux_x2+delta_a(i,j)*sin(theta_p(d))
                            if(k.eq.1) then
                                cr_list(i)%ext_el_num(1) = tnbe+b*n_el+1
                                cr_list(i)%ext_el_num(4) = tnbe+2*b*n_el
                                b = b + 1
                            end if
                    end if
                end if
            end do
            b = b + 1
        end do
        !
        !
        !----------------------------------------------
        ! Creating the lists with the new information:
        !----------------------------------------------
        !
        aux_be = tnbe+2*d*n_el
        aux_bn = tnbn+2*d*n_el*last_elem%nodes_num
        deallocate(be_list)
        deallocate(bn_list)
        allocate(be_list(aux_be))
        allocate(bn_list(aux_bn))
        be_list(1:tnbe) = lbe
        bn_list(1:tnbn) = lbn
        be_num = tnbe
        bn_num = tnbn
        do c = 1,d
            call straight_line_for_propagation(coord(c,1),coord(c,2),coord(c,3),coord(c,4),n_el,&
                & last_elem%nodes_num,last_elem%lamb_1,last_elem%lamb_2)
        end do
        !
        !-----------------------------------------------
        ! Creating a new list with the symmetric nodes:
        !-----------------------------------------------
        !
        call symetric_nodes(bn_num,tnbn,d)
        !
        !-------------------------------------------------------
        ! Solving the fracture problem for the new information:
        !-------------------------------------------------------
        !
        allocate(a(2*bn_num,2*bn_num))
        allocate(f(2*bn_num))
        if(allocated(h)) deallocate(h)
        if(allocated(g)) deallocate(g)
        allocate(h(2*bn_num,2*bn_num))
        allocate(g(2*bn_num,2*bn_num))
        call assembler_for_propagation(a,f,h,g)
        call solve_linear_system(a,f)
        deallocate(a)
        call set_results_on_boundary_nodes(f)
        deallocate(f)
        !
        bound_recovery = .false.
        call evaluate_stress_strain_on_boundary_nodes(bound_recovery)
        call stress_intensity_factor_evaluation()
        call propagation_angle()
        call equivalent_k()
        !
        !-----------------------------------
        ! Calculating the correction angle:
        !-----------------------------------
        !
        !
        d = 0
        do i = 1,cr_num
            do j = 1,cr_list(i)%tips_num
                d = d + 1
                beta_i1(d) = cr_list(i)%theta(j)/2.0_8
            end do
        end do
        !
        !------------------------------------
        ! Checking if |Beta i+1| < |Beta i|:
        !------------------------------------
        !
        if (k .gt. 1) then
            do i = 1,tips_num
                if (abs(beta_i1(i)) .lt. abs(beta_i(i)) .and. propagation_key(i) .eqv. .true.) then
                    propagation_key(i) = .true.
                    else
                        propagation_key(i) = .false.
                end if
            end do
        end if
        !
        !-------------------------------------------------
        ! Checking if the tip will continue to propagate:
        !-------------------------------------------------
        !
        iteration = .false.
        do i = 1,tips_num
            if(propagation_key(i) .eqv. .true.) then
                iteration = .true.
            end if
        end do
        k = k + 1
    end do
    !
    !--------------------------------------------------
    ! Calculating delta_Keff(a) (after the increment):
    !--------------------------------------------------
    !
    do i = 1,cr_num
        do j = 1,cr_list(i)%tips_num
            !
            call get_delta_k_eff(i,j,delta_k_eff_2(i,j))
            !
            if ((mod(n,2) == 1)) then
                y(i,j,2) = 1/(C0*(delta_k_eff_2(i,j))**m)
                else if ((mod(n,2) == 0)) then
                    y(i,j,3) = 1/(C0*(delta_k_eff_2(i,j)**m))
            end if
        end do
    end do
    !
    !
    call set_initial_h(h)
    call set_initial_g(g)
    deallocate(lbe)
    deallocate(lbn)
    deallocate(beta_i1)
    deallocate(beta_i)
    deallocate(theta_p)
    deallocate(coord)
    !
    !
    return
  end subroutine crack_propagation
  !****************************************************
  !****************************************************
  subroutine number_of_cycles()
    !----------------------------------------------------
    ! Subroutine that propagates the crack.
    !
    ! Date: 07 December 2022
    ! Last update: ---
    !----------------------------------------------------
    integer :: i,j
    !
    !
    do i = 1,cr_num
        do j = 1,cr_list(i)%tips_num
          N_step(i,j) = delta_a_i(i,j)*(y(i,j,1) + 4*y(i,j,2) + y(i,j,3))/3
          N_total(i,j) = N_total(i,j) + N_step(i,j)
        end do
    end do
  end subroutine number_of_cycles
  !****************************************************
  !****************************************************
  subroutine get_tips_num(tips_num)
    !----------------------------------------------------
    ! Subroutine that calculates the number of tips.
    !
    ! Date: 20 December 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(out) :: tips_num
    integer :: i,j
    !
    tips_num = 0
    do i = 1,cr_num
        do j = 1,cr_list(i)%tips_num
            tips_num = tips_num+1
        end do
    end do
    !
  end subroutine get_tips_num
  !****************************************************
  !****************************************************
  subroutine get_delta_k_eff(i,j,delta_k_eff)
    !----------------------------------------------------
    ! Subroutine that calculates the number of tips.
    !
    ! Date: 20 December 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: i,j
    real(kind=8),intent(out) :: delta_k_eff
    real(kind=8) :: ang,delta_ki,delta_kii
    !
    delta_ki = cr_list(i)%ki(j)
    delta_kii = cr_list(i)%kii(j)
    if (k_eff_key .eq. 1) then
        delta_k_eff = (delta_ki**4.0_8+8.0_8*delta_kii**4.0_8)**(0.25_8)
        else if (k_eff_key .eq. 2) then
            ang = cr_list(i)%theta(j)
            delta_k_eff = delta_ki*(cos(ang/2.0_8))**3-3._8*delta_kii*sin(ang/2.0_8)*(cos(ang/2.0_8))**2
    end if
    !
  end subroutine get_delta_k_eff
  !****************************************************
  !****************************************************
  subroutine get_delta_a_n_el(i,j,l,n,delta_a,n_el)
    !----------------------------------------------------
    ! Subroutine that calculates delta a and the number.
    ! of elements.
    !
    ! Date: 20 December 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: i,j,n
    real(kind=8),intent(in) :: l
    real(kind=8),intent(out) :: delta_a
    integer,intent(out) :: n_el
    real(kind=8) :: aux_n_el
    !
    if (n .eq. 1 .or. n .eq. 2) then
        delta_a = Da_0
        delta_a_i(i,j) = delta_a
        else
            if ((mod(n,2) == 1)) then
                delta_a = r_Da * delta_a_i(i,j)
                else
                    delta_a = delta_a_i(i,j)
            end if
            if (delta_a .gt. Da_max) then
                delta_a = Da_max
            end if
            delta_a_i(i,j) = delta_a
    end if
    aux_n_el = delta_a/l
    if (aux_n_el .lt. 1.0_8) then
        n_el = 1
        else
            n_el = nint(aux_n_el)
    end if
    !
  end subroutine get_delta_a_n_el
  !****************************************************
  !****************************************************
  subroutine symetric_nodes(bn_num,tnbn,d)
    !----------------------------------------------------
    ! Subroutine that creates the new list of symetric
    ! nodes.
    !
    ! Date: 20 December 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: d,bn_num,tnbn
    integer :: c,i,j
    integer :: aux1,aux2,aux3,aux4,aux5,aux6,aux7,aux8
    integer,allocatable :: nodes_num_temp(:,:)
    !
    c = 0
    do i = 1,cr_num
        aux1 = (bn_num - tnbn)/(2*d)
        aux2 = size(cr_list(i)%nodes_num(:,1))
        aux3 = aux1+aux2
        aux4 = tnbn+2*aux1+c*aux1+1
        aux5 = tnbn+4*aux1+c*aux1
        aux6 = cr_list(i)%tips_num*aux1+aux2
        aux7 = tnbn+c*aux1+1
        aux8 = tnbn+2*aux1+c*aux1
        allocate(nodes_num_temp(aux2,2))
        nodes_num_temp = cr_list(i)%nodes_num
        deallocate(cr_list(i)%nodes_num)
        allocate(cr_list(i)%nodes_num(aux6,2))
        if (cr_list(i)%tips_num .eq. 2) then
            do j = 1,aux1
                cr_list(i)%nodes_num(j,1) = aux4
                cr_list(i)%nodes_num(j,2) = aux5
                aux4 = aux4+1
                aux5 = aux5-1
            end do
            cr_list(i)%nodes_num(aux1+1:aux3,:) = nodes_num_temp
            deallocate(nodes_num_temp)
            do j = aux3+1,aux6
                cr_list(i)%nodes_num(j,1) = aux7
                cr_list(i)%nodes_num(j,2) = aux8
                aux7 = aux7+1
                aux8 = aux8-1
            end do
            else if (cr_list(i)%tips_num .eq. 1) then
                cr_list(i)%nodes_num(1:aux2,:) = nodes_num_temp
                deallocate(nodes_num_temp)
                do j = aux2+1,aux6
                    cr_list(i)%nodes_num(j,1) = aux7
                    cr_list(i)%nodes_num(j,2) = aux8
                    aux7 = aux7+1
                    aux8 = aux8-1
                end do
        end if
        c = c + 2*cr_list(i)%tips_num
    end do
    !
  end subroutine symetric_nodes
  !****************************************************
  !****************************************************
  subroutine propagation_angle()
    !----------------------------------------------------
    ! Calculates the angle of crack propagation direction.
    !
    ! Date: 20 October 2022
    ! Last update: ---
    !---------------------------------------------------
    integer :: i,j
    real(kind=8) :: theta_p
    real(kind=8) :: k_i,k_ii
    !
    do i = 1,cr_num
        do j = 1,cr_list(i)%tips_num
            k_i = cr_list(i)%ki(j)
            k_ii = cr_list(i)%kii(j)
            if (k_ii .gt. 0.0_8) then
                theta_p = 2*atan(0.25_8*((k_i/k_ii)-sqrt((k_i/k_ii)**2+8.0_8)))
                else
                    theta_p = 2*atan(0.25_8*((k_i/k_ii)+sqrt((k_i/k_ii)**2+8.0_8)))
            end if
            cr_list(i)%theta(j) = theta_p
        end do
    end do
    !
    return
  end subroutine propagation_angle
  !****************************************************
  !****************************************************
  subroutine equivalent_k()
    !----------------------------------------------------
    ! Calculates the equivalent k, to determine if the
    ! crack is going to propagate.
    !
    ! Date: 20 October 2022
    ! Last update: ---
    !---------------------------------------------------
    real(kind=8) :: k_i,k_ii,theta_p
    integer :: i,j
    !
    do i = 1,cr_num
        do j = 1,cr_list(i)%tips_num
            k_i = cr_list(i)%ki(j)
            k_ii = cr_list(i)%kii(j)
            theta_p = cr_list(i)%theta(j)
     !Heider!cr_list(i)%keq(j) = cos(theta_p/2.0_8)*(k_i*(cos(theta_p/2.0_8))**2-1.5_8*k_ii*sin(theta_p/2.0_8))
            cr_list(i)%keq(j) = k_i*(cos(theta_p/2.0_8))**3-3._8*k_ii*sin(theta_p/2.0_8)*(cos(theta_p/2.0_8))**2
            !Portela
        end do
    end do
    !
    return
  end subroutine equivalent_k
  !****************************************************
  !****************************************************
  subroutine straight_line_for_propagation(x1,y1,x2,y2,nbe,nbn,lb1,lb2)
    !----------------------------------------------------
    ! Generate elements and boundary nodes for the crack
    ! propagation, which occurs in a straight line.
    !
    !
    ! Date: 28 October 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: nbe,nbn
    real(kind=8),intent(in) :: x1,y1,x2,y2,lb1,lb2
    integer :: bck(2),bcdf(2),bc_dir
    real(kind=8) :: bi(2,4)
    type(point_2d) :: p1,p2
    type(point_2d),allocatable :: pts(:),pts2(:)
    integer :: j,n,bie,ck
    !
    !
    !-------------------------------------------------------------
    ! Evaluating a list of points to generate elements and nodes:
    !-------------------------------------------------------------
    !
    p1%x = x1; p1%y = y1
    p2%x = x2; p2%y = y2
    !
    bie = 2
    ck = 2
    !
    bc_dir = 0
    bck = 1
    bcdf = 0
    bi = 0._8
    !
    if (lb1.eq.1.0_8 .or. lb2.eq.1.0_8) then
        call error_message_stop("Invalid lambda 1 or lambda 2, generating a continuous boundary elements.")
    end if
    !
    if (nbn.lt.3) then
      call divide_straight_line_2d_scaling(p1,p2,nbe,1.0_8,pts)
    else
      call divide_straight_line_2d_semi_scaling(p1,p2,nbe,1.0_8,nbn-1,pts)
    end if
    !
    !
    !--------------------------------
    ! Generating elements and nodes:
    !--------------------------------
    !
    select case (nbn)
    case(1)
      call constant_elements_for_propagation(pts,bc_dir,bck,bcdf,bi,bie)
    case(2)
      call linear_elements_for_propagation(pts,lb1,lb2,bc_dir,bck,bcdf,bi,bie,ck)
    case(3)
      call quadratic_elements_for_propagation(pts,lb1,lb2,bc_dir,bck,bcdf,bi,bie,ck)
    case default
      call error_message_stop("Undefined element order")
    end select
    !
    !
    !-------------------------------------------
    ! Generating elements and nodes for cracks:
    !-------------------------------------------
    !
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
            call constant_elements_for_propagation(pts2,bc_dir,bck,bcdf,bi,bie)
        case(2)
            call linear_elements_for_propagation(pts2,lb2,lb1,bc_dir,bck,bcdf,bi,bie,ck)
        case(3)
            call quadratic_elements_for_propagation(pts2,lb2,lb1,bc_dir,bck,bcdf,bi,bie,ck)
        case default
            call error_message_stop("Undefined element order.")
    end select
    !
    return
  end subroutine straight_line_for_propagation
  !****************************************************
  !****************************************************
  subroutine constant_elements_for_propagation(pts,bc_dir,bck,bcdf,bi,bie)
    !----------------------------------------------------
    ! Generate constant elements for a given list of
    ! points.
    !
    ! Date: 28 October 2022
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
    call get_parametric_coords_for_propagation(pts,1.0_8,1.0_8,1,s)
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
      be_num = be_num + 1 !Adding a boundary element to the model
      !
      be_list(be_num)%bp_1 = p1
      be_list(be_num)%bp_2 = p2
      allocate(be_list(be_num)%incid(1))
      be_list(be_num)%incid(1) = be_num
      be_list(be_num)%lamb_1 = 1.0_8
      be_list(be_num)%lamb_2 = 1.0_8
      be_list(be_num)%nodes_num = 1
      !
      bn_num = bn_num + 1 !Adding a boundary node to the model
      !
      bn_list(bn_num)%point = midpoint_2d(p1,p2)
      bn_list(bn_num)%bie = bie
      call set_boundary_condition_on_node_for_propagation(bn_list(bn_num),bc_dir,bck,bcdf,s(i),bi)
      !
    end do
    !
    return
  end subroutine constant_elements_for_propagation
  !****************************************************
  !****************************************************
  subroutine linear_elements_for_propagation(pts,lb1,lb2,bc_dir,bck,bcdf,bi,bie,ck)
    !----------------------------------------------------
    ! Generate linear elements for a given list of
    ! points.
    !
    ! Date: 28 October 2022
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
    call get_parametric_coords_for_propagation(pts,lb1,lb2,2,s)
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
      be_num = be_num + 1 !Adding a boundary element to the model
      !
      be_list(be_num)%bp_1 = p1
      be_list(be_num)%bp_2 = p2
      be_list(be_num)%nodes_num = 2
      allocate(be_list(be_num)%incid(2))
      if (ck.eq.1) then
        if (i.eq.1) then
            be_list(be_num)%lamb_1 = lb1
        else
            be_list(be_num)%lamb_1 = 1.0_8
        end if
        if (i.eq.size(pts)-1) then
            be_list(be_num)%lamb_2 = lb2
        else
            be_list(be_num)%lamb_2 = 1.0_8
        end if
      end if
      if (ck.eq.2) then
        be_list(be_num)%lamb_1 = lb1
        be_list(be_num)%lamb_2 = lb2
      end if
      !
      if (ck.eq.2) then
        if (i.eq.1) then
            bn_num = bn_num + 1 !Adding a boundary node to the model
            el_pts(1) = p1
            el_pts(2) = p2
            bn_list(bn_num)%point = cartesian_point_from_intrp_func_1d(el_pts,-lb1)
            bn_list(bn_num)%bie = bie
            call set_boundary_condition_on_node_for_propagation(bn_list(bn_num),bc_dir,bck,bcdf,s(i),bi)
        end if
        else if (ck.eq.2) then
            bn_num = bn_num + 1 !Adding a boundary node to the model
            el_pts(1) = p1
            el_pts(2) = p2
            bn_list(bn_num)%point = cartesian_point_from_intrp_func_1d(el_pts,-lb1)
            bn_list(bn_num)%bie = bie
            call set_boundary_condition_on_node_for_propagation(bn_list(bn_num),bc_dir,bck,bcdf,s(i),bi)
      end if
      !
      be_list(be_num)%incid(1) = bn_num
      !
      bn_num = bn_num + 1 !Adding a boundary node to the model
      if (ck.eq.1) then
        if (i.eq.size(pts)-1 .and. ck.eq.1) then
            el_pts(1) = p1
            el_pts(2) = p2
            bn_list(bn_num)%point = cartesian_point_from_intrp_func_1d(el_pts,+lb2)
        else
            bn_list(bn_num)%point = p2
        end if
      else if (ck.eq.2) then
            el_pts(1) = p1
            el_pts(2) = p2
            bn_list(bn_num)%point = cartesian_point_from_intrp_func_1d(el_pts,+lb2)
      end if
      be_list(be_num)%incid(2) = bn_num
      bn_list(bn_num)%bie = bie
      call set_boundary_condition_on_node_for_propagation(bn_list(bn_num),bc_dir,bck,bcdf,s(i+1),bi)
      !
    end do
    !
    return
  end subroutine linear_elements_for_propagation
  !****************************************************
  !****************************************************
  subroutine quadratic_elements_for_propagation(pts,lb1,lb2,bc_dir,bck,bcdf,bi,bie,ck)
    !----------------------------------------------------
    ! Generate quadratic elements for a given list of
    ! points.
    !
    ! Date: 28 October 2022
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
    call get_parametric_coords_for_propagation(pts,lb1,lb2,3,s)
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
      be_num = be_num + 1 !Adding a boundary element to the model
      !
      be_list(be_num)%bp_1 = p1
      be_list(be_num)%bp_2 = p3
      be_list(be_num)%nodes_num = 3
      allocate(be_list(be_num)%incid(3))
      if (ck.eq.1) then
        if (i.eq.1) then
            be_list(be_num)%lamb_1 = lb1
        else
            be_list(be_num)%lamb_1 = 1.0_8
        end if
        if (i.eq.n) then
            be_list(be_num)%lamb_2 = lb2
        else
            be_list(be_num)%lamb_2 = 1.0_8
        end if
      end if
      if (ck.eq.2) then
        be_list(be_num)%lamb_1 = lb1
        be_list(be_num)%lamb_2 = lb2
      end if
      !
      if (i.eq.1 .and. ck.eq.1) then
        bn_num = bn_num + 1 !Adding a boundary node to the model
        el_pts(1) = p1
        el_pts(2) = p2
        el_pts(3) = p3
        bn_list(bn_num)%point = cartesian_point_from_intrp_func_1d(el_pts,-lb1)
        bn_list(bn_num)%bie = bie
        call set_boundary_condition_on_node_for_propagation(bn_list(bn_num),bc_dir,bck,bcdf,s(i),bi)
        else if (ck.eq.2) then
            bn_num = bn_num + 1 !Adding a boundary node to the model
            el_pts(1) = p1
            el_pts(2) = p2
            el_pts(3) = p3
            bn_list(bn_num)%point = cartesian_point_from_intrp_func_1d(el_pts,-lb1)
            bn_list(bn_num)%bie = bie
            call set_boundary_condition_on_node_for_propagation(bn_list(bn_num),bc_dir,bck,bcdf,s(i),bi)
      end if
      !
      be_list(be_num)%incid(1) = bn_num
      !
      bn_num = bn_num + 1 !Adding a boundary node to the model
      be_list(be_num)%incid(2) = bn_num
      bn_list(bn_num)%point = p2
      bn_list(bn_num)%bie = bie
      call set_boundary_condition_on_node_for_propagation(bn_list(bn_num),bc_dir,bck,bcdf,s(2*i),bi)
      !
      bn_num = bn_num + 1 !Adding a boundary node to the model
      if (ck.eq.1) then
        if (i.eq.n) then
            el_pts(1) = p1
            el_pts(2) = p2
            el_pts(3) = p3
            bn_list(bn_num)%point = cartesian_point_from_intrp_func_1d(el_pts,+lb2)
        else
            bn_list(bn_num)%point = p3
        end if
        be_list(be_num)%incid(3) = bn_num
        bn_list(bn_num)%bie = bie
        call set_boundary_condition_on_node_for_propagation(bn_list(bn_num),bc_dir,bck,bcdf,s(2*i+1),bi)
      else if (ck.eq.2) then
        bn_list(bn_num)%point = cartesian_point_from_intrp_func_1d(el_pts,+lb2)
        be_list(be_num)%incid(3) = bn_num
        bn_list(bn_num)%bie = bie
        call set_boundary_condition_on_node_for_propagation(bn_list(bn_num),bc_dir,bck,bcdf,s(2*i+1),bi)
      end if
      !
    end do
    !
    return
  end subroutine quadratic_elements_for_propagation
  !****************************************************
  !****************************************************
  subroutine get_parametric_coords_for_propagation(pts,lb1,lb2,nbn,s)
    !----------------------------------------------------
    ! Return a vector with approximate parametric
    ! coordinates, s = [0,1], for a given set of points.
    !
    ! Date: 28 October 2022
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
  end subroutine get_parametric_coords_for_propagation
  !****************************************************
  !****************************************************
  subroutine set_boundary_condition_on_node_for_propagation(bn,bc_dir,bck,bcdf,s,bi)
    !----------------------------------------------------
    ! Set boundary condition parameters on a given node.
    !
    ! Date: 28 October 2022
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
        bn%disp(i) = prescribed_value_for_propagation(bcdf(i),s,bi(i,:))
        bn%trac(i) = 0.0
        bn%h1(i) = 0.0
        bn%h2(i) = 0.0
      case(1)
        !Traction prescribed:
        bn%trac(i) = prescribed_value_for_propagation(bcdf(i),s,bi(i,:))
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
  end subroutine set_boundary_condition_on_node_for_propagation
  !****************************************************
  !****************************************************
  function prescribed_value_for_propagation(bcdf,s,bi)
    !----------------------------------------------------
    ! Evaluate the prescribed boundary value for a given
    ! parametric coordinate s = [0,1].
    !
    ! Date: 28 October 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8) :: prescribed_value_for_propagation
    integer,intent(in) :: bcdf
    real(kind=8),intent(in) :: s
    real(kind=8),intent(in) :: bi(4)
    !
    if (s.lt.0.0 .or. s.gt.1.0) then
      call error_message_stop("Undefined parametric coordinate to interpolate boundary condition.")
    end if
    !
    prescribed_value_for_propagation = 0.0_8 !MUDEI ISSO PORQUE ESTAVA DANDO WARNING QUE ESSE VALOR NÃO TINHA SIDO INICIALIZADO
    !
    select case (bcdf)
    case(0)
      prescribed_value_for_propagation = bi(1) + bi(2)*s + bi(3)*s*s + bi(4)*s*s*s
    case default
      call error_message_stop("Undefined boundary condition distribution function key.")
    end select
    !
    return
  end function prescribed_value_for_propagation
  !****************************************************
  !****************************************************
  subroutine assembler_for_propagation(a,f,h_initial,g_initial)
    !----------------------------------------------------
    ! Firstly mounts matrices [H] and [G] and, then,
    ! reorganizes according the boundary conditions to
    ! generate [A] and {f}.
    !
    ! Corners and traction discontinuities are treated by
    ! discontinuous elements or double nodes.
    !
    ! Date: 02 November 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(out) :: a(2*bn_num,2*bn_num)
    real(kind=8),intent(out) :: f(2*bn_num)
    real(kind=8),intent(out) :: h_initial(2*bn_num,2*bn_num),g_initial(2*bn_num,2*bn_num)
    real(kind=8) :: h(2*bn_num,2*bn_num),g(2*bn_num,2*bn_num),f_vec(2*bn_num)
    real(kind=8) :: aux(2*bn_num)
    real(kind=8) :: nv(2)
    type(boundary_node) :: cp
    integer :: j,k
    !
    !
    call get_matrices_H_and_G_for_propagation(h,g)
    !
    h_initial = h
    g_initial = g
    !
    !----------------------------------
    ! Applying the boundary conditions:
    !----------------------------------
    !
    do j = 1,bn_num
      !
      cp = bn_list(j)
      !
      ! Modifying matrices for normal/tangential b.c.:
      !
      if (cp%bc_dir.ne.0) then
        !
        nv = get_normal_to_a_boundary_point(j)
        k = 2*j-1
        !
        aux = h(:,k)*nv(1) + h(:,k+1)*nv(2)
        h(:,k+1) = h(:,k+1)*nv(1) - h(:,k)*nv(2)
        h(:,k) = aux
        !
        aux = g(:,k)*nv(1) + g(:,k+1)*nv(2)
        g(:,k+1) = g(:,k+1)*nv(1) - g(:,k)*nv(2)
        g(:,k) = aux
        !
      end if
      !
      ! Reorganizing matrices H and G according the b.c.:
      !
      do k = 1,2
        !
        select case (cp%bc(k))
        case(0)
          !Displacement component prescribed:
          f_vec(2*j+k-2) = cp%disp(k)
          aux = g(:,2*j+k-2)
          g(:,2*j+k-2) = -h(:,2*j+k-2)
          h(:,2*j+k-2) = -aux
        case(1)
          !Traction component prescribed:
          f_vec(2*j+k-2) = cp%trac(k)
        case(2)
          !Elastic support (Robin b.c.):
          f_vec(2*j+k-2) = cp%h2(k)
          h(:,2*j+k-2) = h(:,2*j+k-2) - g(:,2*j+k-2)*cp%h1(k)
        case default
          !Invalid boundary condition
          call error_message_int_stop("Undefined boundary condition at node: ",j)
        end select
        !
      end do
      !
    end do
    !
    !----------------------------------
    ! Mounting matrices [A] and {f}:
    !----------------------------------
    !
    f = matmul(g,f_vec)
    a = h
    !
    return
  end subroutine assembler_for_propagation
  !****************************************************
  !****************************************************
  subroutine get_matrices_H_and_G_for_propagation(h,g)
    !----------------------------------------------------
    ! Verify if the problem will be solved by displacement
    ! or traction BIE and chose the appropriate code.
    !
    ! Date: 12 August 2022
    ! Last update: ---
    !---------------------------------------------------
    real(kind=8),intent(inout) :: h(2*bn_num,2*bn_num)
    real(kind=8),intent(inout) :: g(2*bn_num,2*bn_num)
    real(kind=8),allocatable :: initial_h(:,:),initial_g(:,:)
    real(kind=8),allocatable :: gi(:,:),hi(:,:)
    type(boundary_node) :: cp
    type(boundary_elem) :: be
    real(kind=8) :: normal_cp(2)
    integer :: b,c,d,i,j,k,l,ii,kk,ll,bie
    !
    allocate(initial_h(2*tnbn,2*tnbn))
    allocate(initial_g(2*tnbn,2*tnbn))
    call get_initial_h(initial_h)
    call get_initial_g(initial_g)
    !
    h = 0.0_8
    g = 0.0_8
    h(1:2*tnbn,1:2*tnbn) = initial_h(1:2*tnbn,1:2*tnbn)
    g(1:2*tnbn,1:2*tnbn) = initial_g(1:2*tnbn,1:2*tnbn)
    deallocate(initial_h)
    deallocate(initial_g)
    !
    !---------------------------------------------------
    !Loop of new nodes on all the boundary elements:
    !---------------------------------------------------
    !
    do i = tnbn+1,bn_num
      !
      cp = bn_list(i)
      bie = cp%bie
      ii = 2*i-1
      normal_cp = get_normal_to_a_boundary_point(i)
      !
      !----------------------------
      ! Loop on boundary elements:
      !----------------------------
      !
      do j = 1,be_num
        !
        be = be_list(j)
        !
        allocate(hi(2,2*be%nodes_num))
        allocate(gi(2,2*be%nodes_num))
        !
        hi = 0.0_8
        gi = 0.0_8
        if (bie.eq.2) then
            call displ_bie_matrices_H_and_G(cp,be,hi,gi)
            else if (bie.eq.3) then
                call tract_bie_matrices_H_and_G(cp,be,normal_cp,hi,gi)
            else
                call error_message_stop("Incorrect Boundary Integral Equation for node.")
        end if
        !
        do k = 1,be%nodes_num
          kk = 2*be%incid(k)-1
          h(ii:ii+1,kk:kk+1) = h(ii:ii+1,kk:kk+1) + hi(:,2*k-1:2*k)
          g(ii:ii+1,kk:kk+1) = g(ii:ii+1,kk:kk+1) + gi(:,2*k-1:2*k)
        end do
        !
        deallocate(hi)
        deallocate(gi)
        !
    end do
    !
    ll = 0.0_8
    if (bie.eq.1) then
        !
        ! Diagonal terms of matrix H
        !
        h(ii,ii) = h(ii,ii) + 0.5_8
        h(ii+1,ii+1) = h(ii+1,ii+1) + 0.5_8
        else if (bie.eq.2) then
            !
            ! Diagonal terms of matrix H on the source point
            !
            h(ii,ii) = h(ii,ii) + 0.5_8
            h(ii+1,ii+1) = h(ii+1,ii+1) + 0.5_8
            !
            ! Diagonal terms of matrix H on the coincidental node
            !
            do c = 1,cr_num
                d = size(cr_list(c)%nodes_num(:,1))
                do b = 1,d
                    if (cr_list(c)%nodes_num(b,1) .eq. i) then
                        l = cr_list(c)%nodes_num(b,2)
                        ll = 2*l-1
                    end if
                end do
            end do
            h(ii,ll) = h(ii,ll) + 0.5_8
            h(ii+1,ll+1) = h(ii+1,ll+1) + 0.5_8
            !
            else if(bie.eq.3) then
                !
                !
                ! Diagonal terms of matrix G on the source point
                !
                g(ii,ii) = g(ii,ii) - 0.5_8
                g(ii+1,ii+1) = g(ii+1,ii+1) - 0.5_8
                !
                ! Diagonal terms of matrix H on the coincidental node
                !
                do c = 1,cr_num
                    d = size(cr_list(c)%nodes_num(:,1))
                    do b = 1,d
                        if (cr_list(c)%nodes_num(b,2) .eq. i) then
                            l = cr_list(c)%nodes_num(b,1)
                            ll = 2*l-1
                        end if
                    end do
                end do
                g(ii,ll) = g(ii,ll) + 0.5_8
                g(ii+1,ll+1) = g(ii+1,ll+1) + 0.5_8
                !
                else
                call error_message_stop("Incorrect Boundary Integral Equation for node.")
    end if
    end do
    !
    !---------------------------------------------------
    !Loop of old nodes on new boundary elements:
    !---------------------------------------------------
    !
    do i = 1,tnbn
      !
      cp = bn_list(i)
      bie = cp%bie
      ii = 2*i-1
      normal_cp = get_normal_to_a_boundary_point(i)
      !
      !----------------------------
      ! Loop on boundary elements:
      !----------------------------
      !
      do j = tnbe+1,be_num
        !
        be = be_list(j)
        !
        allocate(hi(2,2*be%nodes_num))
        allocate(gi(2,2*be%nodes_num))
        !
        hi = 0.0_8
        gi = 0.0_8
        if (bie.eq.1 .or. bie.eq.2) then
            call displ_bie_matrices_H_and_G(cp,be,hi,gi)
            else if (bie.eq.3) then
                call tract_bie_matrices_H_and_G(cp,be,normal_cp,hi,gi)
            else
                call error_message_stop("Incorrect Boundary Integral Equation for node.")
        end if
        !
        do k = 1,be%nodes_num
          kk = 2*be%incid(k)-1
          h(ii:ii+1,kk:kk+1) = h(ii:ii+1,kk:kk+1) + hi(:,2*k-1:2*k)
          g(ii:ii+1,kk:kk+1) = g(ii:ii+1,kk:kk+1) + gi(:,2*k-1:2*k)
        end do
        !
        deallocate(hi)
        deallocate(gi)
        !
    end do
    !
    end do
    !
    return
  end subroutine get_matrices_H_and_G_for_propagation
  !****************************************************
  !****************************************************
  subroutine set_initial_delta_a(p)
    !----------------------------------------------------
    ! Set the initial delta_a.
    !
    !
    ! Date: 06 December 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: p
    !
    Da_0 = p
    !
    return
  end subroutine set_initial_delta_a
  !****************************************************
  !****************************************************
  subroutine set_max_delta_a(p)
    !----------------------------------------------------
    ! Set the maximum delta_a.
    !
    !
    ! Date: 07 December 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: p
    !
    Da_max = p
    !
    return
  end subroutine set_max_delta_a
  !****************************************************
  !****************************************************
    subroutine set_r_Da(p)
    !----------------------------------------------------
    ! Set the scaling factor of Da.
    !
    !
    ! Date: 07 December 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: p
    !
    r_Da = p
    !
    return
  end subroutine set_r_Da
  !****************************************************
  !****************************************************
  subroutine set_parameter_C(p)
    !----------------------------------------------------
    ! Set the parameter C.
    !
    !
    ! Date: 18 November 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: p
    !
    C0 = p
    !
    return
  end subroutine set_parameter_C
  !****************************************************
  !****************************************************
  subroutine set_parameter_m(p)
    !----------------------------------------------------
    ! Set the parameter m.
    !
    !
    ! Date: 18 November 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: p
    !
    m = p
    !
    return
  end subroutine set_parameter_m
  !****************************************************
  !****************************************************
  subroutine get_parameter_m(p)
    !----------------------------------------------------
    ! Get the parameter m.
    !
    !
    ! Date: 19 December 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(out) :: p
    !
    p = m
    !
    return
  end subroutine get_parameter_m
  !****************************************************
  !****************************************************
  subroutine set_parameter_R(p)
    !----------------------------------------------------
    ! Set the parameter R.
    !
    !
    ! Date: 18 November 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: p
    !
    R = p
    !
    return
  end subroutine set_parameter_R
  !****************************************************
  !****************************************************
  subroutine get_parameter_R(p)
    !----------------------------------------------------
    ! Get the parameter R.
    !
    !
    ! Date: 19 December 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(out) :: p
    !
    p = R
    !
    return
  end subroutine get_parameter_R
  !****************************************************
  !****************************************************
  subroutine set_parameter_gamma(p)
    !----------------------------------------------------
    ! Set the parameter Gamma.
    !
    !
    ! Date: 07 December 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: p
    !
    gama = p
    !
    return
  end subroutine set_parameter_gamma
  !****************************************************
  !****************************************************
  subroutine get_parameter_gamma(p)
    !----------------------------------------------------
    ! Get the parameter Gamma.
    !
    !
    ! Date: 19 December 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(out) :: p
    !
    p = gama
    !
    return
  end subroutine get_parameter_gamma
  !****************************************************
  !****************************************************
  subroutine set_K_th(p)
    !----------------------------------------------------
    ! Set the parameter K_th.
    !
    !
    ! Date: 07 December 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: p
    !
    K_th = p
    !
    return
  end subroutine set_K_th
  !****************************************************
  !****************************************************
  subroutine set_K_IC(p)
    !----------------------------------------------------
    ! Set the parameter K_IC.
    !
    !
    ! Date: 07 December 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(in) :: p
    !
    K_IC = p
    !
    return
  end subroutine set_K_IC
  !****************************************************
  !****************************************************
  subroutine get_K_IC(p)
    !----------------------------------------------------
    ! Get the parameter K_IC.
    !
    !
    ! Date: 07 December 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(out) :: p
    !
    p = K_IC
    !
    return
  end subroutine get_K_IC
  !****************************************************
  !****************************************************
  subroutine get_number_of_cycles_per_step(i,j,number_of_cycles)
    !----------------------------------------------------
    ! Get the total number of cycles per step of
    ! propagation.
    !
    !
    ! Date: 18 November 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: i,j
    real(kind=8),intent(out) :: number_of_cycles
    !
    number_of_cycles = N_step(i,j)
    !
    return
  end subroutine get_number_of_cycles_per_step
  !****************************************************
  !****************************************************
  subroutine get_number_of_cycles_total(i,j,number_of_cycles)
    !----------------------------------------------------
    ! Get the total number of cycles during propagation.
    !
    !
    ! Date: 18 November 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: i,j
    real(kind=8),intent(out) :: number_of_cycles
    !
    number_of_cycles = N_total(i,j)
    !
    return
  end subroutine get_number_of_cycles_total
  !****************************************************
  !****************************************************
  subroutine set_k_eff_key(p)
    !----------------------------------------------------
    ! Set the key to evaluate effective K.
    !
    !
    ! Date: 15 December 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: p
    !
    k_eff_key = p
    !
    return
  end subroutine set_k_eff_key
  !****************************************************
  !****************************************************
  subroutine get_aq(aq)
    !----------------------------------------------------
    ! Get the multiplier factor aq.
    !
    !
    ! Date: 19 December 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(out) :: aq
    !
    aq = 1.0_8*(1-R)**gama
    !
    return
  end subroutine get_aq
  !****************************************************
  !****************************************************
end module fracture_2d_crack_propagation
