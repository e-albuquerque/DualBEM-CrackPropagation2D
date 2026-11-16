module fracture_2d_output
!----------------------------------------------------
! This module has routines to generate printout
! files for a two-dimensional elastic analysis.
!
! Date: August 2022
! Last update: ---
!----------------------------------------------------
use error_messages
use interpolation_functions
use fracture_2d_model
use fracture_2d_stress_intens_factor
use fracture_2d_crack_propagation
!
implicit none
!
real(kind=8),private :: pi = 3.1415926535897932_8
!
contains
  subroutine initial_data_output()
    !----------------------------------------------------
    ! This routine generates a file with the initial
    ! data.
    !
    ! Date: 09 August 2022
    ! Last update: ---
    !----------------------------------------------------
    character(len=200) :: file_name
    character(len=200) :: analysis_id
    character(len=200) :: analysis_name
    character(len=100) :: e_msg
    character(len=6) :: bc1,bc2
    character(len=2) :: dir1,dir2
    real(kind=8) :: pv1,pv2
    integer :: psk
    real(kind=8) :: el,nu
    integer :: i,j
    !
    !----------------------------------------
    ! Analysis type and material properties:
    !----------------------------------------
    !
    call get_plane_state_key(psk)
    call get_young_modulus(el)
    call get_poisson_ratio(nu)
    !
    !--------------------
    ! Creating the file:
    !--------------------
    !
    call get_model_id(analysis_id)
    call get_model_name(analysis_name)
    !
    file_name = trim(analysis_id) // "_input_data.txt"
    !
    open(10,file=file_name) !Opening the file
    !
    write(10,100)
    write(10,101)
    write(10,102) analysis_name
    !
    if (psk.eq.1) then
      write(10,104)
    else
      write(10,105)
    end if
    !
    write(10,106) el
    write(10,107) nu
    !
    write(10,111) be_num
    write(10,112) bn_num
    if (in_num.gt.0) write(10,113) in_num
    !
    !
    !-------------------------
    ! Printing elements data:
    !-------------------------
    !
    write(10,120)
    do i = 1,be_num
      write(10,121) i,be_list(i)%bp_1%x,be_list(i)%bp_1%y, &
                  & be_list(i)%bp_2%x,be_list(i)%bp_2%y
    end do
    !
    write(10,122)
    do i = 1,be_num
      write(10,123) i,be_list(i)%lamb_1,be_list(i)%lamb_2, &
                  & (be_list(i)%incid(j),j=1,size(be_list(i)%incid))
    end do
    !
    !-------------------------------
    ! Printing boundary nodes data:
    !-------------------------------
    !
    write(10,130)
    do i = 1,bn_num
      !
      if (bn_list(i)%bc_dir.eq.0) then
        dir1 = "_x"
        dir2 = "_y"
      else
        dir1 = "_n"
        dir2 = "_t"
      end if
      !
      select case (bn_list(i)%bc(1))
      case(0)
        bc1 = "disp" // dir1
        pv1 = bn_list(i)%disp(1)
      case(1)
        bc1 = "trac" // dir1
        pv1 = bn_list(i)%trac(1)
      case(2)
        bc1 = "sprg" // dir1
        pv1 = bn_list(i)%h1(1)
      case default
        write(e_msg,"(a49,i0)") "Invalid boundary condition key for boundary node ",i
        call error_message_stop(e_msg)
      end select
      !
      select case (bn_list(i)%bc(2))
      case(0)
        bc2 = "disp" // dir2
        pv2 = bn_list(i)%disp(2)
      case(1)
        bc2 = "trac" // dir2
        pv2 = bn_list(i)%trac(2)
      case(2)
        bc2 = "sprg" // dir2
        pv2 = bn_list(i)%h1(2)
      case default
        write(e_msg,"(a49,i0)") "Invalid boundary condition key for boundary node ",i
        call error_message_stop(e_msg)
      end select
      !
      write(10,131) i, bn_list(i)%point%x, bn_list(i)%point%y, bc1, pv1, bc2, pv2
      !
    end do
    !
    !--------------------------------
    ! Printing internal points data:
    !--------------------------------
    !
    if (in_num.gt.0) then
      !
      write(10,140)
      do i = 1,in_num
        write(10,141) i, in_list(i)%point%x, in_list(i)%point%y
      end do
      !
    end if
    !
    write(10,100)
    !
    close(10) !Closing the file
    !
    !----------
    ! Formats:
    !----------
    !
    100 format(82("*"))
    101 format("Fracture 2D",//,"A Dual-BEM based program - Version 1.0 - August 2022",/)
    102 format("Analysis name: ",(a),/)
    104 format(1x,"Plane stress problem.",/)
    105 format(1x,"Plane strain problem.",/)
    106 format(1x,"Elasticity modulus = ",e20.12)
    107 format(1x,"Poisson ratio = ",f0.8,/)
    !
    111 format(1x,"Number of boundary elements = ",i0)
    112 format(1x,"Number of boundary points = ",i0)
    113 format(1x,"Number of internal points = ",i0)
    !
    120 format(/,"Boundary Elements:",//,"Elem.",5x,"X1",12x,"Y1",12x,"X2",12x,"Y2")
    121 format(i5,4(2x,e12.5))
    122 format(/,"Elem.",5x,"Lambda 1",6x,"Lambda 2",6x,"Incidences")
    123 format(i5,2(2x,e12.5),3(2x,i7))
    !
    130 format(/,"Boundary Points:",//,"Point",5x,"X",13x,"Y",10x,"B.C. 1",4x,"Value",&
              &8x,"B.C. 2",4x,"Value")
    131 format(i5,2(2x,e12.5),2(2x,(a),1x,e14.7))
    !
    140 format(/,"Internal Points:",//,"Point",5x,"X",13x,"Y")
    141 format(i5,2(2x,e12.5))
    !
    return
  end subroutine initial_data_output
  !****************************************************
  !****************************************************
  subroutine results_output(t_sf,t_ss,t_ip,t_cr,n)
    !----------------------------------------------------
    ! This routine generates a text file with results
    ! for a two-dimensional isotropic linear elastic
    ! analysis
    !
    ! Date: 09 August 2022
    ! Last update: 28 October 2022
    !----------------------------------------------------
    real,intent(in) :: t_sf,t_ss,t_ip,t_cr
    integer,intent(in) :: n
    character(len=200) :: file_name
    character(len=200) :: analysis_id
    character(len=200) :: analysis_name
    integer :: psk
    real(kind=8) :: el,nu
    real(kind=8) :: s33,vms
    real(kind=8) :: ps(2)
    integer :: i,j,k,K_key
    !
    character(len=3) :: num
    !
    !----------------------------------------
    ! Analysis type and material properties:
    !----------------------------------------
    !
    call get_plane_state_key(psk)
    call get_young_modulus(el)
    call get_poisson_ratio(nu)
    !
    !--------------------
    ! Creating the file:
    !--------------------
    !
    call get_model_id(analysis_id)
    call get_model_name(analysis_name)
    !
    write(num,'(i3.3)') n
    !
    file_name = trim(analysis_id) // "_results_step_" // trim(num) // ".txt"
    !
    open(10,file=file_name) !Opening the file
    !
    write(10,100)
    write(10,101)
    write(10,102) analysis_name
    !
    if (psk.eq.1) then
      write(10,104)
    else
      write(10,105)
    end if
    !
    call get_K_key(K_key)
    call get_method_key(k)
    if (K_key.eq.1) then
        write (10,143) k
        else if (K_key.eq.2) then
            write(10,144) k
            else if (K_key.eq.3) then
                write(10,145) k
    end if
    !
    write(10,106) el
    write(10,107) nu
    !
    write(10,108) be_num
    write(10,109) bn_num
    if (in_num.gt.0) write(10,110) in_num
    !
    !
    write(10,113) t_sf,t_ss,t_ip,t_cr
    !
    !-------------------------------
    ! Printing boundary nodes data:
    !-------------------------------
    !
    write(10,120)
    do i = 1,bn_num
      write(10,121) i, bn_list(i)%point%x, bn_list(i)%point%y, &
                    & bn_list(i)%disp(1), bn_list(i)%disp(2), &
                    & bn_list(i)%trac(1), bn_list(i)%trac(2)
    end do
    !
    if (psk.eq.1) then
      !Plane stress:
      write(10,132)
      do i = 1,bn_num
        s33 = out_of_plane_strain(bn_list(i)%stress)
        write(10,121) i, bn_list(i)%point%x, bn_list(i)%point%y, &
                    & bn_list(i)%strain(1), bn_list(i)%strain(2), bn_list(i)%strain(3), s33
      end do
      !
      write(10,133)
      do i = 1,bn_num
        write(10,121) i, bn_list(i)%point%x, bn_list(i)%point%y, &
                    & bn_list(i)%stress(1), bn_list(i)%stress(2), bn_list(i)%stress(3)
      end do
      !
    else
      !Plane strain:
      write(10,131)
      do i = 1,bn_num
        write(10,121) i, bn_list(i)%point%x, bn_list(i)%point%y, &
                    & bn_list(i)%strain(1), bn_list(i)%strain(2), bn_list(i)%strain(3)
      end do
      !
      write(10,134)
      do i = 1,bn_num
        s33 = out_of_plane_stress(bn_list(i)%stress)
        write(10,121) i, bn_list(i)%point%x, bn_list(i)%point%y, &
                    & bn_list(i)%stress(1), bn_list(i)%stress(2), bn_list(i)%stress(3), s33
      end do
      !
    end if
    !
    write(10,135)
    do i = 1,bn_num
      ps = eigenvalue(bn_list(i)%strain)
      write(10,121) i, bn_list(i)%point%x, bn_list(i)%point%y, ps(1), ps(2)
    end do
    !
    write(10,136)
    do i = 1,bn_num
      ps = eigenvalue(bn_list(i)%stress)
      vms = von_mises_stress(bn_list(i)%stress)
      write(10,121) i, bn_list(i)%point%x, bn_list(i)%point%y, ps(1), ps(2), vms
    end do
    !
    !--------------------------------
    ! Printing internal points data:
    !--------------------------------
    !
    if (in_num.gt.0) then
      !
      write(10,130)
      do i =1,in_num
        write(10,121) i, in_list(i)%point%x, in_list(i)%point%y, &
                      & in_list(i)%disp(1), in_list(i)%disp(2)
      end do
      !
      if (psk.eq.1) then
        !Plane stress:
        write(10,132)
        do i = 1,in_num
          s33 = out_of_plane_strain(in_list(i)%stress)
          write(10,121) i, in_list(i)%point%x, in_list(i)%point%y, &
                      & in_list(i)%strain(1), in_list(i)%strain(2), in_list(i)%strain(3), s33
        end do
        !
        write(10,133)
        do i = 1,in_num
          write(10,121) i, in_list(i)%point%x, in_list(i)%point%y, &
                      & in_list(i)%stress(1), in_list(i)%stress(2), in_list(i)%stress(3)
        end do
        !
      else
        !Plane strain:
        write(10,131)
        do i = 1,in_num
          write(10,121) i, in_list(i)%point%x, in_list(i)%point%y, &
                      & in_list(i)%strain(1), in_list(i)%strain(2), in_list(i)%strain(3)
        end do
        !
        write(10,134)
        do i = 1,in_num
          s33 = out_of_plane_stress(in_list(i)%stress)
          write(10,121) i, in_list(i)%point%x, in_list(i)%point%y, &
                      & in_list(i)%stress(1), in_list(i)%stress(2), in_list(i)%stress(3), s33
        end do
        !
      end if
      !
      write(10,135)
      do i = 1,in_num
        ps = eigenvalue(in_list(i)%strain)
        write(10,121) i, in_list(i)%point%x, in_list(i)%point%y, ps(1), ps(2)
      end do
      !
      write(10,136)
      do i = 1,in_num
        ps = eigenvalue(in_list(i)%stress)
        vms = von_mises_stress(in_list(i)%stress)
        write(10,121) i, in_list(i)%point%x, in_list(i)%point%y, ps(1), ps(2), vms
      end do
      !
    end if
    !
    !----------------------------------------
    ! Printing stress intensity factor data:
    !----------------------------------------
    !
    write(10,140)
    write(10,141)
    do i = 1,cr_num
        do j = 1,cr_list(i)%tips_num
            write(10,142) i, j, cr_list(i)%ki(j), cr_list(i)%kii(j), cr_list(i)%theta(j)*180/pi,cr_list(i)%keq(j)
        end do
    end do
    write(10,100)
    close(10) !Closing the file
    !
    !----------
    ! Formats:
    !----------
    !
    100 format(100("*"))
    101 format("Fracture 2D",//,"A Dual-BEM based program - Version 2.0 - October 2022",/)
    102 format("Analysis name: ",(a),/)
    104 format(1x,"Plane stress problem.",/)
    105 format(1x,"Plane strain problem.",/)
    106 format(1x,"Elasticity modulus = ",e20.12)
    107 format(1x,"Poisson ratio = ",f0.8,/)
    108 format(1x,"Number of boundary elements = ",i0)
    109 format(1x,"Number of boundary points = ",i0)
    110 format(1x,"Number of internal points = ",i0)
    113 format(/,1x,"Time elapsed to assemble matrices = ",f0.3," sec",&
              &/,1x,"Time elapsed to solve the system = ",f0.3," sec",&
              &/,1x,"Time elapsed to evaluate internal points = ",f0.3," sec",&
              &/,1x,"Time elapsed to evaluate cracks = ",f0.3," sec")
    !
    120 format(/,"Boundary Points:",//,"Point",5x,"X",13x,"Y",13x, &
              &  "Disp. X",9x,"Disp. Y",9x,"Trac. X",9x,"Trac. Y")
    121 format(i5,2(2x,e12.5),4(2x,e14.7))
    130 format(/,"Internal Points:",//,"Point",5x,"X",13x,"Y",13x,"Disp. X",9x,"Disp. Y")
    131 format(/,"Point",5x,"X",13x,"Y",13x,"Strain XX",7x,"Strain YY",7x,"Strain XY")
    132 format(/,"Point",5x,"X",13x,"Y",13x,"Strain XX",7x,"Strain YY",7x,"Strain XY",7x,"Strain ZZ")
    133 format(/,"Point",5x,"X",13x,"Y",13x,"Stress XX",7x,"Stress YY",7x,"Stress XY")
    134 format(/,"Point",5x,"X",13x,"Y",13x,"Stress XX",7x,"Stress YY",7x,"Stress XY",7x,"Stress ZZ")
    135 format(/,"Point",5x,"X",13x,"Y",11x,"Pr. Strain 1",4x,"Pr. Strain 2")
    136 format(/,"Point",5x,"X",13x,"Y",11x,"Pr. Stress 1",4x,"Pr. Stress 2",6x,"von Mises")
    140 format(/,"Cracks:",//,"Stress Intensity Factors:")
    141 format(/,"Crack",3x,"Tip",7x,"KI",14x,"KII",12x,"Theta",12x,"Keq")
    142 format(i5,i5,4(2x,e14.7))
    143 format(1x,"Stress intensity factor evaluated by the displacement correlation method: ",i0,/)
    144 format(/,"Stress intensity factor evaluated by the J-Integral method: ",i0,/)
    145 format(/,"Stress intensity factor evaluated by the displacement fitting method: ",i0,/)
    !
    return
  end subroutine results_output
  !****************************************************
  !****************************************************
  subroutine results_per_crack_output(n,amax)
    !----------------------------------------------------
    ! This routine generates a text file with results
    ! for the stress intensity factors.
    !
    ! Date: 09 August 2022
    ! Last update: 20 December 2022
    !----------------------------------------------------
    integer,intent(in) :: n
    real(kind=8),intent(in) :: amax
    character(len=200) :: analysis_id
    character(len=200) :: file_name
    real(kind=8) :: length,ki,kii,ki_2,kii_2,keq,keq_2,K_IC
    integer :: i
    !
    character(len=3) :: num
    !
    !
    call get_model_id(analysis_id)
    !
    if (n .eq. 0) then
        do i = 1,cr_num
            write(num,'(i3.3)') i
            !
            file_name = trim(analysis_id) // "_results_crack_" // trim(num) // ".txt"
            !
            open(i,file=file_name) !Opening the file
            !
            if (cr_list(i)%tips_num .eq. 1) then
                write(i,101)
                else
                    write(i,102)
            end if
            length = cr_list(i)%length
            ki = amax*cr_list(i)%ki(1)
            kii = amax*cr_list(i)%kii(1)
            ki_2 = amax*cr_list(i)%ki(2)
            kii_2 = amax*cr_list(i)%kii(2)
            keq = amax*cr_list(i)%keq(1)
            keq_2 = amax*cr_list(i)%keq(2)
            if (cr_list(i)%tips_num .eq. 1) then
                write(i,103) n, length, ki, kii, keq
                else
                    write(i,104) n, length, ki, kii, keq, ki_2, kii_2, keq_2
            end if
            !
            close(i)
            !
        end do
        else if(n .gt. 0) then
            do i = 1,cr_num
                write(num,'(i3.3)') i
                !
                file_name = trim(analysis_id) // "_results_crack_" // trim(num) // ".txt"
                !
                open(i,file=file_name,status='OLD',access='append') !Opening the file
                !
                length = cr_list(i)%length
                ki = amax*cr_list(i)%ki(1)
                kii = amax*cr_list(i)%kii(1)
                ki_2 = amax*cr_list(i)%ki(2)
                kii_2 = amax*cr_list(i)%kii(2)
                keq = amax*cr_list(i)%keq(1)
                keq_2 = amax*cr_list(i)%keq(2)
                call get_K_IC(K_IC)
                if (keq .gt. K_IC .or. keq_2 .gt. K_IC) then
                    if (cr_list(i)%tips_num .eq. 1) then
                        write(i,105) n, length, ki, kii, keq
                        else
                            write(i,106) n, length, ki, kii, keq, ki_2, kii_2, keq_2
                    end if
                    else
                        if (cr_list(i)%tips_num .eq. 1) then
                            write(i,103) n, length, ki, kii, keq
                            else
                                write(i,104) n, length, ki, kii, keq, ki_2, kii_2, keq_2
                    end if
                end if
                close(i)
                !
            end do
            !
    end if
    !
    !----------
    ! Formats:
    !----------
    !
    101 format(2x,"Step",5x,"Length",11x,"K_I_1",11x,"K_II_1",11x, "Keq_1")
    102 format(2x,"Step",5x,"Length",11x,"K_I_1",11x,"K_II_1",11x, "Keq_1", 11x,"K_I_2",10x,"K_II_2", 11x, "Keq_2")
    103 format(i5,4(2x,e14.7))
    104 format(i5,7(2x,e14.7))
    105 format(i5,4(2x,e14.7),5x,"**** Keq greater than K_IC!****")
    106 format(i5,7(2x,e14.7),5x,"**** Keq greater than K_IC!****")
    !
    return
  end subroutine results_per_crack_output
  !****************************************************
  !****************************************************
  subroutine crack_geometry_output(n)
    !----------------------------------------------------
    ! This routine generates a text file with the
    ! geometry of the crack.
    !
    ! Date: 20 December 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: n
    character(len=200) :: analysis_id
    character(len=200) :: file_name
    real(kind=8) :: x_initial,y_initial,x_final,y_final,length
    integer :: i
    !
    character(len=3) :: num
    !
    !
    call get_model_id(analysis_id)
    !
    if (n .eq. 0) then
        do i = 1,cr_num
            write(num,'(i3.3)') i
            !
            file_name = trim(analysis_id) // "_geometry_crack_" // trim(num) // ".txt"
            !
            open(i,file=file_name) !Opening the file
            !
            write(i,101)
            !
            x_initial = be_list(cr_list(i)%ext_el_num(4))%bp_2%x
            y_initial = be_list(cr_list(i)%ext_el_num(4))%bp_2%y
            x_final = be_list(cr_list(i)%ext_el_num(2))%bp_2%x
            y_final = be_list(cr_list(i)%ext_el_num(2))%bp_2%y
            length = cr_list(i)%length
            !
            write(i,102) n, x_initial, y_initial, x_final, y_final, length
            !
            close(i)
            !
        end do
        else if(n .gt. 0) then
            do i = 1,cr_num
                write(num,'(i3.3)') i
                !
                file_name = trim(analysis_id) // "_geometry_crack_" // trim(num) // ".txt"
                !
                open(i,file=file_name,status='OLD',access='append') !Opening the file
                !
                x_initial = be_list(cr_list(i)%ext_el_num(4))%bp_2%x
                y_initial = be_list(cr_list(i)%ext_el_num(4))%bp_2%y
                x_final = be_list(cr_list(i)%ext_el_num(2))%bp_2%x
                y_final = be_list(cr_list(i)%ext_el_num(2))%bp_2%y
                length = cr_list(i)%length
                !
                write(i,102) n, x_initial, y_initial, x_final, y_final, length
                !
                close(i)
                !
            end do
            !
    end if
    !
    !----------
    ! Formats:
    !----------
    !
    101 format(2x,"Step",3x,"Initial_x",5x,"Initial_y",6x, "Final_x",7x, "Final_y",8x,"Length")
    102 format(i5,5(e14.5))
    !
    return
  end subroutine crack_geometry_output
  !****************************************************
  !****************************************************
  subroutine fatigue_results_per_crack_output(n)
    !----------------------------------------------------
    ! This routine generates a text file with results
    ! for a two-dimensional isotropic linear elastic
    ! analysis
    !
    ! Date: 09 August 2022
    ! Last update: 28 October 2022
    !----------------------------------------------------
    integer,intent(in) :: n
    character(len=200) :: analysis_id
    character(len=200) :: file_name
    real(kind=8) :: length,N_1_step,N_2_step,N_1_total,N_2_total
    integer :: i

    !
    character(len=3) :: num
    !
    !
    call get_model_id(analysis_id)
    !
    if (n .eq. 0) then
        do i = 1,cr_num
            write(num,'(i3.3)') i
            !
            file_name = trim(analysis_id) // "_fatigue_results_crack_" // trim(num) // ".txt"
            !
            open(i,file=file_name) !Opening the file
            !
            if (cr_list(i)%tips_num .eq. 1) then
                write(i,101)
                else
                    write(i,102)
            end if
            length = cr_list(i)%length
            N_1_step = 0.0
            N_1_total = N_1_step
            N_2_step = 0.0
            N_2_total = N_2_step
            !
            if (cr_list(i)%tips_num .eq. 1) then
                write(i,103) length,N_1_step,N_1_total
                else
                    write(i,104) length,N_1_step,N_1_total,N_2_step,N_2_total
            end if
            !
            close(i)
            !
        end do
        else if(n .gt. 0) then
            do i = 1,cr_num
                write(num,'(i3.3)') i
                !
                file_name = trim(analysis_id) // "_fatigue_results_crack_" // trim(num) // ".txt"
                !
                open(i,file=file_name,status='OLD',access='append') !Opening the file
                !
                length = cr_list(i)%length
                call get_number_of_cycles_per_step(i,1,N_1_step)
                call get_number_of_cycles_total(i,1,N_1_total)
                call get_number_of_cycles_per_step(i,2,N_2_step)
                call get_number_of_cycles_total(i,2,N_2_total)
                if (cr_list(i)%tips_num .eq. 1) then
                    write(i,103) length,N_1_step,N_1_total
                    else
                        write(i,104) length,N_1_step,N_1_total,N_2_step,N_2_total
                end if
                close(i)
                !
            end do
            !
    end if
    !
    !----------
    ! Formats:
    !----------
    !
    101 format(5x,"Length",10x,"N_1_step",8x,"N_1_total")
    102 format(5x,"Length",10x,"N_1_step",8x,"N_1_total",8x,"N_2_step",8x,"N_2_total")
    103 format(3(2x,e14.7))
    104 format(5(2x,e14.7))
    !
    return
  end subroutine fatigue_results_per_crack_output
  !****************************************************
  !****************************************************
  subroutine cycles_output(cycles,alpha_max,R)
    !----------------------------------------------------
    ! This routine generates a text file with the cycles
    ! counted by the rainflow cycle counting.
    !
    ! Date: 19 December 2022
    ! Last update: ---
    !----------------------------------------------------
    integer,intent(in) :: cycles
    real(kind=8),intent(in) :: alpha_max(:),R(:)
    character(len=200) :: analysis_id
    character(len=200) :: file_name
    integer :: i
    !
    !
    call get_model_id(analysis_id)
    !
    file_name = trim(analysis_id) // "_cycles.txt"
    !
    open(10,file=file_name) !Opening the file
    !
    write(10,101)
    do i = 1,cycles
        write(10,102) i,alpha_max(i),R(i)
    end do
    !
    close(10)
    !
    !----------
    ! Formats:
    !----------
    !
    101 format("Cycle",5x,"alpha_max",8x,"R")
    102 format(i5,2(2x,e14.7))
    !
    return
  end subroutine cycles_output
  !****************************************************
  !****************************************************
  subroutine output_file_for_insane(n)
    !----------------------------------------------------
    ! Write a file for INSANE interface with mechanical
    ! results.
    !
    ! Date: 09 August 2022
    ! Last update: 28 October 2022
    !----------------------------------------------------
    integer,intent(in) :: n
    type(boundary_node),allocatable :: geom_nodes_list(:) !A local list of geometric nodes to treat discontinuous elements.
    integer :: pseudo_elem_num !The number of pseudo-elements to plot discontinuous elements results.
    real(kind=8),allocatable :: ni(:)
    real(kind=8) :: el,nu
    real(kind=8) :: eps_33,sig_33,vms
    real(kind=8) :: ps(2)
    integer :: psk
    integer :: i,j,k
    !
    character(len=200) :: file_name
    character(len=200) :: analysis_id
    !
    character(len=38) :: header_1 = "<?xml version='1.1' encoding='UTF-8'?>"
    character(len=40) :: header_2 = "<Insane xmlns='http://www.dees.ufmg.br' "
    character(len=54) :: header_3 = "xmlns:xsi='http://www.w3.org/2001/XMLSchema-instance' "
    character(len=56) :: header_4 = "xsi:schemaLocation='http://www.dees.ufmg.br insane.xsd'>"
    !
    character(len=14) :: an_type
    !
    character(len=5) :: bound_cond(2)
    !
    character(len=32) :: linear_bound_elem = "ParametricElement.BemElement.BL2"
    character(len=32) :: quadratic_bound_elem = "ParametricElement.BemElement.BL3"
    !
    character(len=3) :: num
    !
    !-----------------------------------------
    ! Allocating the geometric nodes list
    ! to plot discontinuous elements results:
    !-----------------------------------------
    !
    pseudo_elem_num = 0
    !
    do i = 1,be_num
      select case (be_list(i)%nodes_num)
      case (1)
        pseudo_elem_num = pseudo_elem_num + 2
      case default
        if (be_list(i)%lamb_1.lt.1._8) pseudo_elem_num = pseudo_elem_num + 1
        if (be_list(i)%lamb_2.lt.1._8) pseudo_elem_num = pseudo_elem_num + 1
      end select
    end do
    !
    allocate(geom_nodes_list(pseudo_elem_num))
    !
    do i = 1,pseudo_elem_num
      geom_nodes_list(i)%disp = 0.0
      geom_nodes_list(i)%trac = 0.0
      geom_nodes_list(i)%strain = 0.0
      geom_nodes_list(i)%stress = 0.0
    end do
    !
    !-----------------------
    ! Opening the xml file:
    !-----------------------
    !
    call get_model_id(analysis_id)
    !
    write(num,'(i3.3)') n
    !
    file_name = trim(analysis_id) // "-STEP-" // trim(num) // ".xml"
    !
    open(10,file=file_name)
    !
    !---------------------
    ! Writing the header:
    !---------------------
    !
    write(10,100) header_1
    write(10,101) header_2,header_3,header_4
    !
    !-----------------------------
    ! Solution class information:
    !-----------------------------
    !
    write(10,110)
    write(10,111)
    write(10,112)
    !
    !--------------------------
    ! Model class information:
    !--------------------------
    !
    write(10,120)
    !
    write(10,130)
    !
    call get_plane_state_key(psk)
    !
    select case (psk)
    case(1)
      an_type = "BemPlaneStress"
    case(2)
      an_type = "BemPlaneStrain"
    case default
      call error_message_stop("Undefined analysis type for xml file generation.")
    end select
    !
    write(10,140) an_type
    !
    write(10,150)
    write(10,151)
    write(10,152)
    write(10,153)
    write(10,154)
    write(10,155)
    write(10,156)
    !
    ! - Material information:
    !
    call get_young_modulus(el)
    call get_poisson_ratio(nu)
    !
    write(10,160)
    write(10,161)
    write(10,162) el
    write(10,163) nu
    write(10,164)
    write(10,165)
    !
    ! - Elements information:
    !
    write(10,210)
    !
    k = 0
    !
    do i = 1,be_num
      !
      ! Creating a pseudo-element "before" a discontinuous or constant element:
      !
      if (be_list(i)%nodes_num.eq.1 .or. be_list(i)%lamb_1.lt.1._8) then
        !
        k = k + 1
        !
        geom_nodes_list(k)%point = be_list(i)%bp_1
        geom_nodes_list(k)%bc_dir = bn_list(be_list(i)%incid(1))%bc_dir
        geom_nodes_list(k)%bc = bn_list(be_list(i)%incid(1))%bc
        !
        allocate(ni(be_list(i)%nodes_num))
        !
        call intrp_func_1d(ni,-1._8,be_list(i)%lamb_1,be_list(i)%lamb_2,be_list(i)%nodes_num-1)
        !
        do j = 1,be_list(i)%nodes_num
          geom_nodes_list(k)%disp = geom_nodes_list(k)%disp + ni(j)*bn_list(be_list(i)%incid(j))%disp
          geom_nodes_list(k)%trac = geom_nodes_list(k)%trac + ni(j)*bn_list(be_list(i)%incid(j))%trac
          geom_nodes_list(k)%strain = geom_nodes_list(k)%strain + ni(j)*bn_list(be_list(i)%incid(j))%strain !Verify this strategy
          geom_nodes_list(k)%stress = geom_nodes_list(k)%stress + ni(j)*bn_list(be_list(i)%incid(j))%stress !Verify this strategy
        end do
        !
        deallocate(ni)
        !
        write(10,224) linear_bound_elem
        write(10,222) k+bn_num+in_num, be_list(i)%incid(1)
        !
        write(10,230) an_type
        !
        write(10,221)
        !
      end if
      !
      ! Real element information:
      !
      select case (be_list(i)%nodes_num)
      case (1)
        ! Constant elements results are plotted by pseudo-elements only
      case (2)
        write(10,220) linear_bound_elem, i
        write(10,222) (be_list(i)%incid(j), j=1,be_list(i)%nodes_num)
        write(10,230) an_type
        write(10,221)
      case (3)
        write(10,220) quadratic_bound_elem, i
        write(10,223) be_list(i)%incid(1), be_list(i)%incid(3), be_list(i)%incid(2)
        write(10,230) an_type
        write(10,221)
      case default
        call warning_message("Undefined boundary element type for *.xml file generation.")
        write(10,*)
        write(10,*) "*** Undefined boundary element type"
        write(10,*) "    Number of nodes = ",be_list(i)%nodes_num
        write(10,*)
      end select
      !
      ! Creating a pseudo-element "after" a discontinuous or constant element:
      !
      if (be_list(i)%nodes_num.eq.1 .or. be_list(i)%lamb_2.lt.1._8) then
        !
        k = k + 1
        !
        geom_nodes_list(k)%point = be_list(i)%bp_2
        geom_nodes_list(k)%bc_dir = bn_list(be_list(i)%incid(be_list(i)%nodes_num))%bc_dir
        geom_nodes_list(k)%bc = bn_list(be_list(i)%incid(be_list(i)%nodes_num))%bc
        !
        allocate(ni(be_list(i)%nodes_num))
        !
        call intrp_func_1d(ni,+1._8,be_list(i)%lamb_1,be_list(i)%lamb_2,be_list(i)%nodes_num-1)
        !
        do j = 1,be_list(i)%nodes_num
          geom_nodes_list(k)%disp = geom_nodes_list(k)%disp + ni(j)*bn_list(be_list(i)%incid(j))%disp
          geom_nodes_list(k)%trac = geom_nodes_list(k)%trac + ni(j)*bn_list(be_list(i)%incid(j))%trac
          geom_nodes_list(k)%strain = geom_nodes_list(k)%strain + ni(j)*bn_list(be_list(i)%incid(j))%strain !Verify this strategy
          geom_nodes_list(k)%stress = geom_nodes_list(k)%stress + ni(j)*bn_list(be_list(i)%incid(j))%stress !Verify this strategy
        end do
        !
        deallocate(ni)
        !
        write(10,224) linear_bound_elem
        write(10,222) be_list(i)%incid(be_list(i)%nodes_num), k+bn_num+in_num
        !
        write(10,230) an_type
        !
        write(10,221)
        !
      end if
      !
    end do
    !
    write(10,211)
    !
    ! - Nodes information:
    !
    write(10,170)
    !
    ! Boundary nodes:
    !
    do i = 1,bn_num
      !
      do j = 1,2
        !
        select case (bn_list(i)%bc(j))
        case(0)
          bound_cond(j) = "true"
        case default
          bound_cond(j) = "false"
        end select
        !
      end do
      !
      write(10,180) i
      write(10,181) bn_list(i)%point%x, bn_list(i)%point%y
      write(10,182)
      write(10,183)
      write(10,184)
      write(10,185)
      write(10,186) bound_cond(1), bound_cond(2)
      write(10,187) bn_list(i)%disp(1), bn_list(i)%disp(2)
      write(10,188) bn_list(i)%trac(1), bn_list(i)%trac(2)
      !
      if (psk.eq.1) then
        !Plane stress:
        sig_33 = 0.0
        eps_33 = out_of_plane_strain(bn_list(i)%stress)
      else
        !Plane strain:
        eps_33 = 0.0
        sig_33 = out_of_plane_stress(bn_list(i)%stress)
      end if
      !
      write(10,189) bn_list(i)%strain(1), bn_list(i)%strain(2), bn_list(i)%strain(3), eps_33
      write(10,190) bn_list(i)%stress(1), bn_list(i)%stress(2), bn_list(i)%stress(3), sig_33
      write(10,193)
      !
      ps = eigenvalue(bn_list(i)%strain)
      write(10,194) ps(1), ps(2)
      !
      ps = eigenvalue(bn_list(i)%stress)
      write(10,195) ps(1), ps(2)
      !
      vms = von_mises_stress(bn_list(i)%stress)
      write(10,196) vms
      !
      write(10,191)
      write(10,192)
      !
    end do
    !
    ! Internal points:
    !
    do i = 1,in_num
      !
      j = i + bn_num
      write(10,180) j
      write(10,181) in_list(i)%point%x, in_list(i)%point%y
      write(10,182)
      write(10,200)
      write(10,184)
      write(10,185)
      write(10,187) in_list(i)%disp(1), in_list(i)%disp(2)
      !
      if (psk.eq.1) then
        !Plane stress:
        sig_33 = 0.0
        eps_33 = out_of_plane_strain(in_list(i)%stress)
      else
        !Plane strain:
        eps_33 = 0.0
        sig_33 = out_of_plane_stress(in_list(i)%stress)
      end if
      !
      write(10,189) in_list(i)%strain(1), in_list(i)%strain(2), in_list(i)%strain(3), eps_33
      write(10,190) in_list(i)%stress(1), in_list(i)%stress(2), in_list(i)%stress(3), sig_33
      write(10,193)
      !
      ps = eigenvalue(in_list(i)%strain)
      write(10,194) ps(1), ps(2)
      !
      ps = eigenvalue(in_list(i)%stress)
      write(10,195) ps(1), ps(2)
      !
      vms = von_mises_stress(in_list(i)%stress)
      write(10,196) vms
      !
      write(10,191)
      write(10,192)
      !
    end do
    !
    ! Geometric points for discontinuous and constant elements plot:
    !
    do i = 1,pseudo_elem_num
      !
      do j = 1,2
        !
        select case (geom_nodes_list(i)%bc(j))
        case(0)
          bound_cond(j) = "true"
        case default
          bound_cond(j) = "false"
        end select
      end do
      !
      write(10,180) i+bn_num+in_num
      write(10,181) geom_nodes_list(i)%point%x, geom_nodes_list(i)%point%y
      write(10,182)
      write(10,183)
      write(10,184)
      write(10,185)
      write(10,186) bound_cond(1), bound_cond(2)
      write(10,187) geom_nodes_list(i)%disp(1), geom_nodes_list(i)%disp(2)
      write(10,188) geom_nodes_list(i)%trac(1), geom_nodes_list(i)%trac(2)
      !
      if (psk.eq.1) then
        !Plane stress:
        sig_33 = 0.0
        eps_33 = out_of_plane_strain(geom_nodes_list(i)%stress)
      else
        !Plane strain:
        eps_33 = 0.0
        sig_33 = out_of_plane_stress(geom_nodes_list(i)%stress)
      end if
      !
      write(10,189) geom_nodes_list(i)%strain(1), geom_nodes_list(i)%strain(2), geom_nodes_list(i)%strain(3), eps_33
      write(10,190) geom_nodes_list(i)%stress(1), geom_nodes_list(i)%stress(2), geom_nodes_list(i)%stress(3), sig_33
      write(10,193)
      !
      ps = eigenvalue(geom_nodes_list(i)%strain)
      write(10,194) ps(1), ps(2)
      !
      ps = eigenvalue(geom_nodes_list(i)%stress)
      write(10,195) ps(1), ps(2)
      !
      vms = von_mises_stress(geom_nodes_list(i)%stress)
      write(10,196) vms
      !
      write(10,191)
      write(10,192)
      !
    end do
    !
    write(10,171)
    !
    write(10,121)
    !
    !-----------------------
    ! Closing the xml file:
    !-----------------------
    !
    write(10,900)
    !
    close(10)
    !
    deallocate(geom_nodes_list)
    !
    !----------
    ! Formats:
    !----------
    !
    100 format ((a))
    101 format ((a),(a),(a))
    !
    110 format (2x,"<Solution class='BemSteadyState'>")
    111 format (4x,"<SolverType>1</SolverType>")
    112 format (2x,"</Solution>")
    !
    120 format (2x,"<Model class='BemModel'>")
    121 format (2x,"</Model>")
    !
    130 format (4x,"<ProblemDriver>BemProblemDriver</ProblemDriver>")
    !
    140 format (4x,"<GlobalAnalysisModel>",(a),"</GlobalAnalysisModel>")
    !
    150 format (4x,"<DegenerationList>")
    151 format (6x,"<Degeneration class='PrescribedDegeneration' label='PrescribedDegeneration'>")
    152 format (8x,"<Height>1.0000000</Height>")
    153 format (8x,"<CSMaterial>LinearElasticIsotropic</CSMaterial>")
    154 format (8x,"<Thickness>1.0000000</Thickness>")
    155 format (6x,"</Degeneration>")
    156 format (4x,"</DegenerationList>")
    !
    160 format (4x,"<MaterialList>")
    161 format (6x,"<Material class='LinearElasticIsotropic' label='LinearElasticIsotropic'>")
    162 format (8x,"<Elasticity>",e20.10,"</Elasticity>")
    163 format (8x,"<Poisson>",e20.10,"</Poisson>")
    164 format (6x,"</Material>")
    165 format (4x,"</MaterialList>")
    !
    170 format (4x,"<BemNodeList>")
    171 format (4x,"</BemNodeList>")
    !
    180 format (6x,"<BemNode label='",i0,"'>")
    181 format (8x,"<Coord>",e20.10,1x,e20.10,1x,"0.0</Coord>")
    182 format (8x,"<BemNodeValues>")
    183 format (10x,"<Internal>false</Internal>")
    184 format (10x,"<DOFLabels>Dx Dy</DOFLabels>")
    185 format (10x,"<StrainStressLabels>XX YY XY ZZ</StrainStressLabels>")
    186 format (10x,"<Restraints>",(a),1x,(a),"</Restraints>")
    187 format (10x,"<Displ>",e20.10,1x,e20.10,"</Displ>")
    188 format (10x,"<Tract>",e20.10,1x,e20.10,"</Tract>")
    189 format (10x,"<Strain>",e20.10,1x,e20.10,1x,e20.10,1x,e20.10,"</Strain>")
    190 format (10x,"<Stress>",e20.10,1x,e20.10,1x,e20.10,1x,e20.10,"</Stress>")
    191 format (8x,"</BemNodeValues>")
    192 format (6x,"</BemNode>")
    193 format (10x,"<PrStrainStressLabels>1 2</PrStrainStressLabels>")
    194 format (10x,"<PrStrain>",e20.10,1x,e20.10,"</PrStrain>")
    195 format (10x,"<PrStress>",e20.10,1x,e20.10,"</PrStress>")
    196 format (10x,"<VonMisesStress>",e20.10,"</VonMisesStress>")
    !
    200 format (10x,"<Internal>true</Internal>")
    !
    210 format (4x,"<ElementList>")
    211 format (4x,"</ElementList>")
    !
    220 format (6x,"<Element class='",(a),"' label='E",i0,"'>")
    221 format (6x,"</Element>")
    222 format (8x,"<Incidence>",i0,1x,i0,"</Incidence>")
    223 format (8x,"<Incidence>",i0,1x,i0,1x,i0,"</Incidence>")
    224 format (6x,"<Element class='",(a),"' label=' '>")
    !
    230 format (8x,"<AnalysisModel>",(a),"</AnalysisModel>")
    !
    900 format ("</Insane>")
    !
    return
  end subroutine output_file_for_insane
  !****************************************************
  !****************************************************
  function von_mises_stress(sig)
    !----------------------------------------------------
    ! Evaluate von Mises equivalent stress for a given
    ! strain/stress state.
    !
    ! Date: 27 June 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8) :: von_mises_stress
    real(kind=8),intent(in) :: sig(3)
    real(kind=8) :: s33
    real(kind=8) :: s1,s2,s3
    integer :: psk
    !
    call get_plane_state_key(psk)
    s33 = 0.0
    !
    select case (psk)
    case(1)
      s33 = 0.0_8
    case(2)
      s33 = out_of_plane_stress(sig)
    case default
      call error_message_stop("Undefined plane state key for von Mises stress evaluation.")
    end select
    !
    s1 = (sig(1)-sig(2))*(sig(1)-sig(2))
    s2 = (sig(2)-s33)*(sig(2)-s33)
    s3 = (s33-sig(1))*(s33-sig(1))
    !
    von_mises_stress = sqrt(0.5_8*(s1+s2+s3)+3.0_8*sig(3)*sig(3))
    !
    return
  end function von_mises_stress
  !****************************************************
  !****************************************************
  function out_of_plane_strain(sig)
    !----------------------------------------------------
    ! Evaluate the out-of-plane strain component for
    ! a plane stress analysis.
    !
    ! Date: 27 June 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8) :: out_of_plane_strain
    real(kind=8),intent(in) :: sig(3)
    real(kind=8) :: el,nu
    !
    call get_young_modulus(el)
    call get_poisson_ratio(nu)
    !
    out_of_plane_strain = -nu*(sig(1)+sig(2))/el
    !
    return
  end function out_of_plane_strain
  !****************************************************
  !****************************************************
  function out_of_plane_stress(sig)
    !----------------------------------------------------
    ! Evaluate the out-of-plane stress component for
    ! a plane strain analysis.
    !
    ! Date: 27 June 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8) :: out_of_plane_stress
    real(kind=8),intent(in) :: sig(3)
    real(kind=8) :: el,nu
    !
    call get_young_modulus(el)
    call get_poisson_ratio(nu)
    !
    out_of_plane_stress = nu*(sig(1)+sig(2))
    !
    return
  end function out_of_plane_stress
  !****************************************************
  !****************************************************
  function eigenvalue(ss)
    !----------------------------------------------------
    ! Evaluate principal strains or stresses.
    !
    ! Date: 27 June 2022
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8) :: eigenvalue(2)
    real(kind=8),intent(in) :: ss(3)
    real(kind=8) :: b,c,d,r1,r2
    real(kind=8) :: tol = 1.0e-8_8
    !
    b = ss(1)+ss(2)
    c = ss(1)*ss(2)-ss(3)*ss(3)
    !
    d = b*b-4*c
    !
    if (d.lt.0.0) then
      if (abs(d).lt.tol) then
        d = 0.0
      else
        call error_message_stop("Complex root during eigenvalue evaluation.")
      end if
    end if
    !
    r1 = 0.5_8*(b+sqrt(d))
    r2 = 0.5_8*(b-sqrt(d))
    !
    if (r1.ge.r2) then
      eigenvalue(1) = r1
      eigenvalue(2) = r2
    else
      eigenvalue(1) = r2
      eigenvalue(2) = r1
    end if
    !
    return
  end function eigenvalue
  !****************************************************
  !****************************************************
end module fracture_2d_output
