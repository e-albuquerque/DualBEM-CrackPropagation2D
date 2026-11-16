program Fracture_2D_v3_0
  !----------------------------------------------------
  ! A Program for two-dimensional fracture analysis
  ! of isotropic media by the boundary element method.
  !
  ! Version 3.0
  !
  ! Date: October 2022
  ! Last update: ---
  !----------------------------------------------------
  use fracture_2d_input
  use fracture_2d_assembler
  use fracture_2d_set_bound_results
  use fracture_2d_internal_points
  use fracture_2d_output
  use fracture_2d_model
  use fracture_2d_stress_intens_factor
  use fracture_2d_crack_propagation
  use fracture_2d_create_internal_points
  use matrix_algebra
  !
  implicit none
  !
  character(len=256) :: input_file_name
  !
  real(kind=8),allocatable :: a(:,:),f(:) !The system of equations
  real :: t1,t2,t_sf,t_ss,t_ip,t_cr
  integer :: n
  logical :: bound_recovery,variable_ampl
  integer :: steps
  real(kind=8) :: aq,amax
  !
  t_sf = 0.0; t_ss = 0.0; t_ip = 0.0; t_cr = 0.0
  !
  print 100
  100 format(78("*"))
  !
  print 101
  101 format("Fracture 2D - Version 2.0 - October 2022",/)
  print 102
  102 format("A Dual-BEM based program for isotropic materials and single region.",/)
  !
  !-------------------------------------------
  ! Reading data file and mounting the model:
  !-------------------------------------------
  !
  if (command_argument_count().ge.1) then
    call get_command_argument(1,input_file_name)
  else
    print 110
    110 format("Input File Name:")
    read(*,'(A)') input_file_name
  end if
  input_file_name = adjustl(input_file_name)
  !
  call input_from_text_file(trim(input_file_name),steps,variable_ampl)
  !
  print*
  !
  !-----------------------------------------
  ! Printing output file with initial data:
  !-----------------------------------------
  !
  print*,"Printing file with initial parameters... "
  !
  call initial_data_output()
  !
  print*,"                                         ...done!"
  print*
  !
  !
  !----------------------
  ! Counting the cycles:
  !----------------------
  !
  if (variable_ampl) then
    print*,"Counting the cycles... "
    !
    call rainflow(aq,amax)
    call set_parameter_R(0.0_8)
    !
    print*,"                         ...done!"
    print*
    else
        call get_aq(aq)
        amax = 1.0_8/aq
        call set_parameter_R(0.0_8)
  end if
  !
  !-----------------------------------
  ! Mounting the system of equations:
  !-----------------------------------
  !
  print*,"Mounting the system of equations... "
  !
  call cpu_time(t1)
  !
  allocate(a(2*bn_num,2*bn_num))
  allocate(f(2*bn_num))
  !
  call assembler(aq,a,f)
  !
  !
  call cpu_time(t2)
  t_sf = t2 - t1
  !
  print*,"                                    ...done!"
  print*
  !
  !----------------------------------
  ! Solving the system of equations:
  !----------------------------------
  !
  print*,"Solving the system of equations..."
  !
  call cpu_time(t1)
  !
  call solve_linear_system(a,f)
  !
  call cpu_time(t2)
  t_ss = t2 - t1
  !
  print*,"                                   ...done!"
  print*
  !
  deallocate(a)
  !
  !----------------------------------------
  ! Setting boundary results in the model:
  !----------------------------------------
  !
  print*,"Setting results on boundary nodes... "
  !
  call set_results_on_boundary_nodes(f)
  !
  bound_recovery = .false.
  call evaluate_stress_strain_on_boundary_nodes(bound_recovery)
  !
  print*,"                                     ...done!"
  print*
  !
  deallocate(f)
  !
  !------------------------------------
  ! Evaluating internal points fields:
  !------------------------------------
  !
  call create_internal_points(20,2)
  !
  if (in_num.gt.0) then
    !
    print*,"Evaluating internal points... "
    !
    call cpu_time(t1)
    !
    call evaluate_internal_points()
    !
    call cpu_time(t2)
    t_ip = t2 - t1
    !
    print*,"                              ...done!"
    print*
    !
  end if
  !
  !------------------------------------
  ! Evaluating cracks:
  !------------------------------------
  !
  if (cr_num.gt.0) then
    !
    print*,"Evaluating cracks... "
    !
    call cpu_time(t1)
    !
    call stress_intensity_factor_evaluation()
    call propagation_angle()
    call equivalent_k()
    call cpu_time(t2)
    t_cr = t2 - t1
    print*,"                     ...done!"
    print*
    !
    !
    !------------------------
    ! Printing output files:
    !------------------------
    !
    ! A file with the results per step:
    !
    call results_output(t_sf,t_ss,t_ip,t_cr,0)
    !
    ! Files with the results per crack:
    !
    call crack_geometry_output(0)
    call results_per_crack_output(0,amax) !amax
    call fatigue_results_per_crack_output(0) !aq
    !
    ! A file with the results for INSANE Graphical interface:
    !
    call output_file_for_insane(0) !amax
    !
    call original_internal_points()
    !
    !
    print*,"Propagating cracks..."
    do n = 1,steps
        !
        call crack_propagation(n)
        !
        !
        !------------------------------------
        ! Evaluating internal points fields:
        !------------------------------------
        !
        call create_internal_points(20,2)
        !
        if (in_num.gt.0) then
            !
            call cpu_time(t1)
            !
            call evaluate_internal_points()
            !
            call cpu_time(t2)
            t_ip = t2 - t1
            !
        end if
        !
        if ((mod(n,2) == 0)) then
            call number_of_cycles()
            call fatigue_results_per_crack_output(n)
        end if
        !
        !------------------------
        ! Printing output files:
        !------------------------
        !
        ! A file with the results:
        !
        call results_output(t_sf,t_ss,t_ip,t_cr,n)
        !
        ! Files with the results per crack:
        !
        call results_per_crack_output(n,amax)
        call crack_geometry_output(n)
        !
        ! A file with the results for INSANE Graphical interface:
        !
        call output_file_for_insane(n)
        !
        call original_internal_points()
        !
    end do
    print*,"                      ...done!"
  end if
  print 100
end
