!-----------------------------------------------------------------------
! MODULE: io_mod
!> @brief Rotinas de leitura do arquivo de dados pré-processado (.dat) e
!>        geração dos arquivos de saída (.txt e .xml para o INSANE / Python).
!>
!> O arquivo de entrada (.dat) é gerado pelo script Python `preproc.py` a
!> partir da malha do Gmsh (.msh) lida com `meshio`, contendo apenas
!> elementos de contorno quadráticos descontínuos de 3 nós.
!-----------------------------------------------------------------------
module io_mod
  use constants_mod
  use quadrature_mod
  use model_mod
  use integration_mod
  use sif_mod
  use fatigue_prop_mod
  implicit none

contains

  !---------------------------------------------------------------------
  !> @brief Lê a próxima linha útil do arquivo (ignorando linhas em branco
  !>        e linhas de comentário iniciadas por '#', '*' ou '!').
  !---------------------------------------------------------------------
  subroutine next_data_line(unit_num, line)
    integer, intent(in)           :: unit_num
    character(len=*), intent(out) :: line
    integer :: ios
    character(len=512) :: buf

    do
      read(unit_num, '(A)', iostat=ios) buf
      if (ios /= 0) then
        call error_message_stop("Fim inesperado do arquivo de entrada (.dat).")
      end if
      buf = adjustl(buf)
      if (len_trim(buf) == 0) cycle
      if (buf(1:1) == '#' .or. buf(1:1) == '*' .or. buf(1:1) == '!') cycle
      line = trim(buf)
      return
    end do
  end subroutine next_data_line

  !---------------------------------------------------------------------
  !> @brief Lê o arquivo `.dat` pré-processado (com elementos quadráticos
  !>        descontínuos) e inicializa todas as estruturas do modelo Dual-BEM.
  !---------------------------------------------------------------------
  subroutine input_from_text_file(file_name, steps, variable_ampl)
    character(len=*), intent(in) :: file_name
    integer, intent(out)         :: steps
    logical, intent(out)         :: variable_ampl

    character(len=512) :: line
    character(len=200) :: analysis_name, analysis_id
    integer  :: psk, npow, k_key_in, method_key_in, k_eff_in, vari_ampl_key
    integer  :: nbe_j, nbn_j, np_j, init_node_j, nd_correl
    integer  :: i, k, id_dummy, n_pairs, tips_num, n_alpha
    real(dp) :: el, nu, lb1_j, lb2_j
    real(dp) :: x1, y1, x2, y2, lb1, lb2, val1_a, val1_b, val2_a, val2_b
    real(dp) :: da0, r_da, da_max, R_in, C0_in, m_in, gama_in, kth_in, kic_in, a_val
    integer  :: bie_in, bc_dir_in, bck1, bck2, inc(3), ext_el(4), nf1, nf2
    real(dp) :: cr_len
    integer, parameter :: u_in = 10

    write(*, '(1X,A,A)') "Lendo arquivo de dados pre-processado: ", trim(file_name)
    open(unit=u_in, file=trim(file_name), status='old', action='read')

    ! 1. Cabeçalho e propriedades do material
    call next_data_line(u_in, line)
    analysis_name = trim(line)
    call set_model_name(analysis_name)

    call next_data_line(u_in, line)
    read(line, *) analysis_id
    call set_model_id(analysis_id)

    call next_data_line(u_in, line)
    read(line, *) el, nu, psk, npow
    call set_young_modulus(el)
    call set_poisson_ratio(nu)
    call set_plane_state_key(psk)
    call set_npow(npow)

    ! 2. Dimensões da malha descontínua
    call next_data_line(u_in, line)
    read(line, *) be_num, bn_num, cr_num, in_num

    if (allocated(be_list)) deallocate(be_list)
    if (allocated(bn_list)) deallocate(bn_list)
    if (allocated(cr_list)) deallocate(cr_list)
    if (allocated(in_list)) deallocate(in_list)

    allocate(be_list(be_num))
    allocate(bn_list(bn_num))
    allocate(cr_list(cr_num))
    allocate(in_list(in_num))

    ! 3. Tabela de elementos de contorno quadráticos descontínuos
    do i = 1, be_num
      call next_data_line(u_in, line)
      read(line, *) id_dummy, x1, y1, x2, y2, lb1, lb2, inc(1), inc(2), inc(3)
      be_list(i)%bp_1%x    = x1
      be_list(i)%bp_1%y    = y1
      be_list(i)%bp_2%x    = x2
      be_list(i)%bp_2%y    = y2
      be_list(i)%nodes_num = 3
      be_list(i)%lamb_1    = lb1
      be_list(i)%lamb_2    = lb2
      be_list(i)%incid     = inc
    end do

    ! 4. Tabela de nós funcionais de colocação e condições de contorno
    do i = 1, bn_num
      call next_data_line(u_in, line)
      read(line, *) id_dummy, x1, y1, bie_in, bc_dir_in, bck1, val1_a, val1_b, bck2, val2_a, val2_b
      bn_list(i)%point%x = x1
      bn_list(i)%point%y = y1
      bn_list(i)%bie     = bie_in
      bn_list(i)%bc_dir  = bc_dir_in
      bn_list(i)%bc(1)   = bck1
      bn_list(i)%bc(2)   = bck2
      bn_list(i)%disp    = 0.0_dp
      bn_list(i)%trac    = 0.0_dp
      bn_list(i)%h1      = 0.0_dp
      bn_list(i)%h2      = 0.0_dp
      bn_list(i)%strain  = 0.0_dp
      bn_list(i)%stress  = 0.0_dp

      select case (bck1)
      case (0)
        bn_list(i)%disp(1) = val1_a
      case (1)
        bn_list(i)%trac(1) = val1_a
      case (2)
        bn_list(i)%h1(1)   = val1_a
        bn_list(i)%h2(1)   = val1_b
      case default
        call error_message_int_stop("Chave de CC invalida na direcao 1 do no ", i)
      end select

      select case (bck2)
      case (0)
        bn_list(i)%disp(2) = val2_a
      case (1)
        bn_list(i)%trac(2) = val2_a
      case (2)
        bn_list(i)%h1(2)   = val2_a
        bn_list(i)%h2(2)   = val2_b
      case default
        call error_message_int_stop("Chave de CC invalida na direcao 2 do no ", i)
      end select
    end do

    ! 5. Topologia das trincas (Dual-BEM)
    do i = 1, cr_num
      call next_data_line(u_in, line)
      read(line, *) id_dummy, tips_num, n_pairs, ext_el(1), ext_el(2), ext_el(3), ext_el(4), cr_len
      cr_list(i)%tips_num   = tips_num
      cr_list(i)%ext_el_num = ext_el
      cr_list(i)%length     = cr_len
      cr_list(i)%ki         = 0.0_dp
      cr_list(i)%kii        = 0.0_dp
      cr_list(i)%theta      = 0.0_dp
      cr_list(i)%keq        = 0.0_dp
      allocate(cr_list(i)%nodes_num(n_pairs, 2))
      do k = 1, n_pairs
        call next_data_line(u_in, line)
        read(line, *) id_dummy, nf1, nf2
        cr_list(i)%nodes_num(k, 1) = nf1
        cr_list(i)%nodes_num(k, 2) = nf2
      end do
    end do

    ! 6. Pontos internos
    do i = 1, in_num
      call next_data_line(u_in, line)
      read(line, *) id_dummy, x1, y1
      in_list(i)%point%x = x1
      in_list(i)%point%y = y1
      in_list(i)%disp    = 0.0_dp
      in_list(i)%stress  = 0.0_dp
      in_list(i)%strain  = 0.0_dp
    end do

    ! 7. Parâmetros de cálculo dos Fatores de Intensidade de Tensão (FIT)
    call next_data_line(u_in, line)
    read(line, *) k_key_in, method_key_in
    call set_K_key(k_key_in)
    call set_method_key(method_key_in)

    call next_data_line(u_in, line)
    if (k_key_in == 1) then
      if (method_key_in == 7) then
        read(line, *) nd_correl
        call set_node(nd_correl)
      end if
    else if (k_key_in == 2) then
      if (method_key_in == 1 .or. method_key_in == 2) then
        read(line, *) np_j, init_node_j
        call set_number_of_points(np_j)
        call set_initial_node_num(init_node_j)
      else if (method_key_in == 3) then
        read(line, *) nbe_j, nbn_j, lb1_j, lb2_j, np_j, init_node_j
        call set_element_information(nbe_j, nbn_j, lb1_j, lb2_j)
        call set_number_of_points(np_j)
        call set_initial_node_num(init_node_j)
      end if
    end if

    ! 8. Parâmetros de propagação de trinca e fadiga (Paris-Walker)
    call next_data_line(u_in, line)
    read(line, *) k_eff_in
    call set_k_eff_key(k_eff_in)

    call next_data_line(u_in, line)
    read(line, *) da0, r_da, da_max, steps
    if (steps > 0 .and. mod(steps, 2) == 1) steps = steps + 1
    call set_initial_delta_a(da0)
    call set_r_Da(r_da)
    call set_max_delta_a(da_max)

    call next_data_line(u_in, line)
    read(line, *) R_in, C0_in, m_in, gama_in, kth_in, kic_in
    call set_parameter_R(R_in)
    call set_parameter_C(C0_in)
    call set_parameter_m(m_in)
    call set_parameter_gamma(gama_in)
    call set_K_th(kth_in)
    call set_K_IC(kic_in)

    ! 9. Carregamento de amplitude variável (Rainflow)
    call next_data_line(u_in, line)
    read(line, *) vari_ampl_key
    if (vari_ampl_key == 0) then
      variable_ampl = .false.
    else if (vari_ampl_key == 1) then
      variable_ampl = .true.
      call next_data_line(u_in, line)
      read(line, *) n_alpha
      call allocate_alpha(n_alpha)
      do i = 1, n_alpha
        call next_data_line(u_in, line)
        read(line, *) a_val
        call set_alpha(i, a_val)
      end do
    else
      call error_message_int_stop("Chave vari_ampl_key invalida: ", vari_ampl_key)
    end if

    close(u_in)

    write(*, '(1X,A,I0)') "Numero de elementos de contorno (descontinuos) = ", be_num
    write(*, '(1X,A,I0)') "Numero de nos funcionais de colocacao          = ", bn_num
    if (in_num > 0) write(*, '(1X,A,I0)') "Numero de pontos internos                      = ", in_num
    write(*, *)
  end subroutine input_from_text_file

  !=====================================================================
  ! ARQUIVOS DE SAÍDA (.txt e .xml)
  !=====================================================================

  subroutine initial_data_output()
    character(len=200) :: file_name, analysis_id, analysis_name
    character(len=6)   :: bc1, bc2
    character(len=2)   :: dir1, dir2
    real(dp) :: pv1, pv2, el, nu
    integer  :: psk, i, j

    call get_plane_state_key(psk)
    call get_young_modulus(el)
    call get_poisson_ratio(nu)
    call get_model_id(analysis_id)
    call get_model_name(analysis_name)

    file_name = trim(analysis_id) // "_input_data.txt"
    open(10, file=file_name, status='replace', action='write')

    write(10, 100)
    write(10, 101)
    write(10, 102) trim(analysis_name)
    if (psk == 1) then
      write(10, 104)
    else
      write(10, 105)
    end if
    write(10, 106) el
    write(10, 107) nu
    write(10, 111) be_num
    write(10, 112) bn_num
    if (in_num > 0) write(10, 113) in_num

    write(10, 120)
    do i = 1, be_num
      write(10, 121) i, be_list(i)%bp_1%x, be_list(i)%bp_1%y, be_list(i)%bp_2%x, be_list(i)%bp_2%y
    end do

    write(10, 122)
    do i = 1, be_num
      write(10, 123) i, be_list(i)%lamb_1, be_list(i)%lamb_2, (be_list(i)%incid(j), j = 1, 3)
    end do

    write(10, 130)
    do i = 1, bn_num
      if (bn_list(i)%bc_dir == 0) then
        dir1 = "_x"; dir2 = "_y"
      else
        dir1 = "_n"; dir2 = "_t"
      end if

      select case (bn_list(i)%bc(1))
      case (0)
        bc1 = "disp" // dir1; pv1 = bn_list(i)%disp(1)
      case (1)
        bc1 = "trac" // dir1; pv1 = bn_list(i)%trac(1)
      case (2)
        bc1 = "sprg" // dir1; pv1 = bn_list(i)%h1(1)
      end select

      select case (bn_list(i)%bc(2))
      case (0)
        bc2 = "disp" // dir2; pv2 = bn_list(i)%disp(2)
      case (1)
        bc2 = "trac" // dir2; pv2 = bn_list(i)%trac(2)
      case (2)
        bc2 = "sprg" // dir2; pv2 = bn_list(i)%h1(2)
      end select

      write(10, 131) i, bn_list(i)%point%x, bn_list(i)%point%y, bc1, pv1, bc2, pv2
    end do

    if (in_num > 0) then
      write(10, 140)
      do i = 1, in_num
        write(10, 141) i, in_list(i)%point%x, in_list(i)%point%y
      end do
    end if

    write(10, 100)
    close(10)

100 format(82("*"))
101 format("Fracture 2D",//,"A Dual-BEM based program - Discontinuous Elements Edition",/)
102 format("Analysis name: ",(a),/)
104 format(1x,"Plane stress problem.",/)
105 format(1x,"Plane strain problem.",/)
106 format(1x,"Elasticity modulus = ",e20.12)
107 format(1x,"Poisson ratio = ",f0.8,/)
111 format(1x,"Number of boundary elements = ",i0)
112 format(1x,"Number of boundary points = ",i0)
113 format(1x,"Number of internal points = ",i0)
120 format(/,"Boundary Elements:",//,"Elem.",5x,"X1",12x,"Y1",12x,"X2",12x,"Y2")
121 format(i5,4(2x,e12.5))
122 format(/,"Elem.",5x,"Lambda 1",6x,"Lambda 2",6x,"Incidences")
123 format(i5,2(2x,e12.5),3(2x,i7))
130 format(/,"Boundary Points:",//,"Point",5x,"X",13x,"Y",10x,"B.C. 1",4x,"Value",8x,"B.C. 2",4x,"Value")
131 format(i5,2(2x,e12.5),2(2x,(a),1x,e14.7))
140 format(/,"Internal Points:",//,"Point",5x,"X",13x,"Y")
141 format(i5,2(2x,e12.5))
  end subroutine initial_data_output

  subroutine results_output(t_sf, t_ss, t_ip, t_cr, n)
    real, intent(in)    :: t_sf, t_ss, t_ip, t_cr
    integer, intent(in) :: n
    character(len=200)  :: file_name, analysis_id, analysis_name
    character(len=3)    :: num
    integer  :: psk, i, j, k, K_key_cur
    real(dp) :: el, nu, s33, vms, ps(2)

    call get_plane_state_key(psk)
    call get_young_modulus(el)
    call get_poisson_ratio(nu)
    call get_model_id(analysis_id)
    call get_model_name(analysis_name)

    write(num, '(i3.3)') n
    file_name = trim(analysis_id) // "_results_step_" // trim(num) // ".txt"
    open(10, file=file_name, status='replace', action='write')

    write(10, 100)
    write(10, 101)
    write(10, 102) trim(analysis_name)
    if (psk == 1) then
      write(10, 104)
    else
      write(10, 105)
    end if

    call get_K_key(K_key_cur)
    call get_method_key(k)
    if (K_key_cur == 1) then
      write(10, 143) k
    else if (K_key_cur == 2) then
      write(10, 144) k
    end if

    write(10, 106) el
    write(10, 107) nu
    write(10, 108) be_num
    write(10, 109) bn_num
    if (in_num > 0) write(10, 110) in_num
    write(10, 113) t_sf, t_ss, t_ip, t_cr

    write(10, 120)
    do i = 1, bn_num
      write(10, 121) i, bn_list(i)%point%x, bn_list(i)%point%y, &
                     bn_list(i)%disp(1), bn_list(i)%disp(2), &
                     bn_list(i)%trac(1), bn_list(i)%trac(2)
    end do

    if (psk == 1) then
      write(10, 132)
      do i = 1, bn_num
        s33 = out_of_plane_strain(bn_list(i)%stress)
        write(10, 121) i, bn_list(i)%point%x, bn_list(i)%point%y, &
                       bn_list(i)%strain(1), bn_list(i)%strain(2), bn_list(i)%strain(3), s33
      end do
      write(10, 133)
      do i = 1, bn_num
        write(10, 121) i, bn_list(i)%point%x, bn_list(i)%point%y, &
                       bn_list(i)%stress(1), bn_list(i)%stress(2), bn_list(i)%stress(3)
      end do
    else
      write(10, 131)
      do i = 1, bn_num
        write(10, 121) i, bn_list(i)%point%x, bn_list(i)%point%y, &
                       bn_list(i)%strain(1), bn_list(i)%strain(2), bn_list(i)%strain(3)
      end do
      write(10, 134)
      do i = 1, bn_num
        s33 = out_of_plane_stress(bn_list(i)%stress)
        write(10, 121) i, bn_list(i)%point%x, bn_list(i)%point%y, &
                       bn_list(i)%stress(1), bn_list(i)%stress(2), bn_list(i)%stress(3), s33
      end do
    end if

    write(10, 135)
    do i = 1, bn_num
      ps = eigenvalue(bn_list(i)%strain)
      write(10, 121) i, bn_list(i)%point%x, bn_list(i)%point%y, ps(1), ps(2)
    end do

    write(10, 136)
    do i = 1, bn_num
      ps  = eigenvalue(bn_list(i)%stress)
      vms = von_mises_stress(bn_list(i)%stress)
      write(10, 121) i, bn_list(i)%point%x, bn_list(i)%point%y, ps(1), ps(2), vms
    end do

    if (in_num > 0) then
      write(10, 130)
      do i = 1, in_num
        write(10, 121) i, in_list(i)%point%x, in_list(i)%point%y, in_list(i)%disp(1), in_list(i)%disp(2)
      end do

      if (psk == 1) then
        write(10, 132)
        do i = 1, in_num
          s33 = out_of_plane_strain(in_list(i)%stress)
          write(10, 121) i, in_list(i)%point%x, in_list(i)%point%y, &
                         in_list(i)%strain(1), in_list(i)%strain(2), in_list(i)%strain(3), s33
        end do
        write(10, 133)
        do i = 1, in_num
          write(10, 121) i, in_list(i)%point%x, in_list(i)%point%y, &
                         in_list(i)%stress(1), in_list(i)%stress(2), in_list(i)%stress(3)
        end do
      else
        write(10, 131)
        do i = 1, in_num
          write(10, 121) i, in_list(i)%point%x, in_list(i)%point%y, &
                         in_list(i)%strain(1), in_list(i)%strain(2), in_list(i)%strain(3)
        end do
        write(10, 134)
        do i = 1, in_num
          s33 = out_of_plane_stress(in_list(i)%stress)
          write(10, 121) i, in_list(i)%point%x, in_list(i)%point%y, &
                         in_list(i)%stress(1), in_list(i)%stress(2), in_list(i)%stress(3), s33
        end do
      end if

      write(10, 135)
      do i = 1, in_num
        ps = eigenvalue(in_list(i)%strain)
        write(10, 121) i, in_list(i)%point%x, in_list(i)%point%y, ps(1), ps(2)
      end do

      write(10, 136)
      do i = 1, in_num
        ps  = eigenvalue(in_list(i)%stress)
        vms = von_mises_stress(in_list(i)%stress)
        write(10, 121) i, in_list(i)%point%x, in_list(i)%point%y, ps(1), ps(2), vms
      end do
    end if

    write(10, 140)
    write(10, 141)
    do i = 1, cr_num
      do j = 1, cr_list(i)%tips_num
        write(10, 142) i, j, cr_list(i)%ki(j), cr_list(i)%kii(j), cr_list(i)%theta(j)*180.0_dp/pi, cr_list(i)%keq(j)
      end do
    end do
    write(10, 100)
    close(10)

100 format(100("*"))
101 format("Fracture 2D",//,"A Dual-BEM based program - Discontinuous Elements Edition",/)
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
120 format(/,"Boundary Points:",//,"Point",5x,"X",13x,"Y",13x,"Disp. X",9x,"Disp. Y",9x,"Trac. X",9x,"Trac. Y")
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
  end subroutine results_output

  subroutine results_per_crack_output(n, amax)
    integer, intent(in)  :: n
    real(dp), intent(in) :: amax
    character(len=200) :: analysis_id, file_name
    character(len=3)   :: num
    real(dp) :: length, ki, kii, ki_2, kii_2, keq, keq_2, K_IC_val
    integer  :: i, u_cr

    call get_model_id(analysis_id)
    call get_K_IC(K_IC_val)

    do i = 1, cr_num
      write(num, '(i3.3)') i
      file_name = trim(analysis_id) // "_results_crack_" // trim(num) // ".txt"
      u_cr      = 20 + i

      if (n == 0) then
        open(u_cr, file=file_name, status='replace', action='write')
        if (cr_list(i)%tips_num == 1) then
          write(u_cr, 101)
        else
          write(u_cr, 102)
        end if
      else
        open(u_cr, file=file_name, status='old', position='append', action='write')
      end if

      length = cr_list(i)%length
      ki     = amax * cr_list(i)%ki(1)
      kii    = amax * cr_list(i)%kii(1)
      ki_2   = amax * cr_list(i)%ki(2)
      kii_2  = amax * cr_list(i)%kii(2)
      keq    = amax * cr_list(i)%keq(1)
      keq_2  = amax * cr_list(i)%keq(2)

      if (n > 0 .and. (keq > K_IC_val .or. keq_2 > K_IC_val)) then
        if (cr_list(i)%tips_num == 1) then
          write(u_cr, 105) n, length, ki, kii, keq
        else
          write(u_cr, 106) n, length, ki, kii, keq, ki_2, kii_2, keq_2
        end if
      else
        if (cr_list(i)%tips_num == 1) then
          write(u_cr, 103) n, length, ki, kii, keq
        else
          write(u_cr, 104) n, length, ki, kii, keq, ki_2, kii_2, keq_2
        end if
      end if
      close(u_cr)
    end do

101 format(2x,"Step",5x,"Length",11x,"K_I_1",11x,"K_II_1",11x,"Keq_1")
102 format(2x,"Step",5x,"Length",11x,"K_I_1",11x,"K_II_1",11x,"Keq_1",11x,"K_I_2",10x,"K_II_2",11x,"Keq_2")
103 format(i5,4(2x,e14.7))
104 format(i5,7(2x,e14.7))
105 format(i5,4(2x,e14.7),5x,"**** Keq greater than K_IC!****")
106 format(i5,7(2x,e14.7),5x,"**** Keq greater than K_IC!****")
  end subroutine results_per_crack_output

  subroutine crack_geometry_output(n)
    integer, intent(in) :: n
    character(len=200) :: analysis_id, file_name
    character(len=3)   :: num
    real(dp) :: x_initial, y_initial, x_final, y_final, length
    integer  :: i, u_cr

    call get_model_id(analysis_id)

    do i = 1, cr_num
      write(num, '(i3.3)') i
      file_name = trim(analysis_id) // "_geometry_crack_" // trim(num) // ".txt"
      u_cr      = 30 + i

      if (n == 0) then
        open(u_cr, file=file_name, status='replace', action='write')
        write(u_cr, 101)
      else
        open(u_cr, file=file_name, status='old', position='append', action='write')
      end if

      x_initial = be_list(cr_list(i)%ext_el_num(4))%bp_2%x
      y_initial = be_list(cr_list(i)%ext_el_num(4))%bp_2%y
      x_final   = be_list(cr_list(i)%ext_el_num(2))%bp_2%x
      y_final   = be_list(cr_list(i)%ext_el_num(2))%bp_2%y
      length    = cr_list(i)%length

      write(u_cr, 102) n, x_initial, y_initial, x_final, y_final, length
      close(u_cr)
    end do

101 format(2x,"Step",3x,"Initial_x",5x,"Initial_y",6x,"Final_x",7x,"Final_y",8x,"Length")
102 format(i5,5(e14.5))
  end subroutine crack_geometry_output

  subroutine fatigue_results_per_crack_output(n)
    integer, intent(in) :: n
    character(len=200) :: analysis_id, file_name
    character(len=3)   :: num
    real(dp) :: length, N_1_step, N_2_step, N_1_total, N_2_total
    integer  :: i, u_cr

    call get_model_id(analysis_id)

    do i = 1, cr_num
      write(num, '(i3.3)') i
      file_name = trim(analysis_id) // "_fatigue_results_crack_" // trim(num) // ".txt"
      u_cr      = 40 + i

      if (n == 0) then
        open(u_cr, file=file_name, status='replace', action='write')
        if (cr_list(i)%tips_num == 1) then
          write(u_cr, 101)
        else
          write(u_cr, 102)
        end if
        length    = cr_list(i)%length
        N_1_step  = 0.0_dp; N_1_total = 0.0_dp
        N_2_step  = 0.0_dp; N_2_total = 0.0_dp
      else
        open(u_cr, file=file_name, status='old', position='append', action='write')
        length = cr_list(i)%length
        call get_number_of_cycles_per_step(i, 1, N_1_step)
        call get_number_of_cycles_total(i, 1, N_1_total)
        call get_number_of_cycles_per_step(i, 2, N_2_step)
        call get_number_of_cycles_total(i, 2, N_2_total)
      end if

      if (cr_list(i)%tips_num == 1) then
        write(u_cr, 103) length, N_1_step, N_1_total
      else
        write(u_cr, 104) length, N_1_step, N_1_total, N_2_step, N_2_total
      end if
      close(u_cr)
    end do

101 format(5x,"Length",10x,"N_1_step",8x,"N_1_total")
102 format(5x,"Length",10x,"N_1_step",8x,"N_1_total",8x,"N_2_step",8x,"N_2_total")
103 format(3(2x,e14.7))
104 format(5(2x,e14.7))
  end subroutine fatigue_results_per_crack_output

  subroutine cycles_output(cycles, alpha_max, R_vec)
    integer, intent(in)  :: cycles
    real(dp), intent(in) :: alpha_max(:), R_vec(:)
    character(len=200) :: analysis_id, file_name
    integer :: i

    call get_model_id(analysis_id)
    file_name = trim(analysis_id) // "_cycles.txt"
    open(10, file=file_name, status='replace', action='write')
    write(10, 101)
    do i = 1, cycles
      write(10, 102) i, alpha_max(i), R_vec(i)
    end do
    close(10)

101 format("Cycle",5x,"alpha_max",8x,"R")
102 format(i5,2(2x,e14.7))
  end subroutine cycles_output

  !---------------------------------------------------------------------
  !> @brief Gera o arquivo XML de cada passo para visualização no INSANE
  !>        e no script Python `plot_results.py`.
  !---------------------------------------------------------------------
  subroutine output_file_for_insane(n)
    integer, intent(in) :: n
    type(boundary_node), allocatable :: geom_nodes_list(:)
    integer  :: pseudo_elem_num, psk, i, j, k
    real(dp) :: ni(3), el, nu, eps_33, sig_33, vms, ps(2)
    character(len=200) :: file_name, analysis_id
    character(len=38)  :: header_1 = "<?xml version='1.1' encoding='UTF-8'?>"
    character(len=40)  :: header_2 = "<Insane xmlns='http://www.dees.ufmg.br' "
    character(len=54)  :: header_3 = "xmlns:xsi='http://www.w3.org/2001/XMLSchema-instance' "
    character(len=56)  :: header_4 = "xsi:schemaLocation='http://www.dees.ufmg.br insane.xsd'>"
    character(len=14)  :: an_type
    character(len=5)   :: bound_cond(2)
    character(len=32)  :: linear_bound_elem    = "ParametricElement.BemElement.BL2"
    character(len=32)  :: quadratic_bound_elem = "ParametricElement.BemElement.BL3"
    character(len=3)   :: num

    pseudo_elem_num = 0
    do i = 1, be_num
      if (be_list(i)%lamb_1 < 1.0_dp) pseudo_elem_num = pseudo_elem_num + 1
      if (be_list(i)%lamb_2 < 1.0_dp) pseudo_elem_num = pseudo_elem_num + 1
    end do

    allocate(geom_nodes_list(pseudo_elem_num))
    do i = 1, pseudo_elem_num
      geom_nodes_list(i)%disp   = 0.0_dp
      geom_nodes_list(i)%trac   = 0.0_dp
      geom_nodes_list(i)%strain = 0.0_dp
      geom_nodes_list(i)%stress = 0.0_dp
    end do

    call get_model_id(analysis_id)
    write(num, '(i3.3)') n
    file_name = trim(analysis_id) // "-STEP-" // trim(num) // ".xml"
    open(10, file=file_name, status='replace', action='write')

    write(10, 100) header_1
    write(10, 101) header_2, header_3, header_4
    write(10, 110)
    write(10, 111)
    write(10, 112)
    write(10, 120)
    write(10, 130)

    call get_plane_state_key(psk)
    an_type = merge("BemPlaneStress", "BemPlaneStrain", psk == 1)
    write(10, 140) an_type

    write(10, 150)
    write(10, 151)
    write(10, 152)
    write(10, 153)
    write(10, 154)
    write(10, 155)
    write(10, 156)

    call get_young_modulus(el)
    call get_poisson_ratio(nu)
    write(10, 160)
    write(10, 161)
    write(10, 162) el
    write(10, 163) nu
    write(10, 164)
    write(10, 165)

    write(10, 210)
    k = 0
    do i = 1, be_num
      ! Pseudo-elemento antes do 1º nó descontínuo
      if (be_list(i)%lamb_1 < 1.0_dp) then
        k = k + 1
        geom_nodes_list(k)%point  = be_list(i)%bp_1
        geom_nodes_list(k)%bc_dir = bn_list(be_list(i)%incid(1))%bc_dir
        geom_nodes_list(k)%bc     = bn_list(be_list(i)%incid(1))%bc
        call intrp_func_1d(ni, -1.0_dp, be_list(i)%lamb_1, be_list(i)%lamb_2)
        do j = 1, 3
          geom_nodes_list(k)%disp   = geom_nodes_list(k)%disp   + ni(j)*bn_list(be_list(i)%incid(j))%disp
          geom_nodes_list(k)%trac   = geom_nodes_list(k)%trac   + ni(j)*bn_list(be_list(i)%incid(j))%trac
          geom_nodes_list(k)%strain = geom_nodes_list(k)%strain + ni(j)*bn_list(be_list(i)%incid(j))%strain
          geom_nodes_list(k)%stress = geom_nodes_list(k)%stress + ni(j)*bn_list(be_list(i)%incid(j))%stress
        end do
        write(10, 224) linear_bound_elem
        write(10, 222) k + bn_num + in_num, be_list(i)%incid(1)
        write(10, 230) an_type
        write(10, 221)
      end if

      ! Elemento quadrático entre os nós funcionais 1, 2 e 3
      write(10, 220) quadratic_bound_elem, i
      write(10, 223) be_list(i)%incid(1), be_list(i)%incid(3), be_list(i)%incid(2)
      write(10, 230) an_type
      write(10, 221)

      ! Pseudo-elemento após o 3º nó descontínuo
      if (be_list(i)%lamb_2 < 1.0_dp) then
        k = k + 1
        geom_nodes_list(k)%point  = be_list(i)%bp_2
        geom_nodes_list(k)%bc_dir = bn_list(be_list(i)%incid(3))%bc_dir
        geom_nodes_list(k)%bc     = bn_list(be_list(i)%incid(3))%bc
        call intrp_func_1d(ni, +1.0_dp, be_list(i)%lamb_1, be_list(i)%lamb_2)
        do j = 1, 3
          geom_nodes_list(k)%disp   = geom_nodes_list(k)%disp   + ni(j)*bn_list(be_list(i)%incid(j))%disp
          geom_nodes_list(k)%trac   = geom_nodes_list(k)%trac   + ni(j)*bn_list(be_list(i)%incid(j))%trac
          geom_nodes_list(k)%strain = geom_nodes_list(k)%strain + ni(j)*bn_list(be_list(i)%incid(j))%strain
          geom_nodes_list(k)%stress = geom_nodes_list(k)%stress + ni(j)*bn_list(be_list(i)%incid(j))%stress
        end do
        write(10, 224) linear_bound_elem
        write(10, 222) be_list(i)%incid(3), k + bn_num + in_num
        write(10, 230) an_type
        write(10, 221)
      end if
    end do
    write(10, 211)

    write(10, 170)
    do i = 1, bn_num
      do j = 1, 2
        bound_cond(j) = merge("true ", "false", bn_list(i)%bc(j) == 0)
      end do
      write(10, 180) i
      write(10, 181) bn_list(i)%point%x, bn_list(i)%point%y
      write(10, 182)
      write(10, 183)
      write(10, 184)
      write(10, 185)
      write(10, 186) trim(bound_cond(1)), trim(bound_cond(2))
      write(10, 187) bn_list(i)%disp(1), bn_list(i)%disp(2)
      write(10, 188) bn_list(i)%trac(1), bn_list(i)%trac(2)
      if (psk == 1) then
        sig_33 = 0.0_dp
        eps_33 = out_of_plane_strain(bn_list(i)%stress)
      else
        eps_33 = 0.0_dp
        sig_33 = out_of_plane_stress(bn_list(i)%stress)
      end if
      write(10, 189) bn_list(i)%strain(1), bn_list(i)%strain(2), bn_list(i)%strain(3), eps_33
      write(10, 190) bn_list(i)%stress(1), bn_list(i)%stress(2), bn_list(i)%stress(3), sig_33
      write(10, 193)
      ps = eigenvalue(bn_list(i)%strain)
      write(10, 194) ps(1), ps(2)
      ps = eigenvalue(bn_list(i)%stress)
      write(10, 195) ps(1), ps(2)
      vms = von_mises_stress(bn_list(i)%stress)
      write(10, 196) vms
      write(10, 191)
      write(10, 192)
    end do

    do i = 1, in_num
      j = i + bn_num
      write(10, 180) j
      write(10, 181) in_list(i)%point%x, in_list(i)%point%y
      write(10, 182)
      write(10, 200)
      write(10, 184)
      write(10, 185)
      write(10, 187) in_list(i)%disp(1), in_list(i)%disp(2)
      if (psk == 1) then
        sig_33 = 0.0_dp
        eps_33 = out_of_plane_strain(in_list(i)%stress)
      else
        eps_33 = 0.0_dp
        sig_33 = out_of_plane_stress(in_list(i)%stress)
      end if
      write(10, 189) in_list(i)%strain(1), in_list(i)%strain(2), in_list(i)%strain(3), eps_33
      write(10, 190) in_list(i)%stress(1), in_list(i)%stress(2), in_list(i)%stress(3), sig_33
      write(10, 193)
      ps = eigenvalue(in_list(i)%strain)
      write(10, 194) ps(1), ps(2)
      ps = eigenvalue(in_list(i)%stress)
      write(10, 195) ps(1), ps(2)
      vms = von_mises_stress(in_list(i)%stress)
      write(10, 196) vms
      write(10, 191)
      write(10, 192)
    end do

    do i = 1, pseudo_elem_num
      do j = 1, 2
        bound_cond(j) = merge("true ", "false", geom_nodes_list(i)%bc(j) == 0)
      end do
      write(10, 180) i + bn_num + in_num
      write(10, 181) geom_nodes_list(i)%point%x, geom_nodes_list(i)%point%y
      write(10, 182)
      write(10, 183)
      write(10, 184)
      write(10, 185)
      write(10, 186) trim(bound_cond(1)), trim(bound_cond(2))
      write(10, 187) geom_nodes_list(i)%disp(1), geom_nodes_list(i)%disp(2)
      write(10, 188) geom_nodes_list(i)%trac(1), geom_nodes_list(i)%trac(2)
      if (psk == 1) then
        sig_33 = 0.0_dp
        eps_33 = out_of_plane_strain(geom_nodes_list(i)%stress)
      else
        eps_33 = 0.0_dp
        sig_33 = out_of_plane_stress(geom_nodes_list(i)%stress)
      end if
      write(10, 189) geom_nodes_list(i)%strain(1), geom_nodes_list(i)%strain(2), geom_nodes_list(i)%strain(3), eps_33
      write(10, 190) geom_nodes_list(i)%stress(1), geom_nodes_list(i)%stress(2), geom_nodes_list(i)%stress(3), sig_33
      write(10, 193)
      ps = eigenvalue(geom_nodes_list(i)%strain)
      write(10, 194) ps(1), ps(2)
      ps = eigenvalue(geom_nodes_list(i)%stress)
      write(10, 195) ps(1), ps(2)
      vms = von_mises_stress(geom_nodes_list(i)%stress)
      write(10, 196) vms
      write(10, 191)
      write(10, 192)
    end do

    write(10, 171)
    write(10, 121)
    write(10, 900)
    close(10)
    deallocate(geom_nodes_list)

100 format((a))
101 format((a),(a),(a))
110 format(2x,"<Solution class='BemSteadyState'>")
111 format(4x,"<SolverType>1</SolverType>")
112 format(2x,"</Solution>")
120 format(2x,"<Model class='BemModel'>")
121 format(2x,"</Model>")
130 format(4x,"<ProblemDriver>BemProblemDriver</ProblemDriver>")
140 format(4x,"<GlobalAnalysisModel>",(a),"</GlobalAnalysisModel>")
150 format(4x,"<DegenerationList>")
151 format(6x,"<Degeneration class='PrescribedDegeneration' label='PrescribedDegeneration'>")
152 format(8x,"<Height>1.0000000</Height>")
153 format(8x,"<CSMaterial>LinearElasticIsotropic</CSMaterial>")
154 format(8x,"<Thickness>1.0000000</Thickness>")
155 format(6x,"</Degeneration>")
156 format(4x,"</DegenerationList>")
160 format(4x,"<MaterialList>")
161 format(6x,"<Material class='LinearElasticIsotropic' label='LinearElasticIsotropic'>")
162 format(8x,"<Elasticity>",e20.10,"</Elasticity>")
163 format(8x,"<Poisson>",e20.10,"</Poisson>")
164 format(6x,"</Material>")
165 format(4x,"</MaterialList>")
170 format(4x,"<BemNodeList>")
171 format(4x,"</BemNodeList>")
180 format(6x,"<BemNode label='",i0,"'>")
181 format(8x,"<Coord>",e20.10,1x,e20.10,1x,"0.0</Coord>")
182 format(8x,"<BemNodeValues>")
183 format(10x,"<Internal>false</Internal>")
184 format(10x,"<DOFLabels>Dx Dy</DOFLabels>")
185 format(10x,"<StrainStressLabels>XX YY XY ZZ</StrainStressLabels>")
186 format(10x,"<Restraints>",(a),1x,(a),"</Restraints>")
187 format(10x,"<Displ>",e20.10,1x,e20.10,"</Displ>")
188 format(10x,"<Tract>",e20.10,1x,e20.10,"</Tract>")
189 format(10x,"<Strain>",e20.10,1x,e20.10,1x,e20.10,1x,e20.10,"</Strain>")
190 format(10x,"<Stress>",e20.10,1x,e20.10,1x,e20.10,1x,e20.10,"</Stress>")
191 format(8x,"</BemNodeValues>")
192 format(6x,"</BemNode>")
193 format(10x,"<PrStrainStressLabels>1 2</PrStrainStressLabels>")
194 format(10x,"<PrStrain>",e20.10,1x,e20.10,"</PrStrain>")
195 format(10x,"<PrStress>",e20.10,1x,e20.10,"</PrStress>")
196 format(10x,"<VonMisesStress>",e20.10,"</VonMisesStress>")
200 format(10x,"<Internal>true</Internal>")
210 format(4x,"<ElementList>")
211 format(4x,"</ElementList>")
220 format(6x,"<Element class='",(a),"' label='E",i0,"'>")
221 format(6x,"</Element>")
222 format(8x,"<Incidence>",i0,1x,i0,"</Incidence>")
223 format(8x,"<Incidence>",i0,1x,i0,1x,i0,"</Incidence>")
224 format(6x,"<Element class='",(a),"' label=' '>")
230 format(8x,"<AnalysisModel>",(a),"</AnalysisModel>")
900 format("</Insane>")
  end subroutine output_file_for_insane

  function von_mises_stress(sig) result(vms)
    real(dp), intent(in) :: sig(3)
    real(dp) :: vms, s33, s1, s2, s3
    integer  :: psk

    call get_plane_state_key(psk)
    s33 = merge(0.0_dp, out_of_plane_stress(sig), psk == 1)
    s1  = (sig(1) - sig(2))**2
    s2  = (sig(2) - s33)**2
    s3  = (s33 - sig(1))**2
    vms = sqrt(0.5_dp * (s1 + s2 + s3) + 3.0_dp * sig(3)*sig(3))
  end function von_mises_stress

  function out_of_plane_strain(sig) result(eps33)
    real(dp), intent(in) :: sig(3)
    real(dp) :: eps33, el, nu
    call get_young_modulus(el)
    call get_poisson_ratio(nu)
    eps33 = -nu * (sig(1) + sig(2)) / el
  end function out_of_plane_strain

  function out_of_plane_stress(sig) result(sig33)
    real(dp), intent(in) :: sig(3)
    real(dp) :: sig33, nu
    call get_poisson_ratio(nu)
    sig33 = nu * (sig(1) + sig(2))
  end function out_of_plane_stress

  function eigenvalue(ss) result(eig)
    real(dp), intent(in) :: ss(3)
    real(dp) :: eig(2), b, c, d, r1, r2
    b = ss(1) + ss(2)
    c = ss(1)*ss(2) - ss(3)*ss(3)
    d = max(0.0_dp, b*b - 4.0_dp*c)
    r1 = 0.5_dp * (b + sqrt(d))
    r2 = 0.5_dp * (b - sqrt(d))
    eig(1) = max(r1, r2)
    eig(2) = min(r1, r2)
  end function eigenvalue

end module io_mod
