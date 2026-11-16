!-----------------------------------------------------------------------
! MODULE: model_mod
!> @brief Estruturas de dados e variáveis globais do modelo Dual-BEM 2D.
!>
!> Simplificado para utilizar exclusivamente elementos de contorno
!> quadráticos descontínuos de 3 nós (nodes_num = 3, incid(3)).
!> Como nenhum nó funcional é compartilhado entre elementos adjacentes,
!> o cálculo do vetor normal e a aplicação das condições de contorno em
!> quinas ou encontros de trinca tornam-se diretos e livres de ambiguidade.
!-----------------------------------------------------------------------
module model_mod
  use constants_mod
  use quadrature_mod
  implicit none

  ! Identificação do modelo e parâmetros elásticos
  character(len=200), private :: model_name = ''
  character(len=200), private :: model_id   = ''
  integer, private            :: plane_state_key = 1  ! 1 = Estado Plano de Tensões (EPT), 2 = Deformações (EPD)
  real(dp), private           :: young_modulus   = 0.0_dp
  real(dp), private           :: poisson_ratio   = 0.0_dp

  !---------------------------------------------------------------------
  !> @brief Nó funcional (ponto de colocação) de elemento descontínuo.
  !---------------------------------------------------------------------
  type :: boundary_node
    integer        :: bie       = 1       ! 1 = EIC desloc. (contorno), 2 = EIC desloc. (trinca), 3 = EIC tração (trinca)
    type(point_2d) :: point               ! Coordenadas cartesianas (x, y) do nó de colocação
    integer        :: bc_dir    = 0       ! 0 = Direções cartesianas (x, y); 1 = Direções locais (normal, tangencial)
    integer        :: bc(2)     = 1       ! 0 = Deslocamento prescrito; 1 = Tração prescrita; 2 = Apoio elástico (mola)
    real(dp)       :: h1(2)     = 0.0_dp  ! Rigidez da mola: t_i = h1(i)*u_i + h2(i)
    real(dp)       :: h2(2)     = 0.0_dp  ! Termo livre da mola
    real(dp)       :: disp(2)   = 0.0_dp  ! Componentes de deslocamento {u_x, u_y}
    real(dp)       :: trac(2)   = 0.0_dp  ! Componentes de força de superfície (tração) {t_x, t_y}
    real(dp)       :: stress(3) = 0.0_dp  ! Componentes de tensão {sigma_xx, sigma_yy, tau_xy}
    real(dp)       :: strain(3) = 0.0_dp  ! Componentes de deformação {eps_xx, eps_yy, eps_xy}
  end type boundary_node

  !---------------------------------------------------------------------
  !> @brief Elemento de contorno quadrático descontínuo de 3 nós.
  !>        Os pontos geométricos extremos são bp_1 (eta=-1) e bp_2 (eta=+1),
  !>        e o ponto geométrico central (eta=0) coincide com o nó incid(2).
  !---------------------------------------------------------------------
  type :: boundary_elem
    integer        :: nodes_num = 3               ! Sempre 3 nós por elemento quadrático descontínuo
    integer        :: incid(3)  = [0, 0, 0]       ! Índices dos 3 nós de colocação em bn_list
    real(dp)       :: lamb_1    = 2.0_dp / 3.0_dp ! Posição paramétrica |eta_1| do 1º nó funcional
    real(dp)       :: lamb_2    = 2.0_dp / 3.0_dp ! Posição paramétrica |eta_3| do 3º nó funcional
    type(point_2d) :: bp_1                        ! Ponto geométrico inicial (eta = -1)
    type(point_2d) :: bp_2                        ! Ponto geométrico final   (eta = +1)
  end type boundary_elem

  !---------------------------------------------------------------------
  !> @brief Ponto interno para pós-processamento de deslocamentos e tensões.
  !---------------------------------------------------------------------
  type :: internal_node
    type(point_2d) :: point
    real(dp)       :: disp(2)   = 0.0_dp
    real(dp)       :: stress(3) = 0.0_dp
    real(dp)       :: strain(3) = 0.0_dp
  end type internal_node

  !---------------------------------------------------------------------
  !> @brief Estrutura de dados de uma trinca no método Dual-BEM.
  !---------------------------------------------------------------------
  type :: crack
    integer               :: tips_num      = 1       ! Número de pontas ativas (1 = trinca de borda, 2 = interna)
    integer               :: ext_el_num(4) = 0       ! [elem_ini_face1, elem_fim_face1, elem_ini_face2, elem_fim_face2]
    real(dp)              :: ki(2)         = 0.0_dp  ! Fator de intensidade de tensão Modo I  (ponta 1 e ponta 2)
    real(dp)              :: kii(2)        = 0.0_dp  ! Fator de intensidade de tensão Modo II (ponta 1 e ponta 2)
    real(dp)              :: theta(2)      = 0.0_dp  ! Ângulo local de deflexão da trinca (radianos)
    real(dp)              :: keq(2)        = 0.0_dp  ! Fator de intensidade de tensão equivalente
    real(dp)              :: length        = 0.0_dp  ! Comprimento total acumulado da trinca
    integer, allocatable  :: nodes_num(:,:)          ! Pares de nós coincidentes (:,1)=Face 1 (bie=2), (:,2)=Face 2 (bie=3)
  end type crack

  ! Contadores e listas globais do modelo
  integer :: bn_num = 0
  integer :: be_num = 0
  integer :: in_num = 0
  integer :: cr_num = 0

  type(boundary_node), allocatable :: bn_list(:)
  type(boundary_elem), allocatable :: be_list(:)
  type(internal_node), allocatable :: in_list(:)
  type(crack),         allocatable :: cr_list(:)

contains

  !---------------------------------------------------------------------
  !> @brief Verifica se o elemento de contorno `be` contém espacialmente
  !>        o ponto fonte `cp` (incluindo nós coincidentes na face oposta
  !>        da trinca no Dual-BEM) e retorna o índice local `lcn` (1, 2 ou 3).
  !---------------------------------------------------------------------
  pure subroutine be_contains_node(be, cp, sing, lcn)
    type(boundary_elem), intent(in) :: be
    type(boundary_node), intent(in) :: cp
    logical, intent(out)            :: sing
    integer, intent(out)            :: lcn
    type(point_2d) :: pt
    real(dp) :: tol
    integer :: i

    sing = .false.
    lcn  = 0
    tol  = distance_only_2d(be%bp_1, be%bp_2) * 1.0e-8_dp

    do i = 1, 3
      pt = bn_list(be%incid(i))%point
      if (distance_only_2d(cp%point, pt) < tol) then
        sing = .true.
        lcn  = i
        return
      end if
    end do
  end subroutine be_contains_node

  !---------------------------------------------------------------------
  !> @brief Comprimento característico do elemento quadrático de contorno.
  !---------------------------------------------------------------------
  pure function be_size(be) result(le)
    type(boundary_elem), intent(in) :: be
    real(dp) :: le
    le = distance_only_2d(be%bp_1, bn_list(be%incid(2))%point) + &
         distance_only_2d(bn_list(be%incid(2))%point, be%bp_2)
  end function be_size

  !---------------------------------------------------------------------
  !> @brief Retorna os 3 pontos geométricos (eta = -1, 0, +1) do elemento.
  !---------------------------------------------------------------------
  pure subroutine get_bound_elem_geom_points(be, be_pts)
    type(boundary_elem), intent(in) :: be
    type(point_2d), intent(out)     :: be_pts(3)

    be_pts(1) = be%bp_1
    be_pts(2) = bn_list(be%incid(2))%point
    be_pts(3) = be%bp_2
  end subroutine get_bound_elem_geom_points

  !---------------------------------------------------------------------
  !> @brief Obtém o vetor normal unitário externo em um nó de colocação.
  !>        Como todos os elementos são descontínuos, cada nó pertence a
  !>        um único elemento e possui normal única e bem definida.
  !---------------------------------------------------------------------
  function get_normal_to_a_boundary_point(inc_num) result(nv)
    integer, intent(in) :: inc_num
    real(dp) :: nv(2)
    type(point_2d) :: gi_pts(3)
    real(dp) :: eta
    integer :: elem_id, loc_id, i, j
    logical :: found

    if (inc_num < 1 .or. inc_num > bn_num) then
      call error_message_int_stop("Indice de no invalido em get_normal_to_a_boundary_point: ", inc_num)
    end if

    ! Busca direta O(1) quando a numeração segue 3*e-2, 3*e-1, 3*e
    elem_id = (inc_num - 1) / 3 + 1
    loc_id  = mod(inc_num - 1, 3) + 1
    found   = .false.

    if (elem_id >= 1 .and. elem_id <= be_num) then
      if (be_list(elem_id)%incid(loc_id) == inc_num) then
        found = .true.
      end if
    end if

    ! Caso geral (se a numeração não for sequencial)
    if (.not. found) then
      outer: do i = 1, be_num
        do j = 1, 3
          if (be_list(i)%incid(j) == inc_num) then
            elem_id = i
            loc_id  = j
            found   = .true.
            exit outer
          end if
        end do
      end do outer
    end if

    if (.not. found) then
      call error_message_int_stop("No nao encontrado na malha de elementos: ", inc_num)
    end if

    call get_bound_elem_geom_points(be_list(elem_id), gi_pts)

    select case (loc_id)
    case (1)
      eta = -be_list(elem_id)%lamb_1
    case (2)
      eta = 0.0_dp
    case (3)
      eta = +be_list(elem_id)%lamb_2
    end select

    nv = normal_1d(gi_pts, eta)
  end function get_normal_to_a_boundary_point

  !---------------------------------------------------------------------
  !> @brief Calcula o vetor de deformações {eps_xx, eps_yy, eps_xy} a partir
  !>        do vetor de tensões {sig_xx, sig_yy, tau_xy} pela Lei de Hooke.
  !---------------------------------------------------------------------
  pure subroutine get_strain_from_stress_elast(psk, nu, el, sig, eps)
    integer, intent(in)   :: psk
    real(dp), intent(in)  :: nu, el
    real(dp), intent(in)  :: sig(3)
    real(dp), intent(out) :: eps(3)
    real(dp) :: nub, eb

    if (psk == 1) then
      ! Estado Plano de Tensões (EPT)
      nub = nu
      eb  = el
    else
      ! Estado Plano de Deformações (EPD)
      nub = nu / (1.0_dp - nu)
      eb  = el / (1.0_dp - nu * nu)
    end if

    eps(1) = (sig(1) - nub * sig(2)) / eb
    eps(2) = (sig(2) - nub * sig(1)) / eb
    eps(3) = ((1.0_dp + nu) / el) * sig(3)
  end subroutine get_strain_from_stress_elast

  !---------------------------------------------------------------------
  !> @brief Constante elástica de Kolosov (kappa).
  !---------------------------------------------------------------------
  subroutine get_kappa(kappa)
    real(dp), intent(out) :: kappa
    if (plane_state_key == 1) then
      kappa = 3.0_dp - 4.0_dp * (poisson_ratio / (1.0_dp + poisson_ratio))
    else
      kappa = 3.0_dp - 4.0_dp * poisson_ratio
    end if
  end subroutine get_kappa

  !---------------------------------------------------------------------
  !> @brief Módulo de elasticidade transversal (cisalhamento) mu = E / [2(1+nu)].
  !---------------------------------------------------------------------
  subroutine get_shear_modulus(g)
    real(dp), intent(out) :: g
    g = 0.5_dp * young_modulus / (1.0_dp + poisson_ratio)
  end subroutine get_shear_modulus

  subroutine set_plane_state_key(psk)
    integer, intent(in) :: psk
    if (psk < 1 .or. psk > 2) call error_message_stop("Chave de estado plano invalida (use 1=EPT ou 2=EPD).")
    plane_state_key = psk
  end subroutine set_plane_state_key

  subroutine get_plane_state_key(psk)
    integer, intent(out) :: psk
    psk = plane_state_key
  end subroutine get_plane_state_key

  subroutine set_young_modulus(el)
    real(dp), intent(in) :: el
    if (el <= 0.0_dp) call error_message_stop("Modulo de Young deve ser positivo.")
    young_modulus = el
  end subroutine set_young_modulus

  subroutine get_young_modulus(el)
    real(dp), intent(out) :: el
    el = young_modulus
  end subroutine get_young_modulus

  subroutine set_poisson_ratio(nu)
    real(dp), intent(in) :: nu
    if (nu < 0.0_dp .or. nu >= 0.5_dp) call error_message_stop("Coeficiente de Poisson fora do intervalo [0, 0.5).")
    poisson_ratio = nu
  end subroutine set_poisson_ratio

  subroutine get_poisson_ratio(nu)
    real(dp), intent(out) :: nu
    nu = poisson_ratio
  end subroutine get_poisson_ratio

  subroutine set_model_name(m_name)
    character(len=*), intent(in) :: m_name
    model_name = trim(m_name)
  end subroutine set_model_name

  subroutine get_model_name(m_name)
    character(len=200), intent(out) :: m_name
    m_name = model_name
  end subroutine get_model_name

  subroutine set_model_id(m_id)
    character(len=*), intent(in) :: m_id
    model_id = trim(m_id)
  end subroutine set_model_id

  subroutine get_model_id(m_id)
    character(len=200), intent(out) :: m_id
    m_id = model_id
  end subroutine get_model_id

end module model_mod
