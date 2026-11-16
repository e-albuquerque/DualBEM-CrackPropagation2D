!-----------------------------------------------------------------------
! MODULE: assembler_mod
!> @brief Montagem das matrizes globais [H] e [G] do Método dos Elementos de
!>        Contorno Dual (Dual-BEM), aplicação das condições de contorno para
!>        formar o sistema algébrico [A]{x} = {f}, recuperação dos campos no
!>        contorno e avaliação nos pontos internos.
!>
!> Como todos os elementos são quadráticos descontínuos (3 nós por elemento),
!> cada grau de liberdade pertence a um único elemento, eliminando conflitos
!> de condições de contorno em quinas ou pontas de trinca.
!-----------------------------------------------------------------------
module assembler_mod
  use constants_mod
  use quadrature_mod
  use model_mod
  use kelvin_sol_mod
  use integration_mod
  implicit none

  real(dp), allocatable, private :: initial_h(:,:), initial_g(:,:)
  type(internal_node), allocatable, private :: in_list_original(:)

contains

  subroutine set_initial_h(h_in)
    real(dp), intent(in) :: h_in(:,:)
    if (allocated(initial_h)) deallocate(initial_h)
    allocate(initial_h(size(h_in, 1), size(h_in, 2)))
    initial_h = h_in
  end subroutine set_initial_h

  subroutine set_initial_g(g_in)
    real(dp), intent(in) :: g_in(:,:)
    if (allocated(initial_g)) deallocate(initial_g)
    allocate(initial_g(size(g_in, 1), size(g_in, 2)))
    initial_g = g_in
  end subroutine set_initial_g

  subroutine get_initial_h(h_out)
    real(dp), allocatable, intent(inout) :: h_out(:,:)
    h_out = initial_h
  end subroutine get_initial_h

  subroutine get_initial_g(g_out)
    real(dp), allocatable, intent(inout) :: g_out(:,:)
    g_out = initial_g
  end subroutine get_initial_g

  !---------------------------------------------------------------------
  !> @brief Monta as matrizes [H] e [G] iniciais, armazena cópias para os
  !>        incrementos futuros de propagação e aplica as condições de
  !>        contorno para obter [A]{x} = {f}.
  !---------------------------------------------------------------------
  subroutine assembler(alpha, a, f)
    real(dp), intent(in)  :: alpha
    real(dp), intent(out) :: a(2*bn_num, 2*bn_num)
    real(dp), intent(out) :: f(2*bn_num)
    real(dp), allocatable :: h(:,:), g(:,:)
    integer :: i

    ! Escala as trações prescritas pelo fator de carregamento equivalente alpha
    do i = 1, bn_num
      bn_list(i)%trac = bn_list(i)%trac * alpha
    end do

    allocate(h(2*bn_num, 2*bn_num), g(2*bn_num, 2*bn_num))
    call get_matrices_H_and_G(h, g)

    call set_initial_h(h)
    call set_initial_g(g)

    call apply_boundary_conditions(h, g, a, f)
    deallocate(h, g)
  end subroutine assembler

  !---------------------------------------------------------------------
  !> @brief Aplica as condições de contorno (deslocamento, tração ou mola,
  !>        em coordenadas cartesianas ou locais normal/tangencial) às
  !>        matrizes [H] e [G] para montar [A]{x} = {f}.
  !---------------------------------------------------------------------
  subroutine apply_boundary_conditions(h, g, a, f)
    real(dp), intent(inout) :: h(2*bn_num, 2*bn_num), g(2*bn_num, 2*bn_num)
    real(dp), intent(out)   :: a(2*bn_num, 2*bn_num), f(2*bn_num)
    real(dp) :: f_vec(2*bn_num), aux(2*bn_num), nv(2)
    type(boundary_node) :: cp
    integer :: j, k, col

    do j = 1, bn_num
      cp = bn_list(j)

      ! Rotação de colunas quando a CC é dada nas direções normal/tangencial (bc_dir /= 0)
      if (cp%bc_dir /= 0) then
        nv  = get_normal_to_a_boundary_point(j)
        col = 2*j - 1

        aux          = h(:, col)*nv(1)   + h(:, col+1)*nv(2)
        h(:, col+1)  = h(:, col+1)*nv(1) - h(:, col)*nv(2)
        h(:, col)    = aux

        aux          = g(:, col)*nv(1)   + g(:, col+1)*nv(2)
        g(:, col+1)  = g(:, col+1)*nv(1) - g(:, col)*nv(2)
        g(:, col)    = aux
      end if

      ! Troca de colunas conforme o tipo de condição de contorno em cada direção (k = 1, 2)
      do k = 1, 2
        col = 2*j + k - 2
        select case (cp%bc(k))
        case (0)
          ! Deslocamento prescrito (Dirichlet): incógnita é a tração
          f_vec(col) = cp%disp(k)
          aux        = g(:, col)
          g(:, col)  = -h(:, col)
          h(:, col)  = -aux
        case (1)
          ! Tração prescrita (Neumann): incógnita é o deslocamento
          f_vec(col) = cp%trac(k)
        case (2)
          ! Apoio elástico / mola (Robin): t_k = h1(k)*u_k + h2(k)
          f_vec(col) = cp%h2(k)
          h(:, col)  = h(:, col) - g(:, col) * cp%h1(k)
        case default
          call error_message_int_stop("Condicao de contorno indefinida no no: ", j)
        end select
      end do
    end do

    f = matmul(g, f_vec)
    a = h
  end subroutine apply_boundary_conditions

  !---------------------------------------------------------------------
  !> @brief Calcula os coeficientes das matrizes globais [H] e [G] usando:
  !>        - EIC de Deslocamento para nós do contorno externo (bie = 1) e
  !>          da primeira face da trinca (bie = 2);
  !>        - EIC Hipersingular de Tração para nós da segunda face da trinca (bie = 3).
  !---------------------------------------------------------------------
  subroutine get_matrices_H_and_G(h, g)
    real(dp), intent(out) :: h(2*bn_num, 2*bn_num), g(2*bn_num, 2*bn_num)
    real(dp) :: hi(2, 6), gi(2, 6), normal_cp(2)
    type(boundary_node) :: cp
    type(boundary_elem) :: be
    integer :: b, c, d, i, j, k, l, ii, kk, ll, bie

    h = 0.0_dp
    g = 0.0_dp

    !$omp parallel do default(shared) &
    !$omp private(i, cp, bie, ii, normal_cp, j, be, hi, gi, k, kk, c, d, b, l, ll) &
    !$omp schedule(dynamic)
    do i = 1, bn_num
      cp        = bn_list(i)
      bie       = cp%bie
      ii        = 2*i - 1
      normal_cp = get_normal_to_a_boundary_point(i)

      do j = 1, be_num
        be = be_list(j)
        hi = 0.0_dp
        gi = 0.0_dp

        if (bie == 1 .or. bie == 2) then
          call displ_bie_matrices_H_and_G(cp, be, hi, gi)
        else if (bie == 3) then
          call tract_bie_matrices_H_and_G(cp, be, normal_cp, hi, gi)
        else
          call error_message_stop("Chave BIE invalida para o no.")
        end if

        do k = 1, 3
          kk = 2 * be%incid(k) - 1
          h(ii:ii+1, kk:kk+1) = h(ii:ii+1, kk:kk+1) + hi(:, 2*k-1:2*k)
          g(ii:ii+1, kk:kk+1) = g(ii:ii+1, kk:kk+1) + gi(:, 2*k-1:2*k)
        end do
      end do

      ! Termos livres de salto de Cauchy (1/2) para contorno suave (elementos descontínuos)
      ll = 1
      if (bie == 1) then
        h(ii,   ii)   = h(ii,   ii)   + 0.5_dp
        h(ii+1, ii+1) = h(ii+1, ii+1) + 0.5_dp
      else if (bie == 2) then
        h(ii,   ii)   = h(ii,   ii)   + 0.5_dp
        h(ii+1, ii+1) = h(ii+1, ii+1) + 0.5_dp
        do c = 1, cr_num
          d = size(cr_list(c)%nodes_num, 1)
          do b = 1, d
            if (cr_list(c)%nodes_num(b, 1) == i) then
              l  = cr_list(c)%nodes_num(b, 2)
              ll = 2*l - 1
            end if
          end do
        end do
        h(ii,   ll)   = h(ii,   ll)   + 0.5_dp
        h(ii+1, ll+1) = h(ii+1, ll+1) + 0.5_dp
      else if (bie == 3) then
        g(ii,   ii)   = g(ii,   ii)   - 0.5_dp
        g(ii+1, ii+1) = g(ii+1, ii+1) - 0.5_dp
        do c = 1, cr_num
          d = size(cr_list(c)%nodes_num, 1)
          do b = 1, d
            if (cr_list(c)%nodes_num(b, 2) == i) then
              l  = cr_list(c)%nodes_num(b, 1)
              ll = 2*l - 1
            end if
          end do
        end do
        g(ii,   ll)   = g(ii,   ll)   + 0.5_dp
        g(ii+1, ll+1) = g(ii+1, ll+1) + 0.5_dp
      end if
    end do
    !$omp end parallel do
  end subroutine get_matrices_H_and_G

  !---------------------------------------------------------------------
  !> @brief Montagem incremental de [H], [G], [A] e {f} durante a propagação
  !>        da trinca (reaproveitando o bloco superior esquerdo 2*tnbn x 2*tnbn).
  !---------------------------------------------------------------------
  subroutine assembler_for_propagation(tnbe, tnbn, a, f, h_initial, g_initial)
    integer, intent(in)   :: tnbe, tnbn
    real(dp), intent(out) :: a(2*bn_num, 2*bn_num), f(2*bn_num)
    real(dp), intent(out) :: h_initial(2*bn_num, 2*bn_num), g_initial(2*bn_num, 2*bn_num)
    real(dp), allocatable :: h(:,:), g(:,:)

    allocate(h(2*bn_num, 2*bn_num), g(2*bn_num, 2*bn_num))
    call get_matrices_H_and_G_for_propagation(tnbe, tnbn, h, g)

    h_initial = h
    g_initial = g

    call apply_boundary_conditions(h, g, a, f)
    deallocate(h, g)
  end subroutine assembler_for_propagation

  subroutine get_matrices_H_and_G_for_propagation(tnbe, tnbn, h, g)
    integer, intent(in)   :: tnbe, tnbn
    real(dp), intent(out) :: h(2*bn_num, 2*bn_num), g(2*bn_num, 2*bn_num)
    real(dp) :: hi(2, 6), gi(2, 6), normal_cp(2)
    type(boundary_node) :: cp
    type(boundary_elem) :: be
    integer :: b, c, d, i, j, k, l, ii, kk, ll, bie

    h = 0.0_dp
    g = 0.0_dp
    h(1:2*tnbn, 1:2*tnbn) = initial_h(1:2*tnbn, 1:2*tnbn)
    g(1:2*tnbn, 1:2*tnbn) = initial_g(1:2*tnbn, 1:2*tnbn)

    ! 1) Novos nós de colocação (tnbn+1 .. bn_num) integrados sobre todos os elementos (1 .. be_num)
    !$omp parallel do default(shared) &
    !$omp private(i, cp, bie, ii, normal_cp, j, be, hi, gi, k, kk, c, d, b, l, ll) &
    !$omp schedule(dynamic)
    do i = tnbn + 1, bn_num
      cp        = bn_list(i)
      bie       = cp%bie
      ii        = 2*i - 1
      normal_cp = get_normal_to_a_boundary_point(i)

      do j = 1, be_num
        be = be_list(j)
        hi = 0.0_dp
        gi = 0.0_dp

        if (bie == 2) then
          call displ_bie_matrices_H_and_G(cp, be, hi, gi)
        else if (bie == 3) then
          call tract_bie_matrices_H_and_G(cp, be, normal_cp, hi, gi)
        else
          call error_message_stop("Chave BIE invalida em propagacao.")
        end if

        do k = 1, 3
          kk = 2 * be%incid(k) - 1
          h(ii:ii+1, kk:kk+1) = h(ii:ii+1, kk:kk+1) + hi(:, 2*k-1:2*k)
          g(ii:ii+1, kk:kk+1) = g(ii:ii+1, kk:kk+1) + gi(:, 2*k-1:2*k)
        end do
      end do

      ll = 1
      if (bie == 2) then
        h(ii,   ii)   = h(ii,   ii)   + 0.5_dp
        h(ii+1, ii+1) = h(ii+1, ii+1) + 0.5_dp
        do c = 1, cr_num
          d = size(cr_list(c)%nodes_num, 1)
          do b = 1, d
            if (cr_list(c)%nodes_num(b, 1) == i) then
              l  = cr_list(c)%nodes_num(b, 2)
              ll = 2*l - 1
            end if
          end do
        end do
        h(ii,   ll)   = h(ii,   ll)   + 0.5_dp
        h(ii+1, ll+1) = h(ii+1, ll+1) + 0.5_dp
      else if (bie == 3) then
        g(ii,   ii)   = g(ii,   ii)   - 0.5_dp
        g(ii+1, ii+1) = g(ii+1, ii+1) - 0.5_dp
        do c = 1, cr_num
          d = size(cr_list(c)%nodes_num, 1)
          do b = 1, d
            if (cr_list(c)%nodes_num(b, 2) == i) then
              l  = cr_list(c)%nodes_num(b, 1)
              ll = 2*l - 1
            end if
          end do
        end do
        g(ii,   ll)   = g(ii,   ll)   + 0.5_dp
        g(ii+1, ll+1) = g(ii+1, ll+1) + 0.5_dp
      end if
    end do
    !$omp end parallel do

    ! 2) Nós antigos (1 .. tnbn) integrados sobre os novos elementos da trinca (tnbe+1 .. be_num)
    !$omp parallel do default(shared) &
    !$omp private(i, cp, bie, ii, normal_cp, j, be, hi, gi, k, kk) &
    !$omp schedule(dynamic)
    do i = 1, tnbn
      cp        = bn_list(i)
      bie       = cp%bie
      ii        = 2*i - 1
      normal_cp = get_normal_to_a_boundary_point(i)

      do j = tnbe + 1, be_num
        be = be_list(j)
        hi = 0.0_dp
        gi = 0.0_dp

        if (bie == 1 .or. bie == 2) then
          call displ_bie_matrices_H_and_G(cp, be, hi, gi)
        else if (bie == 3) then
          call tract_bie_matrices_H_and_G(cp, be, normal_cp, hi, gi)
        end if

        do k = 1, 3
          kk = 2 * be%incid(k) - 1
          h(ii:ii+1, kk:kk+1) = h(ii:ii+1, kk:kk+1) + hi(:, 2*k-1:2*k)
          g(ii:ii+1, kk:kk+1) = g(ii:ii+1, kk:kk+1) + gi(:, 2*k-1:2*k)
        end do
      end do
    end do
    !$omp end parallel do
  end subroutine get_matrices_H_and_G_for_propagation

  !---------------------------------------------------------------------
  !> @brief Submatrizes locais (2x6) da EIC de Deslocamento para um par (cp, be).
  !---------------------------------------------------------------------
  subroutine displ_bie_matrices_H_and_G(cp, be, hi, gi)
    type(boundary_node), intent(in) :: cp
    type(boundary_elem), intent(in) :: be
    real(dp), intent(inout)         :: hi(2, 6), gi(2, 6)
    real(dp) :: mu, nu
    integer  :: psk, lcn
    logical  :: sing

    call get_plane_state_key(psk)
    call get_shear_modulus(mu)
    call get_poisson_ratio(nu)

    call be_contains_node(be, cp, sing, lcn)
    if (sing) then
      call displ_bie_sing_integr(cp%point, be, lcn, hi)
      call displ_bie_weakly_sing_integr(psk, mu, nu, lcn, be, gi)
    else
      call displ_bie_reg_integr(psk, mu, nu, cp%point, be, hi, gi)
    end if
  end subroutine displ_bie_matrices_H_and_G

  !---------------------------------------------------------------------
  !> @brief Submatrizes locais (2x6) da EIC de Tração para um par (cp, be).
  !---------------------------------------------------------------------
  subroutine tract_bie_matrices_H_and_G(cp, be, normal_cp, hi, gi)
    type(boundary_node), intent(in) :: cp
    type(boundary_elem), intent(in) :: be
    real(dp), intent(in)            :: normal_cp(2)
    real(dp), intent(inout)         :: hi(2, 6), gi(2, 6)
    real(dp) :: mu, nu
    integer  :: psk, lcn
    logical  :: sing

    call get_plane_state_key(psk)
    call get_shear_modulus(mu)
    call get_poisson_ratio(nu)

    call be_contains_node(be, cp, sing, lcn)
    if (sing) then
      call traction_bie_sing_intgr(cp%point, be, lcn, normal_cp, hi, gi)
    else
      call traction_bie_reg_integr(psk, mu, nu, cp%point, be, normal_cp, hi, gi)
    end if
  end subroutine tract_bie_matrices_H_and_G

  !=====================================================================
  ! RECUPERAÇÃO DOS CAMPOS DE CONTORNO E PONTOS INTERNOS
  !=====================================================================

  !---------------------------------------------------------------------
  !> @brief Distribui o vetor solução {u} nos campos disp e trac de cada nó.
  !---------------------------------------------------------------------
  subroutine set_results_on_boundary_nodes(u)
    real(dp), intent(in) :: u(:)
    real(dp) :: nv(2), aux
    integer  :: i, j, k

    if (size(u) /= 2*bn_num) then
      call error_message_stop("Tamanho do vetor solucao incompativel com 2*bn_num.")
    end if

    do i = 1, bn_num
      if (bn_list(i)%bc_dir /= 0) then
        nv  = get_normal_to_a_boundary_point(i)
        aux = bn_list(i)%disp(1)*nv(1) - bn_list(i)%disp(2)*nv(2)
        bn_list(i)%disp(2) = bn_list(i)%disp(1)*nv(2) + bn_list(i)%disp(2)*nv(1)
        bn_list(i)%disp(1) = aux

        aux = bn_list(i)%trac(1)*nv(1) - bn_list(i)%trac(2)*nv(2)
        bn_list(i)%trac(2) = bn_list(i)%trac(1)*nv(2) + bn_list(i)%trac(2)*nv(1)
        bn_list(i)%trac(1) = aux

        do k = 1, 2
          if (k == 2) then
            aux   = nv(1)
            nv(1) = -nv(2)
            nv(2) = aux
          end if
          select case (bn_list(i)%bc(k))
          case (0)
            do j = 1, 2
              bn_list(i)%trac(j) = bn_list(i)%trac(j) + u(2*i+k-2)*nv(j)
            end do
          case (1)
            do j = 1, 2
              bn_list(i)%disp(j) = bn_list(i)%disp(j) + u(2*i+k-2)*nv(j)
            end do
          case (2)
            aux = u(2*i+k-2)*bn_list(i)%h1(k) + bn_list(i)%h2(k)
            do j = 1, 2
              bn_list(i)%disp(j) = bn_list(i)%disp(j) + u(2*i+k-2)*nv(j)
              bn_list(i)%trac(j) = bn_list(i)%trac(j) + aux*nv(j)
            end do
          end select
        end do
      else
        do k = 1, 2
          select case (bn_list(i)%bc(k))
          case (0)
            bn_list(i)%trac(k) = u(2*i+k-2)
          case (1)
            bn_list(i)%disp(k) = u(2*i+k-2)
          case (2)
            bn_list(i)%disp(k) = u(2*i+k-2)
            bn_list(i)%trac(k) = u(2*i+k-2)*bn_list(i)%h1(k) + bn_list(i)%h2(k)
          end select
        end do
      end if
    end do
  end subroutine set_results_on_boundary_nodes

  !---------------------------------------------------------------------
  !> @brief Avalia tensões e deformações nos nós do contorno (por recuperação
  !>        de derivadas de funções de forma ou por equações integrais singulares).
  !---------------------------------------------------------------------
  subroutine evaluate_stress_strain_on_boundary_nodes(bound_recovery)
    logical, intent(in) :: bound_recovery

    if (bound_recovery) then
      call stress_strain_using_boundary_recovery()
    else
      call stress_strain_using_singular_integration()
    end if
  end subroutine evaluate_stress_strain_on_boundary_nodes

  subroutine stress_strain_using_boundary_recovery()
    real(dp)       :: grad_u(2, 2), l1(2), l2(2), aa, bb, u(2, 3), t(2)
    real(dp)       :: mu, nu, nub, lb, c1, c2, c3, cc, eta, jac, dni(3)
    type(point_2d) :: gi_pts(3)
    integer        :: psk, inum, i, j, k, l, m, n

    call get_plane_state_key(psk)
    call get_shear_modulus(mu)
    call get_poisson_ratio(nu)

    nub = merge(nu / (1.0_dp + nu), nu, psk == 1)
    lb  = 2.0_dp * mu * nub / (1.0_dp - 2.0_dp * nub)
    c1  = 1.0_dp / mu
    c2  = lb / (2.0_dp * mu + lb)
    c3  = 1.0_dp / (2.0_dp * mu + lb)

    do l = 1, be_num
      call get_bound_elem_geom_points(be_list(l), gi_pts)
      do k = 1, 3
        u(:, k) = bn_list(be_list(l)%incid(k))%disp
      end do

      do n = 1, 3
        inum = be_list(l)%incid(n)
        if (n == 1) then
          eta = -be_list(l)%lamb_1
        else if (n == 3) then
          eta = +be_list(l)%lamb_2
        else
          eta = 0.0_dp
        end if

        call deriv_intrp_func_1d(dni, eta, be_list(l)%lamb_1, be_list(l)%lamb_2)
        call jac_only_1d(jac, gi_pts, eta)

        grad_u = 0.0_dp
        l1     = get_normal_to_a_boundary_point(inum)
        l2(1)  = -l1(2); l2(2) = l1(1)
        t      = bn_list(inum)%trac

        do k = 1, 2
          do m = 1, 2
            do j = 1, 2
              bb = c1*l2(k)*l1(m)*l2(j) + c3*l1(k)*l1(m)*l1(j)
              grad_u(k, m) = grad_u(k, m) + bb * t(j)
              do i = 1, 3
                aa = (l1(k)*l2(m)*l1(j) - l2(k)*l1(m)*l1(j) + l2(k)*l2(m)*l2(j) &
                     - c2*l1(k)*l1(m)*l2(j)) * (1.0_dp / jac) * dni(i)
                grad_u(k, m) = grad_u(k, m) + aa * u(j, i)
              end do
            end do
          end do
        end do

        bn_list(inum)%strain(1) = grad_u(1, 1)
        bn_list(inum)%strain(2) = grad_u(2, 2)
        bn_list(inum)%strain(3) = 0.5_dp * (grad_u(1, 2) + grad_u(2, 1))

        cc = lb * (grad_u(1, 1) + grad_u(2, 2))
        bn_list(inum)%stress(1) = 2.0_dp * mu * grad_u(1, 1) + cc
        bn_list(inum)%stress(2) = 2.0_dp * mu * grad_u(2, 2) + cc
        bn_list(inum)%stress(3) = mu * (grad_u(1, 2) + grad_u(2, 1))
      end do
    end do
  end subroutine stress_strain_using_boundary_recovery

  subroutine stress_strain_using_singular_integration()
    real(dp) :: eps(3), sig(3), hs(3, 6), gs(3, 6), u(6), t(6)
    real(dp) :: mu, nu, el
    integer  :: psk, i, j, k, lcn
    logical  :: sing
    type(boundary_node) :: cp
    type(boundary_elem) :: be

    call get_plane_state_key(psk)
    call get_young_modulus(el)
    call get_shear_modulus(mu)
    call get_poisson_ratio(nu)

    !$omp parallel do default(shared) &
    !$omp private(i, cp, eps, sig, j, be, hs, gs, sing, lcn, k, u, t) &
    !$omp schedule(dynamic)
    do i = 1, bn_num
      eps = 0.0_dp
      sig = 0.0_dp
      cp  = bn_list(i)

      do j = 1, be_num
        be = be_list(j)
        call be_contains_node(be, cp, sing, lcn)
        if (sing) then
          call stress_sing_intgr(cp%point, be, lcn, hs, gs)
        else
          call stress_reg_integr(psk, mu, nu, cp%point, be, hs, gs)
        end if

        do k = 1, 3
          u(2*k-1:2*k) = bn_list(be%incid(k))%disp
          t(2*k-1:2*k) = bn_list(be%incid(k))%trac
        end do
        sig = sig - matmul(hs, u) + matmul(gs, t)
      end do

      sig = 2.0_dp * sig
      bn_list(i)%stress = sig
      call get_strain_from_stress_elast(psk, nu, el, sig, eps)
      bn_list(i)%strain = eps
    end do
    !$omp end parallel do
  end subroutine stress_strain_using_singular_integration

  !---------------------------------------------------------------------
  !> @brief Avalia deslocamentos, tensões e deformações nos pontos internos.
  !---------------------------------------------------------------------
  subroutine evaluate_internal_points()
    real(dp) :: u(2), eps(3), sig(3), h(2, 6), g(2, 6), hs(3, 6), gs(3, 6)
    real(dp) :: u_temp(6), t_temp(6), mu, nu, el
    integer  :: psk, i, j, k
    type(internal_node) :: cp
    type(boundary_elem) :: be

    call get_plane_state_key(psk)
    call get_young_modulus(el)
    call get_shear_modulus(mu)
    call get_poisson_ratio(nu)

    !$omp parallel do default(shared) &
    !$omp private(i, cp, u, eps, sig, j, be, h, g, hs, gs, k, u_temp, t_temp) &
    !$omp schedule(dynamic)
    do i = 1, in_num
      u   = 0.0_dp
      eps = 0.0_dp
      sig = 0.0_dp
      cp  = in_list(i)

      do j = 1, be_num
        be = be_list(j)
        call reg_integr_int_point_to_be(psk, mu, nu, cp%point, be, h, g, hs, gs)
        do k = 1, 3
          u_temp(2*k-1:2*k) = bn_list(be%incid(k))%disp
          t_temp(2*k-1:2*k) = bn_list(be%incid(k))%trac
        end do
        u   = u   - matmul(h,  u_temp) + matmul(g,  t_temp)
        sig = sig - matmul(hs, u_temp) + matmul(gs, t_temp)
      end do

      in_list(i)%disp   = u
      in_list(i)%stress = sig
      call get_strain_from_stress_elast(psk, nu, el, sig, eps)
      in_list(i)%strain = eps
    end do
    !$omp end parallel do
  end subroutine evaluate_internal_points

  !---------------------------------------------------------------------
  !> @brief Cria temporariamente pontos internos circulares ao redor das
  !>        pontas das trincas para exportação detalhada no arquivo XML.
  !---------------------------------------------------------------------
  subroutine create_internal_points(nip, initial_node_num)
    integer, intent(in) :: nip, initial_node_num
    type(point_2d) :: p1, p2, pts(nip)
    type(internal_node), allocatable :: in_list_temp(:)
    real(dp) :: eta(nip), ome(nip), rad(3), theta, theta_initial
    integer  :: i, m, n, j, k, l, aux1, initial_node_pos
    type(boundary_elem) :: last_elem

    if (allocated(in_list_original)) deallocate(in_list_original)
    allocate(in_list_original(in_num))
    if (in_num > 0) in_list_original = in_list

    call gauss_legendre_points_1d(nip, eta, ome)

    do k = 1, cr_num
      do j = 1, cr_list(k)%tips_num
        if (j == 1) then
          aux1             = size(cr_list(k)%nodes_num, 1) - initial_node_num + 1
          initial_node_pos = cr_list(k)%nodes_num(aux1, 2)
          last_elem        = be_list(cr_list(k)%ext_el_num(2))
        else
          initial_node_pos = cr_list(k)%nodes_num(initial_node_num, 1)
          last_elem        = be_list(cr_list(k)%ext_el_num(4))
        end if

        p1            = last_elem%bp_2
        p2            = bn_list(initial_node_pos)%point
        rad           = distance_2d(p1, p2)
        theta_initial = atan2(rad(3), rad(2))

        do i = 1, nip
          theta    = theta_initial + pi * (eta(i) + 1.0_dp)
          pts(i)%x = p1%x + rad(1) * cos(theta)
          pts(i)%y = p1%y + rad(1) * sin(theta)
        end do

        m = size(in_list)
        allocate(in_list_temp(m))
        if (m > 0) in_list_temp = in_list
        deallocate(in_list)
        n = m + nip
        allocate(in_list(n))
        if (m > 0) in_list(1:m) = in_list_temp
        l = 1
        do i = m + 1, n
          in_list(i)%point = pts(l)
          l = l + 1
        end do
        deallocate(in_list_temp)
      end do
    end do
    in_num = size(in_list)
  end subroutine create_internal_points

  !---------------------------------------------------------------------
  !> @brief Restaura a lista original de pontos internos após a escrita XML.
  !---------------------------------------------------------------------
  subroutine original_internal_points()
    integer :: n_orig
    n_orig = size(in_list_original)
    if (allocated(in_list)) deallocate(in_list)
    allocate(in_list(n_orig))
    if (n_orig > 0) in_list = in_list_original
    deallocate(in_list_original)
    in_num = n_orig
  end subroutine original_internal_points

end module assembler_mod
