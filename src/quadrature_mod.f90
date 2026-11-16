!-----------------------------------------------------------------------
! MODULE: quadrature_mod
!> @brief Pontos e pesos de quadratura (Gauss-Legendre e Gauss-Laguerre) e
!>        funções de forma para elementos de contorno quadráticos descontínuos.
!>
!> Neste programa, TODOS os elementos de contorno são quadráticos (3 nós):
!> - A geometria do elemento é interpolada pelos 3 pontos geométricos em
!>   eta = -1 (bp_1), eta = 0 (nó central) e eta = +1 (bp_2).
!> - Os campos mecânicos (deslocamentos e trações) são interpolados pelos
!>   3 nós funcionais de colocação descontínuos em:
!>   eta_1 = -lamb, eta_2 = 0, eta_3 = +lamb  (tipicamente lamb = 2/3).
!-----------------------------------------------------------------------
module quadrature_mod
  use constants_mod
  implicit none

contains

  !---------------------------------------------------------------------
  !> @brief Funções de forma quadráticas (3 nós) para interpolação:
  !>        - Para geometria contínua: usar lamb_1 = 1.0, lamb_2 = 1.0
  !>        - Para elemento descontínuo simétrico: usar lamb_1 = lamb, lamb_2 = lamb
  !---------------------------------------------------------------------
  pure subroutine intrp_func_1d(ni, eta, lamb_1, lamb_2)
    real(dp), intent(out) :: ni(3)
    real(dp), intent(in)  :: eta, lamb_1, lamb_2
    real(dp) :: a

    a = 1.0_dp / (lamb_1 + lamb_2)
    ni(1) = a * eta * (eta - lamb_2) / lamb_1
    ni(2) = 1.0_dp - (eta * (eta + lamb_1 - lamb_2) / (lamb_1 * lamb_2))
    ni(3) = a * eta * (eta + lamb_1) / lamb_2
  end subroutine intrp_func_1d

  !---------------------------------------------------------------------
  !> @brief Primeira derivada das funções de forma quadráticas em relação a eta.
  !---------------------------------------------------------------------
  pure subroutine deriv_intrp_func_1d(dni, eta, lamb_1, lamb_2)
    real(dp), intent(out) :: dni(3)
    real(dp), intent(in)  :: eta, lamb_1, lamb_2
    real(dp) :: a

    a = 1.0_dp / (lamb_1 + lamb_2)
    dni(1) = a * (2.0_dp * eta - lamb_2) / lamb_1
    dni(2) = -(2.0_dp * eta + lamb_1 - lamb_2) / (lamb_1 * lamb_2)
    dni(3) = a * (2.0_dp * eta + lamb_1) / lamb_2
  end subroutine deriv_intrp_func_1d

  !---------------------------------------------------------------------
  !> @brief Segunda derivada das funções de forma quadráticas em relação a eta.
  !---------------------------------------------------------------------
  pure subroutine second_deriv_intrp_func_1d(d2ni, lamb_1, lamb_2)
    real(dp), intent(out) :: d2ni(3)
    real(dp), intent(in)  :: lamb_1, lamb_2
    real(dp) :: a

    a = 1.0_dp / (lamb_1 + lamb_2)
    d2ni(1) =  2.0_dp * a / lamb_1
    d2ni(2) =  2.0_dp / (lamb_1 * lamb_2)
    d2ni(3) =  2.0_dp * a / lamb_2
  end subroutine second_deriv_intrp_func_1d

  !---------------------------------------------------------------------
  !> @brief Calcula o Jacobiano e o vetor normal unitário externo em eta
  !>        para um elemento quadrático definido por gi_pts(1:3).
  !---------------------------------------------------------------------
  pure subroutine jac_normal_1d(jac, v_normal, gi_pts, eta)
    real(dp), intent(out)      :: jac
    real(dp), intent(out)      :: v_normal(2)
    type(point_2d), intent(in) :: gi_pts(3)
    real(dp), intent(in)       :: eta
    real(dp) :: dni(3), v_tang(2)

    call deriv_intrp_func_1d(dni, eta, 1.0_dp, 1.0_dp)
    v_tang(1) = dni(1)*gi_pts(1)%x + dni(2)*gi_pts(2)%x + dni(3)*gi_pts(3)%x
    v_tang(2) = dni(1)*gi_pts(1)%y + dni(2)*gi_pts(2)%y + dni(3)*gi_pts(3)%y

    jac = sqrt(v_tang(1)*v_tang(1) + v_tang(2)*v_tang(2))
    v_normal(1) =  v_tang(2) / jac
    v_normal(2) = -v_tang(1) / jac
  end subroutine jac_normal_1d

  !---------------------------------------------------------------------
  !> @brief Calcula apenas o Jacobiano geométrico em eta para gi_pts(1:3).
  !---------------------------------------------------------------------
  pure subroutine jac_only_1d(jac, gi_pts, eta)
    real(dp), intent(out)      :: jac
    type(point_2d), intent(in) :: gi_pts(3)
    real(dp), intent(in)       :: eta
    real(dp) :: dni(3), v_tang(2)

    call deriv_intrp_func_1d(dni, eta, 1.0_dp, 1.0_dp)
    v_tang(1) = dni(1)*gi_pts(1)%x + dni(2)*gi_pts(2)%x + dni(3)*gi_pts(3)%x
    v_tang(2) = dni(1)*gi_pts(1)%y + dni(2)*gi_pts(2)%y + dni(3)*gi_pts(3)%y
    jac = sqrt(v_tang(1)*v_tang(1) + v_tang(2)*v_tang(2))
  end subroutine jac_only_1d

  !---------------------------------------------------------------------
  !> @brief Calcula o vetor normal unitário externo em eta para gi_pts(1:3).
  !---------------------------------------------------------------------
  pure function normal_1d(gi_pts, eta) result(v_normal)
    type(point_2d), intent(in) :: gi_pts(3)
    real(dp), intent(in)       :: eta
    real(dp) :: v_normal(2), jac

    call jac_normal_1d(jac, v_normal, gi_pts, eta)
  end function normal_1d

  !---------------------------------------------------------------------
  !> @brief Retorna as coordenadas cartesianas (x, y) correspondentes à
  !>        coordenada intrínseca eta em [-1, 1] a partir de pts(1:3).
  !---------------------------------------------------------------------
  pure function cartesian_point_from_intrp_func_1d(pts, eta) result(pt)
    type(point_2d), intent(in) :: pts(3)
    real(dp), intent(in)       :: eta
    type(point_2d) :: pt
    real(dp) :: ni(3)

    call intrp_func_1d(ni, eta, 1.0_dp, 1.0_dp)
    pt%x = ni(1)*pts(1)%x + ni(2)*pts(2)%x + ni(3)*pts(3)%x
    pt%y = ni(1)*pts(1)%y + ni(2)*pts(2)%y + ni(3)*pts(3)%y
  end function cartesian_point_from_intrp_func_1d

  !---------------------------------------------------------------------
  !> @brief Pontos e pesos da quadratura de Gauss-Legendre 1D em [-1, 1].
  !---------------------------------------------------------------------
  subroutine gauss_legendre_points_1d(io, coord, weight)
    integer, intent(in)   :: io
    real(dp), intent(out) :: coord(:), weight(:)

    select case (io)
    case (1)
      coord(1) = 0.0_dp
      weight(1) = 2.0_dp
    case (2)
      coord(1) = 0.577350269189626_dp; coord(2) = -coord(1)
      weight(1) = 1.0_dp; weight(2) = weight(1)
    case (3)
      coord(1) = 0.0_dp
      coord(2) = 0.774596669241483_dp; coord(3) = -coord(2)
      weight(1) = 0.888888888888889_dp
      weight(2) = 0.555555555555556_dp; weight(3) = weight(2)
    case (4)
      coord(1) = 0.339981043584856_dp; coord(2) = -coord(1)
      coord(3) = 0.861136311594053_dp; coord(4) = -coord(3)
      weight(1) = 0.652145154862546_dp; weight(2) = weight(1)
      weight(3) = 0.347854845137454_dp; weight(4) = weight(3)
    case (5)
      coord(1) = 0.0_dp
      coord(2) = 0.538469310105683_dp; coord(3) = -coord(2)
      coord(4) = 0.906179845938664_dp; coord(5) = -coord(4)
      weight(1) = 0.568888888888889_dp
      weight(2) = 0.478628670499366_dp; weight(3) = weight(2)
      weight(4) = 0.236926885056189_dp; weight(5) = weight(4)
    case (6)
      coord(1) = 0.238619186083197_dp; coord(2) = -coord(1)
      coord(3) = 0.661209386466265_dp; coord(4) = -coord(3)
      coord(5) = 0.932469514203152_dp; coord(6) = -coord(5)
      weight(1) = 0.467913934572691_dp; weight(2) = weight(1)
      weight(3) = 0.360761573048139_dp; weight(4) = weight(3)
      weight(5) = 0.171324492379170_dp; weight(6) = weight(5)
    case (7)
      coord(1) = 0.0_dp
      coord(2) = 0.405845151377397_dp; coord(3) = -coord(2)
      coord(4) = 0.741531185599394_dp; coord(5) = -coord(4)
      coord(6) = 0.949107912342759_dp; coord(7) = -coord(6)
      weight(1) = 0.417959183673469_dp
      weight(2) = 0.381830050505119_dp; weight(3) = weight(2)
      weight(4) = 0.279705391489277_dp; weight(5) = weight(4)
      weight(6) = 0.129484966168870_dp; weight(7) = weight(6)
    case (8)
      coord(1) = 0.183434642495650_dp; coord(2) = -coord(1)
      coord(3) = 0.525532409916329_dp; coord(4) = -coord(3)
      coord(5) = 0.796666477413627_dp; coord(6) = -coord(5)
      coord(7) = 0.960289856497536_dp; coord(8) = -coord(7)
      weight(1) = 0.362683783378362_dp; weight(2) = weight(1)
      weight(3) = 0.313706645877887_dp; weight(4) = weight(3)
      weight(5) = 0.222381034453374_dp; weight(6) = weight(5)
      weight(7) = 0.101228536290376_dp; weight(8) = weight(7)
    case (9)
      coord(1) = 0.0_dp
      coord(2) = 0.324253423403809_dp; coord(3) = -coord(2)
      coord(4) = 0.613371432700590_dp; coord(5) = -coord(4)
      coord(6) = 0.836031107326636_dp; coord(7) = -coord(6)
      coord(8) = 0.968160239507626_dp; coord(9) = -coord(8)
      weight(1) = 0.330239355001260_dp
      weight(2) = 0.312347077040003_dp; weight(3) = weight(2)
      weight(4) = 0.260610696402935_dp; weight(5) = weight(4)
      weight(6) = 0.180648160694857_dp; weight(7) = weight(6)
      weight(8) = 0.081274388361574_dp; weight(9) = weight(8)
    case (10)
      coord(1) = 0.148874338981631_dp; coord(2) = -coord(1)
      coord(3) = 0.433395394129247_dp; coord(4) = -coord(3)
      coord(5) = 0.679409568299024_dp; coord(6) = -coord(5)
      coord(7) = 0.865063366688985_dp; coord(8) = -coord(7)
      coord(9) = 0.973906528517172_dp; coord(10) = -coord(9)
      weight(1) = 0.295524224714753_dp; weight(2) = weight(1)
      weight(3) = 0.269266719309996_dp; weight(4) = weight(3)
      weight(5) = 0.219086362515982_dp; weight(6) = weight(5)
      weight(7) = 0.149451349150581_dp; weight(8) = weight(7)
      weight(9) = 0.066671344308688_dp; weight(10) = weight(9)
    case (12)
      coord(1) = 0.125233408511469_dp; coord(2) = -coord(1)
      coord(3) = 0.367831498998180_dp; coord(4) = -coord(3)
      coord(5) = 0.587317954286617_dp; coord(6) = -coord(5)
      coord(7) = 0.769902674194305_dp; coord(8) = -coord(7)
      coord(9) = 0.904117256370475_dp; coord(10) = -coord(9)
      coord(11) = 0.981560634246719_dp; coord(12) = -coord(11)
      weight(1) = 0.249147045813403_dp; weight(2) = weight(1)
      weight(3) = 0.233492536538355_dp; weight(4) = weight(3)
      weight(5) = 0.203167426723066_dp; weight(6) = weight(5)
      weight(7) = 0.160078328543346_dp; weight(8) = weight(7)
      weight(9) = 0.106939325995318_dp; weight(10) = weight(9)
      weight(11) = 0.047175336386512_dp; weight(12) = weight(11)
    case (20)
      coord(1) = 0.076526521133497_dp; coord(2) = -coord(1)
      coord(3) = 0.227785851141645_dp; coord(4) = -coord(3)
      coord(5) = 0.373706088715420_dp; coord(6) = -coord(5)
      coord(7) = 0.510867001950827_dp; coord(8) = -coord(7)
      coord(9) = 0.636053680726515_dp; coord(10) = -coord(9)
      coord(11) = 0.746331906460150_dp; coord(12) = -coord(11)
      coord(13) = 0.839116971822218_dp; coord(14) = -coord(13)
      coord(15) = 0.912234428251325_dp; coord(16) = -coord(15)
      coord(17) = 0.963971927277913_dp; coord(18) = -coord(17)
      coord(19) = 0.993128599185094_dp; coord(20) = -coord(19)
      weight(1) = 0.152753387130726_dp; weight(2) = weight(1)
      weight(3) = 0.149172986472603_dp; weight(4) = weight(3)
      weight(5) = 0.142096109318382_dp; weight(6) = weight(5)
      weight(7) = 0.131688638449176_dp; weight(8) = weight(7)
      weight(9) = 0.118194531961518_dp; weight(10) = weight(9)
      weight(11) = 0.101930119817240_dp; weight(12) = weight(11)
      weight(13) = 0.083276741576704_dp; weight(14) = weight(13)
      weight(15) = 0.062672048334109_dp; weight(16) = weight(15)
      weight(17) = 0.040601429800386_dp; weight(18) = weight(17)
      weight(19) = 0.017614007139152_dp; weight(20) = weight(19)
    case default
      call error_message_int_stop("Ordem nao definida para Gauss-Legendre: ", io)
    end select
  end subroutine gauss_legendre_points_1d

  !---------------------------------------------------------------------
  !> @brief Pontos e pesos da quadratura logarítmica de Gauss-Laguerre
  !>        em [0, 1] para integração de núcleos fracamente singulares ln(1/r).
  !---------------------------------------------------------------------
  subroutine gauss_laguerre_points_1d(io, coord, weight)
    integer, intent(in)   :: io
    real(dp), intent(out) :: coord(:), weight(:)

    select case (io)
    case (8)
      coord(1) = 0.013320244160892_dp; coord(2) = 0.079750429013895_dp
      coord(3) = 0.197871029326188_dp; coord(4) = 0.354153994351909_dp
      coord(5) = 0.529458575234917_dp; coord(6) = 0.701814529939099_dp
      coord(7) = 0.849379320441106_dp; coord(8) = 0.953326450056359_dp
      weight(1) = 0.164416604728002_dp; weight(2) = 0.237525610023306_dp
      weight(3) = 0.226841984431919_dp; weight(4) = 0.175754079006070_dp
      weight(5) = 0.112924030246759_dp; weight(6) = 0.057872210717782_dp
      weight(7) = 0.020979073742133_dp; weight(8) = 0.003686407104028_dp
    case (12)
      coord(1) = 0.006548722279080_dp; coord(2) = 0.038946809560450_dp
      coord(3) = 0.098150263106005_dp; coord(4) = 0.181138581590633_dp
      coord(5) = 0.283220067667315_dp; coord(6) = 0.398434435164983_dp
      coord(7) = 0.519952626779130_dp; coord(8) = 0.640510916775482_dp
      coord(9) = 0.752865011892611_dp; coord(10) = 0.850240024421055_dp
      coord(11) = 0.926749682988251_dp; coord(12) = 0.977756129778486_dp
      weight(1) = 0.093192691443930_dp; weight(2) = 0.149751827576329_dp
      weight(3) = 0.166557454364573_dp; weight(4) = 0.159633559436994_dp
      weight(5) = 0.138424831864748_dp; weight(6) = 0.110016570636057_dp
      weight(7) = 0.079961821776732_dp; weight(8) = 0.052406954780971_dp
      weight(9) = 0.030071089000749_dp; weight(10) = 0.014249245402529_dp
      weight(11) = 0.004899924710876_dp; weight(12) = 0.000834029009810_dp
    case default
      call error_message_int_stop("Ordem nao definida para Gauss-Laguerre: ", io)
    end select
  end subroutine gauss_laguerre_points_1d

end module quadrature_mod
