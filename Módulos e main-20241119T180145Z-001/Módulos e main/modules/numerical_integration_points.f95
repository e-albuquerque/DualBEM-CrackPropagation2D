module numerical_integration_points
!----------------------------------------------------
! A module which gives points and weights for
! numerical integrations.
!
! Date: July 2018
! Last update: September 2018
!----------------------------------------------------
use error_messages
implicit none
!
contains
  subroutine gauss_legendre_points_1d(io,coord,weight)
    !----------------------------------------------------
    ! Returns the coordinates and weights for numerical
    ! integration by standard Gaussian quadrature, given
    ! the integration order.
    !
    ! io = integration order.
    ! coord = an array with the coordinates.
    ! weight = an array with the corresponding weights.
    !
    ! Date: 25 July 2018
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(out) :: coord(:),weight(:)
    integer,intent(in) :: io
    !
    select case (io)
    case(1)
      coord(1) = 0.
      weight(1) = 2.
    case(2)
      coord(1) = 0.577350269189626_8; coord(2) = -coord(1)
      weight(1) = 1.; weight(2) = weight(1)
    case(3)
      coord(1) = 0.
      coord(2) = 0.774596669241483_8; coord(3) = -coord(2)
      weight(1) = 0.888888888888889_8
      weight(2) = 0.555555555555556_8; weight(3) = weight(2)
    case(4)
      coord(1) = 0.339981043584856_8; coord(2) = -coord(1)
      coord(3) = 0.861136311594053_8; coord(4) = -coord(3)
      weight(1) = 0.652145154862546_8; weight(2) = weight(1)
      weight(3) = 0.347854845137454_8; weight(4) = weight(3)
    case(5)
      coord(1) = 0.
      coord(2) = 0.538469310105683_8; coord(3) = -coord(2)
      coord(4) = 0.906179845938664_8; coord(5) = -coord(4)
      weight(1) = 0.568888888888889_8
      weight(2) = 0.478628670499366_8; weight(3) = weight(2)
      weight(4) = 0.236926885056189_8; weight(5) = weight(4)
    case(6)
      coord(1) = 0.238619186083197_8; coord(2) = -coord(1)
      coord(3) = 0.661209386466265_8; coord(4) = -coord(3)
      coord(5) = 0.932469514203152_8; coord(6) = -coord(5)
      weight(1) = 0.467913934572691_8; weight(2) = weight(1)
      weight(3) = 0.360761573048139_8; weight(4) = weight(3)
      weight(5) = 0.171324492379170_8; weight(6) = weight(5)
    case(7)
      coord(1) = 0.
      coord(2) = 0.405845151377397_8; coord(3) = -coord(2)
      coord(4) = 0.741531185599394_8; coord(5) = -coord(4)
      coord(6) = 0.949107912342759_8; coord(7) = -coord(6)
      weight(1) = 0.417959183673469_8
      weight(2) = 0.381830050505119_8; weight(3) = weight(2)
      weight(4) = 0.279705391489277_8; weight(5) = weight(4)
      weight(6) = 0.129484966168870_8; weight(7) = weight(6)
    case(8)
      coord(1) = 0.183434642495650_8; coord(2) = -coord(1)
      coord(3) = 0.525532409916329_8; coord(4) = -coord(3)
      coord(5) = 0.796666477413627_8; coord(6) = -coord(5)
      coord(7) = 0.960289856497536_8; coord(8) = -coord(7)
      weight(1) = 0.362683783378362_8; weight(2) = weight(1)
      weight(3) = 0.313706645877887_8; weight(4) = weight(3)
      weight(5) = 0.222381034453374_8; weight(6) = weight(5)
      weight(7) = 0.101228536290376_8; weight(8) = weight(7)
    case(9)
      coord(1) = 0.
      coord(2) = 0.324253423403809_8; coord(3) = -coord(2)
      coord(4) = 0.613371432700590_8; coord(5) = -coord(4)
      coord(6) = 0.836031107326636_8; coord(7) = -coord(6)
      coord(8) = 0.968160239507626_8; coord(9) = -coord(8)
      weight(1) = 0.330239355001260_8
      weight(2) = 0.312347077040003_8; weight(3) = weight(2)
      weight(4) = 0.260610696402935_8; weight(5) = weight(4)
      weight(6) = 0.180648160694857_8; weight(7) = weight(6)
      weight(8) = 0.081274388361574_8; weight(9) = weight(8)
    case(10)
      coord(1) = 0.148874338981631_8; coord(2) = -coord(1)
      coord(3) = 0.433395394129247_8; coord(4) = -coord(3)
      coord(5) = 0.679409568299024_8; coord(6) = -coord(5)
      coord(7) = 0.865063366688985_8; coord(8) = -coord(7)
      coord(9) = 0.973906528517172_8; coord(10) = -coord(9)
      weight(1) = 0.295524224714753_8; weight(2) = weight(1)
      weight(3) = 0.269266719309996_8; weight(4) = weight(3)
      weight(5) = 0.219086362515982_8; weight(6) = weight(5)
      weight(7) = 0.149451349150581_8; weight(8) = weight(7)
      weight(9) = 0.066671344308688_8; weight(10) = weight(9)
    case(12)
      coord(1) = 0.125233408511469_8; coord(2) = -coord(1)
      coord(3) = 0.367831498998180_8; coord(4) = -coord(3)
      coord(5) = 0.587317954286617_8; coord(6) = -coord(5)
      coord(7) = 0.769902674194305_8; coord(8) = -coord(7)
      coord(9) = 0.904117256370475_8; coord(10) = -coord(9)
      coord(11) = 0.981560634246719_8; coord(12) = -coord(11)
      weight(1) = 0.249147045813403_8; weight(2) = weight(1)
      weight(3) = 0.233492536538355_8; weight(4) = weight(3)
      weight(5) = 0.203167426723066_8; weight(6) = weight(5)
      weight(7) = 0.160078328543346_8; weight(8) = weight(7)
      weight(9) = 0.106939325995318_8; weight(10) = weight(9)
      weight(11) = 0.047175336386512_8; weight(12) = weight(11)
    case(20)
      coord(1) = 0.076526521133497_8; coord(2) = -coord(1)
      coord(3) = 0.227785851141645_8; coord(4) = -coord(3)
      coord(5) = 0.373706088715420_8; coord(6) = -coord(5)
      coord(7) = 0.510867001950827_8; coord(8) = -coord(7)
      coord(9) = 0.636053680726515_8; coord(10) = -coord(9)
      coord(11) = 0.746331906460150_8; coord(12) = -coord(11)
      coord(13) = 0.839116971822218_8; coord(14) = -coord(13)
      coord(15) = 0.912234428251325_8; coord(16) = -coord(15)
      coord(17) = 0.963971927277913_8; coord(18) = -coord(17)
      coord(19) = 0.993128599185094_8; coord(20) = -coord(19)
      weight(1) = 0.152753387130726_8; weight(2) = weight(1)
      weight(3) = 0.149172986472603_8; weight(4) = weight(3)
      weight(5) = 0.142096109318382_8; weight(6) = weight(5)
      weight(7) = 0.131688638449176_8; weight(8) = weight(7)
      weight(9) = 0.118194531961518_8; weight(10) = weight(9)
      weight(11) = 0.101930119817240_8; weight(12) = weight(11)
      weight(13) = 0.083276741576704_8; weight(14) = weight(13)
      weight(15) = 0.062672048334109_8; weight(16) = weight(15)
      weight(17) = 0.040601429800386_8; weight(18) = weight(17)
      weight(19) = 0.017614007139152_8; weight(20) = weight(19)
    case default !error message
        call error_message_int_stop("Undefined integration order for Gauss-Legendre method: ",io)
    end select
    return
  end subroutine gauss_legendre_points_1d
  !****************************************************
  !****************************************************
  subroutine gauss_laguerre_points_1d(io,coord,weight)
    !----------------------------------------------------
    ! Returns the coordinates and weights for numerical
    ! integration by Gauss-Laguerre quadrature of
    ! logarithmic kernels, given the integration order.
    !
    ! io = integration order.
    ! coord = an array with the coordinates.
    ! weight = an array with the corresponding weights.
    !
    ! Date: 21 September 2018
    ! Last update: ---
    !----------------------------------------------------
    real(kind=8),intent(out) :: coord(:),weight(:)
    integer,intent(in) :: io
    !
    select case (io)
    case(8)
      coord(1) = 0.013320244160892_8; coord(2) = 0.079750429013895_8
      coord(3) = 0.197871029326188_8; coord(4) = 0.354153994351909_8
      coord(5) = 0.529458575234917_8; coord(6) = 0.701814529939099_8
      coord(7) = 0.849379320441106_8; coord(8) = 0.953326450056359_8
      weight(1) = 0.164416604728002_8; weight(2) = 0.237525610023306_8
      weight(3) = 0.226841984431919_8; weight(4) = 0.175754079006070_8
      weight(5) = 0.112924030246759_8; weight(6) = 0.057872210717782_8
      weight(7) = 0.020979073742133_8; weight(8) = 0.003686407104028_8
    case(12)
      coord(1) = 0.006548722279080_8; coord(2) = 0.038946809560450_8
      coord(3) = 0.098150263106005_8; coord(4) = 0.181138581590633_8
      coord(5) = 0.283220067667315_8; coord(6) = 0.398434435164983_8
      coord(7) = 0.519952626779130_8; coord(8) = 0.640510916775482_8
      coord(9) = 0.752865011892611_8; coord(10) = 0.850240024421055_8
      coord(11) = 0.926749682988251_8; coord(12) = 0.977756129778486_8
      weight(1) = 0.093192691443930_8; weight(2) = 0.149751827576329_8
      weight(3) = 0.166557454364573_8; weight(4) = 0.159633559436994_8
      weight(5) = 0.138424831864748_8; weight(6) = 0.110016570636057_8
      weight(7) = 0.079961821776732_8; weight(8) = 0.052406954780971_8
      weight(9) = 0.030071089000749_8; weight(10) = 0.014249245402529_8
      weight(11) = 0.004899924710876_8; weight(12) = 0.000834029009810_8
    case default !error message
        call error_message_int_stop("Undefined integration order for Gauss-Laguerre method: ",io)
    end select
    return
  end subroutine gauss_laguerre_points_1d
  !****************************************************
  !****************************************************
end module numerical_integration_points
