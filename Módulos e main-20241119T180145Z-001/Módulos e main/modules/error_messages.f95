module error_messages
!----------------------------------------------------
! A module for for manipulating errors.
!
! Date: May 2018
! Last update: November 2021
!----------------------------------------------------
implicit none
!
contains
  subroutine error_message_stop(e_msg)
  !----------------------------------------------------
  ! Print error message and stops the program.
  !
  ! Date: 29 May 2018
  ! Last update: ---
  !----------------------------------------------------
  character(len = *),intent(in) :: e_msg
  print*
  print*,"*** Error: ",e_msg
  print*
  stop 1
  end subroutine error_message_stop
  !****************************************************
  !****************************************************
  subroutine error_2_messages_stop(e_msg1,e_msg2)
  !----------------------------------------------------
  ! Print error message and stops the program.
  !
  ! Date: 27 July 2018
  ! Last update: ---
  !----------------------------------------------------
  character(len = *),intent(in) :: e_msg1,e_msg2
  print*
  print*,"*** Error: ",e_msg1
  print*
  print*,"           ",e_msg2
  print*
  stop 1
  end subroutine error_2_messages_stop
  !****************************************************
  !****************************************************
  subroutine error_msg_2real_stop(e_msg,r1,r2)
  !----------------------------------------------------
  ! Print error message and stops the program.
  !
  ! Date: 12 November 2021
  ! Last update: ---
  !----------------------------------------------------
  character(len = *),intent(in) :: e_msg
  real(kind=8),intent(in) :: r1,r2
  print*
  print*,"*** Error: ",e_msg
  print*
  print*,"           (",r1,", ",r2,")"
  print*
  stop 1
  end subroutine error_msg_2real_stop
  !****************************************************
  !****************************************************
  subroutine error_message_int_stop(e_msg,n)
  !----------------------------------------------------
  ! Print error message and stops the program.
  !
  ! Date: 25 July 2018
  ! Last update: ---
  !----------------------------------------------------
  character(len = *),intent(in) :: e_msg
  integer,intent(in) :: n
  print*
  print*,"*** Error: ",e_msg,n
  print*
  stop 1
  end subroutine error_message_int_stop
  !****************************************************
  !****************************************************
  subroutine error_message_real_stop(e_msg,r)
  !----------------------------------------------------
  ! Print error message and stops the program.
  !
  ! Date: 26 July 2018
  ! Last update: ---
  !----------------------------------------------------
  character(len = *),intent(in) :: e_msg
  real,intent(in) :: r
  print*
  print*,"*** Error: ",e_msg,r
  print*
  stop 1
  end subroutine error_message_real_stop
  !****************************************************
  !****************************************************
  subroutine error_message_real8_stop(e_msg,r)
  !----------------------------------------------------
  ! Print error message and stops the program.
  !
  ! Date: 26 July 2018
  ! Last update: ---
  !----------------------------------------------------
  character(len = *),intent(in) :: e_msg
  real(kind=8),intent(in) :: r
  print*
  print*,"*** Error: ",e_msg,r
  print*
  stop 1
  end subroutine error_message_real8_stop
  !****************************************************
  !****************************************************
  subroutine error_message_2int_stop(e_msg1,n1,e_msg2,n2)
  !----------------------------------------------------
  ! Print error message and stops the program.
  !
  ! Date: 17 June 2020
  ! Last update: ---
  !----------------------------------------------------
  character(len = *),intent(in) :: e_msg1,e_msg2
  integer,intent(in) :: n1,n2
  print*
  print*,"*** Error: ",e_msg1,n1
  print*
  print*,"           ",e_msg2,n2
  print*
  stop 1
  end subroutine error_message_2int_stop
  !****************************************************
  !****************************************************
  subroutine error_message_real8_int_stop(e_msg1,n,e_msg2,r)
  !----------------------------------------------------
  ! Print error message and stops the program.
  !
  ! Date: 25 January 2019
  ! Last update: ---
  !----------------------------------------------------
  character(len = *),intent(in) :: e_msg1,e_msg2
  integer,intent(in) :: n
  real(kind=8),intent(in) :: r
  print*
  print*,"*** Error: ",e_msg1,n
  print*
  print*,"           ",e_msg2,r
  print*
  stop 1
  end subroutine error_message_real8_int_stop
  !****************************************************
  !****************************************************
  subroutine error_message_real8_2int_stop(e_msg1,n1,n2,e_msg2,r)
  !----------------------------------------------------
  ! Print error message and stops the program.
  !
  ! Date: 02 April 2020
  ! Last update: ---
  !----------------------------------------------------
  character(len = *),intent(in) :: e_msg1,e_msg2
  integer,intent(in) :: n1,n2
  real(kind=8),intent(in) :: r
  print*
  print*,"*** Error: ",e_msg1," ",n1," ",n2
  print*
  print*,"           ",e_msg2,r
  print*
  stop 1
  end subroutine error_message_real8_2int_stop
  !****************************************************
  !****************************************************
  subroutine warning_message(w_msg)
  !----------------------------------------------------
  ! Print warning message and pauses the program.
  !
  ! Date: 30 July 2018
  ! Last update: ---
  !----------------------------------------------------
  character(len = *),intent(in) :: w_msg
  print*
  print*,"*** Warning: ",w_msg
  print*
  end subroutine warning_message
  !****************************************************
  !****************************************************
  subroutine warning_2_messages(w_msg1,w_msg2)
  !----------------------------------------------------
  ! Print warning message and pauses the program.
  !
  ! Date: 11 December 2019
  ! Last update: ---
  !----------------------------------------------------
  character(len = *),intent(in) :: w_msg1,w_msg2
  print*
  print*,"*** Warning: ",w_msg1
  print*,"             ",w_msg2
  print*
  end subroutine warning_2_messages
  !****************************************************
  !****************************************************
end module error_messages
