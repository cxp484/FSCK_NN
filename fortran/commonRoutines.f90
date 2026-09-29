!!$**********************************************************************************
!!$                   SPECTRAL RADIATION CALCULATION SOFTWARE ~ SRCS ~
!!$**********************************************************************************
!!$
!!$                     COPYRIGHT (C) 2013 MICHAEL F. MODEST
!!$
!!$                           SCHOOL OF ENGINEERING                      
!!$                      UNIVERSITY OF CALIFORNIA AT MERCED
!!$
!!$	License: 
!!$		This file is a part of SRCS.
!!$		Used for non-commercial purposes only.
!!$		Not to be distributed without the permission of the copyright holder
!!$
!!$**********************************************************************************
!!$ common routines needed by all spectral modules

! Standalone subset: original strpad routine, unchanged.
module commonRoutines
implicit none
contains
subroutine strpad(st)
   integer i
   character,intent(inout)::st*(*)
   i = len_trim(st)
   do while (i .gt. 0)
    if (st(i:i) .eq. ' ') st(i:i)='0'
    i = i - 1
   enddo
   return
end subroutine strpad
end module commonRoutines
