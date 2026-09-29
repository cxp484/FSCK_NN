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
module FSKDBPath
character(512) :: basedbpath, nbpath, kppath, ktpath, oktpath, v4path

public:: set_fsk_db_paths

contains

subroutine set_fsk_db_paths(basepath)

implicit none
character(len=*), intent(in) :: basepath


! Remove any trailing newline from basepath
basedbpath = trim(basepath)

! Set paths using trimmed basepath
nbpath = trim(basedbpath) // '/nbkgDB/'
kppath = trim(basedbpath) // '/kpDB/'
ktpath = trim(basedbpath) // '/FSKTable/'
oktpath = trim(basedbpath) // '/FSKTableOpt/'
v4path = trim(basedbpath) // '/FSKTableV4/'

! Display other paths for debugging
!write(*,*) "nbpath=", nbpath
!write(*,*) "kppath=", kppath
!write(*,*) "ktpath=", ktpath
!write(*,*) "oktpath=", oktpath



end subroutine set_fsk_db_paths


end module 
