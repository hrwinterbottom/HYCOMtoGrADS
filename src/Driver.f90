!    A suite of routines to interpolate regional HYbrid Coordinate
!    Ocean Model (HYCOM) history (e.g., archv) direct-access binary
!    files to a user specified rectilinear grid covering the same
!    geographical region as the regional HYCOM experiment.

!    Copyright (C) 2013 Henry R. Winterbottom

!    Email: Henry.Winterbottom@noaa.gov

!    Snail-mail:

!    Henry R. Winterbottom
!    NOAA/OAR/PSD R/PSD1
!    325 Broadway
!    Boulder, CO 80303-3328

!    This program is free software; you can redistribute it and/or
!    modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation; either version 2 of
!    the License, or (at your option) any later version.

!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
!    General Public License for more details.

!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!    02110-1301 USA.

program hycomtograds_main

  !=====================================================================

  !$$$ PROGRAM DOCUMENTATION BLOCK
  !
  ! PRGMMR: Winterbottom        ORG: ESRL/PSD1       DATE: 2013-10-10
  !
  ! PROGRAM HISTORY LOG:
  !
  !  2013-10-10 Initial version. Henry R. Winterbottom
  !
  ! ATTRIBUTES:
  !
  !  language: f95
  !
  ! EXTERNAL I/O ROUTINES:
  !
  ! EXTERNAL MODULES:
  !$$$

  !=====================================================================

  ! Define associated modules and subroutines

  !---------------------------------------------------------------------

  use kinds

  !---------------------------------------------------------------------

  use mpi_interface
  use namelist
  use hycomtograds_interface

  !---------------------------------------------------------------------

  implicit none

  !=====================================================================

  ! Define variables computed within routine

  real(r_kind)                                             :: exectime_start
  real(r_kind)                                             :: exectime_finish

  !=====================================================================

  ! Initialize MPI session

  call mpi_interface_initialize()

  !---------------------------------------------------------------------

  ! If on root (master) task, perform all necessary tasks

  if(mpi_procid .eq. mpi_masternode) then

     ! Define local variable

     call cpu_time(exectime_start)

  end if ! if(mpi_procid .eq. mpi_masternode)

  ! Enable the root task to catch up from I/O and calculations

  call mpi_barrier(mpi_comm_world,mpi_ierror)

  !---------------------------------------------------------------------

  ! Compute local variables

  call hycomtograds()

  !---------------------------------------------------------------------

  ! If on root (master) task, perform all necessary tasks

  if(mpi_procid .eq. 0) then

     ! Define local variable

     call cpu_time(exectime_finish)

     ! Print message to user

     write(6,500) exectime_finish - exectime_start
  
  end if ! if(mpi_procid .eq. mpi_masternode)

  call mpi_finalize(mpi_ierror)

  !=====================================================================

  ! Define format statements

500 format('MAIN: Execution time: ', f13.5, ' seconds.')

  !=====================================================================

end program hycomtograds_main
