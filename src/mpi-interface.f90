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

module mpi_interface

  !=======================================================================

  ! Define associated modules and subroutines

  !-----------------------------------------------------------------------

  use kinds

  !-------------------------------------------------------------------

  implicit none

  !-------------------------------------------------------------------

  ! Define necessary include files

  include "mpif.h"

  !-------------------------------------------------------------------

  ! Define global variables

  character                                              :: mpi_nodename(mpi_max_processor_name)
  logical                                                :: mpi_abort
  integer(kind=4),           dimension(:),   allocatable :: mpi_ranks
  integer(kind=4)                                        :: mpi_errorstatus(mpi_status_size)
  integer(kind=4)                                        :: mpi_masternode
  integer(kind=4)                                        :: mpi_ierror
  integer(kind=4)                                        :: mpi_procid
  integer(kind=4)                                        :: mpi_nprocs

  !-------------------------------------------------------------------

contains

  !===================================================================
  
  ! mpi_interface_initialize.f90:

  !-------------------------------------------------------------------

  subroutine mpi_interface_initialize()

    ! Initialize MPI session

    call mpi_init(mpi_ierror)

    ! Define rank for all nodes requested by user

    call mpi_comm_rank(mpi_comm_world,mpi_procid,mpi_ierror)

    ! Define the total number of nodes requested by user

    call mpi_comm_size(mpi_comm_world,mpi_nprocs,mpi_ierror)

    ! Define global variables

    mpi_masternode = 0

    ! Initialize global variables

    mpi_abort = .false.

  end subroutine mpi_interface_initialize

  !===================================================================

  ! mpi_interface_terminate.f90:

  !-------------------------------------------------------------------

  subroutine mpi_interface_terminate()

    ! Terminate MPI session

    call mpi_finalize(mpi_ierror)
    
  end subroutine mpi_interface_terminate

  !===================================================================

end module mpi_interface
