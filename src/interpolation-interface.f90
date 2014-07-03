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

module interpolation_interface

  !=======================================================================

  ! Define associated modules and subroutines

  !-----------------------------------------------------------------------

  use constants
  use kinds

  !-----------------------------------------------------------------------

  use kdtree2_module
  use mpi_interface
  use namelist

  !-----------------------------------------------------------------------

  implicit none

  !-----------------------------------------------------------------------

  ! Define all data and structure types for routine; these variables
  ! are variables required by the subroutines within this module

  type interpgrid
     type(kdtree2),        pointer                               :: kdtree_grid
     real(r_kind),                 dimension(:,:,:), allocatable :: weights
     real(r_kind),                 dimension(:,:),   allocatable :: anlysvar
     real(r_kind),                 dimension(:,:),   allocatable :: depth_profile
     real(r_kind),                 dimension(:),     allocatable :: xlong
     real(r_kind),                 dimension(:),     allocatable :: xlat
     real(r_kind),                 dimension(:),     allocatable :: rotang
     real(r_kind),                 dimension(:),     allocatable :: depth
     real(r_kind),                 dimension(:),     allocatable :: temperature
     real(r_kind),                 dimension(:),     allocatable :: scutoff
     real(r_kind)                                                :: distthresh
     real(r_kind)                                                :: dx
     real(r_kind)                                                :: dy
     real(r_kind),                 dimension(:,:),   allocatable :: grdloc
     integer,                      dimension(:,:),   allocatable :: grdnbors
     integer,                      dimension(:),     allocatable :: mpi_count_begin
     integer,                      dimension(:),     allocatable :: mpi_count_end
     integer                                                     :: mpi_maxprocid
     integer                                                     :: ncoords
     integer                                                     :: nvertlevs
     integer                                                     :: npasses
     integer                                                     :: neighbors
  end type interpgrid

  ! Define global variables

  real(r_kind),                   dimension(:),      allocatable :: srcgrid_var
  real(r_kind),                   dimension(:),      allocatable :: dstgrid_var
  real(r_kind)                                                   :: dx
  real(r_kind)                                                   :: dy
  real(r_kind)                                                   :: rlon_min
  real(r_kind)                                                   :: rlon_max
  real(r_kind)                                                   :: rlat_min
  real(r_kind)                                                   :: rlat_max
  integer                                                        :: nlon
  integer                                                        :: nlat

  !-----------------------------------------------------------------------

contains

  !=======================================================================

  ! interpolation_initialize_task_balance.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_initialize_task_balance(grid)

    ! Define variables passed to routine

    type(interpgrid)                                                     :: grid

    ! Define variables computed within routine

    integer                                                              :: mpi_count_interval

    ! Define counting variables

    integer                                                              :: i, j, k, l
    integer                                                              :: count

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(grid%mpi_count_begin))                              &
         & allocate(grid%mpi_count_begin((mpi_nprocs)))
    if(.not. allocated(grid%mpi_count_end))                                &
         & allocate(grid%mpi_count_end((mpi_nprocs)))

    !---------------------------------------------------------------------

    ! If on master (root) node (task), define problem and broadcast
    ! variables to each slave (compute) node (task)

    if(mpi_procid .eq. mpi_masternode) then

       ! Initialize counting variable

       count = 1

       ! Compute local variable

       mpi_count_interval = grid%ncoords/(mpi_nprocs - 1)

       ! Initialize local variables

       grid%mpi_count_begin    = 0
       grid%mpi_count_end      = 0
       grid%mpi_count_begin(1) = 1
       grid%mpi_count_end(1)   = grid%mpi_count_begin(1) +                 &
            & mpi_count_interval

       ! Loop through total number of processors

       do l = 2, mpi_nprocs - 1

          ! Define local variables

          grid%mpi_count_begin(l) = grid%mpi_count_end(l-1) + 1
          grid%mpi_count_end(l)   = grid%mpi_count_begin(l) +              &
               & mpi_count_interval

          ! Check local variables and proceed accordingly

          if(grid%mpi_count_begin(l) .gt. grid%ncoords) then

             ! Define local variables

             grid%mpi_count_begin(l) = grid%ncoords
             grid%mpi_count_end(l)   = grid%ncoords
             grid%mpi_maxprocid      = l

             ! Define exit from loop

             goto 1000

          end if ! if(grid%mpi_count_begin(l) .gt. grid%ncoords)

          ! Check local variables and proceed accordingly

          if(grid%mpi_count_end(l) .gt. grid%ncoords) then

             ! Define local variables

             grid%mpi_count_end(l) = grid%ncoords
             grid%mpi_maxprocid    = l

             ! Define exit from loop

             goto 1000

          end if ! if(grid%mpi_count_end(l) .gt. grid%ncoords)

       end do ! do l = 2, (mpi_nprocs - 1)

       ! Define exit from loop

1000   continue

       ! Loop through local variable and proceed accordingly

       do l = 1, grid%mpi_maxprocid

          ! Print message to user

          if(debug) write(6,500) grid%ncoords, l, grid%mpi_count_begin(l), &
               & grid%mpi_count_end(l)

       end do ! do l = 1, grid%mpi_maxprocid

    end if ! if(mpi_procid .eq. mpi_masternode)

    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(grid%mpi_count_begin,mpi_nprocs,mpi_integer,            &
         & mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_bcast(grid%mpi_count_end,mpi_nprocs,mpi_integer,              &
         & mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_bcast(grid%mpi_maxprocid,1,mpi_integer,mpi_masternode,        &
         & mpi_comm_world,mpi_ierror)

    !=====================================================================

    ! Define format statements

500 format('INTERPOLATION_INITIALIZE_TASK_BALANCE: (grid size/',           &
         & 'task ID/tile min/tile max) : ', i9, i6, i9, i9)

    !=====================================================================

  end subroutine interpolation_initialize_task_balance

  !=======================================================================

  ! interpolation_cleanup_task_balance.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_cleanup_task_balance(grid)

    ! Define variables passed to routine

    type(interpgrid)                                                     :: grid

    !=====================================================================

    ! Deallocate memory for local variable

    if(allocated(grid%mpi_count_begin))                                    &
         & deallocate(grid%mpi_count_begin)
    if(allocated(grid%mpi_count_end))                                      &
         & deallocate(grid%mpi_count_end)

    !=====================================================================

  end subroutine interpolation_cleanup_task_balance

  !=======================================================================

  ! interpolation_initialize_grid.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_initialize_grid(grid)

    ! Define variables passed to routine

    type(interpgrid)                                         :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    call interpolation_cleanup_grid(grid)

    !---------------------------------------------------------------------

    ! Allocate memory for local variables

    if(.not. allocated(grid%grdloc))                                       &
         & allocate(grid%grdloc(3,grid%ncoords))
    if(.not. allocated(grid%weights))                                      &
         & allocate(grid%weights(grid%ncoords,grid%neighbors,              &
         & grid%npasses))
    if(.not. allocated(grid%xlong))                                        &
         & allocate(grid%xlong(grid%ncoords))
    if(.not. allocated(grid%xlat))                                         &
         & allocate(grid%xlat(grid%ncoords))
    if(.not. allocated(grid%rotang))                                       &
         & allocate(grid%rotang(grid%ncoords))
    if(.not. allocated(grid%grdnbors))                                     &
         & allocate(grid%grdnbors(grid%ncoords,grid%neighbors))
    if(.not. allocated(grid%scutoff))                                      &
         & allocate(grid%scutoff(grid%npasses))
    if(.not. allocated(grid%depth))                                        &
         & allocate(grid%depth(grid%ncoords))
    if(.not. allocated(grid%temperature))                                  &
         & allocate(grid%temperature(grid%ncoords))
    if(.not. allocated(grid%depth_profile))                                &
         & allocate(grid%depth_profile(grid%nvertlevs,grid%ncoords))

    !=====================================================================

  end subroutine interpolation_initialize_grid

  !=======================================================================

  ! interpolation_cleanup_grid.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_cleanup_grid(grid)

    ! Define variables passed to routine

    type(interpgrid)                                         :: grid

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(grid%grdloc))                                             &
         & deallocate(grid%grdloc)
    if(allocated(grid%weights))                                            &
         & deallocate(grid%weights)
    if(allocated(grid%xlong))                                              &
         & deallocate(grid%xlong)
    if(allocated(grid%xlat))                                               &
         & deallocate(grid%xlat)
    if(allocated(grid%rotang))                                             &
         & deallocate(grid%rotang)
    if(allocated(grid%grdnbors))                                           &
         & deallocate(grid%grdnbors)
    if(allocated(grid%scutoff))                                            &
         & deallocate(grid%scutoff)
    if(allocated(grid%depth))                                              &
         & deallocate(grid%depth)
    if(allocated(grid%temperature))                                        &
         & deallocate(grid%temperature)
    if(allocated(grid%depth_profile))                                      &
         & deallocate(grid%depth_profile)
    if(mpi_procid .eq. mpi_masternode .and. associated(grid%kdtree_grid))  &
         & call kdtree2_destroy(grid%kdtree_grid)

    !=====================================================================

  end subroutine interpolation_cleanup_grid

  !=======================================================================

  ! interpolation_define_kdtree.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_define_kdtree(grid)

    ! Define variables passed to routine

    type(interpgrid)                                         :: grid

    ! Define counting variables

    integer                                                  :: i, j, k

    !=====================================================================

    ! Initialize local variables

    call init_constants_derived()

    !---------------------------------------------------------------------

    ! Loop through total number of grid coordinates

    do j = 1, grid%ncoords

       ! Compute local variables

       grid%grdloc(1,j) = rearth_equator*cos(grid%xlat(j))*                &
            & cos(grid%xlong(j))
       grid%grdloc(2,j) = rearth_equator*cos(grid%xlat(j))*                &
            & sin(grid%xlong(j))
       grid%grdloc(3,j) = rearth_equator*sin(grid%xlat(j))

    end do ! do j = 1, grid%ncoords

    ! Initialize local variable

    grid%kdtree_grid => kdtree2_create(grid%grdloc,sort=.true.,            &
         & rearrange=.true.)

    !=====================================================================

  end subroutine interpolation_define_kdtree

  !=======================================================================

  ! interpolation_define_weights_mpi.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_define_weights_mpi(srcgrid,dstgrid)

    ! Define variables passed to routine

    type(interpgrid)                                                     :: srcgrid
    type(interpgrid)                                                     :: dstgrid

    ! Define variables computed within routine

    type(kdtree2_result),       dimension(srcgrid%neighbors)             :: sresults
    real(r_kind),               dimension(:,:,:),            allocatable :: mpi_weights
    real(r_kind),               dimension(:,:),              allocatable :: mpi_sresults_dis
    real(r_kind),               dimension(:,:),              allocatable :: mpi_scutoff
    real(r_kind)                                                         :: mpi_distthresh

    ! Define counting variables

    integer                                                              :: i, j, k, l 

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(mpi_weights))                                       &     
         & allocate(mpi_weights(dstgrid%ncoords,dstgrid%neighbors,         &
         & dstgrid%npasses))
    if(.not. allocated(mpi_scutoff))                                       &
         & allocate(mpi_scutoff(dstgrid%ncoords,dstgrid%npasses))
    if(.not. allocated(mpi_sresults_dis))                                  &
         & allocate(mpi_sresults_dis(dstgrid%ncoords,dstgrid%neighbors))

    !---------------------------------------------------------------------

    ! Initialize local variables

    dstgrid%weights = 0.0
    mpi_weights     = 0.0

    !---------------------------------------------------------------------

    ! If on master (root) node (task), define problem and broadcast
    ! variables to each slave (compute) node (task)
       
    if(mpi_procid .eq. mpi_masternode) then
          
       ! Loop through local variable and proceed accordingly

       do i = 1, dstgrid%ncoords

          ! Define local variable
             
          call kdtree2_n_nearest(tp=srcgrid%kdtree_grid,                   &
               & qv=dstgrid%grdloc(:,i),nn=dstgrid%neighbors,              &
               & results=sresults)
             
          ! Define local variables
             
          mpi_distthresh                          =                        &
               & dstgrid%distthresh
          mpi_scutoff(i,1:dstgrid%npasses)        =                        &
               & dstgrid%scutoff(1:dstgrid%npasses)
          mpi_sresults_dis(i,1:dstgrid%neighbors) =                        &
               & sresults(1:dstgrid%neighbors)%dis
          dstgrid%grdnbors(i,1:dstgrid%neighbors) =                        &
               & sresults(1:dstgrid%neighbors)%idx

       end do ! do i = 1, dstgrid%ncoords
       
    end if ! if(mpi_procid .eq. mpi_masternode)
       
    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(mpi_distthresh,1,mpi_real,mpi_masternode,               &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(mpi_scutoff,(dstgrid%ncoords*dstgrid%npasses),          &
         & mpi_real,mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_bcast(mpi_sresults_dis,(dstgrid%ncoords*dstgrid%neighbors),   &
         & mpi_real,mpi_masternode,mpi_comm_world,mpi_ierror)
    call mpi_bcast(dstgrid%ncoords,1,mpi_integer,mpi_masternode,           &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(dstgrid%neighbors,1,mpi_integer,mpi_masternode,         &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(dstgrid%grdnbors,(dstgrid%ncoords*dstgrid%neighbors),   &
         & mpi_integer,mpi_masternode,mpi_comm_world,mpi_ierror)

    ! If on slave (compute) node (task), receive variables, compute
    ! variables, and send variables to master (root) node (task)
    
    if(mpi_procid .ne. mpi_masternode) then

       ! Check local variable and proceed accordingly

       if(dstgrid%mpi_count_begin(mpi_procid) .ne. 0 .and.                 &
            & dstgrid%mpi_count_end(mpi_procid) .ne. 0) then

          ! Loop through local variable and proceed accordingly
       
          do k = dstgrid%mpi_count_begin(mpi_procid),                      &
               & dstgrid%mpi_count_end(mpi_procid)
          
             ! Loop through total number of neighboring points on
             ! destination grid

             do i = 1, dstgrid%neighbors

                ! Loop through the total number of analysis passes to
                ! perform

                do j = 1, dstgrid%npasses

                   ! Compute local variable

                   mpi_weights(k,i,j) = exp((-1.0)*(                       &
                        & sqrt(mpi_sresults_dis(k,i))*                     &
                        & sqrt(mpi_sresults_dis(k,i)))/(4.0*               &
                        & mpi_scutoff(k,j)*mpi_distthresh*                 &
                        & mpi_distthresh))

                end do ! do j = 1, dstgrid%npasses

             end do ! do i = 1, dstgrid%neighbors

          end do ! do k = dstgrid%mpi_count_begin(mpi_procid),             &
                 ! dstgrid%mpi_count_end(mpi_procid)

       endif ! if(dstgrid%mpi_count_begin(mpi_procid) .ne. 0 .and.         &
             ! dstgrid%mpi_count_end(mpi_procid) .ne. 0)

    end if ! if(mpi_procid .ne. mpi_masternode)

    !---------------------------------------------------------------------

    ! Define local variable

    call mpi_reduce(mpi_weights(1:dstgrid%ncoords,1:dstgrid%neighbors,1:   &
         & dstgrid%npasses),dstgrid%weights(1:dstgrid%ncoords,1:           &
         & dstgrid%neighbors,1:dstgrid%npasses),(dstgrid%ncoords*          &
         & dstgrid%neighbors*dstgrid%npasses),mpi_real,mpi_sum,            &
         & mpi_masternode,mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! Deallocate memory for local variable on all compute tasks

    if(allocated(mpi_weights))      deallocate(mpi_weights)
    if(allocated(mpi_sresults_dis)) deallocate(mpi_sresults_dis)
    if(allocated(mpi_scutoff))      deallocate(mpi_scutoff)

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

  end subroutine interpolation_define_weights_mpi

  !=======================================================================

  ! interpolation_barnes_analysis_mpi.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

    ! Define variable passed to routine

    type(interpgrid)                                                     :: srcgrid

    ! Define variable returned by routine

    type(interpgrid)                                                     :: dstgrid

    ! Define variables computed within routine

    real(r_kind),               dimension(:),                allocatable :: mpi_dstgrid_var
    real(r_kind),               dimension(:),                allocatable :: workgrid
    real(r_kind)                                                         :: weights_sum

    ! Define counting variables

    integer                                                              :: i, j, k, l

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(mpi_dstgrid_var))                                   &
         & allocate(mpi_dstgrid_var(dstgrid%ncoords))
    if(.not. allocated(workgrid))                                          &
         & allocate(workgrid(dstgrid%ncoords))

    !---------------------------------------------------------------------

    ! Initialize local variables

    mpi_dstgrid_var = 0.0

    ! Check local variable and proceed accordingly

    if(mpi_procid .ne. mpi_masternode .and. mpi_procid .le.                &
         & dstgrid%mpi_maxprocid) then

       ! Loop through local variable and proceed accordingly
       
       do k = dstgrid%mpi_count_begin(mpi_procid),                         &
            & dstgrid%mpi_count_end(mpi_procid)

          ! Loop through the total number of analysis passes
          
          do j = 1, dstgrid%npasses
             
             ! Initialize local variable
             
             weights_sum = 0.0
             
             ! Loop through total number of neighboring points on
             ! source grid

             do i = 1, dstgrid%neighbors
             
                ! Define local variable
                
                weights_sum = weights_sum + dstgrid%weights(k,i,j)
             
                ! Define local variable
             
                mpi_dstgrid_var(k) = mpi_dstgrid_var(k) +                  &
                     & (srcgrid_var(dstgrid%grdnbors(k,i)) -               &
                     & workgrid(k))*dstgrid%weights(k,i,j)
                
             end do ! do i = 1, dstgrid%neighbors

             ! Define local variable accordingly

             if(weights_sum .gt. barnes_weights_threshold .and. j .eq.     &
                  & 1) then

                ! Define local variable

                mpi_dstgrid_var(k) = mpi_dstgrid_var(k)/weights_sum

             else  ! if(weights_sum .gt. barnes_weights_threshold          &
                   ! .and. k .eq. 1)
                
                ! Define local variable

                mpi_dstgrid_var(k) = srcgrid_var(dstgrid%grdnbors(k,1))

             end if ! if(weights_sum .gt. barnes_weights_threshold         &
                    ! .and. k .eq. 1)
             
             ! Define local variable
             
             workgrid(k) = mpi_dstgrid_var(k)
                
          end do ! do j = 1, dstgrid%npasses

       end do ! do k = dstgrid%mpi_count_begin(mpi_procid),                &
              ! dstgrid%mpi_count_end(mpi_procid)

    endif ! if(mpi_procid .ne. mpi_masternode .and. mpi_procid             &
          ! .le. dstgrid%mpi_maxprocid)

    !---------------------------------------------------------------------

    ! Define local variable

    call mpi_reduce(mpi_dstgrid_var(1:dstgrid%ncoords),                    &
         & dstgrid_var(1:dstgrid%ncoords),dstgrid%ncoords,mpi_real,        &
         & mpi_sum,mpi_masternode,mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! Deallocate memory for local variable

    if(allocated(mpi_dstgrid_var)) deallocate(mpi_dstgrid_var)
    if(allocated(workgrid))        deallocate(workgrid)

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

  end subroutine interpolation_barnes_analysis_mpi

  !=======================================================================

  ! interpolation_nearest_neighbor_mpi.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_nearest_neighbor_mpi(srcgrid,dstgrid)

    ! Define variable passed to routine

    type(interpgrid)                                                     :: srcgrid

    ! Define variable returned by routine

    type(interpgrid)                                                     :: dstgrid

    ! Define variables computed within routine

    type(kdtree2_result),       dimension(srcgrid%neighbors)             :: sresults
    real(r_kind),               dimension(:),                allocatable :: mpi_dstgrid_var
    integer,                    dimension(:,:),              allocatable :: mpi_srcgrid_idx

    ! Define counting variables

    integer                                                              :: i, j, k, l

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(mpi_srcgrid_idx))                                   &
         & allocate(mpi_srcgrid_idx(dstgrid%ncoords,dstgrid%neighbors))
    if(.not. allocated(mpi_dstgrid_var))                                   &
         & allocate(mpi_dstgrid_var(dstgrid%ncoords))

    !---------------------------------------------------------------------

    ! Initialize local variables

    mpi_dstgrid_var = 0.0

    !---------------------------------------------------------------------

    ! If on master (root) node (task), define problem and broadcast
    ! variables to each slave (compute) node (task)

    if(mpi_procid .eq. mpi_masternode) then

       ! Loop through local variable and proceed accordingly

       do i = 1, dstgrid%ncoords

          ! Define local variable
             
          call kdtree2_n_nearest(tp=srcgrid%kdtree_grid,                   &
               & qv=dstgrid%grdloc(:,i),nn=dstgrid%neighbors,              &
               & results=sresults)
             
             ! Define local variable

          mpi_srcgrid_idx(i,1:dstgrid%neighbors) =                         &
               & sresults(1:dstgrid%neighbors)%idx

       end do ! i = 1, dstgrid%ncoords

    end if ! if(mpi_procid .eq. mpi_masternode)

    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(mpi_srcgrid_idx,(dstgrid%ncoords*dstgrid%neighbors),    &
         & mpi_real,mpi_masternode,mpi_comm_world,mpi_ierror)

    ! Check local variable and proceed accordingly
       
    if(mpi_procid .ne. mpi_masternode .and. mpi_procid .le.                &
         & dstgrid%mpi_maxprocid) then

       ! Loop through local variable and proceed accordingly
       
       do k = dstgrid%mpi_count_begin(mpi_procid),                         &
            & dstgrid%mpi_count_end(mpi_procid)
          
          ! Define local variable
          
          mpi_dstgrid_var(k) = srcgrid_var(mpi_srcgrid_idx(k,1))

       end do ! do k = dstgrid%mpi_count_begin(mpi_procid),                &
              ! dstgrid%mpi_count_end(mpi_procid)

    endif ! if(mpi_procid .ne. mpi_masternode .and. mpi_procid .le.        &
          ! dstgrid%mpi_maxprocid)

    !---------------------------------------------------------------------

    ! Define local variable

    call mpi_reduce(mpi_dstgrid_var(1:dstgrid%ncoords),                    &
         & dstgrid_var(1:dstgrid%ncoords),dstgrid%ncoords,mpi_real,        &
         & mpi_sum,mpi_masternode,mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! Deallocate memory for local variable

    if(allocated(mpi_srcgrid_idx)) deallocate(mpi_srcgrid_idx)
    if(allocated(mpi_dstgrid_var)) deallocate(mpi_dstgrid_var)

    !=====================================================================

    ! Return calculated values

    return

    !=====================================================================

  end subroutine interpolation_nearest_neighbor_mpi

  !=======================================================================

  ! interpolation_define_kdtree_mpi.f90:

  !-----------------------------------------------------------------------

  subroutine interpolation_define_kdtree_mpi(grid)

    ! Define variables passed to routine

    type(interpgrid)                                         :: grid

    ! Define variables computed within routine

    real(r_kind),              dimension(:,:),   allocatable :: mpi_grdloc

    ! Define counting variables

    integer                                                  :: i, j, k

    !=====================================================================

    ! Initialize local variables

    call init_constants_derived()

    !---------------------------------------------------------------------

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(grid%xlong,grid%ncoords,mpi_real,mpi_masternode,        &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(grid%xlat,grid%ncoords,mpi_real,mpi_masternode,         &
         & mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! Allocate memory for local variable

    if(.not. allocated(mpi_grdloc)) allocate(mpi_grdloc(3,grid%ncoords))

    ! Initialize local variables

    mpi_grdloc  = 0.0
    grid%grdloc = 0.0

    ! If on slave (compute) node (task), receive variables, compute
    ! variables, and send variables to master (root) node (task)
    
    if(mpi_procid .ne. mpi_masternode) then

       ! Check local variable and proceed accordingly

       if(grid%mpi_count_begin(mpi_procid) .ne. 0 .and.                    &
            & grid%mpi_count_end(mpi_procid) .ne. 0) then

          ! Loop through local variable and proceed accordingly
       
          do k = grid%mpi_count_begin(mpi_procid),                         &
               & grid%mpi_count_end(mpi_procid)
          
             ! Compute local variables

             mpi_grdloc(1,k) = rearth_equator*cos(grid%xlat(k))*           &
                  & cos(grid%xlong(k))
             mpi_grdloc(2,k) = rearth_equator*cos(grid%xlat(k))*           &
                  & sin(grid%xlong(k))
             mpi_grdloc(3,k) = rearth_equator*sin(grid%xlat(k))

          end do ! do k = grid%mpi_count_begin(mpi_procid),                &
                 ! grid%mpi_count_end(mpi_procid)

       end if ! if(grid%mpi_count_begin(mpi_procid) .ne. 0 .and.           &
              ! grid%mpi_count_end(mpi_procid) .ne. 0)

    end if ! if(mpi_procid .ne. mpi_masternode)

    ! Define local variable

    call mpi_reduce(mpi_grdloc(1:3,1:grid%ncoords),grid%grdloc(1:3,1:      &
         & grid%ncoords),(3*grid%ncoords),mpi_real,mpi_sum,                &
         & mpi_masternode,mpi_comm_world,mpi_ierror)  

    ! Deallocate memory for local variable

    if(allocated(mpi_grdloc)) deallocate(mpi_grdloc)

    !---------------------------------------------------------------------

    ! If on master (root) node (task), define problem and broadcast
    ! variables to each slave (compute) node (task)

    if(mpi_procid .eq. mpi_masternode) then

       ! Initialize local variable

       grid%kdtree_grid => kdtree2_create(grid%grdloc,sort=.true.,         &
            & rearrange=.true.)

    end if ! if(mpi_procid .eq. mpi_masternode)

    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    !=====================================================================

  end subroutine interpolation_define_kdtree_mpi

  !=======================================================================

  ! define_scaling_coefficients.f90:

  !-----------------------------------------------------------------------

  subroutine define_scaling_coefficients(grid)

    ! Define variable passed to routine

    type(interpgrid)                                         :: grid

    ! Define counting variables

    integer                                                  :: i, j, k

    !=====================================================================

    ! Initialize local variable

    grid%scutoff(1) = 1.0

    ! Loop through total number of threshold cutoff values

    do k = 1, grid%npasses

       ! Compute local variable

       grid%scutoff(k) = grid%scutoff(1)/(10**(real(k-1)))

    end do ! do k = 1, grid%npasses

    !=====================================================================

  end subroutine define_scaling_coefficients

  !=======================================================================

end module interpolation_interface
