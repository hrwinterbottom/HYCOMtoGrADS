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

module hycomtograds_interface

  !=======================================================================

  ! Define associated modules and subroutines

  !-----------------------------------------------------------------------

  use constants
  use kinds

  !-----------------------------------------------------------------------

  use interpolation_interface
  use mpi_interface
  use namelist
  use variable_interface

  !-----------------------------------------------------------------------

  implicit none

  !-----------------------------------------------------------------------

contains

  !=======================================================================

  ! hycomtograds.f90:

  !-----------------------------------------------------------------------

  subroutine hycomtograds()

    ! Define variables computed within routine
  
    type(hycomvariables)                                     :: hycom

    !=====================================================================

    ! Ingest external hycomtograds.input

    call namelistparams()

    ! Check local variable and proceed accordingly

    if(.not. is_isointerp .and. .not. is_zinterp) nlev = ocean_hycom_zdim

    !---------------------------------------------------------------------

    ! Define local variables
    
    call hycomtograds_interface_initialize_variables(hycom)
    call hycomtograds_interface_initialize(hycom)

    ! Compute local variables

    call variable_interface_process(hycom)
    
    ! Deallocate memory for local variables

    call hycomtograds_interface_cleanup()
    call hycomtograds_interface_cleanup_variables(hycom)

    !=====================================================================

  end subroutine hycomtograds

  !=======================================================================

  ! hycomtograds_interface_initialize_variables.f90:

  !-----------------------------------------------------------------------

  subroutine hycomtograds_interface_initialize_variables(grid)

    ! Define variables passed to routine

    type(hycomvariables)                                                  :: grid

    ! Define variables computed within routine

    real(r_kind),                          dimension(:),      allocatable :: w
    integer                                                               :: n2drec
    integer                                                               :: nrecl
    integer                                                               :: kcoord
    integer                                                               :: numrec
    integer                                                               :: recnum

    !=====================================================================

    ! Deallocate memory for local variables

    call hycomtograds_interface_cleanup_variables(grid)

    !---------------------------------------------------------------------

    ! Allocate memory for local variables

    if(.not. allocated(grid%uvel))                                         &
         & allocate(grid%uvel((ocean_hycom_xdim*ocean_hycom_ydim),         &
         & ocean_hycom_zdim))
    if(.not. allocated(grid%vvel))                                         &
         & allocate(grid%vvel((ocean_hycom_xdim*ocean_hycom_ydim),         &
         & ocean_hycom_zdim))
    if(.not. allocated(grid%thknss))                                       &
         & allocate(grid%thknss((ocean_hycom_xdim*ocean_hycom_ydim),       &
         & ocean_hycom_zdim))
    if(.not. allocated(grid%temp))                                         &
         & allocate(grid%temp((ocean_hycom_xdim*ocean_hycom_ydim),         &
         & ocean_hycom_zdim))
    if(.not. allocated(grid%salin))                                        &
         & allocate(grid%salin((ocean_hycom_xdim*ocean_hycom_ydim),        &
         & ocean_hycom_zdim))
    if(.not. allocated(grid%depth))                                        &
         & allocate(grid%depth((ocean_hycom_xdim*ocean_hycom_ydim),        &
         & ocean_hycom_zdim))
    if(.not. allocated(grid%montg1))                                       &
         & allocate(grid%montg1((ocean_hycom_xdim*ocean_hycom_ydim)))
    if(.not. allocated(grid%srfhgt))                                       &
         & allocate(grid%srfhgt((ocean_hycom_xdim*ocean_hycom_ydim)))
    if(.not. allocated(grid%steric))                                       &
         & allocate(grid%steric((ocean_hycom_xdim*ocean_hycom_ydim)))
    if(.not. allocated(grid%surflx))                                       &
         & allocate(grid%surflx((ocean_hycom_xdim*ocean_hycom_ydim)))
    if(.not. allocated(grid%salflx))                                       &
         & allocate(grid%salflx((ocean_hycom_xdim*ocean_hycom_ydim)))
    if(.not. allocated(grid%bl_dpth))                                      &
         & allocate(grid%bl_dpth((ocean_hycom_xdim*ocean_hycom_ydim)))
    if(.not. allocated(grid%mix_dpth))                                     &
         & allocate(grid%mix_dpth((ocean_hycom_xdim*ocean_hycom_ydim))) 
    if(.not. allocated(grid%covice))                                       &
         & allocate(grid%covice((ocean_hycom_xdim*ocean_hycom_ydim)))   
    if(.not. allocated(grid%thkice))                                       &
         & allocate(grid%thkice((ocean_hycom_xdim*ocean_hycom_ydim)))
    if(.not. allocated(grid%temice))                                       &
         & allocate(grid%temice((ocean_hycom_xdim*ocean_hycom_ydim)))
    if(.not. allocated(grid%u_btrop))                                      &
         & allocate(grid%u_btrop((ocean_hycom_xdim*ocean_hycom_ydim)))
    if(.not. allocated(grid%v_btrop))                                      &
         & allocate(grid%v_btrop((ocean_hycom_xdim*ocean_hycom_ydim)))

    !---------------------------------------------------------------------

    ! Compute local variables
       
    n2drec = ((ocean_hycom_xdim*ocean_hycom_ydim+4095)/4096)*4096

    ! Allocate memory for local variables
 
    if(.not. allocated(w)) allocate(w(n2drec)) 

    ! Define local variables

    inquire(iolength=nrecl) w
    numrec = nint(11 + (ocean_hycom_zdim)*5.0)
    open(11,file=trim(hycom_filename),form='unformatted',                   &
         & access='direct',recl=nrecl,action='read',convert='big_endian')

    ! Check local variable and proceed accordingly

    if(is_steric) then

       ! Define local variables

       read(11,rec=1)  grid%montg1
       read(11,rec=2)  grid%srfhgt
       read(11,rec=3)  grid%steric
       read(11,rec=4)  grid%surflx
       read(11,rec=5)  grid%salflx
       read(11,rec=6)  grid%bl_dpth
       read(11,rec=7)  grid%mix_dpth
       read(11,rec=8)  grid%covice
       read(11,rec=9)  grid%thkice
       read(11,rec=10) grid%temice
       read(11,rec=11) grid%u_btrop
       read(11,rec=12) grid%v_btrop

       ! Initialize local variables

       kcoord = 1
       recnum = 13
       
       ! Loop through local variable
       
       do kcoord = 1, ocean_hycom_zdim
          
          ! Define local variables
          
          read(11,rec=recnum) grid%uvel(:,kcoord)
          recnum = recnum + 1
          read(11,rec=recnum) grid%vvel(:,kcoord)
          recnum = recnum + 1
          read(11,rec=recnum) grid%thknss(:,kcoord)
          recnum = recnum + 1
          read(11,rec=recnum) grid%temp(:,kcoord)
          recnum = recnum + 1
          read(11,rec=recnum) grid%salin(:,kcoord)
          recnum = recnum + 1
          
       end do ! do kcoord = 1, ocean_hycom_zdim

    end if ! if(is_steric)

    ! Check local variable and proceed accordingly

    if(.not. is_steric) then

       ! Define local variables

       read(11,rec=1)  grid%montg1
       read(11,rec=2)  grid%srfhgt
       read(11,rec=3)  grid%surflx
       read(11,rec=4)  grid%salflx
       read(11,rec=5)  grid%bl_dpth
       read(11,rec=6)  grid%mix_dpth
       read(11,rec=7)  grid%covice
       read(11,rec=8)  grid%thkice
       read(11,rec=9) grid%temice
       read(11,rec=10) grid%u_btrop
       read(11,rec=11) grid%v_btrop

       ! Initialize local variables

       kcoord = 1
       recnum = 12
       
       ! Loop through local variable
       
       do kcoord = 1, ocean_hycom_zdim
          
          ! Define local variables
          
          read(11,rec=recnum) grid%uvel(:,kcoord)
          recnum = recnum + 1
          read(11,rec=recnum) grid%vvel(:,kcoord)
          recnum = recnum + 1
          read(11,rec=recnum) grid%thknss(:,kcoord)
          recnum = recnum + 1
          read(11,rec=recnum) grid%temp(:,kcoord)
          recnum = recnum + 1
          read(11,rec=recnum) grid%salin(:,kcoord)
          recnum = recnum + 1
          
       end do ! do kcoord = 1, ocean_hycom_zdim

    end if ! if(.not. is_steric)

    ! Define local variables
       
    close(11)

    ! Deallocate memory for local variables

    if(allocated(w)) deallocate(w)

    !---------------------------------------------------------------------

    ! Initialize local variable

    grid%depth = 0.0
    
    ! Loop through local variable

    do kcoord = 1, ocean_hycom_zdim - 1

       ! Compute local variables

       where(grid%thknss(:,kcoord) .lt. 1.e20) grid%depth(:,kcoord+1) =    &
            & grid%depth(:,kcoord) + grid%thknss(:,kcoord)/9806.0

    end do ! do kcoord = 1, ocean_hycom_zdim - 1

    !=====================================================================

  end subroutine hycomtograds_interface_initialize_variables

  !=======================================================================

  ! hycomtograds_interface_cleanup_variables.f90:

  !-----------------------------------------------------------------------

  subroutine hycomtograds_interface_cleanup_variables(grid)

    ! Define variables passed to routine

    type(hycomvariables)                                                  :: grid

    !=====================================================================

    ! Allocate memory for local variables

    if(allocated(grid%uvel))     deallocate(grid%uvel)
    if(allocated(grid%vvel))     deallocate(grid%vvel)
    if(allocated(grid%thknss))   deallocate(grid%thknss)
    if(allocated(grid%temp))     deallocate(grid%temp)
    if(allocated(grid%salin))    deallocate(grid%salin)
    if(allocated(grid%depth))    deallocate(grid%depth)
    if(allocated(grid%montg1))   deallocate(grid%montg1)
    if(allocated(grid%srfhgt))   deallocate(grid%srfhgt)
    if(allocated(grid%steric))   deallocate(grid%steric)
    if(allocated(grid%surflx))   deallocate(grid%surflx)
    if(allocated(grid%salflx))   deallocate(grid%salflx)
    if(allocated(grid%bl_dpth))  deallocate(grid%bl_dpth)
    if(allocated(grid%mix_dpth)) deallocate(grid%mix_dpth) 
    if(allocated(grid%covice))   deallocate(grid%covice)   
    if(allocated(grid%thkice))   deallocate(grid%thkice)
    if(allocated(grid%temice))   deallocate(grid%temice)
    if(allocated(grid%u_btrop))  deallocate(grid%u_btrop)
    if(allocated(grid%v_btrop))  deallocate(grid%v_btrop)

    !=====================================================================

  end subroutine hycomtograds_interface_cleanup_variables

  !=======================================================================

  ! hycomtograds_interface_initialize.f90:

  !-----------------------------------------------------------------------

  subroutine hycomtograds_interface_initialize(hycom)

    ! Define variables passed to routine

    type(hycomvariables)                                                  :: hycom

    ! Define variables computed within routine

    real(r_kind),                          dimension(:,:),    allocatable :: xlong
    real(r_kind),                          dimension(:,:),    allocatable :: xlat   
    real(r_kind),                          dimension(:,:),    allocatable :: rotang 
    real(r_kind),                          dimension(:,:),    allocatable :: depth
    real(r_kind),                          dimension(:),      allocatable :: w
    integer                                                               :: n2drec
    integer                                                               :: nrecl

    ! Define counting variables

    integer                                                               :: i, j, k
    integer                                                               :: count

    !=====================================================================

    ! Initialize local variables

    call init_constants_derived()

    !---------------------------------------------------------------------

    ! If on master (root) node (task), define problem and broadcast
    ! variables to each slave (compute) node (task)
       
    if(mpi_procid .eq. mpi_masternode) then

       ! Compute local variable
       
       n2drec = ((ocean_hycom_xdim*ocean_hycom_ydim+4095)/4096)*4096

       ! Allocate memory for local variables

       if(.not. allocated(xlong))                                          &
            & allocate(xlong(ocean_hycom_xdim,ocean_hycom_ydim))
       if(.not. allocated(xlat))                                           &
            & allocate(xlat(ocean_hycom_xdim,ocean_hycom_ydim))
       if(.not. allocated(rotang))                                         &
            & allocate(rotang(ocean_hycom_xdim,ocean_hycom_ydim))
       if(.not. allocated(depth))                                          &
            & allocate(depth(ocean_hycom_xdim,ocean_hycom_ydim))
       if(.not. allocated(w))                                              &
            & allocate(w(n2drec)) 

       ! Define local variables

       inquire(iolength=nrecl) w
       open(11,file=trim(regional_grid_filename),form='unformatted',       &
            & access='direct',recl=nrecl,action='read',                    &
            & convert='big_endian')
       read(11,rec=1) xlong
       read(11,rec=2) xlat
       read(11,rec=9) rotang
       close(11)
       open(11,file=trim(regional_depth_filename),form='unformatted',      &
            & access='direct',recl=nrecl,action='read',                    &
            & convert='big_endian')
       read(11,rec=1) depth
       close(11)

    end if ! if(mpi_procid .eq. mpi_masternode)
       
    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! Define local variables

    srcgrid%ncoords    = ocean_hycom_xdim*ocean_hycom_ydim
    srcgrid%npasses    = barnes_npasses
    srcgrid%neighbors  = barnes_nneighbors
    srcgrid%distthresh = barnes_distance_threshold

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(srcgrid%ncoords,1,mpi_integer,mpi_masternode,           &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(srcgrid%npasses,1,mpi_integer,mpi_masternode,           &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(srcgrid%neighbors,1,mpi_integer,mpi_masternode,         &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(srcgrid%distthresh,1,mpi_real,mpi_masternode,           &
         & mpi_comm_world,mpi_ierror)

    ! Define local variables

    call interpolation_initialize_grid(srcgrid)
    call interpolation_initialize_task_balance(srcgrid)

    ! If on master (root) node (task), define problem and broadcast
    ! variables to each slave (compute) node (task)
    
    if(mpi_procid .eq. mpi_masternode) then
       
       ! Define local variables

       rlon_min = minval(xlong)
       rlon_max = maxval(xlong)
       rlat_min = minval(xlat)
       rlat_max = maxval(xlat)

       ! Initialize local variables
       
       count = 1
       
       ! Loop through local variable
       
       do j = 1, ocean_hycom_ydim
          
          ! Loop through local variable

          do i = 1, ocean_hycom_xdim

             ! Define local variables

             srcgrid%xlong(count)  = xlong(i,j)*deg2rad
             srcgrid%xlat(count)   = xlat(i,j)*deg2rad
             srcgrid%depth(count)  = depth(i,j)
             srcgrid%rotang(count) = rotang(i,j)

             ! Update local variables
             
             count = count + 1
             
          end do ! do i = 1, ocean_hycom_xdim
          
       end do ! do j = 1, ocean_hycom_ydim

       ! Deallocate memory for local variables

       if(allocated(xlong))  deallocate(xlong)
       if(allocated(xlat))   deallocate(xlat)
       if(allocated(rotang)) deallocate(rotang)
       if(allocated(depth))  deallocate(depth)
       if(allocated(w))      deallocate(w)

    end if ! if(mpi_procid .eq. mpi_masternode)
       
    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(srcgrid%xlong,srcgrid%ncoords,mpi_real,mpi_masternode,  &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(srcgrid%xlat,srcgrid%ncoords,mpi_real,mpi_masternode,   &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(srcgrid%depth,srcgrid%ncoords,mpi_real,mpi_masternode,  &
         & mpi_comm_world,mpi_ierror)

    ! Compute local variables

    call interpolation_define_kdtree_mpi(srcgrid)
    call define_scaling_coefficients(srcgrid)

    !---------------------------------------------------------------------

    ! Compute local variables

    nlon = nint((rlon_max - rlon_min)/interp_dx) + 1
    nlat = nint((rlat_max - rlat_min)/interp_dy) + 1

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(nlon,1,mpi_integer,mpi_masternode,mpi_comm_world,       &
         & mpi_ierror)
    call mpi_bcast(nlat,1,mpi_integer,mpi_masternode,mpi_comm_world,       &
         & mpi_ierror)

    ! Define local variables

    dstgrid%ncoords    = nlon*nlat
    dstgrid%npasses    = barnes_npasses
    dstgrid%neighbors  = barnes_nneighbors
    dstgrid%distthresh = barnes_distance_threshold

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(dstgrid%ncoords,1,mpi_integer,mpi_masternode,           &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(dstgrid%npasses,1,mpi_integer,mpi_masternode,           &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(dstgrid%neighbors,1,mpi_integer,mpi_masternode,         &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(dstgrid%distthresh,1,mpi_real,mpi_masternode,           &
         & mpi_comm_world,mpi_ierror)

    ! Define local variables

    call interpolation_initialize_grid(dstgrid)
    call interpolation_initialize_task_balance(dstgrid)

    ! Allocate memory for local variables

    if(.not. allocated(static_weights))                                    &
         & allocate(static_weights(dstgrid%ncoords,dstgrid%neighbors,      &
         & dstgrid%npasses))

    ! If on master (root) node (task), define problem and broadcast
    ! variables to each slave (compute) node (task)
    
    if(mpi_procid .eq. mpi_masternode) then

       ! Initialize local variable
       
       count = 1
       
       ! Loop through local variable
       
       do j = 1, nlat
          
          ! Loop through local variable

          do i = 1, nlon

             ! Compute local variables

             dstgrid%xlong(count) = (rlon_min + (i-1)*interp_dx)*deg2rad
             dstgrid%xlat(count)  = (rlat_min + (j-1)*interp_dy)*deg2rad

             ! Update local variables

             count = count + 1

          end do ! do i = 1, nlon

       end do ! do j = 1, nlat

    end if ! if(mpi_procid .eq. mpi_masternode)
       
    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(dstgrid%xlong,dstgrid%ncoords,mpi_real,mpi_masternode,  &
         & mpi_comm_world,mpi_ierror)
    call mpi_bcast(dstgrid%xlat,dstgrid%ncoords,mpi_real,mpi_masternode,   &
         & mpi_comm_world,mpi_ierror)

    ! Compute local variables

    call interpolation_define_kdtree_mpi(dstgrid)
    call define_scaling_coefficients(dstgrid)
    call interpolation_define_weights_mpi(srcgrid,dstgrid)
    
    ! Allocate memory for local variables

    if(.not. allocated(srcgrid_var))                                       &
         & allocate(srcgrid_var(srcgrid%ncoords))
    if(.not. allocated(dstgrid_var))                                       &
         & allocate(dstgrid_var(dstgrid%ncoords))   

    ! Define local variable

    srcgrid_var = srcgrid%depth

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,mpi_masternode,    &
         & mpi_comm_world,mpi_ierror)

    ! Compute local variable

    call interpolation_nearest_neighbor_mpi(srcgrid,dstgrid)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,mpi_masternode,    &
         & mpi_comm_world,mpi_ierror)

    ! Define local variable

    dstgrid%depth = dstgrid_var

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(dstgrid%depth,dstgrid%ncoords,mpi_real,mpi_masternode,  &
         & mpi_comm_world,mpi_ierror)

    ! Loop through local variable

    do i = 1, srcgrid%ncoords

       ! Define local variable

       srcgrid_var(i) = maxval(hycom%temp(i,1:ocean_hycom_zdim))

    end do ! do i = 1, srcgrid%ncoords

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,mpi_masternode,    &
         & mpi_comm_world,mpi_ierror)

    ! Compute local variable

    call interpolation_nearest_neighbor_mpi(srcgrid,dstgrid)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,mpi_masternode,    &
         & mpi_comm_world,mpi_ierror)

    ! Define local variable

    dstgrid%temperature = dstgrid_var
    where(dstgrid%temperature .gt. 1.e20) dstgrid%temperature = spval

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(dstgrid%temperature,dstgrid%ncoords,mpi_real,           &
         & mpi_masternode,mpi_comm_world,mpi_ierror)

    ! Deallocate memory for local variables

    if(allocated(srcgrid_var)) deallocate(srcgrid_var)
    if(allocated(dstgrid_var)) deallocate(dstgrid_var)

    ! Loop through local variable

    do i = 1, dstgrid%neighbors

       ! Loop through local variable

       do j = 1, dstgrid%npasses

          ! Define local variable

          where(dstgrid%depth .gt. 1.e20) dstgrid%weights(:,i,j) = 0.0

       end do ! do j = 1, dstgrid%npasses

    end do ! do i = 1, dstgrid%neighbors

    ! Define local variable

    static_weights = dstgrid%weights

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)

    call mpi_bcast(dstgrid%weights,(dstgrid%ncoords*dstgrid%neighbors*     &
         & dstgrid%npasses),mpi_real,mpi_masternode,mpi_comm_world,        &
         & mpi_ierror)
    call mpi_bcast(static_weights,(dstgrid%ncoords*dstgrid%neighbors*      &
         & dstgrid%npasses),mpi_real,mpi_masternode,mpi_comm_world,        &
         & mpi_ierror)

    !=====================================================================

  end subroutine hycomtograds_interface_initialize

  !=======================================================================

  ! hycomtograds_interface_cleanup.f90:

  !-----------------------------------------------------------------------

  subroutine hycomtograds_interface_cleanup()

    !=====================================================================

    ! Deallocate memory for local variables
    
    call interpolation_cleanup_grid(srcgrid)
    call interpolation_cleanup_task_balance(srcgrid)
    call interpolation_cleanup_grid(dstgrid)
    call interpolation_cleanup_task_balance(dstgrid)
    if(allocated(static_weights)) deallocate(static_weights)

    !=====================================================================

  end subroutine hycomtograds_interface_cleanup
  
  !=======================================================================

end module hycomtograds_interface
