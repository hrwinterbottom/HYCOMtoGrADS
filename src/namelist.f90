module namelist

  !=======================================================================

  use kinds

  !-----------------------------------------------------------------------

  use mpi_interface

  !-----------------------------------------------------------------------

  implicit none
  
  !-----------------------------------------------------------------------

  ! Define global variables

  character(len=500)                             :: hycom_filename            = 'NOT USED'
  character(len=500)                             :: regional_grid_filename    = 'NOT USED'
  character(len=500)                             :: regional_depth_filename   = 'NOT USED'
  logical                                        :: debug                     = .false.
  logical                                        :: is_steric                 = .false.
  logical                                        :: is_sst                    = .false.
  logical                                        :: is_ssh                    = .false.
  logical                                        :: is_tchp                   = .false.
  logical                                        :: is_ohc                    = .false.
  logical                                        :: is_omld                   = .false.
  logical                                        :: is_omlt                   = .false.
  logical                                        :: is_omls                   = .false.
  logical                                        :: is_omlu                   = .false.
  logical                                        :: is_omlv                   = .false.
  logical                                        :: is_omlke                  = .false.
  logical                                        :: is_uvel                   = .false.
  logical                                        :: is_vvel                   = .false.
  logical                                        :: is_temp                   = .false.
  logical                                        :: is_salin                  = .false.
  logical                                        :: is_dens                   = .false.
  logical                                        :: is_dpth                   = .false.
  logical                                        :: is_isotdpth               = .false.
  logical                                        :: is_zinterp                = .false.
  logical                                        :: is_isointerp              = .false.
  real(r_kind)                                   :: zlevs(100)                = -999.0
  real(r_kind)                                   :: tlevs(100)                = -999.0
  real(r_kind)                                   :: hycom_time                = 10000.0 
  real(r_kind)                                   :: interp_dx                 = 1.0
  real(r_kind)                                   :: interp_dy                 = 1.0
  real(r_kind)                                   :: barnes_weights_threshold  = 1.e-5
  real(r_kind)                                   :: barnes_distance_threshold = 300000.0
  real(r_kind)                                   :: depth_integral_dz         = 10.0
  integer                                        :: ocean_hycom_xdim
  integer                                        :: ocean_hycom_ydim
  integer                                        :: ocean_hycom_zdim
  integer                                        :: barnes_nneighbors         = 10
  integer                                        :: barnes_npasses            = 2
  integer                                        :: nlev                      = 0
  integer                                        :: namelist_io_error         = 0
  namelist /fileio/      debug, hycom_filename, is_steric
  namelist /variableio/  is_sst, is_ssh, is_tchp, is_ohc, is_omld,       &
       & is_omlt, is_omls, is_omlu, is_omlv, is_omlke, is_uvel, is_vvel, &
       & is_temp, is_salin, is_dens, is_dpth, is_isotdpth,               &
       & depth_integral_dz
  namelist /interpio/    interp_dx, interp_dy, barnes_nneighbors,        &
       & barnes_npasses, barnes_distance_threshold,                      &
       & barnes_weights_threshold
  namelist /hycomio/     regional_grid_filename,                         &
       & regional_depth_filename, hycom_time, ocean_hycom_xdim,          &
       & ocean_hycom_ydim, ocean_hycom_zdim
  namelist /zinterpio/   is_zinterp, zlevs
  namelist /isointerpio/ is_isointerp, tlevs

  !---------------------------------------------------------------------

contains

  !=======================================================================

  ! namelistparams.f90:

  !-----------------------------------------------------------------------

  subroutine namelistparams()
    
    ! Define variables computed within routine

    logical                                        :: is_it_there
    integer                                        :: unit_nml

    ! Define counting variables

    integer                                        :: i, j, k
    
    !=====================================================================

    ! Define local variables

    unit_nml    = 9
    is_it_there = .false.
    inquire(file='hycomtograds.input',exist = is_it_there)

    ! Check local variable and proceed accordingly

    if(is_it_there) then

       ! Define local variables

       open(file = 'hycomtograds.input',                                  &
            unit = unit_nml        ,                                      &
            status = 'old'         ,                                      &
            form = 'formatted'     ,                                      &
            action = 'read'        ,                                      &
            access = 'sequential'  )
       read(unit_nml,NML = fileio)
       read(unit_nml,NML = variableio)
       read(unit_nml,NML = interpio)
       read(unit_nml,NML = hycomio)
       read(unit_nml,NML = zinterpio)
       read(unit_nml,NML = isointerpio)
       close(unit_nml)

       ! Check local variable and proceed accordingly

       if(is_zinterp) then

          ! Initialize local variables

          nlev = 0

          ! Loop through all user specified vertical levels

          do k = 1, 100

             ! Update local variable accordingly

             if(zlevs(k) .ne. -999.0) then
                
                ! Update local variables

                nlev = nlev + 1
          
             end if ! if(zlevs(k) .ne. -999.0)

          end do ! do k = 1, 100

       end if ! if(is_zinterp)

       ! Check local variable and proceed accordingly

       if(is_isointerp) then

          ! Initialize local variables

          nlev = 0

          ! Loop through all user specified vertical levels

          do k = 1, 100

             ! Update local variable accordingly

             if(tlevs(k) .ne. -999.0) then
                
                ! Update local variables

                nlev = nlev + 1
          
             end if ! if(tlevs(k) .ne. -999.0)

          end do ! do k = 1, 100

       end if ! if(is_isointerp)

    end if ! if(is_it_there)

    ! Check local variable and proceed accordingly

    if(.not. is_it_there) then 

       ! Print message to user

       if(mpi_procid .eq. mpi_masternode) write(6,500)

       ! Enable the root task to catch up from I/O and calculations
    
       call mpi_barrier(mpi_comm_world,mpi_ierror)

       ! Update local variables

       namelist_io_error = 1.0

    end if ! if(.not. is_it_there)

    !---------------------------------------------------------------------

    ! Check local variable and proceed accordingly

    if(is_ohc) then

       ! Check local variable and proceed accordingly

       if(.not. is_zinterp) then

          ! Print message to user

          if(mpi_procid .eq. mpi_masternode) write(6,501)

          ! Enable the root task to catch up from I/O and calculations
    
          call mpi_barrier(mpi_comm_world,mpi_ierror)

          ! Abort routine

          call mpi_interface_terminate()
          stop

       end if ! if(.not. is_zinterp)

    end if ! if(is_ohc)

    ! Check local variable and proceed accordingly

    if(is_isotdpth) then

       ! Check local variable and proceed accordingly

       if(.not. is_isointerp) then

          ! Print message to user

          if(mpi_procid .eq. mpi_masternode) write(6,502)

          ! Enable the root task to catch up from I/O and calculations
    
          call mpi_barrier(mpi_comm_world,mpi_ierror)

          ! Abort routine

          call mpi_interface_terminate()
          stop

       end if ! if(.not. is_zinterp)

    end if ! if(is_isotdpth)

    !---------------------------------------------------------------------

    ! If on master (root) node (task), define problem and broadcast
    ! variables to each slave (compute) node (task)

    if(mpi_procid .eq. mpi_masternode) then

       ! Print message to user

       write(6,*) '&FILEIO'
       write(6,*) 'DEBUG             = ', debug
       write(6,*) 'HYCOM_FILENAME    = ',                                 &
            & trim(adjustl(hycom_filename))
       write(6,*) 'IS_STERIC         = ', is_steric
       write(6,*) '/'
       write(6,*) '&VARIABLEIO'
       if(is_sst)      write(6,*) 'IS_SST                    = ',         &
            & is_sst
       if(is_ssh)      write(6,*) 'IS_SSH                    = ',         &
            & is_ssh
       if(is_tchp)     write(6,*) 'IS_TCHP                   = ',         &
            & is_tchp
       if(is_ohc)      write(6,*) 'IS_OHC                    = ',         &
            & is_ohc
       if(is_omld)     write(6,*) 'IS_OMLD                   = ',         &
            & is_omld
       if(is_omlt)     write(6,*) 'IS_OMLT                   = ',         &
            & is_omlt
       if(is_omls)     write(6,*) 'IS_OMLS                   = ',         &
            & is_omls
       if(is_omlu)     write(6,*) 'IS_OMLU                   = ',         &
            & is_omlu
       if(is_omlv)     write(6,*) 'IS_OMLV                   = ',         &
            & is_omlv
       if(is_omlke)    write(6,*) 'IS_OMLKE                  = ',         &
            & is_omlke
       if(is_uvel)     write(6,*) 'IS_UVEL                   = ',         &
            & is_uvel
       if(is_vvel)     write(6,*) 'IS_VVEL                   = ',         &
            & is_vvel
       if(is_temp)     write(6,*) 'IS_TEMP                   = ',         &
            & is_temp
       if(is_salin)    write(6,*) 'IS_SALIN                  = ',         &
            & is_salin
       if(is_dens)     write(6,*) 'IS_DENS                   = ',         &
	    & is_dens
       if(is_dpth)     write(6,*) 'IS_DPTH                   = ',         &
	    & is_dpth
       if(is_isotdpth) write(6,*) 'IS_ISOTDPTH               = ',         &
            & is_isotdpth
       write(6,*) 'DEPTH_INTEGRAL_DZ         = ', depth_integral_dz
       write(6,*) '/'
       write(6,*) '&INTERPIO'
       write(6,*) 'INTERP_DX                 = ', interp_dx
       write(6,*) 'INTERP_DX                 = ', interp_dy
       write(6,*) 'BARNES_NNEIGHBORS         = ', barnes_nneighbors
       write(6,*) 'BARNES_NPASSES            = ', barnes_npasses
       write(6,*) 'BARNES_DISTANCE_THRESHOLD = ',                         &
            & barnes_distance_threshold
       write(6,*) 'BARNES_WEIGHTS_THRESHOLD  = ',                         &
            & barnes_weights_threshold
       write(6,*) '/'
       write(6,*) '&HYCOMIO'
       write(6,*) 'REGIONAL_GRID_FILENAME    = ',                         &
            & trim(adjustl(regional_grid_filename))
       write(6,*) 'REGIONAL_DEPTH_FILENAME   = ',                         &
            & trim(adjustl(regional_depth_filename))
       write(6,*) 'OCEAN_HYCOM_XDIM          = ', ocean_hycom_xdim
       write(6,*) 'OCEAN_HYCOM_YDIM          = ', ocean_hycom_ydim
       write(6,*) 'OCEAN_HYCOM_ZDIM          = ', ocean_hycom_zdim
       write(6,*) 'HYCOM_TIME                = ', hycom_time
       write(6,*) '/'
       if(is_zinterp) then
          write(6,*) '&ZINTERPIO'
          write(6,*) 'IS_ZINTERP   = ', is_zinterp
          write(6,*) 'ZLEVS        = ', zlevs(1:nlev)
          write(6,*) '/'
       end if ! if(is_zinterp)
       if(is_isointerp) then
          write(6,*) '&ISOINTERPIO'
          write(6,*) 'IS_ISOINTERP = ', is_isointerp
          write(6,*) 'TLEVS        = ', tlevs(1:nlev)
          write(6,*) '/'
       end if ! if(is_isointerp)
       write(6,*) ' '

    end if ! if(mpi_procid .eq. mpi_masternode)

    ! Enable the root task to catch up from I/O and calculations
    
    call mpi_barrier(mpi_comm_world,mpi_ierror)

    !=====================================================================

    ! Define format statements

500 format('NAMELISTPARAMS: hycomtograds.input not found in current ',    &
         & 'working directory. ABORTING!!!!')
501 format('NAMELISTPARAMS: In order to compute the ocean heat content ', &
         & '(OHC), the user must set is_zinterp = .true. and provide ',   &
         & 'the isobathes for interpolation. ABORTING!!!!')
502 format('NAMELISTPARAMS: In order to compute the depth for ',          &
         & 'an isotherm(s), the user must set is_isointerp = .true. ',    &
         & 'and provide the isotherms for interpolation. ABORTING!!!!')

    !=====================================================================

  end subroutine namelistparams
  
  !=======================================================================

end module namelist
