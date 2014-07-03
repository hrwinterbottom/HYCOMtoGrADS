module variable_interface

  !=======================================================================
  
  ! Define associated modules and subroutines
  
  !-----------------------------------------------------------------------
  
  use constants
  use kinds

  !-----------------------------------------------------------------------

  use interpolation_interface
  use mpi_interface
  use namelist

  !-----------------------------------------------------------------------
  
  implicit none

  !-----------------------------------------------------------------------

  ! Define all data and structure types for routine; these variables
  ! are variables required by the subroutines within this module

  type hycomvariables
     real(r_kind),                 dimension(:,:),   allocatable :: uvel
     real(r_kind),                 dimension(:,:),   allocatable :: vvel
     real(r_kind),                 dimension(:,:),   allocatable :: thknss
     real(r_kind),                 dimension(:,:),   allocatable :: temp
     real(r_kind),                 dimension(:,:),   allocatable :: salin
     real(r_kind),                 dimension(:,:),   allocatable :: depth
     real(r_kind),                 dimension(:),     allocatable :: montg1
     real(r_kind),                 dimension(:),     allocatable :: srfhgt
     real(r_kind),                 dimension(:),     allocatable :: steric
     real(r_kind),                 dimension(:),     allocatable :: surflx
     real(r_kind),                 dimension(:),     allocatable :: salflx
     real(r_kind),                 dimension(:),     allocatable :: bl_dpth
     real(r_kind),                 dimension(:),     allocatable :: mix_dpth
     real(r_kind),                 dimension(:),     allocatable :: covice
     real(r_kind),                 dimension(:),     allocatable :: thkice
     real(r_kind),                 dimension(:),     allocatable :: temice
     real(r_kind),                 dimension(:),     allocatable :: u_btrop
     real(r_kind),                 dimension(:),     allocatable :: v_btrop
  end type hycomvariables

  type variable_info
     character(len=50)                                           :: grads_string
     character(len=10)                                           :: grads_id
     real(r_kind),                dimension(:),      allocatable :: array
     integer                                                     :: zdim
  end type variable_info

  ! Define global variables

  type(interpgrid)                                               :: srcgrid
  type(interpgrid)                                               :: dstgrid
  real(r_kind),                   dimension(:,:,:),  allocatable :: static_weights
  real(r_kind)                                                   :: spval             = 1.e30
  integer                                                        :: var_counter       = 0
  integer                                                        :: number_of_var

  !-----------------------------------------------------------------------

contains

  !=======================================================================

  ! variable_interface_process.f90:

  !-----------------------------------------------------------------------

  subroutine variable_interface_process(hycom)

    ! Define variables passed to routine

    type(hycomvariables)                                           :: hycom

    ! Define variables computed within routine

    type(variable_info),           dimension(:),       allocatable :: variable_info_grid

    ! Define counting variables

    integer                                                        :: i, j, k, l

    !=====================================================================

    ! Define local variables

    call variable_interface_num_var()

    ! Allocate memory for local variables

    call variable_interface_initialize(variable_info_grid)

    ! Define local variables

    call variable_interface_define(variable_info_grid)

    !---------------------------------------------------------------------

    ! Loop through local variable

    do l = 1, number_of_var

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)

       call mpi_bcast(variable_info_grid(l)%grads_id,10,mpi_character,     &
            & mpi_masternode,mpi_comm_world,mpi_ierror)
       call mpi_bcast(variable_info_grid(l)%zdim,1,mpi_integer,            &
            & mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Define local variable accordingly

       nlev = max(variable_info_grid(l)%zdim,nlev)

       ! Compute local variable
                             
       call variable_interface_compute(hycom,variable_info_grid(l))

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)

       call mpi_bcast(nlev,1,mpi_integer,mpi_masternode,mpi_comm_world,    &
           & mpi_ierror)

    end do ! do l = 1, number_of_var

    !---------------------------------------------------------------------

    ! If on master (root) node (task), define problem and broadcast
    ! variables to each slave (compute) node (task)
    
    if(mpi_procid .eq. mpi_masternode) then

       ! Create external files

       call grads_interface_descriptor(variable_info_grid)
       call grads_interface_binary(variable_info_grid)

    end if ! if(mpi_procid .eq. mpi_masternode)
       
    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! Deallocate memory for local variables

    call variable_interface_cleanup(variable_info_grid)

    !=====================================================================

  end subroutine variable_interface_process

  !=======================================================================

  ! variable_interface_compute.f90:

  !-----------------------------------------------------------------------

  subroutine variable_interface_compute(hycom,grid)

    ! Define variables passed to routine

    type(hycomvariables)                                           :: hycom
    type(variable_info)                                            :: grid

    !=====================================================================

    ! Check local variable and proceed accordingly

    if(grid%grads_id .eq. 'sst') then

       ! Allocate memory for local variable

       if(.not. allocated(grid%array))                                     &
            & allocate(grid%array(dstgrid%ncoords))

       ! Compute local variable

       call calculate_seasurfacetemperature(hycom,grid)

    end if ! if(grid%grads_id .eq. 'sst')

    ! Check local variable and proceed accordingly

    if(grid%grads_id .eq. 'ssh') then

       ! Allocate memory for local variable

       if(.not. allocated(grid%array))                                     &
            & allocate(grid%array(dstgrid%ncoords))

       ! Compute local variable

       call calculate_seasurfaceheight(hycom,grid)

    end if ! if(grid%grads_id .eq. 'ssh')

    ! Check local variable and proceed accordingly

    if(grid%grads_id .eq. 'tchp') then

       ! Allocate memory for local variable

       if(.not. allocated(grid%array))                                     &
            & allocate(grid%array(dstgrid%ncoords))

       ! Compute local variable

       call calculate_tropicalcycloneheatpotential(hycom,grid)

    end if ! if(grid%grads_id .eq. 'tchp')

    ! Check local variable and proceed accordingly

    if(grid%grads_id .eq. 'omld') then

       ! Allocate memory for local variable

       if(.not. allocated(grid%array))                                     &
            & allocate(grid%array(dstgrid%ncoords))

       ! Compute local variable

       call calculate_oceanmixedlayerdepth(hycom,grid)

    end if ! if(grid%grads_id .eq. 'omld')

    ! Check local variable and proceed accordingly

    if(grid%grads_id .eq. 'omlt') then

       ! Allocate memory for local variable

       if(.not. allocated(grid%array))                                     &
            & allocate(grid%array(dstgrid%ncoords))

       ! Compute local variable

       call calculate_oceanmixedlayertemperature(hycom,grid)

    end if ! if(grid%grads_id .eq. 'omlt')

    ! Check local variable and proceed accordingly

    if(grid%grads_id .eq. 'omls') then

       ! Allocate memory for local variable

       if(.not. allocated(grid%array))                                     &
            & allocate(grid%array(dstgrid%ncoords))

       ! Compute local variable

       call calculate_oceanmixedlayersalinity(hycom,grid)

    end if ! if(grid%grads_id .eq. 'omls')

    ! Check local variable and proceed accordingly

    if(grid%grads_id .eq. 'omlu') then

       ! Allocate memory for local variable

       if(.not. allocated(grid%array))                                     &
            & allocate(grid%array(dstgrid%ncoords))

       ! Compute local variable

       call calculate_oceanmixedlayeruvelocity(hycom,grid)

    end if ! if(grid%grads_id .eq. 'omlu')

    ! Check local variable and proceed accordingly

    if(grid%grads_id .eq. 'omlv') then

       ! Allocate memory for local variable

       if(.not. allocated(grid%array))                                     &
            & allocate(grid%array(dstgrid%ncoords))

       ! Compute local variable

       call calculate_oceanmixedlayervvelocity(hycom,grid)

    end if ! if(grid%grads_id .eq. 'omlv')

    ! Check local variable and proceed accordingly

    if(grid%grads_id .eq. 'omlke') then

       ! Allocate memory for local variable

       if(.not. allocated(grid%array))                                     &
            & allocate(grid%array(dstgrid%ncoords))

       ! Compute local variable

       call calculate_oceanmixedlayerkineticenergy(hycom,grid)

    end if ! if(grid%grads_id .eq. 'omlke')

    ! Check local variable and proceed accordingly

    if(grid%grads_id .eq. 'ohc') then

       ! Allocate memory for local variable

       if(.not. allocated(grid%array))                                     &
            & allocate(grid%array(dstgrid%ncoords*grid%zdim))

       ! Compute local variable

       call calculate_oceanheatcontent(hycom,grid)

    end if ! if(grid%grads_id .eq. 'ohc')

    ! Check local variable and proceed accordingly

    if(grid%grads_id .eq. 'isotdpth') then

       ! Allocate memory for local variable

       if(.not. allocated(grid%array))                                     &
            & allocate(grid%array(dstgrid%ncoords*grid%zdim))

       ! Compute local variable

       call calculate_isotdpth(hycom,grid)

    end if ! if(grid%grads_id .eq. 'isotdpth')

    ! Check local variable and proceed accordingly

    if(grid%grads_id .eq. 'uvel') then

       ! Allocate memory for local variable

       if(.not. allocated(grid%array))                                     &
            & allocate(grid%array(dstgrid%ncoords*grid%zdim))

       ! Compute local variable

       call calculate_oceanuvelocity(hycom,grid)

    end if ! if(grid%grads_id .eq. 'uvel')

    ! Check local variable and proceed accordingly

    if(grid%grads_id .eq. 'vvel') then

       ! Allocate memory for local variable

       if(.not. allocated(grid%array))                                     &
            & allocate(grid%array(dstgrid%ncoords*grid%zdim))

       ! Compute local variable

       call calculate_oceanvvelocity(hycom,grid)

    end if ! if(grid%grads_id .eq. 'vvel')

    ! Check local variable and proceed accordingly

    if(grid%grads_id .eq. 'temp') then

       ! Allocate memory for local variable

       if(.not. allocated(grid%array))                                     &
            & allocate(grid%array(dstgrid%ncoords*grid%zdim))

       ! Compute local variable

       call calculate_oceantemperature(hycom,grid)

    end if ! if(grid%grads_id .eq. 'temp')

    ! Check local variable and proceed accordingly

    if(grid%grads_id .eq. 'salin') then

       ! Allocate memory for local variable

       if(.not. allocated(grid%array))                                     &
            & allocate(grid%array(dstgrid%ncoords*grid%zdim))

       ! Compute local variable

       call calculate_oceansalinity(hycom,grid)

    end if ! if(grid%grads_id .eq. 'salin')

    ! Check local variable and proceed accordingly

    if(grid%grads_id .eq. 'dens') then

       ! Allocate memory for local variable

       if(.not. allocated(grid%array))                                     &
            & allocate(grid%array(dstgrid%ncoords*grid%zdim))

       ! Compute local variable

       call calculate_oceandensity(hycom,grid)

    end if ! if(grid%grads_id .eq. 'dens')

    ! Check local variable and proceed accordingly

    if(grid%grads_id .eq. 'dpth') then

       ! Allocate memory for local variable

       if(.not. allocated(grid%array))                                     &
            & allocate(grid%array(dstgrid%ncoords*grid%zdim))

       ! Compute local variable

       call calculate_oceandepth(hycom,grid)

    end if ! if(grid%grads_id .eq. 'dpth')

    !=====================================================================

  end subroutine variable_interface_compute

  !=======================================================================

  ! variable_interface_define.f90:

  !-----------------------------------------------------------------------

  subroutine variable_interface_define(grid)

    ! Define variable returned by routine

    type(variable_info),     dimension(number_of_var)              :: grid

    !=====================================================================

    ! Define local variables accordingly

    if(is_sst) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%grads_string = 'sea-surface temperature (K)'
       grid(var_counter)%grads_id     = 'sst'
       grid(var_counter)%zdim         = 1

    end if ! if(is_sst)

    ! Define local variables accordingly

    if(is_ssh) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%grads_string = 'sea-surface height (m)'
       grid(var_counter)%grads_id     = 'ssh'
       grid(var_counter)%zdim         = 1

    end if ! if(is_ssh)

    ! Define local variables accordingly

    if(is_tchp) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%grads_string = 'tropical cyclone heat potential (J/m^2)'
       grid(var_counter)%grads_id     = 'tchp'
       grid(var_counter)%zdim         = 1

    end if ! if(is_tchp)

    ! Define local variables accordingly

    if(is_omld) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%grads_string = 'ocean mixed-layer depth (m)'
       grid(var_counter)%grads_id     = 'omld'
       grid(var_counter)%zdim         = 1

    end if ! if(is_omld)

    ! Define local variables accordingly

    if(is_omlt) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%grads_string = 'ocean mixed-layer mean potential temperature (K)'
       grid(var_counter)%grads_id     = 'omlt'
       grid(var_counter)%zdim         = 1

    end if ! if(is_omlt)

    ! Define local variables accordingly

    if(is_omls) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%grads_string = 'ocean mixed-layer mean salinity (PSU)'
       grid(var_counter)%grads_id     = 'omls'
       grid(var_counter)%zdim         = 1

    end if ! if(is_omls)

    ! Define local variables accordingly

    if(is_omlu) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%grads_string = 'ocean mixed-layer mean zonal velocity (m/s)'
       grid(var_counter)%grads_id     = 'omlu'
       grid(var_counter)%zdim         = 1

    end if ! if(is_omlu)

    ! Define local variables accordingly

    if(is_omlv) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%grads_string = 'ocean mixed-layer mean meridional velocity (m/s)'
       grid(var_counter)%grads_id     = 'omlv'
       grid(var_counter)%zdim         = 1

    end if ! if(is_omlv)

    ! Define local variables accordingly

    if(is_omlke) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%grads_string = 'ocean mixed-layer mean kinetic energy (m^2/s^2)'
       grid(var_counter)%grads_id     = 'omlke'
       grid(var_counter)%zdim         = 1

    end if ! if(is_omlke)

    ! Define local variables accordingly

    if(is_ohc) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%grads_string = 'ocean heat content (J/m^2)'
       grid(var_counter)%grads_id     = 'ohc'
       grid(var_counter)%zdim         = nlev

    end if ! if(is_ohc)

    ! Define local variables accordingly

    if(is_isotdpth) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%grads_string = 'isotherm depth (m)'
       grid(var_counter)%grads_id     = 'isotdpth'
       grid(var_counter)%zdim         = nlev

    end if ! if(is_isotdpth)

    ! Define local variables accordingly

    if(is_uvel) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%grads_string = 'zonal velocity (m/s)'
       grid(var_counter)%grads_id     = 'uvel'
       grid(var_counter)%zdim         = nlev

    end if ! if(is_uvel)

    ! Define local variables accordingly

    if(is_vvel) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%grads_string = 'meridional velocity (m/s)'
       grid(var_counter)%grads_id     = 'vvel'
       grid(var_counter)%zdim         = nlev

    end if ! if(is_vvel)

    ! Define local variables accordingly

    if(is_temp) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%grads_string = 'sea-water potential temperature (K)'
       grid(var_counter)%grads_id     = 'temp'
       grid(var_counter)%zdim         = nlev

    end if ! if(is_temp)

    ! Define local variables accordingly

    if(is_salin) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%grads_string = 'sea-water salinity (PSU)'
       grid(var_counter)%grads_id     = 'salin'
       grid(var_counter)%zdim         = nlev

    end if ! if(is_salin)

    ! Define local variables accordingly

    if(is_dens) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%grads_string = 'sea-water density (kg/m^3)'
       grid(var_counter)%grads_id     = 'dens'
       grid(var_counter)%zdim         = nlev

    end if ! if(is_dens)

    ! Define local variables accordingly

    if(is_dpth) then

       ! Update counting variable

       var_counter = var_counter + 1

       ! Define local variables

       grid(var_counter)%grads_string = 'depth (m)'
       grid(var_counter)%grads_id     = 'dpth'
       grid(var_counter)%zdim         = nlev

    end if ! if(is_dpth)

    !=====================================================================

  end subroutine variable_interface_define

  !=======================================================================

  ! variable_interface_num_var.f90:

  !-----------------------------------------------------------------------

  subroutine variable_interface_num_var()

    !=====================================================================

    ! Initialize local variable

    number_of_var = 0

    !---------------------------------------------------------------------

    ! Update local variable accordingly

    if(is_sst)      number_of_var = number_of_var + 1
    if(is_ssh)      number_of_var = number_of_var + 1
    if(is_tchp)     number_of_var = number_of_var + 1
    if(is_ohc)      number_of_var = number_of_var + 1
    if(is_isotdpth) number_of_var = number_of_var + 1
    if(is_omld)     number_of_var = number_of_var + 1
    if(is_omlt)     number_of_var = number_of_var + 1
    if(is_omls)     number_of_var = number_of_var + 1
    if(is_omlu)     number_of_var = number_of_var + 1
    if(is_omlv)     number_of_var = number_of_var + 1
    if(is_omlke)    number_of_var = number_of_var + 1
    if(is_uvel)     number_of_var = number_of_var + 1
    if(is_vvel)     number_of_var = number_of_var + 1
    if(is_temp)     number_of_var = number_of_var + 1
    if(is_salin)    number_of_var = number_of_var + 1
    if(is_dens)     number_of_var = number_of_var + 1
    if(is_dpth)     number_of_var = number_of_var + 1

    !=====================================================================

  end subroutine variable_interface_num_var

  !=======================================================================

  ! variable_interface_initialize.f90:

  !-----------------------------------------------------------------------

  subroutine variable_interface_initialize(grid)

    ! Define variables computed within routine

    type(variable_info),     dimension(:),       allocatable :: grid

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(grid)) allocate(grid(number_of_var))

    !=====================================================================

  end subroutine variable_interface_initialize

  !=======================================================================

  ! variable_interface_cleanup.f90:

  !-----------------------------------------------------------------------

  subroutine variable_interface_cleanup(grid)

    ! Define variables computed within routine

    type(variable_info),     dimension(:),       allocatable :: grid

    !=====================================================================

    ! Allocate memory for local variables

    if(allocated(grid)) deallocate(grid)

    !=====================================================================

  end subroutine variable_interface_cleanup

  !=======================================================================

  ! calculate_isotherminterpolation.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_isotherminterpolation(hycom,grid,isotherm,          &
       & grid_isotherm)

    ! Define variables passed to routine

    type(hycomvariables)                                                  :: hycom
    real(r_kind),             dimension(srcgrid%ncoords,ocean_hycom_zdim) :: grid
    real(r_kind)                                                          :: isotherm

    ! Define variables returned by routine

    real(r_kind),             dimension(srcgrid%ncoords)                  :: grid_isotherm

    ! Define variables computed within routine

    real(r_kind),             dimension(srcgrid%ncoords)                  :: workgrid
    real(r_kind)                                                          :: isotherm_depth

    ! Define counting variables

    integer                                                               :: i, j, k

    !=====================================================================

    ! Initialize local variable

    grid_isotherm = 0.0
    workgrid      = 0.0

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(grid,(srcgrid%ncoords*ocean_hycom_zdim),mpi_real,       &
         & mpi_masternode,mpi_comm_world,mpi_ierror)

    ! If on slave (compute) node (task), receive variables, compute
    ! variables, and send variables to master (root) node (task)

    if(mpi_procid .ne. mpi_masternode) then

       ! Loop through local variable and proceed accordingly
       
       do i = srcgrid%mpi_count_begin(mpi_procid),                         &
            & srcgrid%mpi_count_end(mpi_procid)
   
          ! Compute local variable

          call calculate_isothermdepth(ocean_hycom_zdim,                   &
               & hycom%temp(i,1:ocean_hycom_zdim),                         &
               & hycom%depth(i,1:ocean_hycom_zdim),isotherm,               &
               & isotherm_depth)

          ! Check local variable and proceed accordingly

          if(isotherm_depth .le. 0.0 .or.                                  &
               & maxval(hycom%temp(i,1:ocean_hycom_zdim)) .gt. 1.e20       &
               & .or. minval(hycom%temp(i,1:ocean_hycom_zdim)) .ge.        &
               & isotherm) isotherm_depth = spval

          ! Check local variable and proceed accordingly

          if(isotherm_depth .ne. spval) then

             ! Compute local variable

             call calculate_variableatdepth(ocean_hycom_zdim,              &
                  & hycom%depth(i,1:ocean_hycom_zdim),                     &
                  & grid(i,1:ocean_hycom_zdim),isotherm_depth,             &
                  & workgrid(i))

          else   ! if(isotherm_depth .ne. spval)  

             ! Define local variable

             workgrid(i) = spval

          end if ! if(isotherm_depth .ne. spval)

       end do ! do i = srcgrid%mpi_count_begin(mpi_procid),                &
              ! srcgrid%mpi_count_end(mpi_procid)

    end if ! if(mpi_procid .ne. mpi_masternode)

    ! Define local variable

    call mpi_reduce(workgrid(1:srcgrid%ncoords),                           &
         & grid_isotherm(1:srcgrid%ncoords),srcgrid%ncoords,mpi_real,      &
         & mpi_sum,mpi_masternode,mpi_comm_world,mpi_ierror)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(grid_isotherm,srcgrid%ncoords,mpi_real,                 &
         & mpi_masternode,mpi_comm_world,mpi_ierror)

    !=====================================================================

  end subroutine calculate_isotherminterpolation

  !=======================================================================

  ! calculate_isobathyinterpolation.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_isobathyinterpolation(hycom,grid,isobathy,          &
       & grid_isobathy)

    ! Define variables passed to routine

    type(hycomvariables)                                                  :: hycom
    real(r_kind),             dimension(srcgrid%ncoords,ocean_hycom_zdim) :: grid
    real(r_kind)                                                          :: isobathy

    ! Define variables returned by routine

    real(r_kind),             dimension(srcgrid%ncoords)                  :: grid_isobathy

    ! Define variables computed within routine

    real(r_kind),             dimension(srcgrid%ncoords)                  :: workgrid
    real(r_kind)                                                          :: w1
    real(r_kind)                                                          :: w2
    integer                                                               :: level

    ! Define counting variables

    integer                                                               :: i, j, k

    !=====================================================================

    ! Initialize local variable

    grid_isobathy = 0.0
    workgrid      = 0.0

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(grid,(srcgrid%ncoords*ocean_hycom_zdim),mpi_real,       &
         & mpi_masternode,mpi_comm_world,mpi_ierror)

    ! If on slave (compute) node (task), receive variables, compute
    ! variables, and send variables to master (root) node (task)

    if(mpi_procid .ne. mpi_masternode) then

       ! Loop through local variable and proceed accordingly
       
       do i = srcgrid%mpi_count_begin(mpi_procid),                         &
            & srcgrid%mpi_count_end(mpi_procid)
    
          ! Check local variable

          if(hycom%depth(i,ocean_hycom_zdim) .gt. isobathy) then

             ! Initialize local variable

             level = 1

             ! Loop through local variable

             do k = 1, ocean_hycom_zdim - 1

                ! Check local variable and proceed accordingly

                if(hycom%depth(i,k) .le. isobathy .and.                    &
                     & hycom%depth(i,k+1) .ge. isobathy) then

                   ! Update local variable

                   level = k

                end if ! if(hycom%depth(i,k) .le. isobathy .and.           &
                       ! hycom%depth(i,k+1) .ge. isobathy)

             end do ! do k = 1, ocean_hycom_zdim - 1

             ! Compute local variables

             w1          = (hycom%depth(i,level+1) - isobathy)/            &
                  & (hycom%depth(i,level+1) - hycom%depth(i,level))
             w2          = 1.0 - w1
             workgrid(i) = w1*grid(i,level) + w2*grid(i,level+1)

          else   ! if(maxval(hycom%depth(i,k)) .le. isobathy)

             ! Define local variable

             workgrid(i) = spval

          end if ! if(maxval(hycom%depth(i,:)) .le. isobathy)

       end do ! do i = srcgrid%mpi_count_begin(mpi_procid),                &
              ! srcgrid%mpi_count_end(mpi_procid)

    end if ! if(mpi_procid .ne. mpi_masternode)

    ! Define local variable

    call mpi_reduce(workgrid(1:srcgrid%ncoords),                           &
         & grid_isobathy(1:srcgrid%ncoords),srcgrid%ncoords,mpi_real,      &
         & mpi_sum,mpi_masternode,mpi_comm_world,mpi_ierror)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(grid_isobathy,srcgrid%ncoords,mpi_real,                 &
         & mpi_masternode,mpi_comm_world,mpi_ierror)

    !=====================================================================

  end subroutine calculate_isobathyinterpolation

  !=======================================================================

  ! calculate_isothermdepth.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_isotdpth(hycom,grid)

    ! Define variables passed to routine

    type(hycomvariables)                                           :: hycom
    type(variable_info)                                            :: grid

    ! Define variables computed within routine

    real(r_kind),              dimension(:),           allocatable :: mpi_isotdpth
    real(r_kind)                                                   :: isotdpth

    ! Define counting variables

    integer                                                        :: i, j, k

    !=====================================================================

    ! Initialize local variables

    call init_constants_derived()

    !---------------------------------------------------------------------

    ! Allocate memory for local variables

    if(.not. allocated(srcgrid_var))                                       &
         & allocate(srcgrid_var(srcgrid%ncoords))
    if(.not. allocated(dstgrid_var))                                       &
         & allocate(dstgrid_var(dstgrid%ncoords))
    if(.not. allocated(mpi_isotdpth))                                      &
         & allocate(mpi_isotdpth(srcgrid%ncoords))

    !---------------------------------------------------------------------

    ! Loop through local variable

    do k = 1, nlev

       ! Initialize local variables

       srcgrid_var  = 0.0
       mpi_isotdpth = 0.0

       ! If on slave (compute) node (task), receive variables, compute
       ! variables, and send variables to master (root) node (task)
    
       if(mpi_procid .ne. mpi_masternode) then

          ! Loop through local variable and proceed accordingly
       
          do i = srcgrid%mpi_count_begin(mpi_procid),                       &
               & srcgrid%mpi_count_end(mpi_procid)

             ! Compute local variable

             call calculate_isothermdepth(ocean_hycom_zdim,                 &
                  & hycom%temp(i,1:ocean_hycom_zdim),                       &
                  & hycom%depth(i,1:ocean_hycom_zdim),tlevs(k),isotdpth)

             ! Define local variable

             if(isotdpth .le. 0.0 .or.                                      &
                  & maxval(hycom%temp(i,1:ocean_hycom_zdim)) .gt. 1.e20     &
                  & .or. maxval(hycom%temp(i,1:ocean_hycom_zdim)) .lt.      &
                  & tlevs(k)) isotdpth = spval
             mpi_isotdpth(i) = isotdpth

          end do ! do i = srcgrid%mpi_count_begin(mpi_procid),              &
                 ! srcgrid%mpi_count_end(mpi_procid)

       end if ! if(mpi_procid .ne. mpi_masternode)

       ! Define local variable

       call mpi_reduce(mpi_isotdpth(1:srcgrid%ncoords),                     &
            & srcgrid_var(1:srcgrid%ncoords),srcgrid%ncoords,mpi_real,      &
            & mpi_sum,mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
       
       call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,mpi_masternode,  &
            & mpi_comm_world,mpi_ierror)

       ! Compute local variable

       call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
       
       call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,mpi_masternode,  &
            & mpi_comm_world,mpi_ierror)    

       ! Define local variable
       
       grid%array((k-1)*dstgrid%ncoords+1:(k*dstgrid%ncoords)) =            &
            & dstgrid_var(1:dstgrid%ncoords)

    end do ! k = 1, nlev

    !---------------------------------------------------------------------

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(grid%array,(dstgrid%ncoords*nlev),mpi_real,              &
         & mpi_masternode,mpi_comm_world,mpi_ierror) 

    !---------------------------------------------------------------------

    ! Deallocate memory for local variables

    if(allocated(srcgrid_var))  deallocate(srcgrid_var)
    if(allocated(dstgrid_var))  deallocate(dstgrid_var)
    if(allocated(mpi_isotdpth)) deallocate(mpi_isotdpth)

    !=====================================================================

  end subroutine calculate_isotdpth

  !=======================================================================

  ! calculate_oceanuvelocity.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_oceanuvelocity(hycom,grid)

    ! Define variables passed to routine

    type(hycomvariables)                                           :: hycom
    type(variable_info)                                            :: grid

    ! Define variables computed within routine

    real(r_kind),              dimension(:,:),         allocatable :: uvel_grid
    real(r_kind),              dimension(:),           allocatable :: workgrid
    real(r_kind),              dimension(:),           allocatable :: uvel_isobathy
    real(r_kind),              dimension(:),           allocatable :: uvel_isotherm
    real(r_kind),              dimension(:),           allocatable :: mpi_uvel
    real(r_kind)                                                   :: uvel
    real(r_kind)                                                   :: vvel

    ! Define counting variables

    integer                                                        :: i, j, k

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(srcgrid_var)) deallocate(srcgrid_var)
    if(allocated(dstgrid_var)) deallocate(dstgrid_var)

    ! Allocate memory for local variable
    
    if(.not. allocated(uvel_grid))                                         &
         & allocate(uvel_grid(srcgrid%ncoords,ocean_hycom_zdim))
    if(.not. allocated(workgrid))                                          &
         & allocate(workgrid(dstgrid%ncoords*ocean_hycom_zdim))
    if(.not. allocated(srcgrid_var))                                       &
         & allocate(srcgrid_var(srcgrid%ncoords))
    if(.not. allocated(dstgrid_var))                                       &
         & allocate(dstgrid_var(dstgrid%ncoords))
    if(.not. allocated(uvel_isobathy))                                     &
         & allocate(uvel_isobathy(srcgrid%ncoords))
    if(.not. allocated(uvel_isotherm))                                     &
         & allocate(uvel_isotherm(srcgrid%ncoords))
    if(.not. allocated(mpi_uvel))                                          &
         & allocate(mpi_uvel(srcgrid%ncoords))

    !---------------------------------------------------------------------

    ! Initialize local variable

    workgrid = spval

    ! Loop through local variable

    do k = 1, ocean_hycom_zdim

       ! Initialize local variables
       
       srcgrid_var = 0.0
       mpi_uvel    = 0.0

       ! If on slave (compute) node (task), receive variables, compute
       ! variables, and send variables to master (root) node (task)
    
       if(mpi_procid .ne. mpi_masternode) then

          ! Loop through local variable and proceed accordingly
          
          do i = srcgrid%mpi_count_begin(mpi_procid),                      &
               & srcgrid%mpi_count_end(mpi_procid)
             
             ! Check local variable and proceed accordingly
             
             if(hycom%uvel(i,k) .lt. 1.e20 .and. hycom%vvel(i,k) .lt.      &
                  & 1.e20) then
                
                ! Define local variables
                
                uvel = hycom%uvel(i,k)
                vvel = hycom%vvel(i,k)
                
                ! Compute local variable
                
                mpi_uvel(i) = hycom%uvel(i,k)*cos(srcgrid%rotang(i)) +     &
                     & hycom%vvel(i,k)*sin(-srcgrid%rotang(i))
                
             else   ! if(hycom%uvel(i,k) .lt. 1.e20 .and. hycom%vvel(i,k)  &
                    ! .lt. 1.e20)
                
                ! Define local variable

                mpi_uvel(i) = spval

             end if ! if(hycom%uvel(i,k) .lt. 1.e20 .and. hycom%vvel(i,k)  &
                    ! .lt. 1.e20)

          end do ! do i = srcgrid%mpi_count_begin(mpi_procid),             &
                 ! srcgrid%mpi_count_end(mpi_procid)

       end if ! if(mpi_procid .ne. mpi_masternode)

       ! Define local variable

       call mpi_reduce(mpi_uvel(1:srcgrid%ncoords),                        &
            & uvel_grid(1:srcgrid%ncoords,k),srcgrid%ncoords,mpi_real,     &
            & mpi_sum,mpi_masternode,mpi_comm_world,mpi_ierror)

    end do ! do k = 1, ocean_hycom_zdim

    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    ! Check local variable and proceed accordingly

    if(.not. is_zinterp .and. .not. is_isointerp) then

       ! Loop through local variable

       do k = 1, nlev

          ! Define local variable

          srcgrid_var = uvel_grid(:,k)

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)
       
          call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Compute local variable

          call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Define local variable

          workgrid((k-1)*dstgrid%ncoords+1:(k*dstgrid%ncoords)) =          &
               & dstgrid_var(1:dstgrid%ncoords)

       end do !  do k = 1, nlev

    else if(is_zinterp) then   ! grid%array = workgrid

       ! Loop through local variable

       do k = 1, nlev

          ! Compute local variable

          call calculate_isobathyinterpolation(hycom,uvel_grid,zlevs(k),   &
               & uvel_isobathy)

          ! Define local variable

          srcgrid_var = uvel_isobathy

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)
       
          call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Define local variable

          dstgrid%weights = static_weights

          ! Loop through local variable

          do i = 1, dstgrid%neighbors

             ! Loop through local variable

             do j = 1, dstgrid%npasses

                where(dstgrid%depth .gt. zlevs(k))                         &
                     & dstgrid%weights(:,i,j) = 0.0

             end do !  do j = 1, dstgrid%npasses

          end do ! do i = 1, dstgrid%neighbors

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid%weights,(dstgrid%ncoords*                 &
               & dstgrid%neighbors*dstgrid%npasses),mpi_real,              &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Compute local variable

          call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Define local variable

          workgrid((k-1)*dstgrid%ncoords+1:(k*dstgrid%ncoords)) =          &
               & dstgrid_var(1:dstgrid%ncoords)
       
          ! Define local variable

          dstgrid%weights = static_weights

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid%weights,(dstgrid%ncoords*                 &
               & dstgrid%neighbors*dstgrid%npasses),mpi_real,              &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

       end do ! do k = 1, nlev

    else if(is_isointerp) then ! grid%array = workgrid

       ! Loop through local variable

       do k = 1, nlev

          ! Compute local variable

          call calculate_isotherminterpolation(hycom,uvel_grid,tlevs(k),   &
               & uvel_isotherm)

          ! Define local variable

          srcgrid_var = uvel_isotherm

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)
       
          call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Define local variable

          dstgrid%weights = static_weights

          ! Loop through local variable

          do i = 1, dstgrid%neighbors

             ! Loop through local variable

             do j = 1, dstgrid%npasses

                ! Check local variable and proceed accordingly

                where(dstgrid%temperature .lt. tlevs(k) .and.              &
                     & dstgrid%temperature .eq. spval)                     &
                     & dstgrid%weights(:,i,j) = 0.0

             end do !  do j = 1, dstgrid%npasses

          end do ! do i = 1, dstgrid%neighbors

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid%weights,(dstgrid%ncoords*                 &
               & dstgrid%neighbors*dstgrid%npasses),mpi_real,              &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Compute local variable

          call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Define local variable

          workgrid((k-1)*dstgrid%ncoords+1:(k*dstgrid%ncoords)) =          &
               & dstgrid_var(1:dstgrid%ncoords)
       
          ! Define local variable

          dstgrid%weights = static_weights

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid%weights,(dstgrid%ncoords*                 &
               & dstgrid%neighbors*dstgrid%npasses),mpi_real,              &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

       end do ! do k = 1, nlev

    end if                     ! grid%array = workgrid

    ! Define local variable

    grid%array = workgrid

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(grid%array,(nlev*dstgrid%ncoords),mpi_real,              &
         & mpi_masternode,mpi_comm_world,mpi_ierror) 

    !---------------------------------------------------------------------

    ! Deallocate memory for local variables

    if(allocated(uvel_grid))     deallocate(uvel_grid)
    if(allocated(workgrid))      deallocate(workgrid)
    if(allocated(srcgrid_var))   deallocate(srcgrid_var)
    if(allocated(dstgrid_var))   deallocate(dstgrid_var)
    if(allocated(uvel_isobathy)) deallocate(uvel_isobathy)
    if(allocated(mpi_uvel))      deallocate(mpi_uvel)

    !=====================================================================

  end subroutine calculate_oceanuvelocity

  !=======================================================================

  ! calculate_oceanvvelocity.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_oceanvvelocity(hycom,grid)

    ! Define variables passed to routine

    type(hycomvariables)                                           :: hycom
    type(variable_info)                                            :: grid

    ! Define variables computed within routine

    real(r_kind),              dimension(:,:),         allocatable :: vvel_grid
    real(r_kind),              dimension(:),           allocatable :: workgrid
    real(r_kind),              dimension(:),           allocatable :: vvel_isobathy
    real(r_kind),              dimension(:),           allocatable :: vvel_isotherm
    real(r_kind),              dimension(:),           allocatable :: mpi_vvel
    real(r_kind)                                                   :: uvel
    real(r_kind)                                                   :: vvel

    ! Define counting variables

    integer                                                        :: i, j, k

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(srcgrid_var)) deallocate(srcgrid_var)
    if(allocated(dstgrid_var)) deallocate(dstgrid_var)

    ! Allocate memory for local variable
    
    if(.not. allocated(vvel_grid))                                         &
         & allocate(vvel_grid(srcgrid%ncoords,ocean_hycom_zdim))
    if(.not. allocated(workgrid))                                          &
         & allocate(workgrid(dstgrid%ncoords*ocean_hycom_zdim))
    if(.not. allocated(srcgrid_var))                                       &
         & allocate(srcgrid_var(srcgrid%ncoords))
    if(.not. allocated(dstgrid_var))                                       &
         & allocate(dstgrid_var(dstgrid%ncoords))
    if(.not. allocated(vvel_isobathy))                                     &
         & allocate(vvel_isobathy(srcgrid%ncoords))
    if(.not. allocated(vvel_isotherm))                                     &
         & allocate(vvel_isotherm(srcgrid%ncoords))
    if(.not. allocated(mpi_vvel))                                          &
         & allocate(mpi_vvel(srcgrid%ncoords))

    !---------------------------------------------------------------------

    ! Initialize local variable

    workgrid = spval

    ! Loop through local variable

    do k = 1, ocean_hycom_zdim

       ! Initialize local variables
       
       srcgrid_var = 0.0
       mpi_vvel    = 0.0

       ! If on slave (compute) node (task), receive variables, compute
       ! variables, and send variables to master (root) node (task)
    
       if(mpi_procid .ne. mpi_masternode) then

          ! Loop through local variable and proceed accordingly
          
          do i = srcgrid%mpi_count_begin(mpi_procid),                      &
               & srcgrid%mpi_count_end(mpi_procid)
             
             ! Check local variable and proceed accordingly
             
             if(hycom%uvel(i,k) .lt. 1.e20 .and. hycom%vvel(i,k) .lt.      &
                  & 1.e20) then
                
                ! Define local variables
                
                uvel = hycom%uvel(i,k)
                vvel = hycom%vvel(i,k)
                
                ! Compute local variable
                
                mpi_vvel(i) = hycom%vvel(i,k)*cos(srcgrid%rotang(i)) -     &
                     & hycom%uvel(i,k)*sin(-srcgrid%rotang(i))
        
             else   ! if(hycom%uvel(i,k) .lt. 1.e20 .and. hycom%vvel(i,k)  &
                    ! .lt. 1.e20)
                
                ! Define local variable

                mpi_vvel(i) = spval

             end if ! if(hycom%uvel(i,k) .lt. 1.e20 .and. hycom%vvel(i,k)  &
                    ! .lt. 1.e20)

          end do ! do i = srcgrid%mpi_count_begin(mpi_procid),             &
                 ! srcgrid%mpi_count_end(mpi_procid)

       end if ! if(mpi_procid .ne. mpi_masternode)

       ! Define local variable

       call mpi_reduce(mpi_vvel(1:srcgrid%ncoords),                        &
            & vvel_grid(1:srcgrid%ncoords,k),srcgrid%ncoords,mpi_real,     &
            & mpi_sum,mpi_masternode,mpi_comm_world,mpi_ierror)

    end do ! do k = 1, ocean_hycom_zdim

    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    ! Check local variable and proceed accordingly

    if(.not. is_zinterp .and. .not. is_isointerp) then

       ! Loop through local variable

       do k = 1, nlev

          ! Define local variable

          srcgrid_var = vvel_grid(:,k)

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)
       
          call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Compute local variable

          call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Define local variable

          workgrid((k-1)*dstgrid%ncoords+1:(k*dstgrid%ncoords)) =          &
               & dstgrid_var(1:dstgrid%ncoords)

       end do !  do k = 1, nlev

    else if(is_zinterp) then   ! grid%array = workgrid

       ! Loop through local variable

       do k = 1, nlev

          ! Compute local variable

          call calculate_isobathyinterpolation(hycom,vvel_grid,zlevs(k),   &
               & vvel_isobathy)

          ! Define local variable

          srcgrid_var = vvel_isobathy

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)
       
          call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Define local variable

          dstgrid%weights = static_weights

          ! Loop through local variable

          do i = 1, dstgrid%neighbors

             ! Loop through local variable

             do j = 1, dstgrid%npasses

                where(dstgrid%depth .gt. zlevs(k))                         &
                     & dstgrid%weights(:,i,j) = 0.0

             end do !  do j = 1, dstgrid%npasses

          end do ! do i = 1, dstgrid%neighbors

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid%weights,(dstgrid%ncoords*                 &
               & dstgrid%neighbors*dstgrid%npasses),mpi_real,              &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Compute local variable

          call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Define local variable

          workgrid((k-1)*dstgrid%ncoords+1:(k*dstgrid%ncoords)) =          &
               & dstgrid_var(1:dstgrid%ncoords)
       
          ! Define local variable

          dstgrid%weights = static_weights

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid%weights,(dstgrid%ncoords*                 &
               & dstgrid%neighbors*dstgrid%npasses),mpi_real,              &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

       end do ! do k = 1, nlev

    else if(is_isointerp) then ! grid%array = workgrid

       ! Loop through local variable

       do k = 1, nlev

          ! Compute local variable

          call calculate_isotherminterpolation(hycom,vvel_grid,tlevs(k),   &
               & vvel_isotherm)

          ! Define local variable

          srcgrid_var = vvel_isotherm

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)
       
          call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Define local variable

          dstgrid%weights = static_weights

          ! Loop through local variable

          do i = 1, dstgrid%neighbors

             ! Loop through local variable

             do j = 1, dstgrid%npasses

                ! Check local variable and proceed accordingly

                where(dstgrid%temperature .lt. tlevs(k) .and.              &
                     & dstgrid%temperature .eq. spval)                     &
                     & dstgrid%weights(:,i,j) = 0.0

             end do !  do j = 1, dstgrid%npasses

          end do ! do i = 1, dstgrid%neighbors

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid%weights,(dstgrid%ncoords*                 &
               & dstgrid%neighbors*dstgrid%npasses),mpi_real,              &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Compute local variable

          call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Define local variable

          workgrid((k-1)*dstgrid%ncoords+1:(k*dstgrid%ncoords)) =          &
               & dstgrid_var(1:dstgrid%ncoords)
       
          ! Define local variable

          dstgrid%weights = static_weights

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid%weights,(dstgrid%ncoords*                 &
               & dstgrid%neighbors*dstgrid%npasses),mpi_real,              &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

       end do ! do k = 1, nlev

    end if                     ! grid%array = workgrid

    ! Define local variable

    grid%array = workgrid

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(grid%array,(nlev*dstgrid%ncoords),mpi_real,              &
         & mpi_masternode,mpi_comm_world,mpi_ierror) 

    !---------------------------------------------------------------------

    ! Deallocate memory for local variables

    if(allocated(vvel_grid))     deallocate(vvel_grid)
    if(allocated(workgrid))      deallocate(workgrid)
    if(allocated(srcgrid_var))   deallocate(srcgrid_var)
    if(allocated(dstgrid_var))   deallocate(dstgrid_var)
    if(allocated(vvel_isobathy)) deallocate(vvel_isobathy)
    if(allocated(mpi_vvel))      deallocate(mpi_vvel)

    !=====================================================================

  end subroutine calculate_oceanvvelocity

  !=======================================================================

  ! calculate_oceantemperature.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_oceantemperature(hycom,grid)

    ! Define variables passed to routine

    type(hycomvariables)                                           :: hycom
    type(variable_info)                                            :: grid

    ! Define variables computed within routine

    real(r_kind),              dimension(:,:),         allocatable :: temp_grid
    real(r_kind),              dimension(:),           allocatable :: workgrid
    real(r_kind),              dimension(:),           allocatable :: temp_isobathy
    real(r_kind),              dimension(:),           allocatable :: temp_isotherm
    real(r_kind),              dimension(:),           allocatable :: mpi_temp

    ! Define counting variables

    integer                                                        :: i, j, k

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(srcgrid_var)) deallocate(srcgrid_var)
    if(allocated(dstgrid_var)) deallocate(dstgrid_var)

    ! Allocate memory for local variable
    
    if(.not. allocated(temp_grid))                                         &
         & allocate(temp_grid(srcgrid%ncoords,ocean_hycom_zdim))
    if(.not. allocated(workgrid))                                          &
         & allocate(workgrid(dstgrid%ncoords*ocean_hycom_zdim))
    if(.not. allocated(srcgrid_var))                                       &
         & allocate(srcgrid_var(srcgrid%ncoords))
    if(.not. allocated(dstgrid_var))                                       &
         & allocate(dstgrid_var(dstgrid%ncoords))
    if(.not. allocated(temp_isobathy))                                     &
         & allocate(temp_isobathy(srcgrid%ncoords))
    if(.not. allocated(temp_isotherm))                                     &
         & allocate(temp_isotherm(srcgrid%ncoords))
    if(.not. allocated(mpi_temp))                                          &
         & allocate(mpi_temp(srcgrid%ncoords))

    !---------------------------------------------------------------------

    ! Initialize local variable

    workgrid = spval

    ! Loop through local variable

    do k = 1, ocean_hycom_zdim

       ! Initialize local variables
       
       srcgrid_var = 0.0
       mpi_temp    = 0.0

       ! If on slave (compute) node (task), receive variables, compute
       ! variables, and send variables to master (root) node (task)
    
       if(mpi_procid .ne. mpi_masternode) then

          ! Loop through local variable and proceed accordingly
          
          do i = srcgrid%mpi_count_begin(mpi_procid),                      &
               & srcgrid%mpi_count_end(mpi_procid)
             
             ! Check local variable and proceed accordingly
             
             if(hycom%temp(i,k) .lt. 1.e20) then
                
                ! Compute local variable
                
                mpi_temp(i) = hycom%temp(i,k) + t0c
                
             else   ! if(hycom%temp(i,k) .lt. 1.e20)
                
                ! Define local variable

                mpi_temp(i) = spval

             end if ! if(hycom%temp(i,k) .lt. 1.e20)

          end do ! do i = srcgrid%mpi_count_begin(mpi_procid),             &
                 ! srcgrid%mpi_count_end(mpi_procid)

       end if ! if(mpi_procid .ne. mpi_masternode)

       ! Define local variable

       call mpi_reduce(mpi_temp(1:srcgrid%ncoords),                        &
            & temp_grid(1:srcgrid%ncoords,k),srcgrid%ncoords,mpi_real,     &
            & mpi_sum,mpi_masternode,mpi_comm_world,mpi_ierror)

    end do ! do k = 1, ocean_hycom_zdim

    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    ! Check local variable and proceed accordingly

    if(.not. is_zinterp .and. .not. is_isointerp) then

       ! Loop through local variable

       do k = 1, nlev

          ! Define local variable

          srcgrid_var = temp_grid(:,k)

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)
       
          call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Compute local variable

          call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Define local variable

          workgrid((k-1)*dstgrid%ncoords+1:(k*dstgrid%ncoords)) =          &
               & dstgrid_var(1:dstgrid%ncoords)

       end do !  do k = 1, nlev

    else if(is_zinterp) then   ! grid%array = workgrid

       ! Loop through local variable

       do k = 1, nlev

          ! Compute local variable

          call calculate_isobathyinterpolation(hycom,temp_grid,zlevs(k),   &
               & temp_isobathy)

          ! Define local variable

          srcgrid_var = temp_isobathy

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)
       
          call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Define local variable

          dstgrid%weights = static_weights

          ! Loop through local variable

          do i = 1, dstgrid%neighbors

             ! Loop through local variable

             do j = 1, dstgrid%npasses

                where(dstgrid%depth .gt. zlevs(k))                         &
                     & dstgrid%weights(:,i,j) = 0.0

             end do !  do j = 1, dstgrid%npasses

          end do ! do i = 1, dstgrid%neighbors

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid%weights,(dstgrid%ncoords*                 &
               & dstgrid%neighbors*dstgrid%npasses),mpi_real,              &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Compute local variable

          call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Define local variable

          workgrid((k-1)*dstgrid%ncoords+1:(k*dstgrid%ncoords)) =          &
               & dstgrid_var(1:dstgrid%ncoords)
       
          ! Define local variable

          dstgrid%weights = static_weights

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid%weights,(dstgrid%ncoords*                 &
               & dstgrid%neighbors*dstgrid%npasses),mpi_real,              &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

       end do ! do k = 1, nlev

    else if(is_isointerp) then ! grid%array = workgrid

       ! Loop through local variable

       do k = 1, nlev

          ! Check local variable and proceed accordingly

          where(temp_grid .lt. 1.e20) temp_grid = temp_grid - t0c

          ! Compute local variable

          call calculate_isotherminterpolation(hycom,temp_grid,tlevs(k),   &
               & temp_isotherm)

          ! Define local variable

          srcgrid_var = temp_isotherm
          
          ! Check local variable and proceed accordingly

          where(srcgrid_var .lt. 1.e20) srcgrid_var = srcgrid_var + t0c

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)
       
          call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Define local variable

          dstgrid%weights = static_weights

          ! Loop through local variable

          do i = 1, dstgrid%neighbors

             ! Loop through local variable

             do j = 1, dstgrid%npasses

                ! Check local variable and proceed accordingly

                where(dstgrid%temperature .lt. tlevs(k) .and.              &
                     & dstgrid%temperature .eq. spval)                     &
                     & dstgrid%weights(:,i,j) = 0.0

             end do !  do j = 1, dstgrid%npasses

          end do ! do i = 1, dstgrid%neighbors

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid%weights,(dstgrid%ncoords*                 &
               & dstgrid%neighbors*dstgrid%npasses),mpi_real,              &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Compute local variable

          call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Define local variable

          workgrid((k-1)*dstgrid%ncoords+1:(k*dstgrid%ncoords)) =          &
               & dstgrid_var(1:dstgrid%ncoords)
       
          ! Define local variable

          dstgrid%weights = static_weights

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid%weights,(dstgrid%ncoords*                 &
               & dstgrid%neighbors*dstgrid%npasses),mpi_real,              &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

       end do ! do k = 1, nlev

    end if                     ! grid%array = workgrid

    ! Define local variable

    grid%array = workgrid

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(grid%array,(nlev*dstgrid%ncoords),mpi_real,              &
         & mpi_masternode,mpi_comm_world,mpi_ierror) 

    !---------------------------------------------------------------------

    ! Deallocate memory for local variables

    if(allocated(temp_grid))     deallocate(temp_grid)
    if(allocated(workgrid))      deallocate(workgrid)
    if(allocated(srcgrid_var))   deallocate(srcgrid_var)
    if(allocated(dstgrid_var))   deallocate(dstgrid_var)
    if(allocated(temp_isobathy)) deallocate(temp_isobathy)
    if(allocated(mpi_temp))      deallocate(mpi_temp)

    !=====================================================================

  end subroutine calculate_oceantemperature

  !=======================================================================

  ! calculate_oceansalinity.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_oceansalinity(hycom,grid)

    ! Define variables passed to routine

    type(hycomvariables)                                           :: hycom
    type(variable_info)                                            :: grid

    ! Define variables computed within routine

    real(r_kind),              dimension(:,:),         allocatable :: salin_grid
    real(r_kind),              dimension(:),           allocatable :: workgrid
    real(r_kind),              dimension(:),           allocatable :: salin_isobathy
    real(r_kind),              dimension(:),           allocatable :: salin_isotherm
    real(r_kind),              dimension(:),           allocatable :: mpi_salin

    ! Define counting variables

    integer                                                        :: i, j, k

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(srcgrid_var)) deallocate(srcgrid_var)
    if(allocated(dstgrid_var)) deallocate(dstgrid_var)

    ! Allocate memory for local variable
    
    if(.not. allocated(salin_grid))                                        &
         & allocate(salin_grid(srcgrid%ncoords,ocean_hycom_zdim))
    if(.not. allocated(workgrid))                                          &
         & allocate(workgrid(dstgrid%ncoords*ocean_hycom_zdim))
    if(.not. allocated(srcgrid_var))                                       &
         & allocate(srcgrid_var(srcgrid%ncoords))
    if(.not. allocated(dstgrid_var))                                       &
         & allocate(dstgrid_var(dstgrid%ncoords))
    if(.not. allocated(salin_isobathy))                                    &
         & allocate(salin_isobathy(srcgrid%ncoords))
    if(.not. allocated(salin_isotherm))                                    &
         & allocate(salin_isotherm(srcgrid%ncoords))
    if(.not. allocated(mpi_salin))                                         &
         & allocate(mpi_salin(srcgrid%ncoords))

    !---------------------------------------------------------------------

    ! Initialize local variable

    workgrid = spval

    ! Loop through local variable

    do k = 1, ocean_hycom_zdim

       ! Initialize local variables
       
       srcgrid_var = 0.0
       mpi_salin   = 0.0

       ! If on slave (compute) node (task), receive variables, compute
       ! variables, and send variables to master (root) node (task)
    
       if(mpi_procid .ne. mpi_masternode) then

          ! Loop through local variable and proceed accordingly
          
          do i = srcgrid%mpi_count_begin(mpi_procid),                      &
               & srcgrid%mpi_count_end(mpi_procid)
             
             ! Check local variable and proceed accordingly
             
             if(hycom%salin(i,k) .lt. 1.e20) then
                
                ! Compute local variable
                
                mpi_salin(i) = hycom%salin(i,k)
                
             else   ! if(hycom%salin(i,k) .lt. 1.e20)
                
                ! Define local variable

                mpi_salin(i) = spval

             end if ! if(hycom%salin(i,k) .lt. 1.e20)

          end do ! do i = srcgrid%mpi_count_begin(mpi_procid),             &
                 ! srcgrid%mpi_count_end(mpi_procid)

       end if ! if(mpi_procid .ne. mpi_masternode)

       ! Define local variable

       call mpi_reduce(mpi_salin(1:srcgrid%ncoords),                       &
            & salin_grid(1:srcgrid%ncoords,k),srcgrid%ncoords,mpi_real,    &
            & mpi_sum,mpi_masternode,mpi_comm_world,mpi_ierror)

    end do ! do k = 1, ocean_hycom_zdim

    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    ! Check local variable and proceed accordingly

    if(.not. is_zinterp .and. .not. is_isointerp) then

       ! Loop through local variable

       do k = 1, nlev

          ! Define local variable

          srcgrid_var = salin_grid(:,k)

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)
       
          call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Compute local variable

          call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Define local variable

          workgrid((k-1)*dstgrid%ncoords+1:(k*dstgrid%ncoords)) =          &
               & dstgrid_var(1:dstgrid%ncoords)

       end do !  do k = 1, nlev

    else if(is_zinterp) then   ! grid%array = workgrid

       ! Loop through local variable

       do k = 1, nlev

          ! Compute local variable

          call calculate_isobathyinterpolation(hycom,salin_grid,zlevs(k),  &
               & salin_isobathy)

          ! Define local variable

          srcgrid_var = salin_isobathy

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)
       
          call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Define local variable

          dstgrid%weights = static_weights

          ! Loop through local variable

          do i = 1, dstgrid%neighbors

             ! Loop through local variable

             do j = 1, dstgrid%npasses

                where(dstgrid%depth .gt. zlevs(k))                         &
                     & dstgrid%weights(:,i,j) = 0.0

             end do !  do j = 1, dstgrid%npasses

          end do ! do i = 1, dstgrid%neighbors

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid%weights,(dstgrid%ncoords*                 &
               & dstgrid%neighbors*dstgrid%npasses),mpi_real,              &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Compute local variable

          call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Define local variable

          workgrid((k-1)*dstgrid%ncoords+1:(k*dstgrid%ncoords)) =          &
               & dstgrid_var(1:dstgrid%ncoords)
       
          ! Define local variable

          dstgrid%weights = static_weights

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid%weights,(dstgrid%ncoords*                 &
               & dstgrid%neighbors*dstgrid%npasses),mpi_real,              &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

       end do ! do k = 1, nlev

    else if(is_isointerp) then ! grid%array = workgrid

       ! Loop through local variable

       do k = 1, nlev

          ! Compute local variable

          call calculate_isotherminterpolation(hycom,salin_grid,tlevs(k),  &
               & salin_isotherm)

          ! Define local variable

          srcgrid_var = salin_isotherm

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)
       
          call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Define local variable

          dstgrid%weights = static_weights

          ! Loop through local variable

          do i = 1, dstgrid%neighbors

             ! Loop through local variable

             do j = 1, dstgrid%npasses

                ! Check local variable and proceed accordingly

                where(dstgrid%temperature .lt. tlevs(k) .and.              &
                     & dstgrid%temperature .eq. spval)                     &
                     & dstgrid%weights(:,i,j) = 0.0

             end do !  do j = 1, dstgrid%npasses

          end do ! do i = 1, dstgrid%neighbors

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid%weights,(dstgrid%ncoords*                 &
               & dstgrid%neighbors*dstgrid%npasses),mpi_real,              &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Compute local variable

          call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Define local variable

          workgrid((k-1)*dstgrid%ncoords+1:(k*dstgrid%ncoords)) =          &
               & dstgrid_var(1:dstgrid%ncoords)
       
          ! Define local variable

          dstgrid%weights = static_weights

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid%weights,(dstgrid%ncoords*                 &
               & dstgrid%neighbors*dstgrid%npasses),mpi_real,              &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

       end do ! do k = 1, nlev

    end if                     ! grid%array = workgrid

    ! Define local variable

    grid%array = workgrid

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(grid%array,(nlev*dstgrid%ncoords),mpi_real,              &
         & mpi_masternode,mpi_comm_world,mpi_ierror) 

    !---------------------------------------------------------------------

    ! Deallocate memory for local variables

    if(allocated(salin_grid))     deallocate(salin_grid)
    if(allocated(workgrid))       deallocate(workgrid)
    if(allocated(srcgrid_var))    deallocate(srcgrid_var)
    if(allocated(dstgrid_var))    deallocate(dstgrid_var)
    if(allocated(salin_isobathy)) deallocate(salin_isobathy)
    if(allocated(mpi_salin))      deallocate(mpi_salin)

    !=====================================================================

  end subroutine calculate_oceansalinity

  !=======================================================================

  ! calculate_oceandensity.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_oceandensity(hycom,grid)

    ! Define variables passed to routine

    type(hycomvariables)                                           :: hycom
    type(variable_info)                                            :: grid

    ! Define variables computed within routine

    real(r_kind),              dimension(:,:),         allocatable :: dens_grid
    real(r_kind),              dimension(:),           allocatable :: workgrid
    real(r_kind),              dimension(:),           allocatable :: dens_isobathy
    real(r_kind),              dimension(:),           allocatable :: dens_isotherm
    real(r_kind),              dimension(:),           allocatable :: mpi_dens

    ! Define counting variables

    integer                                                        :: i, j, k

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(srcgrid_var)) deallocate(srcgrid_var)
    if(allocated(dstgrid_var)) deallocate(dstgrid_var)

    ! Allocate memory for local variable
    
    if(.not. allocated(dens_grid))                                         &
         & allocate(dens_grid(srcgrid%ncoords,ocean_hycom_zdim))
    if(.not. allocated(workgrid))                                          &
         & allocate(workgrid(dstgrid%ncoords*ocean_hycom_zdim))
    if(.not. allocated(srcgrid_var))                                       &
         & allocate(srcgrid_var(srcgrid%ncoords))
    if(.not. allocated(dstgrid_var))                                       &
         & allocate(dstgrid_var(dstgrid%ncoords))
    if(.not. allocated(dens_isobathy))                                     &
         & allocate(dens_isobathy(srcgrid%ncoords))
    if(.not. allocated(dens_isotherm))                                     &
         & allocate(dens_isotherm(srcgrid%ncoords))
    if(.not. allocated(mpi_dens))                                          &
         & allocate(mpi_dens(srcgrid%ncoords))

    !---------------------------------------------------------------------

    ! Initialize local variable

    workgrid = spval

    ! Loop through local variable

    do k = 1, ocean_hycom_zdim

       ! Initialize local variables
       
       srcgrid_var = 0.0
       mpi_dens    = 0.0

       ! If on slave (compute) node (task), receive variables, compute
       ! variables, and send variables to master (root) node (task)
    
       if(mpi_procid .ne. mpi_masternode) then

          ! Loop through local variable and proceed accordingly
          
          do i = srcgrid%mpi_count_begin(mpi_procid),                      &
               & srcgrid%mpi_count_end(mpi_procid)
             
             ! Check local variable and proceed accordingly
             
             if(hycom%temp(i,k) .lt. 1.e20 .and. hycom%salin(i,k) .lt.     &
                  & 1.e20) then
                
                ! Compute local variable
                
                call calculate_seawaterdensity(hycom%temp(i,k)+t0c,        &
                     & hycom%salin(i,k),mpi_dens(i))
                
             else   ! if(hycom%temp(i,k) .lt. 1.e20 .and.                  &
                    ! hycom%salin(i,k) .lt. 1.e20)
                
                ! Define local variable

                mpi_dens(i) = spval

             end if ! if(hycom%temp(i,k) .lt. 1.e20 .and.                  &
                    ! hycom%salin(i,k) .lt. 1.e20) 

          end do ! do i = srcgrid%mpi_count_begin(mpi_procid),             &
                 ! srcgrid%mpi_count_end(mpi_procid)

       end if ! if(mpi_procid .ne. mpi_masternode)

       ! Define local variable

       call mpi_reduce(mpi_dens(1:srcgrid%ncoords),                        &
            & dens_grid(1:srcgrid%ncoords,k),srcgrid%ncoords,mpi_real,     &
            & mpi_sum,mpi_masternode,mpi_comm_world,mpi_ierror)

    end do ! do k = 1, ocean_hycom_zdim

    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    ! Check local variable and proceed accordingly

    if(.not. is_zinterp .and. .not. is_isointerp) then

       ! Loop through local variable

       do k = 1, nlev

          ! Define local variable

          srcgrid_var = dens_grid(:,k)

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)
       
          call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Compute local variable

          call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Define local variable

          workgrid((k-1)*dstgrid%ncoords+1:(k*dstgrid%ncoords)) =          &
               & dstgrid_var(1:dstgrid%ncoords)

       end do !  do k = 1, nlev

    else if(is_zinterp) then   ! grid%array = workgrid

       ! Loop through local variable

       do k = 1, nlev

          ! Compute local variable

          call calculate_isobathyinterpolation(hycom,dens_grid,zlevs(k),   &
               & dens_isobathy)

          ! Define local variable

          srcgrid_var = dens_isobathy

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)
       
          call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Define local variable

          dstgrid%weights = static_weights

          ! Loop through local variable

          do i = 1, dstgrid%neighbors

             ! Loop through local variable

             do j = 1, dstgrid%npasses

                where(dstgrid%depth .gt. zlevs(k))                         &
                     & dstgrid%weights(:,i,j) = 0.0

             end do !  do j = 1, dstgrid%npasses

          end do ! do i = 1, dstgrid%neighbors

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid%weights,(dstgrid%ncoords*                 &
               & dstgrid%neighbors*dstgrid%npasses),mpi_real,              &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Compute local variable

          call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Define local variable

          workgrid((k-1)*dstgrid%ncoords+1:(k*dstgrid%ncoords)) =          &
               & dstgrid_var(1:dstgrid%ncoords)
       
          ! Define local variable

          dstgrid%weights = static_weights

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid%weights,(dstgrid%ncoords*                 &
               & dstgrid%neighbors*dstgrid%npasses),mpi_real,              &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

       end do ! do k = 1, nlev

    else if(is_isointerp) then ! grid%array = workgrid

       ! Loop through local variable

       do k = 1, nlev

          ! Compute local variable

          call calculate_isotherminterpolation(hycom,dens_grid,tlevs(k),   &
               & dens_isotherm)

          ! Define local variable

          srcgrid_var = dens_isotherm

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)
       
          call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Define local variable

          dstgrid%weights = static_weights

          ! Loop through local variable

          do i = 1, dstgrid%neighbors

             ! Loop through local variable

             do j = 1, dstgrid%npasses

                ! Check local variable and proceed accordingly

                where(dstgrid%temperature .lt. tlevs(k) .and.              &
                     & dstgrid%temperature .eq. spval)                     &
                     & dstgrid%weights(:,i,j) = 0.0

             end do !  do j = 1, dstgrid%npasses

          end do ! do i = 1, dstgrid%neighbors

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid%weights,(dstgrid%ncoords*                 &
               & dstgrid%neighbors*dstgrid%npasses),mpi_real,              &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Compute local variable

          call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Define local variable

          workgrid((k-1)*dstgrid%ncoords+1:(k*dstgrid%ncoords)) =          &
               & dstgrid_var(1:dstgrid%ncoords)
       
          ! Define local variable

          dstgrid%weights = static_weights

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid%weights,(dstgrid%ncoords*                 &
               & dstgrid%neighbors*dstgrid%npasses),mpi_real,              &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

       end do ! do k = 1, nlev

    end if                     ! grid%array = workgrid

    ! Define local variable

    grid%array = workgrid

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(grid%array,(nlev*dstgrid%ncoords),mpi_real,              &
         & mpi_masternode,mpi_comm_world,mpi_ierror) 

    !---------------------------------------------------------------------

    ! Deallocate memory for local variables

    if(allocated(dens_grid))     deallocate(dens_grid)
    if(allocated(workgrid))      deallocate(workgrid)
    if(allocated(srcgrid_var))   deallocate(srcgrid_var)
    if(allocated(dstgrid_var))   deallocate(dstgrid_var)
    if(allocated(dens_isobathy)) deallocate(dens_isobathy)
    if(allocated(mpi_dens))      deallocate(mpi_dens)

    !=====================================================================

  end subroutine calculate_oceandensity

  !=======================================================================

  ! calculate_oceandepth.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_oceandepth(hycom,grid)

    ! Define variables passed to routine

    type(hycomvariables)                                           :: hycom
    type(variable_info)                                            :: grid

    ! Define variables computed within routine

    real(r_kind),              dimension(:,:),         allocatable :: dpth_grid
    real(r_kind),              dimension(:),           allocatable :: workgrid
    real(r_kind),              dimension(:),           allocatable :: dpth_isobathy
    real(r_kind),              dimension(:),           allocatable :: dpth_isotherm
    real(r_kind),              dimension(:),           allocatable :: mpi_dpth

    ! Define counting variables

    integer                                                        :: i, j, k

    !=====================================================================

    ! Deallocate memory for local variables

    if(allocated(srcgrid_var)) deallocate(srcgrid_var)
    if(allocated(dstgrid_var)) deallocate(dstgrid_var)

    ! Allocate memory for local variable
    
    if(.not. allocated(dpth_grid))                                         &
         & allocate(dpth_grid(srcgrid%ncoords,ocean_hycom_zdim))
    if(.not. allocated(workgrid))                                          &
         & allocate(workgrid(dstgrid%ncoords*ocean_hycom_zdim))
    if(.not. allocated(srcgrid_var))                                       &
         & allocate(srcgrid_var(srcgrid%ncoords))
    if(.not. allocated(dstgrid_var))                                       &
         & allocate(dstgrid_var(dstgrid%ncoords))
    if(.not. allocated(dpth_isobathy))                                     &
         & allocate(dpth_isobathy(srcgrid%ncoords))
    if(.not. allocated(dpth_isotherm))                                     &
         & allocate(dpth_isotherm(srcgrid%ncoords))
    if(.not. allocated(mpi_dpth))                                          &
         & allocate(mpi_dpth(srcgrid%ncoords))

    !---------------------------------------------------------------------

    ! Initialize local variable

    workgrid = spval

    ! Loop through local variable

    do k = 1, ocean_hycom_zdim

       ! Initialize local variables
       
       srcgrid_var = 0.0
       mpi_dpth    = 0.0

       ! If on slave (compute) node (task), receive variables, compute
       ! variables, and send variables to master (root) node (task)
    
       if(mpi_procid .ne. mpi_masternode) then

          ! Loop through local variable and proceed accordingly
          
          do i = srcgrid%mpi_count_begin(mpi_procid),                      &
               & srcgrid%mpi_count_end(mpi_procid)
             
             ! Check local variable and proceed accordingly
             
             if(hycom%thknss(i,1) .lt. 1.e20) then
                
                ! Define local variable

                mpi_dpth(i) = hycom%depth(i,k)
                
             else   ! if(hycom%thknss(i,1) .lt. 1.e20)
                
                ! Define local variable

                mpi_dpth(i) = spval

             end if ! if(hycom%thknss(i,1) .lt. 1.e20)

          end do ! do i = srcgrid%mpi_count_begin(mpi_procid),             &
                 ! srcgrid%mpi_count_end(mpi_procid)

       end if ! if(mpi_procid .ne. mpi_masternode)

       ! Define local variable

       call mpi_reduce(mpi_dpth(1:srcgrid%ncoords),                        &
            & dpth_grid(1:srcgrid%ncoords,k),srcgrid%ncoords,mpi_real,     &
            & mpi_sum,mpi_masternode,mpi_comm_world,mpi_ierror)

    end do ! do k = 1, ocean_hycom_zdim

    ! Enable the root task to catch up from I/O and calculations

    call mpi_barrier(mpi_comm_world,mpi_ierror)

    ! Check local variable and proceed accordingly

    if(.not. is_zinterp .and. .not. is_isointerp) then

       ! Loop through local variable

       do k = 1, nlev

          ! Define local variable

          srcgrid_var = dpth_grid(:,k)

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)
       
          call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Compute local variable

          call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Define local variable

          workgrid((k-1)*dstgrid%ncoords+1:(k*dstgrid%ncoords)) =          &
               & dstgrid_var(1:dstgrid%ncoords)

       end do !  do k = 1, nlev

    else if(is_zinterp) then   ! grid%array = workgrid

       ! Loop through local variable

       do k = 1, nlev

          ! Compute local variable

          call calculate_isobathyinterpolation(hycom,dpth_grid,zlevs(k),   &
               & dpth_isobathy)

          ! Define local variable

          srcgrid_var = dpth_isobathy

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)
       
          call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Define local variable

          dstgrid%weights = static_weights

          ! Loop through local variable

          do i = 1, dstgrid%neighbors

             ! Loop through local variable

             do j = 1, dstgrid%npasses

                where(dstgrid%depth .gt. zlevs(k))                         &
                     & dstgrid%weights(:,i,j) = 0.0

             end do !  do j = 1, dstgrid%npasses

          end do ! do i = 1, dstgrid%neighbors

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid%weights,(dstgrid%ncoords*                 &
               & dstgrid%neighbors*dstgrid%npasses),mpi_real,              &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Compute local variable

          call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Define local variable

          workgrid((k-1)*dstgrid%ncoords+1:(k*dstgrid%ncoords)) =          &
               & dstgrid_var(1:dstgrid%ncoords)
       
          ! Define local variable

          dstgrid%weights = static_weights

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid%weights,(dstgrid%ncoords*                 &
               & dstgrid%neighbors*dstgrid%npasses),mpi_real,              &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

       end do ! do k = 1, nlev

    else if(is_isointerp) then ! grid%array = workgrid

       ! Loop through local variable

       do k = 1, nlev

          ! Compute local variable

          call calculate_isotherminterpolation(hycom,dpth_grid,tlevs(k),   &
               & dpth_isotherm)

          ! Define local variable

          srcgrid_var = dpth_isotherm

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)
       
          call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Define local variable

          dstgrid%weights = static_weights

          ! Loop through local variable

          do i = 1, dstgrid%neighbors

             ! Loop through local variable

             do j = 1, dstgrid%npasses

                ! Check local variable and proceed accordingly

                where(dstgrid%temperature .lt. tlevs(k) .and.              &
                     & dstgrid%temperature .eq. spval)                     &
                     & dstgrid%weights(:,i,j) = 0.0

             end do !  do j = 1, dstgrid%npasses

          end do ! do i = 1, dstgrid%neighbors

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid%weights,(dstgrid%ncoords*                 &
               & dstgrid%neighbors*dstgrid%npasses),mpi_real,              &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Compute local variable

          call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,             &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

          ! Define local variable

          workgrid((k-1)*dstgrid%ncoords+1:(k*dstgrid%ncoords)) =          &
               & dstgrid_var(1:dstgrid%ncoords)
       
          ! Define local variable

          dstgrid%weights = static_weights

          ! Broadcast all necessary variables to slave (compute) nodes
          ! (tasks)

          call mpi_bcast(dstgrid%weights,(dstgrid%ncoords*                 &
               & dstgrid%neighbors*dstgrid%npasses),mpi_real,              &
               & mpi_masternode,mpi_comm_world,mpi_ierror)

       end do ! do k = 1, nlev

    end if                     ! grid%array = workgrid

    ! Define local variable

    grid%array = workgrid

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(grid%array,(nlev*dstgrid%ncoords),mpi_real,              &
         & mpi_masternode,mpi_comm_world,mpi_ierror) 

    !---------------------------------------------------------------------

    ! Deallocate memory for local variables

    if(allocated(dpth_grid))     deallocate(dpth_grid)
    if(allocated(workgrid))      deallocate(workgrid)
    if(allocated(srcgrid_var))   deallocate(srcgrid_var)
    if(allocated(dstgrid_var))   deallocate(dstgrid_var)
    if(allocated(dpth_isobathy)) deallocate(dpth_isobathy)
    if(allocated(mpi_dpth))      deallocate(mpi_dpth)

    !=====================================================================

  end subroutine calculate_oceandepth

  !=======================================================================

  ! calculate_oceanmixedlayerkineticenergy.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_oceanmixedlayerkineticenergy(hycom,grid)

    ! Define variables passed to routine

    type(hycomvariables)                                           :: hycom
    type(variable_info)                                            :: grid

    ! Define variables computed within routine

    real(r_kind),              dimension(:),           allocatable :: mpi_omlke
    real(r_kind)                                                   :: uvel
    real(r_kind)                                                   :: vvel
    real(r_kind)                                                   :: z
    real(r_kind)                                                   :: omlu
    real(r_kind)                                                   :: omlv
    real(r_kind)                                                   :: omld
    real(r_kind)                                                   :: dz

    ! Define counting variables

    integer                                                        :: i, j, k
    integer                                                        :: count

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(srcgrid_var))                                       &
         & allocate(srcgrid_var(srcgrid%ncoords))
    if(.not. allocated(dstgrid_var))                                       &
         & allocate(dstgrid_var(dstgrid%ncoords))
    if(.not. allocated(mpi_omlke))                                         &
         & allocate(mpi_omlke(srcgrid%ncoords))

    !---------------------------------------------------------------------

    ! Initialize local variables

    srcgrid_var = 0.0
    mpi_omlke   = 0.0

    ! If on slave (compute) node (task), receive variables, compute
    ! variables, and send variables to master (root) node (task)
    
    if(mpi_procid .ne. mpi_masternode) then

       ! Loop through local variable and proceed accordingly
       
       do i = srcgrid%mpi_count_begin(mpi_procid),                         &
            & srcgrid%mpi_count_end(mpi_procid)

          ! Initialize local variable

          count = 0

          ! Check local variable and proceed accordingly

          if(hycom%mix_dpth(i) .lt. 1.e20) then

             ! Define local variables

             omld = hycom%mix_dpth(i)/9860.0
             omlu = 0.0
             omlv = 0.0
             dz   = depth_integral_dz

             ! Initialize local variable

             z = omld

             ! Loop through local variable

             do while(z .ge. 0.0)

                ! Compute local variables

                call calculate_variableatdepth(ocean_hycom_zdim,            &
                     & hycom%depth(i,1:ocean_hycom_zdim),                   &
                     & hycom%uvel(i,1:ocean_hycom_zdim),z,uvel)
                omlu = omlu + uvel
                call calculate_variableatdepth(ocean_hycom_zdim,            &
                     & hycom%depth(i,1:ocean_hycom_zdim),                   &
                     & hycom%vvel(i,1:ocean_hycom_zdim),z,vvel)
                omlv = omlv + vvel
                
                ! Update local variables

                z     = z - dz
                count = count + 1
                
             end do ! do while(z .ge. 0.0)

          end if ! if(hycom%mix_dpth(i) .lt. 1.e20)

          ! Define local variable

          if(count .ge. 1) omlu = omlu/real(count)
          if(count .ge. 1) omlv = omlv/real(count)
          if(count .eq. 0) omlu = spval
          if(count .eq. 0) omlv = spval
          
          ! Check local variable and proceed accordingly

          if(omlu .ne. spval .and. omlv .ne. spval) then

             ! Compute local variable

             mpi_omlke(i) = 0.5*(((omlu*cos(srcgrid%rotang(i)) +            &
                  & omlv*sin(-srcgrid%rotang(i)))*                          &
                  & (omlu*cos(srcgrid%rotang(i)) +                          &
                  & omlv*sin(-srcgrid%rotang(i)))) +                        &
                  & ((omlv*cos(srcgrid%rotang(i)) -                         &
                  & omlu*sin(-srcgrid%rotang(i)))*                          &
                  & omlv*cos(srcgrid%rotang(i)) -                           &
                  & omlu*sin(-srcgrid%rotang(i))))

          else   ! if(omlu .ne. spval .and. omlv .ne. spval)

             ! Define local variable

             mpi_omlke(i) = spval

          end if ! if(omlu .ne. spval .and. omlv .ne. spval)

       end do ! do i = srcgrid%mpi_count_begin(mpi_procid),                 &
              ! srcgrid%mpi_count_end(mpi_procid)

    end if ! if(mpi_procid .ne. mpi_masternode)

    !---------------------------------------------------------------------

    ! Define local variable

    call mpi_reduce(mpi_omlke(1:srcgrid%ncoords),                           &
         & srcgrid_var(1:srcgrid%ncoords),srcgrid%ncoords,mpi_real,         &
         & mpi_sum,mpi_masternode,mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,mpi_masternode,     &
         & mpi_comm_world,mpi_ierror)

    ! Compute local variable

    call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,mpi_masternode,     &
         & mpi_comm_world,mpi_ierror)    

    ! Define local variable

    grid%array = dstgrid_var

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(grid%array,dstgrid%ncoords,mpi_real,mpi_masternode,      &
         & mpi_comm_world,mpi_ierror) 

    !---------------------------------------------------------------------

    ! Deallocate memory for local variables

    if(allocated(srcgrid_var)) deallocate(srcgrid_var)
    if(allocated(dstgrid_var)) deallocate(dstgrid_var)
    if(allocated(mpi_omlke))   deallocate(mpi_omlke)

    !=====================================================================
    
  end subroutine calculate_oceanmixedlayerkineticenergy

  !=======================================================================

  ! calculate_oceanmixedlayeruvelocity.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_oceanmixedlayeruvelocity(hycom,grid)

    ! Define variables passed to routine

    type(hycomvariables)                                           :: hycom
    type(variable_info)                                            :: grid

    ! Define variables computed within routine

    real(r_kind),              dimension(:),           allocatable :: mpi_omlu
    real(r_kind)                                                   :: uvel
    real(r_kind)                                                   :: vvel
    real(r_kind)                                                   :: z
    real(r_kind)                                                   :: omlu
    real(r_kind)                                                   :: omlv
    real(r_kind)                                                   :: omld
    real(r_kind)                                                   :: dz

    ! Define counting variables

    integer                                                        :: i, j, k
    integer                                                        :: count

    !=====================================================================

    ! Initialize local variables

    call init_constants_derived()

    !---------------------------------------------------------------------

    ! Allocate memory for local variables

    if(.not. allocated(srcgrid_var))                                       &
         & allocate(srcgrid_var(srcgrid%ncoords))
    if(.not. allocated(dstgrid_var))                                       &
         & allocate(dstgrid_var(dstgrid%ncoords))
    if(.not. allocated(mpi_omlu))                                          &
         & allocate(mpi_omlu(srcgrid%ncoords))

    !---------------------------------------------------------------------

    ! Initialize local variables

    srcgrid_var = 0.0
    mpi_omlu    = 0.0

    ! If on slave (compute) node (task), receive variables, compute
    ! variables, and send variables to master (root) node (task)
    
    if(mpi_procid .ne. mpi_masternode) then

       ! Loop through local variable and proceed accordingly
       
       do i = srcgrid%mpi_count_begin(mpi_procid),                         &
            & srcgrid%mpi_count_end(mpi_procid)

          ! Initialize local variable

          count = 0

          ! Check local variable and proceed accordingly

          if(hycom%mix_dpth(i) .lt. 1.e20) then

             ! Define local variables

             omld = hycom%mix_dpth(i)/9860.0
             omlu = 0.0
             omlv = 0.0
             dz   = depth_integral_dz

             ! Initialize local variable

             z = omld

             ! Loop through local variable

             do while(z .ge. 0.0)

                ! Compute local variables

                call calculate_variableatdepth(ocean_hycom_zdim,            &
                     & hycom%depth(i,1:ocean_hycom_zdim),                   &
                     & hycom%uvel(i,1:ocean_hycom_zdim),z,uvel)
                omlu = omlu + uvel
                call calculate_variableatdepth(ocean_hycom_zdim,            &
                     & hycom%depth(i,1:ocean_hycom_zdim),                   &
                     & hycom%vvel(i,1:ocean_hycom_zdim),z,vvel)
                omlv = omlv + vvel
                
                ! Update local variables

                z     = z - dz
                count = count + 1
                
             end do ! do while(z .ge. 0.0)

          end if ! if(hycom%mix_dpth(i) .lt. 1.e20)

          ! Define local variable

          if(count .ge. 1) omlu = omlu/real(count)
          if(count .ge. 1) omlv = omlv/real(count)
          if(count .eq. 0) omlu = spval
          if(count .eq. 0) omlv = spval
          
          ! Check local variable and proceed accordingly

          if(omlu .ne. spval .and. omlv .ne. spval) then

             ! Compute local variable

             mpi_omlu(i) = omlu*cos(srcgrid%rotang(i)) +                    &
                  & omlv*sin(-srcgrid%rotang(i))

          else   ! if(omlu .ne. spval .and. omlv .ne. spval)

             ! Define local variable

             mpi_omlu(i) = spval

          end if ! if(omlu .ne. spval .and. omlv .ne. spval)

       end do ! do i = srcgrid%mpi_count_begin(mpi_procid),                 &
              ! srcgrid%mpi_count_end(mpi_procid)

    end if ! if(mpi_procid .ne. mpi_masternode)

    !---------------------------------------------------------------------

    ! Define local variable

    call mpi_reduce(mpi_omlu(1:srcgrid%ncoords),                            &
         & srcgrid_var(1:srcgrid%ncoords),srcgrid%ncoords,mpi_real,         &
         & mpi_sum,mpi_masternode,mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,mpi_masternode,     &
         & mpi_comm_world,mpi_ierror)

    ! Compute local variable

    call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,mpi_masternode,     &
         & mpi_comm_world,mpi_ierror)    

    ! Define local variable

    grid%array = dstgrid_var

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(grid%array,dstgrid%ncoords,mpi_real,mpi_masternode,      &
         & mpi_comm_world,mpi_ierror) 

    !---------------------------------------------------------------------

    ! Deallocate memory for local variables

    if(allocated(srcgrid_var)) deallocate(srcgrid_var)
    if(allocated(dstgrid_var)) deallocate(dstgrid_var)
    if(allocated(mpi_omlu))    deallocate(mpi_omlu)

    !=====================================================================
    
  end subroutine calculate_oceanmixedlayeruvelocity

  !=======================================================================

  ! calculate_oceanmixedlayervvelocity.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_oceanmixedlayervvelocity(hycom,grid)

    ! Define variables passed to routine

    type(hycomvariables)                                           :: hycom
    type(variable_info)                                            :: grid

    ! Define variables computed within routine

    real(r_kind),              dimension(:),           allocatable :: mpi_omlv
    real(r_kind)                                                   :: uvel
    real(r_kind)                                                   :: vvel
    real(r_kind)                                                   :: z
    real(r_kind)                                                   :: omlu
    real(r_kind)                                                   :: omlv
    real(r_kind)                                                   :: omld
    real(r_kind)                                                   :: dz

    ! Define counting variables

    integer                                                        :: i, j, k
    integer                                                        :: count

    !=====================================================================

    ! Initialize local variables

    call init_constants_derived()

    !---------------------------------------------------------------------

    ! Allocate memory for local variables

    if(.not. allocated(srcgrid_var))                                       &
         & allocate(srcgrid_var(srcgrid%ncoords))
    if(.not. allocated(dstgrid_var))                                       &
         & allocate(dstgrid_var(dstgrid%ncoords))
    if(.not. allocated(mpi_omlv))                                          &
         & allocate(mpi_omlv(srcgrid%ncoords))

    !---------------------------------------------------------------------

    ! Initialize local variables

    srcgrid_var = 0.0
    mpi_omlv    = 0.0

    ! If on slave (compute) node (task), receive variables, compute
    ! variables, and send variables to master (root) node (task)
    
    if(mpi_procid .ne. mpi_masternode) then

       ! Loop through local variable and proceed accordingly
       
       do i = srcgrid%mpi_count_begin(mpi_procid),                         &
            & srcgrid%mpi_count_end(mpi_procid)

          ! Initialize local variable

          count = 0

          ! Check local variable and proceed accordingly

          if(hycom%mix_dpth(i) .lt. 1.e20) then

             ! Define local variables

             omld = hycom%mix_dpth(i)/9860.0
             omlu = 0.0
             omlv = 0.0
             dz   = depth_integral_dz

             ! Initialize local variable

             z = omld

             ! Loop through local variable

             do while(z .ge. 0.0)

                ! Compute local variables

                call calculate_variableatdepth(ocean_hycom_zdim,            &
                     & hycom%depth(i,1:ocean_hycom_zdim),                   &
                     & hycom%uvel(i,1:ocean_hycom_zdim),z,uvel)
                omlu = omlu + uvel
                call calculate_variableatdepth(ocean_hycom_zdim,            &
                     & hycom%depth(i,1:ocean_hycom_zdim),                   &
                     & hycom%vvel(i,1:ocean_hycom_zdim),z,vvel)
                omlv = omlv + vvel
                
                ! Update local variables

                z     = z - dz
                count = count + 1
                
             end do ! do while(z .ge. 0.0)

          end if ! if(hycom%mix_dpth(i) .lt. 1.e20)

          ! Define local variable

          if(count .ge. 1) omlu = omlu/real(count)
          if(count .ge. 1) omlv = omlv/real(count)
          if(count .eq. 0) omlu = spval
          if(count .eq. 0) omlv = spval
          
          ! Check local variable and proceed accordingly

          if(omlu .ne. spval .and. omlv .ne. spval) then

             ! Compute local variable

             mpi_omlv(i) = omlv*cos(srcgrid%rotang(i)) -                    &
                  & omlu*sin(-srcgrid%rotang(i))

          else   ! if(omlu .ne. spval .and. omlv .ne. spval)

             ! Define local variable

             mpi_omlv(i) = spval

          end if ! if(omlu .ne. spval .and. omlv .ne. spval)

       end do ! do i = srcgrid%mpi_count_begin(mpi_procid),                 &
              ! srcgrid%mpi_count_end(mpi_procid)

    end if ! if(mpi_procid .ne. mpi_masternode)

    !---------------------------------------------------------------------

    ! Define local variable

    call mpi_reduce(mpi_omlv(1:srcgrid%ncoords),                            &
         & srcgrid_var(1:srcgrid%ncoords),srcgrid%ncoords,mpi_real,         &
         & mpi_sum,mpi_masternode,mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,mpi_masternode,     &
         & mpi_comm_world,mpi_ierror)

    ! Compute local variable

    call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,mpi_masternode,     &
         & mpi_comm_world,mpi_ierror)    

    ! Define local variable

    grid%array = dstgrid_var

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(grid%array,dstgrid%ncoords,mpi_real,mpi_masternode,      &
         & mpi_comm_world,mpi_ierror) 

    !---------------------------------------------------------------------

    ! Deallocate memory for local variables

    if(allocated(srcgrid_var)) deallocate(srcgrid_var)
    if(allocated(dstgrid_var)) deallocate(dstgrid_var)
    if(allocated(mpi_omlv))    deallocate(mpi_omlv)

    !=====================================================================
    
  end subroutine calculate_oceanmixedlayervvelocity

  !=======================================================================

  ! calculate_oceanmixedlayersalinity.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_oceanmixedlayersalinity(hycom,grid)

    ! Define variables passed to routine

    type(hycomvariables)                                           :: hycom
    type(variable_info)                                            :: grid

    ! Define variables computed within routine

    real(r_kind),              dimension(:),           allocatable :: mpi_omls
    real(r_kind)                                                   :: salin
    real(r_kind)                                                   :: z
    real(r_kind)                                                   :: omls
    real(r_kind)                                                   :: omld
    real(r_kind)                                                   :: dz

    ! Define counting variables

    integer                                                        :: i, j, k
    integer                                                        :: count

    !=====================================================================

    ! Initialize local variables

    call init_constants_derived()

    !---------------------------------------------------------------------

    ! Allocate memory for local variables

    if(.not. allocated(srcgrid_var))                                       &
         & allocate(srcgrid_var(srcgrid%ncoords))
    if(.not. allocated(dstgrid_var))                                       &
         & allocate(dstgrid_var(dstgrid%ncoords))
    if(.not. allocated(mpi_omls))                                          &
         & allocate(mpi_omls(srcgrid%ncoords))

    !---------------------------------------------------------------------

    ! Initialize local variables

    srcgrid_var = 0.0
    mpi_omls    = 0.0

    ! If on slave (compute) node (task), receive variables, compute
    ! variables, and send variables to master (root) node (task)
    
    if(mpi_procid .ne. mpi_masternode) then

       ! Loop through local variable and proceed accordingly
       
       do i = srcgrid%mpi_count_begin(mpi_procid),                         &
            & srcgrid%mpi_count_end(mpi_procid)

          ! Initialize local variable

          count = 0

          ! Check local variable and proceed accordingly

          if(hycom%mix_dpth(i) .lt. 1.e20) then

             ! Define local variables

             omld = hycom%mix_dpth(i)/9860.0
             omls = 0.0
             dz   = depth_integral_dz

             ! Initialize local variable

             z = omld

             ! Loop through local variable

             do while(z .ge. 0.0)

                ! Compute local variables

                call calculate_variableatdepth(ocean_hycom_zdim,            &
                     & hycom%depth(i,1:ocean_hycom_zdim),                   &
                     & hycom%salin(i,1:ocean_hycom_zdim),z,salin)
                omls = omls + salin
                
                ! Update local variables

                z     = z - dz
                count = count + 1
                
             end do ! do while(z .ge. 0.0)

          end if ! if(hycom%mix_dpth(i) .lt. 1.e20)

          ! Define local variable

          if(count .ge. 1) omls = omls/real(count)
          if(count .eq. 0) omls = spval
          mpi_omls(i) = omls

       end do ! do i = srcgrid%mpi_count_begin(mpi_procid),                 &
              ! srcgrid%mpi_count_end(mpi_procid)

    end if ! if(mpi_procid .ne. mpi_masternode)

    !---------------------------------------------------------------------

    ! Define local variable

    call mpi_reduce(mpi_omls(1:srcgrid%ncoords),                            &
         & srcgrid_var(1:srcgrid%ncoords),srcgrid%ncoords,mpi_real,         &
         & mpi_sum,mpi_masternode,mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,mpi_masternode,     &
         & mpi_comm_world,mpi_ierror)

    ! Compute local variable

    call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,mpi_masternode,     &
         & mpi_comm_world,mpi_ierror)    

    ! Define local variable

    grid%array = dstgrid_var

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(grid%array,dstgrid%ncoords,mpi_real,mpi_masternode,      &
         & mpi_comm_world,mpi_ierror) 

    !---------------------------------------------------------------------

    ! Deallocate memory for local variables

    if(allocated(srcgrid_var)) deallocate(srcgrid_var)
    if(allocated(dstgrid_var)) deallocate(dstgrid_var)
    if(allocated(mpi_omls))    deallocate(mpi_omls)

    !=====================================================================
    
  end subroutine calculate_oceanmixedlayersalinity

  !=======================================================================

  ! calculate_oceanmixedlayertemperature.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_oceanmixedlayertemperature(hycom,grid)

    ! Define variables passed to routine

    type(hycomvariables)                                           :: hycom
    type(variable_info)                                            :: grid

    ! Define variables computed within routine

    real(r_kind),              dimension(:),           allocatable :: mpi_omlt
    real(r_kind)                                                   :: temp
    real(r_kind)                                                   :: z
    real(r_kind)                                                   :: omlt
    real(r_kind)                                                   :: omld
    real(r_kind)                                                   :: dz

    ! Define counting variables

    integer                                                        :: i, j, k
    integer                                                        :: count

    !=====================================================================

    ! Initialize local variables

    call init_constants_derived()

    !---------------------------------------------------------------------

    ! Allocate memory for local variables

    if(.not. allocated(srcgrid_var))                                       &
         & allocate(srcgrid_var(srcgrid%ncoords))
    if(.not. allocated(dstgrid_var))                                       &
         & allocate(dstgrid_var(dstgrid%ncoords))
    if(.not. allocated(mpi_omlt))                                          &
         & allocate(mpi_omlt(srcgrid%ncoords))

    !---------------------------------------------------------------------

    ! Initialize local variables

    srcgrid_var = 0.0
    mpi_omlt    = 0.0

    ! If on slave (compute) node (task), receive variables, compute
    ! variables, and send variables to master (root) node (task)
    
    if(mpi_procid .ne. mpi_masternode) then

       ! Loop through local variable and proceed accordingly
       
       do i = srcgrid%mpi_count_begin(mpi_procid),                         &
            & srcgrid%mpi_count_end(mpi_procid)

          ! Initialize local variable

          count = 0

          ! Check local variable and proceed accordingly

          if(hycom%mix_dpth(i) .lt. 1.e20) then

             ! Define local variables

             omld = hycom%mix_dpth(i)/9860.0
             omlt = 0.0
             dz   = depth_integral_dz

             ! Initialize local variable

             z = omld

             ! Loop through local variable

             do while(z .ge. 0.0)

                ! Compute local variables

                call calculate_variableatdepth(ocean_hycom_zdim,            &
                     & hycom%depth(i,1:ocean_hycom_zdim),                   &
                     & hycom%temp(i,1:ocean_hycom_zdim),z,temp)
                omlt = omlt + temp
                
                ! Update local variables

                z     = z - dz
                count = count + 1
                
             end do ! do while(z .ge. 0.0)

          end if ! if(hycom%mix_dpth(i) .lt. 1.e20)

          ! Define local variable

          if(count .ge. 1) omlt = (omlt/real(count)) + t0c
          if(count .eq. 0) omlt = spval
          mpi_omlt(i) = omlt

       end do ! do i = srcgrid%mpi_count_begin(mpi_procid),                 &
              ! srcgrid%mpi_count_end(mpi_procid)

    end if ! if(mpi_procid .ne. mpi_masternode)

    !---------------------------------------------------------------------

    ! Define local variable

    call mpi_reduce(mpi_omlt(1:srcgrid%ncoords),                            &
         & srcgrid_var(1:srcgrid%ncoords),srcgrid%ncoords,mpi_real,         &
         & mpi_sum,mpi_masternode,mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,mpi_masternode,     &
         & mpi_comm_world,mpi_ierror)

    ! Compute local variable

    call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,mpi_masternode,     &
         & mpi_comm_world,mpi_ierror)    

    ! Define local variable

    grid%array = dstgrid_var

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(grid%array,dstgrid%ncoords,mpi_real,mpi_masternode,      &
         & mpi_comm_world,mpi_ierror) 

    !---------------------------------------------------------------------

    ! Deallocate memory for local variables

    if(allocated(srcgrid_var)) deallocate(srcgrid_var)
    if(allocated(dstgrid_var)) deallocate(dstgrid_var)
    if(allocated(mpi_omlt))    deallocate(mpi_omlt)

    !=====================================================================

  end subroutine calculate_oceanmixedlayertemperature

  !=======================================================================

  ! calculate_oceanmixedlayerdepth.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_oceanmixedlayerdepth(hycom,grid)

    ! Define variables passed to routine

    type(hycomvariables)                                           :: hycom
    type(variable_info)                                            :: grid

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(srcgrid_var))                                       &
         & allocate(srcgrid_var(srcgrid%ncoords))
    if(.not. allocated(dstgrid_var))                                       &
         & allocate(dstgrid_var(dstgrid%ncoords))

    !---------------------------------------------------------------------

    ! Define local variable

    srcgrid_var(1:srcgrid%ncoords) = hycom%mix_dpth(1:srcgrid%ncoords)

    ! Check local variable and proceed accordingly

    where(srcgrid_var(1:srcgrid%ncoords) .lt. 1.e20)                       &
         & srcgrid_var(1:srcgrid%ncoords) = srcgrid_var(1:srcgrid%ncoords) &
         & /9860.0

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,mpi_masternode,    &
         & mpi_comm_world,mpi_ierror)

    ! Compute local variable

    call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,mpi_masternode,    &
         & mpi_comm_world,mpi_ierror)    

    ! Define local variable

    grid%array = dstgrid_var

    ! Update local variable accordingly

    where(grid%array .gt. 1e20) grid%array = spval

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(grid%array,dstgrid%ncoords,mpi_real,mpi_masternode,     &
         & mpi_comm_world,mpi_ierror) 

    !---------------------------------------------------------------------

    ! Deallocate memory for local variables

    if(allocated(srcgrid_var)) deallocate(srcgrid_var)
    if(allocated(dstgrid_var)) deallocate(dstgrid_var)

    !=====================================================================

  end subroutine calculate_oceanmixedlayerdepth

  !=======================================================================

  ! calculate_seasurfacetemperature.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_seasurfacetemperature(hycom,grid)

    ! Define variables passed to routine

    type(hycomvariables)                                           :: hycom
    type(variable_info)                                            :: grid

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(srcgrid_var))                                       &
         & allocate(srcgrid_var(srcgrid%ncoords))
    if(.not. allocated(dstgrid_var))                                       &
         & allocate(dstgrid_var(dstgrid%ncoords))

    !---------------------------------------------------------------------

    ! Define local variable

    srcgrid_var(1:srcgrid%ncoords) = hycom%temp(1:srcgrid%ncoords,1)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,mpi_masternode,    &
         & mpi_comm_world,mpi_ierror)

    ! Compute local variable

    call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,mpi_masternode,    &
         & mpi_comm_world,mpi_ierror)    

    ! Define local variable

    grid%array = dstgrid_var

    ! Update local variable accordingly

    where(grid%array .gt. 1e20)  grid%array = spval
    where(grid%array .ne. spval) grid%array = grid%array + t0c

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(grid%array,dstgrid%ncoords,mpi_real,mpi_masternode,     &
         & mpi_comm_world,mpi_ierror) 

    !---------------------------------------------------------------------

    ! Deallocate memory for local variables

    if(allocated(srcgrid_var)) deallocate(srcgrid_var)
    if(allocated(dstgrid_var)) deallocate(dstgrid_var)

    !=====================================================================

  end subroutine calculate_seasurfacetemperature

  !=======================================================================

  ! calculate_seasurfaceheight.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_seasurfaceheight(hycom,grid)

    ! Define variables passed to routine

    type(hycomvariables)                                           :: hycom
    type(variable_info)                                            :: grid

    !=====================================================================

    ! Allocate memory for local variables

    if(.not. allocated(srcgrid_var))                                       &
         & allocate(srcgrid_var(srcgrid%ncoords))
    if(.not. allocated(dstgrid_var))                                       &
         & allocate(dstgrid_var(dstgrid%ncoords))

    !---------------------------------------------------------------------

    ! Define local variable

    srcgrid_var(1:srcgrid%ncoords) = hycom%srfhgt(1:srcgrid%ncoords)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,mpi_masternode,    &
         & mpi_comm_world,mpi_ierror)

    ! Compute local variable

    call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,mpi_masternode,    &
         & mpi_comm_world,mpi_ierror)    

    ! Define local variable

    grid%array = dstgrid_var

    ! Update local variable accordingly

    where(grid%array .gt. 1e20)  grid%array = spval
    where(grid%array .ne. spval) grid%array = grid%array/100.0

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(grid%array,dstgrid%ncoords,mpi_real,mpi_masternode,     &
         & mpi_comm_world,mpi_ierror) 

    !---------------------------------------------------------------------

    ! Deallocate memory for local variables

    if(allocated(srcgrid_var)) deallocate(srcgrid_var)
    if(allocated(dstgrid_var)) deallocate(dstgrid_var)

    !=====================================================================

  end subroutine calculate_seasurfaceheight

  !=======================================================================

  ! calculate_oceanheatcontent.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_oceanheatcontent(hycom,grid)

    ! Define variables passed to routine

    type(hycomvariables)                                           :: hycom
    type(variable_info)                                            :: grid

    ! Define variables computed within routine

    real(r_kind),              dimension(:),           allocatable :: mpi_ohc
    real(r_kind)                                                   :: isobathy
    real(r_kind)                                                   :: temp
    real(r_kind)                                                   :: salin
    real(r_kind)                                                   :: z
    real(r_kind)                                                   :: rho
    real(r_kind)                                                   :: cp
    real(r_kind)                                                   :: ohc
    real(r_kind)                                                   :: dz

    ! Define counting variables

    integer                                                        :: i, j, k

    !=====================================================================

    ! Initialize local variables

    call init_constants_derived()

    !---------------------------------------------------------------------

    ! Allocate memory for local variables

    if(.not. allocated(srcgrid_var))                                       &
         & allocate(srcgrid_var(srcgrid%ncoords))
    if(.not. allocated(dstgrid_var))                                       &
         & allocate(dstgrid_var(dstgrid%ncoords))
    if(.not. allocated(mpi_ohc))                                           &
         & allocate(mpi_ohc(srcgrid%ncoords))

    !---------------------------------------------------------------------

    ! Loop through local variable

    do k = 1, nlev

       ! Initialize local variables

       srcgrid_var = 0.0
       mpi_ohc     = 0.0

       ! If on slave (compute) node (task), receive variables, compute
       ! variables, and send variables to master (root) node (task)
    
       if(mpi_procid .ne. mpi_masternode) then

          ! Loop through local variable and proceed accordingly
       
          do i = srcgrid%mpi_count_begin(mpi_procid),                      &
               & srcgrid%mpi_count_end(mpi_procid)

             ! Define local variable

             ohc = 0.0
             dz  = depth_integral_dz

             ! Check local variable and proceed accordingly

             if(hycom%depth(i,ocean_hycom_zdim) .ge. zlevs(k) .and.         &
                  & hycom%depth(i,ocean_hycom_zdim) .ne. 0.0) then

                ! Initialize local variable

                z = zlevs(k)

                ! Loop through local variable

                do while(z .ge. 0.0)

                   ! Compute local variables

                   call calculate_variableatdepth(ocean_hycom_zdim,         &
                        & hycom%depth(i,1:ocean_hycom_zdim),                &
                        & hycom%temp(i,1:ocean_hycom_zdim),z,temp)
                   call calculate_variableatdepth(ocean_hycom_zdim,         &
                        & hycom%depth(i,1:ocean_hycom_zdim),                &
                        & hycom%salin(i,1:ocean_hycom_zdim),z,salin)
                   call calculate_seawaterspecificheatcapacity(temp+t0c,    &
                        & salin,cp)
                   call calculate_seawaterdensity(temp+t0c,salin,rho)
                   ohc = ohc + rho*cp*temp
                
                   ! Update local variable

                   z = z - dz
                
                end do ! do while(z .ge. 0.0)

             end if ! if(hycom%depth(i,ocean_hycom_zdim) .ge. zlevs(k)      &
                    ! .and. hycom%depth(i,ocean_hycom_zdim) .ne. 0.0)

             ! Define local variable

             if(ohc .le. 0.0) ohc = spval
             mpi_ohc(i) = ohc

          end do ! do i = srcgrid%mpi_count_begin(mpi_procid),              &
                 ! srcgrid%mpi_count_end(mpi_procid)

       end if ! if(mpi_procid .ne. mpi_masternode)

       ! Define local variable

       call mpi_reduce(mpi_ohc(1:srcgrid%ncoords),                          &
            & srcgrid_var(1:srcgrid%ncoords),srcgrid%ncoords,mpi_real,      &
            & mpi_sum,mpi_masternode,mpi_comm_world,mpi_ierror)

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
       
       call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,mpi_masternode,  &
            & mpi_comm_world,mpi_ierror)

       ! Compute local variable

       call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

       ! Broadcast all necessary variables to slave (compute) nodes
       ! (tasks)
       
       call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,mpi_masternode,  &
            & mpi_comm_world,mpi_ierror)    

       ! Define local variable
       
       grid%array((k-1)*dstgrid%ncoords+1:(k*dstgrid%ncoords)) =            &
            & dstgrid_var(1:dstgrid%ncoords)

    end do ! k = 1, nlev

    !---------------------------------------------------------------------

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(grid%array,(dstgrid%ncoords*nlev),mpi_real,              &
         & mpi_masternode,mpi_comm_world,mpi_ierror) 

    !---------------------------------------------------------------------

    ! Deallocate memory for local variables

    if(allocated(srcgrid_var)) deallocate(srcgrid_var)
    if(allocated(dstgrid_var)) deallocate(dstgrid_var)
    if(allocated(mpi_ohc))     deallocate(mpi_ohc)

    !=====================================================================

  end subroutine calculate_oceanheatcontent

  !=======================================================================

  ! calculate_tropicalcycloneheatpotential.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_tropicalcycloneheatpotential(hycom,grid)

    ! Define variables passed to routine

    type(hycomvariables)                                           :: hycom
    type(variable_info)                                            :: grid

    ! Define variables computed within routine

    real(r_kind),              dimension(:),           allocatable :: mpi_tchp
    real(r_kind)                                                   :: t26c
    real(r_kind)                                                   :: z26c
    real(r_kind)                                                   :: temp
    real(r_kind)                                                   :: salin
    real(r_kind)                                                   :: z
    real(r_kind)                                                   :: rho
    real(r_kind)                                                   :: cp
    real(r_kind)                                                   :: tchp
    real(r_kind)                                                   :: dz

    ! Define counting variables

    integer                                                        :: i, j, k

    !=====================================================================

    ! Initialize local variables

    call init_constants_derived()

    !---------------------------------------------------------------------

    ! Allocate memory for local variables

    if(.not. allocated(srcgrid_var))                                       &
         & allocate(srcgrid_var(srcgrid%ncoords))
    if(.not. allocated(dstgrid_var))                                       &
         & allocate(dstgrid_var(dstgrid%ncoords))
    if(.not. allocated(mpi_tchp))                                          &
         & allocate(mpi_tchp(srcgrid%ncoords))

    !---------------------------------------------------------------------

    ! Initialize local variables

    srcgrid_var = 0.0
    mpi_tchp    = 0.0

    ! If on slave (compute) node (task), receive variables, compute
    ! variables, and send variables to master (root) node (task)
    
    if(mpi_procid .ne. mpi_masternode) then

       ! Loop through local variable and proceed accordingly
       
       do i = srcgrid%mpi_count_begin(mpi_procid),                         &
            & srcgrid%mpi_count_end(mpi_procid)

          ! Define local variable

          t26c = 26.0
          tchp = 0.0
          dz   = depth_integral_dz

          ! Compute local variable

          call calculate_isothermdepth(ocean_hycom_zdim,                   &
               & hycom%temp(i,1:ocean_hycom_zdim),                         &
               & hycom%depth(i,1:ocean_hycom_zdim),t26c,z26c)

          ! Check local variable and proceed accordingly

          if(z26c .ne. spval) then

             ! Initialize local variable

             z = z26c

             ! Loop through local variable

             do while(z .ge. 0.0)

                ! Compute local variables

                call calculate_variableatdepth(ocean_hycom_zdim,            &
                     & hycom%depth(i,1:ocean_hycom_zdim),                   &
                     & hycom%temp(i,1:ocean_hycom_zdim),z,temp)
                call calculate_variableatdepth(ocean_hycom_zdim,            &
                     & hycom%depth(i,1:ocean_hycom_zdim),                   &
                     & hycom%salin(i,1:ocean_hycom_zdim),z,salin)
                call calculate_seawaterspecificheatcapacity(temp+t0c,       &
                     & salin,cp)
                call calculate_seawaterdensity(temp+t0c,salin,rho)
                tchp = tchp + rho*cp*(temp - t26c)
                
                ! Update local variable

                z = z - dz
                
             end do ! do while(z .ge. 0.0)

          end if ! if(z26c .ne. spval)

          ! Define local variable

          if(tchp .le. 0.0) tchp = spval
          mpi_tchp(i) = tchp

       end do ! do i = srcgrid%mpi_count_begin(mpi_procid),                 &
              ! srcgrid%mpi_count_end(mpi_procid)

    end if ! if(mpi_procid .ne. mpi_masternode)

    !---------------------------------------------------------------------

    ! Define local variable

    call mpi_reduce(mpi_tchp(1:srcgrid%ncoords),                            &
         & srcgrid_var(1:srcgrid%ncoords),srcgrid%ncoords,mpi_real,         &
         & mpi_sum,mpi_masternode,mpi_comm_world,mpi_ierror)

    !---------------------------------------------------------------------

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(srcgrid_var,srcgrid%ncoords,mpi_real,mpi_masternode,     &
         & mpi_comm_world,mpi_ierror)

    ! Compute local variable

    call interpolation_barnes_analysis_mpi(srcgrid,dstgrid)

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(dstgrid_var,dstgrid%ncoords,mpi_real,mpi_masternode,     &
         & mpi_comm_world,mpi_ierror)    

    ! Define local variable

    grid%array = dstgrid_var

    ! Broadcast all necessary variables to slave (compute) nodes
    ! (tasks)
       
    call mpi_bcast(grid%array,dstgrid%ncoords,mpi_real,mpi_masternode,      &
         & mpi_comm_world,mpi_ierror) 

    !---------------------------------------------------------------------

    ! Deallocate memory for local variables

    if(allocated(srcgrid_var)) deallocate(srcgrid_var)
    if(allocated(dstgrid_var)) deallocate(dstgrid_var)
    if(allocated(mpi_tchp))    deallocate(mpi_tchp)

    !=====================================================================

  end subroutine calculate_tropicalcycloneheatpotential

  !=======================================================================

  ! calculate_seawaterdensity.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_seawaterdensity(temp,salin,rho)

    ! Define variables passed to routine

    real(r_kind)                                             :: temp
    real(r_kind)                                             :: salin

    ! Define variables returned by routine

    real(r_kind)                                             :: rho

    ! Define variables computed within routine

    real(r_kind)                                             :: a(5)
    real(r_kind)                                             :: b(5)
    real(r_kind)                                             :: temp_c
    real(r_kind)                                             :: salin_s
    real(r_kind)                                             :: rho_w
    real(r_kind)                                             :: d_rho

    !=====================================================================

    ! Initialize local variables

    call init_constants_derived()

    !---------------------------------------------------------------------

    ! Define local variables

    a(1) = 9.9992293295e+02  ; b(1) = 8.0200240891e+02 
    a(2) = 2.0341179217e-02  ; b(2) = -2.0005183488e+00  
    a(3) = -6.1624591598e-03 ; b(3) = 1.6771024982e-02  
    a(4) = 2.2614664708e-05  ; b(4) = -3.0600536746e-05 
    a(5) = -4.6570659168e-08 ; b(5) = -1.6132224742e-05

    ! Compute local variables

    temp_c  = temp - t0c
    salin_s = salin/1000.0

    !---------------------------------------------------------------------

    ! Compute local variables

    rho_w = a(1) + a(2)*temp_c + a(3)*temp_c**2.0 + a(4)*temp_c**3.0 +     &
         & a(5)*temp_c**4.0
    d_rho = b(1)*salin_s + b(2)*salin_s*temp_c +                           &
         & b(3)*salin_s*(temp_c**2.0) + b(4)*salin_s*(temp_c**3.0)         &
         & + b(5)*(salin_s**2.0)*(temp_c**2.0)
    rho   = rho_w + d_rho

    !=====================================================================

  end subroutine calculate_seawaterdensity

  !=======================================================================

  ! calculate_isothermdepth.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_isothermdepth(zdim,temp,depth,isotherm,isodepth)

    ! Define array dimension variables

    integer                                                 :: zdim

    ! Define variables passed to routine

    real(r_kind),                  dimension(zdim)          :: depth
    real(r_kind),                  dimension(zdim)          :: temp
    real(r_kind)                                            :: isotherm

    ! Define variables returned by routine

    real(r_kind)                                            :: isodepth

    ! Define variables computed within routine

    logical                                                 :: is_depth
    logical                                                 :: is_temp
    real(r_kind)                                            :: w1
    real(r_kind)                                            :: w2

    ! Define counting variables

    integer                                                 :: i, j, k

    !=====================================================================

    ! Initialize local variables

    is_depth = .false.
    is_temp  = .false.
    isodepth = spval

    ! Check local variables and proceed accordingly

    if(maxval(depth(1:zdim)) .gt. 0.0) is_depth = .true.
    if(temp(1) .ge. isotherm)          is_temp  = .true.

    !---------------------------------------------------------------------

    ! Check local variables and proceed accordingly

    if(is_depth .and. is_temp) then

       ! Loop through local variable

       do k = 1, zdim - 1

          ! Check local variable and proceed accordingly

          if(temp(k) .ge. isotherm .and. temp(k+1) .le. isotherm) then

             ! Compute local variables

             w1 = (temp(k+1) - isotherm)/(temp(k+1) - temp(k))
             w2 = 1.0 - w1

             ! Define exit from loop

             goto 1000

          end if ! if(temp(k) .ge. isotherm .and. temp(k+1)               &
                 ! .le. isotherm)

       end do ! do k = 1, zdim

       ! Define exit from loop

1000   continue

       ! Compute local variable

       isodepth = w1*depth(k) + w2*depth(k+1)

    end if ! if(is_depth .and. is_temp)

    !=====================================================================

  end subroutine calculate_isothermdepth

  !=======================================================================

  ! calculate_variableatdepth.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_variableatdepth(zdim,depth,vari,isodepth,isovari)

    ! Define array dimension variables

    integer                                                 :: zdim

    ! Define variables passed to routine

    real(r_kind),                  dimension(zdim)          :: depth
    real(r_kind),                  dimension(zdim)          :: vari
    real(r_kind)                                            :: isodepth

    ! Define variables returned by routine

    real(r_kind)                                            :: isovari

    ! Define variables computed within routine

    logical                                                 :: is_depth
    logical                                                 :: is_vari
    real(r_kind)                                            :: w1
    real(r_kind)                                            :: w2

    ! Define counting variables

    integer                                                 :: i, j, k

    !=====================================================================

    ! Initialize local variables

    is_depth = .false.
    is_vari  = .false.
    isovari  = spval

    ! Check local variables and proceed accordingly

    if(maxval(depth(1:zdim)) .gt. 0.0) is_depth = .true.
    if(vari(1) .lt. 1.e20)             is_vari  = .true.

    !---------------------------------------------------------------------

    ! Check local variables and proceed accordingly

    if(is_depth .and. is_vari) then

       ! Loop through local variable

       do k = 1, zdim - 1

          ! Check local variable and proceed accordingly

          if(depth(k) .le. isodepth .and. depth(k+1) .ge. isodepth) then

             ! Compute local variables

             w1 = (depth(k+1) - isodepth)/(depth(k+1) - depth(k))
             w2 = 1.0 - w1

             ! Define exit from loop

             goto 1000

          end if ! if(depth(k) .le. isodepth .and. depth(k+1)               &
                 ! .ge. isodepth)

       end do ! do k = 1, zdim

       ! Define exit from loop

1000   continue

       ! Compute local variable

       isovari = w1*vari(k) + w2*vari(k+1)

    end if ! if(is_depth .and. is_vari)

    !=====================================================================

  end subroutine calculate_variableatdepth

  !=======================================================================

  ! calculate_seawaterspecificheatcapacity.f90:

  !-----------------------------------------------------------------------

  subroutine calculate_seawaterspecificheatcapacity(temp,salin,cp)

    ! Define variables passed to routine

    real(r_kind)                                             :: temp
    real(r_kind)                                             :: salin

    ! Define variables returned by routine

    real(r_kind)                                             :: cp

    ! Define variables computed within routine

    real(r_kind)                                             :: temp_c
    real(r_kind)                                             :: A
    real(r_kind)                                             :: B
    real(r_kind)                                             :: C
    real(r_kind)                                             :: D

    !=====================================================================

    ! Initialize local variables

    call init_constants_derived()

    !---------------------------------------------------------------------

    ! Compute local variables

    temp_c = temp - t0c
    A      = 4206.8 - 6.6197*salin + 1.2288e-2*(salin**2.0)
    B      = -1.1262 + 5.4178e-2*salin - 2.2719e-4*(salin**2.0)
    C      = 1.2026e-2 - 5.3566e-4*salin + 1.8906e-6*(salin**2.0)
    D      = 6.8777e-7 + 1.517e-6*salin - 4.4268e-9*(salin**2.0)
    cp     = A + B*temp_c + C*(temp_c**2.0) + D*(temp_c**3.0)

    !=====================================================================

  end subroutine calculate_seawaterspecificheatcapacity

  !=======================================================================

  ! grads_interface_define_timestamp.f90:

  !-----------------------------------------------------------------------

  subroutine grads_interface_define_timestamp(grads_timestamp)

    ! Define variable returned by routine

    character(len=15)                                        :: grads_timestamp

    ! Define variables computed within routine

    character(len=3)                                         :: month_char(12)
    logical                                                  :: is_leap
    real(r_kind)                                             :: wnfracday
    real(r_kind)                                             :: month_days(12)
    real(r_kind)                                             :: ryear
    real(r_kind)                                             :: rday
    integer                                                  :: imonth
    integer                                                  :: iday
    integer                                                  :: ihour
    integer                                                  :: iminute

    ! Define counting variables

    integer                                                  :: i, j, k

    !===================================================================== 

    ! Compute local variables

    call wnday(hycom_time,ryear,rday)

    !---------------------------------------------------------------------

    ! Define local variable

    month_char(1)  = 'JAN'
    month_char(2)  = 'FEB'
    month_char(3)  = 'MAR'
    month_char(4)  = 'APR'
    month_char(5)  = 'MAY'
    month_char(6)  = 'JUN'
    month_char(7)  = 'JUL'
    month_char(8)  = 'AUG'
    month_char(9)  = 'SEP'
    month_char(10) = 'OCT'
    month_char(11) = 'NOV'
    month_char(12) = 'DEC'

    ! Initialize local variable

    is_leap = .false.

    ! Check local variable and proceed accordingly

    if(mod(ryear,400.) .eq. 0) is_leap = .true.
    if(mod(ryear,100.) .eq. 0) is_leap = .false.
    if(mod(ryear,4.)   .eq. 0) is_leap = .true.

    ! Check local variable and proceed accordingly

    if(is_leap) then

       ! Define local variable

       month_days(1)  = 1
       month_days(2)  = 32
       month_days(3)  = 61
       month_days(4)  = 92 
       month_days(5)  = 122
       month_days(6)  = 153
       month_days(7)  = 183
       month_days(8)  = 214
       month_days(9)  = 245
       month_days(10) = 275
       month_days(11) = 306
       month_days(12) = 336

    else   ! if(is_leap)

       ! Define local variable

       month_days(1)  = 1
       month_days(2)  = 32
       month_days(3)  = 60
       month_days(4)  = 91 
       month_days(5)  = 121
       month_days(6)  = 152
       month_days(7)  = 182
       month_days(8)  = 213
       month_days(9)  = 244
       month_days(10) = 274
       month_days(11) = 305
       month_days(12) = 335

    end if ! if(is_leap)

    ! Loop through local variable

    do i = 1, 11

       ! Check local variable and proceed accordingly

       if(rday .ge. month_days(i) .and. rday .lt. month_days(i+1)) then

          ! Define local variable

          imonth = i

       end if ! if(rday .ge. month_days(i) .and. rday
              ! .le. month_days(i+1))

    end do ! do i = 1, 11

    ! Compute local variable

    iday = int(rday) - month_days(imonth) + 1
    
    !---------------------------------------------------------------------

    ! Compute local variables
    
    wnfracday = hycom_time - int(hycom_time)
    ihour     = int(wnfracday*24.0)
    iminute   = int(((wnfracday*24.0) - int(wnfracday*24.0))*60.0)

    !---------------------------------------------------------------------

    ! Define local variable

    write(grads_timestamp,500) ihour, iminute, iday, month_char(imonth),   &
         & int(ryear)

    !===================================================================== 

    ! Define format statement
    
500 format(i2.2,":",i2.2,"Z",i2.2,a3,i4)

    !===================================================================== 

  end subroutine grads_interface_define_timestamp

  !=======================================================================

  ! grads_interface_descriptor.f90:

  !-----------------------------------------------------------------------

  subroutine grads_interface_descriptor(grid)

    ! Define variables passed to routine

    type(variable_info),     dimension(number_of_var)        :: grid

    ! Define variables computed within routine

    character(len=50)                                        :: fmt
    character(len=19)                                        :: grid_time
    character(len=15)                                        :: grads_timestamp
    character(len=8)                                         :: grid_output_interval

    ! Define counting variables

    integer                                                  :: i, j, k

    !===================================================================== 

    ! Define local variable

    call grads_interface_define_timestamp(grads_timestamp)

    ! Open external file

    open(999,file='hycomtograds.ctl',form='formatted')

    ! Write values to external file

    write(999,500) 'hycomtograds.bin'
    write(999,501)

    ! Write values to external file

    write(999,502) nlon, rlon_min, interp_dx
    write(999,503) nlat, rlat_min, interp_dy

    ! Define local variable

    write(fmt,'("("i"(f13.5))")') nlev

    ! Write values to external file

    if(.not. is_zinterp .and. .not. is_isointerp)                         &
         & write(999,504) nlev
    if((is_zinterp .or. is_isointerp).and. nlev .gt. 1)                   &
         & write(999,506) nlev

    ! Check local variable and proceed accordingly

    if(nlev .gt. 1) then

       ! Write values to external file

       if(is_zinterp)   write(999,(adjustl(fmt)))                         &
            & (zlevs(k),k = 1, nlev)
       if(is_isointerp) write(999,(adjustl(fmt)))                         &
            & (tlevs(k),k = 1, nlev)

    else   ! if(nlev .gt. 1) 

       ! Write values to external file

       if(is_zinterp)   write(999,507) zlevs(1)
       if(is_isointerp) write(999,507) tlevs(1) 

    end if ! if(nlev .gt. 1) 

    ! Write values to external file

    write(999,509) grads_timestamp

    ! Write values to external file

    write(999,510) number_of_var
    do i = 1, number_of_var
       write(999,511) grid(i)%grads_id, grid(i)%zdim, 99,                  &
            & grid(i)%grads_string
    end do ! do i = 1, number_of_var
    write(999,512)

    ! Close external file

    close(999)

    !---------------------------------------------------------------------

    ! Define format statements

500 format('dset ', a)
501 format('undef 1.e30')
502 format('xdef ', i6, ' linear ', f13.5,1x,f13.5)
503 format('ydef ', i6, ' linear ', f13.5,1x,f13.5)
504 format('zdef ', i6, ' linear   0.0  1.0')
505 format(i6)
506 format('zdef ', i6, ' levels')
507 format('zdef 1 levels ',f13.5)
509 format('tdef 1 linear ', a15, ' 1hr')
510 format('vars ', i4)
511 format(a,1x,i2,1x,i6,1x,a)
512 format('endvars')

    !===================================================================== 

  end subroutine grads_interface_descriptor

  !=======================================================================

  ! grads_interface_binary.f90:

  !-----------------------------------------------------------------------

  subroutine grads_interface_binary(grid)

    ! Define variables passed to routine

    type(variable_info),     dimension(number_of_var)        :: grid

    ! Define counting variables

    integer                                                  :: i, j, k

    !===================================================================== 
 
    ! Open external file

    open(999,file='hycomtograds.bin',form='unformatted',status='unknown',  &
         & recordtype='stream')

    ! Loop through each variable and proceed accordingly

    do i = 1, number_of_var

       ! Write variable to external file

       write(999) grid(i)%array
    
       ! Deallocate memory for local variable

       if(allocated(grid(i)%array)) deallocate(grid(i)%array)

    end do !  do i = 1, number_of_var

    ! Close external file

    close(999)

    !===================================================================== 

  end subroutine grads_interface_binary

  !=======================================================================

end module variable_interface
