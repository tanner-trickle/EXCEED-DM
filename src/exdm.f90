program exdm
    ! EXCEED-DM : EXtended Calculation of Electronic Excitations for Direct detection of Dark Matter
    ! 
    ! The main program.

    use iso_fortran_env
    use mpi

    use hdf5_utils

    use logger_util, only: logger_t
    use timer_util, only: timer_t

    use exdm_inputs_type, only: exdm_inputs_t
    use exdm_elec_config_type, only: exdm_elec_config_t

    ! calculations
    use exdm_calc_binned_scatter_rate, only: exdm_binned_scatter_rate
    use exdm_calc_absorption_rate, only: exdm_absorption_rate
    use exdm_calc_dielectric, only: exdm_dielectric

    implicit none

    ! MPI variables
    integer :: proc_id
        ! processor ID
    integer :: root_proc_id = 0
        ! root process ID
    integer :: n_proc
        ! number of processors
    integer :: mpi_err
        ! MPI error code

    ! hdf5_utils
    integer(HID_T) :: file_id

    ! utilities
    type(logger_t) :: logger
        ! Handles logging messages to user
    type(timer_t) :: timer
        ! Handles timing the program

    ! EXDM specific types
    type(exdm_inputs_t) :: exdm_inputs 
    type(exdm_elec_config_t) :: exdm_elec_config

    ! initialize MPI
    call MPI_INIT(mpi_err)
    call MPI_COMM_RANK(MPI_COMM_WORLD, proc_id, mpi_err)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, n_proc, mpi_err)

    if ( proc_id == root_proc_id ) then

        ! start up message
        call logger%startup_message('EXCEED-DM', _cmake_version, n_proc)

        ! start timer
        call timer%start()

    end if

    ! load inputs
    call exdm_inputs%load(proc_id, root_proc_id, logger)

    ! create output file
    if ( exdm_inputs%control%verbose ) then
        print*, 'Creating output file...'
        print*
    end if

    if ( proc_id == root_proc_id ) then
        ! call hdf_set_print_messages(exdm_inputs%control%verbose)
        call hdf_open_file(file_id, exdm_inputs%control%out_filename, status='NEW')
        call hdf_create_group(file_id, 'timing')
        call hdf_close_file(file_id)

        ! save inputs
        call exdm_inputs%save()
    end if

    if ( exdm_inputs%control%verbose ) then
        print*, 'Done creating output file!'
        print*
    end if

    ! initialize electronic configuration
    call exdm_elec_config%load(exdm_inputs)

    ! perform calculation
    select case ( trim(adjustl(exdm_inputs%control%calculation)) )

        case ( 'binned_scatter_rate' )

            call exdm_binned_scatter_rate(n_proc, proc_id, root_proc_id, &
                exdm_inputs, exdm_elec_config)

        case ( 'absorption_rate' )

            call exdm_absorption_rate(n_proc, proc_id, root_proc_id, &
                exdm_inputs, exdm_elec_config)

        case ( 'dielectric' )

            call exdm_dielectric(n_proc, proc_id, root_proc_id, &
                exdm_inputs, exdm_elec_config)

        case default

            print*, 'Please specify a calculation type in the [control] input.'
            print*
            
    end select

    if ( proc_id == root_proc_id ) then

        ! end timer
        call timer%end()

        call hdf_open_file(file_id, exdm_inputs%control%out_filename, &
            status='OLD', action='WRITE')

        ! save timing information
        call hdf_write_dataset(file_id, 'timing/dt_total', timer%dt)
        call hdf_write_dataset(file_id, 'timing/start_date', timer%start_date)
        call hdf_write_dataset(file_id, 'timing/end_date', timer%end_date)

        ! save version information
        call hdf_write_dataset(file_id, 'exdm_version', _cmake_version) 

        call hdf_close_file(file_id)

        ! shut down message
        call logger%shutdown_message(timer%pretty_dt_str())

    end if

    ! Finalize MPI
    call MPI_FINALIZE(mpi_err)

end program
