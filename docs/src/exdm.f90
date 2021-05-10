program exdm
    !! EXCEED-DM : EXtended Calculation of Electronic Excitations for Direct
    !! detection of Dark Matter

    use iso_fortran_env
    use mpi 

    use version_control
    use control_input
    use timing
    use io_input
    use material_input

    use dme_scatter

    implicit none

    integer :: proc_id
        !! Open MPI, processor ID
    integer :: n_proc
        !! Open MPI, number of processors
    integer :: root_process = 0
        !! Open MPI, root processor ID 
    integer :: err
        !! Open MPI error code

    logical :: verbose = .FALSE.
        !! If verbose = .TRUE., print output

    call MPI_INIT(err)
    call MPI_COMM_RANK(MPI_COMM_WORLD, proc_id, err)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, n_proc, err)

    if (proc_id .eq. root_process) then 

        print*
        print*, '--------------------'
        print*
        print*, '   EXCEED-DM - v', version
        print*
        print*, '--------------------'
        print*
        print*, 'Running on ', n_proc, 'processors'
        print*, 'Compiled with ', compiler_version()
        print*
        print*, '----------'
        print*

        ! prints output information
        if ( .not. quiet ) then
            verbose = .TRUE.
        else
            verbose = .FALSE.
        end if

        time(1) = MPI_Wtime()

    end if

    call get_command_argument(1, nml_input_filename)

    ! load inputs
    call load_control(nml_input_filename, verbose=verbose)
    call load_io(nml_input_filename, verbose=verbose)
    call load_material(nml_input_filename, verbose=verbose)

    ! create the output file
    if ( proc_id == root_process ) then
        call create_output_file(out_filename, overwrite_output, verbose=verbose)
    end if

    ! compute and save data
    if ( trim(process) == 'scatter' ) then

        call run_dme_scatter(proc_id, root_process, out_filename, &
            nml_input_filename, DFT_input_filename, sto_wf_filename, &
            core_elec_config_filename, n_proc, save_binned_rate_if, verbose = verbose)

    else

        print*, '!!! ERROR !!!'
        print*
        print*, '    Process : ', trim(process), ' is not implemented.'
        print*

        stop

    end if

    ! save input data common to all processes
    if ( proc_id .eq. root_process ) then

        if ( verbose ) then
            print*, 'Saving input data...'
            print*
        end if

        call save_material(out_filename, verbose=verbose)

        if ( verbose ) then
            print*, '----------'
            print*
        end if

    end if

    ! Time program
    if ( proc_id .eq. root_process ) then

        time(2) = MPI_Wtime()

        if ( verbose ) then

            print*, 'Run time : '
            print*, trim(pretty_time_format(time(2) - time(1)))
            print*
            print*, '----------'
            print*

        end if

    end if

    call MPI_FINALIZE(err)

end program
