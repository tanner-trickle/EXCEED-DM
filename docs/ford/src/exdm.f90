program exdm
    !! EXCEED-DM : EXtended Calculation of Electronic Excitations for Direct detection of Dark Matter
    !!
    !! Github    : https://github.com/tanner-trickle/EXCEED-DM 

    use iso_fortran_env
    use mpi 

    use info_messages
    use timing

    use version_control

    use control_input
    use io_input

    use material_type
    use dm_model_type
    use expt_type

    use exdm_absorption
    use exdm_scatter
    use exdm_dielectric

    implicit none

    integer :: proc_id
        !! Processor ID
    integer :: n_proc
        !! Number of processors
    integer :: root_process = 0
        !! Root processor ID 
    integer :: err
        !! Error code
    character(len=512) :: nml_input_filename = ''
        !! Namelist input filename.

    logical :: verbose = .FALSE.
        !! If verbose = .TRUE., print output

    type(control_t) :: main_control
        !! Control parameters
    type(io_files_t) :: io_files
        !! Input/output filenames
    type(material_t) :: target_mat
        !! The target material
    type(dm_model_t) :: dm_model
        !! Dark matter model parameters
    type(expt_t) :: expt
        !! Experimental parameters

    ! Initialize MPI
    call MPI_INIT(err)
    call MPI_COMM_RANK(MPI_COMM_WORLD, proc_id, err)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, n_proc, err)

    call exdm_startup_message(proc_id, root_process, n_proc, version)

    call get_command_argument(1, nml_input_filename)

    ! start timing the program
    if ( proc_id == root_process ) then 
        time(1) = MPI_Wtime()
    end if

    ! load inputs
    if ( proc_id == root_process ) then

        call main_control%load(nml_input_filename, verbose = .TRUE.)

        if ( main_control%quiet ) then
            verbose = .FALSE.
        else
            verbose = .TRUE.
        end if

    else

        call main_control%load(nml_input_filename, verbose = .FALSE.)

    end if

    call io_files%load(nml_input_filename, verbose = verbose)
    call target_mat%load(io_files%nml_input_filename, verbose = verbose)
    call dm_model%load(io_files%nml_input_filename, verbose = verbose)
    call expt%load(io_files%nml_input_filename, verbose = verbose)

    ! create the output file
    if ( proc_id == root_process ) then
        call create_output_file(&
            io_files,& 
            main_control%overwrite_output,&
            verbose = verbose)
    end if

    ! compute and save data
    if ( trim(main_control%process) == 'scatter' ) then

        call run_exdm_scatter(proc_id, root_process, n_proc, &
            io_files, main_control, target_mat, expt, &
            dm_model, verbose = verbose)

    else if ( trim(main_control%process) == 'absorption' ) then

        call run_exdm_absorption(proc_id, root_process, n_proc, &
            io_files, main_control, target_mat, expt, &
            dm_model, verbose = verbose)

    else if ( trim(main_control%process) ==  'dielectric' ) then

        call run_exdm_dielectric(proc_id, root_process, n_proc, &
            io_files, main_control, target_mat, verbose = verbose)

    else

        call print_error_message(&
            'Process : '//trim(main_control%process)//' is not implemented.', &
            verbose = verbose)
        stop

    end if

    ! save input data common to all processes
    if ( proc_id == root_process ) then
        call target_mat%save(io_files%out_filename, verbose = verbose)

        call io_files%save(io_files%out_filename, verbose = verbose)
        call main_control%save(io_files%out_filename, verbose = verbose)

        if ( trim(main_control%process) /= 'dielectric' ) then

            call dm_model%save(io_files%out_filename, verbose = verbose)
            call expt%save(io_files%out_filename, verbose = verbose)

        end if
        call save_version(io_files%out_filename, verbose = verbose)
    end if

    ! Time program
    if ( proc_id == root_process ) then
        time(2) = MPI_Wtime()
        call print_timing_info(time(2) - time(1), verbose = verbose)
    end if

    call MPI_FINALIZE(err)

end program
