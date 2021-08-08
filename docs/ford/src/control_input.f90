module control_input
    !! Collection of variables controlling how `EXCEED-DM` is run, i.e. what is computed and how it is computed.

    use info_messages

    implicit none

    type control_t
        !! Stores all variables associated with the program control.

        character(len=64) :: process = ''
        character(len=64) :: calc_mode = ''

        logical :: timer = .TRUE.
            !! Optional
            !!
            !! If .TRUE. the program will output timing information

        logical :: quiet = .FALSE.
            !! Don't print any output

        logical :: overwrite_output = .FALSE.
            !! If .TRUE. the output file will be overwritten

        contains

            procedure :: print => print_control
            procedure :: load => load_control_nml
            procedure :: save => save_control

    end type

contains

    subroutine save_control(self, filename, verbose)
        !! Saves `control`.
        use hdf5
        use h5lt

        use info_messages
        
        implicit none

        class(control_t) :: self
        character(len=*) :: filename
        logical, optional :: verbose

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id
        logical :: file_exists
        integer(HSIZE_T) :: dims1(1) = [1]
        integer :: error

        if ( verbose ) then
            print*, 'Saving control parameters...'
            print*
        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'control', group_id, error)

            ! write data

            call h5ltmake_dataset_string_f(file_id, &
                'control/process',                  &
                trim(self%process),                       &
                error)
            call h5ltmake_dataset_string_f(file_id, &
                'control/calc_mode',                  &
                trim(self%calc_mode),                       &
                error)

            call h5fclose_f(file_id, error)
            call h5close_f(error)

        else

            call print_error_message(&
                'Output file : '//trim(filename)//' does NOT exist.', &
                verbose = verbose)
            stop

        end if

    end subroutine

    subroutine print_control(self, verbose)
        !! Prints `control` components.
        implicit none

        class(control_t) :: self

        logical, optional :: verbose

        if ( verbose ) then
            call print_section_seperator()
            print*, '    -------'
            print*, '    Control'
            print*, '    -------'
            print*
            print*, '        Physics process  : ', trim(self%process)
            print*, '        Calculation mode : ', trim(self%calc_mode)
            print*, '        Timing?          : ', self%timer
            print*
            call print_section_seperator()
            print*
        end if

    end subroutine

    subroutine load_control_nml(self, filename, verbose)
        !! Loads `control` parameters from a namelist.
        implicit none

        class(control_t) :: self
        character(len=*) :: filename
        logical, optional :: verbose

        logical :: file_exists
        integer :: error

        ! namelist
        logical :: timer = .TRUE.
        logical :: quiet = .FALSE.
        character(len=64) :: process = ''
        character(len=64) :: calc_mode = ''
        logical :: overwrite_output = .FALSE.

        NAMELIST /control/ timer,                &
                            quiet,               &
                            process,             &
                            calc_mode,           &
                            overwrite_output

        if ( verbose ) then
            print*, 'Loading control parameters...'
            print*
        end if

        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            open(100, file=trim(filename), iostat=error)
            read(100, nml=control, iostat=error)
            close(100)

            if ( error /= 0 ) then
                call print_error_message(&
                    'Problem reading control namelist.',&
                    verbose = verbose)
                stop
            end if

            self%timer            = timer
            self%quiet            = quiet
            self%process          = process
            self%calc_mode        = calc_mode
            self%overwrite_output = overwrite_output

            call self%print(verbose = verbose)

        else
            call print_error_message(&
                'Input file for control parameters : '//trim(filename)//' does NOT exist.', &
                verbose = verbose)
            stop
        end if

    end subroutine

end module
