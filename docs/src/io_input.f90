module io_input
    !! Input and output filenames.

    use info_messages

    implicit none

    type io_files_t
        !! All input and output filenames.
        character(len=512) :: nml_input_filename = ''
            !! Namelist input filename
        character(len=512) :: PW_data_filename = ''
            !! Bloch wave function coefficient file name
        character(len=64) :: run_description = ''
            !! description of the calculation
        character(len=512) :: out_folder = '.'
            !! Ouput folder
        character(len=512) :: out_filename = ''
            !! Output filename
            !!
            !! Setting this variable overrides the settings of out_folder
            !! and run_description

        character(len = 512) :: sto_data_filename = ''
            !! Input file specifying the Slater type orbital (sto)
            !! wave function coefficients for the core electron
            !! states

        character(len = 512) :: core_elec_config_filename = ''
            !! File specifying the core electron configuration

        character(len = 512) :: dielectric_filename = ''
            !! File specifying the dielectric function. If the dielectric
            !! will be computed this will be where the computed values are
            !! stored. If this file already exists, the dielectric function
            !! will be loaded from this file.

        contains

            procedure :: print => print_io
            procedure :: load => load_io_nml
            procedure :: save => save_io

    end type

contains

    subroutine save_io(self, filename, verbose)
        !! Saves `io_files`.
        use hdf5
        use h5lt

        use info_messages
        
        implicit none

        class(io_files_t) :: self
        character(len=*) :: filename
        logical, optional :: verbose

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id
        logical :: file_exists
        integer(HSIZE_T) :: dims1(1) = [1]
        integer :: error

        if ( verbose ) then
            print*, 'Saving I/O parameters...'
            print*
        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'io', group_id, error)

            ! write data

            call h5ltmake_dataset_string_f(file_id, &
                'io/nml_input_filename',                  &
                trim(self%nml_input_filename),                       &
                error)
            call h5ltmake_dataset_string_f(file_id, &
                'io/PW_data_filename',                  &
                trim(self%PW_data_filename),                       &
                error)
            call h5ltmake_dataset_string_f(file_id, &
                'io/sto_data_filename',                  &
                trim(self%sto_data_filename),                       &
                error)
            call h5ltmake_dataset_string_f(file_id, &
                'io/core_elec_config_filename',                  &
                trim(self%core_elec_config_filename),                       &
                error)
            call h5ltmake_dataset_string_f(file_id, &
                'io/dielectric_filename',                  &
                trim(self%dielectric_filename),                       &
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

    subroutine print_io(self, verbose)
        !! Prints `io_files` components.

        implicit none

        class(io_files_t) :: self
        logical, optional :: verbose

        if ( verbose ) then

            call print_section_seperator()
            print*, '    ---'
            print*, '    I/O'
            print*, '    ---'
            print*
            print*, '        PW data file         : ', trim(self%PW_data_filename)
            print*, '        STO data file        : ', trim(self%sto_data_filename)
            print*, '        Core electron config : ', trim(self%core_elec_config_filename)
            print*, '        Dielectric           : ', trim(self%dielectric_filename)
            print*
            call print_section_seperator()
            print*

        end if

    end subroutine

    subroutine load_io_nml(self, filename, verbose)
        !! Loads `io` parameters from a namelist.
        implicit none

        class(io_files_t) :: self
        character(len=*) :: filename
        logical, optional :: verbose

        logical :: file_exists
        integer :: error

        ! namelist
        character(len=512) :: PW_data_filename = ''
            !! DFT calculations input file
        character(len=64) :: run_description = ''
            !! description of the calculation
        character(len=512) :: out_folder = '.'
            !! Ouput folder
        character(len=512) :: out_filename = ''
            !! Output filename
            !!
            !! Setting this variable overrides the settings of out_folder
            !! and run_description

        character(len = 512) :: sto_data_filename = ''
            !! Input file specifying the Slater type orbital (sto)
            !! wave function coefficients for the core electron
            !! states

        character(len = 512) :: core_elec_config_filename = ''
            !! File specifying the core electron configuration

        character(len = 512) :: dielectric_filename = ''
            !! File specifying the dielectric function. If the dielectric
            !! will be computed this will be where the computed values are
            !! stored. If this file already exists, the dielectric function
            !! will be loaded from this file.

        NAMELIST /io/ PW_data_filename,          &
                        run_description,           &
                        out_folder,                &
                        out_filename,              &
                        sto_data_filename,           &
                        core_elec_config_filename, &
                        dielectric_filename

        if ( verbose ) then
            print*, 'Loading IO parameters...'
            print*
        end if

        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            open(100, file = trim(filename), iostat = error)
            read(100, nml=io, iostat=error)
            close(100)

            if ( error /= 0 ) then
                call print_error_message(&
                    'Problem reading IO namelist.', &
                    verbose = verbose)
                stop
            end if

            self%nml_input_filename        = filename
            self%PW_data_filename          = PW_data_filename
            self%sto_data_filename         = sto_data_filename
            self%core_elec_config_filename = core_elec_config_filename
            self%dielectric_filename       = dielectric_filename

            if ( trim(out_filename) .eq. '') then
                if ( trim(run_description) .eq. '' ) then
                    out_filename = trim(out_folder)//'/EXDMout.hdf5'
                else
                    out_filename = trim(out_folder)//'/EXDMout_'//trim(run_description)//'.hdf5'
                end if
            end if

            self%run_description = run_description
            self%out_folder = out_folder
            self%out_filename = out_filename

            call self%print(verbose = verbose)

        else

            call print_error_message(&
                'Input file for IO parameters : '//trim(filename)//' does NOT exist.',&
                verbose = verbose)
            stop

        end if

    end subroutine

    recursive subroutine create_output_file(io_files, overwrite_output, verbose)
        !! Creates the output file, if the file already exists, tries to create a new file with a random number
        !! attached recursively until a number is found.

        use hdf5
        use h5lt

        use prec

        implicit none

        type(io_files_t) :: io_files

        character(len=512) :: out_filename
        character(len=512) :: new_filename

        logical, optional :: verbose

        integer(HID_T) :: file_id

        logical :: file_exists
        integer :: error

        logical :: overwrite_output

        real(dp) :: r
        character(len=64) :: rand_num_str

        if ( verbose ) then

            print*, 'Creating output file...'
            print*

        end if

        out_filename = io_files%out_filename

        inquire(file = out_filename, exist = file_exists)

        if ( ( .not. file_exists ) .or. ( overwrite_output ) ) then

            call h5open_f(error)
            call h5fcreate_f(trim(out_filename), H5F_ACC_TRUNC_F, file_id, error)

            if ( error /= 0 ) then

                call print_error_message(&
                    'Could not create output file : '//trim(out_filename),&
                    verbose = verbose)
                stop

            end if

            call h5fclose_f(file_id, error)
            call h5close_f(error)

        else

            call random_number(r)
            write(rand_num_str, *) int(10**3*r)
            new_filename = './EXDMout_'//trim(adjustl(rand_num_str))//'.hdf5'

            call print_warning_message(&
                'Output file : '//trim(out_filename)//&
                ' already exists. Attempting to set output filename to :'//&
                trim(new_filename),&
                verbose = verbose)

            io_files%out_filename = new_filename

            call create_output_file(io_files, .FALSE., verbose = verbose)

        end if

    end subroutine

end module
