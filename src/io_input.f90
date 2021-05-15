module io_input
    !! Input and output file names
    implicit none

    character(len=512) :: nml_input_filename = ''
        !! Namelist input filename, 
    character(len=512) :: DFT_input_filename = ''
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

    character(len = 512) :: sto_wf_filename = ''
        !! Input file specifying the Slater type orbital (sto)
        !! wave function coefficients for the core electron
        !! states

    character(len = 512) :: core_elec_config_filename = ''
        !! File specifying the core electron configuration

    NAMELIST /io/ DFT_input_filename, &
                    run_description, &
                    out_folder, &
                    out_filename, &
                    sto_wf_filename, &
                    core_elec_config_filename
contains

    subroutine print_io(verbose)

        implicit none

        logical, optional :: verbose

        if ( verbose ) then

            print*, '----------------------------------------'
            print*, '    ---'
            print*, '    I/O'
            print*, '    ---'
            print*
            print*, '        DFT input            : ', trim(DFT_input_filename)
            print*, '        STO coefficients     : ', trim(sto_wf_filename)
            print*, '        Core electron config : ', trim(core_elec_config_filename)
            print*

        end if

    end subroutine

    subroutine load_io(filename, verbose)
        !! Loads the io namelist
        implicit none

        character(len=*) :: filename

        logical :: file_exists

        logical, optional :: verbose

        integer :: error

        if ( verbose ) then
            print*, 'Loading IO parameters...'
            print*
        end if

        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            open(100, file = trim(filename), iostat = error)
            read(100, nml=io, iostat=error)
            close(100)

            if ( error .ne. 0 ) then

                if ( verbose ) then

                    print*, '!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                    print*
                    print*, '    Problem reading IO namelist.'
                    print*
                    print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                    print*

                end if

                stop

            end if

            if ( trim(out_filename) .eq. '') then
                if ( trim(run_description) .eq. '' ) then
                    out_filename = trim(out_folder)//'/EXDMout.hdf5'
                else
                    out_filename = trim(out_folder)//'/EXDMout_'//trim(run_description)//'.hdf5'
                end if
            end if

            call print_io(verbose = verbose)

            if ( verbose ) then

                print*, '----------------------------------------'
                print*

            end if

        else

            if ( verbose ) then

                print*, '!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*
                print*, '    Input file for IO parameters : ', trim(filename), ' does NOT exist.'
                print*
                print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*

            end if

            stop

        end if

    end subroutine

    recursive subroutine create_output_file(filename, overwrite_output, verbose)
        !! Creates the output file

        use hdf5
        use h5lt

        use prec

        implicit none

        character(len=*) :: filename
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

        inquire(file = filename, exist = file_exists)

        if ( ( .not. file_exists ) .or. ( overwrite_output ) ) then

            call h5open_f(error)
            call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, error)

            if ( (error .ne. 0) .and. verbose ) then

                print*, '!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*
                print*, '   Could not create output file : ', trim(filename)
                print*
                print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*

                stop

            end if

            call h5fclose_f(file_id, error)
            call h5close_f(error)

        else

            if ( verbose ) then

                print*, '~~~ WARNING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
                print*
                print*, '   Output file : ', trim(filename), ' already exists.'
                print*

            end if

            call random_number(r)
            write(rand_num_str, *) int(10**3*r)
            out_filename = './EXDMout_'//trim(adjustl(rand_num_str))//'.hdf5'

            if ( verbose ) then

                print*, '   Attempting to set output filename to : ', trim(out_filename)
                print*
                print*, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
                print*

            end if

            call create_output_file(trim(out_filename), .FALSE., verbose = verbose)

        end if

    end subroutine

end module
