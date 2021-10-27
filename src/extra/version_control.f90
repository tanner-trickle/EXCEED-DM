module version_control
    !! Keeps track of the current version of the program.
    implicit none

    character(len=64) :: version = "0.2.7"

contains

    subroutine save_version(filename, verbose)
        !! Saves the version number.

        use hdf5
        use h5lt

        use info_messages

        implicit none

        character(len=*) :: filename
        logical, optional :: verbose

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id
        logical :: file_exists
        integer :: error

        if ( verbose ) then
            print*, 'Saving version parameters...'
            print*
        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'version', group_id, error)

            ! write data
            call h5ltmake_dataset_string_f(file_id, &
                'version/number', &
                version, &
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

end module
