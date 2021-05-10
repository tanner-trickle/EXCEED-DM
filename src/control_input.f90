module control_input
    !! Collection of variables which specify how the program should run
    implicit none

    character(len=64) :: process = ''
    character(len=64) :: calc_mode = ''

    logical :: timer = .TRUE.
        !! Optional
        !!
        !! If .TRUE. the program will output timing information

    logical :: quiet = .FALSE.
        !! Don't print any output

    logical :: save_binned_rate_if = .FALSE.
        !! save the 2d differential rate data for 
        !! every i -> f transition
        !!
        !! memory intensive

    logical :: overwrite_output = .FALSE.
        !! if True the output file will be overwritten

    NAMELIST /control/ timer, &
                        quiet, &
                        process, &
                        calc_mode, &
                        save_binned_rate_if, &
                        overwrite_output

contains

    subroutine load_control(filename, verbose)
        !! Loads the control variables
        implicit none

        character(len=*) :: filename

        logical :: file_exists

        logical, optional :: verbose

        integer :: error

        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            open(100, file = trim(filename), iostat = error)

            read(100, nml=control, iostat=error)

            if ( error .ne. 0 ) then

                if ( verbose ) then

                    print*, '!! ERROR !!'
                    print*
                    print*, '    Problem with control namelist.'
                    print*
                    print*, '!!!!!!!!!!!'
                    print*

                end if

                stop

            end if

            close(100)

        else

            if ( verbose ) then

                print*, '!! ERROR !!'
                print*
                print*, '   Input file for control parameters : ', trim(filename), ' does NOT exist.'
                print*
                print*, '!!!!!!!!!!!'
                print*

            end if

            stop

        end if

    end subroutine

end module
