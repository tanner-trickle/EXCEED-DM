module info_messages
    !! Procedures for printing program information and errrors/warnings.

    implicit none

    integer :: sec_length = 80

contains

    subroutine exdm_shutdown_message(proc_id, root_process)
        !! Message to display when EXCEED-DM starts.

        use iso_fortran_env

        implicit none

        integer :: proc_id
        integer :: root_process
        integer :: n_proc
        character(len=64) :: version

        character(len=64) :: n_proc_str

        if ( proc_id == root_process ) then

            call print_section_seperator()
            print*
            print*, '    Ended at', trim(pretty_date_and_time())
            print*
            call print_section_seperator()
            print*

        end if

    end subroutine

    subroutine print_warning_message(message, verbose)
        !! Prints a warning message.

        implicit none

        character(len=*) :: message
        logical, optional :: verbose

        if ( verbose ) then

            call print_warning_msg_header()
            print*
            print*, '    ', message
            print*
            call print_warning_msg_seperator()
            print*

        end if

    end subroutine

    subroutine print_error_message(message, verbose)
        !! Prints an error message.

        implicit none

        character(len=*) :: message
        logical, optional :: verbose

        if ( verbose ) then

            call print_error_msg_header()
            print*
            print*, '    ', message
            print*
            call print_error_msg_seperator()
            print*

        end if

    end subroutine

    subroutine exdm_startup_message(proc_id, root_process, n_proc, version)
        !! Message to display when EXCEED-DM starts.

        use iso_fortran_env

        implicit none

        integer :: proc_id
        integer :: root_process
        integer :: n_proc
        character(len=64) :: version

        character(len=64) :: n_proc_str

        if ( proc_id == root_process ) then

            n_proc_str = int_to_str(n_proc)

            print*
            call print_section_seperator()
            print*
            print*, '    EXCEED-DM - v', trim(version)
            print*
            print*, '    Running on ', trim(adjustl(n_proc_str)), ' processors'
            print*, '    Compiled with ', compiler_version()
            print*
            print*, '    Started at', trim(pretty_date_and_time())
            print*
            call print_section_seperator()
            print*

        end if

    end subroutine

    subroutine print_section_seperator()
        !! Prints the section seperator.

        implicit none

        print*, repeat('-', sec_length)

    end subroutine

    subroutine print_error_msg_header()
        !! Prints the error message header.

        implicit none

        print*, repeat('!', 3)//' ERROR '//repeat('!', sec_length - 10)

    end subroutine

    subroutine print_error_msg_seperator()
        !! Prints the error message seperator.

        implicit none

        print*, repeat('!', sec_length)

    end subroutine

    subroutine print_warning_msg_header()
        !! Prints the warning message header.

        implicit none

        print*, repeat('~', 3)//' WARNING '//repeat('~', sec_length - 12)

    end subroutine

    subroutine print_warning_msg_seperator()
        !! Prints the warning message seperator.

        implicit none

        print*, repeat('~', sec_length)

    end subroutine

    function int_to_str(num) result(string)
        !! Converts an integer to a string.

        implicit none

        integer :: num
        character(len=64) :: string

        write(string, *) num

    end function

    function pretty_date_and_time() result( dt_str )
        !! Returns a nicely formatted string of the current date and time.

        implicit none

        integer :: dt_values(8)

        character(len=512) :: dt_str

        character(len=30) :: hr_str
        character(len=30) :: min_str
        character(len=30) :: s_str
        character(len=30) :: ms_str

        character(len=30) :: month_str 
        character(len=30) :: day_str 
        character(len=30) :: year_str 

        call date_and_time(values=dt_values)

        write(year_str, *) dt_values(1)
        write(month_str, *) dt_values(2)
        write(day_str, *) dt_values(3)

        ! ! append zeros for pretty printing
        if ( dt_values(5) < 10 ) then

            write(hr_str, *) dt_values(5)
            write(hr_str, *) '0'//trim(adjustl(hr_str))

        else

            write(hr_str, *) dt_values(5)

        end if

        if ( dt_values(6) < 10 ) then

            write(min_str, *) dt_values(6)
            write(min_str, *) '0'//trim(adjustl(min_str))

        else

            write(min_str, *) dt_values(6)

        end if

        if ( dt_values(7) < 10 ) then

            write(s_str, *) dt_values(7)
            write(s_str, *) '0'//trim(adjustl(s_str))

        else

            write(s_str, *) dt_values(7)

        end if

        if ( dt_values(8) < 10 ) then

            write(ms_str, *) dt_values(8)
            write(ms_str, *) '00'//trim(adjustl(ms_str))

        end if

        if ( ( dt_values(8) >= 10 ) .and. ( dt_values(8) < 100 ) ) then

            write(ms_str, *) dt_values(8)
            write(ms_str, *) '0'//trim(adjustl(ms_str))

        end if

        if ( dt_values(8) >= 100 ) then

            write(ms_str, *) dt_values(8)

        end if

        write(dt_str, *) trim(adjustl(hr_str)),&
                         ':', trim(adjustl(min_str)),&
                         ':', trim(adjustl(s_str)),&
                         '.', trim(adjustl(ms_str)), &
                         ' ', trim(adjustl(month_str)),&
                         '/', trim(adjustl(day_str)),&
                         '/', trim(adjustl(year_str))

    end function

end module
