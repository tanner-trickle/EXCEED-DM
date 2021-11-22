module info_messages
    !! Procedures for printing program information and errrors/warnings.

    implicit none

    integer :: sec_length = 80

contains

    subroutine exdm_shutdown_message(proc_id, root_process)
        !! Message to display when EXCEED-DM starts.

        use iso_fortran_env
        use timing

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

    subroutine print_timing_info(dt, verbose)
        !! Prints a program timing info.

        use timing

        implicit none

        real(dp) :: dt
        logical, optional :: verbose

        if ( verbose ) then

            call print_section_seperator()
            print*, '    ------'
            print*, '    Timing'
            print*, '    ------'
            print*
            print*, '        Run time : '
            print*, '        ', trim(pretty_time_format(dt))
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
        use timing

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

end module
