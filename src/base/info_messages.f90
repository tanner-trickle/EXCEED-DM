module info_messages
    !! Procedures for printing program information and errrors/warnings.

    implicit none

contains

    subroutine print_warning_message(message, verbose)
        !! Prints a warning message.

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
            call print_section_seperator()
            print*

        end if

    end subroutine

    subroutine print_section_seperator()
        !! Prints the section seperator.

        implicit none

        print*, '----------------------------------------------------------------------'

    end subroutine

    subroutine print_error_msg_header()
        !! Prints the error message header.

        implicit none

        print*, '!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

    end subroutine

    subroutine print_error_msg_seperator()
        !! Prints the error message seperator.

        implicit none

        print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

    end subroutine

    subroutine print_warning_msg_header()
        !! Prints the warning message header.

        implicit none

        print*, '~~~ WARNING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

    end subroutine

    subroutine print_warning_msg_seperator()
        !! Prints the warning message seperator.

        implicit none

        print*, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

    end subroutine

    function int_to_str(num) result(string)
        !! Converts an integer to a string.

        implicit none

        integer :: num
        character(len=64) :: string

        write(string, *) num

    end function

end module
