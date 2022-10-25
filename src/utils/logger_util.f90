module logger_util
    ! A utility for logging messages to the user through a new type, :code:`logger_t`
    ! 
    ! .. note:: No external dependancies.

    implicit none

    type :: logger_t
        ! Collection of variables and procedures for easily printing nice log messages to the user.
        integer :: section_length = 80
            ! Length of a section header/footer

        contains

            procedure :: startup_message => logger_util_startup_message
            procedure :: shutdown_message => logger_util_shutdown_message

    end type

contains

    subroutine logger_util_shutdown_message(self, exec_time_str)
        ! Prints a shutdown message about the program.

        implicit none

        class(logger_t) :: self

        character(len=*), intent(in) :: exec_time_str

        print*, repeat('-', self%section_length)
        print*
        print*, '    Run time: '//trim(adjustl(exec_time_str))
        print*
        print*, '    Ended at '//trim(adjustl(pretty_date_and_time()))
        print*
        print*, repeat('-', self%section_length)
        print*

    end subroutine

    subroutine logger_util_startup_message(self, prog_name, version, n_proc)
        ! Prints a startup message about the program.

        use iso_fortran_env

        implicit none

        class(logger_t) :: self

        character(len=*), intent(in) :: prog_name
            ! Program name
        character(len=*), intent(in) :: version
            ! Program version

        integer, intent(in) :: n_proc
            ! Number of processors
        character(len=64) :: n_proc_str
            ! Number of processors, string

        write(n_proc_str, *) n_proc

        print*
        print*, repeat('-', self%section_length)
        print*
        print*, '    '//trim(prog_name)//' - v'//trim(version)
        print*
        print*, '    Running on '//trim(adjustl(n_proc_str))//' processors'
        print*, '    Compiled with '//compiler_version()
        print*
        print*, '    Started at '//trim(adjustl(pretty_date_and_time()))
        print*
        print*, repeat('-', self%section_length)
        print*

    end subroutine

    function pretty_date_and_time() result( dt_str )
        ! Returns a nicely formatted string of the current date and time.

        implicit none

        integer :: dt_values(8)

        character(len=:), allocatable :: dt_str

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

        dt_str = trim(adjustl(hr_str))//&
                         ':'//trim(adjustl(min_str))//&
                         ':'//trim(adjustl(s_str))//&
                         '.'//trim(adjustl(ms_str))// &
                         ' '//trim(adjustl(month_str))//&
                         '/'//trim(adjustl(day_str))//&
                         '/'//trim(adjustl(year_str))

    end function

end module
