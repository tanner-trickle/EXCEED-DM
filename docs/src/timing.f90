module timing
    !! Useful functions for timing a program.
    use prec

    implicit none

    real(dp) :: time(100)
        !! holds raw timing variables
    real(dp) :: delta_t(100)
        !! holds difference in timing variables

contains

    subroutine print_timing_info(dt, verbose)
        !! Prints a program timing info.

        use info_messages

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

    function pretty_time_format(t) result(time_str)
        !! Returns a nicely formatted string of the time.

        implicit none

        real(dp) :: t

        integer :: time_hr
        integer :: time_min
        integer :: time_sec

        character(len = 512) :: time_str

        if ( t .lt. 1.0_dp ) then

            write(time_str, *) t*10**3, 'ms'

        else

            time_hr = floor(t/3600.0_dp)
            time_min = floor((t - time_hr*3600.0_dp)/60.0_dp)
            time_sec = nint(t - time_hr*3600.0_dp - time_min*60.0_dp)

            write(time_str, *) time_hr, 'hr ', time_min, 'min', time_sec, 's'

        end if

    end function

end module
