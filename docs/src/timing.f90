module timing
    !! Useful functions for timing a program
    use prec

    implicit none

    real(dp) :: time(100)
        !! holds raw timing variables
    real(dp) :: delta_t(100)
        !! holds difference in timing variables

contains

    function pretty_time_format(t) result(time_str)
        !! Fills time_str with a nicely formatted version of time 

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
