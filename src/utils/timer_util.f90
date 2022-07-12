module timer_util
    ! A utility for timing an MPI program.
    !
    ! ..note: Depends on MPI.

    use iso_fortran_env
    use mpi

    implicit none

    type :: timer_t
        ! Collection of variables and procedures for easily keeping track of a programs run time

        real(real64) :: start_time
        real(real64) :: end_time
        real(real64) :: dt
            ! Units : s

        integer :: start_date(3)
            ! year, month, day
        integer :: end_date(3)
            ! year, month, day

        contains

            procedure :: start => timer_util_start
            procedure :: end => timer_util_end
            procedure :: pretty_dt_str => timer_util_pretty_dt_str

    end type

contains

    subroutine timer_util_start(self)
        ! Start the timer.

        implicit none

        class(timer_t) :: self

        integer :: dt_values(8)

        self%start_time = MPI_Wtime()
        self%end_time = MPI_Wtime()
        self%dt = self%end_time - self%start_time

        call date_and_time(values=dt_values)

        self%start_date = dt_values(:3)
        self%end_date = dt_values(:3)

    end subroutine

    subroutine timer_util_end(self)
        ! End the timer.

        implicit none

        class(timer_t) :: self

        integer :: dt_values(8)

        self%end_time = MPI_Wtime()
        self%dt = self%end_time - self%start_time

        call date_and_time(values=dt_values)

        self%end_date = dt_values(:3)

    end subroutine

    function timer_util_pretty_dt_str(self) result(time_str)
        ! Returns a nicely formatted string of the time.

        implicit none

        class(timer_t) :: self

        integer :: time_hr
        integer :: time_min
        integer :: time_sec

        character(len=30) :: hr_str
        character(len=30) :: min_str
        character(len=30) :: s_str
        character(len=30) :: ms_str

        character(len = :), allocatable :: time_str

        if ( self%dt < 1.0_real64 ) then

            write(ms_str, *) self%dt*10**3
            time_str = trim(adjustl(ms_str))//' ms'

        else

            time_hr = floor(self%dt/3600.0_real64)
            time_min = floor((self%dt - time_hr*3600.0_real64)/60.0_real64)
            time_sec = nint(self%dt - time_hr*3600.0_real64 - time_min*60.0_real64)

            write(hr_str, *) time_hr
            write(min_str, *) time_min
            write(s_str, *) time_sec

            time_str = trim(adjustl(hr_str))//' hr '//&
                               trim(adjustl(min_str))//' min '//&
                               trim(adjustl(s_str))//' s'

        end if

    end function

end module
