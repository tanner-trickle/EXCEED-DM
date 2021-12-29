module timer_util
    !! Defines the 'timer' type to help with timing the program.

    use prec, only: dp
    use mpi

    implicit none

    ! keeps track of all the named timers
    character(len=256) :: timer_id_list(100)

    type timer_t

        real(dp) :: start_time
        real(dp) :: end_time
        real(dp) :: dt

        ! year/month/day
        integer :: start_date(3)
        integer :: end_date(3)

        contains

            procedure :: start => start_timer
            procedure :: end => end_timer
            procedure :: save => save_timer
            ! procedure :: print => print_timer

    end type

contains

    subroutine start_timer(self)

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

    subroutine end_timer(self)

        implicit none

        class(timer_t) :: self

        integer :: dt_values(8)

        self%end_time = MPI_Wtime()
        self%dt = self%end_time - self%start_time

        call date_and_time(values=dt_values)

        self%end_date = dt_values(:3)

    end subroutine

    subroutine save_timer(self, filename, name, verbose)

        use hdf5
        use h5lt

        use info_messages

        implicit none

        class(timer_t) :: self
        character(len=*) :: filename
        character(len=*) :: name
        logical, optional :: verbose

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id
        logical :: file_exists

        logical :: group_exists

        integer(HSIZE_T) :: dims1(1) = [1]

        integer :: error

        character(len=256) :: timing_group

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5lexists_f(file_id, 'timing', group_exists, error)

            ! create the entry if it doesn't exist
            if ( .not. group_exists ) then

                call h5gcreate_f(file_id, 'timing', group_id, error)

            end if

            timing_group = 'timing/'//trim(adjustl(name))

            call h5gcreate_f(file_id, trim(adjustl(timing_group)) , group_id, error)

            dims1 = [3]
            call h5ltmake_dataset_int_f(file_id, &
                trim(adjustl(timing_group))//'/start_date', &
                size(dims1), dims1,&
                self%start_date, &
                error)
            call h5ltmake_dataset_int_f(file_id, &
                trim(adjustl(timing_group))//'/end_date', &
                size(dims1), dims1,&
                self%end_date, &
                error)
            dims1 = [1]
            call h5ltmake_dataset_double_f(file_id, &
                trim(adjustl(timing_group))//'/dt',&
                size(dims1), dims1,&
                self%dt, error)
            call h5ltmake_dataset_string_f(file_id,&
                trim(adjustl(timing_group))//'/pretty_dt',&
                trim(adjustl(pretty_time_format(self%dt))), error)

            call h5fclose_f(file_id, error)
            call h5close_f(error)

        else

            call print_error_message('Output file : '//trim(filename)//' does NOT exist.')
            stop

        end if

    end subroutine

    function pretty_time_format(t) result(time_str)
        !! Returns a nicely formatted string of the time.

        implicit none

        real(dp) :: t

        integer :: time_hr
        integer :: time_min
        integer :: time_sec

        character(len=30) :: hr_str
        character(len=30) :: min_str
        character(len=30) :: s_str
        character(len=30) :: ms_str

        character(len = 512) :: time_str

        if ( t < 1.0_dp ) then

            write(ms_str, *) t*10**3
            write(time_str, *) trim(adjustl(ms_str))//' ms'

        else

            time_hr = floor(t/3600.0_dp)
            time_min = floor((t - time_hr*3600.0_dp)/60.0_dp)
            time_sec = nint(t - time_hr*3600.0_dp - time_min*60.0_dp)

            write(hr_str, *) time_hr
            write(min_str, *) time_min
            write(s_str, *) time_sec

            write(time_str, *) trim(adjustl(hr_str))//' hr '//&
                               trim(adjustl(min_str))//' min '//&
                               trim(adjustl(s_str))//' s'

        end if

    end function

end module
