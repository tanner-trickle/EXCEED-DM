module timing
    !! Useful functions for timing a program.
    use prec

    implicit none

    real(dp) :: time(100)
        !! holds raw timing variables
    real(dp) :: delta_t(100)
        !! holds difference in timing variables

contains

    subroutine save_timing_info(filename, dt, verbose)
        !! Saves the version number.

        use hdf5
        use h5lt

        implicit none

        character(len=*) :: filename
        logical, optional :: verbose

        real(dp) :: dt

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id
        integer(HSIZE_T) :: dims1(1) = [1]
        logical :: file_exists
        integer :: error

        if ( verbose ) then
            print*, 'Saving timing parameters...'
            print*
        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'timing', group_id, error)

            ! write data
            call h5ltmake_dataset_string_f(file_id, &
                'timing/total_run_time_str', &
                trim(adjustl(pretty_time_format(dt))), &
                error)

            call h5ltmake_dataset_string_f(file_id, &
                'timing/start_date_str', &
                trim(adjustl(pretty_date_and_time())), &
                error)

            call h5ltmake_dataset_double_f(file_id, &
                'timing/total_run_time_s', &
                size(dims1), dims1, &
                dt, &
                error)

            call h5fclose_f(file_id, error)
            call h5close_f(error)

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
