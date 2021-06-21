module absorption_input
    !! Load parameters relevant for absorption calculations

    use hdf5
    use h5lt

    use prec

    implicit none

    integer :: n_omega = 1
        !! Number of omega/mass points to compute for

    real(dp) :: log_omega_min = -2.0_dp
    real(dp) :: log_omega_max = 2.0_dp

    integer :: n_width_b = 1
    integer :: n_width_a = 1

    integer :: n_width_max = 0

    integer :: n_widths
        !! n_width_a*n_width_b
    
    real(dp) :: width_a_min = 0.0_dp
    real(dp) :: width_a_max = 0.0_dp

    real(dp) :: log_width_b_min = -6.0_dp
    real(dp) :: log_width_b_max = -1.0_dp

    real(dp) :: width_max_min = 100.0_dp
    real(dp) :: width_max_max = 100.0_dp

    real(dp) :: sigma_gamma = 3

    NAMELIST /absorption/ n_omega, &
                            log_omega_min, &
                            log_omega_max, &
                            n_width_a, &
                            width_a_min, &
                            width_a_max, &
                            n_width_b, &
                            log_width_b_min, &
                            log_width_b_max, &
                            n_width_max, &
                            width_max_min, &
                            width_max_max, &
                            sigma_gamma

    !!!

    real(dp), allocatable :: omega_list(:)
        !! Dim : [n_omega]
    real(dp), allocatable :: delta_list(:, :)
        !! Dim : [n_omega, n_widths]
    real(dp), allocatable :: width_info(:, :)
        !! Dim : [n_widths, 2]
        !!
        !! Width parameters, [a, b]

contains

    subroutine save_absorption_input(filename, verbose)

        implicit none
        character(len=*) :: filename

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id

        logical, optional :: verbose
        logical :: file_exists

        integer(HSIZE_T) :: dims1(1) = [1]
        integer(HSIZE_T) :: dims2(2)

        integer :: error

        if ( verbose ) then

            print*, 'Saving absorption input...'
            print*

        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'abs_input', group_id, error)

            ! ! write data
            dims1 = [n_omega]
            call h5ltmake_dataset_double_f(file_id, 'abs_input/omega_list', size(dims1), dims1,&
                omega_list, error)

            ! dims2 = [n_omega, n_widths]
            ! call h5ltmake_dataset_double_f(file_id, 'abs_input/delta_list', size(dims1), dims1,&
            !     delta_list, error)

            dims2 = [n_widths, 3]
            call h5ltmake_dataset_double_f(file_id, 'abs_input/width_info', size(dims2), dims2,&
                width_info, error)

            call h5fclose_f(file_id, error)
            call h5close_f(error)

            if ( verbose ) then
                print*, '----------'
                print*
            end if

        else

            if ( verbose ) then

                print*, '!! ERROR !!'
                print*
                print*, '   Output file : ', trim(filename), ' does NOT exist.'
                print*
                print*, '!!!!!!!!!!!'
                print*

            end if

            stop

        end if

    end subroutine

    subroutine print_absorption_input(verbose)
        implicit none

        logical, optional :: verbose

        if ( verbose ) then

            print*, 'Absorption parameters : '
            print*
            print*, '   Number of omega            = ', n_omega
            print*, '   Number of width parameters = ', n_widths
            print*
            print*, '   Minimum omega = ', minval(omega_list), 'eV'
            print*, '   Maximum omega = ', maxval(omega_list), 'eV'
            print*
            print*, '----------'
            print*

        end if
    end subroutine

    subroutine load_absorption_input(filename, verbose)
        !! Load the absorption inputs

        implicit none

        logical, optional :: verbose

        character(len=*) :: filename

        logical :: file_exists

        integer :: error

        integer :: i, a, b, j, c
        real(dp) :: log_omega

        real(dp) :: a_param, b_param, width, width_max_param

        if ( verbose ) then

            print*, 'Loading absorption parameters...'
            print*

        end if

        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            open(100, file = trim(filename), iostat = error)
            read(100, nml=absorption, iostat=error)
            close(100)

            n_widths = n_width_a*n_width_b*(n_width_max + 1)

            allocate(omega_list(n_omega))
            omega_list = 10.0_dp**log_omega_min

            allocate(delta_list(n_omega, n_widths))
            delta_list = 0.0_dp

            allocate(width_info(n_widths, 3))
            width_info = 0.0_dp

            do i = 1, n_omega

                log_omega = log_omega_min + (log_omega_max - log_omega_min)*&
                    (i - 1.0_dp)/max(1.0_dp, n_omega - 1.0_dp)

                omega_list(i) = 10.0_dp**log_omega

                j = 1

                do a = 1, n_width_a

                    a_param = width_a_min + (width_a_max - width_a_min)*(a - 1.0_dp)/&
                       max(1.0_dp, n_width_a - 1.0_dp) 

                    do b = 1, n_width_b

                        b_param = log_width_b_min + (log_width_b_max - log_width_b_min)*(b - 1.0_dp)/&
                           max(1.0_dp, n_width_b - 1.0_dp) 

                        do c = 1, n_width_max + 1

                            if ( c == n_width_max + 1 ) then

                                width_max_param = 1.0e3_dp

                            else

                                width_max_param = width_max_min + (width_max_max - width_max_min)*(c - 1.0_dp)/&
                                   max(1.0_dp, n_width_max - 1.0_dp) 

                            end if

                            if ( i == 1 ) then

                                width_info(j, 1) = a_param
                                width_info(j, 2) = 10.0_dp**b_param
                                width_info(j, 3) = width_max_param

                            end if

                            delta_list(i, j) = min( width_max_param, a_param + (10.0_dp**b_param)*omega_list(i) )

                            j = j + 1

                        end do

                    end do

                end do

            end do

            call print_absorption_input(verbose=verbose)

        else

            if ( verbose ) then

                print*, '!! ERROR !!'
                print*
                print*, '   Input file for absorption parameters : ', trim(filename), ' does NOT exist.'
                print*
                print*, '!!!!!!!!!!!'
                print*

            end if

            stop

        end if

    end subroutine


end module
