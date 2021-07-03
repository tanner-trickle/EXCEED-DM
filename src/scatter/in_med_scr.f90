module in_med_scr
    !! Handle analytic/numeric dielectric or other in medium effects which screen the interaction
    !! 
    !! rate ~ 1/screen^2

    use hdf5
    use h5lt

    use prec
    use constants
    use math_mod

    use dielectric_input
    use calc_dielectric

    implicit none

    character(len=64) :: screen_type = ''

    !!!
    !! Eq 6 from 
    !!      https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892
    real(dp) :: di_e0
    real(dp) :: di_q_tf
    real(dp) :: di_omega_p
    real(dp) :: di_alpha
    !!!

    logical :: include_screen = .TRUE.
        !! toggle whether screening effects are included or not

    !! For numerical calculation of dielectric
    complex(dp), allocatable :: screen_mat(:, :, :, :)
        !! Numerically computed screening matrix, binned in [ omega, q, q_theta, q_phi ]

    integer :: scr_n_w_bins
    integer :: scr_n_q_bins
    integer :: scr_n_q_theta_bins
    integer :: scr_n_q_phi_bins

    real(dp) :: scr_w_bin_width
    real(dp) :: scr_q_bin_width

    NAMELIST /in_medium/ screen_type,&
                            di_e0, &
                            di_q_tf, &
                            di_omega_p, &
                            di_alpha, &
                            include_screen

contains

    subroutine load_scr_dielectric(filename, verbose)

        implicit none

        character(len=*) :: filename

        logical, optional :: verbose

        integer(HID_T) :: file_id

        logical :: file_exists

        integer(HSIZE_T) :: dims1(1) = [1]
        integer(HSIZE_T) :: dims4(4)

        integer :: error

        real(dp), allocatable :: dielectric_buff(:, :, :, :)

        if ( verbose ) then

            print*, 'Loading dielectric from file... '
            print*

        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)

            call h5ltread_dataset_int_f(file_id, 'dielectric/n_w_bins',&
                scr_n_w_bins, dims1, error)
            call h5ltread_dataset_int_f(file_id, 'dielectric/n_q_bins',&
                scr_n_q_bins, dims1, error)
            call h5ltread_dataset_int_f(file_id, 'dielectric/n_q_theta_bins',&
                scr_n_q_theta_bins, dims1, error)
            call h5ltread_dataset_int_f(file_id, 'dielectric/n_q_phi_bins',&
                scr_n_q_phi_bins, dims1, error)

            call h5ltread_dataset_double_f(file_id, 'dielectric/w_bin_width',&
                scr_w_bin_width, dims1, error)
            call h5ltread_dataset_double_f(file_id, 'dielectric/q_bin_width',&
                scr_q_bin_width, dims1, error)

            dims4 = [scr_n_w_bins, scr_n_q_bins, scr_n_q_theta_bins, scr_n_q_phi_bins]

            allocate(dielectric_buff(scr_n_w_bins, scr_n_q_bins, scr_n_q_theta_bins, &
                scr_n_q_phi_bins))
            allocate(screen_mat(scr_n_w_bins, scr_n_q_bins, scr_n_q_theta_bins, &
                scr_n_q_phi_bins))

            screen_mat = (0.0_dp, 0.0_dp)

            call h5ltread_dataset_double_f(file_id, 'dielectric/dielectric_r',&
                dielectric_buff, dims4, error)

            screen_mat = screen_mat + dielectric_buff

            call h5ltread_dataset_double_f(file_id, 'dielectric/dielectric_c',&
                dielectric_buff, dims4, error)

            screen_mat = screen_mat + ii*dielectric_buff

            call h5fclose_f(file_id, error)
            call h5close_f(error)

            if ( verbose ) then
                print*, '----------------------------------------'
                print*
            end if

        else

            if ( verbose ) then

                print*, '!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*
                print*, '    Input file for dielectric : ', trim(filename), ' does NOT exist.'
                print*
                print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*

            end if

            stop

        end if

    end subroutine

    subroutine print_in_med(verbose)
        implicit none

        logical, optional :: verbose

        if ( verbose ) then

            print*, '----------------------------------------'
            print*, '    ---------'
            print*, '    In-Medium'
            print*, '    ---------'
            print*
            print*, '        Include screen? : ', include_screen
            print*, '        Screen type     : ', trim(screen_type)
            print*

            if ( trim(screen_type) == 'analytic' ) then

                print*, '        e_0      : ', di_e0
                print*, '        alpha    : ', di_alpha
                print*, '        q_TF     : ', di_q_tf/1.0e3_dp, 'keV'
                print*, '        omega_p  : ', di_omega_p, 'eV'
                print*

            end if

        end if

    end subroutine

    subroutine load_in_med_scr(filename, dielectric_filename, DFT_input_filename, &
           proc_id, root_process, n_proc, verbose)

        implicit none

        character(len=*) :: filename
        character(len=*) :: dielectric_filename
        character(len=*) :: DFT_input_filename

        integer :: proc_id, root_process, n_proc

        logical :: file_exists

        logical, optional :: verbose

        integer :: error, err

        if ( verbose ) then

            print*, 'Loading in-medium screening parameters...'
            print*

        end if

        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            open(100, file = trim(filename), iostat = error)
            read(100, nml=in_medium, iostat=error)
            close(100)

            call print_in_med(verbose = verbose)

            if ( verbose ) then

                print*, '----------------------------------------'
                print*

            end if

            if ( trim(screen_type) == 'numeric' ) then

                call load_dielectric_input(filename, verbose = verbose)

                if ( .not. load_dielectric_from_file ) then

                    ! compute and save the dielectric
                    call run_dielectric_calc(filename, dielectric_filename, &
                        DFT_input_filename, proc_id, root_process, n_proc, verbose = verbose)

                end if

                if ( proc_id == root_process ) then

                    ! load the dielectric file to one processor
                    call load_scr_dielectric(dielectric_filename, verbose = verbose)

                end if

                ! send loaded data to all other processors
                call MPI_Bcast(scr_n_w_bins, 1, MPI_INTEGER, root_process,&
                  MPI_COMM_WORLD, err)
                call MPI_Bcast(scr_n_q_bins, 1, MPI_INTEGER, root_process,&
                  MPI_COMM_WORLD, err)
                call MPI_Bcast(scr_n_q_theta_bins, 1, MPI_INTEGER, root_process,&
                  MPI_COMM_WORLD, err)
                call MPI_Bcast(scr_n_q_phi_bins, 1, MPI_INTEGER, root_process,&
                  MPI_COMM_WORLD, err)

                call MPI_Bcast(scr_w_bin_width, 1, MPI_DOUBLE, root_process,&
                  MPI_COMM_WORLD, err)
                call MPI_Bcast(scr_q_bin_width, 1, MPI_DOUBLE, root_process,&
                  MPI_COMM_WORLD, err)

                if ( proc_id /= root_process ) then

                    !! allocate screen_mat
                    allocate(screen_mat(scr_n_w_bins, scr_n_q_bins, scr_n_q_theta_bins, &
                        scr_n_q_phi_bins))

                    screen_mat = (0.0_dp, 0.0_dp)

                end if

                call MPI_Bcast(screen_mat, size(screen_mat), MPI_DOUBLE_COMPLEX, &
                    root_process, MPI_COMM_WORLD, err)

                !! For some reason this does not work when load_dielectric_from_file is FALSE.
                !! When load_dielectric_from_file = TRUE all the processes
                !! open and read the file as expected. Seems like a file 
                !! permission issue with HDF5.
                !!
                !! Workaround : read by the main processor then broadcast
                !! the data to the other ones

                ! call MPI_Barrier(MPI_COMM_WORLD, ierr)
                ! call load_scr_dielectric(dielectric_filename, verbose = verbose)

            end if

        else

            if ( verbose ) then

                print*, '!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*
                print*, '    Input file for in medium parameters : ', trim(filename), ' does NOT exist.'
                print*
                print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*

            end if

            stop

        end if

    end subroutine

    subroutine save_in_med_scr(filename, verbose)

        implicit none
        character(len=*) :: filename

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id

        logical, optional :: verbose
        logical :: file_exists

        integer(HSIZE_T) :: dims1(1) = [1]

        integer :: error

        if ( verbose ) then

            print*, 'Saving screening parameters...'
            print*

        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'in_med_screening', group_id, error)

            ! ! write data
            call h5ltmake_dataset_string_f(file_id, 'in_med_screening/screen_type',&
                screen_type, error)
            if ( trim(screen_type) == 'analytic' ) then
                call h5ltmake_dataset_double_f(file_id, 'in_med_screening/di_e0', size(dims1), dims1,&
                    di_e0, error)
                call h5ltmake_dataset_double_f(file_id, 'in_med_screening/di_alpha', size(dims1), dims1,&
                    di_alpha, error)
                call h5ltmake_dataset_double_f(file_id, 'in_med_screening/di_q_tf', size(dims1), dims1,&
                    di_q_tf, error)
                call h5ltmake_dataset_double_f(file_id, 'in_med_screening/di_omega_p', size(dims1), dims1,&
                    di_omega_p, error)
            end if

            call h5fclose_f(file_id, error)
            call h5close_f(error)

        else

            if ( verbose ) then

                print*, '!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*
                print*, '    Output file : ', trim(filename), ' does NOT exist.'
                print*
                print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*
                print*

            end if

            stop

        end if

    end subroutine

    function screening(q_vec, omega) result(scr)

        implicit none

        real(dp) :: q_vec(3)
        real(dp) :: omega

        real(dp) :: q_mag
        real(dp) :: q_hat(3)
        real(dp) :: q_theta, q_phi

        integer :: w_bin_num, q_bin_num, q_theta_bin_num, q_phi_bin_num

        real(dp) :: scr

        if ( include_screen ) then

            if ( trim(screen_type) .eq. 'analytic' ) then
                
                scr = model_dielectric(q_vec, omega, &
                                        di_e0, &
                                        di_q_tf, &
                                        di_alpha, &
                                        di_omega_p)

            else if ( trim(screen_type) == 'numeric' ) then

                q_mag = norm2(q_vec)

                if ( q_mag >= 1.0e-8_dp ) then

                    q_hat = q_vec/q_mag

                    q_theta = get_theta(q_hat)
                    q_phi = get_phi(q_hat)

                    w_bin_num = 1 + floor(omega/scr_w_bin_width)
                    q_bin_num = 1 + floor(q_mag/scr_q_bin_width)
                    q_theta_bin_num = 1 + floor(q_theta/(pi/max(1.0_dp, 1.0_dp*scr_n_q_theta_bins)))
                    q_phi_bin_num = 1 + floor(q_phi/(2.0_dp*pi/max(1.0_dp, 1.0_dp*scr_n_q_phi_bins)))

                    if ( ( w_bin_num > scr_n_w_bins ) .or. &
                         ( q_bin_num > scr_n_q_bins) .or. & 
                         ( q_theta_bin_num > scr_n_q_theta_bins) .or. &
                         ( q_phi_bin_num > scr_n_q_phi_bins) ) then

                         scr = 1.0_dp

                    else

                        scr = abs( screen_mat(w_bin_num, q_bin_num, &
                            q_theta_bin_num, q_phi_bin_num) )

                    end if

                end if

            else

                scr = 1.0_dp

            end if

        else 

            scr = 1.0_dp

        end if

    end function

    function model_dielectric(q_vec, omega,&
            e0, q_tf, alpha, omega_p) result(di)
        !! Analytic form of the dielectric function
        !! Eq 6 from 
        !!      https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892
        implicit none

        real(dp) :: di

        real(dp) :: e0, q_tf, alpha, omega_p

        real(dp) :: q_vec(3)
        real(dp) :: q_mag
        real(dp) :: omega

        q_mag = norm2(q_vec)

        di = 1.0_dp + ( (e0 - 1.0_dp)**(-1) + &
            alpha*(q_mag/q_tf)**2 + &
            q_mag**4/(4*m_elec**2*omega_p**2) - &
            (omega/omega_p)**2)**(-1) 

    end function

end module
