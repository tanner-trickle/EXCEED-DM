module numerics_dielectric
    !! Numerics parameters specific to the dielectric calculation.

    use prec

    use info_messages

    implicit none

    type numerics_dielectric_t
        !! Numerics parameters specific to the dielectric calculation.

        integer :: n_val_max = 0
            !! Maximum number of valence bands
        integer :: n_cond_max = 0
            !! Maximum number of conduction bands
        integer, allocatable :: val_id_list(:)
            !! Dim : [n_val_max]
            !!
            !! List of valence ID's to compute for.
        integer, allocatable :: k_id_list(:)
            !! Dim : [n_k]
            !!
            !! List of k's to compute for.
        integer :: n_k_vec(3)
            !! Number of \( \mathbf{k} \) vectors in each direction.
        real(dp) :: q_grid_min(3)
            !! Offset parameter to map \( \mathbf{q} \) coordinates to
            !! indicies.
        integer :: n_q_grid(3)
            !! Number of points in the q grid with \( q < bins%n_q*bins%q_width \) 
        integer, allocatable :: n_q_bin(:, :, :)
            !! Dim : [n_q, n_q_theta, n_q_phi]
            !!
            !! Number of points in the dielectric \( \mathbf{q} \) grid in each \( \{q, \theta_q, \phi_q \} \) bin.

        contains

            procedure :: print => numerics_dielectric_print
            procedure :: load => numerics_dielectric_load_nml
            procedure :: save => numerics_dielectric_save
            procedure :: create_val_id_list => dielectric_create_val_id_list
            procedure :: create_k_id_list => dielectric_create_k_id_list
            procedure :: compute_n_q_bin

    end type

contains

    subroutine compute_n_q_bin(self, bins, PW_dataset, FFT_grid, proc_id, root_process, verbose)
        !! Find the number of \( \mathbf{q} \) points from the uniformly spaced grid that are in each
        !! in a \( \{q, \theta_q, \phi_q \} \) bin.

        use MPI_util
        use PW_dataset_type
        use FFT_util
        use math_mod
        use bins_dielectric_type

        class(numerics_dielectric_t) :: self
        type(bins_dielectric_t) :: bins
        type(PW_dataset_t) :: PW_dataset
        type(FFT_grid_t) :: FFT_grid
        integer :: proc_id, root_process

        type(parallel_manager_t) :: ki_kf_manager
        integer :: k_id_list(PW_dataset%n_k)
        integer :: k, j, ki_id, kf_id, g1, g2, g3
        integer :: job_id

        integer, allocatable :: job_id_to_ki_kf(:, :)

        integer, allocatable :: total_n_q_bin(:, :, :)

        logical, optional :: verbose

        integer :: err

        real(dp) :: q_vec(3), q_hat(3)
        real(dp) :: q_mag, q_theta, q_phi

        integer :: q_theta_bin, q_phi_bin, q_bin

        real(dp) :: dk_grid_min(3)
        real(dp) :: dk_grid_max(3)

        real(dp) :: G_grid_min(3)
        real(dp) :: G_grid_max(3)

        real(dp) :: q_grid_min(3)
        real(dp) :: q_grid_max(3)

        integer :: n_q_grid(3)

        integer :: q1, q2, q3

        integer :: G_red(3)

        real(dp) :: q_red(3)

        if ( verbose ) then

            print*, 'Starting calculation of number of q points in each bin...'
            print*

        end if

        allocate(self%n_q_bin(bins%n_q, bins%n_q_theta, bins%n_q_phi))
        self%n_q_bin = 0

        dk_grid_min(1) = minval(PW_dataset%k_grid_red(:, 1)) &
            - maxval(PW_dataset%k_grid_red(:, 1))
        dk_grid_min(2) = minval(PW_dataset%k_grid_red(:, 2)) &
            - maxval(PW_dataset%k_grid_red(:, 2))
        dk_grid_min(3) = minval(PW_dataset%k_grid_red(:, 3)) &
            - maxval(PW_dataset%k_grid_red(:, 3))

        dk_grid_max(1) = maxval(PW_dataset%k_grid_red(:, 1)) &
            - minval(PW_dataset%k_grid_red(:, 1))
        dk_grid_max(2) = maxval(PW_dataset%k_grid_red(:, 2)) &
            - minval(PW_dataset%k_grid_red(:, 2))
        dk_grid_max(3) = maxval(PW_dataset%k_grid_red(:, 3)) &
            - minval(PW_dataset%k_grid_red(:, 3))

        G_grid_min = 0.0_dp
        G_grid_max = 0.0_dp

        do g3 = 1, FFT_grid%n_grid(3)
            do g2 = 1, FFT_grid%n_grid(2)
                do g1 = 1, FFT_grid%n_grid(1)

                    G_red = FFT_grid%sym_G_grid_red(g1, g2, g3, :)
                    
                    G_grid_min(1) = min(1.0_dp*G_red(1), G_grid_min(1))
                    G_grid_min(2) = min(1.0_dp*G_red(2), G_grid_min(2))
                    G_grid_min(3) = min(1.0_dp*G_red(3), G_grid_min(3))

                    G_grid_max(1) = max(1.0_dp*G_red(1), G_grid_max(1))
                    G_grid_max(2) = max(1.0_dp*G_red(2), G_grid_max(2))
                    G_grid_max(3) = max(1.0_dp*G_red(3), G_grid_max(3))

                end do
            end do
        end do

        q_grid_min = dk_grid_min + G_grid_min
        q_grid_max = dk_grid_min + G_grid_max

        n_q_grid(1) = 1 + int( self%n_k_vec(1)*( q_grid_max(1) - q_grid_min(1) ) )
        n_q_grid(2) = 1 + int( self%n_k_vec(2)*( q_grid_max(2) - q_grid_min(2) ) )
        n_q_grid(3) = 1 + int( self%n_k_vec(3)*( q_grid_max(3) - q_grid_min(3) ) )

        do q1 = 1, n_q_grid(1)
            do q2 = 1, n_q_grid(2)
                do q3 = 1, n_q_grid(3)

                    q_red(1) = q_grid_min(1) + ( q1 - 1.0_dp )/self%n_k_vec(1) 
                    q_red(2) = q_grid_min(2) + ( q2 - 1.0_dp )/self%n_k_vec(2) 
                    q_red(3) = q_grid_min(3) + ( q3 - 1.0_dp )/self%n_k_vec(3) 

                    q_vec = matmul(PW_dataset%k_red_to_xyz, q_red)

                    q_mag = norm2(q_vec)

                    if ( ( q_mag > 1.0e-8_dp ) .and. &
                        ( q_mag < bins%n_q*bins%q_width ) ) then

                        q_hat = q_vec/q_mag

                        q_theta = get_theta(q_hat)
                        q_phi = get_phi(q_hat)

                        q_theta_bin = Q_func(q_theta, 0.0_dp,&
                            pi/max(1.0_dp, 1.0_dp*bins%n_q_theta), bins%n_q_theta)

                        q_phi_bin = Q_func(q_phi, 0.0_dp,&
                            2.0_dp*pi/max(1.0_dp, 1.0_dp*bins%n_q_phi), bins%n_q_phi)

                        q_bin = 1 + floor(q_mag/bins%q_width)

                        self%n_q_bin(q_bin, q_theta_bin, q_phi_bin) = &
                            self%n_q_bin(q_bin, q_theta_bin, q_phi_bin) + 1

                    end if

                end do
            end do
        end do

        if ( verbose ) then

            print*, 'Done computing number of q points in each bin!'
            print*

        end if

    end subroutine

    subroutine dielectric_create_k_id_list(self, n_k)
        !! Specify the indicies for each \( \mathbf{k} \) point that should be included. 
        !! Specific to the dielectric calculation. 

        implicit none

        class(numerics_dielectric_t) :: self
        integer :: n_k

        integer :: k

        allocate(self%k_id_list(n_k))

        do k = 1, n_k
            self%k_id_list(k) = k
        end do

    end subroutine

    subroutine dielectric_create_val_id_list(self, n_val)
        !! Specify the indicies for each valence band that should be included. 
        !! Specific to the dielectric calculation. 

        implicit none
        class(numerics_dielectric_t) :: self

        integer :: n_val

        integer :: j

        allocate(self%val_id_list(self%n_val_max))

        do j = 1, self%n_val_max
            self%val_id_list(j) = n_val - j + 1
        end do

    end subroutine

    subroutine numerics_dielectric_save(self, filename, verbose)
        !! Saves `numerics_dielectric`.

        use hdf5
        use h5lt

        implicit none

        class(numerics_dielectric_t) :: self
        character(len=*) :: filename
        logical, optional :: verbose

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id
        logical :: file_exists
        integer(HSIZE_T) :: dims1(1) = [1]
        integer :: error

        if ( verbose ) then
            print*, 'Saving numerics - dielectric parameters...'
            print*
        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'numerics_dielectric', group_id, error)

            ! write data
            call h5ltmake_dataset_int_f(file_id, &
                'numerics_dielectric/n_val_max', &
                size(dims1), dims1, &
                self%n_val_max, &
                error)
            call h5ltmake_dataset_int_f(file_id, &
                'numerics_dielectric/n_cond_max', &
                size(dims1), dims1, &
                self%n_cond_max, &
                error)

            dims1 = [3]

            call h5ltmake_dataset_int_f(file_id, &
                'numerics_dielectric/n_k_vec', &
                size(dims1), dims1, &
                self%n_k_vec, &
                error)

            call h5fclose_f(file_id, error)
            call h5close_f(error)

        else

            call print_error_message(&
                'Output file : '//trim(filename)//' does NOT exist.', &
                verbose = verbose)
            stop

        end if

    end subroutine

    subroutine numerics_dielectric_load_nml(self, filename, n_val, n_cond, verbose)
        !! Loads `numerics_dielectric` parameters from a namelist.

        implicit none

        class(numerics_dielectric_t) :: self
        character(len=*) :: filename
        integer :: n_val, n_cond
        logical, optional :: verbose

        logical :: file_exists
        integer :: error

        integer :: n_val_max = 0
            !! Maximum number of valence bands
        integer :: n_cond_max = 0
            !! Maximum number of conduction bands
        integer :: n_k_vec(3)

        NAMELIST /numerics_dielectric/ n_val_max, &
                                 n_cond_max, &
                                 n_k_vec

        if ( verbose ) then
            print*, 'Loading numerics parameters for dielectric calculation...'
            print*
        end if

        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            open(100, file = trim(filename), iostat = error)
            read(100, nml = numerics_dielectric, iostat = error)
            close(100)

            if ( n_val_max == 0 ) then
                self%n_val_max = n_val
            else
                self%n_val_max = min(n_val, n_val_max)
            end if

            if ( n_cond_max == 0 ) then
                self%n_cond_max = n_cond
            else
                self%n_cond_max = min(n_cond, n_cond_max)
            end if

            self%n_k_vec = n_k_vec

            call self%print(verbose = verbose)

        else

            call print_error_message(&
                'Input file for numerics-dielectric parameters : '//trim(filename)//' does NOT exist.', &
                verbose = verbose)
            stop

        end if

    end subroutine

    subroutine numerics_dielectric_print(self, verbose)
        !! Prints `numerics_dielectric` components.

        implicit none

        class(numerics_dielectric_t) :: self
        logical, optional :: verbose

        if ( verbose ) then

            call print_section_seperator()
            print*, '    ---------------------'
            print*, '    Numerics - Dielectric'
            print*, '    ---------------------'
            print*
            print*, '        Maximum number of valence bands    : ', trim(adjustl(int_to_str(self%n_val_max)))
            print*, '        Maximum number of conduction bands : ', trim(adjustl(int_to_str(self%n_cond_max)))
            print*
            print*, '        Number of k points in each direction : ', self%n_k_vec
            print*
            call print_section_seperator()
            print*

        end if

    end subroutine

    subroutine check_dielectric_memory(n, verbose)
        !! Checks to see if the dielectric is going to take up too
        !! much memory.
        
        use info_messages

        implicit none
        
        integer :: n
        logical, optional :: verbose

        if ( 16.0_dp*n >= 1.0e10_dp ) then

            call print_warning_message('Attempting to store more than 10 GB of dielectric data on a single processor.'//&
                    'Try increasing the number of processors or decreasing the number '//&
                    'of energy/momentum bins in the dielectric namelist.', verbose = verbose)

        end if

    end subroutine

end module
