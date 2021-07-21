module calc_dielectric
    !! Compute the dielectric tensor from Pi_11 from valence -> conduction transitions.

    use mpi
    use hdf5
    use h5lt

    use prec
    use constants
    use math_mod

    use control_input
    use material_input
    use DFT_parameters
    use FFT_util
    use io_input
    use transition_form_factor

    use di_transition
    use dielectric_input
    use di_grid

    implicit none

    complex(dp), allocatable :: dielec(:, :, :, :)
        !! Dim : [ n_omega_bins + 1, n_q_bins + 1, n_q_theta_bins + 1, n_q_phi_bins + 1 ]
        !!
        !! eps(omega, q) = 1 - \frac{e^2}{q^2} Pi_11(q, omega)
        !!
        !! Pi_11 = (1/V) sum_{II'} \frac{1}{omega - omega_{I'} - omega_{I} + i delta} 
        !!      |<I'| e^{i q x} |I>|^2
        !!
        !! Units : None

    complex(dp), allocatable :: dielec_t(:, :, :, :, :)
        !! Dim : [ n_tran_per_proc, n_omega_bins + 1, n_q_bins + 1, n_q_theta_bins + 1, n_q_phi_bins + 1 ]
        !!
        !! Each processors contribution to the dielectric
        !!
        !! Units : None

    interface compute_dielectric
        module procedure compute_dielectric_no_spin
        module procedure compute_dielectric_spin
    end interface

contains

    subroutine check_dielectric_memory(verbose)
        !! Checks to see if the dielectric is going to take up too
        !! much memory.

        implicit none
        
        logical, optional :: verbose

        if ( 16.0_dp*(di_n_q_bins)*&
            (di_n_q_phi_bins)*&
            (di_n_q_theta_bins)*&
            (di_n_omega_bins)*&
            (di_n_tran_per_proc) >= 1.0e10_dp ) then

            if ( verbose ) then

                print*, '~~~ WARNING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
                print*
                print*, '    Attempting to store more than 10 GB of dielectric data on a single processor.'
                print*, '    Try increasing the number of processors or decreasing the number '
                print*, '    of energy/momentum bins in the dielectric namelist.'
                print*
                print*, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
                print*

            end if

        end if

    end subroutine

    subroutine run_dielectric_calc(nml_filename, filename, DFT_input_filename,&
            proc_id, root_process, n_proc, verbose)
        !! Compute the dielectric binned in energy and momentum transfer.

        implicit none

        character(len=*) :: filename
        character(len=*) :: nml_filename
        character(len=*) :: DFT_input_filename

        logical, optional :: verbose

        integer :: proc_id, root_process
        integer :: n_proc

        integer :: tran_id
        integer :: val_id, cond_id

        integer :: n_FFT_grid(3)
        integer :: n_FFT

        integer :: wfc_fft_plan(8)
        integer :: Tif_fft_plan(8)

        complex(dp), allocatable :: wfc_FT_i(:, :)
        complex(dp), allocatable :: wfc_FT_f(:, :)

        complex(dp), allocatable :: wfc_FT_i_s(:, :, :)
        complex(dp), allocatable :: wfc_FT_f_s(:, :, :)

        integer :: t, i

        integer :: n_q_grid(3)
        real(dp) :: q_grid_min(3)

        integer :: status(MPI_STATUS_SIZE)
        integer :: tag = 0
        integer :: err

        if ( verbose ) then
            print*, 'Starting dielectric calculation...'
            print*
        end if

        !! load inputs
        ! call load_dielectric_input(nml_filename, verbose = verbose)
        call load_DFT_parameters(DFT_input_filename, verbose = verbose)
        call do_scissor_correction(band_gap, verbose = verbose)

        call di_set_job_table(n_proc, di_n_init, di_n_fin, verbose = verbose)

        ! doubled to avoid wrapping problems
        n_FFT_grid(1) = 2*(maxval(in_G_grid_red(:, 1)) - minval(in_G_grid_red(:, 1))) + 1
        n_FFT_grid(2) = 2*(maxval(in_G_grid_red(:, 2)) - minval(in_G_grid_red(:, 2))) + 1
        n_FFT_grid(3) = 2*(maxval(in_G_grid_red(:, 3)) - minval(in_G_grid_red(:, 3))) + 1

        n_FFT = n_FFT_grid(1)*n_FFT_grid(2)*n_FFT_grid(3)

        if ( verbose ) then

            print*, '----------------------------------------'
            print*, '    --------------'
            print*, '    Dielectric FFT'
            print*, '    --------------'
            print*
            print*, '        Dimensions of FFT grid : '
            print*, '        ', n_FFT_grid
            print*
            print*, '        Number of G points in FFT grid = ', n_FFT
            print*
            print*, '----------------------------------------'
            print*

        end if

        call set_fft_plan_backward_3d(n_FFT_grid, wfc_fft_plan)
        call set_fft_plan_backward_3d(n_FFT_grid, Tif_fft_plan)

        call define_q_grid(di_q_bin_width*di_n_q_bins, &
            k_red_to_xyz, n_q_grid, q_grid_min, n_k_vec, n_FFT_grid, verbose = .TRUE.)

        !! Don't need this here, get G vector from function in FFT_util.
        !! Less than ideal because we don't precompute the G vectors, but 
        !! we avoid reallocating the array in FFT_util twice.

        ! call set_sym_FFT_G_grid_xyz(n_FFT_grid, k_red_to_xyz, verbose = verbose)

        if ( .not. load_dielectric_from_file ) then

            !! hard part of calculation starts here

            call check_dielectric_memory(verbose = verbose)

            if ( ( proc_id == root_process ) .and. ( timer ) ) then

                call time_dielectric_calc(DFT_input_filename, 1, &
                                          n_FFT_grid, 1, wfc_fft_plan, &
                                          Tif_fft_plan, n_q_grid, q_grid_min, n_k_vec, &
                                          verbose = verbose)

            end if

            !! allocate memory

            if ( proc_id == root_process ) then

                allocate(dielec(di_n_omega_bins, di_n_q_bins, &
                    di_n_q_theta_bins, di_n_q_phi_bins))

                !! factor of 1 in the dielectric formula
                dielec = (1.0_dp, 0.0_dp)

            end if

            !! initialize variables
            allocate(dielec_t(di_n_tran_per_proc, di_n_omega_bins, di_n_q_bins, &
                di_n_q_theta_bins, di_n_q_phi_bins))

            dielec_t = (0.0_dp, 0.0_dp)

            if ( include_spin ) then

                allocate(wfc_FT_i_s(n_k, n_in_G, 2))
                allocate(wfc_FT_f_s(n_k, n_in_G, 2))

            else

                allocate(wfc_FT_i(n_k, n_in_G))
                allocate(wfc_FT_f(n_k, n_in_G))

            end if

            !! do calculation for each job

            do t = 1, di_n_tran_per_proc

                tran_id = di_job_table(proc_id + 1, t)

                !! only compute what the processor should be
                !! computing
                if ( tran_id > 0 ) then

                    val_id = di_tran_to_init_fin_id(tran_id, 1)
                    cond_id = di_tran_to_init_fin_id(tran_id, 2) + n_val

                    if ( include_spin ) then

                        call get_in_wfc_FT(DFT_input_filename, val_id, wfc_FT_i_s)
                        call get_in_wfc_FT(DFT_input_filename, cond_id, wfc_FT_f_s)

                        call compute_dielectric(dielec_t(t, :, :, :, :), &
                        val_id, cond_id, wfc_FT_i_s, wfc_FT_f_s, n_FFT_grid, n_k, &
                            wfc_FFT_plan, Tif_FFT_plan, n_q_grid, q_grid_min, &
                            n_k_vec, verbose = verbose)

                    else

                        call get_in_wfc_FT(DFT_input_filename, val_id, wfc_FT_i)
                        call get_in_wfc_FT(DFT_input_filename, cond_id, wfc_FT_f)

                        !! do the calculation for an individual job

                        call compute_dielectric(dielec_t(t, :, :, :, :), &
                        val_id, cond_id, wfc_FT_i, wfc_FT_f, n_FFT_grid, n_k, &
                            wfc_FFT_plan, Tif_FFT_plan, n_q_grid, q_grid_min, &
                            n_k_vec, verbose = verbose)

                    end if

                end if

            end do

            !! collect and sum results from all the processors, MPI part

            if ( proc_id /= root_process ) then

                !! send to main processor
                call MPI_SEND(dielec_t, &
                   size(dielec_t), MPI_DOUBLE_COMPLEX, &
                   root_process, tag, MPI_COMM_WORLD, err)

            end if

            if ( proc_id == root_process ) then

                !! add the main processors contribution
                dielec = dielec + sum(dielec_t, 1)

                !! receive the other processors contributions
                do i = 1, n_proc

                    if ( i - 1 /= root_process ) then

                        call MPI_RECV(dielec_t, &
                           size(dielec_t), MPI_DOUBLE_COMPLEX, &
                           i - 1, MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                        dielec = dielec + sum(dielec_t, 1)

                    end if

                end do

            end if

            !! save dielectric data

            if ( proc_id == root_process ) then

                call save_dielectric(filename, verbose = verbose)

            end if

        end if

    end subroutine

    subroutine time_dielectric_calc(DFT_input_filename, tran_id, &
            n_FFT_grid, k_cut, wfc_fft_plan, Tif_fft_plan, n_q_grid, q_grid_min, n_k_vec, verbose)

        use timing

        implicit none

        character(len=*) :: DFT_input_filename

        integer :: tran_id

        integer :: n_FFT_grid(3)

        integer :: k_cut

        integer :: wfc_fft_plan(8)
        integer :: Tif_fft_plan(8)

        integer :: n_q_grid(3)

        real(dp) :: q_grid_min(3)

        integer :: n_k_vec(3)

        complex(dp), allocatable :: wfc_FT_i(:, :)
        complex(dp), allocatable :: wfc_FT_f(:, :)

        complex(dp), allocatable :: wfc_FT_i_s(:, :, :)
        complex(dp), allocatable :: wfc_FT_f_s(:, :, :)

        logical, optional :: verbose

        integer :: val_id, cond_id

        complex(dp) :: dielec_test(di_n_omega_bins, di_n_q_bins, di_n_q_theta_bins, di_n_q_phi_bins)

        if ( verbose ) then

            print*, 'Timing dielectric calculation...'
            print*

        end if

        if ( include_spin ) then

            allocate(wfc_FT_i_s(n_k, n_in_G, 2))
            allocate(wfc_FT_f_s(n_k, n_in_G, 2))

        else

            allocate(wfc_FT_i(n_k, n_in_G))
            allocate(wfc_FT_f(n_k, n_in_G))

        end if

        val_id = di_tran_to_init_fin_id(tran_id, 1)
        cond_id = di_tran_to_init_fin_id(tran_id, 2) + n_val

        if ( include_spin ) then

            call get_in_wfc_FT(DFT_input_filename, val_id, wfc_FT_i_s)
            call get_in_wfc_FT(DFT_input_filename, cond_id, wfc_FT_f_s)

            time(3) = MPI_Wtime()

            call compute_dielectric(dielec_test, &
            val_id, cond_id, wfc_FT_i_s, wfc_FT_f_s, n_FFT_grid, k_cut, &
                wfc_FFT_plan, Tif_FFT_plan, n_q_grid, q_grid_min, &
                n_k_vec, verbose = verbose)

            time(4) = MPI_Wtime()

        else

            call get_in_wfc_FT(DFT_input_filename, val_id, wfc_FT_i)
            call get_in_wfc_FT(DFT_input_filename, cond_id, wfc_FT_f)

            !! do the calculation for an individual job

            time(3) = MPI_Wtime()

            call compute_dielectric(dielec_test, &
            val_id, cond_id, wfc_FT_i, wfc_FT_f, n_FFT_grid, k_cut, &
                wfc_FFT_plan, Tif_FFT_plan, n_q_grid, q_grid_min, &
                n_k_vec, verbose = verbose)

            time(4) = MPI_Wtime()

        end if

        if ( verbose ) then

            print*, '----------------------------------------'
            print*, '    --------------------------'
            print*, '    Timing (TEST - dielectric)'
            print*, '    --------------------------'
            print*
            print*, '        (TEST) Run time: '
            print*, '            ', trim(pretty_time_format(time(4) - time(3)))
            print*
            print*, '        Expected run time for dielectric calculation : '
            print*, '            ', trim(pretty_time_format(&
                di_n_tran_per_proc*n_k**2*(time(4) - time(3))&
                ))
            print*
            print*, '----------------------------------------'
            print*

        end if

    end subroutine

    subroutine compute_dielectric_spin(di, val_id, cond_id, wfc_FT_i, wfc_FT_f, n_FFT_grid, &
            k_cut, wfc_FFT_plan, Tif_FFT_plan, n_q_grid, q_grid_min, &
            n_k_vec, verbose)
        !! Compute the contribution to the dimensionless dielectric
        !! from a given i -> f transition.
        
        implicit none

        integer :: n_k_vec(3)
        integer :: n_q_grid(3)

        real(dp) :: q_grid_min(3)

        complex(dp) :: di(:, :, :, :)

        complex(dp) :: di_unbinned(di_n_omega_bins, n_q_grid(1), n_q_grid(2), n_q_grid(3))

        integer :: ki, kf

        integer :: k_cut

        integer :: val_id, cond_id

        integer :: n_FFT_grid(3)

        complex(dp) :: wfc_FT_i(:, :, :)
        complex(dp) :: wfc_FT_f(:, :, :)

        complex(dp) :: wfc_FT_i_exp(2, n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3))
        complex(dp) :: wfc_FT_f_exp(2, n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3))

        complex(dp) :: wfc_i(2, n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3))
        complex(dp) :: wfc_f(2, n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3))

        integer :: wfc_fft_plan(8)
        integer :: Tif_fft_plan(8)

        logical, optional :: verbose

        real(dp) :: omega
        real(dp) :: q_vec(3)

        integer :: g1, g2, g3

        integer :: n_FFT

        real(dp) :: f_sq(n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3))

        integer :: w
        complex(dp) :: elec_props(di_n_omega_bins)

        real(dp) :: d_omega

        real(dp) :: delta

        integer :: q_bin, q_theta_bin, q_phi_bin

        real(dp) :: q_hat(3)
        real(dp) :: q_mag, q_phi, q_theta

        real(dp) :: Ei, Ef

        integer :: num_q_in_bins(di_n_q_bins, di_n_q_theta_bins, di_n_q_phi_bins)

        integer :: q1, q2, q3

        real(dp) :: q_red(3)

        integer :: q_ind(3)

        di = (0.0_dp, 0.0_dp)
        di_unbinned = (0.0_dp, 0.0_dp)

        wfc_FT_i_exp = (0.0_dp, 0.0_dp)
        wfc_FT_f_exp = (0.0_dp, 0.0_dp)
        wfc_i = (0.0_dp, 0.0_dp)
        wfc_f = (0.0_dp, 0.0_dp)

        n_FFT = n_FFT_grid(1)*n_FFT_grid(2)*n_FFT_grid(3)

        num_q_in_bins = 0

        do ki = 1, k_cut

            ! expand both components
            call expand_wfc_FT_for_FFT(n_FFT_grid, wfc_FT_i(ki, :, 1), wfc_FT_i_exp(1, :, :, :),&
                verbose = verbose)
            call expand_wfc_FT_for_FFT(n_FFT_grid, wfc_FT_i(ki, :, 2), wfc_FT_i_exp(2, :, :, :),&
                verbose = verbose)

            ! Fourier transform both components
            call dfftw_execute_dft(wfc_fft_plan, wfc_FT_i_exp(1, :, :, :), wfc_i(1, :, :, :)) 
            call dfftw_execute_dft(wfc_fft_plan, wfc_FT_i_exp(2, :, :, :), wfc_i(2, :, :, :)) 

            Ei = energy_bands(ki, val_id)

            do kf = 1, k_cut

                Ef = energy_bands(kf, cond_id)

                call expand_wfc_FT_for_FFT(n_FFT_grid, wfc_FT_f(kf, :, 1), wfc_FT_f_exp(1, :, :, :),&
                    verbose = verbose)
                call expand_wfc_FT_for_FFT(n_FFT_grid, wfc_FT_f(kf, :, 2), wfc_FT_f_exp(2, :, :, :),&
                    verbose = verbose)

                call dfftw_execute_dft(wfc_fft_plan, wfc_FT_f_exp(1, :, :, :), wfc_f(1, :, :, :)) 
                call dfftw_execute_dft(wfc_fft_plan, wfc_FT_f_exp(2, :, :, :), wfc_f(2, :, :, :)) 

                call calc_tff_vc(f_sq, wfc_i, wfc_f, n_FFT_grid, Tif_fft_plan, verbose = .FALSE.)

                do w = 1, di_n_omega_bins

                    omega = di_omega_bin_width*(w - 0.5_dp)

                    d_omega = Ef - Ei

                    delta = di_width_func(omega)

                    elec_props(w) = ( omega - d_omega + ii*delta )**(-1) - &
                        ( omega + d_omega - ii*delta )**(-1)

                end do

                do g3 = 1, n_FFT_grid(3)
                    do g2 = 1, n_FFT_grid(2)
                        do g1 = 1, n_FFT_grid(1)

                            q_vec = k_grid_xyz(kf, :) - k_grid_xyz(ki, :) + &
                                get_sym_FFT_G_grid_xyz(n_FFT_grid, [g1, g2, g3], &
                                    k_red_to_xyz, verbose = .FALSE.)

                            q_red = k_grid_red(kf, :) - k_grid_red(ki, :) + &
                                get_sym_FFT_G_grid_red(n_FFT_grid, [g1, g2, g3], &
                                    verbose = .FALSE.)

                            q_mag = norm2(q_vec)

                            if ( ( q_mag > 1.0e-8_dp ) .and. &
                                ( q_mag < di_n_q_bins*di_q_bin_width ) ) then

                                q_ind(1) = 1 + int(n_k_vec(1)*(q_red(1) - q_grid_min(1)))
                                q_ind(2) = 1 + int(n_k_vec(2)*(q_red(2) - q_grid_min(2)))
                                q_ind(3) = 1 + int(n_k_vec(3)*(q_red(3) - q_grid_min(3)))

                                di_unbinned(:, q_ind(1), q_ind(2), q_ind(3)) = &
                                    di_unbinned(:, q_ind(1), q_ind(2), q_ind(3)) + &
                                    (-1.0_dp)*(spin_degen/2.0_dp)*(e_EM**2/q_mag**2)*(pc_vol)**(-1)*&
                                    k_weight(ki)*(1.0_dp*n_FFT)**(-2)*&
                                    elec_props(:)*f_sq(g1, g2, g3)

                            end if

                        end do
                    end do
                end do

            end do

        end do

        ! now bin the unbinned dielectric
        do q3 = 1, n_q_grid(3)
            do q2 = 1, n_q_grid(2)
                do q1 = 1, n_q_grid(1) 

                    q_red(1) = (1.0_dp/n_k_vec(1))*( q1 - 1 ) + q_grid_min(1)
                    q_red(2) = (1.0_dp/n_k_vec(2))*( q2 - 1 ) + q_grid_min(2)
                    q_red(3) = (1.0_dp/n_k_vec(3))*( q3 - 1 ) + q_grid_min(3)

                    q_vec = matmul(k_red_to_xyz, q_red)

                    q_mag = norm2(q_vec)

                    if ( ( q_mag > 1.0e-8_dp ) .and. &
                        ( q_mag < di_n_q_bins*di_q_bin_width ) ) then

                        q_hat = q_vec/q_mag

                        q_theta = get_theta(q_hat)
                        q_phi = get_phi(q_hat)

                        q_theta_bin = Q_func(q_theta, 0.0_dp,&
                            pi/max(1.0_dp, 1.0_dp*di_n_q_theta_bins), di_n_q_theta_bins)

                        q_phi_bin = Q_func(q_phi, 0.0_dp,&
                            2.0_dp*pi/max(1.0_dp, 1.0_dp*di_n_q_phi_bins), di_n_q_phi_bins)

                        q_bin = 1 + floor(q_mag/di_q_bin_width)

                        num_q_in_bins(q_bin, q_theta_bin, q_phi_bin) = &
                            num_q_in_bins(q_bin, q_theta_bin, q_phi_bin) + 1

                        di(:, q_bin, q_theta_bin, q_phi_bin) = &
                            di(:, q_bin, q_theta_bin, q_phi_bin) + &
                            di_unbinned(:, q1, q2, q3)

                    end if

                end do
            end do
        end do

        ! divide by the number of elements in each bin
        do q1 = 1, di_n_q_bins
            do q2 = 1, di_n_q_theta_bins
                do q3 = 1, di_n_q_phi_bins 

                    di(:, q1, q2, q3) = di(:, q1, q2, q3)/max(1, num_q_in_bins(q1, q2, q3))

                end do
            end do
        end do

    end subroutine 

    subroutine compute_dielectric_no_spin(di, val_id, cond_id, wfc_FT_i, wfc_FT_f, n_FFT_grid, &
            k_cut, wfc_FFT_plan, Tif_FFT_plan, n_q_grid, q_grid_min, &
            n_k_vec, verbose)
        !! Compute the contribution to the dimensionless dielectric
        !! from a given i -> f transition.
        
        implicit none

        integer :: n_k_vec(3)
        integer :: n_q_grid(3)

        real(dp) :: q_grid_min(3)

        complex(dp) :: di(:, :, :, :)

        complex(dp) :: di_unbinned(di_n_omega_bins, n_q_grid(1), n_q_grid(2), n_q_grid(3))

        integer :: ki, kf

        integer :: k_cut

        integer :: val_id, cond_id

        integer :: n_FFT_grid(3)

        complex(dp) :: wfc_FT_i(:, :)
        complex(dp) :: wfc_FT_f(:, :)

        complex(dp) :: wfc_FT_i_exp(n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3))
        complex(dp) :: wfc_FT_f_exp(n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3))

        complex(dp) :: wfc_i(n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3))
        complex(dp) :: wfc_f(n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3))

        integer :: wfc_fft_plan(8)
        integer :: Tif_fft_plan(8)

        logical, optional :: verbose

        real(dp) :: omega
        real(dp) :: q_vec(3)

        integer :: g1, g2, g3

        integer :: n_FFT

        real(dp) :: f_sq(n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3))

        integer :: w
        complex(dp) :: elec_props(di_n_omega_bins)

        real(dp) :: d_omega

        real(dp) :: delta

        integer :: q_bin, q_theta_bin, q_phi_bin

        real(dp) :: q_hat(3)
        real(dp) :: q_mag, q_phi, q_theta

        real(dp) :: Ei, Ef

        integer :: num_q_in_bins(di_n_q_bins, di_n_q_theta_bins, di_n_q_phi_bins)

        integer :: q1, q2, q3

        real(dp) :: q_red(3)

        integer :: q_ind(3)

        di = (0.0_dp, 0.0_dp)
        di_unbinned = (0.0_dp, 0.0_dp)

        wfc_FT_i_exp = (0.0_dp, 0.0_dp)
        wfc_FT_f_exp = (0.0_dp, 0.0_dp)
        wfc_i = (0.0_dp, 0.0_dp)
        wfc_f = (0.0_dp, 0.0_dp)

        n_FFT = n_FFT_grid(1)*n_FFT_grid(2)*n_FFT_grid(3)

        num_q_in_bins = 0

        do ki = 1, k_cut

            call expand_wfc_FT_for_FFT(n_FFT_grid, wfc_FT_i(ki, :), wfc_FT_i_exp,&
                verbose = verbose)

            call dfftw_execute_dft(wfc_fft_plan, wfc_FT_i_exp, wfc_i) 

            Ei = energy_bands(ki, val_id)

            do kf = 1, k_cut

                Ef = energy_bands(kf, cond_id)

                call expand_wfc_FT_for_FFT(n_FFT_grid, wfc_FT_f(kf, :), wfc_FT_f_exp,&
                    verbose = verbose)

                call dfftw_execute_dft(wfc_fft_plan, wfc_FT_f_exp, wfc_f) 

                call calc_tff_vc(f_sq, wfc_i, wfc_f, n_FFT_grid, Tif_fft_plan, verbose = .FALSE.)

                do w = 1, di_n_omega_bins

                    omega = di_omega_bin_width*(w - 0.5_dp)

                    d_omega = Ef - Ei

                    delta = di_width_func(omega)

                    elec_props(w) = ( omega - d_omega + ii*delta )**(-1) - &
                        ( omega + d_omega - ii*delta )**(-1)

                end do

                do g3 = 1, n_FFT_grid(3)
                    do g2 = 1, n_FFT_grid(2)
                        do g1 = 1, n_FFT_grid(1)

                            q_vec = k_grid_xyz(kf, :) - k_grid_xyz(ki, :) + &
                                get_sym_FFT_G_grid_xyz(n_FFT_grid, [g1, g2, g3], &
                                    k_red_to_xyz, verbose = .FALSE.)

                            q_red = k_grid_red(kf, :) - k_grid_red(ki, :) + &
                                get_sym_FFT_G_grid_red(n_FFT_grid, [g1, g2, g3], &
                                    verbose = .FALSE.)

                            q_mag = norm2(q_vec)

                            if ( ( q_mag > 1.0e-8_dp ) .and. &
                                ( q_mag < di_n_q_bins*di_q_bin_width ) ) then

                                q_ind(1) = 1 + int(n_k_vec(1)*(q_red(1) - q_grid_min(1)))
                                q_ind(2) = 1 + int(n_k_vec(2)*(q_red(2) - q_grid_min(2)))
                                q_ind(3) = 1 + int(n_k_vec(3)*(q_red(3) - q_grid_min(3)))

                                di_unbinned(:, q_ind(1), q_ind(2), q_ind(3)) = &
                                    di_unbinned(:, q_ind(1), q_ind(2), q_ind(3)) + &
                                    (-1.0_dp)*(spin_degen/2.0_dp)*(e_EM**2/q_mag**2)*(pc_vol)**(-1)*&
                                    k_weight(ki)*(1.0_dp*n_FFT)**(-2)*&
                                    elec_props(:)*f_sq(g1, g2, g3)

                            end if

                        end do
                    end do
                end do

            end do

        end do

        ! now bin the unbinned dielectric
        do q3 = 1, n_q_grid(3)
            do q2 = 1, n_q_grid(2)
                do q1 = 1, n_q_grid(1) 

                    q_red(1) = (1.0_dp/n_k_vec(1))*( q1 - 1 ) + q_grid_min(1)
                    q_red(2) = (1.0_dp/n_k_vec(2))*( q2 - 1 ) + q_grid_min(2)
                    q_red(3) = (1.0_dp/n_k_vec(3))*( q3 - 1 ) + q_grid_min(3)

                    q_vec = matmul(k_red_to_xyz, q_red)

                    q_mag = norm2(q_vec)

                    if ( ( q_mag > 1.0e-8_dp ) .and. &
                        ( q_mag < di_n_q_bins*di_q_bin_width ) ) then

                        q_hat = q_vec/q_mag

                        q_theta = get_theta(q_hat)
                        q_phi = get_phi(q_hat)

                        q_theta_bin = Q_func(q_theta, 0.0_dp,&
                            pi/max(1.0_dp, 1.0_dp*di_n_q_theta_bins), di_n_q_theta_bins)

                        q_phi_bin = Q_func(q_phi, 0.0_dp,&
                            2.0_dp*pi/max(1.0_dp, 1.0_dp*di_n_q_phi_bins), di_n_q_phi_bins)

                        q_bin = 1 + floor(q_mag/di_q_bin_width)

                        num_q_in_bins(q_bin, q_theta_bin, q_phi_bin) = &
                            num_q_in_bins(q_bin, q_theta_bin, q_phi_bin) + 1

                        di(:, q_bin, q_theta_bin, q_phi_bin) = &
                            di(:, q_bin, q_theta_bin, q_phi_bin) + &
                            di_unbinned(:, q1, q2, q3)

                    end if

                end do
            end do
        end do

        ! divide by the number of elements in each bin
        do q1 = 1, di_n_q_bins
            do q2 = 1, di_n_q_theta_bins
                do q3 = 1, di_n_q_phi_bins 

                    di(:, q1, q2, q3) = di(:, q1, q2, q3)/max(1, num_q_in_bins(q1, q2, q3))

                end do
            end do
        end do

    end subroutine 

    subroutine save_dielectric(filename, verbose)
        !! Save the computed dielectric to 'filename'.

        implicit none

        character(len=*) :: filename

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id

        logical, optional :: verbose
        logical :: file_exists

        integer(HSIZE_T) :: dims1(1) = [1]
        ! integer(HSIZE_T) :: dims2(2)
        integer(HSIZE_T) :: dims4(4)

        integer :: error

        if ( verbose ) then

            print*, 'Saving dielectric...'
            print*

        end if

        ! create the file
        call h5open_f(error)

        call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

        call h5fclose_f(file_id, error)
        call h5close_f(error)

        ! open the file to write
        call h5open_f(error)
        call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

        call h5gcreate_f(file_id, 'dielectric', group_id, error)

        dims4 = [di_n_omega_bins, di_n_q_bins, di_n_q_theta_bins, di_n_q_phi_bins]

        call h5ltmake_dataset_double_f(file_id, 'dielectric/dielectric_r', &
            size(dims4), dims4,&
            real(dielec), error)

        call h5ltmake_dataset_double_f(file_id, 'dielectric/dielectric_c', &
            size(dims4), dims4,&
            aimag(dielec), error)

        call h5ltmake_dataset_double_f(file_id, 'dielectric/w_bin_width', &
            size(dims1), dims1,&
            di_omega_bin_width, error)

        call h5ltmake_dataset_double_f(file_id, 'dielectric/q_bin_width', &
            size(dims1), dims1,&
            di_q_bin_width, error)

        call h5ltmake_dataset_double_f(file_id, 'dielectric/q_theta_bin_width', &
            size(dims1), dims1,&
            pi/max(1.0_dp, 1.0_dp*di_n_q_theta_bins), error)

        call h5ltmake_dataset_double_f(file_id, 'dielectric/q_phi_bin_width', &
            size(dims1), dims1,&
            2.0_dp*pi/max(1.0_dp, 1.0_dp*di_n_q_phi_bins), error)

        call h5ltmake_dataset_int_f(file_id, 'dielectric/n_w_bins', &
            size(dims1), dims1,&
            di_n_omega_bins, error)

        call h5ltmake_dataset_int_f(file_id, 'dielectric/n_q_bins', &
            size(dims1), dims1,&
            di_n_q_bins, error)

        call h5ltmake_dataset_int_f(file_id, 'dielectric/n_q_theta_bins', &
            size(dims1), dims1,&
            di_n_q_theta_bins, error)

        call h5ltmake_dataset_int_f(file_id, 'dielectric/n_q_phi_bins', &
            size(dims1), dims1,&
            di_n_q_phi_bins, error)

        call h5fclose_f(file_id, error)
        call h5close_f(error)

    end subroutine

end module
