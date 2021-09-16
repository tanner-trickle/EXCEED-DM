module calc_dielectric_vc
    !! Compute contributions to the dielectric, \( \varepsilon \), from valence \( \rightarrow \) conduction
    !! transitions.

    use prec
    use constants
    use units
    use math_mod

    use width_parameters_type
    use bins_dielectric_type
    use expt_type
    use material_type
    use dm_model_type
    use PW_dataset_type

    use numerics_dielectric

    use transition_form_factor

    use FFT_util

    implicit none

    interface dielectric_calc_vc
        module procedure dielectric_calc_vc_no_spin
        module procedure dielectric_calc_vc_spin
    end interface

contains

    subroutine dielectric_calc_vc_spin(di, &
            FFT_grid, PW_dataset, target_mat, bins, &
            widths, numerics, &
            wfc_ik, wfc_fkf, &
            val_id, cond_id, k, kf, &
            verbose)
        !! Compute contributions to the dielectric, \( \varepsilon \), from valence \( \rightarrow \) conduction
        !! transitions with spin dependent wave functions.
        
        implicit none

        complex(dp) :: di(:, :, :, :)
        type(FFT_grid_t) :: FFT_grid
        type(PW_dataset_t) :: PW_dataset
        type(material_t) :: target_mat
        type(bins_dielectric_t) :: bins
        type(width_parameters_t) :: widths
        type(numerics_dielectric_t) :: numerics
        complex(dp) :: wfc_ik(:, :, :, :)
        complex(dp) :: wfc_fkf(:, :, :, :)
        integer :: val_id
        integer :: cond_id
        integer :: k
        integer :: kf
        logical, optional :: verbose

        real(dp) :: omega
        real(dp) :: q_vec(3)

        integer :: g1, g2, g3

        real(dp) :: f_sq(FFT_grid%n_grid(1), FFT_grid%n_grid(2), FFT_grid%n_grid(3))

        integer :: w
        complex(dp) :: elec_props(bins%n_E)

        real(dp) :: d_omega

        real(dp) :: delta

        integer :: q_bin, q_theta_bin, q_phi_bin

        real(dp) :: q_hat(3)
        real(dp) :: q_mag, q_phi, q_theta

        real(dp) :: Ei, Ef

        integer :: q1, q2, q3

        real(dp) :: q_red(3)

        integer :: q_ind(3)

        real(dp) :: ddi(bins%n_E, bins%n_q, bins%n_q_theta, bins%n_q_phi)

        ddi = (0.0_dp, 0.0_dp)

        call calc_tff_pw_pw([1, 1], f_sq, wfc_ik, wfc_fkf, &
            FFT_grid%n_grid, FFT_grid%plan, verbose = .FALSE.)

        Ei = PW_dataset%energy_bands(k, val_id)
        Ef = PW_dataset%energy_bands(kf, cond_id)

        do w = 1, bins%n_E

            omega = bins%E_width*(w - 0.5_dp)

            d_omega = Ef - Ei

            ! currently, only compute with one \( \delta \) parameterization at a time.
            delta = widths%get_width(1, omega)

            elec_props(w) = ( omega - d_omega + ii*delta )**(-1) - &
                ( omega + d_omega - ii*delta )**(-1)

        end do

        do g3 = 1, FFT_grid%n_grid(3)
            do g2 = 1, FFT_grid%n_grid(2)
                do g1 = 1, FFT_grid%n_grid(1)

                    q_vec = PW_dataset%k_grid_xyz(kf, :) &
                        - PW_dataset%k_grid_xyz(k, :) + &
                        FFT_grid%sym_G_grid_xyz(g1, g2, g3, :)

                    q_red = PW_dataset%k_grid_red(kf, :) &
                        - PW_dataset%k_grid_red(k, :) + &
                        FFT_grid%sym_G_grid_red(g1, g2, g3, :)

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

                        ddi(:, q_bin, q_theta_bin, q_phi_bin) = &
                            ddi(:, q_bin, q_theta_bin, q_phi_bin) + &
                                (-1.0_dp)*(PW_dataset%spin_degen/2.0_dp)*&
                                (e_EM**2/q_mag**2)*(target_mat%pc_vol)**(-1)*&
                                PW_dataset%k_weight(k)*(1.0_dp*FFT_grid%n)**(-2)*&
                                elec_props(:)*f_sq(g1, g2, g3)

                    end if

                end do
            end do
        end do

        ! divide by the number of elements in each bin
        do q1 = 1, bins%n_q
            do q2 = 1, bins%n_q_theta
                do q3 = 1, bins%n_q_phi

                    ddi(:, q1, q2, q3) = ddi(:, q1, q2, q3)/max(1, numerics%n_q_bin(q1, q2, q3))

                end do
            end do
        end do

        di = di + ddi

    end subroutine 

    subroutine dielectric_calc_vc_no_spin(di, &
            FFT_grid, PW_dataset, target_mat, bins, &
            widths, numerics, &
            wfc_ik, wfc_fkf, &
            val_id, cond_id, k, kf, &
            verbose)
        !! Compute contributions to the dielectric, \( \varepsilon \), from valence \( \rightarrow \) conduction
        !! transitions with spin independent wave functions.
        
        implicit none

        complex(dp) :: di(:, :, :, :)
        type(FFT_grid_t) :: FFT_grid
        type(PW_dataset_t) :: PW_dataset
        type(material_t) :: target_mat
        type(bins_dielectric_t) :: bins
        type(width_parameters_t) :: widths
        type(numerics_dielectric_t) :: numerics
        complex(dp) :: wfc_ik(:, :, :)
        complex(dp) :: wfc_fkf(:, :, :)
        integer :: val_id
        integer :: cond_id
        integer :: k
        integer :: kf
        logical, optional :: verbose

        complex(dp) :: di_unbinned(bins%n_E, &
            numerics%n_q_grid(1), numerics%n_q_grid(2), numerics%n_q_grid(3))

        real(dp) :: omega
        real(dp) :: q_vec(3)

        integer :: g1, g2, g3

        real(dp) :: f_sq(FFT_grid%n_grid(1), FFT_grid%n_grid(2), FFT_grid%n_grid(3))

        integer :: w
        complex(dp) :: elec_props(bins%n_E)

        real(dp) :: d_omega

        real(dp) :: delta

        integer :: q_bin, q_theta_bin, q_phi_bin

        real(dp) :: q_hat(3)
        real(dp) :: q_mag, q_phi, q_theta

        real(dp) :: Ei, Ef

        integer :: num_q_in_bins(bins%n_q, bins%n_q_theta, bins%n_q_phi)

        integer :: q1, q2, q3

        real(dp) :: q_red(3)

        integer :: q_ind(3)

        complex(dp) :: ddi(bins%n_E, bins%n_q, bins%n_q_theta, &
            bins%n_q_phi)

        ddi = (0.0_dp, 0.0_dp)
        di_unbinned = (0.0_dp, 0.0_dp)

        num_q_in_bins = 0

        call calc_tff_pw_pw([1, 1], f_sq, wfc_ik, wfc_fkf, &
            FFT_grid%n_grid, FFT_grid%plan, verbose = .FALSE.)

        Ei = PW_dataset%energy_bands(k, val_id)
        Ef = PW_dataset%energy_bands(kf, cond_id)

        do w = 1, bins%n_E

            omega = bins%E_width*(w - 0.5_dp)

            d_omega = Ef - Ei

            ! currently, only compute with one \( \delta \) parameterization at a time.
            delta = widths%get_width(1, omega)

            elec_props(w) = ( omega - d_omega + ii*delta )**(-1) - &
                ( omega + d_omega - ii*delta )**(-1)

        end do

        do g3 = 1, FFT_grid%n_grid(3)
            do g2 = 1, FFT_grid%n_grid(2)
                do g1 = 1, FFT_grid%n_grid(1)

                    q_vec = PW_dataset%k_grid_xyz(kf, :) &
                        - PW_dataset%k_grid_xyz(k, :) + &
                        FFT_grid%sym_G_grid_xyz(g1, g2, g3, :)

                    q_red = PW_dataset%k_grid_red(kf, :) &
                        - PW_dataset%k_grid_red(k, :) + &
                        FFT_grid%sym_G_grid_red(g1, g2, g3, :)

                    q_mag = norm2(q_vec)

                    if ( ( q_mag > 1.0e-8_dp ) .and. &
                        ( q_mag < bins%n_q*bins%q_width ) ) then

                        q_ind(1) = 1 + int(numerics%n_k_vec(1)*(q_red(1) - numerics%q_grid_min(1)))
                        q_ind(2) = 1 + int(numerics%n_k_vec(2)*(q_red(2) - numerics%q_grid_min(2)))
                        q_ind(3) = 1 + int(numerics%n_k_vec(3)*(q_red(3) - numerics%q_grid_min(3)))

                        di_unbinned(:, q_ind(1), q_ind(2), q_ind(3)) = &
                            di_unbinned(:, q_ind(1), q_ind(2), q_ind(3)) + &
                            (-1.0_dp)*(PW_dataset%spin_degen/2.0_dp)*&
                            (e_EM**2/q_mag**2)*(target_mat%pc_vol)**(-1)*&
                            PW_dataset%k_weight(k)*(1.0_dp*FFT_grid%n)**(-2)*&
                            elec_props(:)*f_sq(g1, g2, g3)

                    end if

                end do
            end do
        end do

        ! now bin the unbinned dielectric
        do q3 = 1, numerics%n_q_grid(3)
            do q2 = 1, numerics%n_q_grid(2)
                do q1 = 1, numerics%n_q_grid(1) 

                    q_red(1) = (1.0_dp/numerics%n_k_vec(1))*( q1 - 1 ) + numerics%q_grid_min(1)
                    q_red(2) = (1.0_dp/numerics%n_k_vec(2))*( q2 - 1 ) + numerics%q_grid_min(2)
                    q_red(3) = (1.0_dp/numerics%n_k_vec(3))*( q3 - 1 ) + numerics%q_grid_min(3)

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

                        num_q_in_bins(q_bin, q_theta_bin, q_phi_bin) = &
                            num_q_in_bins(q_bin, q_theta_bin, q_phi_bin) + 1

                        ddi(:, q_bin, q_theta_bin, q_phi_bin) = &
                            ddi(:, q_bin, q_theta_bin, q_phi_bin) + &
                            di_unbinned(:, q1, q2, q3)

                    end if

                end do
            end do
        end do

        ! divide by the number of elements in each bin
        do q1 = 1, bins%n_q
            do q2 = 1, bins%n_q_theta
                do q3 = 1, bins%n_q_phi

                    ddi(:, q1, q2, q3) = ddi(:, q1, q2, q3)/max(1, num_q_in_bins(q1, q2, q3))

                end do
            end do
        end do

        di = di + ddi

    end subroutine 

end module
