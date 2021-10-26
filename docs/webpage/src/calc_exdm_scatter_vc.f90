module calc_exdm_scatter_vc
    !! Compute \(i, f, \mathbf{k}_i, \mathbf{k}_f \) contributions to the valence \( \rightarrow \) conduction DM-electron scattering rate.

    use prec
    use constants
    use units
    use math_mod

    use binned_scatter_rate_type
    use bins_scatter_type
    use expt_type
    use material_type
    use dm_model_type
    use PW_dataset_type
    use in_med_scr_type

    use physics_scatter_functions
    use transition_form_factor

    use FFT_util

    implicit none

    interface exdm_scatter_vc_calc
        module procedure exdm_scatter_vc_calc_no_spin 
        module procedure exdm_scatter_vc_calc_spin
    end interface

contains

    subroutine exdm_scatter_vc_calc_spin(binned_rate,&
            FFT_grid, PW_dataset, target_mat, &
            bins, dm_model, expt, in_med_scr, &
            wfc_ik, wfc_fkf, &
            val_id, cond_id, k, kf, verbose)
        !! Compute \(i, f, \mathbf{k}_i, \mathbf{k}_f \) contributions to the valence \( \rightarrow \) conduction 
        !! DM-electron scattering rate from spin dependent wave functions.

        implicit none

        type(binned_scatter_rate_t) :: binned_rate
        type(FFT_grid_t) :: FFT_grid
        type(PW_dataset_t) :: PW_dataset
        type(material_t) :: target_mat
        type(bins_scatter_t) :: bins
        type(dm_model_t) :: dm_model
        type(expt_t) :: expt
        type(in_med_scr_t) :: in_med_scr
        integer :: val_id, cond_id, k, kf
        logical, optional :: verbose

        complex(dp) :: wfc_ik(:, :, :, :)
        complex(dp) :: wfc_fkf(:, :, :, :)
        real(dp) :: f_sq(FFT_grid%n_grid(1), FFT_grid%n_grid(2), FFT_grid%n_grid(3))

        real(dp) :: omega, omega_max

        integer :: m, f, t
        integer :: g1, g2, g3

        real(dp) :: v_m, g_func_val
        real(dp) :: q_vec(3), q_mag, q_min, q_max

        real(dp) :: F_med_sq(dm_model%n_med_FF)

        real(dp) :: scr

        integer :: q_bin, E_bin

        real(dp) :: f_sq_val

        real(dp) :: b_rate(bins%n_q, bins%n_E, dm_model%n_mX, dm_model%n_med_FF, expt%n_time)

        b_rate = 0.0_dp

        f_sq = (0.0_dp, 0.0_dp)

        omega = PW_dataset%energy_bands(kf, cond_id) - PW_dataset%energy_bands(k, val_id)

        if ( ( PW_dataset%energy_bands(kf, cond_id) <= PW_dataset%Ef_max ) &
            .and. ( omega >= expt%E_threshold ) ) then

            E_bin = Q_func(omega, &
                        target_mat%band_gap, &
                        bins%E_width, &
                        bins%n_E)

            call calc_tff_pw_pw(dm_model%tff_id, f_sq, wfc_ik, wfc_fkf, &
                FFT_grid%n_grid, FFT_grid%plan, verbose = .FALSE.)

            ! particle physics
            do m = 1, dm_model%n_mX

                omega_max = 0.5_dp*dm_model%mX(m)*dm_model%vX_max**2

                if ( omega < omega_max ) then

                    q_min = dm_model%mX(m)*dm_model%vX_max*(1.0_dp - sqrt(1.0_dp - (omega/omega_max)))
                    q_max = min(&
                        FFT_grid%q_max, &
                        dm_model%mX(m)*dm_model%vX_max*(1.0_dp + sqrt(1.0_dp - (omega/omega_max))) &
                        )

                    do g3 = 1, FFT_grid%n_grid(3)
                        do g2 = 1, FFT_grid%n_grid(2)
                            do g1 = 1, FFT_grid%n_grid(1)

                                q_vec = PW_dataset%k_grid_xyz(kf, :) &
                                    - PW_dataset%k_grid_xyz(k, :) &
                                    + FFT_grid%sym_G_grid_xyz(g1, g2, g3, :)
                                q_mag = norm2(q_vec)

                                if ( ( q_mag > q_min ) .and. ( q_mag < q_max ) ) then

                                    scr = in_med_scr%screening(q_vec, omega)

                                    q_bin = Q_func(q_mag, 0.0_dp, &
                                                bins%q_width, &
                                                bins%n_q)

                                    do f = 1, dm_model%n_med_FF

                                        F_med_sq(f) = F_med_sq_func(q_mag, &
                                            dm_model%med_FF(f))

                                    end do

                                    f_sq_val = f_sq(g1, g2, g3)

                                    do t = 1, expt%n_time

                                        v_m = v_minus(q_vec, dm_model%mX(m), &
                                            dm_model%vE*expt%vE_direction(t, :), &
                                            omega)

                                        if ( v_m < dm_model%vEsc ) then 

                                            g_func_val = g_func(q_mag, v_m, dm_model)

                                            b_rate(q_bin, E_bin, m, :, t) = &
                                                        b_rate(q_bin, E_bin, m, :, t) + &
                                                        F_med_sq(:)*&
                                                        g_func_val*&
                                                        f_sq_val*&
                                                        scr**(-2)

                                        end if 
                                    end do 

                                end if

                            end do
                        end do
                    end do

                end if
            end do
        end if

        ! multiply overall constants
        do m = 1, dm_model%n_mX

            b_rate(:, :, m, :, :) = (PW_dataset%spin_degen/2.0_dp)*&
                            (2*pi)*(dm_model%rhoX/target_mat%rho_T)*(2*target_mat%pc_vol)**(-2)*&
                            PW_dataset%k_weight(k)*PW_dataset%k_weight(kf)*&
                            red_mass(dm_model%mX(m), m_elec)**(-2)*(dm_model%mX(m))**(-1)*&
                            (1.0_dp*FFT_grid%n)**(-2)*&
                            expt%m_T*expt%exposure*&
                            eV_to_inv_cm**2*&
                            b_rate(:, :, m, :, :)
        end do

        ! add new contribution
        binned_rate%binned_rate = binned_rate%binned_rate + b_rate

    end subroutine

    subroutine exdm_scatter_vc_calc_no_spin(binned_rate,&
            FFT_grid, PW_dataset, target_mat, &
            bins, dm_model, expt, in_med_scr, &
            wfc_ik, wfc_fkf, &
            val_id, cond_id, k, kf, verbose)
        !! Compute \(i, f, \mathbf{k}_i, \mathbf{k}_f \) contributions to the valence \( \rightarrow \) conduction 
        !! DM-electron scattering rate from spin independent wave functions.

        implicit none

        type(binned_scatter_rate_t) :: binned_rate
        type(FFT_grid_t) :: FFT_grid
        type(PW_dataset_t) :: PW_dataset
        type(material_t) :: target_mat
        type(bins_scatter_t) :: bins
        type(dm_model_t) :: dm_model
        type(expt_t) :: expt
        type(in_med_scr_t) :: in_med_scr
        integer :: val_id, cond_id, k, kf
        logical, optional :: verbose

        complex(dp) :: wfc_ik(:, :, :)
        complex(dp) :: wfc_fkf(:, :, :)
        real(dp) :: f_sq(FFT_grid%n_grid(1), FFT_grid%n_grid(2), FFT_grid%n_grid(3))

        real(dp) :: omega, omega_max

        integer :: m, f, t
        integer :: g1, g2, g3

        real(dp) :: v_m, g_func_val
        real(dp) :: q_vec(3), q_mag, q_min, q_max

        real(dp) :: F_med_sq(dm_model%n_med_FF)

        real(dp) :: scr

        integer :: q_bin, E_bin

        real(dp) :: f_sq_val

        real(dp) :: b_rate(bins%n_q, bins%n_E, dm_model%n_mX, dm_model%n_med_FF, expt%n_time)

        b_rate = 0.0_dp

        f_sq = (0.0_dp, 0.0_dp)

        omega = PW_dataset%energy_bands(kf, cond_id) - PW_dataset%energy_bands(k, val_id)

        if ( ( PW_dataset%energy_bands(kf, cond_id) <= PW_dataset%Ef_max ) &
            .and. ( omega >= expt%E_threshold ) ) then

            E_bin = Q_func(omega, &
                        target_mat%band_gap, &
                        bins%E_width, &
                        bins%n_E)

            call calc_tff_pw_pw(dm_model%tff_id, f_sq, wfc_ik, wfc_fkf, &
                FFT_grid%n_grid, FFT_grid%plan, verbose = .FALSE.)

            ! particle physics
            do m = 1, dm_model%n_mX

                omega_max = 0.5_dp*dm_model%mX(m)*dm_model%vX_max**2

                if ( omega < omega_max ) then

                    q_min = dm_model%mX(m)*dm_model%vX_max*(1.0_dp - sqrt(1.0_dp - (omega/omega_max)))
                    q_max = min(&
                        FFT_grid%q_max, &
                        dm_model%mX(m)*dm_model%vX_max*(1.0_dp + sqrt(1.0_dp - (omega/omega_max))) &
                        )

                    do g3 = 1, FFT_grid%n_grid(3)
                        do g2 = 1, FFT_grid%n_grid(2)
                            do g1 = 1, FFT_grid%n_grid(1)

                                q_vec = PW_dataset%k_grid_xyz(kf, :) &
                                    - PW_dataset%k_grid_xyz(k, :) &
                                    + FFT_grid%sym_G_grid_xyz(g1, g2, g3, :)
                                q_mag = norm2(q_vec)

                                if ( ( q_mag > q_min ) .and. ( q_mag < q_max ) ) then

                                    scr = in_med_scr%screening(q_vec, omega)

                                    q_bin = Q_func(q_mag, 0.0_dp, &
                                                bins%q_width, &
                                                bins%n_q)

                                    do f = 1, dm_model%n_med_FF

                                        F_med_sq(f) = F_med_sq_func(q_mag, &
                                            dm_model%med_FF(f))

                                    end do

                                    f_sq_val = f_sq(g1, g2, g3)

                                    do t = 1, expt%n_time

                                        v_m = v_minus(q_vec, dm_model%mX(m), &
                                            dm_model%vE*expt%vE_direction(t, :), &
                                            omega)

                                        if ( v_m < dm_model%vEsc ) then 

                                            g_func_val = g_func(q_mag, v_m, dm_model)

                                            b_rate(q_bin, E_bin, m, :, t) = &
                                                        b_rate(q_bin, E_bin, m, :, t) + &
                                                        F_med_sq(:)*&
                                                        g_func_val*&
                                                        f_sq_val*&
                                                        scr**(-2)

                                        end if 
                                    end do 

                                end if

                            end do
                        end do
                    end do

                end if
            end do
        end if

        ! multiply overall constants
        do m = 1, dm_model%n_mX

            b_rate(:, :, m, :, :) = (PW_dataset%spin_degen/2.0_dp)*&
                            (2*pi)*(dm_model%rhoX/target_mat%rho_T)*(2*target_mat%pc_vol)**(-2)*&
                            PW_dataset%k_weight(k)*PW_dataset%k_weight(kf)*&
                            red_mass(dm_model%mX(m), m_elec)**(-2)*(dm_model%mX(m))**(-1)*&
                            (1.0_dp*FFT_grid%n)**(-2)*&
                            expt%m_T*expt%exposure*&
                            eV_to_inv_cm**2*&
                            b_rate(:, :, m, :, :)
        end do

        ! add new contributions
        binned_rate%binned_rate = binned_rate%binned_rate + b_rate

    end subroutine
    
end module
