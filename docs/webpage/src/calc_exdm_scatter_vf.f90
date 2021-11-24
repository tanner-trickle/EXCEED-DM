module calc_exdm_scatter_vf
    !! Compute \( i \) contributions to the valence \( \rightarrow \) free DM-electron scattering rate.

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

    use numerics_scatter_vf

    implicit none

contains

    subroutine exdm_scatter_vf_calc(binned_rate,&
            PW_dataset, target_mat, &
            bins, dm_model, expt, in_med_scr, &
            numerics, &
            wfc_FT_ik, val_id, k, verbose)
        !! Computes the valence to free scattering rate by first computing \( \frac{dR}{d \omega} \)
        !! for a list of \( \omega \) values, then interpolates and integrates to get the binned rate
        !! in `bins`.

        implicit none

        type(binned_scatter_rate_t) :: binned_rate
        type(PW_dataset_t) :: PW_dataset
        type(material_t) :: target_mat
        type(bins_scatter_t) :: bins
        type(dm_model_t) :: dm_model
        type(expt_t) :: expt
        type(in_med_scr_t) :: in_med_scr
        type(numerics_scatter_vf_t) :: numerics
        complex(dp) :: wfc_FT_ik(:)
        integer :: val_id, k
        logical, optional :: verbose

        real(dp) :: dRdw(bins%n_q, numerics%n_omega, dm_model%n_mX, dm_model%n_med_FF, expt%n_time)

        real(dp) :: b_rate(bins%n_q, bins%n_E, dm_model%n_mX, dm_model%n_med_FF, expt%n_time)

        real(dp) :: kf_theta, kf_phi, kf_mag
        real(dp) :: k_f(3), k_i(3)

        real(dp) :: q_min, q_max, q_mag
        real(dp) :: q_vec(3)
        integer :: q_bin

        real(dp) :: omega, Eik, Ef

        real(dp) :: v_m, g_func_val
        real(dp) :: F_med_sq(dm_model%n_med_FF)
        real(dp) :: f_sq
        real(dp) :: scr

        real(dp) :: fermi_val, fermi_factor

        real(dp) :: omega_kin_max

        integer :: m, w, a, t, q, g, f

        b_rate = 0.0_dp
        dRdw = 1.0e-100_dp

        do m = 1, dm_model%n_mX

            omega_kin_max = 0.5_dp*dm_model%mX(m)*dm_model%vX_max**2

            ! make sure that the kinematically allowed energy deposition is above the lower
            ! threshold energy
            if ( omega_kin_max > numerics%omega_list(1) ) then

                ! log_omega_bounds(1) = log10(numerics%omega_list(1))
                ! log_omega_bounds(2) = min(log10(omega_kin_max), log_omegas(2))

                ! compute dRdw for all w at the endpoints, then interpolate and integrate between to get 
                ! binned rate
                do w = 1, numerics%n_omega

                    omega = numerics%omega_list(w)

                    q_min = dm_model%mX(m)*dm_model%vX_max*(1.0_dp - sqrt(abs(1.0_dp - (omega/omega_kin_max))))
                    q_max = dm_model%mX(m)*dm_model%vX_max*(1.0_dp + sqrt(abs(1.0_dp - (omega/omega_kin_max))))

                    Eik = PW_dataset%energy_bands(k, val_id) - maxval(PW_dataset%energy_bands(:, :PW_dataset%n_val))
                    Ef = omega + Eik

                    if ( ( Ef >= PW_dataset%Ef_max ) .and. ( omega > expt%E_threshold ) ) then

                        kf_mag = sqrt(2*m_elec*Ef)
                            
                        fermi_val = 2.0_dp*pi*numerics%Zeff(val_id, k)*(alpha_EM*m_elec/kf_mag)
                        fermi_factor = fermi_val*(1.0_dp - exp(-fermi_val))**(-1)

                        do a = 1, size(numerics%kf_angular_mesh, 1)

                            kf_theta = numerics%kf_angular_mesh(a, 1)
                            kf_phi = numerics%kf_angular_mesh(a, 2)

                            k_f(1) = kf_mag*sin(kf_theta)*cos(kf_phi)
                            k_f(2) = kf_mag*sin(kf_theta)*sin(kf_phi)
                            k_f(3) = kf_mag*cos(kf_theta)

                            do g = 1, PW_dataset%n_G

                                f_sq = wfc_FT_ik(g)*conjg(wfc_FT_ik(g))

                                k_i = PW_dataset%k_grid_xyz(k, :) + PW_dataset%G_grid_xyz(g, :)

                                q_vec = k_f - k_i
                                q_mag = norm2(q_vec)

                                if ( ( q_mag > q_min ) .and. ( q_mag < q_max ) ) then

                                    scr = in_med_scr%screening(q_vec, omega)

                                    q_bin = Q_func(q_mag, 0.0_dp, bins%q_width, bins%n_q)

                                    do f = 1, dm_model%n_med_FF
                                        F_med_sq(f) = F_med_sq_func(q_mag, dm_model%med_FF(f))
                                    end do

                                    do t = 1, expt%n_time

                                        v_m = v_minus(q_vec, dm_model%mX(m), dm_model%vE*expt%vE_direction(t, :), omega)

                                        if (v_m < dm_model%vEsc) then 

                                            g_func_val = g_func(q_mag, v_m, dm_model)

                                            dRdw(q_bin, w, m, :, t) = &
                                                        dRdw(q_bin, w, m, :, t) + &
                                                        f_sq*&
                                                        fermi_factor*&
                                                        g_func_val*F_med_sq*&
                                                        2.0_dp*m_elec*kf_mag*&
                                                        scr**(-2)

                                        end if

                                    end do

                                end if

                            end do

                        end do

                    end if

                end do

                ! overall constants
                dRdw(:, :, m, :, :) = (2*pi)*PW_dataset%k_weight(k)*&
                                        (dm_model%rhoX/target_mat%rho_T)*(2*target_mat%pc_vol)**(-1)*&
                                        ((4*pi)/(numerics%n_kf_theta*numerics%n_kf_phi*1.0_dp))*&
                                        red_mass(m_elec, dm_model%mX(m))**(-2)*(dm_model%mX(m))**(-1)*&
                                        (2*pi)**(-3)*&
                                        expt%m_T*expt%exposure*&
                                        eV_to_inv_cm**2*&
                                        dRdw(:, :, m, :, :)

                do q = 1, bins%n_q
                    do t = 1, expt%n_time
                        do f = 1, dm_model%n_med_FF

                            call vf_interpolate_dRdw_to_binned_rate(numerics%omega_list, &
                                dRdw(q, :, m, f, t), &
                                b_rate(q, :, m, f, t), &
                                bins%E_width, target_mat%band_gap)

                        end do
                    end do
                end do

            end if

        end do

        binned_rate%binned_rate = binned_rate%binned_rate + b_rate

    end subroutine

    subroutine vf_interpolate_dRdw_to_binned_rate(omega_list, dRdw, b_rate, E_bin_width, band_gap)
        !! Given a list of \( \omega_i \) values, and \( \frac{dR}{d \omega} \) at those values, go through adjacent pairs,
        !! \( [\omega_{i}, \omega_{i + 1}] \) and \( [ \frac{dR}{d \omega}(\omega_{i}), \frac{dR}{d\omega}(\omega_{i + 1})  ] \),
        !! find the power law fit parameter, \( \beta \), such that 
        !! \( \frac{dR}{d \omega}(\omega_{i + 1}) = \frac{dR}{d \omega}(\omega_i) \left( \frac{\omega_{i + 1}}{\omega_i} \right)^\beta \).
        !! Then integrate this function over each bin, \( e \), \( \Delta R_e = \int_{(e - 1)\Delta E}^{e \Delta E} \frac{dR}{d\omega}
        !! d\omega \) for each \( i \). Note that \( \frac{dR}{d\omega} \) is assumed to be zero for \( \omega \notin \omega_i \).
        !!
        !! Specific to valence \( \rightarrow \) free DM-electron scattering rate calculation.

        implicit none

        real(dp) :: omega_list(:)

        real(dp) :: dRdw(:)
        real(dp) :: b_rate(:)
        real(dp) :: E_bin_width

        real(dp) :: omega_bounds(2)
        real(dp) :: dRdw_pair(2)

        real(dp) :: fit_params

        integer :: w, e

        real(dp) :: E1, E2

        real(dp) :: E_bin_L, E_bin_R

        real(dp) :: band_gap

        if ( size(b_rate) > 1 ) then

            do e = 1, size(b_rate) - 1

                E_bin_L = band_gap + (e - 1)*E_bin_width
                E_bin_R = band_gap + e*E_bin_width

                b_rate(e) = b_rate(e) + integrate_log_interpolate(&
                    omega_list,& 
                    dRdw, E_bin_L, E_bin_R)


            end do

        end if

        ! last bin integrates over the rest of the list
        b_rate(size(b_rate)) = b_rate(size(b_rate)) + integrate_log_interpolate(&
                                omega_list, dRdw,&
                                band_gap + (size(b_rate) - 1)*E_bin_width,&
                                omega_list(size(omega_list)))

    end subroutine

end module
