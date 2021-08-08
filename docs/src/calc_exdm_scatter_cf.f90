module calc_exdm_scatter_cf
    !! Compute \( i \) contributions to the core \( \rightarrow \) free DM-electron scattering rate.

    use prec
    use constants
    use units
    use math_mod

    use binned_scatter_rate_type
    use bins_scatter_type
    use expt_type
    use material_type
    use dm_model_type
    use core_electron_type
    use in_med_scr_type

    use physics_scatter_functions
    use transition_form_factor

    use numerics_scatter_cf

    implicit none

contains

    subroutine exdm_scatter_cf_calc(binned_rate,&
            core_electron, target_mat, bins, dm_model, expt, in_med_scr, &
            numerics, &
            init_id, w_id, verbose)
        !! Compute \( i \) contributions to the core \( \rightarrow \) free DM-electron scattering rate.

        implicit none

        type(binned_scatter_rate_t) :: binned_rate
        type(core_electron_t) :: core_electron
        type(material_t) :: target_mat
        type(bins_scatter_t) :: bins
        type(dm_model_t) :: dm_model
        type(expt_t) :: expt
        type(in_med_scr_t) :: in_med_scr
        type(numerics_scatter_cf_t) :: numerics
        integer :: init_id, w_id
        logical, optional :: verbose

        real(dp) :: dRdw(bins%n_q, numerics%n_omega, dm_model%n_mX, dm_model%n_med_FF, expt%n_time)

        real(dp) :: b_rate(bins%n_q, bins%n_E, dm_model%n_mX, dm_model%n_med_FF, expt%n_time)

        real(dp) :: Ei, omega, Ef

        real(dp) :: ki_max, ki_mag, jac_ki, ki_theta, ki_phi
        real(dp) :: ki_vec(3)

        real(dp) :: kf_mag, kf_theta, kf_phi
        real(dp) :: kf_vec(3)

        real(dp) :: omega_kin_max

        real(dp) :: q_min, q_max, q_mag
        real(dp) :: q_vec(3)
        integer :: q_bin

        real(dp) :: v_m, g_func_val
        real(dp) :: F_med_sq(dm_model%n_med_FF)
        real(dp) :: f_sq
        real(dp) :: scr

        real(dp) :: fermi_val, fermi_factor

        integer :: m, w, ki, a, b, t, f, q

        b_rate = 0.0_dp
        dRdw = 0.0_dp

        Ei = core_electron%energy(init_id)

        ki_max = numerics%ki_s*numerics%Zeff(init_id)*alpha_EM*m_elec

        do m = 1, dm_model%n_mX

            omega_kin_max = 0.5_dp*dm_model%mX(m)*dm_model%vX_max**2

            if ( omega_kin_max > numerics%omega_list(w_id) ) then

                do w = 1, 2

                    omega = numerics%omega_list(&
                        min(w_id + (w - 1), numerics%n_omega))

                    q_min = dm_model%mX(m)*dm_model%vX_max*(1.0_dp - sqrt(abs(1.0_dp - (omega/omega_kin_max))))
                    q_max = dm_model%mX(m)*dm_model%vX_max*(1.0_dp + sqrt(abs(1.0_dp - (omega/omega_kin_max))))

                    Ef = Ei + omega

                    if ( ( Ef > 0.0_dp ) &
                        .and. ( omega > expt%E_threshold ) &
                        .and. ( Ef > numerics%Ef_min ) ) then

                        kf_mag = sqrt(2.0_dp*m_elec*Ef)

                        fermi_val = 2.0_dp*pi*numerics%Zeff(init_id)*(alpha_EM*m_elec/kf_mag)
                        fermi_factor = fermi_val*(1.0_dp - exp(-fermi_val))**(-1)

                        do ki = 1, numerics%n_ki

                            ki_mag = (ki_max/numerics%ki_min)**((ki - 1.0_dp)/(numerics%n_ki - 1.0_dp))&
                                *numerics%ki_min
                            jac_ki = (1.0_dp/numerics%n_ki)*ki_mag**3*log10(ki_max/numerics%ki_min)

                            do a = 1, size(numerics%ki_angular_mesh, 1)

                                ki_theta = numerics%ki_angular_mesh(a, 1)
                                ki_phi = numerics%ki_angular_mesh(a, 2)

                                ki_vec(1) = ki_mag*sin(ki_theta)*cos(ki_phi)
                                ki_vec(2) = ki_mag*sin(ki_theta)*sin(ki_phi)
                                ki_vec(3) = ki_mag*cos(ki_theta)

                                f_sq = abs(core_electron%atomic_sto_wf_FT(init_id, ki_vec))**2

                                do b = 1, size(numerics%kf_angular_mesh, 1)

                                    kf_theta = numerics%kf_angular_mesh(b, 1)
                                    kf_phi = numerics%kf_angular_mesh(b, 2)

                                    kf_vec(1) = kf_mag*sin(kf_theta)*cos(kf_phi)
                                    kf_vec(2) = kf_mag*sin(kf_theta)*sin(kf_phi)
                                    kf_vec(3) = kf_mag*cos(kf_theta)

                                    q_vec = kf_vec - ki_vec
                                    q_mag = norm2(q_vec)

                                    if ( ( q_mag > q_min ) .and. ( q_mag < q_max ) ) then

                                        do f = 1, dm_model%n_med_FF
                                            F_med_sq(f) = F_med_sq_func(q_mag, &
                                                dm_model%med_FF(f))
                                        end do

                                        q_bin = Q_func(q_mag, 0.0_dp, &
                                            bins%q_width, bins%n_q)

                                        scr = in_med_scr%screening(q_vec, omega)

                                        do t = 1, expt%n_time

                                            v_m = v_minus(q_vec, dm_model%mX(m), &
                                                dm_model%vE*expt%vE_direction(t, :), &
                                                omega)

                                            if (v_m < dm_model%vEsc) then 

                                                g_func_val = g_func(q_mag, v_m, dm_model)

                                                dRdw(q_bin, w, m, :, t) = dRdw(q_bin, w, m, :, t) + &
                                                    fermi_factor*&
                                                    F_med_sq*&
                                                    g_func_val*&
                                                    f_sq*&
                                                    jac_ki*&
                                                    m_elec*kf_mag*&
                                                    scr**(-2)

                                            end if

                                        end do

                                    end if

                                end do

                            end do

                        end do

                    end if

                end do

                ! add overall constants
                dRdw(:, :, m, :, :) = (core_electron%config(init_id, 5)*pi)*&
                    (dm_model%rhoX/target_mat%rho_T)*&
                    (target_mat%pc_vol)**(-1)*red_mass(dm_model%mX(m), m_elec)**(-2)*&
                    dm_model%mX(m)**(-1)*&
                    (2*pi)**(-6)*(4*pi)/(1.0_dp*numerics%n_ki_theta*numerics%n_ki_phi)*&
                    (4*pi)/(1.0_dp*numerics%n_kf_theta*numerics%n_kf_phi)*&
                    expt%m_T*expt%exposure*&
                    eV_to_inv_cm**2*&
                    dRdw(:, :, m, :, :)

                do q = 1, bins%n_q
                    do t = 1, expt%n_time
                        do f = 1, dm_model%n_med_FF
                            
                            call cf_interpolate_dRdw_to_binned_rate(w_id, &
                                numerics%omega_list, dRdw(q, :, m, f, t), b_rate(q, :, m, f, t), bins%E_width)

                        end do
                    end do
                end do

            end if

        end do

        binned_rate%binned_rate = binned_rate%binned_rate + b_rate

    end subroutine

    subroutine cf_interpolate_dRdw_to_binned_rate(w_id, omega_list, dRdw, b_rate, E_bin_width)
        !! Given \( [\omega_1, \omega_2] \) and \( [ \frac{dR}{d \omega}(\omega_1), \frac{dR}{d\omega}(\omega_2)  ] \),
        !! find the power law fit parameter, \( \beta \), such that 
        !! \( \frac{dR}{d \omega}(\omega_2) = \frac{dR}{d \omega}(\omega_1) \left( \frac{\omega_2}{\omega_1} \right)^\beta \).
        !! Then integrate this function over each bin, \( e \), \( \Delta R_e = \int_{(e - 1)\Delta E}^{e \Delta E} \frac{dR}{d\omega}
        !! d\omega \). Note that \( \frac{dR}{d\omega} \) is assumed to be zero for \( \omega \notin [\omega_1, \omega_2] \).
        !!
        !! Specific to core \( \rightarrow \) free DM-electron scattering rate calculation.

        implicit none

        integer :: w_id
        real(dp) :: omega_list(:)

        real(dp) :: dRdw(:)
        real(dp) :: b_rate(:)
        real(dp) :: E_bin_width

        real(dp) :: omega_bounds(2)
        real(dp) :: dRdw_pair(2)

        real(dp) :: fit_params

        integer :: w, e

        real(dp) :: E1, E2

        if ( w_id < size(omega_list) ) then

            omega_bounds(1) = omega_list(w_id)
            omega_bounds(2) = omega_list(w_id + 1)

            dRdw_pair(1) = dRdw(1)
            dRdw_pair(2) = dRdw(2)

            ! make sure both end points aren't zero
            if ( ( dRdw_pair(1) > 0 ) .and. ( dRdw_pair(2) > 0 ) ) then

                ! find the fit parameters
                fit_params = power_law_fit(&
                    log10(omega_bounds), &
                    log10(dRdw_pair))

                ! now we 'know' dRdw between these two points, lets integrate this function 
                ! over each bin
                do e = 1, size(b_rate)

                    ! edges of the bin
                    E1 = (e - 1)*E_bin_width
                    E2 = e*E_bin_width

                    if ( ( E2 > omega_bounds(1) ) & 
                        .and. ( E1 < omega_bounds(2) ) &
                        .and. ( E1 > 0.0_dp ) ) then

                        ! integrate 
                        b_rate(e) = b_rate(e) + dRdw_pair(1)*&
                            integrate_power_law(fit_params, max(E1, omega_bounds(1)), &
                                            min(E2, omega_bounds(2)), omega_bounds(1))

                    end if

                end do

            end if

        end if

    end subroutine

end module
