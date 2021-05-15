module calc_exdm_scatter_cf
    !! Holds the subroutine which calculates the binned rate for the cf calculation

    use prec
    use constants
    use math_mod

    use numerics_input
    use material_input
    use particle_physics_scatter
    use in_med_scr
    use Zeff_input
    use core_electrons

    implicit none

contains

    subroutine dme_scatter_cf_calc(binned_rate_t, init_id, log_omegas, &
            ki_cut, ki_angular_mesh, kf_angular_mesh, verbose)

        implicit none

        logical, optional :: verbose

        real(dp) :: binned_rate_t(n_q_bins + 1, n_E_bins + 1, n_mX, n_FDM, n_time)
        real(dp) :: b_rate(n_q_bins + 1, n_E_bins + 1, n_mX, n_FDM, n_time)
        real(dp) :: dRdw(n_q_bins + 1, 2, n_mX, n_FDM, n_time)

        complex(dp) :: wfc_FT_i(n_k, n_in_G)

        integer :: init_id

        real(dp) :: log_omegas(2)
        integer :: ki_cut

        real(dp) :: ki_angular_mesh(:, :)
        real(dp) :: kf_angular_mesh(:, :)

        real(dp) :: Ei, Ef, omega, log_omega
        real(dp) :: omega_1, omega_2, E1, E2

        real(dp) :: log_omega_bounds(2)

        real(dp) :: ki_vec(3), ki_mag, ki_theta, ki_phi
        real(dp) :: ki_max
        real(dp) :: jac_ki
        real(dp) :: kf_vec(3), kf_mag, kf_theta, kf_phi

        real(dp) :: q_vec(3), q_mag, q_min, q_max

        integer :: q_bin_num

        real(dp) :: v_m, g_func_val

        real(dp) :: Zeff
        real(dp) :: fermi_val, fermi_factor

        real(dp) :: scr
        real(dp) :: F_med_sq(n_FDM)

        real(dp) :: f_sq

        real(dp) :: omega_kin_max
        real(dp) :: log_dRdw(2)

        real(dp) :: fit_params

        integer :: i, a, b, t, ki, kf, m, q, f, e

        b_rate = 0.0_dp
        dRdw = 0.0_dp

        Ei = core_energy(init_id)

        Zeff = get_val_Z_eff(init_id, 1)

        ki_max = ki_s*Zeff*alpha_EM*m_elec

        do m = 1, n_mX

            omega_kin_max = 0.5_dp*mX(m)*v_max**2

            if ( omega_kin_max > 10.0_dp**log_omegas(1) ) then

                log_omega_bounds(1) = log_omegas(1)
                log_omega_bounds(2) = min(log10(omega_kin_max), log_omegas(2))

                do i = 1, 2

                    omega = 10.0_dp**log_omega_bounds(i)

                    q_min = mX(m)*v_max*(1.0_dp - sqrt(abs(1.0_dp - (omega/omega_kin_max))))
                    q_max = mX(m)*v_max*(1.0_dp + sqrt(abs(1.0_dp - (omega/omega_kin_max))))

                    Ef = Ei + omega

                    if ( ( Ef > 0.0_dp ) .and. ( Ef > Ef_max ) ) then

                        kf_mag = sqrt(2.0_dp*m_elec*Ef)

                        fermi_val = 2.0_dp*pi*Zeff*(alpha_EM*m_elec/kf_mag)
                        fermi_factor = fermi_val*(1.0_dp - exp(-fermi_val))**(-1)

                        do ki = 1 + n_ki - ki_cut, n_ki

                            ki_mag = (ki_max/ki_min)**((ki - 1.0_dp)/(n_ki - 1.0_dp))*ki_min
                            jac_ki = (1.0_dp/n_ki)*ki_mag**3*log10(ki_max/ki_min)

                            do a = 1, size(ki_angular_mesh, 1)

                                ki_theta = ki_angular_mesh(a, 1)
                                ki_phi = ki_angular_mesh(a, 2)

                                ki_vec(1) = ki_mag*sin(ki_theta)*cos(ki_phi)
                                ki_vec(2) = ki_mag*sin(ki_theta)*sin(ki_phi)
                                ki_vec(3) = ki_mag*cos(ki_theta)

                                f_sq = abs(core_sto_wf_FT(init_id, ki_vec))**2

                                do b = 1, size(kf_angular_mesh, 1)

                                    kf_theta = kf_angular_mesh(b, 1)
                                    kf_phi = kf_angular_mesh(b, 2)

                                    kf_vec(1) = kf_mag*sin(kf_theta)*cos(kf_phi)
                                    kf_vec(2) = kf_mag*sin(kf_theta)*sin(kf_phi)
                                    kf_vec(3) = kf_mag*cos(kf_theta)

                                    q_vec = kf_vec - ki_vec
                                    q_mag = norm2(q_vec)

                                    if ( ( q_mag > q_min ) .and. ( q_mag < q_max ) ) then

                                        do f = 1, n_FDM
                                            F_med_sq(f) = F_med_sq_func(q_mag, FDMPowerList(f))
                                        end do

                                        q_bin_num = Q_func(q_mag, 0.0_dp, q_bin_width, n_q_bins + 1)

                                        scr = screening(q_vec, omega)

                                        do t = 1, n_time

                                            v_m = v_minus(q_vec, mX(m), vEVecList(t, :), omega)

                                            if (v_m .lt. vEsc) then 

                                                g_func_val = g_func(q_mag, v_m)

                                                dRdw(q_bin_num, i, m, :, t) = dRdw(q_bin_num, i, m, :, t) + &
                                                    fermi_factor*&
                                                    F_med_sq*&
                                                    g_func_val*&
                                                    f_sq*&
                                                    jac_ki*&
                                                    m_elec*kf_mag*scr**(-2)

                                            end if

                                        end do

                                    end if

                                end do

                            end do

                        end do

                    end if

                end do

                ! add overall constants
                dRdw(:, :, m, :, :) = (core_elec_conf(i, 5)*pi)*(rhoX/rho_T)*&
                    (pc_vol)**(-1)*red_mass(mX(m), m_elec)**(-2)*mX(m)**(-1)*&
                    (2*pi)**(-6)*(4*pi)/(1.0_dp*n_ki_theta*n_ki_phi)*&
                    (4*pi)/(1.0_dp*n_kf_theta*n_kf_phi)*&
                    dRdw(:, :, m, :, :)

                ! now that dRdw has been computed at the end points we need to interpolate and integrate
                ! to get the binned rate

                ! our interpolating function is dRdw = dRdw(1)*(w/w_1)**b

                ! This procedure is equivalent to saving dRdw and interpolating later

                do q = 1, n_q_bins + 1
                    do t = 1, n_time
                        do f = 1, n_FDM

                            if ( ( dRdw(q, 1, m, f, t) > 0.0_dp ) .and. ( dRdw(q, 2, m, f, t) > 0.0_dp ) ) then

                                log_dRdw = log10(dRdw(q, :, m, f, t))

                                fit_params = power_law_fit(log_omega_bounds, log_dRdw)

                                omega_1 = 10.0_dp**log_omega_bounds(1)
                                omega_2 = 10.0_dp**log_omega_bounds(2)

                                do e = E_bin_threshold, n_E_bins + 1

                                    E1 = (e - 1)*E_bin_width
                                    E2 = e*E_bin_width

                                    if ( ( E2 > omega_1 ) & 
                                        .and. ( E1 < omega_2 ) &
                                        .and. ( E1 > 0.0_dp ) ) then

                                        b_rate(q, e, m, f, t) = 10.0_dp**log_dRdw(1)*&
                                            integrate_power_law(fit_params, max(E1, omega_1), &
                                                                min(E2, omega_2), omega_1)

                                    end if

                                end do

                            end if

                        end do
                    end do
                end do

            end if

        end do

        binned_rate_t = b_rate

    end subroutine

end module
