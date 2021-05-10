module calc_dme_scatter_vf

    use prec
    use constants
    use math_mod

    use numerics_input
    use material_input
    use particle_physics_scatter
    use DFT_parameters
    use in_med_scr
    use Zeff_input

    implicit none

contains

    subroutine dme_scatter_vf_calc(binned_rate_t,&
            wfc_FT_i, val_id, log_omegas, k_cut, angular_mesh, &
            verbose)

        implicit none

        real(dp) :: binned_rate_t(n_q_bins + 1, n_E_bins + 1, n_mX, n_FDM, n_time)
        real(dp) :: b_rate(n_q_bins + 1, n_E_bins + 1, n_mX, n_FDM, n_time)
        real(dp) :: dRdw(n_q_bins + 1, 2, n_mX, n_FDM, n_time)

        complex(dp) :: wfc_FT_i(n_k, n_in_G)

        integer :: val_id

        real(dp) :: log_omegas(2)

        real(dp) :: omega_kin_max
        real(dp) :: log_omega_bounds(2)
        real(dp) :: log_dRdw(2)

        integer :: k_cut

        real(dp) :: angular_mesh(n_kf_theta*n_kf_phi, 2)

        logical, optional :: verbose

        integer :: i, t, m, f, q, e, g, k, a

        real(dp) :: Zeff

        real(dp) :: kf_theta, kf_phi, kf_mag
        real(dp) :: k_f(3), k_i(3)

        real(dp) :: q_min, q_max, q_mag
        real(dp) :: q_vec(3)
        integer :: q_bin_num

        real(dp) :: omega, Eik, Ef
        real(dp) :: E1, E2
        real(dp) :: omega_1, omega_2

        real(dp) :: fit_params

        real(dp) :: v_m, g_func_val
        real(dp) :: F_med_sq(n_FDM)
        real(dp) :: f_sq
        real(dp) :: scr

        real(dp) :: fermi_val, fermi_factor

        real(dp) :: log_b_rate

        b_rate = 0.0_dp
        dRdw = 0.0_dp

        do m = 1, n_mX

            omega_kin_max = 0.5_dp*mX(m)*v_max**2

            if ( omega_kin_max > 10.0_dp**log_omegas(1) ) then

                log_omega_bounds(1) = log_omegas(1)
                log_omega_bounds(2) = min(log10(omega_kin_max), log_omegas(2))

                ! compute dRdw at the endpoints, then interpolate and integrate between to get 
                ! binned rate
                do i = 1, 2

                    omega = 10.0_dp**log_omega_bounds(i)

                    q_min = mX(m)*v_max*(1.0_dp - sqrt(abs(1.0_dp - (omega/omega_kin_max))))
                    q_max = mX(m)*v_max*(1.0_dp + sqrt(abs(1.0_dp - (omega/omega_kin_max))))

                    do k = 1, k_cut

                        Eik = energy_bands(k, val_id) - maxval(energy_bands(:, :n_val))
                        Ef = omega + Eik

                        if ( ( Ef .gt. 0.0_dp ) .and. ( Ef .gt. Ef_max ) ) then

                            Zeff = get_val_Z_eff(val_id, k) 

                            kf_mag = sqrt(2*m_elec*Ef)
                            
                            fermi_val = 2.0_dp*pi*Zeff*(alpha_EM*m_elec/kf_mag)
                            fermi_factor = fermi_val*(1.0_dp - exp(-fermi_val))**(-1)

                            do a = 1, n_kf_theta*n_kf_phi

                                kf_theta = angular_mesh(a, 1)
                                kf_phi = angular_mesh(a, 2)

                                k_f(1) = kf_mag*sin(kf_theta)*cos(kf_phi)
                                k_f(2) = kf_mag*sin(kf_theta)*sin(kf_phi)
                                k_f(3) = kf_mag*cos(kf_theta)

                                do g = 1, n_in_G

                                    f_sq = wfc_FT_i(k, g)*conjg(wfc_FT_i(k, g))

                                    k_i = k_grid_xyz(k, :) + in_G_grid_xyz(g, :)

                                    q_vec = k_f - k_i
                                    q_mag = norm2(q_vec)

                                    if ( ( q_mag > q_min ) .and. ( q_mag < q_max ) ) then

                                        scr = screening(q_vec, omega)

                                        q_bin_num = Q_func(q_mag, 0.0_dp, q_bin_width, n_q_bins + 1)

                                        do f = 1, n_FDM
                                            F_med_sq(f) = F_med_sq_func(q_mag, FDMPowerList(f))
                                        end do

                                        do t = 1, n_time

                                            v_m = v_minus(q_vec, mX(m), vEVecList(t, :), omega)

                                            if (v_m .lt. vEsc) then 

                                                g_func_val = g_func(q_mag, v_m)

                                                dRdw(q_bin_num, i, m, :, t) = &
                                                            dRdw(q_bin_num, i, m, :, t) + &
                                                            k_weight(k)*&
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

                end do

                ! overall constants
                dRdw(:, :, m, :, :) = (2*pi)*(rhoX/rho_T)*(2*pc_vol)**(-1)*((4*pi)/(n_kf_theta*n_kf_phi*1.0_dp))*&
                                        red_mass(m_elec, mX(m))**(-2)*(mX(m))**(-1)*(2*pi)**(-3)*dRdw(:, :, m, :, :)

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
