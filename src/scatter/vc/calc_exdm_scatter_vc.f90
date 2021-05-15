module calc_exdm_scatter_vc
    !! Holds the subroutine which calculates the binned rate for the vc calculation

    use prec
    use constants
    use math_mod

    use numerics_input
    use material_input
    use particle_physics_scatter
    use DFT_parameters
    use FFT_util
    use in_med_scr

    implicit none

contains

    subroutine dme_scatter_vc_calc(binned_rate_t,&
            wfc_FT_i, wfc_FT_f,&
            val_id, cond_id, &
            n_FFT_grid,&
            k_cut, &
            wfc_FFT_plan, Tif_FFT_plan,&
            verbose)

        implicit none

        real(dp) :: binned_rate_t(n_q_bins + 1, n_E_bins + 1, n_mX, n_FDM, n_time)
        real(dp) :: b_rate(n_q_bins + 1, n_E_bins + 1, n_mX, n_FDM, n_time)

        integer :: n_FFT_grid(3)

        complex(dp) :: wfc_FT_i(n_k, n_in_G)
        complex(dp) :: wfc_FT_f(n_k, n_in_G)

        complex(dp) :: wfc_i(n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3))
        complex(dp) :: wfc_f(n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3))
        complex(dp) :: wfc_FT_i_exp(n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3))
        complex(dp) :: wfc_FT_f_exp(n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3))

        complex(dp) :: T_if(n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3))
        complex(dp) :: f_G(n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3))

        real(dp) :: f_sq(n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3))

        integer :: k_cut
        integer :: val_id, cond_id

        integer :: wfc_fft_plan(8)
        integer :: Tif_fft_plan(8)

        logical, optional :: verbose

        integer :: ki, kf
        integer :: m, f, t
        integer :: g1, g2, g3

        real(dp) :: v_m, g_func_val
        real(dp) :: q_vec(3), q_mag, q_min, q_max, q_max_FFT

        real(dp) :: omega, omega_max, Ei, Ef
        real(dp) :: F_med_sq(n_FDM)

        real(dp) :: scr

        integer :: q_bin_num, E_bin_num

        integer :: n_FFT

        real(dp) :: f_sq_val

        b_rate = 0.0_dp
        wfc_FT_i_exp = (0.0_dp, 0.0_dp)
        wfc_FT_f_exp = (0.0_dp, 0.0_dp)
        wfc_i = (0.0_dp, 0.0_dp)
        wfc_f = (0.0_dp, 0.0_dp)
        T_if = (0.0_dp, 0.0_dp)
        f_G = (0.0_dp, 0.0_dp)
        f_sq = (0.0_dp, 0.0_dp)

        n_FFT = n_FFT_grid(1)*n_FFT_grid(2)*n_FFT_grid(3)

        q_max_FFT = n_FFT_grid(1)*q_s_FFT

        do ki = 1, k_cut

            call expand_wfc_FT_for_FFT(n_FFT_grid, wfc_FT_i(ki, :), wfc_FT_i_exp,&
                verbose = verbose)

            call dfftw_execute_dft(wfc_fft_plan, wfc_FT_i_exp, wfc_i) 

            Ei = energy_bands(ki, val_id)

            do kf = 1, k_cut
                
                Ef = energy_bands(kf, cond_id)

                omega = Ef - Ei

                E_bin_num = Q_func(omega, band_gap, E_bin_width, n_E_bins + 1)

                if ( ( E_bin_num >= E_bin_threshold ) .and. &
                     ( Ef <= Ef_max ) ) then

                    call expand_wfc_FT_for_FFT(n_FFT_grid, wfc_FT_f(kf, :), wfc_FT_f_exp,&
                        verbose = verbose)

                    call dfftw_execute_dft(wfc_fft_plan, wfc_FT_f_exp, wfc_f) 

                    T_if = conjg(wfc_f)*wfc_i

                    f_G = (0.0_dp, 0.0_dp)
                    call dfftw_execute_dft(Tif_fft_plan, T_if, f_G) 

                    f_sq = abs(f_G)**2

                    ! particle physics
                    do m = 1, n_mX

                        omega_max = 0.5_dp*mX(m)*v_max**2

                        if ( omega < omega_max ) then

                            q_min = mX(m)*v_max*(1.0_dp - sqrt(1.0_dp - (omega/omega_max)))
                            q_max = min(&
                                q_max_FFT, &
                                mX(m)*v_max*(1.0_dp + sqrt(1.0_dp - (omega/omega_max))) &
                                )

                            do g3 = 1, n_FFT_grid(3)
                                do g2 = 1, n_FFT_grid(2)
                                    do g1 = 1, n_FFT_grid(1)

                                        q_vec = k_grid_xyz(kf, :) - k_grid_xyz(ki, :) + sym_FFT_G_grid_xyz(g1, g2, g3, :)
                                        q_mag = norm2(q_vec)

                                        if ( ( q_mag > q_min ) .and. ( q_mag < q_max ) ) then

                                            scr = screening(q_vec, omega)

                                            q_bin_num = Q_func(q_mag, 0.0_dp, q_bin_width, n_q_bins + 1)

                                            do f = 1, n_FDM
                                                F_med_sq(f) = F_med_sq_func(q_mag, FDMPowerList(f))
                                            end do

                                            f_sq_val = f_sq(g1, g2, g3)

                                            do t = 1, n_time

                                                v_m = v_minus(q_vec, mX(m), vEVecList(t, :), omega)

                                                if ( v_m < vEsc ) then 

                                                    g_func_val = g_func(q_mag, v_m)

                                                    b_rate(q_bin_num, E_bin_num, m, :, t) = &
                                                                        b_rate(q_bin_num, E_bin_num, m, :, t) + &
                                                                        k_weight(ki)*k_weight(kf)*&
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
            end do
        end do

        ! add overall constants
        do m = 1, n_mX

            b_rate(:, :, m, :, :) = (2*pi)*(rhoX/rho_T)*(2*pc_vol)**(-2)*&
                            red_mass(mX(m), m_elec)**(-2)*(mX(m))**(-1)*&
                            (1.0_dp*n_FFT)**(-2)*b_rate(:, :, m, :, :)
        end do

        binned_rate_t = b_rate

    end subroutine
    
end module
