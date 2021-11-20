module calc_exdm_scatter_cc_ext
    !! Compute \( i, f, \mathbf{k}_f \) contributions to the core \( \rightarrow \) conduction (extended) DM-electron scattering 
    !! rate.

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
    use core_electron_type
    use in_med_scr_type

    use numerics_scatter_cc_ext

    use physics_scatter_functions
    use transition_form_factor

    implicit none

contains

    subroutine exdm_scatter_cc_ext_calc(binned_rate, &
            core_electron, PW_dataset, target_mat, &
            bins, dm_model, expt, in_med_scr, &
            wfc_FT_fkf, &
            init_id, cond_id, kf, verbose)

        implicit none

        type(binned_scatter_rate_t) :: binned_rate
        type(core_electron_t) :: core_electron
        type(PW_dataset_t) :: PW_dataset
        type(material_t) :: target_mat
        type(bins_scatter_t) :: bins
        type(dm_model_t) :: dm_model
        type(expt_t) :: expt
        type(in_med_scr_t) :: in_med_scr
        type(numerics_scatter_cc_ext_t) :: numerics
        integer :: init_id, cond_id, k, kf
        logical, optional :: verbose

        real(dp) :: dRdq(numerics%n_q, bins%n_E, dm_model%n_mX, dm_model%n_med_FF, expt%n_time)

        complex(dp) :: wfc_FT_fkf(:)

        real(dp) :: f_sq

        real(dp) :: omega, omega_max
        real(dp) :: Ei, Ef

        integer :: m, f, t
        integer :: g1, g2, g3

        real(dp) :: v_m, g_func_val
        real(dp) :: q_vec(3), q_mag, q_min, q_max

        real(dp) :: kf_vec(3), G_vec(3)

        real(dp) :: F_med_sq(dm_model%n_med_FF)

        real(dp) :: scr

        integer :: q_bin, E_bin

        real(dp) :: f_sq_val

        real(dp) :: b_rate(bins%n_q, bins%n_E, dm_model%n_mX, dm_model%n_med_FF, expt%n_time)

        real(dp) :: q_list(numerics%n_q)

        complex(dp) :: form_fac

        complex(dp) :: wfc_i_FT

        real(dp) :: eq_pos(3)

        ! initialize
        b_rate = 0.0_dp
        dRdq = 0.0_dp

        Ei = core_electron%energy(init_id)
        Ef = PW_dataset%energy_bands(kf, cond_id)

        omega = Ef - Ei

        E_bin = Q_func(omega, target_mat%band_gap, bins%E_width, bins%n_E)

        kf_vec = PW_dataset%k_grid_xyz(kf)

        eq_pos = core_electron%eq_pos_xyz(core_electron%config(init_id, 1), :)

        if ( omega >= expt%E_threshold ) then

            do m = 1, dm_model%n_mX

                omega_kin_max = 0.5_dp*dm_model%mX(m)*dm_model%vX_max**2

                q_min = dm_model%mX(m)*dm_model%vX_max*(1.0_dp - sqrt(abs(1.0_dp - (omega/omega_kin_max))))
                q_max = dm_model%mX(m)*dm_model%vX_max*(1.0_dp + sqrt(abs(1.0_dp - (omega/omega_kin_max))))

                q_list = 10.0_dp**uniform_list(numerics%n_q, log10(q_min), log10(q_max))

                do q = 1, numerics%n_q

                    q_mag = q_list(q)

                    do a = 1, size(numerics%q_angular_mesh, 1)

                        q_theta = numerics%q_angular_mesh(a, 1)
                        q_phi = numerics%q_angular_mesh(a, 2)

                        q_vec(1) = q_mag*sin(q_theta)*cos(q_phi)
                        q_vec(2) = q_mag*sin(q_theta)*sin(q_phi)
                        q_vec(3) = q_mag*cos(q_theta)

                        form_fac = (0.0_dp, 0.0_dp)

                        do g = 1, PW_dataset%n_G

                            G_vec = PW_dataset%G_grid_xyz(g, :)

                            wfc_i_FT = core_electron%atomic_sto_wf_FT(init_id, & 
                                G_vec + kf_vec - q_vec)

                            form_fac = form_fac + conjg(wfc_FT_kf)*&
                                                    wfc_i_FT*&
                                                    exp(-ii*dot_product(G_vec, eq_pos))

                        end do

                        f_sq_val = conjg(form_fac)*form_fac

                        scr = in_med_scr%screening(q_vec, omega)

                        do f = 1, dm_model%n_med_FF

                            F_med_sq(f) = F_med_sq_func(q_mag, &
                                dm_model%med_FF(f))

                        end do

                        do t = 1, expt%n_time

                            v_m = v_minus(q_vec, dm_model%mX(m), &
                                dm_model%vE*expt%vE_direction(t, :), &
                                omega)

                            if ( v_m < dm_model%vEsc ) then 

                                g_func_val = g_func(q_mag, v_m, dm_model)

                                dRdq(q, E_bin, m, :, t) = &
                                    dRdq(q, E_bin, m, :, t) + &
                                        q_mag**2*&
                                        scr**(-2)*&
                                        F_med_sq(:)*&
                                        g_func_val*&
                                        f_sq_val

                            end if 
                        end do 

                    end do

                end do

                ! overall constants

                dRdq(:, :, m, :, :) = (2*pi)*PW_dataset%k_weight(kf)*&
                   (dm_model%rhoX/target_mat%rho_T)*(2*target_mat%pc_vol**2)**(-1)*&
                   ((4*pi)/(numerics%n_q_theta*numerics%n_q_phi*1.0_dp))*&
                   (dm_model%mX(m)*red_mass(m_elec, dm_model%mX(m))**2)**(-1)*&
                   (2*pi)**(-3)*dRdq(:, :, m, :, :)

                ! interpolate to get b_rate

                do w = 1, bins%n_E
                    do t = 1, expt%n_time
                        do f = 1, dm_model%n_med_FF

                            call cc_ext_interpolate_dRdq_to_binned_rate(q_list, &
                                dRdq(:, w, m, f, t), &
                                b_rate(:, w, m, f, t), &
                                bins%q_width)

                        end do
                    end do
                end do

            end do

        end if

        binned_rate%binned_rate = binned_rate%binned_rate + b_rate

    end subroutine

    subroutine cc_ext_interpolate_dRdq_to_binned_rate(q_list, dRdq, b_rate, q_bin_width)

        implicit none

        real(dp) :: q_list(:)
        real(dp) :: dRdq(:)
        real(dp) :: b_rate(:)
        real(dp) :: q_bin_width

        ! go through pairs of points
        do q = 1, size(dRdq) - 1

        end do





    end subroutine

end module
