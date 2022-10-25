module binned_scatter_rate_RIF_calculator

    implicit none

contains

    subroutine binned_scatter_rate_RIF_compute(binned_rate, &
            omega, q_vec_list, &
            FIF, jac_q_list, &
            jac_I, jac_F, &
            spin_dof, &
            extra_FF, &
            screen_factor_list, &
            exdm_inputs)

        use prec_util, only: dp
        use math_util

        use constants_util

        use exdm_inputs_type

        use scatter_physics_functions

        implicit none

        real(dp) :: binned_rate(:, :, :, :, :)
        real(dp) :: omega
        real(dp) :: q_vec_list(:, :)
        real(dp) :: FIF(:)
        real(dp) :: jac_q_list(:)
        real(dp) :: jac_I, jac_F
        integer :: spin_dof
        real(dp) :: extra_FF
        real(dp) :: screen_factor_list(:)
        type(exdm_inputs_t) :: exdm_inputs

        real(dp) :: b_rate(size(binned_rate, 1),& 
                           size(binned_rate, 2),&
                           size(binned_rate, 3),&
                           size(binned_rate, 4),&
                           size(binned_rate, 5))
            ! Dim : [ n_q, n_E, n_mX, n_med_FF, n_vE ]

        real(dp) :: E_bin_width
        real(dp) :: q_bin_width

        integer :: E_bin

        real(dp) :: q_mag_list(size(q_vec_list, 1))
        real(dp) :: half_q2_mag_list(size(q_vec_list, 1))
        real(dp) :: qm1_mag_list(size(q_vec_list, 1))
        real(dp) :: screen_factor_m2_list(size(q_vec_list, 1))
        real(dp) :: jac_FIF_screen_list(size(q_vec_list, 1))

        integer :: q_bin_list(size(q_vec_list, 1))

        real(dp) :: F_med_sq_list(size(q_vec_list, 1), size(exdm_inputs%dm_model%med_FF))

        real(dp) :: q_vE_list(size(q_vec_list, 1))
        real(dp) :: v_m_list(size(q_vec_list, 1))
        real(dp) :: kinematic_function_list(size(q_vec_list, 1), & 
            size(exdm_inputs%dm_model%mX), &
            size(exdm_inputs%astroph_model%v_e_km_per_sec, 1))

        real(dp) :: b_rate_q(size(binned_rate, 1))
        real(dp) :: norm

        integer :: m, f, v

        b_rate = 0.0_dp

        E_bin_width = exdm_inputs%numerics_binned_scatter_rate%E_bin_width
        ! keV -> eV
        q_bin_width = 1.0e3_dp*exdm_inputs%numerics_binned_scatter_rate%q_bin_width

        E_bin = clamp_int(&
            1 + floor( (omega - exdm_inputs%material%band_gap)/E_bin_width ), &
            1, exdm_inputs%numerics_binned_scatter_rate%n_E_bins)

        ! variables only depending on q
        q_mag_list = norm2(q_vec_list, 2) + 1.0e-8_dp
        half_q2_mag_list = 0.5_dp*q_mag_list**2
        qm1_mag_list = q_mag_list**(-1)
        screen_factor_m2_list = screen_factor_list**(-2)
        jac_FIF_screen_list = jac_q_list*FIF*screen_factor_m2_list 

        q_bin_list = clamp_int(&
            1 + floor( q_mag_list/q_bin_width ), &
            1, exdm_inputs%numerics_binned_scatter_rate%n_q_bins) 

        ! mediator
        do f = 1, size(F_med_sq_list, 2)

            F_med_sq_list(:, f) = (q_mag_list*(alpha_EM*m_elec)**(-1))&
                **(-2*exdm_inputs%dm_model%med_FF(f))

        end do

        ! kinematics
        do v = 1, size(exdm_inputs%astroph_model%v_e_list, 1) 

            q_vE_list = matmul(q_vec_list, exdm_inputs%astroph_model%v_e_list(v, :))

            do m = 1, size(exdm_inputs%dm_model%mX)

                v_m_list = compute_v_minus(q_vE_list, &
                    half_q2_mag_list, qm1_mag_list, &
                    exdm_inputs%dm_model%mX(m)**(-1), &
                    omega, &
                    exdm_inputs%astroph_model%v_esc)

                kinematic_function_list(:, m, v) = &
                    exdm_inputs%astroph_model%kinematic_function(v_m_list, qm1_mag_list)

            end do

        end do
                
        do v = 1, size(kinematic_function_list, 3)
            do f = 1, size(F_med_sq_list, 2)
                do m = 1, size(kinematic_function_list, 2)

                    b_rate_q = 0.0_dp

                    b_rate_q(q_bin_list) = b_rate_q(q_bin_list) + &
                                                jac_FIF_screen_list*&
                                                kinematic_function_list(:, m, v)*&
                                                F_med_sq_list(:, f)

                    b_rate(:, E_bin, m, f, v) = b_rate_q

                end do

            end do
        end do

        do m = 1, size(exdm_inputs%dm_model%mX)

            norm = (2.0_dp*pi/spin_dof)*&
                       (exdm_inputs%dm_model%rho_X/exdm_inputs%material%rho_T)*&
                       exdm_inputs%material%pc_vol**(-2)*&
                       reduced_mass(exdm_inputs%dm_model%mX(m), m_elec)**(-2)*&
                       exdm_inputs%dm_model%mX(m)**(-1)*&
                       exdm_inputs%experiment%M*exdm_inputs%experiment%T*&
                       eV_to_inv_cm**2

            b_rate(:, :, m, :, :) = jac_I*jac_F*extra_FF*&
                                    norm*b_rate(:, :, m, :, :)
        end do

        binned_rate = binned_rate + b_rate

    end subroutine

end module
