module binned_scatter_rate_RIF_calculator

    implicit none

contains

    subroutine binned_scatter_rate_RIF_compute(binned_rate, &
            omega, q_vec_list, &
            FIF, jac_q_list, &
            jac_I, jac_F, &
            spin_dof, &
            extra_FF, &
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

        type(exdm_inputs_t) :: exdm_inputs

        integer :: E_bin
        real(dp) :: E_bin_width
        real(dp) :: q_bin_width

        real(dp) :: q_mag_list(size(q_vec_list, 1))
        real(dp) :: screen_factor_list(size(q_vec_list, 1))
        real(dp) :: v_m_list(size(q_vec_list, 1))
        logical :: kinematic_mask(size(q_vec_list, 1))

        real(dp), allocatable :: v_m_packed(:)
        real(dp), allocatable :: q_mag_packed(:)
        real(dp), allocatable :: F_med_sq_packed(:)
        real(dp), allocatable :: FIF_packed(:)
        real(dp), allocatable :: jac_q_packed(:)
        integer, allocatable :: q_bin_packed(:)
        real(dp), allocatable :: kinematic_function_packed(:)
        real(dp), allocatable :: screen_factor_list_packed(:)

        real(dp) :: b_rate(size(binned_rate, 1),& 
                           size(binned_rate, 2),&
                           size(binned_rate, 3),&
                           size(binned_rate, 4),&
                           size(binned_rate, 5))
            !* Dim : [ n_q, n_E, n_mX, n_med_FF, n_vE ]
        real(dp) :: b_rate_q(size(binned_rate, 1))

        real(dp) :: norm

        integer :: m, f, v

        b_rate = 0.0_dp

        E_bin_width = exdm_inputs%numerics_binned_scatter_rate%E_bin_width
        q_bin_width = 1.0e3_dp*exdm_inputs%numerics_binned_scatter_rate%q_bin_width

        E_bin = clamp_int(&
            1 + floor( (omega - exdm_inputs%material%band_gap)/E_bin_width ), &
            1, exdm_inputs%numerics_binned_scatter_rate%n_E_bins)

        q_mag_list = norm2(q_vec_list, 2) + 1.0e-8_dp

        if ( trim(adjustl(exdm_inputs%screening%type)) /= '' ) then
            screen_factor_list = exdm_inputs%screening%screening_factor(q_vec_list, omega)
        else 
            screen_factor_list = 1.0_dp
        end if

        do m = 1, size(exdm_inputs%dm_model%mX)
            do v = 1, size(exdm_inputs%astroph_model%v_e_list, 1)

                v_m_list = compute_v_minus(q_vec_list, q_mag_list, &
                    exdm_inputs%dm_model%mX(m), &
                    exdm_inputs%astroph_model%v_e_list(v, :), omega)

                kinematic_mask = ( v_m_list < exdm_inputs%astroph_model%v_esc )

                v_m_packed = pack(v_m_list, kinematic_mask)
                q_mag_packed = pack(q_mag_list, kinematic_mask)
                FIF_packed = pack(FIF, kinematic_mask) 
                jac_q_packed = pack(jac_q_list, kinematic_mask)
                screen_factor_list_packed = pack(screen_factor_list, kinematic_mask)

                q_bin_packed = clamp_int(&
                    1 + floor( q_mag_packed/q_bin_width ), &
                    1, exdm_inputs%numerics_binned_scatter_rate%n_q_bins) 

                kinematic_function_packed = &
                    exdm_inputs%astroph_model%kinematic_function(v_m_packed, q_mag_packed)

                do f = 1, size(binned_rate, 4)

                    F_med_sq_packed = (q_mag_packed*(alpha_EM*m_elec)**(-1))&
                        **(-2*exdm_inputs%dm_model%med_FF(f))

                    b_rate_q = 0.0_dp

                    b_rate_q(q_bin_packed) = b_rate_q(q_bin_packed) + &
                                                jac_q_packed*&
                                                kinematic_function_packed*&
                                                F_med_sq_packed*&
                                                FIF_packed*&
                                                screen_factor_list_packed**(-2)

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
