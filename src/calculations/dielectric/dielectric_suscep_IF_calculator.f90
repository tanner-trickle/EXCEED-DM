module dielectric_suscep_IF_calculator

    use prec_util, only: dp

    use exdm_inputs_type

    implicit none

contains

    subroutine dielectric_suscep_IF_compute(dielectric, omega_IF, &
                                            q_vec_list, jac_I, jac_F, spin_dof, &
                                            exdm_inputs, extra_FF_IF, TIF_1, TIF_v_q0)

        use abs_physics_functions, only: width_func, green_func

        use constants_util 
        use math_util

        implicit none

        complex(dp) :: dielectric(:, :, :, :, :) 

        real(dp) :: omega, omega_IF

        real(dp) :: q_vec_list(:, :)
        complex(dp) :: TIF_1(:)
        complex(dp) :: TIF_v_q0(3)

        real(dp) :: jac_I, jac_F

        real(dp) :: extra_FF_IF

        integer :: spin_dof

        type(exdm_inputs_t) :: exdm_inputs

        integer :: e, w, q

        real(dp) :: delta

        integer :: q_bin, q_theta_bin, q_phi_bin
        real(dp) :: q_mag, q_theta, q_phi

        real(dp) :: q_theta_bin_width, q_phi_bin_width
        real(dp) :: q_bin_width

        complex(dp) :: ddi(size(dielectric, 1), size(dielectric, 2), size(dielectric, 3), size(dielectric, 4), size(dielectric, 5))
        real(dp) :: ddi_q(size(dielectric, 3), size(dielectric, 2), size(dielectric, 1))
        real(dp) :: ddi_q_T(size(dielectric, 1), size(dielectric, 2), size(dielectric, 3))

        integer :: q_t, q_p

        real(dp) :: V_q_bin, q_mag_bin

        real(dp) :: q_hat(3)

        complex(dp) :: TIF_vq

        real(dp) :: TIF_1_sq
        real(dp) :: TIF_v_sq

        ddi = ( 0.0_dp, 0.0_dp )
        ddi_q = 0.0_dp

        q_bin_width = 1.0e3_dp*exdm_inputs%numerics_dielectric%q_bin_width

        q_theta_bin_width = pi/exdm_inputs%numerics_dielectric%n_q_theta
        q_phi_bin_width = 2.0_dp*pi/exdm_inputs%numerics_dielectric%n_q_phi

        do q = 1, size(q_vec_list, 1)

            q_mag = norm2(q_vec_list(q, :)) + 1.0e-8_dp

            if ( q_mag < exdm_inputs%numerics_dielectric%n_q_bins*q_bin_width ) then

                q_hat = q_vec_list(q, :)/q_mag

                q_theta = acos( q_hat(3) )
                q_phi = mod(atan2( q_hat(2), q_hat(1) ) + 2.0_dp*pi, 2.0_dp*pi)

                q_bin = clamp_int(&
                    1 + floor( q_mag/q_bin_width ), &
                    1, exdm_inputs%numerics_dielectric%n_q_bins)

                q_theta_bin = clamp_int( 1 + floor(q_theta/q_theta_bin_width), &
                    1, max( 1, exdm_inputs%numerics_dielectric%n_q_theta ))
                q_phi_bin = clamp_int( 1 + floor(q_phi/q_phi_bin_width), &
                    1, max( 1, exdm_inputs%numerics_dielectric%n_q_phi ))

                V_q_bin = (1.0_dp/3.0_dp)*q_bin_width**3*( q_bin*(3*q_bin - 3) + 1 )*&
                                ( cos( (q_theta_bin - 1)*q_theta_bin_width ) - cos(q_theta_bin*q_theta_bin_width) )*&
                                ( q_phi_bin_width )

                TIF_vq = dot_product(q_hat, TIF_v_q0)

                TIF_1_sq = conjg(TIF_1(q))*TIF_1(q)
                TIF_v_sq = conjg(TIF_vq)*TIF_vq 

                ddi_q(q_phi_bin, q_theta_bin, q_bin) = ddi_q(q_phi_bin, q_theta_bin, q_bin) + &
                                                            (-2.0_dp/spin_dof)*jac_I*jac_F*&
                                                            (e_EM**2)*&
                                                            ( q_mag**2 + omega_IF**2*TIF_1_sq/TIF_v_sq )**(-1)*&
                                                            exdm_inputs%material%pc_vol**(-2)*&
                                                            TIF_1_sq*&
                                                            V_q_bin**(-1)*&
                                                            extra_FF_IF*&
                                                            (2.0_dp*pi)**3

            end if

        end do

        ddi_q_T = reshape(ddi_q, [size(ddi_q_T, 1), size(ddi_q_T, 2), size(ddi_q_T, 3)])

        do w = 1, size(exdm_inputs%numerics_dielectric%widths, 1)

            do e = 1, exdm_inputs%numerics_dielectric%n_E_bins

                omega = exdm_inputs%numerics_dielectric%E_bin_width*(e - 0.5_dp)

                delta = width_func(omega, exdm_inputs%numerics_dielectric%widths(w, 1), &
                                          exdm_inputs%numerics_dielectric%widths(w, 2), &
                                          exdm_inputs%numerics_dielectric%widths(w, 3))

                ddi(:, :, :, e, w) = ddi_q_T*green_func(omega, omega_IF, delta, &
                    exdm_inputs%numerics_dielectric%smear_type)


            end do

        end do

        dielectric = dielectric + ddi

    end subroutine

end module
