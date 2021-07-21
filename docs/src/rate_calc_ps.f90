module rate_calc_ps
    !! Given the self energies, compute the absorption rate of 
    !! pseudoscalar DM

    use prec
    use constants
    use particle_physics_abs
    use absorption_input
    use material_input

    implicit none

contains

    subroutine calc_rate_ps(pi_vi_vj, v_vec, v_max, abs_rate, verbose)

        implicit none

        complex(dp) :: pi_vi_vj(3, 3, n_omega, n_widths)

        real(dp) :: abs_rate(n_omega, n_widths, n_time)

        logical, optional :: verbose

        integer :: w, p, t, v, a

        real(dp) :: gam
        real(dp) :: av_rate, rate

        real(dp) :: omega

        real(dp) :: v_mag, v_max
        real(dp) :: v_vec(3)
        real(dp) :: q_vec(3), q_mag

        real(dp) :: pi_r, pi_c
        real(dp) :: ve_vec(3)

        real(dp) :: mb_val

        complex(dp) :: pi_eigvals(3)
        complex(dp) :: pi_eigvectors(3, 3)

        integer :: i

        v_mag = norm2(v_vec)

        do w = 1, n_omega

            omega = omega_list(w)

            do p = 1, n_widths

                call calc_eig_system_33(e_EM**2*pi_vi_vj(:, :, w, p), pi_eigvals, pi_eigvectors)

                do t = 1, n_time

                    ve_vec = vE_vec_list(t, :) 

                    av_rate = 0.0_dp

                    q_vec = omega*v_vec
                    q_mag = norm2(q_vec)

                    if ( q_mag > 0.0_dp ) then

                        gam = -omega**2*(4.0_dp*m_elec**2*omega*e_EM**2)**(-1)*&
                            aimag(sum(pi_eigvals))

                        rate = (rhoX/rho_T)*(omega)**(-1)*gam

                        mb_val = mb_vel_distribution(v_vec, boost_vec_in = ve_vec)

                        av_rate = av_rate + v_mag**2*(4.0_dp*pi*v_max)*(1.0_dp*n_v_mag*n_v_theta*n_v_phi)**(-1)*&
                                            rate*mb_val

                    end if

                    abs_rate(w, p, t) = abs_rate(w, p, t) + av_rate

                end do
            end do
        end do

    end subroutine

end module
