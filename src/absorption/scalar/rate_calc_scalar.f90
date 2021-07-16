module rate_calc_scalar
    !! Given the self energies, compute the absorption rate of 
    !! scalar DM

    use prec
    use constants
    use particle_physics_abs
    use absorption_input
    use material_input

    implicit none

contains

    subroutine calc_rate_scalar(pi_1_1_mat, pi_v2_v2, v_vec, v_max, abs_rate, verbose)

        implicit none

        complex(dp) :: pi_1_1_mat(3, 3, n_omega, n_widths)
        complex(dp) :: pi_v2_v2(n_omega, n_widths)

        real(dp) :: abs_rate(n_omega, n_widths, n_time)

        logical, optional :: verbose

        integer :: w, p, t, v, a

        real(dp) :: gam
        real(dp) :: av_rate, rate

        real(dp) :: omega
        real(dp) :: v_angular_mesh(n_v_theta*n_v_phi, 2)

        real(dp) :: v_mag, v_max
        real(dp) :: v_vec(3)
        real(dp) :: q_vec(3), q_mag

        real(dp) :: pi_r, pi_c
        real(dp) :: pi_1_1

        real(dp) :: ve_vec(3)

        real(dp) :: mb_val

        v_mag = norm2(v_vec)

        do w = 1, n_omega

            omega = omega_list(w)

            do p = 1, n_widths
                do t = 1, n_time

                    ve_vec = vE_vec_list(t, :) 

                    av_rate = 0.0_dp

                    q_vec = omega*v_vec
                    q_mag = norm2(q_vec)

                    if ( q_mag > 0.0_dp ) then

                        pi_1_1 = dot_product( q_vec/m_elec, matmul( pi_1_1_mat(:, :, w, p), q_vec/m_elec ) )

                        gam = -(omega)**(-1)*aimag(&
                            pi_v2_v2(w, p) + &
                            (e_EM)**(-2)*( q_mag**2 )*&
                            ( q_mag**2 )*&
                            ( q_mag**2 - e_EM**2*pi_1_1 )**(-1)&
                            )

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
