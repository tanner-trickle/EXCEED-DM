module rate_calc_vector
    !! Given the self energies, compute the absorption rate of vector DM

    use prec
    use constants
    use particle_physics_abs
    use absorption_input
    use material_input

    implicit none

contains

    subroutine calc_rate_vector(pi_11_mat, abs_rate, verbose)

        implicit none

        complex(dp) :: pi_11_mat(3, 3, n_omega, n_widths)

        real(dp) :: abs_rate(n_omega, n_widths, n_time)

        logical, optional :: verbose

        integer :: w, p, t, v, a

        real(dp) :: gam
        real(dp) :: av_rate, rate

        real(dp) :: omega
        real(dp) :: v_angular_mesh(n_v_theta*n_v_phi, 2)

        real(dp) :: v_mag, v_theta, v_phi, v_max
        real(dp) :: v_mag_list(n_v_mag)
        real(dp) :: v_vec(3)
        real(dp) :: q_vec(3), q_mag

        real(dp) :: pi_r, pi_c
        real(dp) :: ve_vec(3)

        real(dp) :: mb_val

        v_max = vEsc + vE

        do v = 1, n_v_mag

            v_mag_list(v) = v_max*(v - 1.0_dp)/max(1.0_dp, n_v_mag - 1.0_dp)

        end do

        v_angular_mesh = generate_uniform_points_on_sphere(n_v_theta, n_v_phi)

        do w = 1, n_omega

            omega = omega_list(w)

            do p = 1, n_widths
                do t = 1, n_time

                    ve_vec = vE_vec_list(t, :) 

                    av_rate = 0.0_dp

                    ! velocity integral
                    do v = 1, n_v_mag
                        do a = 1, n_v_theta*n_v_phi

                            v_mag = v_mag_list(v)
                            v_theta = v_angular_mesh(a, 1)
                            v_phi = v_angular_mesh(a, 2)

                            v_vec(1) = v_mag*sin(v_theta)*cos(v_phi)
                            v_vec(2) = v_mag*sin(v_theta)*sin(v_phi)
                            v_vec(3) = v_mag*cos(v_theta)

                            q_vec = omega*v_vec
                            q_mag = norm2(q_vec)

                            if ( q_mag > 0.0_dp ) then
                                pi_r = real(dot_product( q_vec/m_elec, matmul( pi_11_mat(:, :, w, p), q_vec/m_elec ) ))
                                pi_c = aimag(dot_product( q_vec/m_elec, matmul( pi_11_mat(:, :, w, p), q_vec/m_elec ) ))

                                gam = -(omega)**(-1)*(q_mag**2*omega**2)*&
                                    ( ( q_mag**2 - e_EM**2*pi_r )**2 + (e_EM**2*pi_c)**2 )**(-1)*&
                                    pi_c

                                rate = (rhoX/rho_T)*(omega)**(-1)*gam

                                mb_val = mb_vel_distribution(v_vec, boost_vec_in = ve_vec)

                                av_rate = av_rate + v_mag**2*(4.0_dp*pi*v_max)*(1.0_dp*n_v_mag*n_v_theta*n_v_phi)**(-1)*&
                                                    rate*mb_val
                            end if

                        end do
                    end do

                    abs_rate(w, p, t) = av_rate

                end do
            end do
        end do

    end subroutine

end module
