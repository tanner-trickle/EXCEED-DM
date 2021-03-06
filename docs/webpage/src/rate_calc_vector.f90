module rate_calc_vector
    !! Compute the vector DM absorption rate.

    use prec
    use constants
    use math_mod

    use dm_model_type
    use expt_type
    use material_type
    use width_parameters_type

    use numerics_abs
    use physics_abs_functions

    implicit none

contains

    subroutine calc_rate_vector(pi_1_1_mat, v_vec, dm_model, expt, widths, &
            target_mat, numerics, abs_rate, verbose)
        !! Compute the vector DM absorption rate at a given \( \mathbf{v} \).

        implicit none

        complex(dp) :: pi_1_1_mat(:, :, :, :)
        real(dp) :: v_vec(3)
        type(dm_model_t) :: dm_model
        type(expt_t) :: expt
        type(width_parameters_t) :: widths
        type(material_t) :: target_mat
        type(numerics_abs_t) :: numerics
        real(dp) :: abs_rate(:, :, :)
        logical, optional :: verbose

        real(dp) :: v_mag

        integer :: w, p, t, v, a

        real(dp) :: gam
        real(dp) :: av_rate, rate

        real(dp) :: omega

        real(dp) :: q_vec(3), q_mag

        real(dp) :: pi_r, pi_c
        real(dp) :: ve_vec(3)

        real(dp) :: mb_val

        complex(dp) :: pi_eigvals(3)
        complex(dp) :: pi_eigvectors(3, 3)

        complex(dp) :: pi_mat(3, 3)

        integer :: i

        v_mag = norm2(v_vec)

        do w = 1, dm_model%n_mX

            omega = dm_model%mX(w)

            do p = 1, widths%n

                pi_mat = e_EM**2*(omega/m_elec)**2*pi_1_1_mat(:, :, w, p)

                call calc_eigvals_33(pi_mat, pi_eigvals)

                do t = 1, expt%n_time

                    ve_vec = dm_model%vE*expt%vE_direction(t, :)

                    av_rate = 0.0_dp

                    q_vec = omega*v_vec
                    q_mag = norm2(q_vec)

                    if ( q_mag > 0.0_dp ) then

                        gam = 0.0_dp

                        ! sum over (diagonalized) polarizations
                        do i = 1, 3

                            ! pi_r = real(pi_eigvals(i))
                            ! pi_c = aimag(pi_eigvals(i))

                            gam = gam + &
                                (-1.0_dp)*e_EM**(-2)*(3.0_dp*omega)**(-1)*&
                                aimag( omega**2*pi_eigvals(i) / ( omega**2 - pi_eigvals(i) ) )

                            ! gam = gam + &
                            !     (-1.0_dp)*(3.0_dp*omega)**(-1)*&
                            !     omega**4*( ( omega**2 - e_EM**2*pi_r )**2 + ( e_EM**2*pi_c )**2 )**(-1)*&
                            !     pi_c

                        end do

                        rate = (dm_model%rhoX/target_mat%rho_T)*(omega)**(-1)*gam

                        mb_val = mb_vel_distribution(v_vec, dm_model, boost_vec_in = ve_vec)

                        av_rate = av_rate + v_mag**2*(4.0_dp*pi*dm_model%vX_max)*&
                            (1.0_dp*size(numerics%v_mesh, 1))**(-1)*rate*mb_val

                    end if

                    abs_rate(w, p, t) = abs_rate(w, p, t) + av_rate*expt%m_T*expt%exposure

                end do
            end do
        end do

    end subroutine

end module
