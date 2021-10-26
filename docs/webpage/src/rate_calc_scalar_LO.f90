module rate_calc_scalar_LO
    !! Compute the scalar DM absorption rate (leading order contribution, from \( \Pi_{\bar{v}^2, \bar{v}^2} \) ), only).

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

    subroutine calc_rate_scalar_LO(pi_v2_v2, v_vec, dm_model, expt, widths, &
            target_mat, numerics, abs_rate, verbose)
        !! Compute the scalar DM absorption rate from the LO contribution, at a given \( \mathbf{v} \).

        implicit none

        complex(dp) :: pi_v2_v2(:, :)
        real(dp) :: v_vec(3)
        type(dm_model_t) :: dm_model
        type(expt_t) :: expt
        type(width_parameters_t) :: widths
        type(material_t) :: target_mat
        type(numerics_abs_t) :: numerics
        real(dp) :: abs_rate(:, :, :)
        logical, optional :: verbose

        integer :: w, p, t, v, a

        real(dp) :: gam
        real(dp) :: av_rate, rate

        real(dp) :: omega

        real(dp) :: v_mag, v_max
        real(dp) :: q_vec(3), q_mag

        real(dp) :: pi_r, pi_c
        real(dp) :: pi_1_1

        real(dp) :: ve_vec(3)

        real(dp) :: mb_val

        v_mag = norm2(v_vec)

        do w = 1, dm_model%n_mX

            omega = dm_model%mX(w)

            do p = 1, widths%n
                do t = 1, expt%n_time

                    ve_vec = dm_model%vE*expt%vE_direction(t, :)

                    av_rate = 0.0_dp

                    q_vec = omega*v_vec
                    q_mag = norm2(q_vec)

                    if ( q_mag > 0.0_dp ) then

                        gam = -(omega)**(-1)*aimag( pi_v2_v2(w, p) )

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
