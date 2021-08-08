module physics_abs_functions
    !! Physics functions needed for DM absorption rate calcuations.

    use prec
    use constants

    use dm_model_type

    implicit none

contains

    function mb_vel_distribution(v_vec, dm_model, boost_vec_in) result(mb_val)
        !! Maxwell-Boltzmann velocity distribution, boosted with \( \mathbf{v}_\text{boost} \).
        !!
        !! $$\begin{align*}
        !! f(\mathbf{v}) = \frac{1}{N_0} \exp( - |\mathbf{v} + \mathbf{v}_\text{boost}|^2/v_0^2)
        !! \Theta( v_\text{esc} - |\mathbf{v} + \mathbf{v}_\text{boost}|)
        !! \end{align*}$$
        
        implicit none

        real(dp) :: v_vec(3)
        type(dm_model_t) :: dm_model
        real(dp), optional :: boost_vec_in(3)

        real(dp) :: mb_val

        real(dp) :: v_mag

        real(dp) :: boost_vec(3)

        real(dp) :: v_p(3), v_p_mag

        ! boost the distribution
        if ( present(boost_vec_in) ) then
            boost_vec = boost_vec_in
        else
            boost_vec = [0.0_dp, 0.0_dp, 0.0_dp]
        end if

        v_p = v_vec + boost_vec
        v_p_mag = norm2(v_p)

        if ( v_p_mag < dm_model%vEsc ) then

            mb_val = (dm_model%g_func_N0)**(-1)*exp( -(v_p_mag/dm_model%v0)**2 )

        else

            mb_val = 0.0_dp

        end if

    end function

    subroutine check_mb_dist_normalization(v_mesh, dm_model, boost_vec_in, verbose)
        !! Check the integral of the MB velocity distribution. If the samepling is fine
        !! enough, the value of the integral should be 1.

        implicit none

        real(dp) :: v_mesh(:, :)
        type(dm_model_t) :: dm_model
        real(dp), optional :: boost_vec_in(3)
        logical, optional :: verbose

        integer :: v, t, p

        real(dp) :: int_val

        real(dp) :: boost_vec(3)

        real(dp) :: v_vec(3)
        real(dp) :: v_mag, v_theta, v_phi
        real(dp) :: v_max

        if ( present(boost_vec_in) ) then
            boost_vec = boost_vec_in
        else
            boost_vec = [0.0_dp, 0.0_dp, 0.0_dp]
        end if

        int_val = 0.0_dp

        do v = 1, size(v_mesh, 1)

            v_vec = v_mesh(v, :)
            v_mag = norm2(v_vec)

            int_val = int_val + v_mag**2*(4.0_dp*pi*dm_model%vX_max)*&
                (1.0_dp*size(v_mesh, 1))**(-1)*&
                mb_vel_distribution(v_vec, dm_model, boost_vec_in = boost_vec)

        end do

        if ( verbose ) then

            print*, 'Integral of MB velocity distrubution : ', int_val
            print*

        end if

    end subroutine

end module
