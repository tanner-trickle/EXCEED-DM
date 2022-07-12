module elec_state_bloch_single_PW_type

    use prec_util, only: dp

    use elec_state_bloch_type

    implicit none

    type, extends(elec_state_bloch_t) :: elec_state_bloch_single_PW_t

        real(dp) :: p_vec(3)
        real(dp) :: G_vec(3)

        integer :: G_vec_red(3)

        contains

            procedure :: compute_u => elec_state_bloch_single_PW_type_compute_u

    end type

contains

    subroutine elec_state_bloch_single_PW_type_compute_u(self, u)

        use constants_util

        implicit none

        class(elec_state_bloch_single_PW_t) :: self
        complex(dp) :: u(:, :, :, :)

        real(dp) :: x_vec_red_array(self%n_x_grid(1), self%n_x_grid(2), self%n_x_grid(3), 3)
        real(dp) :: x_vec_red_list(product(self%n_x_grid), 3)

        complex(dp) :: u_list(product(self%n_x_grid))

        integer :: n1, n2, n3

        do n3 = 1, self%n_x_grid(3)
            do n2 = 1, self%n_x_grid(2)
                do n1 = 1, self%n_x_grid(1)

                    x_vec_red_array(n1, n2, n3, :) = [ (n1 - 1.0_dp)/self%n_x_grid(1), &
                                                       (n2 - 1.0_dp)/self%n_x_grid(2), &
                                                       (n3 - 1.0_dp)/self%n_x_grid(3) ]

                end do
            end do
        end do

        x_vec_red_list = reshape(x_vec_red_array, [ size(x_vec_red_list, 1), size(x_vec_red_list, 2) ])

        u_list = exp(2.0_dp*pi*ii*matmul(x_vec_red_list, self%G_vec_red))

        ! Note: Assumption: spin independent free states
        u(:, :, :, 1) = reshape(u_list, [ size(u, 1), size(u, 2), size(u, 3) ])

    end subroutine

    subroutine compute_kG_from_p(p_vec, k_vec, G_vec, k_red_to_xyz, k_xyz_to_red)

        implicit none

        real(dp) :: p_vec(3)
        real(dp) :: k_vec(3)
        real(dp) :: G_vec(3)

        real(dp) :: k_red_to_xyz(3, 3)
        real(dp) :: k_xyz_to_red(3, 3)

        real(dp) :: p_red(3)

        integer :: G_red_closest(3, 8)

        integer :: G_red(3)

        real(dp) :: diff_vec_list(3, 8)
        real(dp) :: diff_mag_list(8)
        integer :: diff_mag_min_index(1)

        p_red = matmul(k_xyz_to_red, p_vec)

        G_red_closest(:, 1) = [ floor(p_red(1)), floor(p_red(2)), floor(p_red(3)) ]
        G_red_closest(:, 2) = [ floor(p_red(1)), floor(p_red(2)), ceiling(p_red(3)) ]
        G_red_closest(:, 3) = [ floor(p_red(1)), ceiling(p_red(2)), floor(p_red(3)) ]
        G_red_closest(:, 4) = [ ceiling(p_red(1)), floor(p_red(2)), floor(p_red(3)) ]
        G_red_closest(:, 5) = [ floor(p_red(1)), ceiling(p_red(2)), ceiling(p_red(3)) ]
        G_red_closest(:, 6) = [ ceiling(p_red(1)), floor(p_red(2)), ceiling(p_red(3)) ]
        G_red_closest(:, 7) = [ ceiling(p_red(1)), ceiling(p_red(2)), floor(p_red(3)) ]
        G_red_closest(:, 8) = [ ceiling(p_red(1)), ceiling(p_red(2)), ceiling(p_red(3)) ]

        diff_vec_list = spread(p_vec, 2, 8) - matmul(k_red_to_xyz, G_red_closest)

        diff_mag_list = norm2(diff_vec_list, 1)

        diff_mag_min_index = minloc(diff_mag_list, .TRUE.)

        G_red = G_red_closest(:, diff_mag_min_index(1))

        G_vec = matmul(k_red_to_xyz, G_red)

        k_vec = p_vec - G_vec

    end subroutine

end module 
