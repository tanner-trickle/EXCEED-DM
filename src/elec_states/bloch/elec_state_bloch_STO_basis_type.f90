module elec_state_bloch_STO_basis_type

    use prec_util, only: dp

    use hdf5_utils

    use elec_state_bloch_type

    implicit none

    type, extends(elec_state_bloch_t) :: elec_state_bloch_STO_basis_t

        integer :: n ! Principal quantum number, :math:`n`
        integer :: l ! Principal quantum number, :math:`\ell`
        integer :: m ! Principal quantum number, :math:`m`

        real(dp), allocatable :: nl_list(:)
        real(dp), allocatable :: Zl_list(:)
        real(dp), allocatable :: norml_list(:)
        real(dp), allocatable :: Cnl_list(:)

        real(dp) :: eq_pos_red(3)
        real(dp) :: eq_pos(3)

        integer :: n_r_vec_grid(3)
        real(dp), allocatable :: r_vec_grid(:, :)

        ! These variables do not really belong here 
        real(dp) :: red_to_xyz(3, 3)
        real(dp) :: pc_vol

        contains

            procedure :: compute_u => elec_state_bloch_STO_basis_type_compute_u

    end type

contains

    subroutine elec_state_bloch_STO_basis_type_compute_u(self, u)

        use constants_util
        use math_util

        implicit none

        class(elec_state_bloch_STO_basis_t) :: self
        complex(dp) :: u(:, :, :, :)
        complex(dp) :: u_list(product(self%n_x_grid))

        real(dp) :: x_vec_list(product(self%n_x_grid), 3)
        real(dp) :: y_vec_list(product(self%n_x_grid), 3)
        real(dp) :: r_list(product(self%n_x_grid))
        real(dp) :: theta_list(product(self%n_x_grid))
        real(dp) :: phi_list(product(self%n_x_grid))

        complex(dp) :: phase_fac_list(product(self%n_x_grid))
        complex(dp) :: sph_harm_list(product(self%n_x_grid))

        real(dp) :: radial_wf_list(product(self%n_x_grid))

        integer :: r

        integer :: i

        u = (0.0_dp, 0.0_dp)

        u_list = (0.0_dp, 0.0_dp)

        ! compute UC discretized grid
        call n_x_grid_to_x_vec_list(self%n_x_grid, self%red_to_xyz, x_vec_list)

        do r = 1, size(self%r_vec_grid, 1)

            ! compute relative positions
            y_vec_list = spread(self%r_vec_grid(r, :) - self%eq_pos, 1, size(x_vec_list, 1)) + x_vec_list

            ! compute phase factor list
            phase_fac_list = exp(-ii*matmul(y_vec_list, self%k_vec))

            ! compute r_list, theta_list, phi_list
            call cartesian_to_polar_list(y_vec_list, r_list, theta_list, phi_list)

            ! compute spherical harmonics list
            call compute_sph_harmonics(self%l, self%m, theta_list, phi_list, sph_harm_list)

            ! compute radial wf list
            radial_wf_list = 0.0_dp

            do i = 1, size(self%nl_list)

                radial_wf_list = radial_wf_list + &
                    self%Cnl_list(i)*STO_radial(r_list,&
                                                self%nl_list(i),&
                                                self%norml_list(i), &
                                                self%Zl_list(i))

            end do

            ! sum all pieces
            u_list = u_list + phase_fac_list*radial_wf_list*sph_harm_list

        end do

        ! NOTE: Assumption: core wave functions are not spin dependent
        ! normalize and reshape
        u(:, :, :, 1) = sqrt(self%pc_vol)*reshape(u_list, [size(u, 1), size(u, 2), size(u, 3)])

    end subroutine

    subroutine n_x_grid_to_x_vec_list(n_x_grid, red_to_xyz, x_vec_list)

        implicit none

        integer :: n_x_grid(:)

        real(dp) :: red_to_xyz(3, 3)
        real(dp) :: x_vec_array(n_x_grid(1), n_x_grid(2), n_x_grid(3), 3)

        real(dp) :: x_vec_list(:, :)

        integer :: n1, n2, n3

        do n3 = 1, n_x_grid(3)
            do n2 = 1, n_x_grid(2)
                do n1 = 1, n_x_grid(1)

                    x_vec_array(n1, n2, n3, :) = matmul(red_to_xyz,&
                            [ (n1 - 1.0_dp)/n_x_grid(1), (n2 - 1.0_dp)/n_x_grid(2), (n3 - 1.0_dp)/n_x_grid(3) ])

                end do
            end do
        end do

        x_vec_list = reshape(x_vec_array, [ size(x_vec_list, 1), size(x_vec_list, 2) ])

    end subroutine

end module
