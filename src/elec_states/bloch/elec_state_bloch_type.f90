module elec_state_bloch_type

    use prec_util, only: dp

    use elec_state_type

    implicit none

    type, abstract, extends(elec_state_t) :: elec_state_bloch_t

        real(dp) :: k_vec_red(3) ! :math:`\mathbf{k}` in reduced coordinates
        real(dp) :: k_vec(3) ! :math:`\mathbf{k}` in xyz coordinates, units: eV

        integer :: n_x_grid(3) ! Number of points in :math:`\mathbf{x}` grid.

        integer :: FFT_plans(2, 8) ! FFT plans for bloch states, 1 - forward, 2 - backward

        contains

            procedure (compute_u_procedure), deferred :: compute_u

    end type

    abstract interface

        subroutine compute_u_procedure(self, u)

            use prec_util, only: dp

            import :: elec_state_bloch_t
            class(elec_state_bloch_t) :: self
            complex(dp) :: u(:, :, :, :)

        end subroutine

    end interface

end module
