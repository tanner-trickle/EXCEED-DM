module elec_state_atomic_type

    use prec_util, only: dp

    use elec_state_type 

    implicit none

    type, abstract, extends(elec_state_t) :: elec_state_atomic_t
        ! States representing atomic orbitals, :math:`\psi -> R(r) Y(\theta, \phi)`

        integer :: l ! angular quantum number, l
        integer :: m ! angular quantum number, m

        real(dp), allocatable :: sph_x_list(:, :) ! {r, theta, phi}, [ N, 3 ]

        contains

            procedure (compute_radial_wf_procedure), deferred :: compute_radial_wf
            procedure (compute_wf_procedure), deferred :: compute_wf

    end type

    abstract interface

        subroutine compute_radial_wf_procedure(self, R_wf)

            use prec_util, only: dp

            import :: elec_state_atomic_t
            class(elec_state_atomic_t) :: self
            real(dp) :: R_wf(:)

        end subroutine

        subroutine compute_wf_procedure(self, wf)

            use prec_util, only: dp

            import :: elec_state_atomic_t
            class(elec_state_atomic_t) :: self
            complex(dp) :: wf(:)

        end subroutine

    end interface

end module
