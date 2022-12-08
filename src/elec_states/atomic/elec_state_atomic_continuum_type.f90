module elec_state_atomic_continuum_type

    use prec_util, only: dp

    use elec_state_atomic_type 

    implicit none

    type, extends(elec_state_atomic_t) :: elec_state_atomic_continuum_t

        real(dp) :: k ! momentum, Units : eV
        real(dp) :: Z_eff

        contains

            procedure :: compute_radial_wf => elec_state_atomic_continuum_type_compute_radial_wf
            procedure :: compute_wf => elec_state_atomic_continuum_type_compute_wf

    end type

contains

    subroutine Hydrogenic_Continuum_Wave_Single_L(k, r, Z, l, phi_wf)
        ! All inputs are in natural (eV) units.
        ! computes phi_{kl} (dimensionless)

        use constants_util

        use h_continuum_gautchi, only: Hydrogenic_Continuum_Wave

        implicit none

        real(dp) :: k, r, phi_wf, Z

        integer :: l

        real(dp), allocatable :: phi_wf_l(:)

        allocate(phi_wf_l(l+1), source = 0.0_dp)

        call Hydrogenic_Continuum_Wave(k*a0, r/a0, Z, l, phi_wf_l)

        phi_wf = phi_wf_l(l+1)

    end subroutine

    subroutine elec_state_atomic_continuum_type_compute_radial_wf(self, R_wf)
        ! Units: eV

        use constants_util 
        use math_util

        implicit none

        class(elec_state_atomic_continuum_t) :: self
        real(dp) :: R_wf(:)
        real(dp) :: phi_wf(size(R_wf))

        real(dp) :: conversion_factor

        integer :: r

        R_wf = 0.0_dp
        phi_wf = 0.0_dp

        do r = 1, size(self%sph_x_list, 1)

            call Hydrogenic_Continuum_Wave_Single_L(self%k, self%sph_x_list(r, 1), &
                self%Z_eff, self%l, phi_wf(r))

        end do

        ! R = phi / r
        R_wf = phi_wf/self%sph_x_list(:, 1)

    end subroutine

    subroutine elec_state_atomic_continuum_type_compute_wf(self, wf)

        use constants_util 
        use math_util

        implicit none

        class(elec_state_atomic_continuum_t) :: self
        complex(dp) :: wf(:)

        real(dp) :: R_wf(size(self%sph_x_list, 1))
        complex(dp) :: sph_harm_list(size(self%sph_x_list, 1))

        wf = ( 0.0_dp, 0.0_dp )

        call self%compute_radial_wf(R_wf)

        call compute_sph_harmonics(self%l, self%m, self%sph_x_list(:, 2), self%sph_x_list(:, 3), sph_harm_list)

        wf = R_wf*sph_harm_list

    end subroutine

end module
