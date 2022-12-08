module elec_state_atomic_STO_basis_type

    use prec_util, only: dp

    use elec_state_atomic_type 

    implicit none

    type, extends(elec_state_atomic_t) :: elec_state_atomic_STO_basis_t

        integer :: n 

        ! coefficients in STO basis
        real(dp), allocatable :: nl_list(:)
        real(dp), allocatable :: Zl_list(:)
        real(dp), allocatable :: norml_list(:)
        real(dp), allocatable :: Cnl_list(:)

        contains

            procedure :: compute_radial_wf => elec_state_atomic_STO_basis_type_compute_radial_wf
            procedure :: compute_d_radial_wf => elec_state_atomic_STO_basis_type_compute_d_radial_wf
            procedure :: compute_d2_radial_wf => elec_state_atomic_STO_basis_type_compute_d_radial_wf
            procedure :: compute_wf => elec_state_atomic_STO_basis_type_compute_wf

    end type

contains

    subroutine elec_state_atomic_STO_basis_type_compute_d2_radial_wf(self, d2_R_wf)
        ! Compute first radial derivative

        use constants_util
        use math_util

        implicit none

        class(elec_state_atomic_STO_basis_t) :: self
        real(dp) :: d2_R_wf(:)

        integer :: i

        d2_R_wf = 0.0_dp

        do i = 1, size(self%nl_list)

            d2_R_wf = d2_R_wf + &
                self%Cnl_list(i)*&
                ( -( self%nl_list(i) - 1 )/self%sph_x_list(:, 1)**2 &
                + ( (self%nl_list(i) - 1)/self%sph_x_list(:, 1) - self%Zl_list(i)/a0 )**2 )*&
                STO_radial(self%sph_x_list(:, 1),&
                                            self%nl_list(i),&
                                            self%norml_list(i), &
                                            self%Zl_list(i))

        end do

    end subroutine

    subroutine elec_state_atomic_STO_basis_type_compute_d_radial_wf(self, d_R_wf)
        ! Compute first radial derivative

        use constants_util
        use math_util

        implicit none

        class(elec_state_atomic_STO_basis_t) :: self
        real(dp) :: d_R_wf(:)

        integer :: i

        d_R_wf = 0.0_dp

        do i = 1, size(self%nl_list)

            d_R_wf = d_R_wf + &
                self%Cnl_list(i)*&
                ( ( self%nl_list(i) - 1 )/self%sph_x_list(:, 1) - self%Zl_list(i)/a0 )*&
                STO_radial(self%sph_x_list(:, 1),&
                                            self%nl_list(i),&
                                            self%norml_list(i), &
                                            self%Zl_list(i))

        end do

    end subroutine

    subroutine elec_state_atomic_STO_basis_type_compute_radial_wf(self, R_wf)

        use constants_util 
        use math_util

        implicit none

        class(elec_state_atomic_STO_basis_t) :: self
        real(dp) :: R_wf(:)

        integer :: i

        R_wf = 0.0_dp

        do i = 1, size(self%nl_list)

            R_wf = R_wf + &
                self%Cnl_list(i)*STO_radial(self%sph_x_list(:, 1),&
                                            self%nl_list(i),&
                                            self%norml_list(i), &
                                            self%Zl_list(i))

        end do

    end subroutine

    subroutine elec_state_atomic_STO_basis_type_compute_wf(self, wf)

        use constants_util 
        use math_util

        implicit none

        class(elec_state_atomic_STO_basis_t) :: self
        complex(dp) :: wf(:)

        real(dp) :: R_wf(size(self%sph_x_list, 1))
        complex(dp) :: sph_harm_list(size(self%sph_x_list, 1))

        wf = ( 0.0_dp, 0.0_dp )

        call self%compute_radial_wf(R_wf)

        call compute_sph_harmonics(self%l, self%m, self%sph_x_list(:, 2), self%sph_x_list(:, 3), sph_harm_list)

        wf = R_wf*sph_harm_list

    end subroutine

end module
