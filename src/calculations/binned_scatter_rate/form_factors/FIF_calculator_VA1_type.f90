module FIF_calculator_VA1_type

    implicit none

    type :: FIF_calculator_VA1_t

        contains

            procedure :: set_TIF_mask => FIF_calculator_VA1_type_set_TIF_mask
            procedure :: compute => FIF_calculator_VA1_type_compute

    end type

contains

    subroutine FIF_calculator_VA1_type_compute(self, TIF_calculator, FIF)

        use constants_util

        use TIF_calculator_type

        implicit none

        class(FIF_calculator_VA1_t) :: self

        class(TIF_calculator_t) :: TIF_calculator

        real(dp) :: FIF(:)

        integer :: d

        real(dp) :: k_0

        FIF = 0.0_dp

        k_0 = alpha_EM*m_elec

        do d = 1, 3

            FIF = FIF + 4.0_dp*(m_elec/k_0)**2*abs(TIF_calculator%TIF_v(:, d))**2 &
            + 2.0_dp*(m_elec/k_0**2)*conjg(TIF_calculator%TIF_v(:, d))*TIF_calculator%q_vec_list(:, d)*&
                TIF_calculator%TIF_1 &
            + 2.0_dp*(m_elec/k_0**2)*TIF_calculator%TIF_v(:, d)*TIF_calculator%q_vec_list(:, d)*&
                conjg(TIF_calculator%TIF_1)

        end do

        FIF = FIF + k_0**(-2)*norm2(TIF_calculator%q_vec_list, 2)**2*abs(TIF_calculator%TIF_1)**2

    end subroutine

    subroutine FIF_calculator_VA1_type_set_TIF_mask(self, TIF_mask)

        implicit none

        class(FIF_calculator_VA1_t) :: self

        logical :: TIF_mask(:)

        TIF_mask(1) = .TRUE.
        TIF_mask(2) = .TRUE.

    end subroutine

end module
