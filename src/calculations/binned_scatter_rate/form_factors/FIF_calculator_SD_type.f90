module FIF_calculator_SD_type

    implicit none

    type :: FIF_calculator_SD_t

        contains

            procedure :: set_TIF_mask => FIF_calculator_SD_type_set_TIF_mask
            procedure :: compute => FIF_calculator_SD_type_compute

    end type

contains

    subroutine FIF_calculator_SD_type_compute(self, TIF_calculator, FIF)

        use TIF_calculator_type

        implicit none

        class(FIF_calculator_SD_t) :: self
        class(TIF_calculator_t) :: TIF_calculator

        real(dp) :: FIF(:)

        integer :: d

        FIF = 0.0_dp

        do d = 1, 3

            FIF = FIF + (3.0_dp)**(-1)*abs(TIF_calculator%TIF_s(:, d))**2

        end do

    end subroutine

    subroutine FIF_calculator_SD_type_set_TIF_mask(self, TIF_mask)

        implicit none

        class(FIF_calculator_SD_t) :: self

        logical :: TIF_mask(:)

        TIF_mask(4) = .TRUE.

    end subroutine

end module
