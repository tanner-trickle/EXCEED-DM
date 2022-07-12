module FIF_calculator_SI_type

    implicit none

    type :: FIF_calculator_SI_t

        contains

            procedure :: set_TIF_mask => FIF_calculator_SI_type_set_TIF_mask
            procedure :: compute => FIF_calculator_SI_type_compute

    end type

contains

    subroutine FIF_calculator_SI_type_compute(self, TIF_calculator, FIF)

        use TIF_calculator_type

        implicit none

        class(FIF_calculator_SI_t) :: self
        class(TIF_calculator_t) :: TIF_calculator

        real(dp) :: FIF(:)

        FIF = abs(TIF_calculator%TIF_1)**2

    end subroutine

    subroutine FIF_calculator_SI_type_set_TIF_mask(self, TIF_mask)

        implicit none

        class(FIF_calculator_SI_t) :: self

        logical :: TIF_mask(:)

        TIF_mask(1) = .TRUE.

    end subroutine

end module
