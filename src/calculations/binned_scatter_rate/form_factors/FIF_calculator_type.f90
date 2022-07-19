module FIF_calculator_type

    use prec_util, only: dp

    use FIF_calculator_SI_type
    use FIF_calculator_SD_type

    implicit none

    type :: FIF_calculator_t

        real(dp), allocatable :: FIF(:)

        type(FIF_calculator_SI_t) :: FIF_calculator_SI
        type(FIF_calculator_SD_t) :: FIF_calculator_SD

        contains

            procedure :: ID_to_TIF_mask => FIF_calculator_type_ID_to_TIF_mask
            procedure :: compute => FIF_calculator_type_compute

    end type

contains

    subroutine FIF_calculator_type_compute(self, FIF_id, TIF_calculator)

        use TIF_calculator_type

        implicit none

        class(FIF_calculator_t) :: self
        class(TIF_calculator_t) :: TIF_calculator

        character(len=*) :: FIF_id

        select case ( trim(adjustl(FIF_id)) )

            case ( 'SI' )

                call self%FIF_calculator_SI%compute(TIF_calculator, self%FIF)

            case ( 'SD' )

                call self%FIF_calculator_SD%compute(TIF_calculator, self%FIF)

        end select

    end subroutine

    subroutine FIF_calculator_type_ID_to_TIF_mask(self, FIF_id, TIF_mask)

        implicit none

        class(FIF_calculator_t) :: self

        character(len=*) :: FIF_id
        logical :: TIF_mask(:)

        select case ( trim(adjustl(FIF_id)) )

            case ( 'SI' )

                call self%FIF_calculator_SI%set_TIF_mask(TIF_mask)

            case ( 'SD' )

                call self%FIF_calculator_SD%set_TIF_mask(TIF_mask)

        end select

    end subroutine

end module
