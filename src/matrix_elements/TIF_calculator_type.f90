module TIF_calculator_type

    use prec_util, only: dp

    use exdm_inputs_type

    use elec_state_type
    use elec_state_bloch_type
    use elec_state_atomic_type

    use TIF_calculator_bloch_type
    use TIF_calculator_atomic_type

    implicit none

    type :: TIF_calculator_t

        type(TIF_calculator_bloch_t) :: TIF_calculator_bloch
        type(TIF_calculator_atomic_t) :: TIF_calculator_atomic

        integer :: n_TIF = 7 ! The number of matrix element operators
        logical :: mask(7) ! Which matrix element operators to compute

        real(dp), allocatable :: q_vec_list(:, :)
        real(dp), allocatable :: jac_q_list(:)

        complex(dp), allocatable :: TIF_1(:)
        complex(dp), allocatable :: TIF_v(:, :)
        complex(dp), allocatable :: TIF_v2(:)
        complex(dp), allocatable :: TIF_s(:, :)
        complex(dp), allocatable :: TIF_vds(:)
        complex(dp), allocatable :: TIF_vxs(:, :)
        complex(dp), allocatable :: TIF_vivj(:, :, :)

        contains

            procedure :: setup => TIF_calculator_type_setup
            procedure :: init_initialize => TIF_calculator_type_init_initialize
            procedure :: fin_initialize => TIF_calculator_type_fin_initialize
            procedure :: compute_TIF_1 => TIF_calculator_type_compute_TIF_1
            procedure :: compute_TIF_v => TIF_calculator_type_compute_TIF_v
            procedure :: compute_TIF_v2 => TIF_calculator_type_compute_TIF_v2
            procedure :: compute_TIF_s => TIF_calculator_type_compute_TIF_s
            procedure :: compute_TIF_vds => TIF_calculator_type_compute_TIF_vds
            procedure :: compute_TIF_vxs => TIF_calculator_type_compute_TIF_vxs
            procedure :: compute_TIF_vivj => TIF_calculator_type_compute_TIF_vivj
            procedure :: compute_all => TIF_calculator_type_compute_all

    end type

contains

    subroutine TIF_calculator_type_compute_all(self, init_state, fin_state, q0_limit)

        use timer_util

        implicit none

        class(TIF_calculator_t) :: self

        class(elec_state_t) :: init_state, fin_state

        logical, optional :: q0_limit

        type(timer_t) :: timer

        ! compute all relevant TIFs
        if ( self%mask(1) ) then
            call self%compute_TIF_1(init_state, fin_state, q0_limit)
        end if
        if ( self%mask(2) ) then
            call self%compute_TIF_v(init_state, fin_state, q0_limit)
        end if
        if ( self%mask(3) ) then
            call self%compute_TIF_v2(init_state, fin_state, q0_limit)
        end if
        if ( self%mask(4) ) then
            call self%compute_TIF_s(init_state, fin_state, q0_limit)
        end if
        if ( self%mask(5) ) then
            call self%compute_TIF_vds(init_state, fin_state, q0_limit)
        end if
        if ( self%mask(6) ) then
            call self%compute_TIF_vxs(init_state, fin_state, q0_limit)
        end if
        if ( self%mask(7) ) then
            call self%compute_TIF_vivj(init_state, fin_state, q0_limit)
        end if

        ! set q_vec_list and jac_q_list
        select type ( init_state )
        class is ( elec_state_bloch_t )

            select type ( fin_state )
            class is ( elec_state_bloch_t )

                if ( .not. q0_limit ) then

                    self%q_vec_list = spread(fin_state%k_vec - init_state%k_vec, 1, size(self%q_vec_list, 1)) &
                        + self%TIF_calculator_bloch%G_vec_list

                end if

                self%jac_q_list = 1.0_dp

            end select
        end select

        select type ( init_state )
        class is ( elec_state_atomic_t )

            select type ( fin_state )
            class is ( elec_state_atomic_t )

                if ( .not. q0_limit ) then

                    print*, 'ERROR: q /= 0 not implemented for transitions between atomic states.'

                end if

                self%jac_q_list = 1.0_dp

            end select
        end select

    end subroutine

    subroutine TIF_calculator_type_setup(self, exdm_inputs, init_state, fin_state, q0_limit)

        implicit none

        class(TIF_calculator_t) :: self

        type(exdm_inputs_t) :: exdm_inputs
        class(elec_state_t) :: init_state, fin_state

        logical, optional :: q0_limit

        integer :: n_q

        select type ( init_state )
        class is ( elec_state_bloch_t )

            select type ( fin_state )
            class is ( elec_state_bloch_t )

                call self%TIF_calculator_bloch%setup(exdm_inputs, self%mask, init_state, fin_state, q0_limit = q0_limit)
                n_q = self%TIF_calculator_bloch%n_q

            end select
        end select

        select type ( init_state )
        class is ( elec_state_atomic_t )

            select type ( fin_state )
            class is ( elec_state_atomic_t )

                call self%TIF_calculator_atomic%setup(exdm_inputs, self%mask, init_state, fin_state, q0_limit = q0_limit)
                n_q = self%TIF_calculator_atomic%n_q

            end select
        end select

        if ( q0_limit ) then
            n_q = 1
        end if

        ! allocate q_vec_list
        allocate(self%q_vec_list(n_q, 3), source = 0.0_dp)
        allocate(self%jac_q_list(n_q), source = 1.0_dp)

        ! allocate arrays which are specified by the mask
        if ( self%mask(1) ) then
            allocate(self%TIF_1(n_q), source = ( 0.0_dp, 0.0_dp ))
        end if
        if ( self%mask(2) ) then
            allocate(self%TIF_v(n_q, 3), source = ( 0.0_dp, 0.0_dp ))
        end if
        if ( self%mask(3) ) then
            allocate(self%TIF_v2(n_q), source = ( 0.0_dp, 0.0_dp ))
        end if
        if ( self%mask(4) ) then
            allocate(self%TIF_s(n_q, 3), source = ( 0.0_dp, 0.0_dp ))
        end if
        if ( self%mask(5) ) then
            allocate(self%TIF_vds(n_q), source = ( 0.0_dp, 0.0_dp ))
        end if
        if ( self%mask(6) ) then
            allocate(self%TIF_vxs(n_q, 3), source = ( 0.0_dp, 0.0_dp ))
        end if
        if ( self%mask(7) ) then
            allocate(self%TIF_vivj(n_q, 3, 3), source = ( 0.0_dp, 0.0_dp ))
        end if

    end subroutine

    subroutine TIF_calculator_type_init_initialize(self, init_state)

        implicit none

        class(TIF_calculator_t) :: self

        class(elec_state_t) :: init_state

        select type ( init_state )
        class is ( elec_state_bloch_t )

            call self%TIF_calculator_bloch%init_initialize(init_state, self%mask)

        class is ( elec_state_atomic_t )

            call self%TIF_calculator_atomic%init_initialize(init_state, self%mask)

        end select

    end subroutine

    subroutine TIF_calculator_type_fin_initialize(self, fin_state)

        implicit none

        class(TIF_calculator_t) :: self

        class(elec_state_t) :: fin_state

        select type ( fin_state )
        class is ( elec_state_bloch_t )

            call self%TIF_calculator_bloch%fin_initialize(fin_state)

        class is ( elec_state_atomic_t )

            call self%TIF_calculator_atomic%fin_initialize(fin_state)

        end select

    end subroutine

    subroutine TIF_calculator_type_compute_TIF_1(self, init_state, fin_state, q0_limit)

        implicit none

        class(TIF_calculator_t) :: self

        class(elec_state_t) :: init_state, fin_state

        logical, optional :: q0_limit

        select type ( init_state )
        class is ( elec_state_bloch_t )

            select type ( fin_state )
            class is ( elec_state_bloch_t )

                call self%TIF_calculator_bloch%compute_TIF_1(self%TIF_1, init_state, fin_state, q0_limit = q0_limit)

            end select
        end select

        select type ( init_state )
        class is ( elec_state_atomic_t )

            select type ( fin_state )
            class is ( elec_state_atomic_t )

                call self%TIF_calculator_atomic%compute_TIF_1(self%TIF_1, init_state, fin_state, q0_limit = q0_limit)

            end select
        end select

    end subroutine

    subroutine TIF_calculator_type_compute_TIF_v(self, init_state, fin_state, q0_limit)

        implicit none

        class(TIF_calculator_t) :: self

        class(elec_state_t) :: init_state, fin_state

        logical, optional :: q0_limit

        select type ( init_state )
        class is ( elec_state_bloch_t )

            select type ( fin_state )
            class is ( elec_state_bloch_t )

                call self%TIF_calculator_bloch%compute_TIF_v(self%TIF_v, init_state, fin_state, q0_limit = q0_limit)

            end select
        end select

        select type ( init_state )
        class is ( elec_state_atomic_t )

            select type ( fin_state )
            class is ( elec_state_atomic_t )

                call self%TIF_calculator_atomic%compute_TIF_v(self%TIF_v, init_state, fin_state, q0_limit = q0_limit)

            end select
        end select

    end subroutine

    subroutine TIF_calculator_type_compute_TIF_vivj(self, init_state, fin_state, q0_limit)

        implicit none

        class(TIF_calculator_t) :: self

        class(elec_state_t) :: init_state, fin_state

        logical, optional :: q0_limit

        select type ( init_state )
        class is ( elec_state_bloch_t )

            select type ( fin_state )
            class is ( elec_state_bloch_t )

                call self%TIF_calculator_bloch%compute_TIF_vivj(self%TIF_vivj, init_state, fin_state, q0_limit = q0_limit)

            end select
        end select

        select type ( init_state )
        class is ( elec_state_atomic_t )

            select type ( fin_state )
            class is ( elec_state_atomic_t )

                call self%TIF_calculator_atomic%compute_TIF_vivj(self%TIF_vivj, init_state, fin_state, q0_limit = q0_limit)

            end select
        end select

    end subroutine

    subroutine TIF_calculator_type_compute_TIF_v2(self, init_state, fin_state, q0_limit)

        implicit none

        class(TIF_calculator_t) :: self

        class(elec_state_t) :: init_state, fin_state

        logical, optional :: q0_limit

        select type ( init_state )
        class is ( elec_state_bloch_t )

            select type ( fin_state )
            class is ( elec_state_bloch_t )

                call self%TIF_calculator_bloch%compute_TIF_v2(self%TIF_v2, init_state, fin_state, q0_limit = q0_limit)

            end select
        end select

        select type ( init_state )
        class is ( elec_state_atomic_t )

            select type ( fin_state )
            class is ( elec_state_atomic_t )

                call self%TIF_calculator_atomic%compute_TIF_v2(self%TIF_v2, init_state, fin_state, q0_limit = q0_limit)

            end select
        end select

    end subroutine

    subroutine TIF_calculator_type_compute_TIF_s(self, init_state, fin_state, q0_limit)

        implicit none

        class(TIF_calculator_t) :: self

        class(elec_state_t) :: init_state, fin_state

        logical, optional :: q0_limit

        select type ( init_state )
        class is ( elec_state_bloch_t )

            select type ( fin_state )
            class is ( elec_state_bloch_t )

                call self%TIF_calculator_bloch%compute_TIF_s(self%TIF_s, init_state, fin_state, q0_limit = q0_limit)

            end select
        end select

        select type ( init_state )
        class is ( elec_state_atomic_t )

            select type ( fin_state )
            class is ( elec_state_atomic_t )

                call self%TIF_calculator_atomic%compute_TIF_s(self%TIF_s, init_state, fin_state, q0_limit = q0_limit)

            end select
        end select

    end subroutine

    subroutine TIF_calculator_type_compute_TIF_vds(self, init_state, fin_state, q0_limit)

        implicit none

        class(TIF_calculator_t) :: self

        class(elec_state_t) :: init_state, fin_state

        logical, optional :: q0_limit

        select type ( init_state )
        class is ( elec_state_bloch_t )

            select type ( fin_state )
            class is ( elec_state_bloch_t )

                call self%TIF_calculator_bloch%compute_TIF_vds(self%TIF_vds, init_state, fin_state, q0_limit = q0_limit)

            end select
        end select

        select type ( init_state )
        class is ( elec_state_atomic_t )

            select type ( fin_state )
            class is ( elec_state_atomic_t )

                call self%TIF_calculator_atomic%compute_TIF_vds(self%TIF_vds, init_state, fin_state, q0_limit = q0_limit)

            end select
        end select

    end subroutine

    subroutine TIF_calculator_type_compute_TIF_vxs(self, init_state, fin_state, q0_limit)

        implicit none

        class(TIF_calculator_t) :: self

        class(elec_state_t) :: init_state, fin_state

        logical, optional :: q0_limit

        select type ( init_state )
        class is ( elec_state_bloch_t )

            select type ( fin_state )
            class is ( elec_state_bloch_t )

                call self%TIF_calculator_bloch%compute_TIF_vxs(self%TIF_vxs, init_state, fin_state, q0_limit = q0_limit)

            end select
        end select

        select type ( init_state )
        class is ( elec_state_atomic_t )

            select type ( fin_state )
            class is ( elec_state_atomic_t )

                call self%TIF_calculator_atomic%compute_TIF_vxs(self%TIF_vxs, init_state, fin_state, q0_limit = q0_limit)

            end select
        end select

    end subroutine

end module
