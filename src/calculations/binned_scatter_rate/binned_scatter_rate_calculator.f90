module binned_scatter_rate_calculator

    use exdm_inputs_type

    use elec_state_type

    implicit none

contains

    subroutine binned_scatter_rate_compute(n_proc, proc_id, root_proc_id, &
                                            exdm_inputs, &
                                            init_states, fin_states, &
                                            jac_I_list, jac_F_list, &
                                            extra_FF, &
                                            binned_scatter_rate)

        use parallel_index_util

        use TIF_calculator_type
        use FIF_calculator_type

        use binned_scatter_rate_RIF_calculator

        implicit none

        integer :: n_proc, proc_id, root_proc_id

        type(exdm_inputs_t) :: exdm_inputs

        class(elec_state_t) :: init_states(:)
        class(elec_state_t) :: fin_states(:)

        real(dp) :: jac_I_list(:)
        real(dp) :: jac_F_list(:)

        real(dp) :: extra_FF(:, :)

        real(dp) :: binned_scatter_rate(:, :, :, :, :, :)

        type(TIF_calculator_t) :: TIF_calculator
        type(FIF_calculator_t) :: FIF_calculator

        integer :: n_tran_all, n_tran
        integer :: prev_init_id, prev_fin_id

        integer :: init_id, fin_id
        integer :: tran_id
        integer :: state_ids(2)

        integer :: i

        n_tran_all = size(init_states)*size(fin_states)
        n_tran = get_n_transitions(proc_id, n_proc, n_tran_all)

        prev_init_id = 0
        prev_fin_id = 0

        call FIF_calculator%ID_to_TIF_mask(exdm_inputs%dm_model%FIF_id, TIF_calculator%mask)

        ! For all the states this processor has to compute for
        do i = 1, n_tran

            if ( i == 1 ) then

                ! initialize
                call TIF_calculator%setup(exdm_inputs, init_states(i), fin_states(i), q0_limit = .FALSE.)
                allocate(FIF_calculator%FIF(size(TIF_calculator%q_vec_list, 1)), source = 0.0_dp)

            end if

            tran_id = get_transition_id(i, proc_id, n_proc, n_tran_all)
            state_ids = transition_id_to_state_ids(tran_id, size(fin_states))

            init_id = state_ids(1)
            fin_id = state_ids(2)

            if ( prev_init_id /= init_id ) then
                ! Compute all parts of TIF dependent on the initial state
                call TIF_calculator%init_initialize(init_states(init_id))
                prev_init_id = init_id
            end if

            if ( prev_fin_id /= fin_id ) then
                ! Compute all parts of TIF dependent on the final state
                call TIF_calculator%fin_initialize(fin_states(fin_id))
                prev_fin_id = fin_id
            end if

            ! compute all relevant TIFs
            call TIF_calculator%compute_all(init_states(init_id), fin_states(fin_id), &
                q0_limit = .FALSE.)

            ! compute FIF from TIF
            call FIF_calculator%compute(exdm_inputs%dm_model%FIF_id, TIF_calculator)

            ! compute RIF from FIF
            call binned_scatter_rate_RIF_compute(&
                binned_scatter_rate(:, :, :, :, :, init_states(init_id)%i), &
                fin_states(fin_id)%energy - init_states(init_id)%energy, &
                TIF_calculator%q_vec_list, &
                FIF_calculator%FIF, &
                TIF_calculator%jac_q_list, &
                jac_I_list(init_id), &
                jac_F_list(fin_id), &
                min( init_states(init_id)%spin_dof, fin_states(fin_id)%spin_dof ), &
                extra_FF(init_id, fin_id), &
                exdm_inputs)

        end do

    end subroutine

end module
