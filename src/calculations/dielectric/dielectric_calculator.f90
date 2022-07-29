module dielectric_calculator

    use exdm_inputs_type 
    use elec_state_type

    implicit none

contains

    subroutine dielectric_compute(n_proc, proc_id, root_proc_id, &
                                    exdm_inputs, &
                                    init_states, fin_states, &
                                    jac_I_list, &
                                    jac_F_list, &
                                    Zeff_list, &
                                    dielectric, &
                                    fermi_factor_bool)

        use parallel_index_util 

        use TIF_calculator_type

        use scatter_physics_functions

        use dielectric_suscep_IF_calculator

        use timer_util

        implicit none

        integer :: n_proc, proc_id, root_proc_id

        type(exdm_inputs_t) :: exdm_inputs

        class(elec_state_t) :: init_states(:)
        class(elec_state_t) :: fin_states(:)

        real(dp) :: jac_I_list(:)
        real(dp) :: jac_F_list(:)

        real(dp) :: Zeff_list(:)

        complex(dp) :: dielectric(:, :, :, :, :) 

        logical, optional :: fermi_factor_bool
        logical :: fermi_factor_bool_in

        integer :: n_tran_all, n_tran
        integer :: prev_init_id, prev_fin_id

        integer :: init_id, fin_id
        integer :: tran_id
        integer :: state_ids(2)

        integer :: i

        type(TIF_calculator_t) :: TIF_calculator

        real(dp) :: extra_FF

        type(timer_t) :: timer

        n_tran_all = size(init_states)*size(fin_states)

        n_tran = get_n_transitions(proc_id, n_proc, n_tran_all)

        prev_init_id = 0
        prev_fin_id = 0

        ! initialize mask
        TIF_calculator%mask = .FALSE.

        ! set mask
        TIF_calculator%mask(1) = .TRUE.

        ! For all the states this processor has to compute for
        do i = 1, n_tran

            if ( i == 1 ) then

                ! initialize
                call TIF_calculator%setup(exdm_inputs, init_states(i), fin_states(i), q0_limit = .FALSE.)

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

            ! Optional: Fermi factor 
            if ( present(fermi_factor_bool) ) then
                fermi_factor_bool_in = fermi_factor_bool
            else
                fermi_factor_bool_in = .FALSE.
            end if

            if ( fermi_factor_bool_in ) then
                extra_FF = fermi_factor(Zeff_list(init_id), fin_states(fin_id)%energy)
            else
                extra_FF = 1.0_dp
            end if

            ! compute electric susceptability from TIF
            call dielectric_suscep_IF_compute(&
                dielectric, &
                fin_states(fin_id)%energy - init_states(init_id)%energy, &
                TIF_calculator%q_vec_list, &
                jac_I_list(init_id), &
                jac_F_list(fin_id), &
                min( init_states(init_id)%spin_dof, fin_states(fin_id)%spin_dof ), &
                exdm_inputs, &
                extra_FF, &
                TIF_calculator%TIF_1)

        end do

    end subroutine

end module
