module self_energies_calculator

    use prec_util, only: dp

    use exdm_inputs_type

    use elec_state_type

    implicit none

contains

    subroutine self_energies_compute_continuum_atomic(n_proc, proc_id, root_proc_id, &
                                        exdm_inputs, &
                                        init_states, &
                                        Zeff_list, &
                                        PiIF_calculator)

        use elec_state_atomic_STO_basis_type
        use elec_state_atomic_continuum_type

        use parallel_index_util 

        use TIF_calculator_type 
        use PiIF_calculator_type

        use constants_util

        implicit none

        integer :: n_proc, proc_id, root_proc_id

        type(exdm_inputs_t) :: exdm_inputs
        
        class(elec_state_atomic_STO_basis_t) :: init_states(:)
        type(elec_state_atomic_continuum_t) :: fin_state

        real(dp) :: Zeff_list(:)

        type(TIF_calculator_t) :: TIF_calculator

        type(PiIF_calculator_t) :: PiIF_calculator

        integer :: i

        integer :: n_tran
        integer :: n_tran_all

        integer :: tran_id

        integer :: k_id, init_band_id, fin_band_id

        integer :: state_ids(2)

        integer :: init_id, mX_id
        integer :: prev_init_id

        integer :: lp, mp
        integer :: lp_min, lp_max
        integer :: mp_min, mp_max

        real(dp) :: k_F

        n_tran_all = size(init_states)*size(exdm_inputs%dm_model%mX)
        n_tran = get_n_transitions(proc_id, n_proc, n_tran_all)

        ! initialize mask
        TIF_calculator%mask = .FALSE.

        ! set mask
        call PiIF_mask_to_TIF_mask(PiIF_calculator%mask, TIF_calculator%mask)

        prev_init_id = 0

        ! For all the states this processor has to compute for
        do i = 1, n_tran

            if ( i == 1 ) then

                ! initialize
                call TIF_calculator%setup(exdm_inputs, init_states(i), fin_state, q0_limit = .TRUE.)

            end if

            ! get absolute transition id
            tran_id = get_transition_id(i, proc_id, n_proc, n_tran_all)

            ! get unique (i, omega) from transition id
            state_ids = transition_id_to_state_ids(tran_id, size(exdm_inputs%dm_model%mX))

            init_id = state_ids(1)
            mX_id = state_ids(2)

            ! make sure the final state is kinematically accessible
            if ( exdm_inputs%dm_model%mX(mX_id) + init_states(init_id)%energy > 0.0_dp ) then

                if ( prev_init_id /= init_id ) then
                    ! Compute all parts of TIF dependent on the initial state
                    call TIF_calculator%init_initialize(init_states(init_id))
                    prev_init_id = init_id

                    lp_min = TIF_calculator%TIF_calculator_atomic%lp_min
                    lp_max = TIF_calculator%TIF_calculator_atomic%lp_max

                    mp_min = TIF_calculator%TIF_calculator_atomic%mp_min
                    mp_max = TIF_calculator%TIF_calculator_atomic%mp_max

                end if

                ! set up the final state
                k_F = sqrt(2.0_dp*m_elec*(exdm_inputs%dm_model%mX(mX_id) + init_states(init_id)%energy))
                fin_state%k = k_F
                fin_state%energy = (k_F**2)/(2.0_dp*m_elec)
                fin_state%Z_eff = Zeff_list(init_id)
                ! fin_state%Z_eff = init_states(init_id)%n*sqrt(-init_states(init_id)%energy/13.6)
                fin_state%i = 1
                fin_state%spin_dof = 1
                fin_state%sph_x_list = init_states(init_id)%sph_x_list

                do lp = max( 0, lp_min ), lp_max 

                    fin_state%l = lp

                    call fin_state%compute_radial_wf(TIF_calculator%TIF_calculator_atomic%R_F)

                    do mp = max( -lp, mp_min ), min( mp_max, lp )

                        fin_state%m = mp

                        ! compute all relevant TIFs
                        call TIF_calculator%compute_all(init_states(init_id), &
                                                        fin_state, &
                                                        q0_limit = .TRUE.)

                        ! compute all relevant PiIF's
                        call PiIF_calculator%compute_all_atomic_continuum(TIF_calculator, &
                            exdm_inputs%material%n_T, &
                            exdm_inputs%dm_model%mX(mX_id), &
                            mX_id, &
                            init_states(init_id)%energy, &
                            init_states(init_id)%l)

                    end do

                end do

            end if

        end do

    end subroutine

    subroutine self_energies_compute(n_proc, proc_id, root_proc_id, &
                                        exdm_inputs, &
                                        init_states, fin_states, &
                                        init_kid_list, fin_kid_list, &
                                        jac_I_list, &
                                        PiIF_calculator)

        use parallel_index_util 

        use TIF_calculator_type 
        use PiIF_calculator_type

        implicit none

        integer :: n_proc, proc_id, root_proc_id

        type(exdm_inputs_t) :: exdm_inputs
        
        class(elec_state_t) :: init_states(:)
        class(elec_state_t) :: fin_states(:)

        integer :: init_kid_list(:)
        integer :: fin_kid_list(:)

        real(dp) :: jac_I_list(:)

        type(TIF_calculator_t) :: TIF_calculator

        type(PiIF_calculator_t) :: PiIF_calculator

        integer :: i

        integer :: n_tran
        integer :: n_tran_all

        integer :: tran_id

        integer :: k_id, init_band_id, fin_band_id

        integer :: state_ids(3)

        integer :: init_id, fin_id
        integer :: prev_init_id, prev_fin_id

        integer :: n_I_band, n_F_band, n_k

        integer, allocatable :: init_abs_ik_to_id(:, :)
        integer, allocatable :: fin_abs_ik_to_id(:, :)

        ! number of vertical transition slices
        n_k = maxval(init_kid_list)

        ! number of initial states at each k
        n_I_band = size(init_states)/n_k
        ! number of final states at each k
        n_F_band = size(fin_states)/n_k

        allocate(init_abs_ik_to_id(n_I_band, n_k), source = 1)
        allocate(fin_abs_ik_to_id(n_F_band, n_k), source = 1)

        ! total number of absorption transitions
        n_tran_all = n_k*n_I_band*n_F_band
        n_tran = get_n_transitions(proc_id, n_proc, n_tran_all)

        ! initialize mask
        TIF_calculator%mask = .FALSE.

        ! set mask
        call PiIF_mask_to_TIF_mask(PiIF_calculator%mask, TIF_calculator%mask)

        prev_init_id = 0
        prev_fin_id = 0

        ! For all the states this processor has to compute for
        do i = 1, n_tran

            if ( i == 1 ) then

                ! initialize
                call TIF_calculator%setup(exdm_inputs, init_states(i), fin_states(i), q0_limit = .TRUE.)

                call create_abs_ik_to_id(init_kid_list, init_abs_ik_to_id)
                call create_abs_ik_to_id(fin_kid_list, fin_abs_ik_to_id)

            end if

            ! get absolute transition id
            tran_id = get_transition_id(i, proc_id, n_proc, n_tran_all)

            ! get unique (k, i, f) from transition id
            state_ids = transition_id_to_state_ids_abs(tran_id, n_I_band, n_F_band)

            k_id = state_ids(1)
            init_band_id = state_ids(2)
            fin_band_id = state_ids(3)

            ! indicies in init_states() and fin_states()
            init_id = init_abs_ik_to_id(init_band_id, k_id)
            fin_id = fin_abs_ik_to_id(fin_band_id, k_id) 

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
            call TIF_calculator%compute_all(init_states(init_id), &
                                            fin_states(fin_id), &
                                            q0_limit = .TRUE.)

            ! compute all relevant PiIF's
            call PiIF_calculator%compute_all(TIF_calculator, &
                fin_states(fin_id)%energy - init_states(init_id)%energy, &
                exdm_inputs%dm_model%mX, &
                exdm_inputs%numerics_absorption_rate%widths, &
                jac_I_list(init_id), &
                exdm_inputs%material%pc_vol, &
                min( init_states(init_id)%spin_dof, fin_states(fin_id)%spin_dof ), &
                exdm_inputs%numerics_absorption_rate%smear_type)

        end do

    end subroutine

    subroutine create_abs_ik_to_id(kid_list, abs_ik_to_id)

        implicit none

        integer :: kid_list(:)
        integer :: abs_ik_to_id(:, :)

        integer :: i

        integer :: band_count(size(abs_ik_to_id, 2))

        integer :: k_id, band_id

        band_count = 0

        do i = 1, size(kid_list)

            k_id = kid_list(i)

            band_count(k_id) = band_count(k_id) + 1

            band_id = band_count(k_id)

            abs_ik_to_id(band_id, k_id) = i

        end do

    end subroutine

end module
