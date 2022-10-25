module parallel_index_util
    ! A collection of useful procedures for keeping track of indicies when computing in parallel. 
    ! 
    ! .. note: No dependancies!

contains

    function transition_id_to_state_ids_abs(tran_id, n_init_states, n_final_states) result ( state_ids )

        implicit none

        integer :: tran_id
        integer :: n_all_transitions
        integer :: n_final_states
        integer :: n_init_states

        integer :: state_ids(3)

        state_ids(1) = (tran_id - 1)/(n_init_states*n_final_states) + 1 
        state_ids(2) = mod(tran_id - 1, n_init_states*n_final_states)/n_final_states + 1
        state_ids(3) = mod(tran_id - 1, n_final_states) + 1

    end function

    function get_n_transitions(proc_id, n_proc, n_all_transitions) result ( n_transitions )
        ! Returns the number of transitions that a given processor has to compute for.

        implicit none

        integer :: proc_id
        integer :: n_proc
        integer :: n_all_transitions
        integer :: n_transitions

        integer :: m, remainder

        m = n_all_transitions / n_proc 
        remainder = mod(n_all_transitions, n_proc)

        if ( proc_id < remainder ) then
            n_transitions = m + 1
        else
            n_transitions = m
        end if

    end function

    function get_transition_id(i, proc_id, n_proc, n_all_transitions) result ( transition_id )
        ! Returns the absolute id of the transition for processor proc_id.
        ! 1 <= transition_id <= n_all_transitions

        implicit none

        integer :: i
        integer :: proc_id
        integer :: n_proc
        integer :: n_all_transitions

        integer :: transition_id

        integer :: m, remainder

        integer :: start

        m = n_all_transitions / n_proc 
        remainder = mod(n_all_transitions, n_proc)

        if ( proc_id < remainder ) then
            start = proc_id*(m + 1) + 1
        else
            start = proc_id*m + remainder + 1
        end if

        transition_id = start + (i - 1)

    end function

    function transition_id_to_state_ids(tran_id, n_final_states) result ( state_ids )
        !* Absolute transition id to (init_state_id, fin_state_id) pair.

        implicit none

        integer :: tran_id
        integer :: n_all_transitions
        integer :: n_final_states

        integer :: state_ids(2)

        state_ids(1) = (tran_id - 1)/n_final_states + 1
        state_ids(2) = mod(tran_id - 1, n_final_states) + 1

    end function

end module
