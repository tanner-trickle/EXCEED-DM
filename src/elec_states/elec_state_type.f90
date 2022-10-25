module elec_state_type
    ! Defines the `electronic_state` type; an abstract type representing a quantum state :math:`| I \rangle`. 

    use prec_util, only: dp

    implicit none

    type, abstract :: elec_state_t
        ! An abstract type representing an electronic quantum state :math:`| I \rangle`. Meant to be extended by explicit
        ! realizations.
        real(dp) :: energy ! Energy, eV, :math:`E_I`
        integer :: i ! State index, :math:`I`
        integer :: spin_dof

    end type

end module
