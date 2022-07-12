module exdm_elec_config_type

    use prec_util, only: dp

    use elec_config_bloch_PW_basis_type
    use elec_config_bloch_STO_basis_type
    use elec_config_bloch_single_PW_type

    implicit none

    type :: exdm_elec_config_t

        ! init
        integer :: n_init_groups
        type(elec_config_bloch_PW_basis_t) :: init_bloch_PW_basis_config
        type(elec_config_bloch_STO_basis_t) :: init_bloch_STO_basis_config
        type(elec_config_bloch_single_PW_t) :: init_bloch_single_PW_config

        ! fin
        type(elec_config_bloch_PW_basis_t) :: fin_bloch_PW_basis_config
        type(elec_config_bloch_STO_basis_t) :: fin_bloch_STO_basis_config
        type(elec_config_bloch_single_PW_t) :: fin_bloch_single_PW_config

        contains

            procedure :: load => exdm_elec_config_type_load

    end type

contains

    subroutine exdm_elec_config_type_load(self, exdm_inputs)

        use exdm_inputs_type

        implicit none

        class(exdm_elec_config_t) :: self
        type(exdm_inputs_t) :: exdm_inputs

        if ( exdm_inputs%control%verbose ) then
            print*, 'Initializing electronic configuration...'
            print*
        end if

        ! init
        call self%init_bloch_PW_basis_config%initialize(exdm_inputs, 'init')
        call self%init_bloch_STO_basis_config%initialize(exdm_inputs, 'init')
        call self%init_bloch_single_PW_config%initialize(exdm_inputs, 'init')

        call set_n_init_groups(self)

        ! fin
        call self%fin_bloch_PW_basis_config%initialize(exdm_inputs, 'fin')
        call self%fin_bloch_STO_basis_config%initialize(exdm_inputs, 'fin')
        call self%fin_bloch_single_PW_config%initialize(exdm_inputs, 'fin')

        if ( exdm_inputs%control%verbose ) then
            print*, 'Done initializing electronic configuration!'
            print*
        end if

    end subroutine

    subroutine set_n_init_groups(self)

        implicit none

        class(exdm_elec_config_t) :: self

        integer :: n

        self%n_init_groups = 0

        do n = 1, size(self%init_bloch_PW_basis_config%states)
            self%n_init_groups = max(self%init_bloch_PW_basis_config%states(n)%i, self%n_init_groups)
        end do

        do n = 1, size(self%init_bloch_STO_basis_config%states)
            self%n_init_groups = max(self%init_bloch_STO_basis_config%states(n)%i, self%n_init_groups)
        end do

    end subroutine

end module
