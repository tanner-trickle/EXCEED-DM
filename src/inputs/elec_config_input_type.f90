module elec_config_input_type

    implicit none

    type :: elec_config_input_t

        character(len=512) :: filename 

        contains

            procedure :: set_defaults => elec_config_input_type_set_defaults
            procedure :: get_values => elec_config_input_type_get_values

    end type

contains

    subroutine elec_config_input_type_set_defaults(self, cfg)

        use m_config

        implicit none

        class(elec_config_input_t) :: self
        type(CFG_t) :: cfg

        call CFG_add(cfg, &
                     "elec_config_input%filename",&
                     "", &
                     "File containing the electronic configuration")

    end subroutine

    subroutine elec_config_input_type_get_values(self, cfg)

        use m_config

        implicit none

        class(elec_config_input_t) :: self
        type(CFG_t) :: cfg

        call CFG_get(cfg, "elec_config_input%filename", self%filename)

    end subroutine

end module
