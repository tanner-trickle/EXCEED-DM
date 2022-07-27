module exdm_inputs_type
    ! Module holding the type :code:`exdm_inputs_t` which is a container for all of the input variables. Input variables are handled
    ! by the `m_config` library.
    !
    ! For each section, e.g., [control] in the input file there should be a corresponding type inside exdm_inputs. To add more
    ! inputs simply follow the same layout as the others.
    !
    ! The idea is a) set the defaults in a configuration file (%set_defaults), b) read any overwritten data in the input files, c)
    ! copy the new configuration values to the ones being used in the program (%get_values).

    use control_type
    use dm_model_type
    use astroph_model_type
    use elec_config_input_type
    use material_type
    use experiment_type
    use screening_type

    ! calculation specific
    use numerics_binned_scatter_rate_type
    use numerics_absorption_rate_type
    use numerics_dielectric_type

    implicit none

    type :: exdm_inputs_t
        ! Contains all of the input parameters for EXCEED-DM

        type(control_t) :: control
            ! General program runtime control
        type(dm_model_t) :: dm_model
            ! Contains all of the DM model parameters
        type(elec_config_input_t) :: elec_config_input
            ! Contains information about the input electronic configuration
        type(material_t) :: material
            ! Contains information about the target material
        type(astroph_model_t) :: astroph_model
            ! Contains all of the astrophysica parameters
        type(experiment_t) :: experiment
            ! Contains all of the experimental parameters
        type(screening_t) :: screening
            ! Contains all of the screening parameters

        type(numerics_binned_scatter_rate_t) :: numerics_binned_scatter_rate
        type(numerics_absorption_rate_t) :: numerics_absorption_rate
        type(numerics_dielectric_t) :: numerics_dielectric

        contains

            procedure :: load => exdm_inputs_load
            procedure :: save => exdm_inputs_save

    end type

contains

    subroutine exdm_inputs_save(self)

        implicit none

        class(exdm_inputs_t) :: self

        call self%dm_model%save(self%control%out_filename)
        call self%experiment%save(self%control%out_filename)
        call self%material%save(self%control%out_filename)
        
        call self%numerics_binned_scatter_rate%save(self%control%out_filename)
        call self%numerics_absorption_rate%save(self%control%out_filename)
        call self%numerics_dielectric%save(self%control%out_filename)

    end subroutine

    subroutine exdm_inputs_load(self, proc_id, root_process, logger)
        ! Loads all the inputs for EXCEED-DM.

        use m_config

        use logger_util

        implicit none

        class(exdm_inputs_t) :: self
        integer :: proc_id, root_process

        type(logger_t) :: logger

        type(CFG_t) :: cfg
        type(CFG_t) :: default_cfg

        logical :: save_bool
        character(len=512) :: default_inputs_filename

        ! Initialize defaults for each input section
        call self%control%set_defaults(cfg)
        call self%dm_model%set_defaults(cfg)
        call self%elec_config_input%set_defaults(cfg)
        call self%material%set_defaults(cfg)
        call self%astroph_model%set_defaults(cfg)
        call self%experiment%set_defaults(cfg)
        call self%screening%set_defaults(cfg)

        call self%numerics_binned_scatter_rate%set_defaults(cfg)
        call self%numerics_absorption_rate%set_defaults(cfg)
        call self%numerics_dielectric%set_defaults(cfg)

        ! set the default configuration
        default_cfg = cfg

        ! Load all input files
        if ( command_argument_count() == 0 ) then

            print*, 'No input file specified, aborting.'
            print*

            stop
        end if

        call CFG_update_from_arguments(cfg)

        ! Set all variables
        call self%control%get_values(cfg, proc_id, root_process)
        call self%dm_model%get_values(cfg)
        call self%elec_config_input%get_values(cfg)
        call self%astroph_model%get_values(cfg)
        call self%material%get_values(cfg)
        call self%experiment%get_values(cfg)
        call self%screening%get_values(cfg)

        call self%numerics_binned_scatter_rate%get_values(cfg)
        call self%numerics_absorption_rate%get_values(cfg)
        call self%numerics_dielectric%get_values(cfg)

        ! Print
        if ( self%control%verbose ) then
            print*, repeat('-', logger%section_length)
            print*
            print*, 'Inputs:'
            print*
            call CFG_write(cfg, "stdout")
            print*, repeat('-', logger%section_length)
            print*
        end if

        ! Save the input parameters to a markdown file.
        if ( proc_id == root_process ) then
            if ( self%control%save_inputs_markdown ) then
                call CFG_write_markdown(cfg, self%control%input_markdown_filename, include_title = .FALSE.)
            end if
            if ( self%control%save_default_inputs_markdown ) then
                call CFG_write_markdown(default_cfg, self%control%default_input_markdown_filename, include_title = .FALSE.)
            end if
        end if

    end subroutine

end module
