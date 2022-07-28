module control_type
    ! Defines the type containing all of the control variables.

    implicit none

    type :: control_t

        logical :: verbose

        logical :: save_inputs_markdown
        character(len=512):: input_markdown_filename

        logical :: save_default_inputs_markdown
        character(len=512) :: default_input_markdown_filename

        character(len=512) :: run_description
        character(len=512) :: out_folder

        character(len=512) :: out_filename

        character(len=512) :: calculation

        contains

            procedure :: set_defaults => control_type_set_defaults
            procedure :: get_values => control_type_get_values

    end type

contains

    subroutine control_type_set_defaults(self, cfg)
        ! Defines the default values inside the dm_model

        use m_config

        implicit none

        class(control_t) :: self
        type(CFG_t) :: cfg

        call CFG_add(cfg,                       &
                     "control%calculation",             &
                     "", &
                     "Which calculation to perform")

        call CFG_add(cfg,                       &
                     "control%verbose",             &
                     .TRUE., &
                     "Toggle output printing to the console")

        call CFG_add(cfg,                       &
                     "control%save_inputs_markdown",             &
                     .FALSE., &
                     "Toggle to save the input parameters to a markdown file")

        call CFG_add(cfg,                       &
                     "control%input_markdown_filename",             &
                     "./inputs.md", &
                     "Filename to store input parameters in markdown format to")

        call CFG_add(cfg,                       &
                     "control%save_default_inputs_markdown",             &
                     .FALSE., &
                     "Toggle to save the default input parameters to a markdown file")

        call CFG_add(cfg,                       &
                     "control%default_input_markdown_filename",             &
                     "./default_inputs.md", &
                     "Filename to store default input parameters in markdown format to")

        call CFG_add(cfg, &
                     "control%run_description", &
                     "", &
                     "Small description of calculation which will be appended to 'EXDM_out_' to set the output filename")

        call CFG_add(cfg, &
                     "control%out_folder", &
                     "./", &
                     "Folder to store the ouput data")

    end subroutine

    subroutine control_type_get_values(self, cfg, proc_id, root_process)

        use m_config

        implicit none

        class(control_t) :: self
        type(CFG_t) :: cfg
        integer :: proc_id, root_process

        self%verbose = .FALSE.

        if ( proc_id == root_process ) then
            call CFG_get(cfg, "control%verbose", self%verbose)
        end if

        call CFG_get(cfg, "control%save_inputs_markdown", self%save_inputs_markdown)
        call CFG_get(cfg, "control%input_markdown_filename", self%input_markdown_filename)

        call CFG_get(cfg, "control%save_default_inputs_markdown", self%save_default_inputs_markdown)
        call CFG_get(cfg, "control%default_input_markdown_filename", self%default_input_markdown_filename)

        call CFG_get(cfg, "control%run_description", self%run_description)
        call CFG_get(cfg, "control%out_folder", self%out_folder)

        call CFG_get(cfg, "control%calculation", self%calculation)

        ! set the output filename
        if ( trim(adjustl(self%run_description)) == "" ) then

            self%out_filename = trim(adjustl(self%out_folder))//"EXDM_out.hdf5"

        else

            self%out_filename = trim(adjustl(self%out_folder))//"EXDM_out_"//trim(adjustl(self%run_description))//".hdf5"

        end if

    end subroutine


end module

