module experiment_type

    use prec_util, only: dp

    implicit none

    type :: experiment_t

        real(dp) :: M_kg
        real(dp) :: T_year

        real(dp) :: M 
        real(dp) :: T

        character(len=512) :: name = ''

        contains

            procedure :: set_defaults => experiment_type_set_defaults
            procedure :: get_values => experiment_type_get_values
            procedure :: save => experiment_type_save

    end type

contains

    subroutine experiment_type_save(self, filename)

        use hdf5_utils

        implicit none

        class(experiment_t) :: self
        character(len=*) :: filename

        integer(HID_T) :: file_id

        call hdf_open_file(file_id, filename, status='OLD', action='WRITE')

        call hdf_create_group(file_id, 'experiment')

        call hdf_write_dataset(file_id, 'experiment/M', self%M_kg)
        call hdf_write_dataset(file_id, 'experiment/T', self%T_year)

        call hdf_close_file(file_id)

    end subroutine

    subroutine experiment_type_set_defaults(self, cfg)

        use m_config

        implicit none

        class(experiment_t) :: self

        type(CFG_t) :: cfg

        call CFG_add(cfg, &
                        "experiment%M_kg", &
                        1.0_dp, &
                        "Mass of the experimental target <br />Units: kg")

        call CFG_add(cfg, &
                        "experiment%T_year", &
                        1.0_dp, &
                        "Exposure time of the experimental target <br />Units: yr")

    end subroutine

    subroutine experiment_type_get_values(self, cfg)

        use m_config

        use constants_util

        implicit none

        class(experiment_t) :: self

        type(CFG_t) :: cfg

        call CFG_get(cfg, "experiment%M_kg", self%M_kg)
        call CFG_get(cfg, "experiment%T_year", self%T_year)

        self%M = self%M_kg*kg_to_eV
        self%T = self%T_year*yr_to_inv_eV

    end subroutine

end module
