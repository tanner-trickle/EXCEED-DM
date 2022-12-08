module numerics_TIF_calculator_atomic_type

    use prec_util, only: dp

    implicit none

    type :: numerics_TIF_calculator_atomic_t

        integer :: n_r

        real(dp) :: r_min, r_max

        character(len=512) :: integration_scheme

        contains

            procedure :: set_defaults => numerics_TIF_calculator_atomic_type_set_defaults
            procedure :: get_values => numerics_TIF_calculator_atomic_type_get_values
            procedure :: save => numerics_TIF_calculator_atomic_type_save
            procedure :: generate_integration_grid => numerics_TIF_calculator_atomic_type_generate_integration_grid

    end type

contains

    subroutine numerics_TIF_calculator_atomic_type_generate_integration_grid(self, r_grid, jac)

        use math_util

        implicit none

        class(numerics_TIF_calculator_atomic_t) :: self

        real(dp) :: r_grid(self%n_r)
        real(dp) :: jac(self%n_r)

        integer :: r

        real(dp) :: log_r

        select case ( self%integration_scheme )

        case ( 'log' )

            do r = 1, self%n_r

                log_r = log(self%r_min) + ( log(self%r_max) - log(self%r_min) )*( r - 1 )/max(1.0_dp, self%n_r - 1.0_dp)

                r_grid(r) = exp(log_r) 
                jac(r) = log(self%r_max/self%r_min)*(1.0_dp*self%n_r)**(-1)*r_grid(r)

            end do

        end select

    end subroutine

    subroutine numerics_TIF_calculator_atomic_type_save(self, filename)

        use hdf5_utils

        implicit none

        class(numerics_TIF_calculator_atomic_t) :: self
        character(len=*) :: filename

        integer(HID_T) :: file_id

        call hdf_open_file(file_id, filename, status='OLD', action='WRITE')

        call hdf_create_group(file_id, 'numerics_TIF_calculator_atomic')

        call hdf_write_dataset(file_id, 'numerics_TIF_calculator_atomic/n_r', self%n_r)

        call hdf_write_dataset(file_id, 'numerics_TIF_calculator_atomic/r_min', self%r_min)
        call hdf_write_dataset(file_id, 'numerics_TIF_calculator_atomic/r_max', self%r_min)

        call hdf_write_dataset(file_id, 'numerics_TIF_calculator_atomic/integration_scheme', self%integration_scheme)

        call hdf_close_file(file_id)

    end subroutine

    subroutine numerics_TIF_calculator_atomic_type_set_defaults(self, cfg)

        use m_config

        implicit none

        class(numerics_TIF_calculator_atomic_t) :: self

        type(CFG_t) :: cfg

        call CFG_add(cfg, &
                     "numerics_TIF_calculator_atomic%n_r", &
                     1, &
                     "Number of radial points to integrate with.")

        call CFG_add(cfg, &
                     "numerics_TIF_calculator_atomic%r_min_a0", &
                     1.0e-3_dp, &
                     "Minimum radius to use in the integration, in units of the Bohr radius, $a_0$."//&
                     "<ul>"//&
                     "<li><b>[$a_0$]</b></li>"//&
                     "</ul>")

        call CFG_add(cfg, &
                     "numerics_TIF_calculator_atomic%r_max_a0", &
                     1.0e3_dp, &
                     "Maximum radius to use in the integration, in units of the Bohr radius, $a_0$."//&
                     "<ul>"//&
                     "<li><b>[$a_0$]</b></li>"//&
                     "</ul>")

        call CFG_add(cfg, &
                     "numerics_TIF_calculator_atomic%integration_scheme", &
                     "log", &
                     "Specific method of sampling the radial direction.")

    end subroutine

    subroutine numerics_TIF_calculator_atomic_type_get_values(self, cfg)

        use m_config

        use constants_util

        implicit none

        class(numerics_TIF_calculator_atomic_t) :: self
        type(CFG_t) :: cfg 

        real(dp) :: r_min_a0, r_max_a0

        call CFG_get(cfg, "numerics_TIF_calculator_atomic%integration_scheme", self%integration_scheme)

        call CFG_get(cfg, "numerics_TIF_calculator_atomic%n_r", self%n_r)

        call CFG_get(cfg, "numerics_TIF_calculator_atomic%r_min_a0", r_min_a0)
        call CFG_get(cfg, "numerics_TIF_calculator_atomic%r_max_a0", r_max_a0)

        self%r_min = r_min_a0*a0
        self%r_max = r_max_a0*a0

    end subroutine

end module
