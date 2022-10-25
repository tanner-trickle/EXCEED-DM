module numerics_absorption_rate_type

    use prec_util, only: dp

    implicit none

    type :: numerics_absorption_rate_t

        real(dp), allocatable :: widths(:, :)

        character(len=512) :: smear_type

        contains

            procedure :: set_defaults => numerics_absorption_rate_type_set_defaults
            procedure :: get_values => numerics_absorption_rate_type_get_values
            procedure :: save => numerics_absorption_rate_type_save

    end type

contains

    subroutine numerics_absorption_rate_type_save(self, filename)

        use hdf5_utils

        implicit none

        class(numerics_absorption_rate_t) :: self
        character(len=*) :: filename

        integer(HID_T) :: file_id

        call hdf_open_file(file_id, filename, status='OLD', action='WRITE')

        call hdf_create_group(file_id, 'numerics_absorption_rate')

        call hdf_write_dataset(file_id, 'numerics_absorption_rate/widths', self%widths)
        call hdf_write_dataset(file_id, 'numerics_absorption_rate/smear_type', self%smear_type)

        call hdf_close_file(file_id)

    end subroutine

    subroutine numerics_absorption_rate_type_set_defaults(self, cfg)

        use m_config

        implicit none

        class(numerics_absorption_rate_t) :: self

        type(CFG_t) :: cfg

        call CFG_add(cfg,&
                     "numerics_absorption_rate%widths", &
                     [0.0_dp, 1.0e-1_dp, 100.0_dp], &
                     "List of widths, $\delta \, [\text{eV}]$, to compute for, parameterized as [$a$, $b$, $c$]<br />"//&
                     "<ul>"//&
                     "<li><b>Formula</b>: $\delta = \text{min}(a + b \omega, c)$</li>"//&
                     "<li><b>Units</b>: [$\text{eV}$, -, $\text{eV}$]</li>"//&
                     "<li><b>Dim</b>: [ : , 3]</li>"//&
                     "</ul>", dynamic_size = .TRUE.)

        call CFG_add(cfg,&
                     "numerics_absorption_rate%smear_type", &
                     "lorentz", &
                     "Defines broadening behavior for the imaginary part of the Greens function")

    end subroutine

    subroutine numerics_absorption_rate_type_get_values(self, cfg)

        use m_config

        implicit none

        class(numerics_absorption_rate_t) :: self
        type(CFG_t) :: cfg

        integer :: n 

        real(dp), allocatable :: widths(:)

        call CFG_get_size(cfg, "numerics_absorption_rate%widths", n)

        allocate(widths(n))
        allocate(self%widths(n/3, 3), source = 0.0_dp)

        call CFG_get(cfg, "numerics_absorption_rate%widths", widths)
        self%widths = transpose(reshape( widths, [ 3, n/3 ] ) )

        call CFG_get(cfg, "numerics_absorption_rate%smear_type", self%smear_type)

    end subroutine

end module
