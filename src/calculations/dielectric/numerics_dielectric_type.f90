module numerics_dielectric_type

    use prec_util, only: dp

    implicit none

    type :: numerics_dielectric_t

        integer :: n_q_bins
        integer :: n_E_bins

        real(dp) :: E_bin_width
        real(dp) :: q_bin_width

        integer :: n_q_theta
        integer :: n_q_phi

        real(dp), allocatable :: widths(:, :)

        character(len=512) :: smear_type

        contains

            procedure :: set_defaults => numerics_dielectric_type_set_defaults
            procedure :: get_values => numerics_dielectric_type_get_values
            procedure :: save => numerics_dielectric_type_save

    end type

contains

    subroutine numerics_dielectric_type_set_defaults(self, cfg)

        use m_config

        implicit none

        class(numerics_dielectric_t) :: self
        type(CFG_t) :: cfg

        call CFG_add(cfg, &
                    "numerics_dielectric%n_q_bins", &
                    1, &
                    "Number of bins in $q$ space")

        call CFG_add(cfg, &
                    "numerics_dielectric%n_E_bins", &
                    1, &
                    "Number of bins in $\omega$ space")

        call CFG_add(cfg, &
                    "numerics_dielectric%n_q_theta", &
                    1, &
                    "Number of bins in $\theta_\mathbf{q}$")

        call CFG_add(cfg, &
                    "numerics_dielectric%n_q_phi", &
                    1, &
                    "Number of bins in $\phi_\mathbf{q}$")

        call CFG_add(cfg, &
                    "numerics_dielectric%q_bin_width", &
                    1.0_dp, &
                    "Width of bins in $q$ space<br />"//&
                     "<ul>"//&
                     "<li><b>Units</b>: $\text{keV}$</li>"//&
                     "</ul>")

        call CFG_add(cfg, &
                    "numerics_dielectric%E_bin_width", &
                    1.0_dp, &
                    "Width of bins in $\omega$ space<br />"//&
                     "<ul>"//&
                     "<li><b>Units</b>: $\text{eV}$</li>"//&
                     "</ul>")

        call CFG_add(cfg,&
                     "numerics_dielectric%widths", &
                     [0.0_dp, 1.0e-1_dp, 100.0_dp], &
                     "List of widths, $\delta \, [\text{eV}]$, to compute for, parameterized as [$a$, $b$, $c$]<br />"//&
                     "<ul>"//&
                     "<li><b>Formula</b>: $\delta = \text{min}(a + b \omega, c)$</li>"//&
                     "<li><b>Units</b>: [$\text{eV}$, -, $\text{eV}$]</li>"//&
                     "<li><b>Dim</b>: [ : , 3]</li>"//&
                     "</ul>", dynamic_size = .TRUE.)

        call CFG_add(cfg,&
                     "numerics_dielectric%smear_type", &
                     "lorentz", &
                     "Defines broadening behavior for the imaginary part of the Greens function.")

    end subroutine

    subroutine numerics_dielectric_type_get_values(self, cfg)

        use m_config

        implicit none

        class(numerics_dielectric_t) :: self
        type(CFG_t) :: cfg

        integer :: n 
        real(dp), allocatable :: widths(:)

        call CFG_get(cfg, "numerics_dielectric%n_q_bins", self%n_q_bins)
        call CFG_get(cfg, "numerics_dielectric%n_q_theta", self%n_q_theta)
        call CFG_get(cfg, "numerics_dielectric%n_q_phi", self%n_q_phi)
        call CFG_get(cfg, "numerics_dielectric%n_E_bins", self%n_E_bins)
        call CFG_get(cfg, "numerics_dielectric%q_bin_width", self%q_bin_width)
        call CFG_get(cfg, "numerics_dielectric%E_bin_width", self%E_bin_width)

        call CFG_get_size(cfg, "numerics_dielectric%widths", n)

        allocate(widths(n))
        allocate(self%widths(n/3, 3), source = 0.0_dp)

        call CFG_get(cfg, "numerics_dielectric%widths", widths)
        self%widths = transpose(reshape( widths, [ 3, n/3 ] ) )

        call CFG_get(cfg, "numerics_dielectric%smear_type", self%smear_type)

    end subroutine

    subroutine numerics_dielectric_type_save(self, filename)

        use hdf5_utils

        implicit none

        class(numerics_dielectric_t) :: self
        character(len=*) :: filename

        integer(HID_T) :: file_id

        call hdf_open_file(file_id, filename, status='OLD', action='WRITE')

        call hdf_create_group(file_id, 'numerics_dielectric')

        call hdf_write_dataset(file_id, 'numerics_dielectric/E_bin_width', self%E_bin_width)
        call hdf_write_dataset(file_id, 'numerics_dielectric/q_bin_width', self%q_bin_width)

        call hdf_close_file(file_id)

    end subroutine

end module
