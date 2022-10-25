module numerics_binned_scatter_rate_type

    use prec_util, only: dp

    implicit none

    type :: numerics_binned_scatter_rate_t

        integer :: n_q_bins
        integer :: n_E_bins

        real(dp) :: E_bin_width
        real(dp) :: q_bin_width

        contains

            procedure :: set_defaults => numerics_binned_scatter_rate_type_set_defaults
            procedure :: get_values => numerics_binned_scatter_rate_type_get_values
            procedure :: save => numerics_binned_scatter_rate_type_save

    end type

contains

    subroutine numerics_binned_scatter_rate_type_save(self, filename)

        use hdf5_utils

        implicit none

        class(numerics_binned_scatter_rate_t) :: self
        character(len=*) :: filename

        integer(HID_T) :: file_id

        call hdf_open_file(file_id, filename, status='OLD', action='WRITE')

        call hdf_create_group(file_id, 'numerics_binned_scatter_rate')

        call hdf_write_dataset(file_id, 'numerics_binned_scatter_rate/E_bin_width', self%E_bin_width)
        call hdf_write_dataset(file_id, 'numerics_binned_scatter_rate/q_bin_width', self%q_bin_width)

        call hdf_close_file(file_id)

    end subroutine

    subroutine numerics_binned_scatter_rate_type_set_defaults(self, cfg)

        use m_config

        implicit none

        class(numerics_binned_scatter_rate_t) :: self
        type(CFG_t) :: cfg

        call CFG_add(cfg, &
                    "numerics_binned_scatter_rate%n_q_bins", &
                    1, &
                    "Number of bins in $q$ space")

        call CFG_add(cfg, &
                    "numerics_binned_scatter_rate%n_E_bins", &
                    1, &
                    "Number of bins in $\omega$ space")

        call CFG_add(cfg, &
                    "numerics_binned_scatter_rate%q_bin_width", &
                    1.0_dp, &
                    "Width of bins in $q$ space<br />"//&
                     "<ul>"//&
                     "<li><b>Units</b>: $\text{keV}$</li>"//&
                     "</ul>")

        call CFG_add(cfg, &
                    "numerics_binned_scatter_rate%E_bin_width", &
                    1.0_dp, &
                    "Width of bins in $\omega$ space<br />"//&
                     "<ul>"//&
                     "<li><b>Units</b>: $\text{eV}$</li>"//&
                     "</ul>")

    end subroutine

    subroutine numerics_binned_scatter_rate_type_get_values(self, cfg)

        use m_config

        implicit none

        class(numerics_binned_scatter_rate_t) :: self
        type(CFG_t) :: cfg

        call CFG_get(cfg, "numerics_binned_scatter_rate%n_q_bins", self%n_q_bins)
        call CFG_get(cfg, "numerics_binned_scatter_rate%n_E_bins", self%n_E_bins)
        call CFG_get(cfg, "numerics_binned_scatter_rate%q_bin_width", self%q_bin_width)
        call CFG_get(cfg, "numerics_binned_scatter_rate%E_bin_width", self%E_bin_width)

    end subroutine

end module
