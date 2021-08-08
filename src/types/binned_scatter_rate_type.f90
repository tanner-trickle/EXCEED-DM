module binned_scatter_rate_type
    !! Defines the `binned_scatter_rate` data type.

    use prec

    use dm_model_type
    use expt_type
    use bins_scatter_type

    implicit none

    type binned_scatter_rate_t
        !! Binned scattering rate for each DM mass, mediator form factor, and time of day. Binning is
        !! in momentum deposition, \( q \)  and energy deposition, \( \omega \) (referred to as `E` in the code).
        real(dp), allocatable :: binned_rate(:, :, :, :, :)
            !! Dim : [ bins%n_q, bins%n_E, dm_model%n_mX, dm_model%n_med_FF, expt%n_time ]
            !!
            !! Binned scattering rate per cross section, \( \Delta R_{q, E}(m_\chi, t, \beta) \).
            !!
            !! Units : \( \text{cm}^{-2} \)

        contains

            procedure :: init => binned_scatter_rate_init
            procedure :: compute_rate => binned_scatter_rate_compute_rate
            ! procedure :: save => binned_scatter_rate_save

    end type

contains

    function binned_scatter_rate_compute_rate(self) result( rate )
        !! Sums `binned_rate%binned_rate` over the \( q, E \) bins to get the total rate.
        implicit none

        class(binned_scatter_rate_t) :: self
        integer :: q, E
        real(dp) :: rate(&
            size(self%binned_rate, 3),&
            size(self%binned_rate, 4),&
            size(self%binned_rate, 5))

        rate = 0.0_dp

        do E = 1, size(self%binned_rate, 2)
            do q = 1, size(self%binned_rate, 1)

                rate = rate + self%binned_rate(q, E, :, :, :)

            end do
        end do

    end function

    subroutine binned_scatter_rate_init(self, bins, dm_model, expt)
        !! Creates a `binned_rate` object and sets `binned_rate%binned_rate` = 0.

        implicit none

        class(binned_scatter_rate_t) :: self
        type(bins_scatter_t) :: bins
        type(dm_model_t) :: dm_model
        type(expt_t) :: expt

        allocate(self%binned_rate( bins%n_q, bins%n_E, &
            dm_model%n_mX, dm_model%n_med_FF, expt%n_time))
        self%binned_rate = 0.0_dp

    end subroutine

end module
