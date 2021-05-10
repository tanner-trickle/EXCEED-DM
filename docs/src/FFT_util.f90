module FFT_util
    !! Utilities for computing FFT's

    use prec
    use fftw3

    implicit none

    real(dp), allocatable :: sym_FFT_G_grid_xyz(:, :, :, :)
        !! FFT G grid shifted so that both + and - G's are represented

    contains

    subroutine set_sym_FFT_G_grid_xyz(n_grid, k_red_to_xyz, verbose)
        !! Symmetric FFT grid. 
        !! 
        !! Standard convention for FFTs is to compute for frequencies 0 -> N-1. 
        !! This subroutine creates a map which takes the larger half of the positive
        !! frequencies and maps then back to negative values, for a symmetric G grid.  
        implicit none

        real(dp) :: k_red_to_xyz(3, 3)
        logical, optional :: verbose

        integer :: g1, g2, g3
        integer :: G_red(3)

        integer :: n_grid(3)

        allocate(sym_FFT_G_grid_xyz(n_grid(1), n_grid(2), n_grid(3), 3))

        do g3 = 1, n_grid(3)
            do g2 = 1, n_grid(2)
                do g1 = 1, n_grid(1)

                    G_red = [g1 - 1, g2 - 1, g3 - 1]

                    if ( g1 .gt. n_grid(1)/2 ) then
                        G_red(1) = g1 - n_grid(1) - 1
                    end if

                    if ( g2 .gt. n_grid(2)/2 ) then
                        G_red(2) = g2 - n_grid(2) - 1
                    end if

                    if ( g3 .gt. n_grid(3)/2 ) then
                        G_red(3) = g3 - n_grid(3) - 1
                    end if

                    sym_FFT_G_grid_xyz(g1, g2, g3, :) = matmul(k_red_to_xyz, G_red)

                end do
            end do
        end do

    end subroutine

    subroutine set_fft_plan_forward_3d(n_grid, fft_plan)
        !! sets the FFT plan for an array with dimensinos [ n_grid ]
        implicit none

        integer :: fft_plan(8)
        integer :: n_grid(3)

        complex(dp) :: test_mat_in(n_grid(1), n_grid(2), n_grid(3))
        complex(dp) :: test_mat_out(n_grid(1), n_grid(2), n_grid(3))

        test_mat_in = (0.0_dp, 0.0_dp)
        test_mat_out = (0.0_dp, 0.0_dp)

        call dfftw_plan_dft_3d(fft_plan, n_grid(3), n_grid(2), n_grid(1),&
                                test_mat_in, test_mat_out, FFTW_FORWARD, FFTW_PATIENT)

    end subroutine

    subroutine set_fft_plan_backward_3d(n_grid, fft_plan)
        !! sets the FFT plan for an array with dimensinos [ n_grid ]
        implicit none

        integer :: fft_plan(8)
        integer :: n_grid(3)

        complex(dp) :: test_mat_in(n_grid(1), n_grid(2), n_grid(3))
        complex(dp) :: test_mat_out(n_grid(1), n_grid(2), n_grid(3))

        test_mat_in = (0.0_dp, 0.0_dp)
        test_mat_out = (0.0_dp, 0.0_dp)

        call dfftw_plan_dft_3d(fft_plan, n_grid(3), n_grid(2), n_grid(1),&
                                test_mat_in, test_mat_out, FFTW_BACKWARD, FFTW_PATIENT)

    end subroutine

    subroutine G_red_to_FFT_G_grid_index(n_grid, G_red, idx)
        !! Given a G_red vector returns the index of that vector in the FFT 
        !! grid.
        implicit none

        integer :: idx(3)
        integer :: G_red(3)

        integer :: n_grid(3)

        integer :: i

        idx = G_red + [1, 1, 1]

        ! put negative values at the back of the array in FFTW3 convention
        do i = 1, 3
            if ( G_red(i) .lt. 0 ) then
                idx(i) = G_red(i) + n_grid(i) + 1
            end if
        end do

    end subroutine

end module
