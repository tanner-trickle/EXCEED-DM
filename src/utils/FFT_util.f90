module FFT_util
    ! Collection of useful functions when computing the FFT numerically. 
    !
    ! .. note: Only depends on FFTW3.

    use fftw3

    use prec_util, only: dp

    implicit none

contains

    subroutine create_FFT_G_vec_array(G_vec_array, n_FFT_grid, k_red_to_xyz)
        ! Generates an array of G vectors corresponding to an FFT grid

        implicit none

        integer :: n_FFT_grid(3)
        real(dp) :: k_red_to_xyz(3, 3)

        real(dp) :: G_vec_array(n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3), 3)

        integer :: g1, g2, g3

        integer :: G_red(3)

        do g3 = 1, n_FFT_grid(3)
            do g2 = 1, n_FFT_grid(2)
                do g1 = 1, n_FFT_grid(1)

                    G_red = FFT_idx_to_G_red([g1, g2, g3], n_FFT_grid)

                    G_vec_array(g1, g2, g3, :) = matmul(k_red_to_xyz, G_red)

                end do
            end do
        end do

    end subroutine

    subroutine create_FFT_G_vec_list(G_vec_list, n_FFT_grid, k_red_to_xyz)
        ! Generates the list of G vectors corresponding to an FFT grid

        implicit none

        real(dp) :: G_vec_list(:, :)
        integer :: n_FFT_grid(3)
        real(dp) :: k_red_to_xyz(3, 3)

        real(dp) :: G_vec_array(n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3), 3)

        call create_FFT_G_vec_array(G_vec_array, n_FFT_grid, k_red_to_xyz)

        G_vec_list = reshape( G_vec_array, [ size(G_vec_list, 1), 3 ] )

    end subroutine

    function FFT_idx_to_G_red(FFT_idx, n_grid) result(G_red)
        ! Converts an FFT index to a G vector in reduced coordinates.

        implicit none

        integer :: FFT_idx(3)
        integer :: n_grid(3)
        integer :: G_red(3)

        integer :: d

        do d = 1, 3

            if ( mod(N_grid(d), 2) == 0 ) then

                if ( FFT_idx(d) <= N_grid(d)/2 ) then

                    G_red(d) = FFT_idx(d) - 1

                else

                    G_red(d) = FFT_idx(d) - N_grid(d) - 1

                end if

            else

                if ( FFT_idx(d) <= (N_grid(d) + 1)/2 ) then

                    G_red(d) = FFT_idx(d) - 1

                else

                    G_red(d) = FFT_idx(d) - N_grid(d) - 1

                end if

            end if

        end do

    end function

    function G_red_to_FFT_idx(G_red, n_grid) result ( FFT_idx )

        implicit none

        integer :: G_red(3)
        integer :: n_grid(3)
        integer :: FFT_idx(3)

        integer :: i

        FFT_idx = G_red + 1

        do i = 1, 3
            if ( G_red(i) < 0 ) then
                FFT_idx(i) = G_red(i) + n_grid(i) + 1
            end if
        end do

    end function

    subroutine create_FFT_plan_pair(n_grid, FFT_plan_pair)
        ! Creates a forwards/backwards pair of FFT plans for FFTW3.

        implicit none

        integer :: n_grid(3)

        integer :: FFT_plan_pair(2, 8)

        call set_fft_plan_forward_3d(n_grid, FFT_plan_pair(1, :))
        call set_fft_plan_backward_3d(n_grid, FFT_plan_pair(2, :))

    end subroutine

    subroutine set_fft_plan_forward_3d(n_grid, fft_plan)
        ! Sets the FFTW3 FFT plan with a 'forwards' convention. See
        ! [https://www.fftw.org/fftw3_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html#The-1d-Discrete-Fourier-Transform-_0028DFT_0029](https://www.fftw.org/fftw3_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html#The-1d-Discrete-Fourier-Transform-_0028DFT_0029).
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
        ! Sets the FFTW3 FFT plan with a 'backwards' convention. See 
        ! [https://www.fftw.org/fftw3_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html#The-1d-Discrete-Fourier-Transform-_0028DFT_0029](https://www.fftw.org/fftw3_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html#The-1d-Discrete-Fourier-Transform-_0028DFT_0029).
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

    subroutine zero_pad_FFT_matrix(mat_in, mat_out, FFT_plan_1, FFT_plan_2)

        implicit none

        complex(dp) :: mat_in(:, :, :)
        complex(dp) :: mat_out(:, :, :)

        complex(dp) :: mat_in_FFT(size(mat_in, 1), size(mat_in, 2), size(mat_in, 3))

        complex(dp) :: mat_in_FFT_zero_pad(size(mat_out, 1), size(mat_out, 2), size(mat_out, 3))

        integer :: FFT_idx(3)
        integer :: G_red(3)

        integer :: FFT_plan_1(8)
        integer :: FFT_plan_2(8)

        integer :: g1, g2, g3

        ! FFT mat_in
        call dfftw_execute_dft(FFT_plan_1, mat_in, mat_in_FFT) 

        ! normalize
        mat_in_FFT = (1.0_dp*product(shape(mat_in)))**(-1)*mat_in_FFT

        ! zero pad
        mat_in_FFT_zero_pad = (0.0_dp, 0.0_dp)

        do g3 = 1, size(mat_in_FFT, 3)
            do g2 = 1, size(mat_in_FFT, 2)
                do g1 = 1, size(mat_in_FFT, 1)

                    ! G_red = G_red_to_sym_G_red([ g1 - 1, g2 - 1, g3 - 1 ], [size(mat_in, 1), size(mat_in, 2), size(mat_in, 3)])
                    G_red = FFT_idx_to_G_red([g1, g2, g3], [size(mat_in, 1), size(mat_in, 2), size(mat_in, 3) ])

                    FFT_idx = G_red_to_FFT_idx(G_red, [size(mat_out, 1), size(mat_out, 2), size(mat_out, 3)])

                    mat_in_FFT_zero_pad(FFT_idx(1), FFT_idx(2), FFT_idx(3)) = mat_in_FFT(g1, g2, g3) 

                end do
            end do
        end do

        ! FFT back
        call dfftw_execute_dft(FFT_plan_2, mat_in_FFT_zero_pad, mat_out) 

    end subroutine

end module
