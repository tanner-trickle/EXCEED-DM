module di_grid

    use prec

    use DFT_parameters
    use FFT_util

    implicit none

contains

    subroutine define_q_grid(q_max, k_red_to_xyz, n_q_grid, &
            q_grid_min, n_k_vec, n_FFT_grid, verbose)
        !! q = k' - k + G
        !! 
        !! Goes through all possible q's and finds ones which have |q| < q_max.
        !!
        !! From this list find :
        !!     N_q_grid(i) = N_k(i) x ( max(q^red)_i - min(q^red)_i )
        !!
        !!     q_grid_min = min(q^red)_i
        !!
        !! so that given q_red, the unique index in the q_grid is given by
        !!     index = [1, 1, 1] + n_k x [ q_red - q_grid_min ]

        implicit none

        logical, optional :: verbose

        real(dp) :: q_max
        real(dp) :: k_red_to_xyz(3, 3)

        real(dp) :: q_grid_min(3)
        integer :: n_q_grid(3)

        integer :: k, kp, g1, g2, g3

        integer :: n_k_vec(3)
            !! Might have to be input by hand right now, should add
            !! to DFT_input file
        integer :: n_FFT_grid(3)

        real(dp) :: q_mag
        real(dp) :: q_xyz(3)

        integer :: n_q_less_q_max 
        integer :: ind

        real(dp) :: dk_red(3)
        real(dp) :: dk_red_min(3)
        real(dp) :: dk_red_max(3)

        integer :: G_red(3)
        integer :: G_red_min(3)
        integer :: G_red_max(3)

        real(dp) :: q_red(3)
        real(dp) :: q_red_min(3)
        real(dp) :: q_red_max(3)

        q_red_min = 0.0_dp
        q_red_max = 0.0_dp

        do k = 1, n_k
            do kp = 1, n_k

                dk_red = k_grid_red(kp, :) - k_grid_red(k, :)

                do g3 = 1, n_FFT_grid(3)
                    do g2 = 1, n_FFT_grid(2)
                        do g1 = 1, n_FFT_grid(1)

                            G_red = get_sym_FFT_G_grid_red(n_FFT_grid, [g1, g2, g3], &
                                verbose = .FALSE.)

                            q_red = dk_red + G_red

                            q_xyz = matmul(k_red_to_xyz, q_red)

                            q_mag = norm2(q_xyz)

                            if ( q_mag <= q_max ) then

                                q_red_min(1) = min(q_red(1), q_red_min(1))
                                q_red_min(2) = min(q_red(2), q_red_min(2))
                                q_red_min(3) = min(q_red(3), q_red_min(3))

                                q_red_max(1) = max(q_red(1), q_red_max(1))
                                q_red_max(2) = max(q_red(2), q_red_max(2))
                                q_red_max(3) = max(q_red(3), q_red_max(3))

                            end if

                        end do
                    end do
                end do
            end do
        end do

        n_q_grid(1) = 1 + int( n_k_vec(1)*( q_red_max(1) - q_red_min(1) ) )
        n_q_grid(2) = 1 + int( n_k_vec(2)*( q_red_max(2) - q_red_min(2) ) )
        n_q_grid(3) = 1 + int( n_k_vec(3)*( q_red_max(3) - q_red_min(3) ) )

        q_grid_min = q_red_min

    end subroutine

end module
