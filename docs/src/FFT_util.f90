module FFT_util
    !! Utilities for computing FFT's

    use prec
    use fftw3

    implicit none

    real(dp), allocatable :: sym_FFT_G_grid_xyz(:, :, :, :)
        !! FFT G grid shifted so that both + and - G's are represented

contains

    subroutine find_q_max_FFT(q_max_FFT, k_red_to_xyz, n_grid, verbose)
        !! Finds the maximum \( |q| \) the FFT is consistent for.
        !!
        !! Consider a 2D FFT computed on two square grids in xyz coordinates, \( [N, N] \) and \( [2 N, 2 N] \), ( \( \Delta q = 1
        !! \) ).
        !! Now consider two functions, \( f_1( \mathbf{q} ) \), \( f_2( \mathbf{q} ) \), coming from the FFT of the same function on these grids, binned
        !! in \( |q| \).
        !! \( \bar{f}_1 = \bar{f}_2 \) up to some \( |q| \); in this case \( q = N \), since, e.g., \( |q|= \sqrt{2} N \) components from grid 1 will be
        !! different in grid 2 since grid 2 includes all q with magnitude \( \sqrt(2) N \), whereas grid 1 only contains some of
        !! them (i.e. it's missing \( \mathbf{q} = [\sqrt{2} N, 0] \) ).
        !!
        !! However all grids larger than grid 1 will contain all \( |q| \le N \). This is 'q_max_FFT', or the largest |q| for 
        !! which the FFT results will be identical if the grid size in increased.
        !!
        !! When k_red_to_xyz is not diagonal the problem becomes finding the closest face of the parallelipipid, which is what this
        !! routine computes. 

        implicit none

        real(dp) :: q_max_FFT

        real(dp) :: k_red_to_xyz(3, 3)

        integer :: n_grid(3)

        logical, optional :: verbose

        real(dp) :: corners_xyz(8, 3)

        real(dp) :: faces_xyz(6, 4, 3)

        real(dp) :: dist_to_face(6)

        integer :: f

        real(dp) :: basis_vecs(2, 3)

        real(dp) :: p(3)
        real(dp) :: p_perp(3)

        corners_xyz(1, :) = matmul( k_red_to_xyz, [ -n_grid(1), -n_grid(2), -n_grid(3) ] )/2.0_dp
        corners_xyz(2, :) = matmul( k_red_to_xyz, [ -n_grid(1), -n_grid(2),  n_grid(3) ] )/2.0_dp
        corners_xyz(3, :) = matmul( k_red_to_xyz, [ -n_grid(1),  n_grid(2), -n_grid(3) ] )/2.0_dp
        corners_xyz(4, :) = matmul( k_red_to_xyz, [  n_grid(1), -n_grid(2), -n_grid(3) ] )/2.0_dp
        corners_xyz(5, :) = matmul( k_red_to_xyz, [ -n_grid(1),  n_grid(2),  n_grid(3) ] )/2.0_dp
        corners_xyz(6, :) = matmul( k_red_to_xyz, [  n_grid(1), -n_grid(2),  n_grid(3) ] )/2.0_dp
        corners_xyz(7, :) = matmul( k_red_to_xyz, [  n_grid(1),  n_grid(2), -n_grid(3) ] )/2.0_dp
        corners_xyz(8, :) = matmul( k_red_to_xyz, [  n_grid(1),  n_grid(2),  n_grid(3) ] )/2.0_dp

        ! face 1, [1, 2, 3, 4]
        faces_xyz(1, 1, :) = corners_xyz(1, :)
        faces_xyz(1, 2, :) = corners_xyz(2, :)
        faces_xyz(1, 3, :) = corners_xyz(3, :)
        faces_xyz(1, 4, :) = corners_xyz(4, :)

        ! face 2, [2, 4, 5, 8]
        faces_xyz(2, 1, :) = corners_xyz(2, :)
        faces_xyz(2, 2, :) = corners_xyz(4, :)
        faces_xyz(2, 3, :) = corners_xyz(7, :)
        faces_xyz(2, 4, :) = corners_xyz(8, :)

        ! face 3, [1, 2, 5, 8]
        faces_xyz(3, 1, :) = corners_xyz(1, :)
        faces_xyz(3, 2, :) = corners_xyz(2, :)
        faces_xyz(3, 3, :) = corners_xyz(5, :)
        faces_xyz(3, 4, :) = corners_xyz(8, :)

        ! face 4, [1, 3, 5, 6]
        faces_xyz(4, 1, :) = corners_xyz(1, :)
        faces_xyz(4, 2, :) = corners_xyz(3, :)
        faces_xyz(4, 3, :) = corners_xyz(5, :)
        faces_xyz(4, 4, :) = corners_xyz(6, :)

        ! face 5, [5, 6, 7, 8]
        faces_xyz(5, 1, :) = corners_xyz(5, :)
        faces_xyz(5, 2, :) = corners_xyz(6, :)
        faces_xyz(5, 3, :) = corners_xyz(7, :)
        faces_xyz(5, 4, :) = corners_xyz(8, :)

        ! face 6, [2, 4, 6, 7]
        faces_xyz(6, 1, :) = corners_xyz(2, :)
        faces_xyz(6, 2, :) = corners_xyz(4, :)
        faces_xyz(6, 3, :) = corners_xyz(6, :)
        faces_xyz(6, 4, :) = corners_xyz(7, :)

        do f = 1, 6

            !! shift coordinates to first corner
            p = -faces_xyz(f, 1, :)

            !! define the basis vectors of the plane
            basis_vecs(1, :) = ( faces_xyz(f, 2, :) - faces_xyz(f, 1, :) ) / norm2( faces_xyz(f, 2, :) - faces_xyz(f, 1, :) )
            basis_vecs(2, :) = ( faces_xyz(f, 3, :) - faces_xyz(f, 1, :) ) / norm2( faces_xyz(f, 3, :) - faces_xyz(f, 1, :) )

            !! get magnitude of p perpendicular to face
            p_perp = p - basis_vecs(1, :)*dot_product(basis_vecs(1, :), p) - basis_vecs(2, :)*dot_product(basis_vecs(2, :), p)
            
            dist_to_face(f) = norm2(p_perp)

        end do

        q_max_FFT = minval(dist_to_face)

    end subroutine

    function get_sym_FFT_G_grid_red(n_grid, G_ind_vec, verbose) result(G_red)

        implicit none

        integer :: n_grid(3)

        integer :: G_red(3)
        integer :: G_ind_vec(3)
        integer :: g1, g2, g3

        logical, optional :: verbose

        g1 = G_ind_vec(1)
        g2 = G_ind_vec(2)
        g3 = G_ind_vec(3)

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

    end function

    function get_sym_FFT_G_grid_xyz(n_grid, G_ind_vec, k_red_to_xyz, verbose) result(G_xyz)

        implicit none

        integer :: n_grid(3)

        integer :: G_red(3)
        integer :: G_ind_vec(3)
        integer :: g1, g2, g3

        real(dp) :: k_red_to_xyz(3, 3)

        real(dp) :: G_xyz(3)

        logical, optional :: verbose

        g1 = G_ind_vec(1)
        g2 = G_ind_vec(2)
        g3 = G_ind_vec(3)

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

        G_xyz = matmul(k_red_to_xyz, G_red)

    end function

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
