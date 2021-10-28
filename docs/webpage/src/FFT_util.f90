module FFT_util
    !! Utilities for computing FFT's.

    use fftw3

    use prec
    use info_messages

    implicit none

    type FFT_grid_t
        !! Collection of parameters defining the grid the FFT is performed on, 
        !! and how the FFT should be performed.
        integer :: N
            !! Number of points in the FFT
        integer :: n_grid(3)
            !! Number of points in the FFT in each dimension
        real(dp) :: q_max
            !! Maximum q for which the FFT is consistent for. See `find_q_max`
            !!
            !! Units : eV
        integer :: plan(8)
            !! The FFTW3 plan
        integer, allocatable :: sym_G_grid_red(:, :, :, :)
            !! Dim : [ n_grid(1), n_grid(2), n_grid(3), 3 ]
            !!
            !! Shifted, symmetric list of G points in FFT grid in reduced coordinates
            !!
            !! Units : None
        real(dp), allocatable :: sym_G_grid_xyz(:, :, :, :)
            !! Dim : [ n_grid(1), n_grid(2), n_grid(3), 3 ]
            !!
            !! Shifted, symmetric list of G points in FFT grid in xyz coordinates
            !!
            !! Units : eV

        contains

            procedure :: init => FFT_grid_init
            procedure :: save => FFT_grid_save
            procedure :: print => FFT_grid_print
            procedure :: find_q_max
            procedure :: set_sym_G_grid => FFT_grid_set_sym_G_grid

    end type

contains

    subroutine FFT_grid_print(self, verbose)
        !! Prints `FFT_grid` components.

        implicit none

        class(FFT_grid_t) :: self
        logical, optional :: verbose

        if ( verbose ) then
            call print_section_seperator()
            print*, '    --------'
            print*, '    FFT Grid'
            print*, '    --------'
            print*
            print*, '        Number of points : ', self%n
            print*, '        Dimensions       : '
            print*, '        ', self%n_grid
            print*
            print*, '        q_max = ', self%q_max/1.0e3_dp, 'keV'
            print*
            call print_section_seperator()
            print*
        end if

    end subroutine

    subroutine FFT_grid_save(self, filename, verbose)
        !! Saves `FFT_grid`.

        use hdf5
        use h5lt

        implicit none

        class(FFT_grid_t) :: self
        character(len=*) :: filename
        logical, optional :: verbose

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id
        logical :: file_exists
        integer(HSIZE_T) :: dims1(1) = [1]
        integer(HSIZE_T) :: dims2(1) = [3]
        integer :: error

        if ( verbose ) then
            print*, 'Saving FFT grid parameters...'
            print*
        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'FFT', group_id, error)

            ! write data
            call h5ltmake_dataset_int_f(file_id, &
                'FFT/N', &
                size(dims1), dims1, &
                self%N, &
                error)
            call h5ltmake_dataset_int_f(file_id, &
                'FFT/n_grid', &
                size(dims2), dims2, &
                self%n_grid, &
                error)
            call h5ltmake_dataset_double_f(file_id, &
                'FFT/q_max', &
                size(dims1), dims1, &
                self%q_max, &
                error)

            call h5fclose_f(file_id, error)
            call h5close_f(error)

        else

            call print_error_message(&
                'Output file : '//trim(filename)//' does NOT exist.', &
                verbose = verbose)
            stop

        end if

    end subroutine

    subroutine FFT_grid_init(self, n_grid, k_red_to_xyz, fb_char, verbose)
        !! Initialize an `FFT_grid` instance. 
        !!
        !! 1) Find \( q^\text{FFT}_\text{max} \)
        !! 2) Initialize the FFT plan. (Default: FFTW3 'patient' planning)
        !! 3) Set the symmetric FFT grids.
        implicit none

        class(FFT_grid_t) :: self
        integer :: n_grid(3)
        real(dp) :: k_red_to_xyz(3, 3)
        character(len=1) :: fb_char
            !! 'f' for forward rule, 'b' for backwards
        logical, optional :: verbose

        self%n_grid = n_grid
        self%n = self%n_grid(1)*self%n_grid(2)*self%n_grid(3)

        call self%find_q_max(k_red_to_xyz, verbose = verbose)

        if ( verbose ) then

            print*, 'Finding optimal FFT plan...'
            print*

        end if

        if ( fb_char == 'b' ) then

            call set_fft_plan_backward_3d(self%n_grid, self%plan)

        else

            call set_fft_plan_forward_3d(self%n_grid, self%plan)

        end if

        if ( verbose ) then

            print*, 'Setting symmetric FFT G grid...'
            print*

        end if

        call self%set_sym_G_grid(k_red_to_xyz, verbose = verbose)

    end subroutine

    subroutine FFT_grid_set_sym_G_grid(self, k_red_to_xyz, verbose)
        !! Symmetric FFT grid. 
        !! 
        !! Standard convention for FFTs is to compute for frequencies 0 -> N-1. 
        !! This subroutine creates a map which takes the larger half of the positive
        !! frequencies and maps then back to negative values, for a symmetric G grid.  
        implicit none

        class(FFT_grid_t) :: self
        real(dp) :: k_red_to_xyz(3, 3)
        logical, optional :: verbose

        integer :: g1, g2, g3
        integer :: G_red(3)

        allocate(self%sym_G_grid_red(self%n_grid(1), self%n_grid(2), self%n_grid(3), 3))
        allocate(self%sym_G_grid_xyz(self%n_grid(1), self%n_grid(2), self%n_grid(3), 3))

        do g3 = 1, self%n_grid(3)
            do g2 = 1, self%n_grid(2)
                do g1 = 1, self%n_grid(1)

                    G_red = [g1 - 1, g2 - 1, g3 - 1]

                    if ( g1 > self%n_grid(1)/2 ) then
                        G_red(1) = g1 - self%n_grid(1) - 1
                    end if

                    if ( g2 > self%n_grid(2)/2 ) then
                        G_red(2) = g2 - self%n_grid(2) - 1
                    end if

                    if ( g3 > self%n_grid(3)/2 ) then
                        G_red(3) = g3 - self%n_grid(3) - 1
                    end if

                    self%sym_G_grid_red(g1, g2, g3, :) = G_red 
                    self%sym_G_grid_xyz(g1, g2, g3, :) = matmul(k_red_to_xyz, G_red)

                end do
            end do
        end do

    end subroutine

    subroutine find_q_max(self, k_red_to_xyz, verbose)
        !* Finds the maximum \( q \) the FFT is consistent for, `q_max`.
        !
        ! Consider a 2D FFT computed on two square grids in xyz coordinates of dimension, 1) \( [N, N] \) and 2) \( [2 N, 2 N] \),
        ! so the second grid includes all of the points in the first grid. Now consider two functions, \( f_1( \mathbf{q} ) \), \(
        ! f_2( \mathbf{q} ) \), coming from the FFT of the same function on these grids, binned in \( q \). \( \bar{f}_1 \) will
        ! equal \( \bar{f}_2 \) up to some \( q \); in this case \( q = N \), since, e.g., the \( q = \sqrt{2} N \) value of \(
        ! \bar{f}_1 \) will be different in grid 2 since grid 2 includes all \( \mathbf{q} \) with magnitude \( \sqrt{2} N \),
        ! whereas grid 1 only contains some of them (i.e., it's missing \( \mathbf{q} = [\sqrt{2} N, 0] \) ).
        !
        ! However all grids larger than grid 1 will contain all \( q \le N \). This is the `q_max` computed here, or the largest \(
        ! q \) for which the FFT results will be identical if the grid size in increased.
        !
        ! When `k_red_to_xyz` is not diagonal the problem becomes finding the closest face of the parallelipipid, which is (more
        ! practically) what this routine computes. 

        implicit none

        class(FFT_grid_t) :: self
        real(dp) :: k_red_to_xyz(3, 3)
        logical, optional :: verbose

        real(dp) :: corners_xyz(8, 3)

        real(dp) :: faces_xyz(6, 4, 3)

        real(dp) :: dist_to_face(6)

        integer :: f

        real(dp) :: basis_vecs(2, 3)

        real(dp) :: p(3)
        real(dp) :: p_perp(3)

        integer :: n_grid(3)

        if ( verbose ) then

            print*, 'Finding q_max for FFT...'
            print*

        end if

        n_grid = self%n_grid

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

        self%q_max = minval(dist_to_face)

    end subroutine

    subroutine set_fft_plan_forward_3d(n_grid, fft_plan)
        !! Sets the FFTW3 FFT plan with a 'forwards' convention. See
        !![https://www.fftw.org/fftw3_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html#The-1d-Discrete-Fourier-Transform-_0028DFT_0029](https://www.fftw.org/fftw3_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html#The-1d-Discrete-Fourier-Transform-_0028DFT_0029).
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
        !! Sets the FFTW3 FFT plan with a 'backwards' convention. See
        !![https://www.fftw.org/fftw3_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html#The-1d-Discrete-Fourier-Transform-_0028DFT_0029](https://www.fftw.org/fftw3_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html#The-1d-Discrete-Fourier-Transform-_0028DFT_0029).
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
            if ( G_red(i) < 0 ) then
                idx(i) = G_red(i) + n_grid(i) + 1
            end if
        end do

    end subroutine

end module
