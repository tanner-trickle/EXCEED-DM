module PW_dataset_type
    !! Defines the `PW_dataset` type.
   
    use h5lt
    use hdf5

    use prec
    use constants
    use units
    use info_messages

    implicit none

    type PW_dataset_t
        !! Stores the metadata of the Bloch wave function coefficients in the Fourier basis, 
        !!\( \widetilde{u}_{i, \mathbf{k}, \mathbf{G}} \), e.g. number of \( \mathbf{k} \) points, \( N_\mathbf{k} \)
        !! , number of bands, etc and has routines which load the coefficients from an external data file when needed.
        !!
        !! $$\begin{align}
        !! \psi_{i, \mathbf{k}}(\mathbf{x}) = \frac{1}{\sqrt{V}} 
        !! \sum_\mathbf{G} \widetilde{u}_{i, \mathbf{k}, \mathbf{G}} e^{i (\mathbf{k} + \mathbf{G}) \cdot \mathbf{x}}
        !! \nonumber
        !! \end{align}$$.
        integer :: n_k
            !! Number of k points
        integer :: n_val
            !! Number of valence bands
        integer :: n_cond
            !! Number of conduction bands
        integer :: n_G
            !! Number of G points
        integer, allocatable :: G_grid_red(:, :)
            !! Dim : [n_G, 3]
            !!
            !! G vectors in reduced coordinates.
            !!
            !! Note : This G grid is k-independent and assumes that the \( \widetilde{u} \)'s,
            !! for each \( \mathbf{k} \), are defined on the same grid. 
            !!
            !! Units : None
        real(dp) :: a_vecs_A(3, 3)
            !! Primitive lattice vectors
            !!
            !! (i, :) is the ith primitive lattice vector
            !! 
            !! Units : Ang
        real(dp) :: b_vecs_A(3, 3)
            !! Reciprocal lattice vector
            !! 
            !! (i, :) is the ith reciprocal lattice vector
            !! 
            !! Units : Ang^(-1)
        real(dp), allocatable :: k_weight(:)
            !! Dim : [n_k]
            !! 
            !! Weights of each k point in the k_red grid, must sum to 2
            !!
            !! Units : None
        real(dp), allocatable :: k_grid_red(:, :)
            !! Dim : [n_k, 3]
            !!
            !! List of k vectors in reduced coordinates
            !! 
            !! Units : None
        real(dp), allocatable :: energy_bands_raw(:, :)
            !! Dim : [n_k, n_bands]
            !!
            !! Electron energy eigenvalues taken straight from the input file, no scissor 
            !! correction performed.
            !!
            !! Units : eV 
        real(dp) :: Ef_max = 60.0_dp
            !! Maximum electron energy computed for.
            !!
            !! Units : eV
        ! generated from input
        integer :: n_bands
            !! Total number of bands
            !! 
            !! n_bands = n_cond + n_val
        real(dp) :: q_cut
            !! Plane wave expansion cutoff
            !!
            !! \( \text{max} (| \mathbf{k} + \mathbf{G} |) \) 
            !!
            !! Units : eV
        real(dp) :: E_cut
            !! Energy cooresponding to q_cut
            !!
            !! E_cut = q_cut**2/(2*m_elec)
            !!
            !! Units : eV
        real(dp) :: a_vecs(3, 3)
            !! Primitive lattice vectors
            !!
            !! (i, :) is the ith primitive lattice vector
            !! 
            !! Units : eV^(-1)
        real(dp) :: red_to_xyz(3, 3)
            !! Matrix converting reduced cooredinate positions to physical xyz cooredinates
            !! in eV^(-1) via
            !!
            !!      x_xyz = matmul(red_to_xyz, x_red)
            !!
            !! red_to_xyz = transpose(a_vecs)
            !! 
            !! Units : eV^(-1)
        real(dp) :: b_vecs(3, 3)
            !! Reciprocal lattice vectors
            !! 
            !! (i, :) is the ith reciprocal lattice vector
            !! 
            !! Units : eV
        real(dp) :: k_red_to_xyz(3, 3)
            !! Matrix converting reduced coordinate momentum to physical xyz coordinates
            !! in eV via
            !!
            !!      k_xyz = matmul(k_red_to_xyz, k_red)
            !!
            !! k_red_to_xyz = transpose(b_vecs)
            !! 
            !! Units : eV
        real(dp), allocatable :: energy_bands(:, :)
            !! Dim : [n_k, n_bands]
            !!
            !! Scissor corrected electron energy eigenvalues
            !!
            !! Units : eV 
        real(dp), allocatable :: G_grid_xyz(:, :)
            !! Dim : [n_in_G, 3]
            !!
            !! G vectors in physical xyz coordinates
            !!
            !! Units : eV
        real(dp), allocatable :: k_grid_xyz(:, :)
            !! Dim : [n_k, 3]
            !!
            !! List of k vectors in xyz coordinates
            !! 
            !! Units : eV
        logical :: include_spin = .FALSE.
            !! Flag declaring whether the wave functions are spin dependent or not

        real(dp) :: spin_degen = 2.0_dp
            !! Spin degeneracy factor.
            !!
            !! SI wave functions - 2
            !! SD wave functions - 1
        character(len=512) :: filename
            !! Filename of data
        contains

            procedure :: load => load_PW_dataset_hdf5
            procedure :: save => save_PW_dataset
            procedure :: print => print_PW_dataset
            procedure :: load_wfc_FT_ik_no_spin
            procedure :: load_wfc_ik_expanded_no_spin
            procedure :: load_wfc_FT_ik_spin
            procedure :: load_wfc_ik_expanded_spin
            procedure :: compute_PW_cutoff
            procedure :: do_scissor_correction

    end type

contains

    subroutine load_wfc_ik_expanded_no_spin(self, i, k, FFT_grid, wfc_ik)
        !! Loads the spin independent wave function coefficients, \( \widetilde{u} \) 
        !! for a given \( i, \mathbf{k} \), expands (zero-pads), then Fourier transforms them to get \( u \).

        use FFT_util

        implicit none

        class(PW_dataset_t) :: self
        integer :: i, k
        type(FFT_grid_t) :: FFT_grid

        complex(dp) :: wfc_ik(FFT_grid%n_grid(1), FFT_grid%n_grid(2), FFT_grid%n_grid(3))
        
        complex(dp) :: wfc_FT_ik(self%n_G)
        complex(dp) :: wfc_FT_ik_expanded(FFT_grid%n_grid(1), FFT_grid%n_grid(2), FFT_grid%n_grid(3))

        integer :: idx(3)
        integer :: g

        call self%load_wfc_FT_ik_no_spin(i, k, wfc_FT_ik)

        wfc_FT_ik_expanded = (0.0_dp, 0.0_dp)

        do g = 1, self%n_G
           call G_red_to_FFT_G_grid_index(FFT_grid%n_grid, self%G_grid_red(g, :), idx)
           wfc_FT_ik_expanded(idx(1), idx(2), idx(3)) = wfc_FT_ik(g)
        end do

        ! Fourier transform
        call dfftw_execute_dft(FFT_grid%plan, wfc_FT_ik_expanded, wfc_ik) 

    end subroutine

    subroutine load_wfc_ik_expanded_spin(self, i, k, FFT_grid, wfc_ik)
        !! Loads the spin dependent wave function coefficients, \( \widetilde{u} \) 
        !! for a given \( i, \mathbf{k} \), expands (zero-pads), then Fourier transforms them to get \( u \).

        use FFT_util

        implicit none

        class(PW_dataset_t) :: self
        integer :: i, k
        type(FFT_grid_t) :: FFT_grid

        complex(dp) :: wfc_ik(2, FFT_grid%n_grid(1), FFT_grid%n_grid(2), FFT_grid%n_grid(3))
        
        complex(dp) :: wfc_FT_ik(self%n_G, 2)
        complex(dp) :: wfc_FT_ik_expanded(2, FFT_grid%n_grid(1), FFT_grid%n_grid(2), FFT_grid%n_grid(3))

        integer :: idx(3)
        integer :: g, s

        call self%load_wfc_FT_ik_spin(i, k, wfc_FT_ik)

        wfc_FT_ik_expanded = (0.0_dp, 0.0_dp)

        do g = 1, self%n_G
           call G_red_to_FFT_G_grid_index(FFT_grid%n_grid, self%G_grid_red(g, :), idx)
           wfc_FT_ik_expanded(:, idx(1), idx(2), idx(3)) = wfc_FT_ik(g, :)
        end do

        ! Fourier transform
        do s = 1, 2
            call dfftw_execute_dft(FFT_grid%plan, wfc_FT_ik_expanded(s, :, :, :), wfc_ik(s, :, :, :)) 
        end do

    end subroutine

    subroutine load_wfc_FT_ik_no_spin(self, i, k, wfc_FT_ik)
        !! Loads the spin independent Bloch coefficients, \( \widetilde{u}_{ i, \mathbf{k}, \mathbf{G}} \)
        !! for a given \( i, \mathbf{k} \).
        implicit none

        class(PW_dataset_t) :: self
        integer :: i, k

        complex(dp) :: wfc_FT_ik(:)
            !! Dim : [n_G]

        real(dp) :: wfc_FT_ik_r(self%n_G)
            !! real part of the bloch wave functions in fourier space
        real(dp) :: wfc_FT_ik_c(self%n_G)
            !! complex part of the bloch wave functions in fourier space

        integer(HSIZE_T) :: dims(1)

        integer(HID_T) :: file_id
            !! HDF5 file ID number for DFT input file

        integer :: error

        character(len=64) :: dataset_path_r
        character(len=64) :: dataset_path_c

        dataset_path_r = 'wfc_FT_r'//&
                            '/i_'//trim(adjustl(int_to_str(i)))//&
                            '/k_'//trim(adjustl(int_to_str(k)))
        dataset_path_c = 'wfc_FT_c'//&
                            '/i_'//trim(adjustl(int_to_str(i)))//&
                            '/k_'//trim(adjustl(int_to_str(k)))

        dims = [self%n_G]

        call h5open_f(error)
        call h5fopen_f(self%filename, H5F_ACC_RDONLY_F, file_id, error)

        call h5ltread_dataset_double_f(file_id, dataset_path_r, wfc_FT_ik_r, dims, error)
        call h5ltread_dataset_double_f(file_id, dataset_path_c, wfc_FT_ik_c, dims, error)

        call h5fclose_f(file_id, error)
        call h5close_f(error)

        wfc_FT_ik = wfc_FT_ik_r + ii*wfc_FT_ik_c

    end subroutine

    subroutine load_wfc_FT_ik_spin(self, i, k, wfc_FT_ik)
        !! Loads the spin dependent Bloch coefficients, \( \widetilde{u}_{ i, \mathbf{k}, \mathbf{G}, s} \)
        !! for a given \( i, \mathbf{k} \).
        implicit none

        class(PW_dataset_t) :: self
        integer :: i, k

        complex(dp) :: wfc_FT_ik(:, :)
            !! Dim : [n_G, 2]

        real(dp) :: wfc_FT_ik_r(self%n_G, 2)
            !! real part of the bloch wave functions in fourier space
        real(dp) :: wfc_FT_ik_c(self%n_G, 2)
            !! complex part of the bloch wave functions in fourier space

        integer(HSIZE_T) :: dims(2)

        integer(HID_T) :: file_id
            !! HDF5 file ID number for DFT input file

        integer :: error

        character(len=64) :: dataset_path_r
        character(len=64) :: dataset_path_c

        dataset_path_r = 'wfc_FT_r'//&
                            '/i_'//trim(adjustl(int_to_str(i)))//&
                            '/k_'//trim(adjustl(int_to_str(k)))
        dataset_path_c = 'wfc_FT_c'//&
                            '/i_'//trim(adjustl(int_to_str(i)))//&
                            '/k_'//trim(adjustl(int_to_str(k)))

        dims = [self%n_G, 2]

        call h5open_f(error)
        call h5fopen_f(self%filename, H5F_ACC_RDONLY_F, file_id, error)

        call h5ltread_dataset_double_f(file_id, dataset_path_r, wfc_FT_ik_r, dims, error)
        call h5ltread_dataset_double_f(file_id, dataset_path_c, wfc_FT_ik_c, dims, error)

        call h5fclose_f(file_id, error)
        call h5close_f(error)

        wfc_FT_ik = wfc_FT_ik_r + ii*wfc_FT_ik_c

    end subroutine

    subroutine load_PW_dataset_hdf5(self, filename, verbose)
        !! Loads `PW_dataset` parameters from an hdf5 file.
        implicit none

        class(PW_dataset_t) :: self
        character(len=*) :: filename
        logical, optional :: verbose

        integer(HID_T) :: file_id
        logical :: file_exists
        integer(HSIZE_T) :: dims(1) = [0]
        integer(HSIZE_T) :: dims2(2)
        integer :: error

        integer :: g, k
        integer :: wfc_data_rank

        if ( verbose ) then
            print*, 'Loading PW dataset input file...'
            print*
        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)

            self%filename = filename

            ! read the data
            call h5ltread_dataset_int_f(file_id, 'n_k', self%n_k, dims, error)
            call h5ltread_dataset_int_f(file_id, 'n_val', self%n_val, dims, error)
            call h5ltread_dataset_int_f(file_id, 'n_cond', self%n_cond, dims, error)
            self%n_bands = self%n_val + self%n_cond

            call h5ltread_dataset_double_f(file_id, 'Ef_max', self%Ef_max, dims, error)

            call h5ltread_dataset_int_f(file_id, 'n_G', self%n_G, dims, error)

            allocate(self%G_grid_red(self%n_G, 3))        
            dims2 = [self%n_G, 3]
            call h5ltread_dataset_int_f(file_id, 'G_grid_red', self%G_grid_red, dims2, error)

            dims2 = [3, 3]
            call h5ltread_dataset_double_f(file_id, 'a_vecs_A', self%a_vecs_A, dims2, error)
            self%a_vecs = Ang_to_inv_eV*self%a_vecs_A
            self%red_to_xyz = transpose(self%a_vecs)

            call h5ltread_dataset_double_f(file_id, 'b_vecs_A', self%b_vecs_A, dims2, error)
            self%b_vecs = inv_Ang_to_eV*self%b_vecs_A
            self%k_red_to_xyz = transpose(self%b_vecs)

            ! convert reduced coordinates to xyz
            allocate(self%G_grid_xyz(self%n_G, 3))
            do g = 1, self%n_G
                self%G_grid_xyz(g, :) = matmul(self%k_red_to_xyz, self%G_grid_red(g, :))
            end do

            allocate(self%k_weight(self%n_k))
            dims = [self%n_k]
            call h5ltread_dataset_double_f(file_id, 'k_weight', self%k_weight, dims, error) 

            allocate(self%k_grid_red(self%n_k, 3))
            dims2 = [self%n_k, 3]
            call h5ltread_dataset_double_f(file_id, 'k_grid_red', self%k_grid_red, dims2, error) 

            allocate(self%k_grid_xyz(self%n_k, 3))
            do k = 1, self%n_k
                self%k_grid_xyz(k, :) = matmul(self%k_red_to_xyz, self%k_grid_red(k, :))
            end do

            allocate(self%energy_bands_raw(self%n_k, self%n_bands))
            allocate(self%energy_bands(self%n_k, self%n_bands))
            dims2 = [self%n_k, self%n_bands]
            call h5ltread_dataset_double_f(file_id, 'energy_bands', self%energy_bands_raw, dims2, error) 

            call h5ltget_dataset_ndims_f(file_id, 'wfc_FT_c/i_1/k_1', wfc_data_rank, error)

            if ( wfc_data_rank == 2 ) then
                self%include_spin = .TRUE.
                self%spin_degen = 1.0_dp
            end if

            call h5fclose_f(file_id, error)
            call h5close_f(error)

            ! call check_DFT_parameters(verbose)

            call self%compute_PW_cutoff()

            call self%print(filename, verbose = verbose)

        else

            call print_error_message('Input file for PW dataset : '//trim(filename)//' does NOT exist.', &
                verbose = verbose)
            stop

        end if

    end subroutine

    subroutine save_PW_dataset(self, filename, verbose)
        !! Saves `PW_dataset`.
        implicit none

        class(PW_dataset_t) :: self
        character(len=*) :: filename
        logical, optional :: verbose

        integer(HID_T) :: file_id, group_id
        logical :: file_exists
        integer(HSIZE_T) :: dims1(1) = [1]
        integer(HSIZE_T) :: dims2(2) = [3, 3]
        integer :: error

        if ( verbose ) then
            print*, 'Saving PW dataset parameters...'
            print*
        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'PW_dataset', group_id, error)

            ! write data
            call h5ltmake_dataset_int_f(file_id, &
                'PW_dataset/n_val', &
                size(dims1), dims1,&
                self%n_val, &
                error)
            call h5ltmake_dataset_int_f(file_id, &
                'PW_dataset/n_cond', &
                size(dims1), dims1,&
                self%n_cond, &
                error)
            call h5ltmake_dataset_int_f(file_id, &
                'PW_dataset/n_bands', &
                size(dims1), dims1,&
                self%n_bands, &
                error)
            call h5ltmake_dataset_int_f(file_id, &
                'PW_dataset/n_k', &
                size(dims1), dims1,&
                self%n_k, &
                error)
            call h5ltmake_dataset_int_f(file_id, &
                'PW_dataset/n_G', &
                size(dims1), dims1,&
                self%n_G, &
                error)
            call h5ltmake_dataset_double_f(file_id, &
                'PW_dataset/E_cut', &
                size(dims1), dims1,&
                self%E_cut, &
                error)
            call h5ltmake_dataset_double_f(file_id, &
                'PW_dataset/Ef_max', &
                size(dims1), dims1,&
                self%Ef_max, &
                error)
            call h5ltmake_dataset_double_f(file_id, &
                'PW_dataset/q_cut', &
                size(dims1), dims1,&
                self%q_cut, &
                error)
            call h5ltmake_dataset_double_f(file_id, &
                'PW_dataset/spin_degen', &
                size(dims1), dims1,&
                self%spin_degen, &
                error)
            call h5ltmake_dataset_double_f(file_id, &
                'PW_dataset/k_red_to_xyz', &
                size(dims2), dims2,&
                self%k_red_to_xyz, &
                error)
            call h5ltmake_dataset_double_f(file_id, &
                'PW_dataset/red_to_xyz', &
                size(dims2), dims2,&
                self%red_to_xyz, &
                error)

            call h5fclose_f(file_id, error)
            call h5close_f(error)

        else

            call print_error_message('Output file : '//trim(filename)//' does NOT exist.', &
                verbose = verbose)
            stop

        end if

    end subroutine

    ! subroutine check_DFT_parameters(verbose)
    !     !! Check the values of of the input file to make sure they make sense

    !     implicit none

    !     logical, optional :: verbose
        
    !     real(dp), parameter :: pi = 4.0_dp*datan(1.0_dp)

    !     real(dp) :: i_3_2pi(3, 3) 

    !     real(dp) :: eps_val

    !     if ( verbose ) then
    !         print*, 'Running preliminary checks on DFT data...'
    !         print*
    !     end if

    !     i_3_2pi = 0.0_dp

    !     i_3_2pi(1, 1) = 2*pi
    !     i_3_2pi(2, 2) = 2*pi
    !     i_3_2pi(3, 3) = 2*pi

    !     ! a_i . b_j = 2*pi delta_ij
    !     eps_val = abs(sum(matmul(transpose(k_red_to_xyz), red_to_xyz) - i_3_2pi)/3.0_dp)

    !     if ( eps_val > 1e-3_dp ) then

    !         print*, '!! ERROR !!'
    !         print*
    !         print*, '   Basis vectors are not orthonormalized correctly. a_i . b_j != 2 pi delta_ij'
    !         print*
    !         print*, '!!!!!!!!!!!'
    !         print*

    !         stop

    !     end if

    !     ! sum_k k_weight(k) = 2

    !     eps_val = sum(k_weight) - 2.0_dp

    !     if ( eps_val > 1e-3_dp ) then

    !         print*, '!! ERROR !!'
    !         print*
    !         print*, '   Sum of k weights != 2'
    !         print*
    !         print*, '!!!!!!!!!!!'
    !         print*

    !         stop

    !     end if

    ! end subroutine

    subroutine print_PW_dataset(self, filename, verbose)
        !! Prints `PW_dataset` components.
        implicit none

        class(PW_dataset_t) :: self
        character(len=*) :: filename
        logical, optional :: verbose

        if ( verbose ) then
            call print_section_seperator()
            print*, '    ----------'
            print*, '    PW dataset'
            print*, '    ----------'
            print*
            print*, '        Filename : ', trim(filename)
            print*
            print*, '        Primitive lattice vectors (Ang) : '
            print*, '            a1 = ', self%a_vecs_A(1, :)
            print*, '            a2 = ', self%a_vecs_A(2, :)
            print*, '            a3 = ', self%a_vecs_A(3, :)
            print*
            print*, '        Reciprocal lattice vectors Ang^(-1) : '
            print*, '            b1 = ', self%b_vecs_A(1, :) 
            print*, '            b2 = ', self%b_vecs_A(2, :)
            print*, '            b3 = ', self%b_vecs_A(3, :)
            print*
            print*, '        Number of valence bands     = ', trim(adjustl(int_to_str(self%n_val)))
            print*, '        Number of conduction bands  = ', trim(adjustl(int_to_str(self%n_cond)))
            print*
            print*, '        Number of k points          = ', trim(adjustl(int_to_str(self%n_k)))
            print*, '        Number of G points          = ', trim(adjustl(int_to_str(self%n_G)))
            print*
            print*, '        E_cut  = ', self%E_cut, 'eV'
            print*, '        q_cut  = ', self%q_cut/1.0e3_dp, 'keV'
            print*
            print*, '        Ef_max = ', self%Ef_max, 'eV'
            print*
            print*, '        Include spin?               = ', self%include_spin
            print*
            call print_section_seperator()
            print*
        end if

    end subroutine

    subroutine compute_PW_cutoff(self)
        !! Computes plane wave expansion parameters, \( q_\text{cut} \) and \( E_\text{cut} \).
        implicit none

        class(PW_dataset_t) :: self

        integer :: k, g

        real(dp) :: q_xyz(3)
        real(dp) :: q_mag

        self%q_cut = 0.0_dp
        self%E_cut = 0.0_dp

        do k = 1, self%n_k
            do g = 1, self%n_G

                q_xyz = self%k_grid_xyz(k, :) + self%G_grid_xyz(g, :)

                q_mag = norm2(q_xyz)

                if ( q_mag >= self%q_cut ) then
                    self%q_cut = q_mag
                end if

            end do
        end do

        self%E_cut = self%q_cut**2/(2.0_dp*m_elec)

    end subroutine

    subroutine do_scissor_correction(self, band_gap, verbose)
        !! Performs a scissor correction to the band structure. This shifts the valence and conduction
        !! bands independently such that \( \text{min} E_{\text{cond}, \mathbf{k}} - \text{max} E_{\text{val}, \mathbf{k}} = \)
        !!` band_gap`.
        !!
        !! Additionally, shifts all of the bands by a constant value such that the valence band maximum, \( \text{max}
        !! E_{\text{val}, \mathbf{k}} \), is equal to 0.
        implicit none

        class(PW_dataset_t) :: self
        real(dp) :: band_gap
        logical, optional :: verbose

        real(dp) :: scissor_correct

        ! scissor correct band structure
        scissor_correct = band_gap - &
                    (minval(self%energy_bands_raw(:, self%n_val + 1:)) - &
                     maxval(self%energy_bands_raw(:, :self%n_val)))

        self%energy_bands = 0.0_dp
        self%energy_bands(:, :self%n_val) = self%energy_bands_raw(:, :self%n_val) - scissor_correct/2.0_dp 
        self%energy_bands(:, self%n_val + 1:) = self%energy_bands_raw(:, self%n_val + 1:) + scissor_correct/2.0_dp 

        ! shift valence band maximum to 0 by shifting the whole band structure 
        self%energy_bands = self%energy_bands - maxval(self%energy_bands(:, self%n_val))

        if ( verbose ) then

            print*, 'Performed scissor correction.'
            print*

        end if
    
    end subroutine

end module
