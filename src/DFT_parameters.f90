module DFT_parameters
    !! Handles the results from DFT calculations that compute the Bloch wave 
    !! function coefficients.
    !!
    !! Variables with '_A' are with units of Angstroms
   
    use h5lt
    use hdf5

    use prec
    use constants
    use units

    implicit none

    integer :: n_k
        !! Number of k points
    integer :: n_val
        !! Number of valence bands
    integer :: n_cond
        !! Number of conduction bands
    integer :: n_in_G
        !! Number of G points in in_G_grid
    integer, allocatable :: in_G_grid_red(:, :)
        !! Dim : [n_in_G, 3]
        !!
        !! Input grid of G vectors in reduced coordinates
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
        !! Raw electron energy eigenvalues
        !!
        !! Units : eV 

    !! Generated from input

    integer :: n_bands
        !! Number of bands
        !! 
        !! n_bands = n_cond + n_val
    real(dp) :: q_PW_cut
        !! plane wave cutoff in expansion
        !!
        !! Units : eV
    real(dp) :: E_PW_cut
        !! energy cooresponding to q_PW_cut
        !!
        !! E_PW_cut = q_PW_cut**2/(2*m_elec)
        !!
        !! Units : eV
    real(dp) :: a_vecs(3, 3)
        !! Primitive lattice vectors
        !!
        !! (i, :) is the ith primitive lattice vector
        !! 
        !! Units : eV
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
        !! Reciprocal lattice vector
        !! 
        !! (i, :) is the ith reciprocal lattice vector
        !! 
        !! Units : eV
    real(dp) :: k_red_to_xyz(3, 3)
        !! Matrix converting reduced cooredinate momentum to physical xyz cooredinates
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

    real(dp), allocatable :: in_G_grid_xyz(:, :)
        !! Dim : [n_in_G, 3]
        !!
        !! Input grid of G vectors in physical xyz coordinates
        !!
        !! Units : eV
    real(dp), allocatable :: k_grid_xyz(:, :)
        !! Dim : [n_k, 3]
        !!
        !! List of k vectors in xyz coordinates
        !! 
        !! Units : eV

    !!! experimental

    logical :: include_spin = .FALSE.

    interface get_in_wfc_FT
        module procedure get_in_wfc_FT_no_spin
        module procedure get_in_wfc_FT_spin
    end interface

    contains

    subroutine save_DFT_parameters(filename, verbose)
        !! Saves some of the DFT input variables to filename
        implicit none

        character(len=*) :: filename

        integer(HID_T) :: file_id, group_id

        logical, optional :: verbose
        logical :: file_exists

        integer(HSIZE_T) :: dims1(1) = [1]
        integer(HSIZE_T) :: dims2(2) = [3, 3]

        integer :: error

        if ( verbose ) then

            print*, 'Saving DFT parameters...'
            print*

        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'DFT_parameters', group_id, error)

            ! write data
            call h5ltmake_dataset_int_f(file_id, 'DFT_parameters/n_val', size(dims1), dims1,&
                n_val, error)
            call h5ltmake_dataset_int_f(file_id, 'DFT_parameters/n_cond', size(dims1), dims1,&
                n_cond, error)
            call h5ltmake_dataset_int_f(file_id, 'DFT_parameters/n_bands', size(dims1), dims1,&
                n_bands, error)
            call h5ltmake_dataset_int_f(file_id, 'DFT_parameters/n_k', size(dims1), dims1,&
                n_k, error)
            call h5ltmake_dataset_int_f(file_id, 'DFT_parameters/n_in_G', size(dims1), dims1,&
                n_in_G, error)
            call h5ltmake_dataset_double_f(file_id, 'DFT_parameters/E_PW_cut', size(dims1), dims1,&
                E_PW_cut, error)
            call h5ltmake_dataset_double_f(file_id, 'DFT_parameters/q_PW_cut', size(dims1), dims1,&
                q_PW_cut, error)
            call h5ltmake_dataset_double_f(file_id, 'DFT_parameters/k_red_to_xyz', size(dims2), dims2,&
                k_red_to_xyz, error)
            call h5ltmake_dataset_double_f(file_id, 'DFT_parameters/red_to_xyz', size(dims2), dims2,&
                red_to_xyz, error)

            call h5fclose_f(file_id, error)
            call h5close_f(error)

        else

            if ( verbose ) then

                print*, '!! ERROR !!'
                print*
                print*, '   Output file : ', trim(filename), ' does NOT exist.'
                print*
                print*, '!!!!!!!!!!!'
                print*

            end if

            stop

        end if

    end subroutine

    subroutine load_DFT_parameters(filename, verbose)
        !! Loads the DFT input file
        implicit none

        character(len=*) :: filename

        integer(HID_T) :: file_id

        logical, optional :: verbose
        logical :: file_exists

        integer(HSIZE_T) :: dims(1) = [0]
        integer(HSIZE_T) :: dims2(2)

        integer :: error

        integer :: g, k

        integer :: wfc_data_rank

        if ( verbose ) then

            print*, 'Loading DFT input file...'
            print*

        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)

            ! read the data
            call h5ltread_dataset_int_f(file_id, 'n_k', n_k, dims, error)
            call h5ltread_dataset_int_f(file_id, 'n_val', n_val, dims, error)
            call h5ltread_dataset_int_f(file_id, 'n_cond', n_cond, dims, error)
            n_bands = n_val + n_cond

            call h5ltread_dataset_int_f(file_id, 'n_in_G', n_in_G, dims, error)

            allocate(in_G_grid_red(n_in_G, 3))        
            dims2 = [n_in_G, 3]
            call h5ltread_dataset_int_f(file_id, 'in_G_grid_red', in_G_grid_red, dims2, error)

            dims2 = [3, 3]
            call h5ltread_dataset_double_f(file_id, 'a_vecs_A', a_vecs_A, dims2, error)
            a_vecs = Ang_to_inv_eV*a_vecs_A
            red_to_xyz = transpose(a_vecs)

            call h5ltread_dataset_double_f(file_id, 'b_vecs_A', b_vecs_A, dims2, error)
            b_vecs = inv_Ang_to_eV*b_vecs_A
            k_red_to_xyz = transpose(b_vecs)

            ! convert reduced coordinates to xyz
            allocate(in_G_grid_xyz(n_in_G, 3))
            do g = 1, n_in_G
                in_G_grid_xyz(g, :) = matmul(k_red_to_xyz, in_G_grid_red(g, :))
            end do

            allocate(k_weight(n_k))
            dims = [n_k]
            call h5ltread_dataset_double_f(file_id, 'k_weight', k_weight, dims, error) 

            allocate(k_grid_red(n_k, 3))
            dims2 = [n_k, 3]
            call h5ltread_dataset_double_f(file_id, 'k_red', k_grid_red, dims2, error) 

            allocate(k_grid_xyz(n_k, 3))
            do k = 1, n_k
                k_grid_xyz(k, :) = matmul(k_red_to_xyz, k_grid_red(k, :))
            end do

            allocate(energy_bands_raw(n_k, n_bands))
            allocate(energy_bands(n_k, n_bands))
            dims2 = [n_k, n_bands]
            call h5ltread_dataset_double_f(file_id, 'energy_bands_raw', energy_bands_raw, dims2, error) 

            call h5ltget_dataset_ndims_f(file_id, 'in_wfc_FT_c/1', wfc_data_rank, error)

            if ( wfc_data_rank == 3 ) then
                include_spin = .TRUE.
            end if

            call h5fclose_f(file_id, error)
            call h5close_f(error)

            call check_DFT_parameters(verbose)

            call print_DFT_parameters(filename, verbose=verbose)

            ! find how much the wave function coefficients were expanded
            call get_PW_cutoffs(verbose=verbose)

            if ( verbose ) then
                print*, '----------------------------------------'
                print*
            end if

        else

            if ( verbose ) then

                print*, '!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*
                print*, '    Input file for DFT parameters : ', trim(filename), ' does NOT exist.'
                print*
                print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*

            end if

            stop

        end if

    end subroutine

    subroutine check_DFT_parameters(verbose)
        !! Check the values of of the input file to make sure they make sense

        implicit none

        logical, optional :: verbose
        
        real(dp), parameter :: pi = 4.0_dp*datan(1.0_dp)

        real(dp) :: i_3_2pi(3, 3) 

        real(dp) :: eps_val

        if ( verbose ) then
            print*, 'Running preliminary checks on DFT data...'
            print*
        end if

        i_3_2pi = 0.0_dp

        i_3_2pi(1, 1) = 2*pi
        i_3_2pi(2, 2) = 2*pi
        i_3_2pi(3, 3) = 2*pi

        ! a_i . b_j = 2*pi delta_ij
        eps_val = abs(sum(matmul(transpose(k_red_to_xyz), red_to_xyz) - i_3_2pi)/3.0_dp)

        if ( eps_val > 1e-3_dp ) then

            print*, '!! ERROR !!'
            print*
            print*, '   Basis vectors are not orthonormalized correctly. a_i . b_j != 2 pi delta_ij'
            print*
            print*, '!!!!!!!!!!!'
            print*

            stop

        end if

        ! sum_k k_weight(k) = 2

        eps_val = sum(k_weight) - 2.0_dp

        if ( eps_val > 1e-3_dp ) then

            print*, '!! ERROR !!'
            print*
            print*, '   Sum of k weights != 2'
            print*
            print*, '!!!!!!!!!!!'
            print*

            stop

        end if

    end subroutine

    subroutine print_DFT_parameters(filename, verbose)
        !! Prints variables defined in this moudle
        implicit none

        character(len=*) :: filename

        logical, optional :: verbose

        character(len=64) :: n_k_str
        character(len=64) :: n_in_G_str
        character(len=64) :: n_val_str
        character(len=64) :: n_cond_str

        write(n_k_str, *) n_k
        write(n_val_str, *) n_val
        write(n_cond_str, *) n_cond
        write(n_in_G_str, *) n_in_G

        if ( verbose ) then

            print*, '----------------------------------------'
            print*, '    ---------'
            print*, '    DFT Input'
            print*, '    ---------'
            print*
            print*, '        Filename : ', trim(filename)
            print*
            print*, '        Primitive lattice vectors (Ang) : '
            print*, '            a1 = ', a_vecs_A(1, :)
            print*, '            a2 = ', a_vecs_A(2, :)
            print*, '            a3 = ', a_vecs_A(3, :)
            print*
            print*, '        Reciprocal lattice vectors Ang^(-1) : '
            print*, '            b1 = ', b_vecs_A(1, :) 
            print*, '            b2 = ', b_vecs_A(2, :)
            print*, '            b3 = ', b_vecs_A(3, :)
            print*
            print*, '        Number of valence bands     = ', trim(adjustl(n_val_str))
            print*, '        Number of conduction bands  = ', trim(adjustl(n_cond_str))
            print*
            print*, '        Number of k points          = ', trim(adjustl(n_k_str))
            print*
            print*, '        Number of G points          = ', trim(adjustl(n_in_G_str))
            print*
            print*, '        Include spin?               = ', include_spin
            print*

        end if

    end subroutine

    subroutine get_in_wfc_FT_no_spin(filename, band_num, in_wfc_FT)
        !! Reads the input DFT file and sets the [n_k, n_in_G] complex(dp) array 
        !! with whose elements are the dimensionless block wave function 
        !! coefficients, u_ikG
        implicit none

        character(len=*) :: filename
        integer :: band_num

        complex(dp) :: in_wfc_FT(:, :)
            !! Bloch wave functions in reciprocal space

        integer(HID_T) :: file_id
            !! HDF5 file ID number for DFT input file

        real(dp) :: in_wfc_FT_r(n_k, n_in_G)
            !! real part of the bloch wave functions in fourier space
        real(dp) :: in_wfc_FT_c(n_k, n_in_G)
            !! complex part of the bloch wave functions in fourier space

        integer(HSIZE_T) :: dims2(2)

        integer :: error

        character(len=64) :: dataset_path_r
        character(len=64) :: dataset_path_c
        character(len=64) :: band_num_str

        write(band_num_str, *) band_num

        dataset_path_r = 'in_wfc_FT_r/'//trim(adjustl(band_num_str))
        dataset_path_c = 'in_wfc_FT_c/'//trim(adjustl(band_num_str))

        dims2 = [n_k, n_in_G]

        call h5open_f(error)
        call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)

        call h5ltread_dataset_double_f(file_id, dataset_path_r, in_wfc_FT_r, dims2, error)
        call h5ltread_dataset_double_f(file_id, dataset_path_c, in_wfc_FT_c, dims2, error)

        call h5fclose_f(file_id, error)
        call h5close_f(error)

        in_wfc_FT = in_wfc_FT_r + ii*in_wfc_FT_c

    end subroutine

    subroutine get_in_wfc_FT_spin(filename, band_num, in_wfc_FT)
        !! Reads the input DFT file and sets the [n_k, n_in_G, 2] complex(dp) array 
        !! with whose elements are the dimensionless block wave function 
        !! coefficients, u_{i, k, G, s}, where s is the spin index
        implicit none

        character(len=*) :: filename
        integer :: band_num

        complex(dp) :: in_wfc_FT(:, :, :)
            !! Bloch wave functions in reciprocal space

        integer(HID_T) :: file_id
            !! HDF5 file ID number for DFT input file

        real(dp) :: in_wfc_FT_r(n_k, n_in_G, 2)
            !! real part of the bloch wave functions in fourier space
        real(dp) :: in_wfc_FT_c(n_k, n_in_G, 2)
            !! complex part of the bloch wave functions in fourier space

        integer(HSIZE_T) :: dims3(3)

        integer :: error

        character(len=64) :: dataset_path_r
        character(len=64) :: dataset_path_c
        character(len=64) :: band_num_str

        write(band_num_str, *) band_num

        dataset_path_r = 'in_wfc_FT_r/'//trim(adjustl(band_num_str))
        dataset_path_c = 'in_wfc_FT_c/'//trim(adjustl(band_num_str))

        dims3 = [n_k, n_in_G, 2]

        call h5open_f(error)
        call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)

        call h5ltread_dataset_double_f(file_id, dataset_path_r, in_wfc_FT_r, dims3, error)
        call h5ltread_dataset_double_f(file_id, dataset_path_c, in_wfc_FT_c, dims3, error)

        call h5fclose_f(file_id, error)
        call h5close_f(error)

        in_wfc_FT = in_wfc_FT_r + ii*in_wfc_FT_c

    end subroutine

    subroutine expand_wfc_FT_for_FFT(n_grid,&
            wfc_FT, wfc_FT_exp,&
            verbose)
        !! The Fourier components will generally not be defined on a
        !! uniform grid needed for an FFT. This subroutine puts the 
        !! values of the Bloch coefficients in the correct place 
        use FFT_util

        implicit none

        integer :: n_grid(3)

        complex(dp) :: wfc_FT(n_in_G)
        complex(dp) :: wfc_FT_exp(n_grid(1), n_grid(2), n_grid(3))

        logical, optional :: verbose

        integer :: g
        integer :: FFT_idx(3)

        wfc_FT_exp = (0.0_dp, 0.0_dp)

        do g = 1, n_in_G

            call G_red_to_FFT_G_grid_index(n_grid, in_G_grid_red(g, :), FFT_idx)

            wfc_FT_exp(FFT_idx(1), FFT_idx(2), FFT_idx(3)) = wfc_FT(g)

        end do

    end subroutine

    function in_wfc_at_single_pt(band_num, x0_red, in_wfc_FT) result(wfc)
        !! returns the wave function evaluated at a point in
        !! reduced coordinates at a single point
        !!
        !! av_wfc(x0)_k = sum_G e^(i G.x0) wfc(G)_k
        implicit none

        integer :: band_num

        real(dp) :: x0_red(3)
        real(dp) :: x0(3)
        complex(dp) :: phase_fac(n_in_G)
        complex(dp) :: wfc(n_k)

        complex(dp) :: in_wfc_FT(n_k, n_in_G)

        integer :: g

        x0 = matmul(red_to_xyz, x0_red)

        phase_fac = cmplx(0.0_dp, 0.0_dp, dp)

        do g = 1, n_in_G
            phase_fac(g) = exp(ii*dot_product(in_G_grid_xyz(g, :), x0))
        end do

        wfc = matmul(in_wfc_FT, phase_fac)

    end function

    subroutine get_PW_cutoffs(verbose)
        !! Calculate the plane wave expansion parameters

        implicit none

        logical, optional :: verbose

        integer :: k, g

        real(dp) :: q_xyz(3)
        real(dp) :: q_mag

        q_PW_cut = 0.0_dp
        E_PW_cut = 0.0_dp

        do k = 1, n_k
            do g = 1, n_in_G

                q_xyz = k_grid_xyz(k, :) + in_G_grid_xyz(g, :)

                q_mag = norm2(q_xyz)

                if ( q_mag .ge. q_PW_cut ) then
                    q_PW_cut = q_mag
                end if

            end do
        end do

        E_PW_cut = q_PW_cut**2/(2.0_dp*m_elec)

        if ( verbose ) then

           print*, '        Plane wave expansion parameters : '
           print*, '            E_PW_cut = ', E_PW_cut, 'eV'
           print*, '            q_PW_cut = ', q_PW_cut/1.0e3_dp, 'keV'
           print*

        end if

    end subroutine

    subroutine do_scissor_correction(band_gap, verbose)
        implicit none

        logical, optional :: verbose

        real(dp) :: band_gap

        real(dp) :: scissor_correct

        scissor_correct = band_gap - &
                    (minval(energy_bands_raw(:, n_val + 1:)) - maxval(energy_bands_raw(:, :n_val)))

        energy_bands = 0.0_dp
        energy_bands(:, :n_val) = energy_bands_raw(:, :n_val) - scissor_correct/2.0_dp 
        energy_bands(:, n_val + 1:) = energy_bands_raw(:, n_val + 1:) + scissor_correct/2.0_dp 

        if ( verbose ) then

            print*, 'Performed scissor correction.'
            print*

        end if
    
    end subroutine

    function get_w_max(i) result(w_max)
        !! Returns the maximum energy transfer allowed from valence state i
        implicit none

        integer :: i

        real(dp) :: w_max

        w_max = maxval(energy_bands) - minval(energy_bands(:, i))

    end function

end module
