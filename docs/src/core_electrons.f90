module core_electrons
    !! Handles the wave functions and energy levels for core electrons

    use hdf5
    use h5lt

    use prec
    use constants
    use math_mod

    implicit none

    integer :: n_atoms
        !! Number of atoms in the primitive cell

    integer :: n_core_states
        !! Total number of core states

    integer, allocatable :: Z_list(:)
        !! Dim : [n_atoms]
        !! 
        !! Proton number / number of electrons

    real(dp), allocatable :: eq_pos_red(:, :)
        !! Dim : [n_atoms, 3] 
        !! 
        !! Equilibrium position of the atoms in the primitive cell
        !! in reduced coordinates
        !!
        !! Units : None

    real(dp), allocatable :: core_energy(:)
        !! Dim : [n_core_states]
        !!
        !! Energy of the core states
        !! 
        !! Units : eV

    integer, allocatable :: core_elec_conf(:, :)
        !! Dim : [n_core_states, 5]
        !!
        !! Electron configuration, each element is
        !!  core_elec_conf(i, 1) - atom id
        !!  core_elec_conf(i, 2) - n
        !!  core_elec_conf(i, 3) - l
        !!  core_elec_conf(i, 4) - m
        !!  core_elec_conf(i, 5) - n_s (number of spin states)

    !! Generated STO parameters

    integer, allocatable :: core_sto_nj_list(:)
        !! Dim : [n_core_states]
        !!
        !! Number of terms in the analytic STO expansion

    integer :: core_sto_max_nj
        !! Maximum value of core_sto_nj_list

    real(dp), allocatable :: core_sto_data(:, :, :)
        !! Dim : [n_core_states, 4, max(core_sto_nj_list)] 
        !!
        !! All of the parameters of the STO expansion
        !!  core_sto_data(i, 1, :) - n_l[j]
        !!  core_sto_data(i, 2, :) - Z_l[j]
        !!  core_sto_data(i, 3, :) - N_l[j]
        !!  core_sto_data(i, 4, :) - C_nl[j]

contains

    subroutine save_core_electrons(filename, verbose)
        implicit none

        character(len=*) :: filename

        integer(HID_T) :: file_id, group_id

        logical, optional :: verbose
        logical :: file_exists

        integer(HSIZE_T) :: dims1(1) = [1]
        integer(HSIZE_T) :: dims2(2) 

        integer :: error

        if ( verbose ) then

            print*, 'Saving core electron configuration...'
            print*

        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'core', group_id, error)

            ! write data
            call h5ltmake_dataset_int_f(file_id, 'core/n_core_states', size(dims1), dims1,&
                n_core_states, error)
            call h5ltmake_dataset_int_f(file_id, 'core/n_atoms', size(dims1), dims1,&
                n_atoms, error)
            dims1 = [n_atoms]
            call h5ltmake_dataset_int_f(file_id, 'core/Z_list', size(dims1), dims1,&
                Z_list, error)
            dims2 = [n_atoms, 3]
            call h5ltmake_dataset_double_f(file_id, 'core/eq_pos_red', size(dims2), dims2,&
                eq_pos_red, error)
            dims1 = [n_core_states]
            call h5ltmake_dataset_double_f(file_id, 'core/core_energy', size(dims1), dims1,&
                core_energy, error)
            dims2 = [n_core_states, 5]
            call h5ltmake_dataset_int_f(file_id, 'core/core_elec_config', size(dims2), dims2,&
                core_elec_conf, error)

            call h5fclose_f(file_id, error)
            call h5close_f(error)

        else

            if ( verbose ) then

                print*, '!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*
                print*, '    Output file : ', trim(filename), ' does NOT exist.'
                print*
                print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*

            end if

            stop

        end if

    end subroutine

    subroutine calc_core_sto_wf_grid(core_wfc, core_state_idx, n_grid, red_to_xyz,&
            shift, k_vec_in, s_cut_in, verbose) 
        !! computes the core wave functions on a grid, summed over the closest unit cells
        !!
        !! if shift == T : x_shift = x_pos_red (equilibrium position)
        !! else: x_shift = 0
        !! core_sto_wf_grid(x) = sum_s core_sto_wf(x - x_shift + s)
        !!
        !! where x_i is the atomic equilibrium position
        implicit none

        logical, optional :: verbose
        logical, optional :: shift

        integer :: n_grid(3)
            !! Number of grid points in reduced x coordinates

        complex(dp) :: core_wfc(n_grid(1), n_grid(2), n_grid(3))

        real(dp), optional :: k_vec_in(3)
        real(dp) :: k_vec(3)

        integer, optional :: s_cut_in
        integer :: s_cut

        integer :: core_state_idx

        integer :: n1, n2, n3
        integer :: s1, s2, s3

        real(dp) :: x_red(3)
        real(dp) :: s_red(3)
        real(dp) :: x_vec(3)
        real(dp) :: x_shift_red(3)
        real(dp) :: x_shift(3)

        real(dp) :: red_to_xyz(3, 3)

        complex(dp) :: core_sto_wf_sum

        complex(dp) :: phase_fac

        real(dp) :: y_vec(3)
            !! y = x + r - x_i

        integer :: s, s_count, n_s
        real(dp), allocatable :: s_red_list(:, :)
        real(dp), allocatable :: s_vec_list(:, :)

        if ( shift ) then
            x_shift_red = eq_pos_red(core_elec_conf(core_state_idx, 1), :)
        else
            x_shift_red = [0.0_dp, 0.0_dp, 0.0_dp]
        end if

        x_shift = matmul(red_to_xyz, x_shift_red)

        if ( present(k_vec_in) ) then

            k_vec = k_vec_in

        else

            k_vec = [0.0_dp, 0.0_dp, 0.0_dp]

        end if

        if ( present(s_cut_in) ) then

            s_cut = s_cut_in

        else

            s_cut = 1

        end if

        n_s = (2*s_cut + 1)**3

        allocate(s_red_list(3, n_s))
        allocate(s_vec_list(3, n_s))

        s_count = 0
        do s1 = -s_cut, s_cut
            do s2 = -s_cut, s_cut
                do s3 = -s_cut, s_cut
                    
                    s_count = s_count + 1

                    s_red_list(:, s_count) = [1.0_dp*s1, 1.0_dp*s2, 1.0_dp*s3]
                    s_vec_list(:, s_count) = matmul(red_to_xyz, s_red_list(:, s_count))

                end do
            end do
        end do

        do n3 = 1, n_grid(3)
            do n2 = 1, n_grid(2)
                do n1 = 1, n_grid(1)

                    x_red = [(n1 - 1.0_dp)/n_grid(1),&
                                (n2 - 1.0_dp)/n_grid(2),&
                                (n3 - 1.0_dp)/n_grid(3)]

                    x_vec = matmul(red_to_xyz, x_red)

                    core_sto_wf_sum = cmplx(0.0_dp, 0.0_dp, dp)

                    do s = 1, n_s

                        phase_fac = exp(ii*dot_product(k_vec, s_vec_list(:, s) - x_vec))
                        y_vec = x_vec - x_shift + s_vec_list(:, s)

                        core_sto_wf_sum = core_sto_wf_sum + &
                            phase_fac*core_sto_wf(core_state_idx, y_vec)

                    end do

                    core_wfc(n1, n2, n3) = core_sto_wf_sum 

                end do
            end do
        end do

    end subroutine

    function core_sto_wf(core_state_idx, x_vec) result(val)
        !! Core wave function
        !!
        !! Units : eV^(3/2)
        implicit none

        integer :: core_state_idx
            !! references a specific element in core_elec_conf
        integer :: l, m
        real(dp) :: x_vec(3), x_hat(3)
        real(dp) :: x_mag

        complex(dp) :: val

        x_mag = norm2(x_vec)

        !! avoid |x| = 0 problems
        if ( x_mag .gt. 1e-8_dp ) then
            x_hat = x_vec/x_mag
        else
            x_hat = [0, 0, 1]
        end if

        l = core_elec_conf(core_state_idx, 3)
        m = core_elec_conf(core_state_idx, 4)

        val = core_sto_wf_radial(core_state_idx, x_mag)*sph_harmonic(l, m,&
                                            get_theta(x_hat), get_phi(x_hat))

    end function

    function core_sto_wf_FT(core_state_idx, k_vec) result(val)
        !! Fourier transform of the core wave function
        !!
        !! Units : eV^(-3/2)
        implicit none

        integer :: core_state_idx
            !! references a specific element in core_elec_conf
        integer :: l, m
        real(dp) :: k_vec(3), k_hat(3)
        real(dp) :: k_mag

        complex(dp) :: val

        k_mag = norm2(k_vec)

        !! avoid |k| = 0 problems
        if ( k_mag .gt. 1e-8_dp ) then
            k_hat = k_vec/k_mag
        else
            k_hat = [0, 0, 1]
        end if

        l = core_elec_conf(core_state_idx, 3)
        m = core_elec_conf(core_state_idx, 4)

        val = core_sto_wf_FT_radial(core_state_idx, k_mag)*sph_harmonic(l, m,&
                                            get_theta(k_hat), get_phi(k_hat))

    end function

    function core_sto_wf_radial(core_state_idx, x_mag) result(val)
        !! Radial part of the total core wave function, summed over 
        !! the individual sto_wf_radial
        !!
        !! Units : eV^(3/2)
        implicit none

        integer :: core_state_idx
            !! references a specific element in core_elec_conf
        integer :: atom, n, l, m

        integer :: j
        integer :: n_lj, nj
        real(dp) :: Z_lj, N0_lj, C_lnj

        real(dp) :: x_mag

        real(dp) :: val

        val = 0.0_dp

        l = core_elec_conf(core_state_idx, 3)

        nj = core_sto_nj_list(core_state_idx)
        do j = 1, nj

            n_lj = int(core_sto_data(core_state_idx, 1, j))
            Z_lj = core_sto_data(core_state_idx, 2, j)
            N0_lj = core_sto_data(core_state_idx, 3, j)
            C_lnj = core_sto_data(core_state_idx, 4, j)

            val = val + C_lnj*sto_wf_radial(n_lj, N0_lj, Z_lj, x_mag) 

        end do

    end function

    function core_sto_wf_FT_radial(core_state_idx, k_mag) result(chi)
        !! Radial part of the total Fourier transformed core wave function, summed over 
        !! the individual sto_wf_FT_radial
        !!
        !! Units : eV^(-3/2)
        implicit none

        integer :: core_state_idx
            !! references a specific element in core_elec_conf
        integer :: atom, n, l, m

        integer :: j
        integer :: n_lj, nj

        real(dp) :: Z_lj, N0_lj, C_lnj

        real(dp) :: k_mag

        complex(dp) :: chi

        chi = cmplx(0.0_dp, 0.0_dp, dp)

        l = core_elec_conf(core_state_idx, 3)

        nj = core_sto_nj_list(core_state_idx)
        do j = 1, nj

            n_lj = int(core_sto_data(core_state_idx, 1, j))
            Z_lj = core_sto_data(core_state_idx, 2, j)
            N0_lj = core_sto_data(core_state_idx, 3, j)
            C_lnj = core_sto_data(core_state_idx, 4, j)

            chi = chi + C_lnj*sto_wf_FT_radial(n_lj, l, N0_lj, Z_lj, k_mag) 

        end do

    end function

    function sto_wf_radial(n, norm, Z, x_mag) result(val)
        !! Radial part of the Slater type orbital (STO) wave function
        !!
        !! Units : eV^(3/2)
        implicit none
        
        integer :: n
        real(dp) :: norm, Z, x_mag

        real(dp) :: val

        val = a0**(-1.5)*norm*(x_mag/a0)**(n - 1)*exp(-Z*x_mag/a0)

    end function

    function sto_wf_FT_radial(n, l, norm, Z, k_mag) result(val)
        !! Radial part of the fourier transform of a Slater type orbital (STO) wave function
        !! 
        !! sto_wf_FT = int d^3x sto_wf e^(-ikr)
        !!           = (sto_wf_FT_radial) * sph_harmonic(k hat)
        !!
        !! Reference : https://en.wikipedia.org/wiki/Slater-type_orbital
        !!
        !! Units : eV^(-3/2)
        implicit none

        integer :: n, l
        real(dp) :: norm, Z, k_mag 

        complex(dp) :: val

        integer :: s
        real(dp) :: omega_s, xi

        val = cmplx(0.0_dp, 0.0_dp, dp)

        xi = Z/a0

        do s = 0, floor((n - l)/2.0_dp)

            omega_s = (-4*xi**2)**(-s)*factorial(n - s)*&
              (factorial(s)*factorial(n - l - 2*s))**(-1)

            val = val + omega_s*(k_mag**2 + xi**2)**(s - n - 1)

        end do

        val = a0**(-0.5)*a0**(-n)*norm*(4*pi)*factorial(n - l)*&
            (2*xi)**n*(ii*k_mag/xi)**l*val

    end function

    subroutine load_core_elec_config(filename, verbose)
        !! Reads the core electron configuration file

        implicit none

        logical, optional :: verbose

        character(len=*) :: filename

        integer(HID_T) :: file_id
        logical :: file_exists

        integer(HSIZE_T) :: dims(1) = [0]
        integer(HSIZE_T) :: dims1(1)
        integer(HSIZE_T) :: dims2(2)

        integer :: error

        if ( verbose ) then

            print*, 'Loading core electron configuration...'
            print*

        end if

        inquire(file = filename, exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)

            call h5ltread_dataset_int_f(file_id, 'n_atoms', n_atoms, dims, error)
            call h5ltread_dataset_int_f(file_id, 'n_core_states', n_core_states, dims, error)

            dims2 = [n_atoms, 3]
            allocate(eq_pos_red(n_atoms, 3))
            call h5ltread_dataset_double_f(file_id, 'eq_pos_red', eq_pos_red, dims2, error)

            dims1 = [n_core_states]
            allocate(core_energy(n_core_states))
            call h5ltread_dataset_double_f(file_id, 'core_energy', core_energy, dims1, error)

            dims1 = [n_atoms]
            allocate(Z_list(n_atoms))
            call h5ltread_dataset_int_f(file_id, 'Z_list', Z_list, dims1, error)

            dims2 = [n_core_states, 5]
            allocate(core_elec_conf(n_core_states, 5))
            call h5ltread_dataset_int_f(file_id, 'core_elec_config', core_elec_conf, dims2, error)

            call h5fclose_f(file_id, error)
            call h5close_f(error)

            call print_core_elec_config(filename, verbose=verbose)

            if ( verbose ) then
                print*, '----------------------------------------'
                print*
            end if

        else

            if ( verbose ) then

                print*, '!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*
                print*, '    Core electron configuration file : '
                print*, '    ', trim(filename)
                print*, '    does NOT exist.'
                print*
                print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*

            end if

            stop

        end if

    end subroutine

    subroutine print_core_elec_config(filename, verbose)
        implicit none

        logical, optional :: verbose
        integer :: n

        character(len=*) :: filename

        if ( verbose ) then

            print*, '----------------------------------------'
            print*, '    ---------------------------'
            print*, '    Core Electron Configuration'
            print*, '    ---------------------------'
            print*
            print*, '        Core electron configuration filename : ', trim(filename) 
            print*
            print*, '        Number of atoms = ', n_atoms
            print*, '        Proton numbers = ', Z_list
            print*, '        Number of core states = ', n_core_states
            print*
            print*, '        Equilibrium positions (reduced) : '
            do n = 1, n_atoms
                print*, '            atom # = ', n, ', ', eq_pos_red(n, :)
            end do
            print*
            print*, '        Core electron configuration : '
            do n = 1, n_core_states
                
                print*, '            atom = ', core_elec_conf(n, 1)
                print*, '            n    = ', core_elec_conf(n, 2)
                print*, '            l    = ', core_elec_conf(n, 3)
                print*, '            m    = ', core_elec_conf(n, 4)
                print*, '            n_s  = ', core_elec_conf(n, 5)
                print* 
                print*, '            energy = ', core_energy(n), 'eV'
                print*

            end do

        end if

    end subroutine

    subroutine load_core_sto_data(filename, verbose)
        !! Reads the sto wf file to get the relevant
        !! coefficients for the electron configuration
        !! specified in the core_elec_conf
        
        implicit none

        character(len=*) :: filename

        logical, optional :: verbose

        integer(HID_T) :: file_id
        logical :: file_exists

        integer(HSIZE_T) :: dims(1)
        integer(HSIZE_T) :: dims1(1)
        integer(HSIZE_T) :: dims2(2)

        character(len=64) :: dset_name

        integer :: error

        integer :: n, nj

        real(dp), allocatable :: buf(:)

        dims = [1]

        if ( verbose ) then

            print*, 'Loading STO wf data file...'
            print*

        end if

        inquire(file = filename, exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)

            allocate(core_sto_nj_list(n_core_states))

            !! get the relevant data for each core state
            do n = 1, n_core_states

                dset_name = get_sto_dataset_str(Z_list(core_elec_conf(n, 1)),&
                                                core_elec_conf(n, 2),&
                                                core_elec_conf(n, 3))

                call h5ltread_dataset_int_f(file_id,&
                    trim(dset_name)//'/nj',&
                    core_sto_nj_list(n),&
                    dims, error)

            end do

            core_sto_max_nj = maxval(core_sto_nj_list)

            allocate(core_sto_data(n_core_states, 4, core_sto_max_nj))
            core_sto_data = 0.0_dp

            do n = 1, n_core_states

                nj = core_sto_nj_list(n)

                dset_name = get_sto_dataset_str(Z_list(core_elec_conf(n, 1)),&
                                                core_elec_conf(n, 2),&
                                                core_elec_conf(n, 3))

                dims1 = [nj]
                allocate(buf(nj))

                call h5ltread_dataset_double_f(file_id,&
                    trim(dset_name)//'/n_lj', buf, dims1, error)
                core_sto_data(n, 1, :nj) = buf

                call h5ltread_dataset_double_f(file_id,&
                    trim(dset_name)//'/Z_lj', buf, dims1, error)
                core_sto_data(n, 2, :nj) = buf

                call h5ltread_dataset_double_f(file_id,&
                    trim(dset_name)//'/N_lj', buf, dims1, error)
                core_sto_data(n, 3, :nj) = buf

                call h5ltread_dataset_double_f(file_id,&
                    trim(dset_name)//'/C_lnj', buf, dims1, error)
                core_sto_data(n, 4, :nj) = buf

                deallocate(buf)

            end do

            call h5fclose_f(file_id, error)
            call h5close_f(error)

            call print_sto_data(verbose=verbose)

            if ( verbose ) then
                print*, '----------------------------------------'
                print*
            end if

        else

            if ( verbose ) then

                print*, '!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*
                print*, '    STO wave function coefficient file : '
                print*, '    ', trim(filename)
                print*, '    does NOT exist.'
                print*
                print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*

            end if

            stop

        end if

    end subroutine

    function get_sto_dataset_str(Z, n, l) result(dset_name)
        !! Returns the name of the dataset for nth element in the
        !! core electron configuration
        implicit none

        integer :: Z, n, l
        character(len=64) :: Z_str, n_str, l_str

        character(len = 64) :: dset_name

        write(Z_str, *) Z
        Z_str = adjustl(Z_str)
        write(n_str, *) n
        n_str = adjustl(n_str)
        write(l_str, *) l
        l_str = adjustl(l_str)

        dset_name = 'Z_'//trim(Z_str)//'/n_'//trim(n_str)//'/l_'//trim(l_str)

    end function

    subroutine print_sto_data(verbose)

        implicit none

        logical, optional :: verbose
        integer :: n
        integer :: a

        if ( verbose ) then

            print*, '----------------------------------------'
            print*, '    -------------------'
            print*, '    STO WF Coefficients'
            print*, '    -------------------'
            print*
            do n = 1, n_core_states

                if ( core_elec_conf(n, 4) == 0 ) then
                
                    print*, '        atom = ', core_elec_conf(n, 1)
                    print*, '        n    = ', core_elec_conf(n, 2)
                    print*, '        l    = ', core_elec_conf(n, 3)
                    ! print*, '    m    = ', core_elec_conf(n, 4)
                    print* 
                    print*, '            Number of coefficients = ', core_sto_nj_list(n)
                    print*, '            n_lj = ', core_sto_data(n, 1, :)
                    print*, '            Z_lj = ', core_sto_data(n, 2, :)
                    print*, '            N0_lj = ', core_sto_data(n, 3, :)
                    print*, '            C_lnj = ', core_sto_data(n, 4, :)
                    print*

                end if

            end do

        end if

    end subroutine

end module
