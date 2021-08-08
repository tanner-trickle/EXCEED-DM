module core_electron_type
    !! Defines the `core_electron` type.

    use hdf5
    use h5lt

    use prec
    use constants
    use units
    use math_mod

    use info_messages

    implicit none

    type core_electron_t
        !! Information about the core electrons for each atom in the unit cell.
        !! Contains the configuration information, e.g. number of atoms, number of core
        !! states for each atom, but also the RHF wave function coefficients to get the
        !! atomic wave functions, \( \psi^\text{atom}_{\kappa, n, l, m} \).
        integer :: n_atom
            !! Number of atoms in the primitive cell.
        integer :: n_state
            !! Number of core states.
        real(dp), allocatable :: eq_pos_red(:, :)
            !! Dim : [n_atom, 3]
            !!
            !! Equilibrium positions in reduced coordinates
            !!
            !! Units : None
        real(dp), allocatable :: energy(:)
            !! Dim : [n_state]
            !!
            !! Energy of each core state.
            !!
            !! Units : eV
        integer, allocatable :: Z(:)
            !! Dim : [n_atom]
            !!
            !! Proton number of each atom in the primitive cell
        integer, allocatable :: config(:, :)
            !! Dim : [n_state, 5]
            !!
            !! Electron configuration, each element is
            !!      config(i, 1) - atom id
            !!      config(i, 2) - n
            !!      config(i, 3) - l
            !!      config(i, 4) - m
            !!      config(i, 5) - n_s (number of spin states)
        integer, allocatable :: STO_nj(:)
            !! Dim : [n_state]
            !!
            !! Number of STO expansion parameters for each state
        real(dp), allocatable :: STO_nl(:, :)
            !! Dim : [ n_state, STO_nj(n_state) ]
            !!
            !! STO n_lj parameters
            !!
            !! Units : None
        real(dp), allocatable :: STO_Zl(:, :)
            !! Dim : [ n_state, STO_nj(n_state) ]
            !!
            !! STO Z_lj parameters
            !!
            !! Units : None
        real(dp), allocatable :: STO_norm_l(:, :)
            !! Dim : [ n_state, STO_nj(n_state) ]
            !!
            !! STO N_lj parameters
            !!
            !! Units : None
        real(dp), allocatable :: STO_Cnl(:, :)
            !! Dim : [ n_state, STO_nj(n_state) ]
            !!
            !! STO C_nlj parameters
            !!
            !! Units : None
        ! generated
        real(dp), allocatable :: eq_pos_xyz(:, :)
            !! Dim : [n_atom, 3]
            !!
            !! Equilibrium positions in xyz coordinates.
            !!
            !! Units : eV^(-1)
        integer :: STO_max_nj
            !! Maximum value of nj across all states.
        character(len=512) :: config_filename
            !! Core electron configuration filename
        character(len=512) :: sto_filename
            !! STO wave function coefficient filename

        contains

            procedure :: load => load_core_electron
            procedure :: save => save_core_electron
            procedure :: print => print_core_electron
            procedure :: bloch_wf_on_grid => core_bloch_wf_on_grid
            procedure :: atomic_sto_wf
            procedure :: atomic_sto_wf_radial
            procedure :: atomic_sto_wf_FT
            procedure :: atomic_sto_wf_FT_radial
        
    end type

contains

        function atomic_sto_wf_FT_radial(self, id, k_mag) result(chi)
        !! Radial part of `atomic_sto_wf`.
        !!
        !! Units : \( \text{eV}^{-3/2} \)
        implicit none

        class(core_electron_t) :: self

        integer :: id
            !! references a specific element in self%config
        integer :: atom, n, l, m

        integer :: j
        integer :: n_lj, nj

        real(dp) :: Z_lj, N0_lj, C_lnj

        real(dp) :: k_mag

        complex(dp) :: chi

        chi = cmplx(0.0_dp, 0.0_dp, dp)

        l = self%config(id, 3)

        do j = 1, self%STO_nj(id)

            n_lj = int(self%STO_nl(id, j))
            Z_lj = self%STO_Zl(id, j)
            N0_lj = self%STO_norm_l(id, j)
            C_lnj = self%STO_Cnl(id, j)

            chi = chi + C_lnj*sto_wf_FT_radial(n_lj, l, N0_lj, Z_lj, k_mag)

        end do

    end function

    function atomic_sto_wf_FT(self, id, k) result ( wf_FT )
        !! Fourier transform of the atomic wave function which is a sum of STO's.
        !!
        !! Units : \( \text{eV}^{-3/2} \)
        implicit none

        class(core_electron_t) :: self

        integer :: id
            !! references a specific element in self%config()
        integer :: l, m
        real(dp) :: k(3), k_hat(3)
        real(dp) :: k_mag

        complex(dp) :: wf_FT

        k_mag = norm2(k)

        ! avoid |k| = 0 problems
        if ( k_mag > 1e-8_dp ) then
            k_hat = k/k_mag
        else
            k_hat = [0, 0, 1]
        end if

        l = self%config(id, 3)
        m = self%config(id, 4)

        wf_FT = self%atomic_sto_wf_FT_radial(id, k_mag)*&
            sph_harmonic(l, m, get_theta(k_hat), get_phi(k_hat))

    end function

    function sto_wf_radial(n, N0, Z, r) result ( wf )
        !! Individual Slater Type Orbital (STO)
        !!
        !! Units : eV^(3/2)
        integer :: n
        real(dp) :: N0, Z, r

        real(dp) :: wf

        wf = a0**(-1.5_dp)*N0*(r/a0)**(n - 1)*exp(-Z*r/a0)

    end function

    function sto_wf_FT_radial(n, l, norm, Z, k_mag) result ( wf_FT )
        !! Radial part of the Fourier transform of a Slater type orbital (STO) wave function.
        !!
        !! $$\begin{align*}
        !!      \widetilde{\psi}_\text{STO}( \mathbf{k} ) & = \int d^3x \psi_\text{STO} e^{-i ( \mathbf{k} \cdot \mathbf{r} ) } \\
        !!           & = \widetilde{\psi}_\text{STO, radial}( k ) Y^l_m(\hat{\mathbf{k}})
        !! \end{align*}$$
        !!
        !! Reference : [https://en.wikipedia.org/wiki/Slater-type_orbital](https://en.wikipedia.org/wiki/Slater-type_orbital)
        !!
        !! Units : \( \text{eV}^{-3/2} \)
        implicit none

        integer :: n, l
        real(dp) :: norm, Z, k_mag

        complex(dp) :: wf_FT

        integer :: s
        real(dp) :: omega_s, xi

        wf_FT = cmplx(0.0_dp, 0.0_dp, dp)

        xi = Z/a0

        do s = 0, floor((n - l)/2.0_dp)

            omega_s = (-4*xi**2)**(-s)*factorial(n - s)*&
              (factorial(s)*factorial(n - l - 2*s))**(-1)

            wf_FT = wf_FT + omega_s*(k_mag**2 + xi**2)**(s - n - 1)

        end do

        wf_FT = a0**(-0.5)*a0**(-n)*norm*(4*pi)*factorial(n - l)*&
            (2*xi)**n*(ii*k_mag/xi)**l*wf_FT

    end function

    function atomic_sto_wf_radial(self, id, r) result( wf )
        !! Radial part of `atomic_sto_wf`.
        !!
        !! Units : \( \text{eV}^{3/2} \)
        
        class(core_electron_t) :: self
        integer :: id
        real(dp) :: r

        integer :: l, j  
        integer :: n_lj
        real(dp) :: Z_lj, N0_lj, C_lnj

        real(dp) :: wf

        wf = 0.0_dp

        ! l = self%config(id, 3)

        do j = 1, self%STO_nj(id)

            n_lj = int(self%STO_nl(id, j))
            Z_lj = self%STO_Zl(id, j)
            N0_lj = self%STO_norm_l(id, j)
            C_lnj = self%STO_Cnl(id, j)

            wf = wf + C_lnj*sto_wf_radial(n_lj, N0_lj, Z_lj, r) 

        end do

    end function

    function atomic_sto_wf(self, id, x) result( wf )
        !! Atomic wave function, summed over individual STO's.
        !!
        !! Units : \( \text{eV}^{3/2} \)
        implicit none

        class(core_electron_t) :: self
        integer :: id
        real(dp) :: x(3)

        real(dp) :: x_hat(3)
        real(dp) :: x_mag

        complex(dp) :: wf

        integer :: l, m

        x_mag = norm2(x)

        ! avoid |x| = 0 problems
        if ( x_mag > 1e-8_dp ) then
            x_hat = x/x_mag
        else
            x_hat = [0, 0, 1]
        end if

        l = self%config(id, 3)
        m = self%config(id, 4)

        wf = self%atomic_sto_wf_radial(id, x_mag)&
            *sph_harmonic(l, m, get_theta(x_hat), get_phi(x_hat))

    end function

    subroutine core_bloch_wf_on_grid(self, n_grid, wf, id, pc_vol, &
            red_to_xyz, shift, k_vec_in, r_cut_in, verbose)
        !! Compute the core electron wave functions on a grid inside the primitive cell, summed
        !! over neighboring unit cells.
        !!
        !! \( \text{ID} \rightarrow \{ \kappa, n, l, m \} \)
        !! $$\begin{align} 
        !! u^\text{core}_{\text{ID}, \mathbf{k}}(\mathbf{x}) = \sqrt{\Omega} \sum_{r}
        !! e^{-i \mathbf{k} \left( \mathbf{x} - \mathbf{r} - \mathbf{x}_\kappa \right) }
        !! \psi^\text{atom}( \mathbf{x} - \mathbf{r} - \mathbf{x}_\kappa )
        !! \nonumber \end{align}$$
        !!
        !! `wf` \( \rightarrow u^\text{core}_{\text{ID}, \mathbf{k}} \)
        !!
        !! By default, will compute the above assuming \( \mathbf{k} = 0 \) and sum over only the nearest 
        !! neighbor unit cells.
        !!
        !! Units : None
        implicit none

        class(core_electron_t) :: self
        integer :: n_grid(3)
        complex(dp) :: wf(n_grid(1), n_grid(2), n_grid(3))
        integer :: id
        real(dp) :: pc_vol
        real(dp) :: red_to_xyz(3, 3)
        logical, optional :: shift
        real(dp), optional :: k_vec_in(3)
        integer, optional :: r_cut_in
        logical, optional :: verbose

        real(dp) :: k_vec(3)
        integer :: r_cut

        integer :: n1, n2, n3
        integer :: r1, r2, r3

        real(dp) :: x_red(3)
        real(dp) :: r_red(3)
        real(dp) :: x_vec(3)
        real(dp) :: x_shift_red(3)
        real(dp) :: x_shift(3)

        complex(dp) :: phase_fac

        real(dp) :: y_vec(3)
            !! y = x + r - x_i

        integer :: r, r_count, n_r
        real(dp), allocatable :: r_red_list(:, :)
        real(dp), allocatable :: r_vec_list(:, :)

        complex(dp) :: wf_sum

        if ( shift ) then
            x_shift_red = self%eq_pos_red(self%config(id, 1), :)
        else
            x_shift_red = [0.0_dp, 0.0_dp, 0.0_dp]
        end if

        x_shift = matmul(red_to_xyz, x_shift_red)

        if ( present(k_vec_in) ) then
            k_vec = k_vec_in
        else
            k_vec = [0.0_dp, 0.0_dp, 0.0_dp]
        end if

        if ( present(r_cut_in) ) then
            r_cut = r_cut_in
        else
            r_cut = 1
        end if

        n_r = (2*r_cut + 1)**3

        allocate(r_red_list(3, n_r))
        allocate(r_vec_list(3, n_r))

        r_count = 0
        do r1 = -r_cut, r_cut
            do r2 = -r_cut, r_cut
                do r3 = -r_cut, r_cut
                    
                    r_count = r_count + 1

                    r_red_list(:, r_count) = [1.0_dp*r1, 1.0_dp*r2, 1.0_dp*r3]
                    r_vec_list(:, r_count) = matmul(red_to_xyz, r_red_list(:, r_count))

                end do
            end do
        end do

        if ( present(k_vec_in) ) then

            do n3 = 1, n_grid(3)
                do n2 = 1, n_grid(2)
                    do n1 = 1, n_grid(1)

                        x_red = [&
                                    (n1 - 1.0_dp)/n_grid(1),&
                                    (n2 - 1.0_dp)/n_grid(2),&
                                    (n3 - 1.0_dp)/n_grid(3)&
                                ]

                        x_vec = matmul(red_to_xyz, x_red)

                        wf_sum = cmplx(0.0_dp, 0.0_dp, dp)

                        do r = 1, n_r

                            y_vec = x_vec - x_shift + r_vec_list(:, r)
                            phase_fac = exp(-ii*dot_product(k_vec, y_vec))

                            wf_sum = wf_sum + phase_fac*self%atomic_sto_wf(id, y_vec)

                        end do

                        wf(n1, n2, n3) = wf_sum

                    end do
                end do
            end do

        else

            do n3 = 1, n_grid(3)
                do n2 = 1, n_grid(2)
                    do n1 = 1, n_grid(1)

                        x_red = [&
                                    (n1 - 1.0_dp)/n_grid(1),&
                                    (n2 - 1.0_dp)/n_grid(2),&
                                    (n3 - 1.0_dp)/n_grid(3)&
                                ]

                        x_vec = matmul(red_to_xyz, x_red)

                        wf_sum = cmplx(0.0_dp, 0.0_dp, dp)

                        do r = 1, n_r

                            y_vec = x_vec - x_shift + r_vec_list(:, r)

                            wf_sum = wf_sum + self%atomic_sto_wf(id, y_vec)

                        end do

                        wf(n1, n2, n3) = wf_sum

                    end do
                end do
            end do

        end if

        wf = sqrt(pc_vol)*wf

    end subroutine

    subroutine save_core_electron(self, filename, verbose)
        !! Saves `core_electron`.

        implicit none

        class(core_electron_t) :: self
        character(len=*) :: filename
        logical, optional :: verbose

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id
        logical :: file_exists
        integer(HSIZE_T) :: dims1(1) = [1]
        integer(HSIZE_T) :: dims2(2)
        integer :: error

        if ( verbose ) then
            print*, 'Saving core electron parameters...'
            print*
        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'core_electron', group_id, error)

            ! write data
            call h5ltmake_dataset_int_f(file_id, &
                'core_electron/n_atom', &
                size(dims1), dims1, &
                self%n_atom, &
                error)
            call h5ltmake_dataset_int_f(file_id, &
                'core_electron/n_state', &
                size(dims1), dims1, &
                self%n_state, &
                error)

            dims2 = [self%n_atom, 3]
            call h5ltmake_dataset_double_f(file_id, &
                'core_electron/eq_pos_red', &
                size(dims2), dims2, &
                self%eq_pos_red, &
                error)

            dims1 = [self%n_state]
            call h5ltmake_dataset_double_f(file_id, &
                'core_electron/energy', &
                size(dims1), dims1, &
                self%energy, &
                error)

            dims1 = [self%n_atom]
            call h5ltmake_dataset_int_f(file_id, &
                'core_electron/Z', &
                size(dims1), dims1, &
                self%Z, &
                error)

            dims2 = [self%n_state, 5]
            call h5ltmake_dataset_int_f(file_id, &
                'core_electron/config', &
                size(dims2), dims2, &
                self%config, &
                error)

            ! STO parameters

            call h5gcreate_f(file_id, 'core_electron/STO', group_id, error)

            dims1 = [self%n_state]
            call h5ltmake_dataset_int_f(file_id, &
                'core_electron/STO/n_j', &
                size(dims1), dims1, &
                self%STO_nj, &
                error)

            dims1 = [1]
            call h5ltmake_dataset_int_f(file_id, &
                'core_electron/STO/max_n_j', &
                size(dims1), dims1, &
                self%STO_max_nj, &
                error)

            dims2 = [self%n_state, self%STO_max_nj]
            call h5ltmake_dataset_double_f(file_id, &
                'core_electron/STO/n_lj', &
                size(dims2), dims2, &
                self%STO_nl, &
                error)

            call h5ltmake_dataset_double_f(file_id, &
                'core_electron/STO/Z_lj', &
                size(dims2), dims2, &
                self%STO_Zl, &
                error)

            call h5ltmake_dataset_double_f(file_id, &
                'core_electron/STO/N0_lj', &
                size(dims2), dims2, &
                self%STO_norm_l, &
                error)

            call h5ltmake_dataset_double_f(file_id, &
                'core_electron/STO/C_nlj', &
                size(dims2), dims2, &
                self%STO_Cnl, &
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

    subroutine print_core_electron(self, verbose)
        !! Prints `core_electron` components.
        implicit none

        class(core_electron_t) :: self
        logical, optional :: verbose

        integer :: n

        if ( verbose ) then
            call print_section_seperator()
            print*, '    --------------'
            print*, '    Core Electrons'
            print*, '    --------------'
            print*
            print*, '        Core electron configuration filename : ', trim(self%config_filename) 
            print*, '        STO data filename                    : ', trim(self%STO_filename)
            print*
            print*, '        Number of atoms       : ', trim(adjustl(int_to_str(self%n_atom)))
            print*, '        Proton numbers        : ', self%Z
            print*, '        Number of core states : ', trim(adjustl(int_to_str(self%n_state)))
            print*
            print*, '        Equilibrium positions (reduced) : '
            do n = 1, self%n_atom
                print*, '            atom # : ', n, ', ', self%eq_pos_red(n, :)
            end do
            print*
            print*, '        Configuration and STO parameters : '
            print*
            do n = 1, self%n_state
                print*, '            atom # : ', trim(adjustl(int_to_str(self%config(n, 1))))
                print*, '            n      : ', trim(adjustl(int_to_str(self%config(n, 2))))
                print*, '            l      : ', trim(adjustl(int_to_str(self%config(n, 3))))
                print*, '            m      : ', trim(adjustl(int_to_str(self%config(n, 4))))
                print*, '            n_s    : ', trim(adjustl(int_to_str(self%config(n, 5))))
                print* 
                print*, '            energy : ', self%energy(n), 'eV'
                print*
                print*, '            Number of coefficients : ', trim(adjustl(int_to_str(self%STO_nj(n))))
                print*
                print*, '            n_lj : ', self%STO_nl(n, :)
                print*, '            Z_lj : ', self%STO_Zl(n, :)
                print*, '            N0_lj : ', self%STO_norm_l(n, :)
                print*, '            C_lnj : ', self%STO_Cnl(n, :)
                print*
            end do
            call print_section_seperator()
            print*
        end if

    end subroutine

    subroutine load_core_electron(self, core_elec_config_filename, STO_filename, verbose)
        !! Loads `core_electron` parameters from a namelist.
        implicit none

        class(core_electron_t) :: self
        character(len=*) :: core_elec_config_filename
        character(len=*) :: STO_filename
        logical, optional :: verbose

        integer(HID_T) :: file_id
        logical :: core_file_exists, STO_file_exists

        integer(HSIZE_T) :: dims(1) = [0]
        integer(HSIZE_T) :: dims1(1)
        integer(HSIZE_T) :: dims2(2)

        real(dp), allocatable :: buf(:)
        character(len=64) :: dset_name

        integer :: nj, n

        integer :: error

        if ( verbose ) then
            print*, 'Loading core electron configuration...'
            print*
        end if

        inquire(file = core_elec_config_filename, exist = core_file_exists)

        if ( core_file_exists ) then

            self%config_filename = trim(core_elec_config_filename)

            call h5open_f(error)
            call h5fopen_f(core_elec_config_filename, H5F_ACC_RDONLY_F, file_id, error)

            call h5ltread_dataset_int_f(file_id, 'n_atom', self%n_atom, dims, error)
            call h5ltread_dataset_int_f(file_id, 'n_state', self%n_state, dims, error)

            dims2 = [self%n_atom, 3]
            allocate(self%eq_pos_red(self%n_atom, 3))
            call h5ltread_dataset_double_f(file_id, 'eq_pos_red', self%eq_pos_red, dims2, error)

            dims1 = [self%n_state]
            allocate(self%energy(self%n_state))
            call h5ltread_dataset_double_f(file_id, 'energy', self%energy, dims1, error)

            dims1 = [self%n_atom]
            allocate(self%Z(self%n_atom))
            call h5ltread_dataset_int_f(file_id, 'Z', self%Z, dims1, error)

            dims2 = [self%n_state, 5]
            allocate(self%config(self%n_state, 5))
            call h5ltread_dataset_int_f(file_id, 'config', self%config, dims2, error)

            call h5fclose_f(file_id, error)
            call h5close_f(error)

        else

            call print_error_message('Core electron configuration file'//&
               trim(core_elec_config_filename)//'does NOT exist.', &
                verbose = verbose)
            stop

        end if

        if ( verbose ) then
            print*, 'Loading STO coefficients...'
            print*
        end if

        inquire(file = STO_filename, exist = STO_file_exists)

        if ( STO_file_exists ) then

            self%STO_filename = trim(STO_filename)

            call h5open_f(error)
            call h5fopen_f(STO_filename, H5F_ACC_RDONLY_F, file_id, error)

            allocate(self%STO_nj(self%n_state))

            !! get the relevant data for each core state
            do n = 1, self%n_state

                dset_name = get_sto_dataset_str(self%Z(self%config(n, 1)),&
                                                self%config(n, 2),&
                                                self%config(n, 3))

                call h5ltread_dataset_int_f(file_id,&
                    trim(dset_name)//'/nj',&
                    self%STO_nj(n),&
                    dims, error)

            end do

            self%STO_max_nj = maxval(self%STO_nj)

            allocate(self%STO_nl(self%n_state, self%STO_max_nj))
            allocate(self%STO_Zl(self%n_state, self%STO_max_nj))
            allocate(self%STO_norm_l(self%n_state, self%STO_max_nj))
            allocate(self%STO_Cnl(self%n_state, self%STO_max_nj))

            self%STO_nl = 0.0_dp
            self%STO_Zl = 0.0_dp
            self%STO_norm_l = 0.0_dp
            self%STO_Cnl = 0.0_dp

            do n = 1, self%n_state

                nj = self%STO_nj(n)

                dset_name = get_sto_dataset_str(self%Z(self%config(n, 1)),&
                                                self%config(n, 2),&
                                                self%config(n, 3))

                dims1 = [nj]
                allocate(buf(nj))

                call h5ltread_dataset_double_f(file_id,&
                    trim(dset_name)//'/n_lj', buf, dims1, error)
                self%STO_nl(n, :nj) = buf

                call h5ltread_dataset_double_f(file_id,&
                    trim(dset_name)//'/Z_lj', buf, dims1, error)
                self%STO_Zl(n, :nj) = buf

                call h5ltread_dataset_double_f(file_id,&
                    trim(dset_name)//'/N_lj', buf, dims1, error)
                self%STO_norm_l(n, :nj) = buf

                call h5ltread_dataset_double_f(file_id,&
                    trim(dset_name)//'/C_lnj', buf, dims1, error)
                self%STO_Cnl(n, :nj) = buf

                deallocate(buf)

            end do

            call h5fclose_f(file_id, error)
            call h5close_f(error)

        else

            call print_error_message('STO wave function coefficient file'//&
                trim(STO_filename)//' does NOT exist.', verbose = verbose)
            stop

        end if

        call self%print(verbose = verbose)

    end subroutine

    function get_sto_dataset_str(Z, n, l) result(dset_name)
        !! Returns the (string) name of the dataset where the RHF STO coefficients are for atom with 
        !! proton number, \( Z \), and quantum numbers, \( n, l \) inside the RHF data file.
        implicit none

        integer :: Z, n, l
        character(len = 64) :: dset_name

        dset_name = 'Z_'//trim(adjustl(int_to_str(Z)))&
            //'/n_'//trim(adjustl(int_to_str(n)))&
            //'/l_'//trim(adjustl(int_to_str(l)))

    end function

end module
