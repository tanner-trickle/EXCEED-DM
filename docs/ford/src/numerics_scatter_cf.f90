module numerics_scatter_cf
    !! Numerics parameters specific to the core \( \rightarrow \) free DM-electron scattering rate calculation.

    use prec
    use info_messages

    use core_electron_type

    implicit none

    type numerics_scatter_cf_t
        !! Numerics parameters specific to the core \( \rightarrow \) free DM-electron scattering rate calculation.

        integer :: n_principal_max = 10
            !! Maximum principal quantum number to include.
        integer :: n_principal_min = 1
            !! Minimum principal quantum number to include.
        character(len=64) :: Zeff_type = 'one'
            !! Specify what Zeff to use
            !!
            !! - 'one' - all Zeff = 1
            !!
            !! - 'Eb' - use the binding energy of the core states
        real(dp), allocatable :: Zeff(:)
            !! Dim : [ n_state ]
            !!
            !! \( Z_\text{eff} \) parameter in the Fermi form factor.
        integer :: n_kf_theta = 1
            !! Number of \( \theta \) points in the integration of \( \mathbf{k}_f \)
        integer :: n_kf_phi = 1
            !! Number of \( \phi \) points in the integration of \( \mathbf{k}_f \)
        real(dp), allocatable :: kf_angular_mesh(:, :)
            !! Dim : [ n_kf_theta*n_kf_phi, 2 ]
            !!
            !! Grid of \( \theta, \phi \) points for the integration over \( \mathbf{k}_f \)
            !!
            !! Units : None
        integer :: n_omega
            !! Number of \( \omega \) points to compute for. The core -> free scattering
            !! rate calculation is done by computing `n_omega` values of \( \frac{dR}{d \omega} \) between
            !!
            !! \( \omega_\text{min} = E_{f, \text{max}} \)
            !! \( \omega_\text{max} = \frac{1}{2} m_\chi v_\text{max}^2 \)
            !!
            !! and then interpolated to the user defined grid.

        real(dp), allocatable :: omega_list(:)
            !! List of \( \omega \) points to compute for. The core -> free scattering
            !! rate calculation is done by computing `n_omega` values of \( \frac{dR}{d \omega} \) between
            !!
            !! \( \omega_\text{min} = E_{f, \text{max}} \)
            !! \( \omega_\text{max} = \frac{1}{2} m_\chi v_\text{max}^2 \)
            !!
            !! and then interpolated to the user defined grid.

        integer :: n_ki = 2
            !! Number of \( | \mathbf{k}_i | \) points in integration over \( \mathbf{k}_i \)
        integer :: n_ki_theta = 1
            !! Number of \( \theta \) points in integration over \( \mathbf{k}_i \)
        integer :: n_ki_phi = 1
            !! Number of \( \phi \) points in integration over \( \mathbf{k}_i \)
        real(dp), allocatable :: ki_angular_mesh(:, :)
            !! Dim : [ n_ki_theta*n_ki_phi, 2 ]
            !!
            !! Grid of \( \theta, \phi \) points for the integration over \( \mathbf{k}_i \)
            !!
            !! Units : None
        real(dp) :: ki_s = 100.0_dp
            !! Scale parameter for maximum \( |k_i| \) to integrate over.
            !!
            !! `ki_max` = `ki_s` \( \times  Z \alpha m_e \)
            !!
            !! Generally, this parameters should be \( \gg 1 \) since \( Z \alpha m_e \) is the
            !! typical momentum scale of the electronic wave functions.
            !!
            !! Units : None
        real(dp) :: ki_min = 1.0e3_dp
            !! Minimum \( | \mathbf{k}_i | \) in the integration over \( \mathbf{k}_i \)
            !!
            !! Units : eV
        integer, allocatable :: core_id_list(:)
            !! Dim : [n_init]
            !!
            !! List of core ID's to compute for.
        integer, allocatable :: w_id_list(:)
            !! Dim : [n_omega]
            !!
            !! List of omega ID's to compute for.
        real(dp) :: Ef_min
            !! Smallest final electron energy to compute for.

        contains

            procedure :: print => numerics_scatter_cf_print
            procedure :: load => numerics_scatter_cf_load_nml
            procedure :: save => numerics_scatter_cf_save
            procedure :: create_core_id_list => cf_create_core_id_list
            procedure :: create_w_id_list => cf_create_w_id_list

    end type

contains

    subroutine numerics_scatter_cf_save(self, filename, verbose)
        !! Save `numerics_scatter_cf`.
        use hdf5
        use h5lt

        implicit none

        class(numerics_scatter_cf_t) :: self
        character(len=*) :: filename
        logical, optional :: verbose

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id
        logical :: file_exists
        integer(HSIZE_T) :: dims1(1) = [1]
        integer(HSIZE_T) :: dims2(2)
        integer :: error

        if ( verbose ) then
            print*, 'Saving numerics - scatter-cf parameters...'
            print*
        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'numerics_scatter_cf', group_id, error)

            ! write data
            call h5ltmake_dataset_int_f(file_id, &
                'numerics_scatter_cf/n_principal_max', &
                size(dims1), dims1, &
                self%n_principal_max, &
                error)
            call h5ltmake_dataset_int_f(file_id, &
                'numerics_scatter_cf/n_principal_min', &
                size(dims1), dims1, &
                self%n_principal_min, &
                error)

            call h5ltmake_dataset_int_f(file_id, &
                'numerics_scatter_cf/n_kf_theta', &
                size(dims1), dims1, &
                self%n_kf_theta, &
                error)
            call h5ltmake_dataset_int_f(file_id, &
                'numerics_scatter_cf/n_kf_phi', &
                size(dims1), dims1, &
                self%n_kf_phi, &
                error)

            call h5ltmake_dataset_int_f(file_id, &
                'numerics_scatter_cf/n_ki_theta', &
                size(dims1), dims1, &
                self%n_ki_theta, &
                error)
            call h5ltmake_dataset_int_f(file_id, &
                'numerics_scatter_cf/n_ki_phi', &
                size(dims1), dims1, &
                self%n_ki_phi, &
                error)
            call h5ltmake_dataset_int_f(file_id, &
                'numerics_scatter_cf/n_ki', &
                size(dims1), dims1, &
                self%n_ki, &
                error)
            call h5ltmake_dataset_double_f(file_id, &
                'numerics_scatter_cf/ki_s', &
                size(dims1), dims1, &
                self%ki_s, &
                error)
            call h5ltmake_dataset_double_f(file_id, &
                'numerics_scatter_cf/ki_min', &
                size(dims1), dims1, &
                self%ki_min, &
                error)

            call h5ltmake_dataset_double_f(file_id, &
                'numerics_scatter_cf/Ef_min', &
                size(dims1), dims1, &
                self%Ef_min, &
                error)

            dims1 = [size(self%core_id_list)] 
            call h5ltmake_dataset_int_f(file_id, &
                'numerics_scatter_cf/core_id_list', &
                size(dims1), dims1, &
                self%core_id_list, &
                error)

            dims1 = [size(self%Zeff)] 
            call h5ltmake_dataset_double_f(file_id, &
                'numerics_scatter_cf/Zeff', &
                size(dims1), dims1, &
                self%Zeff, &
                error)

            dims1 = [size(self%omega_list)] 
            call h5ltmake_dataset_double_f(file_id, &
                'numerics_scatter_cf/omega_list', &
                size(dims1), dims1, &
                self%omega_list, &
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

    subroutine numerics_scatter_cf_print(self, verbose)
        !! Prints `numerics_scatter_cf` components.

        implicit none

        class(numerics_scatter_cf_t) :: self
        logical, optional :: verbose

        if ( verbose ) then

            call print_section_seperator()
            print*, '    -----------------------'
            print*, '    Numerics - Scatter - cf'
            print*, '    -----------------------'
            print*
            print*, '        Maximum principal quantum number   : ', trim(adjustl(int_to_str(self%n_principal_max)))
            print*, '        Minimum principal quantum number   : ', trim(adjustl(int_to_str(self%n_principal_min)))
            print*
            print*, '        ki Integration'
            print*, '            Number of ki theta points : ', trim(adjustl(int_to_str(self%n_ki_theta)))
            print*, '            Number of ki phi points   : ', trim(adjustl(int_to_str(self%n_ki_phi)))
            print*, '            Number of |ki| points     : ', trim(adjustl(int_to_str(self%n_ki)))
            print*, '            ki scale parameter,       : ', self%ki_s 
            print*, '            Minimum |ki|              : ', self%ki_min, 'eV'
            print*
            print*, '        kf Integration'
            print*, '            Number of kf theta points : ', trim(adjustl(int_to_str(self%n_kf_theta)))
            print*, '            Number of kf phi points   : ', trim(adjustl(int_to_str(self%n_kf_phi)))
            print*
            print*, '        Number of omega points : ', trim(adjustl(int_to_str(self%n_omega)))
            print*, '        Minimum Ef             : ', self%Ef_min, 'eV'
            print*
            print*, '        Zeff type : ', trim(self%Zeff_type)
            print*
            call print_section_seperator()
            print*

        end if

    end subroutine

    subroutine numerics_scatter_cf_load_nml(self, filename, &
            core_electron, verbose)
        !! Loads `numerics_scatter_cf` parameters from a namelist.
        
        implicit none

        class(numerics_scatter_cf_t) :: self
        character(len=*) :: filename
        type(core_electron_t) :: core_electron
        logical, optional :: verbose

        logical :: file_exists
        integer :: error

        integer :: n_principal_max = 10
            !! Maximum number of valence bands
        integer :: n_principal_min = 1

        character(len=64) :: Zeff_type = 'one'
            !! Specify what Zeff to use
            !!
            !! - 'one' - all Zeff = 1
            !! - 'Eb' - use the binding energy of the core states
        integer :: n_kf_theta = 1
            !! Number of \( \theta \) points in the integration of \( \mathbf{k}_f \)
        integer :: n_kf_phi = 1
            !! Number of \( \phi \) points in the integration of \( \mathbf{k}_f \)
        integer :: n_omega
            !! Number of \( \omega \) points to compute for. The core -> free scattering
            !! rate calculation is done by computing `n_omega` values of \( \frac{dR}{d \omega} \) between
            !!
            !! \( \omega_\text{min} = E_{f, \text{max}} \)
            !! \( \omega_\text{max} = \frac{1}{2} m_\chi v_\text{max}^2 \)
            !!
            !! and then interpolated to the user defined grid.

        integer :: n_ki = 2
            !! Number of \( | \mathbf{k}_i | \) points in integration over \( \mathbf{k}_i \)
        integer :: n_ki_theta = 1
            !! Number of \( \theta \) points in integration over \( \mathbf{k}_i \)
        integer :: n_ki_phi = 1
            !! Number of \( \phi \) points in integration over \( \mathbf{k}_i \)
        real(dp) :: ki_s = 100.0_dp
            !! Scale parameter for maximum \( |k_i| \) to integrate over.
            !!
            !! `ki_max` = `ki_s` \( \times  Z \alpha m_e \)
            !!
            !! Generally, this parameters should be \( \gg 1 \) since \( Z \alpha m_e \) is the
            !! typical momentum scale of the electronic wave functions.
            !!
            !! Units : None
        real(dp) :: ki_min = 1.0e3_dp
            !! Minimum \( | \mathbf{k}_i | \) in the integration over \( \mathbf{k}_i \)
            !!
            !! Units : eV
        real(dp) :: Ef_min = 60.0_dp

        integer :: i

        NAMELIST /numerics_s_cf/ n_principal_max, &
                                 n_principal_min, &
                                 Zeff_type, &
                                 n_kf_theta, &
                                 n_kf_phi, &
                                 n_omega, &
                                 n_ki, &
                                 n_ki_theta, &
                                 n_ki_phi, &
                                 ki_s, &
                                 ki_min, &
                                 Ef_min

        if ( verbose ) then
            print*, 'Loading numerics parameters for scatter-cf calculation...'
            print*
        end if

        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            open(100, file = trim(filename), iostat = error)
            read(100, nml = numerics_s_cf, iostat = error)
            close(100)

            self%n_principal_max = min(n_principal_max, 10)
            self%n_principal_min = max(n_principal_min, 1)

            self%Zeff_type = trim(Zeff_type)
            
            self%n_kf_theta = n_kf_theta
            self%n_kf_phi = n_kf_phi
            self%n_omega = n_omega
            
            self%n_ki = n_ki
            self%n_ki_theta = n_ki_theta
            self%n_ki_phi = n_ki_phi
            self%ki_s = ki_s
            self%ki_min = ki_min

            self%Ef_min = Ef_min

            allocate(self%omega_list(self%n_omega))
            self%omega_list = 0.0_dp
            allocate(self%Zeff(core_electron%n_state))
            self%Zeff = 1.0_dp
            allocate(self%kf_angular_mesh(self%n_kf_theta*self%n_kf_phi, 2))
            allocate(self%ki_angular_mesh(self%n_ki_theta*self%n_ki_phi, 2))

            self%kf_angular_mesh = generate_uniform_points_on_sphere(&
                self%n_kf_theta, self%n_kf_phi)
            self%ki_angular_mesh = generate_uniform_points_on_sphere(&
                self%n_ki_theta, self%n_ki_phi)

            if ( trim(self%Zeff_type) == 'Eb' ) then

                do i = 1, core_electron%n_state
                
                    self%Zeff(i) = core_electron%config(i, 2)*sqrt(&
                        ( maxval(core_electron%energy) &
                            - core_electron%energy(i) )&
                        /13.6_dp) 

                        self%Zeff(i) = max(self%Zeff(i), 1.0_dp)

                end do

            end if

            call self%print(verbose = verbose)

        else

            call print_error_message(&
                'Input file for numerics-scatter-cf parameters : '//trim(filename)//' does NOT exist.', &
                verbose = verbose)
            stop

        end if

    end subroutine

    subroutine cf_create_w_id_list(self)
        !! Specify the indicies for each \( \omega \) point that should be included. 
        !! Specific to the core \( \rightarrow \) free DM-electron scattering rate calculation. 
        implicit none

        class(numerics_scatter_cf_t) :: self

        integer :: i

        allocate(self%w_id_list(self%n_omega))

        do i = 1, self%n_omega
            self%w_id_list(i) = i
        end do

    end subroutine

    subroutine cf_create_core_id_list(self, core_electron)
        !! Create the list of core id's to compute for.

        implicit none

        class(numerics_scatter_cf_t) :: self
        type(core_electron_t) :: core_electron

        integer :: n
        integer :: n_init, init_id

        init_id = 0
        do n = 1, core_electron%n_state

            if ( ( core_electron%config(n, 2) >= self%n_principal_min ) &
           .and. ( core_electron%config(n, 2) <= self%n_principal_max ) ) then
                
                init_id = init_id + 1

            end if

        end do

        n_init = init_id

        allocate(self%core_id_list(n_init))

        init_id = 0
        do n = 1, core_electron%n_state

            if ( ( core_electron%config(n, 2) >= self%n_principal_min ) &
           .and. ( core_electron%config(n, 2) <= self%n_principal_max ) ) then
                
                init_id = init_id + 1
                self%core_id_list(init_id) = n

            end if

        end do

    end subroutine

end module
