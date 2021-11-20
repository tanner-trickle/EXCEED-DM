module numerics_scatter_cc_ext

    use info_messages

    use core_electron_type

    type numerics_scatter_cc_ext_t
        !! Numerics parameters specific to the core \( \rightarrow \) conduction DM-electron scattering rate calculation.

        integer :: n_principal_max = 10
            !! Maximum principal quantum number to include.
        integer :: n_principal_min = 1
            !! Minimum principal quantum number to include.
        integer :: n_cond_max = 0
            !! Maximum number of conduction bands.

        integer, allocatable :: core_id_list(:)
            !! Dim : [n_init]
            !!
            !! List of core ID's to compute for.
        integer, allocatable :: k_id_list(:)
            !! Dim : [n_init]
            !!
            !! List of k points to compute for.

        integer :: n_q

        integer :: n_q_theta = 1

        integer :: n_q_phi = 1

        real(dp), allocatable :: q_angular_mesh(:, :)
            !! Dim : [ `n_q_theta` x `n_kf_phi`, 2 ]
            !!
            !! Grid of \( \theta, \phi \) points for the angular integration over \( \mathbf{q} \)
            !!
            !! Units : None

        contains

            procedure :: print => numerics_scatter_cc_ext_print
            procedure :: load => numerics_scatter_cc_ext_load_nml
            procedure :: save => numerics_scatter_cc_ext_save
            procedure :: create_k_id_list => cc_ext_create_k_id_list
            procedure :: create_core_id_list => cc_ext_create_core_id_list

    end type

contains

    subroutine cc_ext_create_core_id_list(self, core_electron)
        !! Specify the indicies for each core state that should be included. 
        !! Specific to the  core \( \rightarrow \) conduction (extended) DM-electron scattering rate calculation. 

        implicit none

        class(numerics_scatter_cc_ext_t) :: self
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

    subroutine cc_ext_create_k_id_list(self, n_k)
        !! Specify the indicies for each \( \mathbf{k} \) point that should be included. 
        !! Specific to the core \( \rightarrow \) conduction DM-electron scattering rate calculations. 

        implicit none

        class(numerics_scatter_cc_ext_t) :: self
        integer :: n_k

        integer :: k

        allocate(self%k_id_list(n_k))

        do k = 1, n_k
            self%k_id_list(k) = k
        end do

    end subroutine

    subroutine numerics_scatter_cc_ext_save(self, filename, verbose)
        !! Saves `numerics_scatter_cc_ext`.

        use hdf5
        use h5lt

        implicit none

        class(numerics_scatter_cc_ext_t) :: self
        character(len=*) :: filename
        logical, optional :: verbose

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id
        logical :: file_exists
        integer(HSIZE_T) :: dims1(1) = [1]
        integer(HSIZE_T) :: dims2(2)
        integer :: error

        if ( verbose ) then
            print*, 'Saving numerics - scatter-cc_ext parameters...'
            print*
        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'numerics_scatter_cc_ext', group_id, error)

            ! write data
            call h5ltmake_dataset_int_f(file_id, &
                'numerics_scatter_cc_ext/n_principal_max', &
                size(dims1), dims1, &
                self%n_principal_max, &
                error)
            call h5ltmake_dataset_int_f(file_id, &
                'numerics_scatter_cc_ext/n_principal_min', &
                size(dims1), dims1, &
                self%n_principal_min, &
                error)
            call h5ltmake_dataset_int_f(file_id, &
                'numerics_scatter_cc_ext/n_cond_max', &
                size(dims1), dims1, &
                self%n_cond_max, &
                error)

            call h5ltmake_dataset_int_f(file_id, &
                'numerics_scatter_cc_ext/n_q', &
                size(dims1), dims1, &
                self%n_q, &
                error)
            call h5ltmake_dataset_int_f(file_id, &
                'numerics_scatter_cc_ext/n_q_theta', &
                size(dims1), dims1, &
                self%n_q_theta, &
                error)
            call h5ltmake_dataset_int_f(file_id, &
                'numerics_scatter_cc_ext/n_q_phi', &
                size(dims1), dims1, &
                self%n_q_phi, &
                error)

            dims1 = [self%n_q]
            call h5ltmake_dataset_double_f(file_id, &
                'numerics_scatter_cc_ext/q_list', &
                size(dims1), dims1, &
                self%q_list, &
                error)

            dims1 = [size(self%core_id_list)] 
            call h5ltmake_dataset_int_f(file_id, &
                'numerics_scatter_cc_ext/core_id_list', &
                size(dims1), dims1, &
                self%core_id_list, &
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

    subroutine numerics_scatter_cc_ext_load_nml(self, filename, n_cond, verbose)
        !! Loads `numerics_scatter_cc` parameters from a namelist.

        implicit none

        class(numerics_scatter_cc_ext_t) :: self
        character(len=*) :: filename
        integer :: n_cond
        logical, optional :: verbose

        logical :: file_exists
        integer :: error

        integer :: n_principal_max = 10
            !! Maximum number of valence bands
        integer :: n_principal_min = 1
        integer :: n_cond_max = 0
            !! Maximum number of conduction bands

        integer :: n_q, n_q_theta, n_q_phi

        NAMELIST /numerics_s_cc_ext/ n_principal_max, &
                                 n_principal_min, &
                                 n_cond_max, &
                                 n_q, &
                                 n_q_theta, &
                                 n_q_phi

        if ( verbose ) then
            print*, 'Loading numerics parameters for scatter-cc_ext calculation...'
            print*
        end if

        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            open(100, file = trim(filename), iostat = error)
            read(100, nml = numerics_s_cc_ext, iostat = error)
            close(100)

            self%n_principal_max = min(n_principal_max, 10)
            self%n_principal_min = max(n_principal_min, 1)

            if ( n_cond_max == 0 ) then
                self%n_cond_max = n_cond
            else
                self%n_cond_max = min(n_cond, n_cond_max)
            end if

            self%n_q = n_q
            self%n_q_theta = n_q_theta
            self%n_q_phi = n_q_phi

            allocate(self%q_angular_mesh(self%n_q_theta*self%n_q_phi, 2))

            self%q_angular_mesh = generate_uniform_points_on_sphere(self%n_q_theta, self%n_q_phi)

            call self%print(verbose = verbose)

        else

            call print_error_message(&
                'Input file for numerics-scatter-cc_ext parameters : '//trim(filename)//' does NOT exist.', &
                verbose = verbose)
            stop

        end if

    end subroutine

    subroutine numerics_scatter_cc_ext_print(self, verbose)
        !! Prints `numerics_scatter_cc` components.

        implicit none

        class(numerics_scatter_cc_ext_t) :: self
        logical, optional :: verbose

        if ( verbose ) then

            call print_section_seperator()
            print*, '    ---------------------------'
            print*, '    Numerics - Scatter - cc_ext'
            print*, '    ---------------------------'
            print*
            print*, '        Maximum principal quantum number   : ', trim(adjustl(int_to_str(self%n_principal_max)))
            print*, '        Minimum principal quantum number   : ', trim(adjustl(int_to_str(self%n_principal_min)))
            print*
            print*, '        Maximum number of conduction bands : ', trim(adjustl(int_to_str(self%n_cond_max)))
            print*
            print*, '        n_q       = ', trim(adjustl(int_to_str(self%n_q)))
            print*, '        n_q_theta = ', trim(adjustl(int_to_str(self%n_q_theta)))
            print*, '        n_q_phi   = ', trim(adjustl(int_to_str(self%n_q_phi)))
            print*
            call print_section_seperator()
            print*

        end if

    end subroutine


end module
