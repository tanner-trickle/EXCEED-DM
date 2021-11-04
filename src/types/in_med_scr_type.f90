module in_med_scr_type
    !! Defines the `in_med_scr` type.

    use prec
    use constants

    implicit none

    type in_med_scr_t
        !! Collection of parameters defining how the DM-electron scattering rate
        !! should be screened.

        character(len=64) :: type = ''
            !! Type of screening to use
            !! 
            !! Options: 'analytic', 'numeric'
            !!
            !! Will default to including no screening
        logical :: include_screen = .TRUE.
            !! Whether or not to include screening effects at all.
        real(dp) :: e0
            !! Static dielectric parameter for analytic screening
            !!
            !! \( \epsilon(0) \) in Eq 6 from [https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892)
            !!
            !! Units : None
        real(dp) :: q_tf
            !! Thomas Fermi momentum for analytic screening
            !!
            !! \( q_\text{TF} \) in Eq 6 from [https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892)
            !!
            !! Units : eV
        real(dp) :: omega_p
            !! Plasma frequency for analytic screening
            !!
            !! \( \omega_p \) in Eq 6 from [https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892)
            !!
            !! Units : eV
        real(dp) :: alpha
            !! \( \alpha \) shape parameter for analytic screening
            !!
            !! \( \alpha \) in Eq 6 from [https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892)
            !!
            !! Units : None
        integer :: n_E
            !! Number of \( \omega \) bins.
            !! 
            !! Numeric screening only.
        integer :: n_q
            !! Number of \( q \) bins.
            !! 
            !! Numeric screening only.
        integer :: n_q_theta
            !! Number of \( \theta_q \) bins.
            !! 
            !! Numeric screening only.
        integer :: n_q_phi
            !! Number of \( \phi_q \) bins
            !! 
            !! Numeric screening only.
        real(dp) :: E_width
            !! Width of the \( \omega \) bins.
            !!
            !! Numeric screening only.
            !!
            !! Units : eV
        real(dp) :: q_width
            !! Width of the \( q \) bins.
            !!
            !! Numeric screening only.
            !!
            !! Units : eV
        complex(dp), allocatable :: numeric_screen_mat(:, :, :, :)
            !! Dim : [ n_E, n_q, n_q_theta, n_q_phi ]
            !!
            !! Numerically computed screening matrix, binned in [ \( \omega, q, \theta_q, \phi_q \) ]
            !!
            !! Only used when `type = numeric`.
            !!
            !! Units : None
        complex(dp), allocatable :: numeric_anisotropic_screen_mat(:, :, :, :, :, :)
            !* Dim : [ n_E, n_q, n_q_theta, n_q_phi, 3, 3 ]
            !
            ! Numerically computed anisotropic screening matrix, binned in [ \( \omega, q, \theta_q, \phi_q \) ]. This will be used 
            ! if the input dielectric file is a 3 \( \times \) 3 matrix in each bin.
            !
            ! Only used when `type = numeric`.
            !
            ! Units : None
        logical :: anisotropic_screen = .FALSE.
            !!Flag to check whether the loaded dielectric is a matrix or scalar quantity. If matrix -> True, if scalar -> False.

        contains

            procedure :: screening
            procedure :: analytic_screening
            procedure :: numeric_screening
            procedure :: load_in_med_scr_nml
            procedure :: load => load_in_med_scr
            procedure :: save => save_in_med_scr
            procedure :: print => print_in_med_scr
            procedure :: load_computed_dielectric

    end type

contains

    subroutine save_in_med_scr(self, filename, verbose)
        !! Saves `in_med_scr`.
        use hdf5
        use h5lt

        use info_messages
        
        implicit none

        class(in_med_scr_t) :: self
        character(len=*) :: filename
        logical, optional :: verbose

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id
        logical :: file_exists
        integer(HSIZE_T) :: dims1(1) = [1]
        integer :: error

        if ( verbose ) then
            print*, 'Saving in medium screening parameters...'
            print*
        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'screening', group_id, error)

            ! write data

            if ( self%include_screen ) then
                call h5ltmake_dataset_string_f(file_id, &
                    'screening/type',                  &
                    self%type,                       &
                    error)

                if ( trim(self%type) == 'analytic' ) then

                    call h5ltmake_dataset_double_f(file_id, &
                        'screening/e0',                  &
                        size(dims1), dims1,                 &
                        self%e0,                       &
                        error)
                    call h5ltmake_dataset_double_f(file_id, &
                        'screening/q_tf',                  &
                        size(dims1), dims1,                 &
                        self%q_tf,                       &
                        error)
                    call h5ltmake_dataset_double_f(file_id, &
                        'screening/omega_p',                  &
                        size(dims1), dims1,                 &
                        self%omega_p,                       &
                        error)
                    call h5ltmake_dataset_double_f(file_id, &
                        'screening/alpha',                  &
                        size(dims1), dims1,                 &
                        self%alpha,                       &
                        error)
                    
                end if

                if ( trim(self%type) == 'numeric' ) then

                    call h5ltmake_dataset_double_f(file_id, &
                        'screening/q_width',                  &
                        size(dims1), dims1,                 &
                        self%q_width,                       &
                        error)
                    call h5ltmake_dataset_double_f(file_id, &
                        'screening/E_width',                  &
                        size(dims1), dims1,                 &
                        self%E_width,                       &
                        error)

                    call h5ltmake_dataset_int_f(file_id, &
                        'screening/n_E',                  &
                        size(dims1), dims1,                 &
                        self%n_E,                       &
                        error)
                    call h5ltmake_dataset_int_f(file_id, &
                        'screening/n_q',                  &
                        size(dims1), dims1,                 &
                        self%n_q,                       &
                        error)
                    call h5ltmake_dataset_int_f(file_id, &
                        'screening/n_q_theta',                  &
                        size(dims1), dims1,                 &
                        self%n_q_theta,                       &
                        error)
                    call h5ltmake_dataset_int_f(file_id, &
                        'screening/n_q_phi',                  &
                        size(dims1), dims1,                 &
                        self%n_q_phi,                       &
                        error)

                end if

            end if

            call h5fclose_f(file_id, error)
            call h5close_f(error)

        else

            call print_error_message(&
                'Output file : '//trim(filename)//' does NOT exist.', &
                verbose = verbose)
            stop

        end if

    end subroutine

    subroutine print_in_med_scr(self, verbose)
        !! Prints `in_med_scr` components.

        use info_messages

        implicit none

        class(in_med_scr_t) :: self
        logical, optional :: verbose

        if ( verbose ) then

            call print_section_seperator()
            print*, '    -------------------'
            print*, '    In-Medium Screening'
            print*, '    -------------------'
            print*
            print*, '        Screening type : ', trim(self%type)
            print*
            print*, '        Include screening?         : ', self%include_screen
            print*
            print*, '        Analytic screening parameters : '
            print*, '            e0      : ', self%e0
            print*, '            q_tf    : ', self%q_tf/1.0e3_dp, 'keV'
            print*, '            omega_p : ', self%omega_p, 'eV'
            print*, '            alpha   : ', self%alpha
            print*
            call print_section_seperator()
            print*
        end if

    end subroutine

    subroutine load_in_med_scr_nml(self, filename, verbose)
        !! Loads `in_med_scr` parameters from a namelist.
        
        use info_messages

        implicit none

        class(in_med_scr_t) :: self
        character(len=*) :: filename
        logical, optional :: verbose

        logical :: file_exists
        integer :: error

        character(len=64) :: type = ''
            !! Type of screening to use
            !! 
            !! Options: 'analytic', 'numeric'
            !!
            !! Will default to including no screening
        logical :: include_screen = .FALSE.
            !! Whether or not to include screening effects at all.
        real(dp) :: e0 = 1.0_dp
            !! Static dielectric parameter for analytic screening
            !!
            !! \( \epsilon(0) \) in Eq 6 from [https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892)
            !!
            !! Units : None
        real(dp) :: q_tf = 1.0e3_dp
            !! Thomas Fermi momentum for analytic screening
            !!
            !! \( q_\text{TF} \) in Eq 6 from [https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892)
            !!
            !! Units : eV
        real(dp) :: omega_p = 10.0_dp
            !! Plasma frequency for analytic screening
            !!
            !! \( \omega_p \) in Eq 6 from [https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892)
            !!
            !! Units : eV
        real(dp) :: alpha = 1.563_dp
            !! \( \alpha \) shape parameter for analytic screening
            !!
            !! \( \alpha \) in Eq 6 from [https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892)
            !!
            !! Units : None

        NAMELIST /in_med_scr/ include_screen, &
                                type, &
                              e0, &
                              q_tf, &
                              omega_p, &
                              alpha

        if ( verbose ) then
            print*, 'Loading in medium screening parameters from the namelist file...'
            print*
        end if

        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            open(100, file = trim(filename), iostat = error)
            read(100, nml = in_med_scr, iostat = error)
            close(100)

            self%include_screen = include_screen
            self%type           = type
            self%e0             = e0
            self%q_tf           = q_tf
            self%omega_p        = omega_p
            self%alpha          = alpha

            call self%print(verbose = verbose)

        else

            call print_error_message(&
                'Input file for dielectric binning parameters : '//trim(filename)//' does NOT exist.', &
                verbose = verbose)
            stop

        end if

    end subroutine

    subroutine load_computed_dielectric(self, filename, verbose)
        !! Load the pre-computed dielectric.

        use hdf5
        use h5lt

        use info_messages

        implicit none

        class(in_med_scr_t) :: self
        character(len=*) :: filename
        logical, optional :: verbose

        integer(HID_T) :: file_id

        logical :: file_exists

        integer(HSIZE_T) :: dims1(1) = [1]
        integer(HSIZE_T) :: dims4(4)
        integer(HSIZE_T) :: dims6(6)

        integer :: wfc_data_rank

        integer :: error

        real(dp), allocatable :: dielectric_buff(:, :, :, :)
        real(dp), allocatable :: anisotropic_dielectric_buff(:, :, :, :, :, :)

        if ( verbose ) then
            print*, 'Loading the dielectric from file...'
            print*
        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error)

            ! load binning parameters
            call h5ltread_dataset_int_f(file_id, 'bins_dielectric/n_E',&
                self%n_E, dims1, error)
            call h5ltread_dataset_int_f(file_id, 'bins_dielectric/n_q',&
                self%n_q, dims1, error)
            call h5ltread_dataset_int_f(file_id, 'bins_dielectric/n_q_theta',&
                self%n_q_theta, dims1, error)
            call h5ltread_dataset_int_f(file_id, 'bins_dielectric/n_q_phi',&
                self%n_q_phi, dims1, error)

            call h5ltread_dataset_double_f(file_id, 'bins_dielectric/E_width',&
                self%E_width, dims1, error)
            call h5ltread_dataset_double_f(file_id, 'bins_dielectric/q_width',&
                self%q_width, dims1, error)

            ! get the dimension of the dielectric dataset
            call h5ltget_dataset_ndims_f(file_id, 'dielectric/dielectric_r', wfc_data_rank, error)

            if ( wfc_data_rank == 4 ) then

                dims4 = [self%n_E, self%n_q, self%n_q_theta, self%n_q_phi]

                allocate(dielectric_buff(self%n_E, self%n_q, self%n_q_theta, self%n_q_phi))
                allocate(self%numeric_screen_mat(self%n_E, self%n_q, self%n_q_theta, self%n_q_phi))

                self%numeric_screen_mat = (0.0_dp, 0.0_dp)

                call h5ltread_dataset_double_f(file_id, 'dielectric/dielectric_r',&
                    dielectric_buff, dims4, error)

                self%numeric_screen_mat = self%numeric_screen_mat + dielectric_buff

                call h5ltread_dataset_double_f(file_id, 'dielectric/dielectric_c',&
                    dielectric_buff, dims4, error)

                self%numeric_screen_mat = self%numeric_screen_mat + ii*dielectric_buff

            else

                self%anisotropic_screen = .TRUE.

                ! load the anisotropic dielectric

                dims6 = [self%n_E, self%n_q, self%n_q_theta, self%n_q_phi, 3, 3]

                allocate(anisotropic_dielectric_buff(self%n_E, self%n_q, self%n_q_theta, self%n_q_phi, 3, 3))
                allocate(self%numeric_anisotropic_screen_mat(self%n_E, self%n_q, self%n_q_theta, self%n_q_phi, 3, 3))

                self%numeric_anisotropic_screen_mat = (0.0_dp, 0.0_dp)

                call h5ltread_dataset_double_f(file_id, 'dielectric/dielectric_r',&
                    anisotropic_dielectric_buff, dims6, error)

                self%numeric_anisotropic_screen_mat = self%numeric_anisotropic_screen_mat + anisotropic_dielectric_buff

                call h5ltread_dataset_double_f(file_id, 'dielectric/dielectric_c',&
                    anisotropic_dielectric_buff, dims6, error)

                self%numeric_anisotropic_screen_mat = self%numeric_anisotropic_screen_mat + ii*anisotropic_dielectric_buff

            end if

            call h5fclose_f(file_id, error)
            call h5close_f(error)

        else

            call print_error_message('Input file for dielectric : '//trim(filename)//' does NOT exist.', &
                verbose = verbose)
            stop

        end if

        if ( verbose ) then
            print*, 'Done loading the dielectric from file!'
            print*
        end if

    end subroutine

    subroutine load_in_med_scr(self, proc_id, root_process, n_proc, &
            io_files, main_control, target_mat, verbose)
        !! Loads the in medium screening parameters. If the user specifies
        !! that the screening should be computed numerically then `exdm_dielectric` will
        !! be run.

        use mpi
        use hdf5
        use h5lt

        use io_input
        use control_input

        use material_type

        use exdm_dielectric

        implicit none

        class(in_med_scr_t) :: self
        integer :: proc_id, root_process, n_proc
        type(io_files_t) :: io_files
        type(control_t) :: main_control
        type(material_t) :: target_mat
        logical, optional :: verbose

        logical :: file_exists
        integer(HID_T) :: file_id
        integer :: error

        integer :: err

        integer :: n_E, n_q, n_q_theta, n_q_phi
        real(dp) :: E_width, q_width

        complex(dp), allocatable :: screen_mat(:, :, :, :)

        call self%load_in_med_scr_nml(io_files%nml_input_filename, verbose = verbose)

        inquire(file = io_files%dielectric_filename, exist = file_exists)

        ! If the user wants to use the numerically computed screening factor AND include_screen = .TRUE.
        ! AND the file doesn't already exist, compute dielectric.

        if ( ( trim(self%type) == 'numeric' ) .and. ( self%include_screen ) ) then

            if ( .not. file_exists ) then

                if ( verbose ) then

                    print*, 'Dielectric file does NOT exist. Computing from scratch...'
                    print*

                end if

                if ( proc_id == root_process ) then

                    ! create the file
                    call h5open_f(error)

                    call h5fcreate_f(io_files%dielectric_filename, H5F_ACC_TRUNC_F, file_id, error)

                    call h5fclose_f(file_id, error)
                    call h5close_f(error)

                end if

                ! compute the dielectric
                call run_exdm_dielectric(proc_id, root_process, n_proc, io_files, & 
                    main_control, target_mat, verbose = verbose)

            end if

            ! load the dielectric, bcast relevant screen_mat data to all processors
            if ( proc_id == root_process ) then

                call self%load_computed_dielectric(io_files%dielectric_filename, verbose = verbose)

                n_E = self%n_E
                n_q = self%n_q
                n_q_theta = self%n_q_theta
                n_q_phi = self%n_q_phi

                E_width = self%E_width
                q_width = self%q_width

            end if

            call MPI_Bcast(n_E, 1, MPI_INTEGER, root_process,&
              MPI_COMM_WORLD, err)
            call MPI_Bcast(n_q, 1, MPI_INTEGER, root_process,&
              MPI_COMM_WORLD, err)
            call MPI_Bcast(n_q_theta, 1, MPI_INTEGER, root_process,&
              MPI_COMM_WORLD, err)
            call MPI_Bcast(n_q_phi, 1, MPI_INTEGER, root_process,&
              MPI_COMM_WORLD, err)

            call MPI_Bcast(E_width, 1, MPI_DOUBLE, root_process,&
              MPI_COMM_WORLD, err)
            call MPI_Bcast(q_width, 1, MPI_DOUBLE, root_process,&
              MPI_COMM_WORLD, err)

            if ( proc_id /= root_process ) then

                self%n_E = n_E
                self%n_q = n_q
                self%n_q_theta = n_q_theta
                self%n_q_phi = n_q_phi

                self%E_width = E_width
                self%q_width = q_width

                !! allocate screen_mat
                allocate(self%numeric_screen_mat(self%n_E, &
                    self%n_q, self%n_q_theta, self%n_q_phi))

            end if

            allocate(screen_mat(n_E, n_q, n_q_theta, n_q_phi))

            if ( proc_id == root_process ) then

                screen_mat = self%numeric_screen_mat

            end if

            call MPI_Bcast(screen_mat, size(screen_mat), MPI_DOUBLE_COMPLEX, &
                root_process, MPI_COMM_WORLD, err)

            if ( proc_id /= root_process ) then

                self%numeric_screen_mat = screen_mat

            end if

        end if 

    end subroutine

    function numeric_screening(self, q_vec, omega) result( scr )
        !! Absolute value of the numerically computed dielectric.

        use math_mod

        implicit none

        class(in_med_scr_t) :: self
        real(dp) :: q_vec(3)
        real(dp) :: omega

        real(dp) :: q_mag
        real(dp) :: q_hat(3)

        real(dp) :: q_theta, q_phi

        integer :: w_bin, q_bin, q_theta_bin, q_phi_bin

        real(dp) :: scr

        q_mag = norm2(q_vec)

        scr = 1.0_dp

        if ( q_mag >= 1.0e-8_dp ) then

            q_hat = q_vec/q_mag

            q_theta = get_theta(q_hat)
            q_phi = get_phi(q_hat)

            w_bin = 1 + floor(omega/self%E_width)
            q_bin = 1 + floor(q_mag/self%q_width)
            q_theta_bin = 1 + floor(q_theta/(pi/max(1.0_dp, 1.0_dp*self%n_q_theta)))
            q_phi_bin = 1 + floor(q_phi/(2.0_dp*pi/max(1.0_dp, 1.0_dp*self%n_q_phi)))

            if ( ( w_bin <= self%n_E ) .and. &
                 ( q_bin <= self%n_q ) .and. & 
                 ( q_theta_bin <= self%n_q_theta ) .and. &
                 ( q_phi_bin <= self%n_q_phi ) ) then

                if ( self%anisotropic_screen ) then

                    scr = abs(& 
                                dot_product(q_hat,&
                                        matmul(self%numeric_anisotropic_screen_mat(&
                                                    w_bin,& 
                                                    q_bin,& 
                                                    q_theta_bin,& 
                                                    q_phi_bin, :, :), q_hat&
                                                )&
                                        )&
                                )
                else

                    scr = abs( self%numeric_screen_mat(w_bin, &
                        q_bin, q_theta_bin, q_phi_bin) )

                end if


            end if

        end if

    end function

    function analytic_screening(self, q_vec, omega) result( scr )
        !! Analytic form of the dielectric function.
        !!
        !! Eq 6 in [https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892)

        implicit none

        class(in_med_scr_t) :: self
        real(dp) :: q_vec(3)
        real(dp) :: omega

        real(dp) :: q_mag
        real(dp) :: scr

        q_mag = norm2(q_vec)

        scr = 1.0_dp + ( (self%e0 - 1.0_dp)**(-1) + &
            self%alpha*(q_mag/self%q_tf)**2 + &
            q_mag**4/(4*m_elec**2*self%omega_p**2) - &
            (omega/self%omega_p)**2)**(-1) 

    end function

    function screening(self, q_vec, omega) result(scr)
        !! Screening factor in scattering rate calculations.
        !!
        !! \( R \sim 1/text{scr}^2 \)

        implicit none

        class(in_med_scr_t) :: self
        real(dp) :: q_vec(3)
        real(dp) :: omega

        real(dp) :: q_mag
        real(dp) :: q_hat(3)
        real(dp) :: q_theta, q_phi

        integer :: w_bin_num, q_bin_num, q_theta_bin_num, q_phi_bin_num

        real(dp) :: scr

        scr = 1.0_dp

        if ( self%include_screen ) then

            if ( trim(self%type) == 'analytic' ) then
                scr = self%analytic_screening(q_vec, omega)
            end if

            if ( trim(self%type) == 'numeric' ) then
                scr = self%numeric_screening(q_vec, omega)
            end if

        end if

    end function

end module
