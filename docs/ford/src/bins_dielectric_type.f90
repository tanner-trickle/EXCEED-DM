module bins_dielectric_type
    !! Defines the `bins_dielectric` type.

    use prec

    implicit none

    type bins_dielectric_t
        !! Collection of parameters which specify how the dielectric function
        !! is binned in \( q, \theta_q, \phi_q, \omega \).
        integer :: n_q = 1
            !! Number of bins for \( q \).
        integer :: n_q_theta = 1
            !! Number of bins for the \( \theta \) coordinate of \( \mathbf{q} \)
        integer :: n_q_phi = 1
            !! Number of bins for the \( \phi \) coordinate of \( \mathbf{q} \)
        integer :: n_E = 1
            !! Number of bins in omega space
        real(dp) :: q_width = 1.0e3_dp
            !! Width of the q bins
            !!
            !! Units : eV
        real(dp) :: E_width = 0.1_dp
            !! Width of the omega bins.
            !!
            !! Units : eV
        contains

            procedure :: print => print_bins_dielectric
            procedure :: load => load_bins_dielectric_nml
            procedure :: save => save_bins_dielectric

    end type

contains

    subroutine print_bins_dielectric(self, verbose)
        !! Print `bins_dielectric` components.

        use info_messages

        implicit none

        class(bins_dielectric_t) :: self
        logical, optional :: verbose

        if ( verbose ) then
            call print_section_seperator()
            print*, '    -----------------'
            print*, '    Bins - Dielectric'
            print*, '    -----------------'
            print*
            print*, '        Number of q bins : ', trim(adjustl(int_to_str(self%n_q)))
            print*, '        q bin width      : ', self%q_width/1.0e3_dp, 'keV'
            print*
            print*, '        Number of q theta bins : ', trim(adjustl(int_to_str(self%n_q_theta)))
            print*, '        Number of q phi bins   : ', trim(adjustl(int_to_str(self%n_q_phi)))
            print*
            print*, '        Number of E bins : ', trim(adjustl(int_to_str(self%n_E)))
            print*, '        E bin width      : ', self%E_width, 'eV'
            print*
            call print_section_seperator()
            print*
        end if

    end subroutine

    subroutine load_bins_dielectric_nml(self, filename, verbose)
        !! Loads `bins_dielectric` parameters from a namelist.

        use info_messages

        implicit none

        class(bins_dielectric_t) :: self
        character(len=*) :: filename
        logical, optional :: verbose

        logical :: file_exists
        integer :: error

        integer :: n_q = 1 
        integer :: n_q_theta = 1
        integer :: n_q_phi = 1
        integer :: n_E = 1
        real(dp) :: q_width = 1.0e3_dp
        real(dp) :: E_width = 0.1_dp

        NAMELIST /bins_dielectric/ n_q    , &
                                   n_q_theta, &
                                   n_q_phi, &
                                   n_E, &
                                   q_width, &
                                   E_width

        if ( verbose ) then
            print*, 'Loading dielectric binning parameters...'
            print*
        end if

        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            open(100, file = trim(filename), iostat = error)
            read(100, nml = bins_dielectric, iostat = error)
            close(100)

            self%n_q     = n_q
            self%n_q_theta     = n_q_theta
            self%n_q_phi     = n_q_phi
            self%n_E     = n_E
            self%q_width = q_width
            self%E_width = E_width

            call self%print(verbose = verbose)

        else

            call print_error_message(&
                'Input file for dielectric binning parameters : '//trim(filename)//' does NOT exist.', &
                verbose = verbose)
            stop

        end if

    end subroutine

    subroutine save_bins_dielectric(self, filename, verbose)
        !! Saves `bins_dielectric`.
        
        use hdf5
        use h5lt

        use info_messages
        
        implicit none

        class(bins_dielectric_t) :: self
        character(len=*) :: filename
        logical, optional :: verbose

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id
        logical :: file_exists
        integer(HSIZE_T) :: dims1(1) = [1]
        integer :: error

        if ( verbose ) then
            print*, 'Saving dielectric binning parameters...'
            print*
        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'bins_dielectric', group_id, error)

            ! write data
            call h5ltmake_dataset_double_f(file_id, &
                'bins_dielectric/q_width',                  &
                size(dims1), dims1,                 &
                self%q_width,                       &
                error)
            call h5ltmake_dataset_double_f(file_id, &
                'bins_dielectric/E_width',                  &
                size(dims1), dims1,                 &
                self%E_width,                       &
                error)
            call h5ltmake_dataset_int_f(file_id, &
                'bins_dielectric/n_q',                   &
                size(dims1), dims1,              &
                self%n_q,                        &
                error)
            call h5ltmake_dataset_int_f(file_id, &
                'bins_dielectric/n_q_theta',                   &
                size(dims1), dims1,              &
                self%n_q_theta,                        &
                error)
            call h5ltmake_dataset_int_f(file_id, &
                'bins_dielectric/n_q_phi',                   &
                size(dims1), dims1,              &
                self%n_q_phi,                        &
                error)
            call h5ltmake_dataset_int_f(file_id, &
                'bins_dielectric/n_E',                   &
                size(dims1), dims1,              &
                self%n_E,                        &
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

end module
