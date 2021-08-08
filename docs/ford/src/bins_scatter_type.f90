module bins_scatter_type
    !! Defines the `bins_scatter` type.

    use prec

    implicit none

    type bins_scatter_t
        !! Collection of parameters which specify how the DM-electron scattering rate
        !! is binned in \( q, \omega \).
        integer :: n_q = 1
            !! Number of bins in momentum space
        integer :: n_E = 1
            !! Number of bins in omega space
        real(dp) :: q_width = 1.0e3_dp
            !! Width of the q bins
            !!
            !! Units : eV
        real(dp) :: E_width = 0.1_dp
            !! Width of the omega bins.
            !!
            !! To get the event rate per ionization threshold, \( Q \), set this to
            !! set this to the \( \epsilon \) parameter.
            !!
            !! Units : eV
        contains

            procedure :: print => print_bins
            procedure :: load => load_bins_nml
            procedure :: save => save_bins

    end type

contains

    subroutine print_bins(self, verbose)
        !! Print `bins_scatter` components.

        use info_messages

        implicit none

        class(bins_scatter_t) :: self
        logical, optional :: verbose

        if ( verbose ) then
            call print_section_seperator()
            print*, '    --------------'
            print*, '    Bins - Scatter'
            print*, '    --------------'
            print*
            print*, '        Number of q bins : ', trim(adjustl(int_to_str(self%n_q)))
            print*, '        q bin width      : ', self%q_width/1.0e3_dp, 'keV'
            print*
            print*, '        Number of E bins : ', trim(adjustl(int_to_str(self%n_E)))
            print*, '        E bin width      : ', self%E_width, 'eV'
            print*
            call print_section_seperator()
            print*
        end if

    end subroutine

    subroutine load_bins_nml(self, filename, verbose)
        !! Loads `bins_scatter` parameters from a namelist.

        use info_messages

        implicit none

        class(bins_scatter_t) :: self
        character(len=*) :: filename
        logical, optional :: verbose

        logical :: file_exists
        integer :: error

        integer :: n_q = 1 
        integer :: n_E = 1
        real(dp) :: q_width = 1.0e3_dp
        real(dp) :: E_width = 0.1_dp

        NAMELIST /bins_scatter/ n_q    , &
                           n_E    , &
                           q_width , &
                           E_width

        if ( verbose ) then
            print*, 'Loading binning parameters...'
            print*
        end if

        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            open(100, file = trim(filename), iostat = error)
            read(100, nml = bins_scatter, iostat = error)
            close(100)

            self%n_q     = n_q
            self%n_E     = n_E
            self%q_width = q_width
            self%E_width = E_width

            call self%print(verbose = verbose)

        else

            call print_error_message(&
                'Input file for binning parameters : '//trim(filename)//' does NOT exist.', &
                verbose = verbose)
            stop

        end if

    end subroutine

    subroutine save_bins(self, filename, verbose)
        !! Saves `bins_scatter`.
        
        use hdf5
        use h5lt

        use info_messages
        
        implicit none

        class(bins_scatter_t) :: self
        character(len=*) :: filename
        logical, optional :: verbose

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id
        logical :: file_exists
        integer(HSIZE_T) :: dims1(1) = [1]
        integer :: error

        if ( verbose ) then
            print*, 'Saving binning parameters...'
            print*
        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'bins_scatter', group_id, error)

            ! write data
            call h5ltmake_dataset_double_f(file_id, &
                'bins_scatter/q_width',                  &
                size(dims1), dims1,                 &
                self%q_width,                       &
                error)
            call h5ltmake_dataset_double_f(file_id, &
                'bins_scatter/E_width',                  &
                size(dims1), dims1,                 &
                self%E_width,                       &
                error)
            call h5ltmake_dataset_int_f(file_id, &
                'bins_scatter/n_q',                   &
                size(dims1), dims1,              &
                self%n_q,                        &
                error)
            call h5ltmake_dataset_int_f(file_id, &
                'bins_scatter/n_E',                   &
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
