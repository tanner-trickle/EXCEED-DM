module expt_type
    !! Defines `expt` type.

    use prec
    use constants

    implicit none

    type expt_t
        !! Experimental parameters, e.g. threshold energy, exposure, detector orientation.
        real(dp) :: E_threshold = 0.0_dp
            !! Energy threshold. Only events which deposit energy omega >= E_threshold will be considered.
            !!
            !! Units : eV
        real(dp) :: m_T_kg = 1.0_dp
            !! Mass of the target
            !!
            !! Units : kg
        real(dp) :: m_T
            !! Mass of the target
            !!
            !! Units : eV
        real(dp) :: exposure_yr = 1.0_dp
            !! Total exposure of the experiment
            !!
            !! Units : year
        real(dp) :: exposure
            !! Total exposure of the experiment
            !!
            !! Units : eV^(-1)
        integer :: n_time = 1
            !! Number of time of days to compute for
        real(dp), allocatable :: times(:)
            !! Dim : [ n_time ]
            !! 
            !! Times of the day.
            !!
            !! Units : days
        real(dp), allocatable :: vE_direction(:, :)
            !! Dim : [ n_time, 3 ]
            !!
            !! Unit vector of the Earth velocity in the galactic frame, assuming the DM wind
            !! is incoming on the Earth in the direction \( (0, 0, -1) \).
            !!
            !! See Fig. 1 of https://arxiv.org/abs/1909.09170 for an illustration.
            !!
            !! The standard setup assumes that the z-axis of the target is aligned anti-parallel
            !! to the DM wind at \( t = 0 \) (defined to be \( \hat{\mathbf{z}} \) ), the Earth then 
            !! rotates about a vector which is \( \theta_E = 42\deg \) off \( \hat{\mathbf{z}} \).
            !!
            !! Units : None
        real(dp) :: theta_E = 42.0_dp*(pi/180.0_dp)
            !! Angle between the Earth's rotation axis and the DM wind.
            !!
            !! Units : None

        contains

            procedure :: load => load_expt_nml
            procedure :: print => print_expt
            procedure :: save => save_expt

    end type

contains

    subroutine load_expt_nml(self, filename, verbose)
        !! Loads `expt` parameters from a namelist.

        use info_messages
        use units
        use math_mod

        implicit none

        class(expt_t) :: self
        character(len=*) :: filename
        logical, optional :: verbose

        logical :: file_exists
        integer :: error

        integer :: t

        ! namelist
        real(dp) :: E_threshold = 0.0_dp
            !! Energy threshold. Only events which deposit energy omega >= E_threshold will be considered.
            !!
            !! Units : eV
        real(dp) :: m_T_kg = 1.0_dp
            !! Mass of the target
            !!
            !! Units : kg
        real(dp) :: exposure_yr = 1.0_dp
            !! Total exposure of the experiment
            !!
            !! Units : year
        integer :: n_time = 1
            !! Number of time of days to compute for
        real(dp) :: theta_E = 42.0_dp*(pi/180.0_dp)

        NAMELIST /experiment/ E_threshold , &
                              m_T_kg      , &
                              exposure_yr , &
                              n_time, &
                              theta_E

        if ( verbose ) then
            print*, 'Loading the experimental parameters...'
            print*
        end if

        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            open(100, file = trim(filename), iostat = error)
            read(100, nml = experiment, iostat = error)
            rewind(100)

            self%E_threshold = E_threshold

            self%m_T_kg = m_T_kg
            self%m_T = kg_to_eV*m_T_kg

            self%exposure_yr = exposure_yr
            self%exposure = yr_to_inv_eV*exposure_yr

            self%n_time = n_time
            allocate(self%times(n_time))
            allocate(self%vE_direction(n_time, 3))
            self%times = uniform_list(n_time, 0.0_dp, 1.0_dp)

            self%theta_E = theta_E

            do t = 1, n_time

                self%vE_direction(t, 1) = sin(theta_E)*sin(2.0_dp*pi*self%times(t))
                self%vE_direction(t, 2) = cos(theta_E)*sin(theta_E)*(cos(2.0_dp*pi*self%times(t)) - 1.0_dp)
                self%vE_direction(t, 3) = ((sin(theta_E)**2)*cos(2.0_dp*pi*self%times(t)) &
                                    + cos(theta_E)**2)

            end do

            call self%print(verbose = verbose)

        else

            call print_error_message(&
                'Input file for experimental parameters : '//trim(filename)//' does NOT exist.', &
                verbose = verbose)
            stop

        end if

    end subroutine

    subroutine print_expt(self, verbose)
        !! Prints `expt` components.

        use info_messages

        class(expt_t) :: self
        logical, optional :: verbose

        if ( verbose ) then
            call print_section_seperator()
            print*, '    -----------------------'
            print*, '    Experimental Parameters'
            print*, '    -----------------------'
            print*
            print*, '        Energy threshold : ', self%E_threshold, 'eV'
            print*
            print*, '        Mass             : ', self%m_T_kg, 'kg'
            print*, '        Exposure Time    : ', self%exposure_yr, 'year'
            print*
            print*, '        Times of the day : ', trim(adjustl(int_to_str(self%n_time)))
            print* 
            call print_section_seperator()
            print*
        end if

    end subroutine

    subroutine save_expt(self, filename, verbose)
        !! Saves `expt`.

        use hdf5
        use h5lt

        use info_messages

        implicit none

        class(expt_t) :: self
        character(len=*) :: filename
        logical, optional :: verbose

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id
        logical :: file_exists
        integer(HSIZE_T) :: dims1(1) = [1]
        integer(HSIZE_T) :: dims2(2)
        integer :: error

        if ( verbose ) then
            print*, 'Saving experimental parameters...'
            print*
        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'experiment', group_id, error)

            ! ! write data
            call h5ltmake_dataset_int_f(file_id, 'experiment/n_time', size(dims1), dims1,&
                self%n_time, error)
            call h5ltmake_dataset_double_f(file_id, 'experiment/E_threshold', size(dims1), dims1,&
                self%E_threshold, error)
            dims1 = [self%n_time]
            call h5ltmake_dataset_double_f(file_id, 'experiment/times', size(dims1), dims1,&
                self%times, error)
            dims2 = [self%n_time, 3]
            call h5ltmake_dataset_double_f(file_id, 'experiment/vE_direction', size(dims2), dims2,&
                self%vE_direction, error)

            call h5ltmake_dataset_double_f(file_id, 'experiment/m_T', size(dims1), dims1,&
                self%m_T_kg, error)
            call h5ltmake_dataset_double_f(file_id, 'experiment/exposure', size(dims1), dims1,&
                self%exposure_yr, error)

            call h5fclose_f(file_id, error)
            call h5close_f(error)

        else

            call print_error_message('Output file : '//trim(filename)//' does NOT exist.')
            stop

        end if

    end subroutine

end module
