module width_parameters_type
    !! Defines the `width_parameters` type. 

    use prec

    implicit none

    type width_parameters_t
        !! Electron lifetime/widths to use in calculations which require the electronic Green's
        !! functions (e.g. absorption calculations). width's are parameterized as
        !!
        !! $$\begin{align*}
        !!      \delta(\omega) = \text{min}\left( \delta_\text{max}, a + b \omega \right)
        !! \end{align*}$$
        !!
        !! and we require that \( |\omega - \omega_{ii'}| < \sigma \delta(\omega) \).
        integer :: n = 1
            !! Total number of width parameters.
            !!
            !! `n` = `n_a` x `n_b` x `n_m`
        integer :: n_a = 1
            !! Number of \( a \) parameters.
        integer :: n_b = 1
            !! Number of \( b \) parameters.
        integer :: n_m = 1
            !! Number of \( \delta_\text{max} \) parameters.
        real(dp) :: a_min
            !! Minimum \( a \) parameter.
            !!
            !! Units : eV
        real(dp) :: a_max
            !! Minimum \( a \) parameter.
            !!
            !! Units : eV
        real(dp) :: log_b_min
            !! Log10 of the minimum \( b \) parameter.
        real(dp) :: log_b_max
            !! Log10 of the maximum \( b \) parameter.
        real(dp) :: m_max
            !! Maximum \( \delta_\text{max} \) parameter.
            !!
            !! Units : eV
        real(dp) :: m_min
            !! Minimum \( \delta_\text{max} \) parameter.
            !!
            !! Units : eV
        real(dp) :: sigma = 1D10
            !! \( \sigma \) parameter, energy difference must be within \( \sigma \delta \) to be
            !! added to the rate.
        real(dp), allocatable :: info(:, :)
            !! Dim : [`n`, 3]
            !!
            !! Width parameters [a, b, delta_max].

        contains
            
            procedure :: load => load_width_parameters_nml
            procedure :: save => save_width_parameters
            procedure :: print => print_width_parameters
            procedure :: get_width

    end type

contains

    subroutine save_width_parameters(self, filename, verbose)
        !! Saves `width_parameters`.
        
        use hdf5
        use h5lt
        
        use info_messages
        
        implicit none

        class(width_parameters_t) :: self
        character(len=*) :: filename
        logical, optional :: verbose

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id
        logical :: file_exists
        integer(HSIZE_T) :: dims1(1) = [1]
        integer(HSIZE_T) :: dims2(2)

        integer :: error

        if ( verbose ) then
            print*, 'Saving width parameters...'
            print*
        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'widths', group_id, error)

            ! ! write data
            call h5ltmake_dataset_int_f(file_id, 'widths/n_a',&
                size(dims1), dims1,&
                self%n_a, error)
            call h5ltmake_dataset_int_f(file_id, 'widths/n_b',&
                size(dims1), dims1,&
                self%n_b, error)
            call h5ltmake_dataset_int_f(file_id, 'widths/n_m',&
                size(dims1), dims1,&
                self%n_m, error)
            call h5ltmake_dataset_int_f(file_id, 'widths/n',&
                size(dims1), dims1,&
                self%n, error)

            call h5ltmake_dataset_double_f(file_id, 'widths/a_min',&
                size(dims1), dims1,&
                self%a_min, error)
            call h5ltmake_dataset_double_f(file_id, 'widths/a_max',&
                size(dims1), dims1,&
                self%a_min, error)

            call h5ltmake_dataset_double_f(file_id, 'widths/log_b_min',&
                size(dims1), dims1,&
                self%log_b_min, error)
            call h5ltmake_dataset_double_f(file_id, 'widths/log_b_max',&
                size(dims1), dims1,&
                self%log_b_max, error)

            call h5ltmake_dataset_double_f(file_id, 'widths/m_min',&
                size(dims1), dims1,&
                self%m_min, error)
            call h5ltmake_dataset_double_f(file_id, 'widths/m_max',&
                size(dims1), dims1,&
                self%m_max, error)

            call h5ltmake_dataset_double_f(file_id, 'widths/sigma',&
                size(dims1), dims1,&
                self%sigma, error)

            dims2 = [self%n, 3]
            call h5ltmake_dataset_double_f(file_id, 'widths/info',&
                size(dims2), dims2,&
                self%info, error)

            call h5fclose_f(file_id, error)
            call h5close_f(error)

        else

            call print_error_message('Output file : '//trim(filename)//' does NOT exist.')
            stop

        end if

    end subroutine

    subroutine print_width_parameters(self, verbose)
        !! Prints `width_parameters` components.

        use info_messages

        implicit none

        class(width_parameters_t) :: self
        logical, optional :: verbose

        if ( verbose ) then

            call print_section_seperator()
            print*, '    ----------------'
            print*, '    Width Parameters'
            print*, '    ----------------'
            print*
            print*, '        Number of `a` parameters         : ', trim(adjustl(int_to_str(self%n_a))) 
            print*, '        Number of `b` parameters         : ', trim(adjustl(int_to_str(self%n_b))) 
            print*, '        Number of delta max parameters   : ', trim(adjustl(int_to_str(self%n_m))) 
            print*, '        Total number of width parameters : ', trim(adjustl(int_to_str(self%n))) 
            print*
            print*, '        Minimum `a` parameter : ', self%a_min, 'eV'
            print*, '        Maximum `a` parameter : ', self%a_max, 'eV'
            print*
            print*, '        Minimum `b` parameter : ', 10.0_dp**self%log_b_min
            print*, '        Maximum `b` parameter : ', 10.0_dp**self%log_b_max
            print*
            print*, '        Minimum delta max parameter : ', self%m_min, 'eV'
            print*, '        Maximum delta max parameter : ', self%m_max, 'eV'
            print*
            print*, '        sigma : ', self%sigma
            print*
            call print_section_seperator()
            print*
        end if

    end subroutine

    function get_width(self, id, omega) result( width )
        !! Returns \( \delta(\omega) \) for a given id.
        implicit none

        class(width_parameters_t) :: self
        integer :: id
        real(dp) :: omega

        real(dp) :: width

        width = min( self%info(id, 3), self%info(id, 1) + self%info(id, 2)*omega )

    end function

    subroutine load_width_parameters_nml(self, filename, verbose)
        !! Loads `width_parameters` parameters from a namelist.

        use info_messages
        use math_mod

        implicit none

        class(width_parameters_t) :: self
        character(len=*)  :: filename
        logical, optional :: verbose

        logical :: file_exists
        integer :: error

        integer :: n_a = 1
            !! Number of \( a \) parameters.
        integer :: n_b = 1
            !! Number of \( b \) parameters.
        integer :: n_m = 1
            !! Number of \( \delta_\text{max} \) parameters.
        real(dp) :: a_min
            !! Minimum \( a \) parameter.
            !!
            !! Units : eV
        real(dp) :: a_max
            !! Minimum \( a \) parameter.
            !!
            !! Units : eV
        real(dp) :: log_b_min
            !! Log10 of the minimum \( b \) parameter.
        real(dp) :: log_b_max
            !! Log10 of the maximum \( b \) parameter.
        real(dp) :: m_max
            !! Maximum \( \delta_\text{max} \) parameter.
            !!
            !! Units : eV
        real(dp) :: m_min
            !! Minimum \( \delta_\text{max} \) parameter.
            !!
            !! Units : eV
        real(dp) :: sigma = 1D10
            !! \( \sigma \) parameter, energy difference must be within \( \sigma \delta \) to be
            !! added to the rate.

        NAMELIST /widths/ n_a, &
                         n_b, &
                         n_m, &
                         a_min, &
                         a_max, &
                         log_b_min, &
                         log_b_max, &
                         m_max, &
                         m_min, &
                         sigma

        integer :: a, b, m, id

        real(dp), allocatable :: a_list(:)
        real(dp), allocatable :: log_b_list(:)
        real(dp), allocatable :: m_list(:)

        if ( verbose ) then
            print*, 'Loading the width parameters...'
            print*
        end if

        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            open(100, file = trim(filename), iostat = error)
            read(100, nml = widths, iostat = error)
            close(100)

            self%n_a = n_a
            self%n_b = n_b
            self%n_m = n_m

            self%a_min = a_min
            self%a_max = a_max
            self%log_b_min = log_b_min
            self%log_b_max = log_b_max
            self%m_max = m_max
            self%m_min = m_min
            self%sigma = sigma

            self%n = self%n_a*self%n_b*self%n_m
            allocate(self%info(self%n, 3))

            allocate(a_list(self%n_a))
            allocate(log_b_list(self%n_b))
            allocate(m_list(self%n_m))

            a_list = uniform_list(self%n_a, self%a_min, self%a_max)
            log_b_list = uniform_list(self%n_b, self%log_b_min, self%log_b_max)
            m_list = uniform_list(self%n_m, self%m_min, self%m_max)

            id = 0
            do a = 1, self%n_a
                do b = 1, self%n_b
                    do m = 1, self%n_m

                        id = id + 1

                        self%info(id, 1) = a_list(a) 
                        self%info(id, 2) = 10.0_dp**log_b_list(b)
                        self%info(id, 3) = m_list(m)

                    end do
                end do
            end do

            call self%print(verbose = verbose)

        else

            call print_error_message(&
                'Input file for width parameters : '//trim(filename)//' does NOT exist.', &
                verbose = verbose)
            stop

        end if

    end subroutine

end module
