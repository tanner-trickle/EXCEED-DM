module dm_model_type
    !! Defines the `dm_model` type.

    use hdf5
    use h5lt

    use info_messages

    use prec
    use constants
    use units
    use math_mod

    implicit none

    type dm_model_t
        !! Parameters defining the dark matter (DM) model, i.e. the DM particle and mediator.

        integer :: n_mX = 1
            !! Number of masses
        real(dp), allocatable :: mX(:)
            !! Dim : [ n_mX ]
            !!
            !! List of masses
            !!
            !! Units : eV
        integer :: n_med_FF = 1
            !! Number of mediator form factors
        real(dp), allocatable :: med_FF(:)
            !! Dim : [ n_med_FF ]
            !!
            !! List of mediator form factors coefficients, \( \beta \)
            !!
            !! $$\begin{align}
            !!      \mathcal{F}_\text{med} = \left( \frac{ \alpha m_e }{ q } \right)^\beta \nonumber
            !! \end{align}$$
            !!
            !! Examples:
            !! <ul>
            !!     <li>Light mediator : \( \beta = 2 \) </li>
            !!     <li>Heavy mediator : \( \beta = 0 \) </li>
            !! </ul>
            !!
            !! Units : None
        real(dp) :: rhoX_GeV_per_cm3 = 0.4_dp
            !! Dark matter density
            !!
            !! Units : GeV/cm^3
        real(dp) :: rhoX
            !! Dark matter density
            !!
            !! Units : eV^4
        real(dp) :: v0_km_per_sec = 230.0_dp
            !! DM velocity distribution parameter, \( v_0 \)
            !!
            !! $$\begin{align}
            !!      f_\chi(\mathbf{v}) = \frac{1}{N_0} 
            !!          e^{-( \mathbf{v} + \mathbf{v}_e)^2/v_0^2}
            !!          \Theta( v_\text{esc} - |\mathbf{v} + \mathbf{v}_e|) \nonumber
            !! \end{align}$$
            !!
            !! Units : km/s
        real(dp) :: vE_km_per_sec = 240.0_dp
            !! DM velocity distribution parameter, \( |\mathbf{v}_e| \)
            !!
            !! $$\begin{align}
            !!      f_\chi(\mathbf{v}) = \frac{1}{N_0} 
            !!          e^{-( \mathbf{v} + \mathbf{v}_e)^2/v_0^2}
            !!          \Theta( v_\text{esc} - |\mathbf{v} + \mathbf{v}_e|) \nonumber
            !! \end{align}$$
            !!
            !! Units : km/s
        real(dp) :: vEsc_km_per_sec = 600.0_dp
            !! DM velocity distribution parameter, \( v_\text{esc} \)
            !!
            !! $$\begin{align}
            !!      f_\chi(\mathbf{v}) = \frac{1}{N_0} 
            !!          e^{-( \mathbf{v} + \mathbf{v}_e)^2/v_0^2}
            !!          \Theta( v_\text{esc} - |\mathbf{v} + \mathbf{v}_e|) \nonumber
            !! \end{align}$$
            !!
            !! Units : km/s
        real(dp) :: v0
            !! DM velocity distribution parameter, \( v_0 \)
            !!
            !! $$\begin{align}
            !!      f_\chi(\mathbf{v}) = \frac{1}{N_0} 
            !!          e^{-( \mathbf{v} + \mathbf{v}_e)^2/v_0^2}
            !!          \Theta( v_\text{esc} - |\mathbf{v} + \mathbf{v}_e|) \nonumber
            !! \end{align}$$
            !!
            !! Units : None
        real(dp) :: vE
            !! DM velocity distribution parameter, \( |\mathbf{v}_e| \)
            !!
            !! $$\begin{align}
            !!      f_\chi(\mathbf{v}) = \frac{1}{N_0} 
            !!          e^{-( \mathbf{v} + \mathbf{v}_e)^2/v_0^2}
            !!          \Theta( v_\text{esc} - |\mathbf{v} + \mathbf{v}_e|) \nonumber
            !! \end{align}$$
            !!
            !! Units : None
        real(dp) :: vEsc
            !! DM velocity distribution parameter, \( v_\text{esc} \)
            !!
            !! $$\begin{align}
            !!      f_\chi(\mathbf{v}) = \frac{1}{N_0} 
            !!          e^{-( \mathbf{v} + \mathbf{v}_e)^2/v_0^2}
            !!          \Theta( v_\text{esc} - |\mathbf{v} + \mathbf{v}_e|) \nonumber
            !! \end{align}$$
            !!
            !! Units : None
        real(dp) :: g_func_N0
            !! Parameter of
            !! $$\begin{align*}
            !!  g( \mathbf{q}, \omega ) & = 2 \pi \int d^3\mathbf{v} f_\chi(\mathbf{v}) \delta(\omega - \omega_\mathbf{q} ) \\
            !!  & = \frac{2 \pi^2 v_0^2}{N_0 q} \left( e^{-v_-^2/v_0^2} - e^{-v_\text{esc}^2/v_0^2} \right) \\
            !!  \omega_\mathbf{q} & = \mathbf{q} \cdot \mathbf{v} - \frac{q^2}{2 m_\chi}
            !! \end{align*}$$
            !!
            !! `g_func_N0` \( = N_0 = \pi v_0^3 
            !! \left( \sqrt{\pi} \text{erf}(v_\text{esc}/v_0) - 2\frac{v_\text{esc}}{v_0}
            !! e^{-(v_\text{esc}/v_0)^2} \right) \)
            !!
            !! Units : None
        real(dp) :: g_func_c1
            !! Parameter of 
            !! $$\begin{align*}
            !!  g( \mathbf{q}, \omega ) & = 2 \pi \int d^3\mathbf{v} f_\chi(\mathbf{v}) \delta(\omega - \omega_\mathbf{q} ) \\
            !!  & = \frac{2 \pi^2 v_0^2}{N_0 q} \left( e^{-v_-^2/v_0^2} - e^{-v_\text{esc}^2/v_0^2} \right) \\
            !!  \omega_\mathbf{q} & = \mathbf{q} \cdot \mathbf{v} - \frac{q^2}{2 m_\chi}
            !! \end{align*}$$
            !!
            !! `g_func_c1` \( = \frac{2 \pi^2 v_0^2}{N_0}  \)
            !!
            !! Units : None
        real(dp) :: g_func_c2
            !! Parameter of 
            !! $$\begin{align*}
            !!  g( \mathbf{q}, \omega ) & = 2 \pi \int d^3\mathbf{v} f_\chi(\mathbf{v}) \delta(\omega - \omega_\mathbf{q} ) \\
            !!  & = \frac{2 \pi^2 v_0^2}{N_0 q} \left( e^{-v_-^2/v_0^2} - e^{-v_\text{esc}^2/v_0^2} \right) \\
            !!  \omega_\mathbf{q} & = \mathbf{q} \cdot \mathbf{v} - \frac{q^2}{2 m_\chi}
            !! \end{align*}$$
            !!
            !! `g_func_c2` \( = e^{-v_\text{esc}^2/v_0^2} \)
            !!
            !! Units : None
        real(dp) :: vX_max
            !! Maximum speed of the incoming dark matter
            !!
            !! `vX_max` = `vE` + `vEsc`
            !!
            !! Units : None
        integer :: tff_id(2) = [1, 1]
            !! Transition form factor ID. Specifies which matrix element to compute, which
            !! will depend on the scattering potential of the DM model
        character(len=64) :: particle_type = 'fermion'
            !! Dark matter particle type

        contains 

            procedure :: load  => load_dm_model_nml
            procedure :: save  => save_dm_model
            procedure :: print => print_dm_model

    end type

contains

    subroutine print_dm_model(self, verbose)
        !! Prints `dm_model` components.

        implicit none

        class(dm_model_t) :: self
        logical, optional :: verbose

        if ( verbose ) then
            call print_section_seperator()
            print*, '    -----------------'
            print*, '    Dark Matter Model'
            print*, '    -----------------'
            print*
            print*, '        Particle type : ', trim(self%particle_type)
            print*
            print*, '        Density     : ', self%rhoX_GeV_per_cm3, 'GeV/cm^3'
            print*, '        Masses (eV) : ', self%mX
            print* 
            print*, '        Mediator Form Factors (-d log F_DM / d log q) : ', self%med_FF
            print*, '        Transition Form Factor ID                     : ', self%tff_id
            print* 
            print*, '        Velocity Distribution Parameters : '
            print*, '            v0   = ', self%v0_km_per_sec, 'km/sec'
            print*, '            vE   = ', self%vE_km_per_sec, 'km/sec'
            print*, '            vEsc = ', self%vEsc_km_per_sec, 'km/sec'
            print* 
            call print_section_seperator()
            print*
        end if

    end subroutine

    subroutine load_dm_model_nml(self, filename, verbose)
        !! Loads `dm_model` parameters from a namelist.

        implicit none

        class(dm_model_t) :: self
        character(len=*)  :: filename
        logical, optional :: verbose

        logical :: file_exists
        integer :: error

        ! namelist
        integer  :: n_mX = 1
        real(dp) :: log_mX_min = 9.0_dp
        real(dp) :: log_mX_max = 9.0_dp
        
        integer :: n_extra_mX = 0
            !! Optional
            !!
            !! User can specify extra mass points to add in addition to 
            !! the log-uniform ones chosen

        integer :: n_med_FF = 1
        real(dp) :: med_FF_min = 0.0_dp
        real(dp) :: med_FF_max = 0.0_dp

        real(dp) :: rhoX_GeV_per_cm3 = 0.4_dp

        real(dp) :: v0_km_per_sec   = 230.0_dp
        real(dp) :: vE_km_per_sec   = 240.0_dp
        real(dp) :: vEsc_km_per_sec = 600.0_dp

        integer :: tff_id(2) = [1, 1]

        character(len=64) :: particle_type = 'fermion'

        real(dp), allocatable :: mX_2(:)

        NAMELIST /dm_model/ n_mX             , &
                            log_mX_min       , &
                            log_mX_max       , &
                            n_extra_mX       , &
                            n_med_FF         , &
                            med_FF_min       , &
                            med_FF_max       , &
                            rhoX_GeV_per_cm3 , &
                            v0_km_per_sec    , &
                            vE_km_per_sec    , &
                            vEsc_km_per_sec  , &
                            tff_id           , &
                            particle_type

        NAMELIST /extra_mX/ mX_2

        if ( verbose ) then
            print*, 'Loading the dark matter model...'
            print*
        end if

        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            open(100, file = trim(filename), iostat = error)
            read(100, nml = dm_model, iostat = error)
            rewind(100)

            self%n_med_FF = n_med_FF
            allocate(self%med_FF(n_med_FF))
            self%med_FF = uniform_list(n_med_FF, med_FF_min, med_FF_max)
            
            self%rhoX_GeV_per_cm3 = rhoX_GeV_per_cm3
            self%rhoX = inv_cm_to_eV**3*1.0e9_dp*rhoX_GeV_per_cm3

            self%v0_km_per_sec = v0_km_per_sec
            self%vE_km_per_sec = vE_km_per_sec
            self%vEsc_km_per_sec = vEsc_km_per_sec
            
            self%v0        = km_per_sec_to_none*v0_km_per_sec
            self%vE        = km_per_sec_to_none*vE_km_per_sec
            self%vEsc      = km_per_sec_to_none*vEsc_km_per_sec
            self%vX_max    = self%vE + self%vEsc
            self%g_func_N0 = (pi*self%v0**3)*&
                (sqrt(pi)*erf(self%vEsc/self%v0) - 2*(self%vEsc/self%v0)*exp(-(self%vEsc/self%v0)**2))
            self%g_func_c1 = (2*pi**2*self%v0**2/self%g_func_N0)
            self%g_func_c2 = exp(-(self%vEsc/self%v0)**2)

            self%particle_type = particle_type

            if ( n_extra_mX /= 0 ) then

                allocate(mX_2(n_extra_mX))

                mX_2 = 1.0e9_dp

                read(100, nml = extra_mX, iostat = error)
                rewind(100)

            end if

            close(100)

            self%n_mX = n_mX + n_extra_mX
            allocate(self%mX(self%n_mX))

            self%mX(:n_mX) = 10.0_dp**uniform_list(n_mX, log_mX_min, log_mX_max)
            if ( n_extra_mX /= 0 ) then
                self%mX(n_mX + 1:) = mX_2
            end if

            call self%print(verbose = verbose)

        else

            call print_error_message(&
                'Input file for dark matter model : '//trim(filename)//' does NOT exist.', &
                verbose = verbose)
            stop

        end if

    end subroutine

    subroutine save_dm_model(self, filename, verbose)
        !! Saves `dm_model`.
        implicit none

        class(dm_model_t) :: self
        character(len=*) :: filename
        logical, optional :: verbose

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id
        logical :: file_exists
        integer(HSIZE_T) :: dims1(1) = [1]

        integer :: error

        if ( verbose ) then
            print*, 'Saving dark matter model...'
            print*
        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'dm_model', group_id, error)

            ! ! write data
            call h5ltmake_dataset_string_f(file_id, 'dm_model/particle_type', &
                trim(self%particle_type), error)
            call h5ltmake_dataset_int_f(file_id, 'dm_model/n_mX', size(dims1), dims1,&
                self%n_mX, error)
            call h5ltmake_dataset_int_f(file_id, 'dm_model/n_med_FF', size(dims1), dims1,&
                self%n_med_FF, error)
            dims1 = [self%n_mX]
            call h5ltmake_dataset_double_f(file_id, 'dm_model/mX', size(dims1), dims1,&
                self%mX, error)
            dims1 = [self%n_med_FF]
            call h5ltmake_dataset_double_f(file_id, 'dm_model/med_FF', size(dims1), dims1,&
                self%med_FF, error)
            call h5ltmake_dataset_double_f(file_id, 'dm_model/v0', size(dims1), dims1,&
                self%v0_km_per_sec, error)
            call h5ltmake_dataset_double_f(file_id, 'dm_model/vEsc', size(dims1), dims1,&
                self%vEsc_km_per_sec, error)
            call h5ltmake_dataset_double_f(file_id, 'dm_model/vE', size(dims1), dims1,&
                self%vE_km_per_sec, error)

            call h5fclose_f(file_id, error)
            call h5close_f(error)

        else

            call print_error_message('Output file : '//trim(filename)//' does NOT exist.')
            stop

        end if

    end subroutine

end module
