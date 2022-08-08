module astroph_model_type

    use prec_util, only: dp

    implicit none

    type :: astroph_model_t

        real(dp) :: v_0_km_per_sec
        real(dp) :: v_0

        real(dp), allocatable :: v_e_km_per_sec(:, :)
        real(dp), allocatable :: v_e_list(:, :)

        real(dp) :: v_esc_km_per_sec
        real(dp) :: v_esc

        real(dp) :: SHM_kinematic_c1
        real(dp) :: SHM_kinematic_c2

        character(len=512) :: vel_distribution_name

        contains

            procedure :: set_defaults => astroph_model_type_set_defaults
            procedure :: get_values => astroph_model_type_get_values
            procedure :: kinematic_function => astroph_model_type_kinematic_function

    end type

contains

    function astroph_model_type_kinematic_function(self, v_m_list, qm1_mag_list) result ( g0 )

        implicit none

        class(astroph_model_t) :: self

        real(dp), intent(in) :: v_m_list(:)
        real(dp), intent(in) :: qm1_mag_list(:)

        real(dp) :: g0(size(v_m_list))

        select case( trim(adjustl(self%vel_distribution_name)) )

            case ( 'SHM' )

                g0 = kinematic_function_SHM(self, v_m_list, qm1_mag_list)

        end select

    end function

    function kinematic_function_SHM(self, v_m_list, qm1_mag_list) result( g0 )

        implicit none

        class(astroph_model_t) :: self

        real(dp), intent(in) :: v_m_list(:)
        real(dp), intent(in) :: qm1_mag_list(:)

        real(dp) :: g0(size(v_m_list))
        integer :: q

        do q = 1, size(v_m_list)

            g0(q) = 0.0_dp

            if ( v_m_list(q) < self%v_esc )  then

                g0(q) = self%SHM_kinematic_c1*qm1_mag_list(q)*(exp(-(v_m_list(q)/self%v_0)**2) - self%SHM_kinematic_c2)

            end if

        end do

    end function

    subroutine astroph_model_type_set_defaults(self, cfg)

        use m_config

        implicit none

        class(astroph_model_t) :: self

        type(CFG_t) :: cfg

        call CFG_add(cfg,&
                     "astroph_model%v_0_km_per_sec", &
                     230.0_dp, &
                     "Dark matter SHM velocity distribution parameter, $v_0$.<br />"//&
                     "<ul>"//&
                     "<li><b>Units</b>: $\text{km}/\text{s}$</li>"//&
                     "</ul>")

        call CFG_add(cfg,&
                     "astroph_model%v_e_km_per_sec", &
                     [ 0.0_dp, 0.0_dp, 240.0_dp ], &
                     "List of Earth velocity vectors, $\mathbf{v}_e$.<br />"//&
                     "<ul>"//&
                     "<li><b>Units</b>: $\text{km}/\text{s}$</li>"//&
                     "<li><b>Dim</b>: [ : , 3]</li>"//&
                     "</ul>", dynamic_size = .TRUE.)

        call CFG_add(cfg,&
                     "astroph_model%v_esc_km_per_sec", &
                     600.0_dp, &
                     "Dark matter SHM velocity distribution parameter, $v_\mathrm{esc}$.<br />"//&
                     "<ul>"//&
                     "<li><b>Units</b>: $\text{km}/\text{s}$</li>"//&
                     "</ul>")

        call CFG_add(cfg, &
                     "astroph_model%vel_distribution_name", &
                     "SHM", &
                     "Specify the velocity distribution to use in the calculation.")

    end subroutine

    subroutine astroph_model_type_get_values(self, cfg)

        use m_config

        use constants_util

        implicit none

        class(astroph_model_t) :: self
        type(CFG_t) :: cfg

        integer :: n 

        real(dp) :: N0

        real(dp), allocatable :: v_e_list(:)

        call CFG_get(cfg, "astroph_model%v_0_km_per_sec", self%v_0_km_per_sec)
        self%v_0 = self%v_0_km_per_sec*km_per_sec_to_none

        call CFG_get_size(cfg, "astroph_model%v_e_km_per_sec", n)

        allocate(v_e_list(n))
        allocate(self%v_e_km_per_sec(n/3, 3), source = 0.0_dp)
        allocate(self%v_e_list(n/3, 3), source = 0.0_dp)

        call CFG_get(cfg, "astroph_model%v_e_km_per_sec", v_e_list)
        self%v_e_km_per_sec = transpose(reshape( v_e_list, [ 3, n/3 ] ) )
        self%v_e_list = self%v_e_km_per_sec*km_per_sec_to_none

        call CFG_get(cfg, "astroph_model%v_esc_km_per_sec", self%v_esc_km_per_sec)
        self%v_esc = self%v_esc_km_per_sec*km_per_sec_to_none

        N0 = (pi*self%v_0**3)*&
            (sqrt(pi)*erf(self%v_esc/self%v_0) - 2*(self%v_esc/self%v_0)*exp(-(self%v_esc/self%v_0)**2))
        self%SHM_kinematic_c1 = (2*pi**2*self%v_0**2/N0)
        self%SHM_kinematic_c2 = exp(-(self%v_esc/self%v_0)**2)

        call CFG_get(cfg, "astroph_model%vel_distribution_name", self%vel_distribution_name)

    end subroutine

end module
