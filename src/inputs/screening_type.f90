module screening_type

    use prec_util, only: dp

    implicit none

    type :: screening_t

        character(len=512) :: type

        character(len=512) :: dielectric_filename

        ! analytic parameters
        real(dp) :: e0

        real(dp) :: q_tf_keV
        real(dp) :: q_tf

        real(dp) :: omega_p

        real(dp) :: alpha

        ! numerics parameters
        integer :: width_id

        complex(dp), allocatable :: num_dielectric(:, :, :, :)
        complex(dp), allocatable :: num_dielectric_ani(:, :, :, :, :, :)

        real(dp) :: q_bin_width
        real(dp) :: E_bin_width

        contains 

            procedure :: set_defaults => screening_type_set_defaults
            procedure :: get_values => screening_type_get_values
            procedure :: screening_factor => screening_type_screening_factor

    end type

contains

    function screening_type_screening_factor(self, q_vec_list, omega) result ( screen_factor_list )

        use constants_util
        use math_util

        implicit none

        class(screening_t) :: self

        real(dp) :: q_vec_list(:, :)
        real(dp) :: q_hat(3)
        real(dp) :: q_mag_list(size(q_vec_list, 1))
        real(dp) :: omega
        real(dp) :: q_phi, q_theta

        integer :: q_bin, q_theta_bin, q_phi_bin, E_bin

        real(dp) :: screen_factor_list(size(q_vec_list, 1))

        integer :: q

        real(dp) :: q_bin_width

        q_mag_list = norm2(q_vec_list, 2) + 1.0e-8_dp

        select case ( trim(adjustl(self%type)) )

            case ( 'analytic' )

                screen_factor_list = 1.0_dp + ( (self%e0 - 1.0_dp)**(-1) + &
                    self%alpha*(q_mag_list/self%q_tf)**2 + &
                    q_mag_list**4/(4*m_elec**2*self%omega_p**2) - &
                    (omega/self%omega_p)**2)**(-1) 

            case ( 'numeric' )

                q_bin_width = self%q_bin_width*1.0e3_dp

                if ( allocated(self%num_dielectric_ani) ) then

                    do q = 1, size(q_vec_list, 1)

                        q_hat = q_vec_list(q, :)/q_mag_list(q)
                        
                        q_theta = acos(q_hat(3))
                        q_phi = mod(atan2( q_hat(2), q_hat(1) ) + 2.0_dp*pi, 2.0_dp*pi)

                        q_theta_bin = 1 + floor(q_theta/(pi/size(self%num_dielectric_ani, 2)))
                        q_phi_bin = 1 + floor(q_phi/(2.0_dp*pi/size(self%num_dielectric_ani, 3)))
                        q_bin = 1 + floor( q_mag_list(q)/q_bin_width )
                        E_bin = 1 + floor(omega/self%E_bin_width)

                        if ( ( q_theta_bin <= size(self%num_dielectric_ani, 2) ) .and. &
                             ( q_phi_bin <= size(self%num_dielectric_ani, 3) ) .and. &
                             ( q_bin <= size(self%num_dielectric_ani, 1) ) .and. &
                             ( E_bin <= size(self%num_dielectric_ani, 4) ) ) then

                             screen_factor_list(q) = abs(&
                                 dot_product(q_hat, &
                                    matmul( self%num_dielectric_ani(q_bin, q_theta_bin, q_phi_bin, E_bin, :, :), &
                                    q_hat ) ) )

                        else

                            screen_factor_list(q) = 1.0_dp

                        end if

                    end do

                else if ( allocated(self%num_dielectric) ) then

                    do q = 1, size(q_vec_list, 1)

                        q_hat = q_vec_list(q, :)/q_mag_list(q)
                        
                        q_theta = acos(q_hat(3))
                        q_phi = mod(atan2( q_hat(2), q_hat(1) ) + 2.0_dp*pi, 2.0_dp*pi)

                        q_theta_bin = 1 + floor(q_theta/(pi/size(self%num_dielectric, 2)))
                        q_phi_bin = 1 + floor(q_phi/(2.0_dp*pi/size(self%num_dielectric, 3)))
                        q_bin = 1 + floor( q_mag_list(q)/q_bin_width )
                        E_bin = 1 + floor(omega/self%E_bin_width)

                        if ( ( q_theta_bin <= size(self%num_dielectric, 2) ) .and. &
                             ( q_phi_bin <= size(self%num_dielectric, 3) ) .and. &
                             ( q_bin <= size(self%num_dielectric, 1) ) .and. &
                             ( E_bin <= size(self%num_dielectric, 4) ) ) then

                             screen_factor_list(q) = abs(self%num_dielectric(q_bin, q_theta_bin, q_phi_bin, E_bin) )

                        else

                            screen_factor_list(q) = 1.0_dp

                        end if

                    end do

                else

                    screen_factor_list = 1.0_dp

                end if

            case default

                screen_factor_list = 1.0_dp

        end select

    end function

    subroutine screening_type_set_defaults(self, cfg)

        use m_config

        implicit none

        class(screening_t) :: self

        type(CFG_t) :: cfg

        call CFG_add(cfg, &
            "screening%type", &
            "", &
            "Model to use for the screening factor. Default is a screening factor of 1.")

        call CFG_add(cfg, &
            "screening%dielectric_filename", &
            "", &
            "Location of the numerically computed dielectric to use as the screening factor.")

        call CFG_add(cfg, &
            "screening%width_id", &
            1, &
            "Width ID number to use for the numeric dielectric.")

        call CFG_add(cfg, &
            "screening%e0", &
            1.0_dp, &
            "Static dielectric parameter for analytic screening, $\epsilon(0)$ in Eq. (6) of"&
    //"[https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892)"&
           //"<br />Units : None")

        call CFG_add(cfg, &
            "screening%q_tf", &
            1.0_dp, &
            "Thomas-Fermi momentum parameter for analytic screening, $q_\text{TF}$ in Eq. (6) of"&
    //"[https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892)"&
           //"<br />Units : keV")

        call CFG_add(cfg, &
            "screening%omega_p", &
            1.0_dp, &
            "Plasma frequency parameter for analytic screening, $\omega_p$ in Eq. (6) of"&
    //" [https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892)"&
           //"<br />Units : eV")

        call CFG_add(cfg, &
            "screening%alpha", &
            1.0_dp, &
            "Shape parameter for analytic screening, $\alpha$ in Eq. (6) of"&
    //" [https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892)"&
           //"<br />Units : None")

        call CFG_add(cfg, &
            "screening%q_bin_width", &
            1.0_dp, &
            "Width of the bins in $q$. Only used for the numeric screening.<br />Units : keV")

        call CFG_add(cfg, &
            "screening%E_bin_width", &
            1.0_dp, &
            "Width of the bins in $\omega$. Only used for the numeric screening.<br />Units : eV")

    end subroutine

    subroutine screening_type_get_values(self, cfg)

        use m_config

        implicit none

        class(screening_t) :: self
        type(CFG_t) :: cfg

        call CFG_get(cfg, "screening%type", self%type)
        call CFG_get(cfg, "screening%dielectric_filename", self%dielectric_filename)
        call CFG_get(cfg, "screening%width_id", self%width_id)

        if ( ( trim(adjustl(self%type)) == 'numeric' ) .and. &
             ( trim(adjustl(self%dielectric_filename) ) /= "" ) ) then

            call load_numeric_dielectric(self, self%dielectric_filename)

        end if

        call CFG_get(cfg, "screening%e0", self%e0)
        call CFG_get(cfg, "screening%q_tf", self%q_tf_keV)
        self%q_tf = self%q_tf_keV*1.0e3_dp

        call CFG_get(cfg, "screening%omega_p", self%omega_p)
        call CFG_get(cfg, "screening%alpha", self%alpha)

        call CFG_get(cfg, "screening%q_bin_width", self%q_bin_width)
        call CFG_get(cfg, "screening%E_bin_width", self%E_bin_width)

    end subroutine

    subroutine load_numeric_dielectric(screening, filename)

        use hdf5_utils

        use constants_util

        implicit none

        class(screening_t) :: screening

        character(len=*) :: filename 

        integer(HID_T) :: file_id

        integer :: rank

        integer :: dims(6)

        real(dp), allocatable :: num_dielectric_r(:, :, :, :)
        real(dp), allocatable :: num_dielectric_c(:, :, :, :)

        real(dp), allocatable :: num_dielectric_ani_r(:, :, :, :, :, :)
        real(dp), allocatable :: num_dielectric_ani_c(:, :, :, :, :, :)

        character(len=512) :: dielectric_path
        character(len=512) :: w_str

        call hdf_open_file(file_id, filename, status='OLD', action='READ')

        write(w_str, *) screening%width_id

        dielectric_path = "dielectric/width_"//trim(adjustl(w_str))

        call hdf_get_rank(file_id, trim(adjustl(dielectric_path))//"/dielectric_r", rank)
        call hdf_get_dims(file_id, trim(adjustl(dielectric_path))//"/dielectric_r", dims)

        if ( rank == 6 ) then

            allocate(screening%num_dielectric_ani(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)), &
                source = (0.0_dp, 0.0_dp))
            allocate(num_dielectric_ani_r(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)), &
                source = 0.0_dp)
            allocate(num_dielectric_ani_c(dims(1), dims(2), dims(3), dims(4), dims(5), dims(6)), &
                source = 0.0_dp)

            call hdf_read_dataset(file_id, trim(adjustl(dielectric_path))//"/dielectric_r", &
                                    num_dielectric_ani_r)

            call hdf_read_dataset(file_id, trim(adjustl(dielectric_path))//"/dielectric_c", &
                                    num_dielectric_ani_c)

            screening%num_dielectric_ani = num_dielectric_ani_r + ii*num_dielectric_ani_c

        else if ( rank == 4 ) then

            allocate(screening%num_dielectric(dims(1), dims(2), dims(3), dims(4)), &
                source = (0.0_dp, 0.0_dp))
            allocate(num_dielectric_r(dims(1), dims(2), dims(3), dims(4)), &
                source = 0.0_dp)
            allocate(num_dielectric_c(dims(1), dims(2), dims(3), dims(4)), &
                source = 0.0_dp)

            call hdf_read_dataset(file_id, trim(adjustl(dielectric_path))//"/dielectric_r", &
                                    num_dielectric_r)

            call hdf_read_dataset(file_id, trim(adjustl(dielectric_path))//"/dielectric_c", &
                                    num_dielectric_c)

            screening%num_dielectric = num_dielectric_r + ii*num_dielectric_c

        end if

        call hdf_close_file(file_id)

    end subroutine

end module
