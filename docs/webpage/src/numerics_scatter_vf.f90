module numerics_scatter_vf
    !! Numerics parameters specific to the valence \( \rightarrow \) free DM-electron scattering rate calculation.

    use prec
    use info_messages

    implicit none

    type numerics_scatter_vf_t
        !! Numerics parameters specific to the valence \( \rightarrow \) free DM-electron scattering rate calculation.

        integer :: n_val_max = 0
            !! Maximum number of valence bands
        character(len=64) :: Zeff_type = 'one'
            !! Specify what Zeff to use
            !!
            !! - 'one' - all Zeff = 1
            !!
            !! - 'Eb' - use the binding energy of the valence states
        real(dp), allocatable :: Zeff(:, :)
            !! Dim : [ n_val, n_k ]
            !!
            !! \( Z_\text{eff} \) parameter in the Fermi form factor.
        integer, allocatable :: val_id_list(:)
            !! Dim : [n_val_max]
            !!
            !! List of valence ID's to compute for.
        integer, allocatable :: k_id_list(:)
            !! Dim : [n_k]
            !!
            !! List of k's to compute for.
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
            !! Number of \( \omega \) points to compute for. The valence -> free scattering
            !! rate calculation is done by computing `n_omega` values of \( \frac{dR}{d \omega} \) between
            !!
            !! \( \omega_\text{min} = E_{f, \text{max}} \)
            !! \( \omega_\text{max} = \frac{1}{2} m_\chi v_\text{max}^2 \)
            !!
            !! and then interpolated to the user defined grid.

        real(dp), allocatable :: omega_list(:)
            !! List of \( \omega \) points to compute for. The valence -> free scattering
            !! rate calculation is done by computing `n_omega` values of \( \frac{dR}{d \omega} \) between
            !!
            !! \( \omega_\text{min} = E_{f, \text{max}} \)
            !! \( \omega_\text{max} = \frac{1}{2} m_\chi v_\text{max}^2 \)
            !!
            !! and then interpolated to the user defined grid.

        contains

            procedure :: print => numerics_scatter_vf_print
            procedure :: load => numerics_scatter_vf_load_nml
            procedure :: save => numerics_scatter_vf_save
            procedure :: create_val_id_list => vf_create_val_id_list
            procedure :: create_k_id_list => vf_create_k_id_list

    end type

contains

    subroutine vf_create_k_id_list(self, n_k)
        !! Specify the indicies for each \( \mathbf{k} \) point that should be included. 
        !! Specific to the valence \( \rightarrow \) free DM-electron scattering rate calculations. 

        implicit none

        class(numerics_scatter_vf_t) :: self
        integer :: n_k

        integer :: k

        allocate(self%k_id_list(n_k))

        do k = 1, n_k
            self%k_id_list(k) = k
        end do

    end subroutine

    subroutine vf_create_val_id_list(self, n_val)
        !! Specify the indicies for each valence band that should be included. 
        !! Specific to the valence \( \rightarrow \) free DM-electron scattering rate calculation. 

        implicit none
        class(numerics_scatter_vf_t) :: self

        integer :: n_val

        integer :: j

        allocate(self%val_id_list(self%n_val_max))

        do j = 1, self%n_val_max
            self%val_id_list(j) = n_val - j + 1
        end do

    end subroutine

    subroutine numerics_scatter_vf_save(self, filename, verbose)
        !! Saves `numerics_scatter_vf`.

        use hdf5
        use h5lt

        implicit none

        class(numerics_scatter_vf_t) :: self
        character(len=*) :: filename
        logical, optional :: verbose

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id
        logical :: file_exists
        integer(HSIZE_T) :: dims1(1) = [1]
        integer(HSIZE_T) :: dims2(2)
        integer :: error

        if ( verbose ) then
            print*, 'Saving numerics - scatter-vf parameters...'
            print*
        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'numerics_scatter_vf', group_id, error)

            ! write data
            call h5ltmake_dataset_int_f(file_id, &
                'numerics_scatter_vf/n_val_max', &
                size(dims1), dims1, &
                self%n_val_max, &
                error)
            call h5ltmake_dataset_int_f(file_id, &
                'numerics_scatter_vf/n_kf_theta', &
                size(dims1), dims1, &
                self%n_kf_theta, &
                error)
            call h5ltmake_dataset_int_f(file_id, &
                'numerics_scatter_vf/n_kf_phi', &
                size(dims1), dims1, &
                self%n_kf_phi, &
                error)
            call h5ltmake_dataset_int_f(file_id, &
                'numerics_scatter_vf/n_omega', &
                size(dims1), dims1, &
                self%n_omega, &
                error)

            dims1 = [size(self%omega_list)]
            call h5ltmake_dataset_double_f(file_id, &
                'numerics_scatter_vf/omega_list', &
                size(dims1), dims1, &
                self%omega_list, &
                error)

            dims2 = [ size(self%Zeff, 1), size(self%Zeff, 2) ]
            call h5ltmake_dataset_double_f(file_id, &
                'numerics_scatter_vf/Zeff', &
                size(dims2), dims2, &
                self%Zeff, &
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

    subroutine numerics_scatter_vf_load_nml(self, filename, PW_dataset, verbose)
        !! Loads `numerics_scatter_vf` parameters from a namelist.

        use PW_dataset_type
        use math_mod

        implicit none

        class(numerics_scatter_vf_t) :: self
        character(len=*) :: filename
        type(PW_dataset_t) :: PW_dataset
        integer :: n_val, n_cond
        logical, optional :: verbose

        logical :: file_exists
        integer :: error

        integer :: n_val_max = 0
            !! Maximum number of valence bands
        integer :: n_kf_theta = 1
            !! Number of \( \theta \) points in the integration of \( \mathbf{k}_f \)
        integer :: n_kf_phi = 1
            !! Number of \( \phi \) points in the integration of \( \mathbf{k}_f \)
        integer :: n_omega = 1

        character(len=64) :: Zeff_type = 'one'

        real(dp) :: n_Eb(PW_dataset%n_val)

        integer :: i, k

        NAMELIST /numerics_s_vf/ n_val_max, &
                                 n_kf_theta, &
                                 n_kf_phi, &
                                 Zeff_type, &
                                 n_Eb, &
                                 n_omega


        if ( verbose ) then
            print*, 'Loading numerics parameters for scatter-vc calculation...'
            print*
        end if

        allocate(self%Zeff(PW_dataset%n_val, PW_dataset%n_k))
        self%Zeff = 1.0_dp

        n_Eb = 1.0_dp

        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            open(100, file = trim(filename), iostat = error)
            read(100, nml = numerics_s_vf, iostat = error)
            close(100)

            if ( n_val_max == 0 ) then
                self%n_val_max = PW_dataset%n_val
            else
                self%n_val_max = min(PW_dataset%n_val, n_val_max)
            end if

            self%n_kf_theta = n_kf_theta
            self%n_kf_phi = n_kf_phi
            self%n_omega = n_omega

            allocate(self%omega_list(self%n_omega))
            allocate(self%kf_angular_mesh(self%n_kf_theta*self%n_kf_phi, 2))

            self%kf_angular_mesh = generate_uniform_points_on_sphere(&
                self%n_kf_theta, self%n_kf_phi)

            self%Zeff_type = trim(Zeff_type)

            if ( trim(self%Zeff_type) == 'Eb' ) then

                do i = 1, PW_dataset%n_val
                
                    do k = 1, PW_dataset%n_k

                        self%Zeff(i, k) = n_Eb(i)*sqrt(&
                            ( maxval(PW_dataset%energy_bands(:, :PW_dataset%n_val)) &
                                - PW_dataset%energy_bands(k, i) )&
                            /13.6_dp) 

                        self%Zeff(i, k) = max(self%Zeff(i, k), 1.0_dp)

                    end do

                end do

            end if

            call self%print(verbose = verbose)

        else

            call print_error_message(&
                'Input file for numerics-scatter-vc parameters : '//trim(filename)//' does NOT exist.', &
                verbose = verbose)
            stop

        end if

    end subroutine

    subroutine numerics_scatter_vf_print(self, verbose)
        !! Prints `numerics_scatter_vf` components.

        implicit none

        class(numerics_scatter_vf_t) :: self
        logical, optional :: verbose

        if ( verbose ) then

            call print_section_seperator()
            print*, '    -----------------------'
            print*, '    Numerics - scatter - vf'
            print*, '    -----------------------'
            print*
            print*, '        Maximum number of valence bands    : ', trim(adjustl(int_to_str(self%n_val_max)))
            print*
            print*, '        Number of kf theta points : ', trim(adjustl(int_to_str(self%n_kf_theta)))
            print*, '        Number of kf phi points   : ', trim(adjustl(int_to_str(self%n_kf_phi)))
            print*, '        Number of omega points    : ', trim(adjustl(int_to_str(self%n_omega)))
            print*
            print*, '        Zeff type : ', trim(self%Zeff_type)
            print*
            call print_section_seperator()
            print*

        end if

    end subroutine

end module
