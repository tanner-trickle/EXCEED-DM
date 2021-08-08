module numerics_abs
    !! Numerics parameters specific to DM absorption rate calculations.

    use info_messages
    use math_mod
    
    use physics_abs_functions

    implicit none

    type numerics_abs_t
        !! Numerics parameters specific to DM absorption calculations.

        integer :: n_val_max = 0
            !! Maximum number of valence bands
        integer :: n_cond_max = 0
            !! Maximum number of conduction bands
        integer, allocatable :: val_id_list(:)
            !! Dim : [n_val_max]
            !!
            !! List of valence ID's to compute for.
        integer, allocatable :: k_id_list(:)
            !! Dim : [n_k]
            !!
            !! List of k's to compute for.
        logical :: save_tran_form = .FALSE.
            !! If .TRUE. will save the transition form factors to the output file.
        integer :: n_v
            !! Number of \( v \) points in integration over \( \mathbf{v} \)
        integer :: n_v_theta
            !! Number of \( \theta \) points in integration over \( \mathbf{v} \)
        integer :: n_v_phi
            !! Number of \( \phi \) points in integration over \( \mathbf{v} \)

        real(dp), allocatable :: v_angular_mesh(:, :)
            !! Dim : [ n_v_theta*n_v_phi, 2]
            !!
            !! List of \( \theta, \phi \) points in integral over \( \mathbf{v} \)
            !!
            !! Units : None

        real(dp), allocatable :: v_list(:)
            !! Dim : [n_v]
            !!
            !! List of \( |\mathbf{v}| \) points in velocity integral
            !!
            !! Units : None
        real(dp), allocatable :: v_mesh(:, :)
            !! Dim : [ n_v*n_v_theta*n_v_phi, 3]
            !!
            !! All points in the \( \mathbf{v} \) integration.
            !!
            !! Units : None

        contains

            procedure :: print => numerics_abs_print
            procedure :: load => numerics_abs_load_nml
            procedure :: save => numerics_abs_save
            procedure :: create_val_id_list => abs_create_val_id_list
            procedure :: create_k_id_list => abs_create_k_id_list

    end type

contains

    subroutine abs_create_k_id_list(self, n_k)
        !! Specify the indicies for each \( \mathbf{k} \) point that should be included. 
        !! Specific to the DM absorption rate calculation. 

        implicit none

        class(numerics_abs_t) :: self
        integer :: n_k

        integer :: k

        allocate(self%k_id_list(n_k))

        do k = 1, n_k
            self%k_id_list(k) = k
        end do

    end subroutine

    subroutine abs_create_val_id_list(self, n_val)
        !! Specify the indicies for each valence band that should be included. 
        !! Specific to the DM absorption rate calculation. 

        implicit none
        class(numerics_abs_t) :: self

        integer :: n_val

        integer :: j

        allocate(self%val_id_list(self%n_val_max))

        do j = 1, self%n_val_max
            self%val_id_list(j) = n_val - j + 1
        end do

    end subroutine

    subroutine numerics_abs_save(self, filename, verbose)
        !! Saves `numerics_abs`.

        use hdf5
        use h5lt

        implicit none

        class(numerics_abs_t) :: self
        character(len=*) :: filename
        logical, optional :: verbose

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id
        logical :: file_exists
        integer(HSIZE_T) :: dims1(1) = [1]
        integer :: error

        if ( verbose ) then
            print*, 'Saving numerics - abs parameters...'
            print*
        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'numerics_abs', group_id, error)

            ! write data
            call h5ltmake_dataset_int_f(file_id, &
                'numerics_abs/n_val_max', &
                size(dims1), dims1, &
                self%n_val_max, &
                error)
            call h5ltmake_dataset_int_f(file_id, &
                'numerics_abs/n_cond_max', &
                size(dims1), dims1, &
                self%n_cond_max, &
                error)

            call h5ltmake_dataset_int_f(file_id, &
                'numerics_abs/n_v', &
                size(dims1), dims1, &
                self%n_v, &
                error)
            call h5ltmake_dataset_int_f(file_id, &
                'numerics_abs/n_v_theta', &
                size(dims1), dims1, &
                self%n_v_theta, &
                error)
            call h5ltmake_dataset_int_f(file_id, &
                'numerics_abs/n_v_phi', &
                size(dims1), dims1, &
                self%n_v_phi, &
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

    subroutine numerics_abs_load_nml(self, filename, n_val, n_cond, dm_model, verbose)
        !! Loads `numerics_abs` parameters from a namelist.

        implicit none

        class(numerics_abs_t) :: self
        character(len=*) :: filename
        integer :: n_val, n_cond
        type(dm_model_t) :: dm_model
        logical, optional :: verbose

        logical :: file_exists
        integer :: error

        integer :: n_val_max = 0
            !! Maximum number of valence bands
        integer :: n_cond_max = 0
            !! Maximum number of conduction bands

        logical :: save_tran_form = .FALSE.

        integer :: n_v = 1
        integer :: n_v_theta = 1
        integer :: n_v_phi = 1

        integer :: v, a, v_id

        NAMELIST /numerics_abs/ n_val_max, &
                                 n_cond_max, &
                                 save_tran_form, &
                                 n_v, &
                                 n_v_theta, &
                                 n_v_phi

        if ( verbose ) then
            print*, 'Loading numerics parameters for absorption calculation...'
            print*
        end if

        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            open(100, file = trim(filename), iostat = error)
            read(100, nml = numerics_abs, iostat = error)
            close(100)

            if ( n_val_max == 0 ) then
                self%n_val_max = n_val
            else
                self%n_val_max = min(n_val, n_val_max)
            end if

            if ( n_cond_max == 0 ) then
                self%n_cond_max = n_cond
            else
                self%n_cond_max = min(n_cond, n_cond_max)
            end if

            self%save_tran_form = save_tran_form

            self%n_v = n_v
            self%n_v_theta = n_v_theta
            self%n_v_phi = n_v_phi

            allocate(self%v_angular_mesh(self%n_v_theta*self%n_v_phi, 2))
            self%v_angular_mesh = generate_uniform_points_on_sphere(self%n_v_theta, &
                self%n_v_phi) 

            allocate(self%v_list(self%n_v))
            self%v_list = uniform_list(self%n_v, 0.0_dp, dm_model%vX_max)

            ! create the v mesh
            allocate(self%v_mesh(self%n_v*self%n_v_theta*self%n_v_phi, 3))

            v_id = 0
            do v = 1, self%n_v

                do a = 1, size(self%v_angular_mesh, 1)

                    v_id = v_id + 1

                    self%v_mesh(v_id, 1) = self%v_list(v)*sin(self%v_angular_mesh(a, 1))*&
                        cos(self%v_angular_mesh(a, 2))
                    self%v_mesh(v_id, 2) = self%v_list(v)*sin(self%v_angular_mesh(a, 1))*&
                        sin(self%v_angular_mesh(a, 2))
                    self%v_mesh(v_id, 3) = self%v_list(v)*cos(self%v_angular_mesh(a, 1))

                end do

            end do

            ! chech the MB integration.
            call check_mb_dist_normalization(self%v_mesh, dm_model, &
                boost_vec_in = dm_model%vE*[0, 0, 1], verbose = verbose)

            call self%print(verbose = verbose)

        else

            call print_error_message(&
                'Input file for numerics-absorption parameters : '//trim(filename)//' does NOT exist.', &
                verbose = verbose)
            stop

        end if

    end subroutine

    subroutine numerics_abs_print(self, verbose)
        !! Prints `numerics_abs` components.

        implicit none

        class(numerics_abs_t) :: self
        logical, optional :: verbose

        if ( verbose ) then

            call print_section_seperator()
            print*, '    ---------------------'
            print*, '    Numerics - Absorption'
            print*, '    ---------------------'
            print*
            print*, '        Maximum number of valence bands    : ', trim(adjustl(int_to_str(self%n_val_max)))
            print*, '        Maximum number of conduction bands : ', trim(adjustl(int_to_str(self%n_cond_max)))
            print*
            print*, '        v Integration'
            print*, '            Number of |v| points   : ', trim(adjustl(int_to_str(self%n_v)))
            print*, '            Number of theta points : ', trim(adjustl(int_to_str(self%n_v_theta)))
            print*, '            Number of phi points   : ', trim(adjustl(int_to_str(self%n_v_phi)))
            print*
            call print_section_seperator()
            print*

        end if

    end subroutine

end module
