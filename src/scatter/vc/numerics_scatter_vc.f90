module numerics_scatter_vc
    !! Numerics parameters specific to the valence \( \rightarrow \) conduction DM-electron scattering rate calculation.

    use info_messages

    implicit none

    type numerics_scatter_vc_t
        !! Numerics parameters specific to the valence \( \rightarrow \) conduction DM-electron scattering rate calculation.

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

        contains

            procedure :: print => numerics_scatter_vc_print
            procedure :: load => numerics_scatter_vc_load_nml
            procedure :: save => numerics_scatter_vc_save
            procedure :: create_val_id_list => vc_create_val_id_list
            procedure :: create_k_id_list => vc_create_k_id_list

    end type

contains

    subroutine vc_create_k_id_list(self, n_k)
        !! Specify the indicies for each \( \mathbf{k} \) point that should be included. 
        !! Specific to the valence \( \rightarrow \) conduction DM-electron rate calculations. 

        implicit none

        class(numerics_scatter_vc_t) :: self
        integer :: n_k

        integer :: k

        allocate(self%k_id_list(n_k))

        do k = 1, n_k
            self%k_id_list(k) = k
        end do

    end subroutine

    subroutine vc_create_val_id_list(self, n_val)
        !! Specify the indicies for each valence band that should be included. 
        !! Specific to the valence \( \rightarrow \) conduction DM-electron scattering rate calculation. 

        implicit none
        class(numerics_scatter_vc_t) :: self

        integer :: n_val

        integer :: j

        allocate(self%val_id_list(self%n_val_max))

        do j = 1, self%n_val_max
            self%val_id_list(j) = n_val - j + 1
        end do

    end subroutine

    subroutine numerics_scatter_vc_save(self, filename, verbose)
        !! Saves `numerics_scatter_vc`.

        use hdf5
        use h5lt

        implicit none

        class(numerics_scatter_vc_t) :: self
        character(len=*) :: filename
        logical, optional :: verbose

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id
        logical :: file_exists
        integer(HSIZE_T) :: dims1(1) = [1]
        integer :: error

        if ( verbose ) then
            print*, 'Saving numerics - scatter-vc parameters...'
            print*
        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'numerics_scatter_vc', group_id, error)

            ! write data
            call h5ltmake_dataset_int_f(file_id, &
                'numerics_scatter_vc/n_val_max', &
                size(dims1), dims1, &
                self%n_val_max, &
                error)
            call h5ltmake_dataset_int_f(file_id, &
                'numerics_scatter_vc/n_cond_max', &
                size(dims1), dims1, &
                self%n_cond_max, &
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

    subroutine numerics_scatter_vc_load_nml(self, filename, n_val, n_cond, verbose)
        !! Loads `numerics_scatter_vc` parameters from a namelist.

        implicit none

        class(numerics_scatter_vc_t) :: self
        character(len=*) :: filename
        integer :: n_val, n_cond
        logical, optional :: verbose

        logical :: file_exists
        integer :: error

        integer :: n_val_max = 0
            !! Maximum number of valence bands
        integer :: n_cond_max = 0
            !! Maximum number of conduction bands

        NAMELIST /numerics_s_vc/ n_val_max, &
                                 n_cond_max       

        if ( verbose ) then
            print*, 'Loading numerics parameters for scatter-vc calculation...'
            print*
        end if

        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            open(100, file = trim(filename), iostat = error)
            read(100, nml = numerics_s_vc, iostat = error)
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

            call self%print(verbose = verbose)

        else

            call print_error_message(&
                'Input file for numerics-scatter-vc parameters : '//trim(filename)//' does NOT exist.', &
                verbose = verbose)
            stop

        end if

    end subroutine

    subroutine numerics_scatter_vc_print(self, verbose)
        !! Prints `numerics_scatter_vc` components.

        implicit none

        class(numerics_scatter_vc_t) :: self
        logical, optional :: verbose

        if ( verbose ) then

            call print_section_seperator()
            print*, '    -----------------------'
            print*, '    Numerics - Scatter - vc'
            print*, '    -----------------------'
            print*
            print*, '        Maximum number of valence bands    : ', trim(adjustl(int_to_str(self%n_val_max)))
            print*, '        Maximum number of conduction bands : ', trim(adjustl(int_to_str(self%n_cond_max)))
            print*
            call print_section_seperator()
            print*

        end if

    end subroutine

end module
