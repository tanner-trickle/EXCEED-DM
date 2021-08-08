module material_type
    !! Defines the `material` type.

    use hdf5
    use h5lt

    use info_messages
    use prec
    use units

    implicit none

    type material_t
        !! Target material properties.

        character(len=64) :: name = ''
            !! Material name

        real(dp) :: pc_vol_A
            !! Volume of the unit cell 
            !!
            !! pc_vol_A = det(a_vecs)
            !!
            !! Units : Ang^3

        real(dp) :: pc_vol
            !! Volume of the primitive cell 
            !!
            !! Units : eV^(-3)

        real(dp) :: rho_T_g_per_cm3
            !! Target density
            !!
            !! Units : g/cm^3
        real(dp) :: rho_T
            !! Target density
            !!
            !! Units : eV^4
        real(dp) :: band_gap = 0.0_dp
            !! Band gap of the target
            !!
            !! Units : eV

        contains

            procedure :: load => load_material_nml
            procedure :: print => print_material
            procedure :: save => save_material

    end type


contains

    subroutine print_material(self, verbose)
        !! Prints `material` components.

        implicit none

        class(material_t) :: self
        logical, optional :: verbose

        if ( verbose ) then
            call print_section_seperator()
            print*, '    --------'
            print*, '    Material'
            print*, '    --------'
            print*
            print*, '        Name      : ', trim(self%name)
            print*, '        Band gap  : ', self%band_gap, 'eV'
            print*, '        Density   : ', self%rho_T_g_per_cm3, 'g/cm^3'
            print*, '        PC volume : ', self%pc_vol_A, 'Ang^3'
            print*
            call print_section_seperator()
            print*
        end if

    end subroutine

    subroutine load_material_nml(self, filename, verbose)
        !! Loads `material` parameters from a namelist.
        implicit none

        class(material_t) :: self
        character(len=*) :: filename
        logical, optional :: verbose

        logical :: file_exists
        integer :: error

        character(len=64) :: name = ''

        real(dp) :: pc_vol_A
            !! Volume of the unit cell 
            !!
            !! pc_vol_A = det(a_vecs)
            !!
            !! Units : Ang^3

        real(dp) :: rho_T_g_per_cm3
            !! Target density
            !!
            !! Units : g/cm^3

        real(dp) :: band_gap = 0.0_dp
            !! Band gap of the target
            !!
            !! Units : eV

        NAMELIST /material/ pc_vol_A        , &
                            band_gap        , &
                            rho_T_g_per_cm3 , &
                            name

        if ( verbose ) then
            print*, 'Loading material parameters...'
            print*
        end if

        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            open(100, file = trim(filename), iostat = error)
            read(100, nml=material, iostat=error)
            close(100)

            if ( error /= 0 ) then

                call print_error_message(&
                    'Problem reading material namelist.', &
                    verbose = verbose)
                stop

            end if

            self%name            = name
            self%pc_vol_A        = pc_vol_A
            self%band_gap        = band_gap
            self%rho_T_g_per_cm3 = rho_T_g_per_cm3

            self%rho_T  = g_to_eV*inv_cm_to_eV**3*rho_T_g_per_cm3
            self%pc_vol = Ang_to_inv_eV**3*pc_vol_A

            call self%print(verbose = verbose)

        else

            call print_error_message(&
                'Input file for material parameters : '//trim(filename)//' does NOT exist.', &
                verbose = verbose)
            stop

        end if

    end subroutine

    subroutine save_material(self, filename, verbose)
        !! Saves `material`.
        implicit none

        class(material_t) :: self
        character(len=*) :: filename
        logical, optional :: verbose

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id
        logical :: file_exists
        integer(HSIZE_T) :: dims1(1) = [1]
        integer :: error

        if ( verbose ) then
            print*, 'Saving material parameters...'
            print*
        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'material', group_id, error)

            ! write data
            call h5ltmake_dataset_string_f(file_id, &
                'material/name', &
                self%name, &
                error)
            call h5ltmake_dataset_double_f(file_id, &
                'material/band_gap', &
                size(dims1), dims1, &
                self%band_gap, &
                error)
            call h5ltmake_dataset_double_f(file_id, &
                'material/density', &
                size(dims1), dims1, &
                self%rho_T_g_per_cm3, &
                error)
            call h5ltmake_dataset_double_f(file_id, &
                'material/pc_vol', &
                size(dims1), dims1, &
                self%pc_vol_A, &
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
