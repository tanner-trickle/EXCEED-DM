module material_input
    !! Information about the target material.

    use hdf5
    use h5lt

    use prec
    use units

    implicit none

    character(len=64) :: mat_name = ''

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

    real(dp) :: m_T_kg = 1.0_dp
        !! Target mass
        !!
        !! Units : kg 

    real(dp) :: band_gap = 0.0_dp
        !! Band gap of the target
        !!
        !! Units : eV

    !! Generated from input

    real(dp) :: pc_vol
        !! Volume of the primitive cell 
        !!
        !! Units : eV^(-3)

    real(dp) :: rho_T
        !! Target density
        !!
        !! Units : eV^4

    real(dp) :: m_T
        !! Target mass
        !!
        !! Units : eV

    NAMELIST /material/ pc_vol_A, &
                        band_gap, &
                        rho_T_g_per_cm3, &
                        m_T_kg, &
                        mat_name

contains

    subroutine print_material(verbose)

        implicit none

        logical, optional :: verbose

        if ( verbose ) then
            print*, '----------------------------------------'
            print*, '    --------'
            print*, '    Material'
            print*, '    --------'
            print*
            print*, '        Name      : ', trim(mat_name)
            print*, '        Band gap  : ', band_gap, ' eV'
            print*, '        Density   : ', rho_T_g_per_cm3, ' g/cm^3'
            print*, '        PC volume : ', pc_vol_A, ' Ang^3'
            print*
        end if

    end subroutine

    subroutine load_material(filename, verbose)
        implicit none

        character(len=*) :: filename

        logical :: file_exists

        logical, optional :: verbose

        integer :: error

        if ( verbose ) then
            print*, 'Loading material parameters...'
            print*
        end if

        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            open(100, file = trim(filename), iostat = error)

            read(100, nml=material, iostat=error)

            close(100)

            if ( error .ne. 0 ) then

                if ( verbose ) then

                    print*, '!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                    print*
                    print*, '    Problem reading material namelist.'
                    print*
                    print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                    print*

                end if

                stop

            end if

            rho_T = g_to_eV*inv_cm_to_eV**3*rho_T_g_per_cm3  
            m_T = kg_to_eV*m_T_kg 
            pc_vol = Ang_to_inv_eV**3*pc_vol_A 

            call print_material(verbose = verbose)

            if ( verbose ) then

                print*, '----------------------------------------'
                print*

            end if

        else

            if ( verbose ) then

                print*, '!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*
                print*, '    Input file for material parameters : ', trim(filename), ' does NOT exist.'
                print*
                print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*

            end if

            stop

        end if

    end subroutine

    subroutine save_material(filename, verbose)
        !! Saves the material properties
        implicit none

        character(len=*) :: filename

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id

        logical, optional :: verbose
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

            ! ! write data
            call h5ltmake_dataset_string_f(file_id, 'material/mat_name',&
                mat_name, error)
            call h5ltmake_dataset_double_f(file_id, 'material/band_gap', size(dims1), dims1,&
                band_gap, error)
            call h5ltmake_dataset_double_f(file_id, 'material/density', size(dims1), dims1,&
                rho_T_g_per_cm3, error)
            call h5ltmake_dataset_double_f(file_id, 'material/target_mass', size(dims1), dims1,&
                m_T_kg, error)

            call h5fclose_f(file_id, error)
            call h5close_f(error)

        else

            if ( verbose ) then

                print*, '!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*
                print*, '    Output file : ', trim(filename), ' does NOT exist.'
                print*
                print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*

            end if

            stop

        end if

    end subroutine

end module
