module numerics_input
    !! Variables relating to binning and other numerics

    use hdf5
    use h5lt

    use prec

    implicit none

    real(dp) :: q_bin_width = 1.0e3_dp
        !! width of momuntum bins
        !! 
        !! Units : eV

    integer :: n_q_bins = 0
        !! Number of q bins, with bin width q_bin_width, to save
        !!
        !! Note: momentum transfers greater than n_q_bins*q_bin_width
        !! will be stored in the last, extra bin

    real(dp) :: E_bin_width = 1.0_dp
        !! Width of the energy bins
        !!
        !! Units : eV

    integer :: n_E_bins = 0
        !! Number of E bins, with bin width E_bin_width, to save
        !!
        !! Note: energy transfers greater than n_E_bins*E_bin_width
        !! will be stored in the last, extra bin

    integer :: E_bin_threshold = 1
        !! smallest E bin added to total rate calculation
        !! default : 1 - include all bins

    real(dp) :: Ef_max = 60.0_dp
        !! maximum final state energy to include in calculations
        !!
        !! Units : eV
        !!
        !! TODO : E_cut is a more appropriate name

    integer :: n_init
        !! number of initial states
    integer :: n_fin
        !! number of final states

    integer :: n_FFT_grid_input(3) = [0, 0, 0]
        !! Default size of the FFT, used to extend core
        !! calculation by making the size of the FFT larger
        !!
        !! Note : only used in calculations which use the FFT

    !real(dp) :: q_s_FFT = 2151.0_dp/2.0_dp
    !    !! Units : eV
    !    !!
    !    !! q_max_FFT = q_s_FFT*N_FFT
    !    !!
    !    !! Minimum eigenvalue of k_red_to_xyz matrix
    !    !!
    !    !! Note : only used in calculations which use the FFT

    integer :: n_kf_theta = 1
        !! Number of theta points in integration of kf
        !!
        !! Note : only used in free calculatinos
    integer :: n_kf_phi = 1
        !! Number of phi points in integration of kf
        !!
        !! Note : only used in free calculatinos

    integer :: n_ki_theta = 1
        !! Number of theta points in integration of ki
        !!
        !! Note : only used in core -> free calculatinos
    integer :: n_ki_phi = 1
        !! Number of phi points in integration of ki
        !!
        !! Note : only used in core -> free calculatinos

    real(dp) :: ki_s = 100.0_dp
        !! Scale parameter for k_i momentum to integrate over in core -> free calculation
        !!
        !! maximum initial electron momentum = ki_s * Z * alpha * m_e
        !!
        !! Generally want this to be >> Z alpha m_e, the scale factor of the electron wave functions

    integer :: n_ki = 2
        !! Number of radial points in integration of ki

    real(dp) :: ki_min = 1.0e3_dp
        !! Minimum electron momentum 
        !!
        !! Units : eV

    NAMELIST /numerics/ q_bin_width, &
                        n_q_bins, &
                        E_bin_width, &
                        n_E_bins, &
                        E_bin_threshold, &
                        Ef_max, &
                        n_init, &
                        n_fin, &
                        ! q_s_FFT, &
                        n_FFT_grid_input, &
                        n_kf_theta, &
                        n_kf_phi, &
                        n_ki, &
                        n_ki_theta, &
                        n_ki_phi, &
                        ki_s, &
                        ki_min

contains

    subroutine load_numerics(filename, verbose)

        implicit none

        logical, optional :: verbose

        character(len=*) :: filename

        logical :: file_exists

        integer :: error

        if ( verbose ) then

            print*, 'Loading numerics parameters...'
            print*

        end if

        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            open(100, file = trim(filename), iostat = error)
            read(100, nml=numerics, iostat=error)
            close(100)

            if ( error .ne. 0 ) then

                if ( verbose ) then

                    print*, '!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                    print*
                    print*, '    Problem reading numerics namelist.'
                    print*
                    print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                    print*

                end if

                stop

            end if

            call print_numerics(verbose=verbose)

            if ( verbose ) then

                print*, '----------------------------------------'
                print*

            end if

        else

            if ( verbose ) then

                print*, '!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*
                print*, '    Input file for numerics parameters : ', trim(filename), ' does NOT exist.'
                print*
                print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*

            end if

            stop

        end if

    end subroutine

    subroutine print_numerics(verbose)

        implicit none

        logical, optional :: verbose

        character(len=64) :: n_E_bins_str
        character(len=64) :: n_q_bins_str

        if ( verbose ) then

            write(n_E_bins_str, *) n_E_bins
            write(n_q_bins_str, *) n_q_bins

            print*, '----------------------------------------'
            print*, '    --------'
            print*, '    Numerics'
            print*, '    --------'
            print*
            print*, '        Number of initial states = ', n_init
            print*, '        Number of final states   = ', n_fin
            print*
            print*, '        Number of E bins = ', trim(adjustl(n_E_bins_str)), ' + 1'
            print*, '        Energy bin width = ', E_bin_width, ' eV'
            print*
            print*, '        Number of q bins = ', trim(adjustl(n_q_bins_str)), ' + 1'
            print*, '        q bin width      = ', q_bin_width/1.0e3_dp, ' keV'
            print*
            print*, '        Ef_max = ', Ef_max, ' eV'
            print*

        end if

    end subroutine

    subroutine save_numerics(filename, verbose)

        implicit none
        character(len=*) :: filename

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id

        logical, optional :: verbose
        logical :: file_exists

        integer(HSIZE_T) :: dims1(1) = [1]
        integer(HSIZE_T) :: dims2(2)

        integer :: error

        if ( verbose ) then

            print*, 'Saving numerics...'
            print*

        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'numerics', group_id, error)

            ! ! write data
            call h5ltmake_dataset_int_f(file_id, 'numerics/n_E_bins', size(dims1), dims1,&
                n_E_bins + 1, error)
            call h5ltmake_dataset_int_f(file_id, 'numerics/n_q_bins', size(dims1), dims1,&
                n_q_bins + 1, error)
            call h5ltmake_dataset_int_f(file_id, 'numerics/n_init', size(dims1), dims1,&
                n_init, error)
            call h5ltmake_dataset_int_f(file_id, 'numerics/n_fin', size(dims1), dims1,&
                n_fin, error)
            call h5ltmake_dataset_double_f(file_id, 'numerics/E_bin_width', size(dims1), dims1,&
                E_bin_width, error)
            call h5ltmake_dataset_double_f(file_id, 'numerics/q_bin_width', size(dims1), dims1,&
                q_bin_width, error)
            call h5ltmake_dataset_int_f(file_id, 'numerics/E_bin_threshold', size(dims1), dims1,&
                E_bin_threshold, error)

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
