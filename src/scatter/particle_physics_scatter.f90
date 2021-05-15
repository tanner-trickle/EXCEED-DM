module particle_physics_scatter
    !! Particle physics parameters needed for a DM-electron
    !! scattering rate calculation

    use hdf5
    use h5lt

    use prec
    use constants
    use units

    implicit none

    integer :: n_mX = 1
        !! Number of masses to compute for

    real(dp) :: log_mmin = 9.0_dp
    real(dp) :: log_mmax = 9.0_dp
    
    integer :: n_extra_mX = 0
        !! Optional
        !!
        !! User can specify extra mass points to add in addition to 
        !! the log-uniform ones chosen

    integer :: n_FDM = 1
        !! Number of mediator types to compute for

    real(dp) :: FDMPower_min = 0.0_dp
    real(dp) :: FDMPower_max = 0.0_dp

    integer :: n_time = 1
        !! Number of time of days to compute for
    real(dp) :: time_day_min = 0.0_dp
    real(dp) :: time_day_max = 0.0_dp

    real(dp) :: percentile_cut = 3.0_dp

    real(dp) :: rhoX_GeV_per_cm3 = 0.4_dp
        !! Dark matter density
        !!
        !! Units : GeV/cm^3
    real(dp) :: v0_km_per_sec = 230.0_dp
    real(dp) :: vE_km_per_sec = 240.0_dp
    real(dp) :: vEsc_km_per_sec = 600.0_dp

    real(dp) :: thetaE = (42.0_dp/180.0_dp)*pi !rad

    !! Generated !!!!!!!

    real(dp) :: rhoX
        !! Dark matter density
        !!
        !! Units : eV^4
    
    real(dp) :: v0
    real(dp) :: vE
    real(dp) :: vEsc

    real(dp), allocatable :: mX(:)
    real(dp), allocatable :: mX_2(:)
    real(dp), allocatable :: FDMPowerList(:)
    real(dp), allocatable :: timeOfDayList(:)
    real(dp), allocatable :: vEVecList(:, :)
        !! Dim : [n_times, 3]

    real(dp) :: g_func_N0
    real(dp) :: g_func_c1, g_func_c2

    real(dp) :: v_max
        !! = vE + vEsc

    NAMELIST /particle_physics/ n_mX, &
                                log_mmin, &
                                log_mmax, &
                                n_extra_mX, &
                                n_FDM, &
                                FDMPower_min, &
                                FDMPower_max, &
                                n_time, &
                                time_day_min, &
                                time_day_max, &
                                percentile_cut, &
                                rhoX_GeV_per_cm3, &
                                v0_km_per_sec, &
                                vE_km_per_sec, &
                                vEsc_km_per_sec

    NAMELIST /extra_mX/ mX_2

contains

    subroutine print_particle_physics_scatter(verbose)

        implicit none

        logical, optional :: verbose

        if ( verbose ) then
            print*, '----------------------------------------'
            print*, '    -----------'
            print*, '    Dark Matter'
            print*, '    -----------'
            print*
            print*, '        Density : ', rhoX_GeV_per_cm3, ' GeV/cm^3'
            print*, '        Masses : ', mX
            print* 
            print*, '        Mediator Form Factors (-d log F_DM / d log q) : ', FDMPowerList
            print* 
            print*, '        Halo Velocity Distribution Parameters : '
            print*, '            v0   = ', v0_km_per_sec, ' km/sec'
            print*, '            vE   = ', vE_km_per_sec, ' km/sec'
            print*, '            vEsc = ', vEsc_km_per_sec, ' km/sec'
            print* 
            print*, '        Time of day : ', timeOfDayList
            print* 
        end if

    end subroutine

    subroutine load_particle_physics_scatter(filename, verbose)

        implicit none

        logical, optional :: verbose

        character(len=*) :: filename

        logical :: file_exists

        integer :: error

        if ( verbose ) then

            print*, 'Loading particle physics parameters...'
            print*

        end if

        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            open(100, file = trim(filename), iostat = error)

            read(100, nml=particle_physics, iostat=error)
            rewind(100)

            if ( error .ne. 0 ) then

                if ( verbose ) then

                    print*, '!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                    print*
                    print*, '    Problem reading particle physics namelist.'
                    print*
                    print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                    print*

                end if

                stop

            end if

            rhoX = inv_cm_to_eV**3*1.0e9_dp*rhoX_GeV_per_cm3
            v0 = km_per_sec_to_none*v0_km_per_sec 
            vE = km_per_sec_to_none*vE_km_per_sec
            vEsc = km_per_sec_to_none*vEsc_km_per_sec

            v_max = vE + vEsc

            if ( n_extra_mX .ne. 0 ) then
                allocate(mX_2(n_extra_mX))

                mX_2 = 1.0e9_dp

                read(100, nml=extra_mX, iostat=error)
                rewind(100)
            end if

            close(100)

            call set_g_func_parameters()
            call set_mX(verbose)
            call set_FDM_powers(verbose)
            call set_time_vE_vec(verbose)

            call print_particle_physics_scatter(verbose=verbose)

            if ( verbose ) then

                print*, '----------------------------------------'
                print*

            end if

        else

            if ( verbose ) then

                print*, '!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*
                print*, '    Input file for particle physics parameters : ', trim(filename), ' does NOT exist.'
                print*
                print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*

            end if

            stop

        end if

    end subroutine

    subroutine save_particle_physics_scatter(filename, verbose)

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

            print*, 'Saving particle physics parameters...'
            print*

        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'particle_physics', group_id, error)

            ! ! write data
            call h5ltmake_dataset_int_f(file_id, 'particle_physics/n_mX', size(dims1), dims1,&
                n_mX, error)
            call h5ltmake_dataset_int_f(file_id, 'particle_physics/n_FDM', size(dims1), dims1,&
                n_FDM, error)
            call h5ltmake_dataset_int_f(file_id, 'particle_physics/n_time', size(dims1), dims1,&
                n_time, error)
            dims1 = [n_mX]
            call h5ltmake_dataset_double_f(file_id, 'particle_physics/mX', size(dims1), dims1,&
                mX, error)
            dims1 = [n_FDM]
            call h5ltmake_dataset_double_f(file_id, 'particle_physics/FDM_list', size(dims1), dims1,&
                FDMPowerList, error)
            dims1 = [n_time]
            call h5ltmake_dataset_double_f(file_id, 'particle_physics/t_list', size(dims1), dims1,&
                timeOfDayList, error)
            dims2 = [n_time, 3]
            call h5ltmake_dataset_double_f(file_id, 'particle_physics/vE_list', size(dims2), dims2,&
                vEVecList, error)
            call h5ltmake_dataset_double_f(file_id, 'particle_physics/v0', size(dims1), dims1,&
                v0, error)
            call h5ltmake_dataset_double_f(file_id, 'particle_physics/vEsc', size(dims1), dims1,&
                vEsc, error)
            call h5ltmake_dataset_double_f(file_id, 'particle_physics/vE', size(dims1), dims1,&
                vE, error)

            call h5fclose_f(file_id, error)
            call h5close_f(error)

            ! if ( verbose ) then
            !     print*, '----------'
            !     print*
            ! end if

        else

            if ( verbose ) then

                print*, '!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*
                print*, '   Output file : ', trim(filename), ' does NOT exist.'
                print*
                print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*

            end if

            stop

        end if

    end subroutine

    subroutine set_mX(verbose)
        implicit none

        logical :: verbose

        integer :: m

        n_mX = n_mX + n_extra_mX

        allocate(mX(n_mX))

        do m = 1, n_mX - n_extra_mX
            mX(m) = 10.0_dp**(log_mmin + (log_mmax - log_mmin)*&
                        (m - 1)/(max(1, n_mX - n_extra_mX - 1))&
                        )
        end do

        if ( n_extra_mX .ne. 0) then 

            mX(n_mX - n_extra_mX + 1:) = mX_2

        end if 

    end subroutine

    subroutine set_g_func_parameters()

        implicit none

        ! precomputed constants for the g function defined in the physics module
        g_func_N0 = (pi*v0**3)*(sqrt(pi)*erf(vEsc/v0) - 2*(vEsc/v0)*exp(-(vEsc/v0)**2))
        g_func_c1 = (2*pi**2*v0**2/g_func_N0)
        g_func_c2 = exp(-(vEsc/v0)**2)

    end subroutine

    subroutine set_FDM_powers(verbose)
        implicit none

        logical :: verbose
        integer :: f

        allocate(FDMPowerList(n_FDM))

        do f = 1, n_FDM
            FDMPowerList(f) = FDMPower_min + & 
                (FDMPower_max - FDMPower_min)*(f - 1)/(max(1, n_FDM - 1))
        end do

    end subroutine

    subroutine set_time_vE_vec(verbose)
        implicit none

        logical :: verbose
        integer :: t

        allocate(timeOfDayList(n_time))
        allocate(vEVecList(n_time, 3))
        
        do t = 1, n_time

            timeOfDayList(t) = time_day_min + (time_day_max - time_day_min)*&
                            (t - 1.0_dp)/(max(1, n_time - 1))

            vEVecList(t, 1) = vE*sin(thetaE)*sin(2.0_dp*PI*timeOfDayList(t))
            vEVecList(t, 2) = vE*cos(thetaE)*sin(thetaE)*(cos(2.0_dp*PI*timeOfDayList(t)) - 1)
            vEVecList(t, 3) = vE*((sin(thetaE)**2)*cos(2.0_dp*PI*timeOfDayList(t)) &
                                + cos(thetaE)**2)
        end do

    end subroutine

    function v_minus(q_vec, mX, vE_vec, omega) result(v_m)
        !! v_- function
        !! 
        !! v_- = (1/q)| q_vec . vE + q^2/mX + w |
        !!
        !! Note that we will explicitly check that this value is < v_Esc so that 
        !! the g function is always > 0
        !! 
        !! Units : None

        implicit none

        real(dp) :: mX, omega
        real(dp) :: q_mag

        real(dp) :: v_m

        real(dp) :: q_vec(3)
        real(dp) :: vE_vec(3)

        q_mag = norm2(q_vec)
        v_m = (1/q_mag)*abs(dot_product(q_vec, vE_vec) + 0.5_dp*q_mag**2/mX + omega)

    end function

    function g_func(q, v_m) result (g_fun)
        !! Kinematic function : 
        !!
        !! g_func = 2*pi*int d^3v f_chi(v) delta(w - w_q)
        !!
        !! Units : eV^(-1)
        implicit none

        real(dp) :: q, v_m
        real(dp) :: g_fun

        g_fun = (g_func_c1/q)*(exp(-(v_m/v0)**2) - g_func_c2)

    end function

    function red_mass(m1, m2) result(mu)
        !! Reduced mass 
        implicit none

        real(dp) :: m1, m2, mu 

        mu = m1*m2/(m1 + m2)
    end function

    function F_med_sq_func(q_mag, power) result (F_med_sq_val)
        !! Mediator for factor squared
        !!
        !! Units : None
        implicit none

        real(dp) :: q_mag
        real(dp) :: F_med_sq_val
        real(dp) :: power

        f_med_sq_val = (alpha_EM*m_elec/q_mag)**(2*power)

    end function

end module
