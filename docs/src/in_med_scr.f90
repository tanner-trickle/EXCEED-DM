module in_med_scr
    !! Handle analytic/numeric dielectric or other in medium effects which screen the interaction
    !! 
    !! rate ~ 1/screen^2

    use hdf5
    use h5lt

    use prec
    use constants

    implicit none

    character(len=64) :: screen_type = ''

    !!!
    !! Eq 6 from 
    !!      https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892
    real(dp) :: di_e0
    real(dp) :: di_q_tf
    real(dp) :: di_omega_p
    real(dp) :: di_alpha
    !!!

    logical :: include_screen = .TRUE.
        !! toggle whether screening effects are included or not

    NAMELIST /in_medium/ screen_type,&
                            di_e0, &
                            di_q_tf, &
                            di_omega_p, &
                            di_alpha, &
                            include_screen

contains

    subroutine print_in_med(verbose)
        implicit none

        logical, optional :: verbose

        if ( verbose ) then

            print*, '----------------------------------------'
            print*, '    ---------'
            print*, '    In-Medium'
            print*, '    ---------'
            print*
            print*, '        Include screen? : ', include_screen
            print*, '        Screen type     : ', trim(screen_type)
            print*
            if ( trim(screen_type) == 'analytic' ) then

                print*, '        e_0      : ', di_e0
                print*, '        alpha    : ', di_alpha
                print*, '        q_TF     : ', di_q_tf/1.0e3_dp, 'keV'
                print*, '        omega_p  : ', di_omega_p, 'eV'
                print*

            end if

        end if

    end subroutine

    subroutine load_in_med_scr(filename, verbose)

        implicit none

        character(len=*) :: filename

        logical :: file_exists

        logical, optional :: verbose

        integer :: error

        if ( verbose ) then

            print*, 'Loading in-medium screening parameters...'
            print*

        end if

        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            open(100, file = trim(filename), iostat = error)
            read(100, nml=in_medium, iostat=error)
            close(100)

            call print_in_med(verbose = verbose)

            if ( verbose ) then

                print*, '----------------------------------------'
                print*

            end if

        else

            if ( verbose ) then

                print*, '!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*
                print*, '    Input file for in medium parameters : ', trim(filename), ' does NOT exist.'
                print*
                print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*

            end if

            stop

        end if

    end subroutine

    subroutine save_in_med_scr(filename, verbose)

        implicit none
        character(len=*) :: filename

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id

        logical, optional :: verbose
        logical :: file_exists

        integer(HSIZE_T) :: dims1(1) = [1]

        integer :: error

        if ( verbose ) then

            print*, 'Saving screening parameters...'
            print*

        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'in_med_screening', group_id, error)

            ! ! write data
            call h5ltmake_dataset_string_f(file_id, 'in_med_screening/screen_type',&
                screen_type, error)
            if ( trim(screen_type) == 'analytic' ) then
                call h5ltmake_dataset_double_f(file_id, 'in_med_screening/di_e0', size(dims1), dims1,&
                    di_e0, error)
                call h5ltmake_dataset_double_f(file_id, 'in_med_screening/di_alpha', size(dims1), dims1,&
                    di_alpha, error)
                call h5ltmake_dataset_double_f(file_id, 'in_med_screening/di_q_tf', size(dims1), dims1,&
                    di_q_tf, error)
                call h5ltmake_dataset_double_f(file_id, 'in_med_screening/di_omega_p', size(dims1), dims1,&
                    di_omega_p, error)
            end if

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
                print*

            end if

            stop

        end if

    end subroutine

    function screening(q_vec, omega) result(scr)

        implicit none

        real(dp) :: q_vec(3)
        real(dp) :: omega

        real(dp) :: scr

        if ( include_screen ) then

            if ( trim(screen_type) .eq. 'analytic' ) then
                
                scr = model_dielectric(q_vec, omega, &
                                        di_e0, &
                                        di_q_tf, &
                                        di_alpha, &
                                        di_omega_p)

            else

                scr = 1.0_dp

            end if

        else 

            scr = 1.0_dp

        end if

    end function

    function model_dielectric(q_vec, omega,&
            e0, q_tf, alpha, omega_p) result(di)
        !! Analytic form of the dielectric function
        !! Eq 6 from 
        !!      https://journals.aps.org/prb/pdf/10.1103/PhysRevB.47.9892
        implicit none

        real(dp) :: di

        real(dp) :: e0, q_tf, alpha, omega_p

        real(dp) :: q_vec(3)
        real(dp) :: q_mag
        real(dp) :: omega

        q_mag = norm2(q_vec)

        di = 1.0_dp + ( (e0 - 1.0_dp)**(-1) + &
            alpha*(q_mag/q_tf)**2 + &
            q_mag**4/(4*m_elec**2*omega_p**2) - &
            (omega/omega_p)**2)**(-1) 

    end function

end module
