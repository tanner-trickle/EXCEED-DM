module particle_physics_abs
    !! Relevant particle physics parameters for an absorption
    !! calculation

    use prec
    use constants
    use units
    use math_mod

    implicit none

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

    real(dp) :: vel_dist_N0
    integer :: n_v_mag, n_v_theta, n_v_phi

    integer :: n_time = 1
    real(dp) :: time_day_min = 0.0_dp
    real(dp) :: time_day_max = 1.0_dp

    real(dp), allocatable :: time_of_day(:)
    real(dp), allocatable :: vE_vec_list(:, :)

    NAMELIST /particle_physics_a/ v0_km_per_sec, &
                                    vE_km_per_sec, &
                                    vEsc_km_per_sec, &
                                    rhoX_GeV_per_cm3, &
                                    thetaE, &
                                    n_time, &
                                    time_day_min, &
                                    time_day_max, &
                                    n_v_mag, &
                                    n_v_theta, &
                                    n_v_phi

contains

    subroutine print_particle_physics_abs(verbose)

        implicit none

        logical, optional :: verbose

        if ( verbose ) then
            print*, 'Dark Matter properties'
            print*
            print*, '   Velocity distribution parameters : '
            print*, '       v0   = ', v0_km_per_sec, ' km/sec'
            print*, '       vE   = ', vE_km_per_sec, ' km/sec'
            print*, '       vEsc = ', vEsc_km_per_sec, ' km/sec'
            print*
            print*, '   Velocity integral parameters : '
            print*, '       n_v_mag = ', n_v_mag
            print*, '       n_v_theta = ', n_v_theta
            print*, '       n_v_phi = ', n_v_phi
            print* 
            print*, '   Times : ', time_of_day(:)
            print* 
            print*, '----------'
            print*
        end if

    end subroutine

    subroutine load_particle_physics_abs(filename, verbose)

        implicit none

        logical, optional :: verbose

        character(len=*) :: filename

        logical :: file_exists

        integer :: error

        integer :: t

        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            open(100, file = trim(filename), iostat = error)

            read(100, nml=particle_physics_a, iostat=error)

            close(100)

            rhoX = inv_cm_to_eV**3*1.0e9_dp*rhoX_GeV_per_cm3
            v0 = km_per_sec_to_none*v0_km_per_sec 
            vE = km_per_sec_to_none*vE_km_per_sec
            vEsc = km_per_sec_to_none*vEsc_km_per_sec

            allocate(time_of_day(n_time))
            allocate(vE_vec_list(n_time, 3))
            
            do t = 1, n_time

                time_of_day(t) = time_day_min + (time_day_max - time_day_min)*&
                                (t - 1.0_dp)/(max(1, n_time - 1))

                vE_vec_list(t, 1) = vE*sin(thetaE)*sin(2.0_dp*PI*time_of_day(t))
                vE_vec_list(t, 2) = vE*cos(thetaE)*sin(thetaE)*(cos(2.0_dp*pi*time_of_day(t)) - 1)
                vE_vec_list(t, 3) = vE*((sin(thetaE)**2)*cos(2.0_dp*pi*time_of_day(t)) &
                                    + cos(thetaE)**2)
            end do

            call print_particle_physics_abs(verbose = verbose)

            vel_dist_N0 = (pi*v0**3)*(sqrt(pi)*erf(vEsc/v0) - 2*(vEsc/v0)*exp(-(vEsc/v0)**2))

            call check_mb_dist_normalization(n_v_mag, n_v_theta, n_v_phi, &
               boost_vec_in = vE*[0, 0, 1], verbose = verbose)

        else

            if ( verbose ) then

                print*, '!! ERROR !!'
                print*
                print*, '   Input file for particle physics (absorption) parameters : ', trim(filename), ' does NOT exist.'
                print*
                print*, '!!!!!!!!!!!'
                print*

            end if

            stop

        end if

    end subroutine

    function mb_vel_distribution(v_vec, boost_vec_in) result(mb_val)
        !! Maxwell boltzmann distribution, boosted with boost_vec
        
        implicit none

        real(dp) :: v_vec(3)
        real(dp) :: v_mag

        real(dp), optional :: boost_vec_in(3)
        real(dp) :: boost_vec(3)

        real(dp) :: mb_val

        real(dp) :: v_p(3), v_p_mag

        ! boost the distribution
        if ( present(boost_vec_in) ) then
            boost_vec = boost_vec_in
        else
            boost_vec = [0.0_dp, 0.0_dp, 0.0_dp]
        end if

        v_p = v_vec + boost_vec
        v_p_mag = norm2(v_p)

        if ( v_p_mag < vEsc ) then

            mb_val = (vel_dist_N0)**(-1)*exp( -(v_p_mag/v0)**2 )

        else

            mb_val = 0.0_dp

        end if

    end function

    subroutine check_mb_dist_normalization(n_v_mag, n_v_theta, n_v_phi, boost_vec_in, verbose)
        !! make sure the mb distribution integrates to 1
        !!
        !! good indicator of whether the mesh is large enough to integrate over

        implicit none

        integer :: v, t, p

        integer :: n_v_mag, n_v_theta, n_v_phi

        real(dp) :: v_mag_list(n_v_mag)

        real(dp) :: v_angular_mesh(n_v_theta*n_v_phi, 2)

        real(dp) :: int_val

        real(dp), optional :: boost_vec_in(3)
        real(dp) :: boost_vec(3)

        real(dp) :: v_vec(3)
        real(dp) :: v_mag, v_theta, v_phi
        real(dp) :: v_max

        logical, optional :: verbose

        if ( present(boost_vec_in) ) then
            boost_vec = boost_vec_in
        else
            boost_vec = [0.0_dp, 0.0_dp, 0.0_dp]
        end if

        v_angular_mesh = generate_uniform_points_on_sphere(n_v_theta, n_v_phi)

        v_max = vEsc + norm2(boost_vec)

        do v = 1, n_v_mag

            v_mag_list(v) = v_max*(v - 1.0_dp)/max(1.0_dp, n_v_mag - 1.0_dp)

        end do

        int_val = 0.0_dp

        do v = 1, n_v_mag
            do t = 1, n_v_theta*n_v_phi

                v_mag = v_mag_list(v)
                v_theta = v_angular_mesh(t, 1)
                v_phi = v_angular_mesh(t, 2)

                v_vec(1) = v_mag*sin(v_theta)*cos(v_phi)
                v_vec(2) = v_mag*sin(v_theta)*sin(v_phi)
                v_vec(3) = v_mag*cos(v_theta)

                int_val = int_val + v_mag**2*(4.0_dp*pi*v_max)*(1.0_dp*n_v_mag*n_v_theta*n_v_phi)**(-1)*&
                    mb_vel_distribution(v_vec, boost_vec_in = boost_vec)

            end do
        end do

        if ( verbose ) then

            print*, '    Integral of MB velocity distrubution : ', int_val
            print*

        end if

    end subroutine

end module
