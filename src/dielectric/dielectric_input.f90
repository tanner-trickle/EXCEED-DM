module dielectric_input
    !! Handles input variables for the calculation of the dielectric.

    use prec

    implicit none

    integer :: di_n_omega_bins = 1
        !! Number of omega bins
    real(dp) :: di_omega_bin_width = 1.0_dp
        !! Width of the omega bins
        !!
        !! Units : eV

    integer :: di_n_q_bins = 1
        !! Number of |q| bins
    real(dp) :: di_q_bin_width = 1.0e3_dp
        !! Width of the q bins.
        !!
        !! Units : eV

    integer :: di_n_q_phi_bins = 1
        !! Number of phi_q bins
    integer :: di_n_q_theta_bins = 1
        !! Number of theta_q bins

    logical :: load_dielectric_from_file = .FALSE.
        !! Specify whether the dielectric matrix should be loaded
        !! from a file or computed.

    real(dp) :: di_width_a = 0.2_dp
        !! Parameter in the model for the width:
        !!
        !! delta = min( a + b*omega, width_max )
        !!
        !! Units: eV

    real(dp) :: di_log_width_b = -3.0_dp
        !! Parameter in the model for the width:
        !!
        !! delta = min( a + b*omega, width_max )
        !!
        !! Units: None

    real(dp) :: di_width_max = 0.2_dp
        !! Parameter in the model for the width:
        !!
        !! delta = min( a + b*omega, width_max )
        !!
        !! Units: eV

    integer :: di_n_init = 1
        !! number of initial bands to compute for.
    integer :: di_n_fin = 1
        !! number of final bands to compute for.

    integer :: n_k_vec(3)

    NAMELIST /dielectric/ di_n_omega_bins, &
                            di_n_q_bins, &
                            di_n_q_phi_bins, &
                            di_n_q_theta_bins, &
                            load_dielectric_from_file, &
                            di_width_a, &
                            di_log_width_b, &
                            di_width_max, &
                            di_n_init, &
                            di_n_fin, &
                            di_q_bin_width, &
                            di_omega_bin_width, &
                            n_k_vec

contains

    function di_width_func(omega) result(delta)
        !! Parameterization of the electron lifetime/width.

        real(dp) :: omega
        real(dp) :: delta

        delta = min( di_width_a + 10.0_dp**(di_log_width_b)*omega, di_width_max )

    end function

    subroutine print_dielectric_input(verbose)
        implicit none

        logical, optional :: verbose

        if ( verbose ) then
            print*, '----------------------------------------'
            print*, '    ----------'
            print*, '    Dielectric'
            print*, '    ----------'
            print*
            print*, '        Load from file? : ', load_dielectric_from_file
            print*
            print*, '        Number of initial states : ', di_n_init
            print*, '        Number of final states   : ', di_n_fin
            print* 
            print*, '        Binning : '
            print*
            print*, '            Number of w bins : ', di_n_omega_bins
            print*, '            w bin width      : ', di_omega_bin_width
            print*, '            Number of q bins : ', di_n_q_bins
            print*, '            q bin width      : ', di_q_bin_width
            print*, '            Number of q theta bins : ', di_n_q_phi_bins
            print*, '            Number of q phi bins : ', di_n_q_theta_bins
            print*
            print*, '        delta = min( ', di_width_max, ', ', di_width_a, ' + ', 10.0_dp**di_log_width_b, ' x w )' 
            print*
        end if

    end subroutine

    subroutine load_dielectric_input(filename, verbose)
        implicit none

        character(len=*) :: filename
        logical :: file_exists

        logical, optional :: verbose

        integer :: error

        if ( verbose ) then
            print*, 'Loading dielectric parameters...'
            print*
        end if

        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            open(100, file = trim(filename), iostat = error)
            read(100, nml=dielectric, iostat=error)
            close(100)

            call print_dielectric_input(verbose = verbose)

            if ( verbose ) then

                print*, '----------------------------------------'
                print*

            end if

        else

            if ( verbose ) then

                print*, '!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*
                print*, '    Input file for dielectric parameters : ', trim(filename), ' does NOT exist.'
                print*
                print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*

            end if

            stop

        end if

    end subroutine

end module
