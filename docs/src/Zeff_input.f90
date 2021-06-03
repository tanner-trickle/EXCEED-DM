module Zeff_input
    !! Loads the Zeff paramters for the free calculation

    use prec

    use DFT_parameters
    use core_electrons

    implicit none

    character(len=64) :: Zeff_type = 'one'
        !! Specify what Zeff to use
        !!
        !! - 'one' - all Zeff = 1
        !! - 'Eb' - use the binding energy of the (valence) state
        !! - 'Eb_c' - use the binding energy of the (core) state
        !! - 'in' - specified by the input

    real(dp) :: n_Eb = 1.0_dp
        !! n to use when computing Z_eff with the binding energy approximation

    real(dp) :: val_Zeff_in(100) = 1.0_dp

    NAMELIST /Zeff/ Zeff_type, &
                        n_Eb, &
                        val_Zeff_in

contains

    subroutine print_Zeff(verbose)
        implicit none

        logical, optional :: verbose

        if ( verbose ) then
            print*, '----------------------------------------'
            print*, '    ----'
            print*, '    Zeff'
            print*, '    ----'
            print*
            print*, '        Type  : ', trim(Zeff_type)
            print*
        end if

    end subroutine

    subroutine load_Zeff_parameters(filename, verbose)

        implicit none

        character(len=*) :: filename

        logical, optional :: verbose

        logical :: file_exists

        integer :: error

        if ( verbose ) then

            print*, 'Loading Zeff parameters...'
            print*

        end if

        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            open(100, file = trim(filename), iostat = error)
            read(100, nml=Zeff, iostat=error)
            close(100)

            if ( error .ne. 0 ) then

                if ( verbose ) then

                    print*, '!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                    print*
                    print*, '    Problem reading Zeff namelist.'
                    print*
                    print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                    print*

                end if

                stop

            end if

            call print_Zeff(verbose = verbose)

            if ( verbose ) then

                print*, '----------------------------------------'
                print*

            end if

        else

            if ( verbose ) then

                print*, '!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*
                print*, '   Input file for Zeff parameters : ', trim(filename), ' does NOT exist.'
                print*
                print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*

            end if

            stop

        end if

    end subroutine

    function get_val_Z_eff(i, k) result(Zeff)
        !! returns the Z_eff value of the valence state

        implicit none

        real(dp) :: Zeff

        integer :: i, k

        real(dp) :: Eb
            !! binding energy

        if ( trim(Zeff_type) == 'one' ) then 

            Zeff = 1.0_dp

        else if ( trim(Zeff_type) == 'Eb' ) then

            Eb = maxval(energy_bands(:, :n_val)) - energy_bands(k, i)

            Zeff = n_Eb*sqrt(Eb/13.6_dp)

            Zeff = max(Zeff, 1.0_dp)

        else if ( trim(Zeff_type) == 'Eb_c' ) then

            Eb = -core_energy(i)

            n_Eb = core_elec_conf(i, 2)

            Zeff = n_Eb*sqrt(Eb/13.6_dp)

            Zeff = max(Zeff, 1.0_dp)

        else if ( trim(Zeff_type) == 'in' ) then

            Zeff = val_Zeff_in(i)

        end if

    end function

end module
