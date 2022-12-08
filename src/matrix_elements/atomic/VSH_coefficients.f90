module VSH_coefficients

    use prec_util, only: dp
    use constants_util

    implicit none

contains

    function VSH_coefficients_Y(i, dl, dm, l, m) result ( c_i )
        ! Computes the coefficients in the vector spherical harmonic, Y

        implicit none

        integer :: i 
        integer :: dl, dm
        integer :: l, m

        complex(dp) :: c_i

        if ( i == 1 ) then

            c_i = VSH_coefficients_Y_x(dl, dm, l, m)

        else if ( i == 2 ) then

            c_i = VSH_coefficients_Y_y(dl, dm, l, m)

        else if ( i == 3 ) then

            c_i = VSH_coefficients_Y_z(dl, dm, l, m)

        end if

    end function

    function VSH_coefficients_Psi(i, dl, dm, l, m) result ( d_i )
        ! Computes the coefficients in the vector spherical harmonic, \Psi

        implicit none

        integer :: i 
        integer :: dl, dm
        integer :: l, m

        complex(dp) :: d_i

        if ( i == 1 ) then

            d_i = VSH_coefficients_Psi_x(dl, dm, l, m)

        else if ( i == 2 ) then

            d_i = VSH_coefficients_Psi_y(dl, dm, l, m)

        else if ( i == 3 ) then

            d_i = VSH_coefficients_Psi_z(dl, dm, l, m)

        end if

    end function

    function VSH_coefficients_Y_x(dl, dm, l, m) result ( c_x )

        implicit none

        integer :: dl, dm
        integer :: l, m

        complex(dp) :: c_x

        complex(dp) :: c_x_mat(3, 3)

        if ( ( l - 1 >= 0 ) .and. ( m + l - 2 >= 0 ) ) then
            c_x_mat(1, 1) = -(0.5_dp)*Amm_coeff(l, m)*( (2*l - 1)*(2*l + 1) )**(-0.5_dp) 
        else
            c_x_mat(1, 1) = ( 0.0_dp, 0.0_dp )
        end if

        c_x_mat(1, 2) = ( 0.0_dp, 0.0_dp )

        if ( ( l - 1 >= 0 ) .and. ( l - m - 2 >= 0 ) ) then
            c_x_mat(1, 3) = (0.5_dp)*Amp_coeff(l, m)*( (2*l - 1)*(2*l + 1) )**(-0.5_dp) 
        else
            c_x_mat(1, 3) = ( 0.0_dp, 0.0_dp )
        end if

        c_x_mat(2, 1) = ( 0.0_dp, 0.0_dp )
        c_x_mat(2, 2) = ( 0.0_dp, 0.0_dp )
        c_x_mat(2, 3) = ( 0.0_dp, 0.0_dp )

        if ( l - m + 2 >= 0 ) then
            c_x_mat(3, 1) = (0.5_dp)*Apm_coeff(l, m)*( (2*l + 3)*(2*l + 1) )**(-0.5_dp) 
        else
            c_x_mat(3, 1) = ( 0.0_dp, 0.0_dp )
        end if

        c_x_mat(3, 2) = ( 0.0_dp, 0.0_dp )

        if ( l + m + 2 >= 0 ) then
            c_x_mat(3, 3) = -(0.5_dp)*App_coeff(l, m)*( (2*l + 3)*(2*l + 1) )**(-0.5_dp) 
        else
            c_x_mat(3, 3) = ( 0.0_dp, 0.0_dp )
        end if

        c_x = c_x_mat(dl + 2, dm + 2)

    end function

    function VSH_coefficients_Y_y(dl, dm, l, m) result ( c_y )

        implicit none

        integer :: dl, dm
        integer :: l, m

        complex(dp) :: c_y

        complex(dp) :: c_y_mat(3, 3)

        if ( ( l - 1 >= 0 ) .and. ( m + l - 2 >= 0 ) ) then
            c_y_mat(1, 1) = -ii*(0.5_dp)*Amm_coeff(l, m)*( (2*l - 1)*(2*l + 1) )**(-0.5_dp) 
        else
            c_y_mat(1, 1) = ( 0.0_dp, 0.0_dp )
        end if

        c_y_mat(1, 2) = ( 0.0_dp, 0.0_dp )

        if ( ( l - 1 >= 0 ) .and. ( l - m - 2 >= 0 ) ) then
            c_y_mat(1, 3) = -ii*(0.5_dp)*Amp_coeff(l, m)*( (2*l - 1)*(2*l + 1) )**(-0.5_dp) 
        else
            c_y_mat(1, 3) = ( 0.0_dp, 0.0_dp )
        end if

        c_y_mat(2, 1) = ( 0.0_dp, 0.0_dp )
        c_y_mat(2, 2) = ( 0.0_dp, 0.0_dp )
        c_y_mat(2, 3) = ( 0.0_dp, 0.0_dp )

        if ( l - m + 2 >= 0 ) then
            c_y_mat(3, 1) = ii*(0.5_dp)*Apm_coeff(l, m)*( (2*l + 3)*(2*l + 1) )**(-0.5_dp) 
        else
            c_y_mat(3, 1) = ( 0.0_dp, 0.0_dp )
        end if

        c_y_mat(3, 2) = ( 0.0_dp, 0.0_dp )

        if ( l + m + 2 >= 0 ) then
            c_y_mat(3, 3) = ii*(0.5_dp)*App_coeff(l, m)*( (2*l + 3)*(2*l + 1) )**(-0.5_dp) 
        else
            c_y_mat(3, 3) = ( 0.0_dp, 0.0_dp )
        end if

        c_y = c_y_mat(dl + 2, dm + 2)

    end function

    function VSH_coefficients_Y_z(dl, dm, l, m) result ( c_z )

        implicit none

        integer :: dl, dm
        integer :: l, m

        complex(dp) :: c_z

        complex(dp) :: c_z_mat(3, 3)

        c_z_mat(1, 1) = ( 0.0_dp, 0.0_dp )

        if ( l - 1 >= 0 ) then
            c_z_mat(1, 2) = Am0_coeff(l, m)*(2.0_dp*(2*l - 1)*(2*l + 1))**(-0.5_dp)
        else 
            c_z_mat(1, 2) = ( 0.0_dp, 0.0_dp )
        end if

        c_z_mat(1, 3) = ( 0.0_dp, 0.0_dp )

        c_z_mat(2, 1) = ( 0.0_dp, 0.0_dp )
        c_z_mat(2, 2) = ( 0.0_dp, 0.0_dp )
        c_z_mat(2, 3) = ( 0.0_dp, 0.0_dp )

        c_z_mat(3, 1) = ( 0.0_dp, 0.0_dp )
        c_z_mat(3, 2) = -Ap0_coeff(l, m)*(2.0_dp*(2*l + 3)*(2*l + 1))**(-0.5_dp)
        c_z_mat(3, 3) = ( 0.0_dp, 0.0_dp )

        c_z = c_z_mat(dl + 2, dm + 2)

    end function

    function VSH_coefficients_Psi_x(dl, dm, l, m) result ( d_x )

        implicit none

        integer :: dl, dm
        integer :: l, m

        complex(dp) :: d_x

        complex(dp) :: d_x_mat(3, 3)

        if ( ( l - 1 >= 0 ) .and. ( m + l - 2 >= 0 ) ) then
            d_x_mat(1, 1) = -(l + 1)*(0.5_dp)*Amm_coeff(l, m)*( (2*l - 1)*(2*l + 1) )**(-0.5_dp) 
        else
            d_x_mat(1, 1) = ( 0.0_dp, 0.0_dp )
        end if

        d_x_mat(1, 2) = ( 0.0_dp, 0.0_dp )

        if ( ( l - 1 >= 0 ) .and. ( l - m - 2 >= 0 ) ) then
            d_x_mat(1, 3) = (l + 1)*(0.5_dp)*Amp_coeff(l, m)*( (2*l - 1)*(2*l + 1) )**(-0.5_dp) 
        else
            d_x_mat(1, 3) = ( 0.0_dp, 0.0_dp )
        end if

        d_x_mat(2, 1) = ( 0.0_dp, 0.0_dp )
        d_x_mat(2, 2) = ( 0.0_dp, 0.0_dp )
        d_x_mat(2, 3) = ( 0.0_dp, 0.0_dp )

        if ( l - m + 2 >= 0 ) then
            d_x_mat(3, 1) = -l*(0.5_dp)*Apm_coeff(l, m)*( (2*l + 3)*(2*l + 1) )**(-0.5_dp) 
        else
            d_x_mat(3, 1) = ( 0.0_dp, 0.0_dp )
        end if

        d_x_mat(3, 2) = ( 0.0_dp, 0.0_dp )

        if ( l + m + 2 >= 0 ) then
            d_x_mat(3, 3) = l*(0.5_dp)*App_coeff(l, m)*( (2*l + 3)*(2*l + 1) )**(-0.5_dp) 
        else
            d_x_mat(3, 3) = ( 0.0_dp, 0.0_dp )
        end if

        d_x = d_x_mat(dl + 2, dm + 2)

    end function

    function VSH_coefficients_Psi_y(dl, dm, l, m) result ( d_y )

        implicit none

        integer :: dl, dm
        integer :: l, m

        complex(dp) :: d_y

        complex(dp) :: d_y_mat(3, 3)

        if ( ( l - 1 >= 0 ) .and. ( m + l - 2 >= 0 ) ) then
            d_y_mat(1, 1) = -ii*(l + 1)*(0.5_dp)*Amm_coeff(l, m)*( (2*l - 1)*(2*l + 1) )**(-0.5_dp) 
        else
            d_y_mat(1, 1) = ( 0.0_dp, 0.0_dp )
        end if

        d_y_mat(1, 2) = ( 0.0_dp, 0.0_dp )

        if ( ( l - 1 >= 0 ) .and. ( l - m - 2 >= 0 ) ) then
            d_y_mat(1, 3) = -ii*(l + 1)*(0.5_dp)*Amp_coeff(l, m)*( (2*l - 1)*(2*l + 1) )**(-0.5_dp) 
        else
            d_y_mat(1, 3) = ( 0.0_dp, 0.0_dp )
        end if

        d_y_mat(2, 1) = ( 0.0_dp, 0.0_dp )
        d_y_mat(2, 2) = ( 0.0_dp, 0.0_dp )
        d_y_mat(2, 3) = ( 0.0_dp, 0.0_dp )

        if ( l - m + 2 >= 0 ) then
            d_y_mat(3, 1) = -ii*l*(0.5_dp)*Apm_coeff(l, m)*( (2*l + 3)*(2*l + 1) )**(-0.5_dp) 
        else
            d_y_mat(3, 1) = ( 0.0_dp, 0.0_dp )
        end if

        d_y_mat(3, 2) = ( 0.0_dp, 0.0_dp )

        if ( l + m + 2 >= 0 ) then
            d_y_mat(3, 3) = -ii*l*(0.5_dp)*App_coeff(l, m)*( (2*l + 3)*(2*l + 1) )**(-0.5_dp) 
        else
            d_y_mat(3, 3) = ( 0.0_dp, 0.0_dp )
        end if

        d_y = d_y_mat(dl + 2, dm + 2)

    end function

    function VSH_coefficients_Psi_z(dl, dm, l, m) result ( d_z )

        implicit none

        integer :: dl, dm
        integer :: l, m

        complex(dp) :: d_z

        complex(dp) :: d_z_mat(3, 3)

        d_z_mat(1, 1) = ( 0.0_dp, 0.0_dp )

        if ( l - 1 >= 0 ) then
            d_z_mat(1, 2) = (l + 1)*Am0_coeff(l, m)*(2.0_dp*(2*l - 1)*(2*l + 1))**(-0.5_dp)
        else 
            d_z_mat(1, 2) = ( 0.0_dp, 0.0_dp )
        end if

        d_z_mat(1, 3) = ( 0.0_dp, 0.0_dp )

        d_z_mat(2, 1) = ( 0.0_dp, 0.0_dp )
        d_z_mat(2, 2) = ( 0.0_dp, 0.0_dp )
        d_z_mat(2, 3) = ( 0.0_dp, 0.0_dp )

        d_z_mat(3, 1) = ( 0.0_dp, 0.0_dp )
        d_z_mat(3, 2) = l*Ap0_coeff(l, m)*(2.0_dp*(2*l + 3)*(2*l + 1))**(-0.5_dp)
        d_z_mat(3, 3) = ( 0.0_dp, 0.0_dp )

        d_z = d_z_mat(dl + 2, dm + 2)

    end function

    function App_coeff(l, m) result ( App )

        implicit none

        integer :: l, m
        real(dp) :: App

        App = sqrt( 1.0_dp*(l + m + 1) * (l + m + 2) )

    end function

    function Ap0_coeff(l, m) result ( Ap0 )

        implicit none

        integer :: l, m
        real(dp) :: Ap0

        Ap0 = -sqrt( 2.0_dp * (l + m + 1) * (l - m + 1) )

    end function

    function Apm_coeff(l, m) result ( Apm )

        implicit none

        integer :: l, m
        real(dp) :: Apm

        Apm = sqrt( 1.0_dp*(l - m + 1) * (l - m + 2) )

    end function

    function Amp_coeff(l, m) result ( Amp )

        implicit none

        integer :: l, m
        real(dp) :: Amp

        Amp = sqrt( 1.0_dp*(l - m - 1) * (l - m) )

    end function

    function Am0_coeff(l, m) result ( Am0 )

        implicit none

        integer :: l, m
        real(dp) :: Am0

        Am0 = sqrt( 2.0_dp * (l + m) * (l - m) )

    end function

    function Amm_coeff(l, m) result ( Amm )

        implicit none

        integer :: l, m
        real(dp) :: Amm

        Amm = sqrt( 1.0_dp*(l + m - 1) * (l + m) )

    end function

end module
