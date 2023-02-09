module derivative_VSH_coefficients
    ! $ \nabla^i ( f(r) Y_{\ell m} ) = \sum_{k = -1}^1 \sum_{q = -1}^1 
    ! \left( c^i(k, q, \ell, m) \frac{df}{dr} + d^i(k, q, \ell, m) \frac{f}{r} \right)
    ! Y_{\ell + k, m + q}$

    use prec_util, only: dp
    use constants_util

    implicit none

contains

    function derivative_VSH_c(i, k, q, l, m) result ( c_i )
        ! Computes the coefficients in the vector spherical harmonic, Y

        implicit none

        integer :: i 
        integer :: k, q
        integer :: l, m

        integer :: l_hat, m_hat

        complex(dp) :: c_i

        select case ( i )
            case ( 1 )
                c_i = derivative_VSH_c_x(k, q, l, m)
            case ( 2 )
                c_i = derivative_VSH_c_y(k, q, l, m)
            case ( 3 )
                c_i = derivative_VSH_c_z(k, q, l, m)
        end select

    end function

    function derivative_VSH_d(i, k, q, l, m) result ( d_i )
        ! Computes the coefficients in the vector spherical harmonic, \Psi

        implicit none

        integer :: i 
        integer :: k, q
        integer :: l, m

        integer :: l_hat, m_hat

        complex(dp) :: d_i

        select case ( i )
            case ( 1 )
                d_i = derivative_VSH_d_x(k, q, l, m)
            case ( 2 )
                d_i = derivative_VSH_d_y(k, q, l, m)
            case ( 3 )
                d_i = derivative_VSH_d_z(k, q, l, m)
        end select

    end function

    function derivative_VSH_c_x(k, q, l, m) result ( c_x )

        implicit none

        integer :: k, q
        integer :: l, m

        complex(dp) :: c_x

        complex(dp) :: c_x_mat(3, 3)

        c_x_mat(1, 1) = -(0.5_dp)*Amm_coeff(l, m)*( (2*l - 1)*(2*l + 1) )**(-0.5_dp) 
        c_x_mat(1, 2) = (0.0_dp, 0.0_dp)
        c_x_mat(1, 3) = (0.5_dp)*Amp_coeff(l, m)*( (2*l - 1)*(2*l + 1) )**(-0.5_dp)

        c_x_mat(2, 1) = ( 0.0_dp, 0.0_dp )
        c_x_mat(2, 2) = ( 0.0_dp, 0.0_dp )
        c_x_mat(2, 3) = ( 0.0_dp, 0.0_dp )

        c_x_mat(3, 1) = (0.5_dp)*Apm_coeff(l, m)*( (2*l + 3)*(2*l + 1) )**(-0.5_dp)
        c_x_mat(3, 2) = ( 0.0_dp, 0.0_dp )
        c_x_mat(3, 3) = -(0.5_dp)*App_coeff(l, m)*( (2*l + 3)*(2*l + 1) )**(-0.5_dp) 

        c_x = c_x_mat(k + 2, q + 2)

    end function

    function derivative_VSH_c_y(k, q, l, m) result ( c_y )

        implicit none

        integer :: k, q
        integer :: l, m

        complex(dp) :: c_y

        complex(dp) :: c_y_mat(3, 3)

        c_y_mat(1, 1) = -ii*(0.5_dp)*Amm_coeff(l, m)*( (2*l - 1)*(2*l + 1) )**(-0.5_dp)
        c_y_mat(1, 2) = ( 0.0_dp, 0.0_dp )
        c_y_mat(1, 3) = -ii*(0.5_dp)*Amp_coeff(l, m)*( (2*l - 1)*(2*l + 1) )**(-0.5_dp) 

        c_y_mat(2, 1) = ( 0.0_dp, 0.0_dp )
        c_y_mat(2, 2) = ( 0.0_dp, 0.0_dp )
        c_y_mat(2, 3) = ( 0.0_dp, 0.0_dp )

        c_y_mat(3, 1) = ii*(0.5_dp)*Apm_coeff(l, m)*( (2*l + 3)*(2*l + 1) )**(-0.5_dp) 
        c_y_mat(3, 2) = ( 0.0_dp, 0.0_dp )
        c_y_mat(3, 3) = ii*(0.5_dp)*App_coeff(l, m)*( (2*l + 3)*(2*l + 1) )**(-0.5_dp) 

        c_y = c_y_mat(k + 2, q + 2)

    end function

    function derivative_VSH_c_z(k, q, l, m) result ( c_z )

        implicit none

        integer :: k, q
        integer :: l, m

        complex(dp) :: c_z

        complex(dp) :: c_z_mat(3, 3)

        c_z_mat(1, 1) = ( 0.0_dp, 0.0_dp )
        c_z_mat(1, 2) = Am0_coeff(l, m)*(2.0_dp*(2*l - 1)*(2*l + 1))**(-0.5_dp)
        c_z_mat(1, 3) = ( 0.0_dp, 0.0_dp )

        c_z_mat(2, 1) = ( 0.0_dp, 0.0_dp )
        c_z_mat(2, 2) = ( 0.0_dp, 0.0_dp )
        c_z_mat(2, 3) = ( 0.0_dp, 0.0_dp )

        c_z_mat(3, 1) = ( 0.0_dp, 0.0_dp )
        c_z_mat(3, 2) = -Ap0_coeff(l, m)*(2.0_dp*(2*l + 3)*(2*l + 1))**(-0.5_dp)
        c_z_mat(3, 3) = ( 0.0_dp, 0.0_dp )

        c_z = c_z_mat(k + 2, q + 2)

    end function

    function derivative_VSH_d_x(k, q, l, m) result ( d_x )

        implicit none

        integer :: k, q
        integer :: l, m

        complex(dp) :: d_x

        complex(dp) :: d_x_mat(3, 3)

        d_x_mat(1, 1) = -(l + 1)*(0.5_dp)*Amm_coeff(l, m)*( (2*l - 1)*(2*l + 1) )**(-0.5_dp) 
        d_x_mat(1, 2) = ( 0.0_dp, 0.0_dp )
        d_x_mat(1, 3) = (l + 1)*(0.5_dp)*Amp_coeff(l, m)*( (2*l - 1)*(2*l + 1) )**(-0.5_dp) 

        d_x_mat(2, 1) = ( 0.0_dp, 0.0_dp )
        d_x_mat(2, 2) = ( 0.0_dp, 0.0_dp )
        d_x_mat(2, 3) = ( 0.0_dp, 0.0_dp )

        d_x_mat(3, 1) = -l*(0.5_dp)*Apm_coeff(l, m)*( (2*l + 3)*(2*l + 1) )**(-0.5_dp) 
        d_x_mat(3, 2) = ( 0.0_dp, 0.0_dp )
        d_x_mat(3, 3) = l*(0.5_dp)*App_coeff(l, m)*( (2*l + 3)*(2*l + 1) )**(-0.5_dp) 

        d_x = d_x_mat(k + 2, q + 2)

    end function

    function derivative_VSH_d_y(k, q, l, m) result ( d_y )

        implicit none

        integer :: k, q
        integer :: l, m

        complex(dp) :: d_y

        complex(dp) :: d_y_mat(3, 3)

        d_y_mat(1, 1) = -ii*(l + 1)*(0.5_dp)*Amm_coeff(l, m)*( (2*l - 1)*(2*l + 1) )**(-0.5_dp) 
        d_y_mat(1, 2) = ( 0.0_dp, 0.0_dp )
        d_y_mat(1, 3) = -ii*(l + 1)*(0.5_dp)*Amp_coeff(l, m)*( (2*l - 1)*(2*l + 1) )**(-0.5_dp) 

        d_y_mat(2, 1) = ( 0.0_dp, 0.0_dp )
        d_y_mat(2, 2) = ( 0.0_dp, 0.0_dp )
        d_y_mat(2, 3) = ( 0.0_dp, 0.0_dp )

        d_y_mat(3, 1) = -ii*l*(0.5_dp)*Apm_coeff(l, m)*( (2*l + 3)*(2*l + 1) )**(-0.5_dp) 
        d_y_mat(3, 2) = ( 0.0_dp, 0.0_dp )
        d_y_mat(3, 3) = -ii*l*(0.5_dp)*App_coeff(l, m)*( (2*l + 3)*(2*l + 1) )**(-0.5_dp) 

        d_y = d_y_mat(k + 2, q + 2)

    end function

    function derivative_VSH_d_z(k, q, l, m) result ( d_z )

        implicit none

        integer :: k, q
        integer :: l, m

        complex(dp) :: d_z

        complex(dp) :: d_z_mat(3, 3)

        d_z_mat(1, 1) = ( 0.0_dp, 0.0_dp )
        d_z_mat(1, 2) = (l + 1)*Am0_coeff(l, m)*(2.0_dp*(2*l - 1)*(2*l + 1))**(-0.5_dp)
        d_z_mat(1, 3) = ( 0.0_dp, 0.0_dp )

        d_z_mat(2, 1) = ( 0.0_dp, 0.0_dp )
        d_z_mat(2, 2) = ( 0.0_dp, 0.0_dp )
        d_z_mat(2, 3) = ( 0.0_dp, 0.0_dp )

        d_z_mat(3, 1) = ( 0.0_dp, 0.0_dp )
        d_z_mat(3, 2) = l*Ap0_coeff(l, m)*(2.0_dp*(2*l + 3)*(2*l + 1))**(-0.5_dp)
        d_z_mat(3, 3) = ( 0.0_dp, 0.0_dp )

        d_z = d_z_mat(k + 2, q + 2)

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
