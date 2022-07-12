module abs_physics_functions

    use prec_util, only: dp
    use constants_util

    implicit none

contains

    function Pi_scaling_func(omega, omega_IF) result ( scaling )
        ! As $\omega \rightarrow 0$, $\Pi \propto \omega^2$ such that $\Gamma \propto \frac{\Pi}{\omega} \rightarrow 0$. 
        ! This function guarantees that scaling behavior.

        implicit none

        real(dp), intent(in) :: omega, omega_IF
        real(dp) :: scaling

        scaling = ( omega / omega_IF )**2

    end function

    function width_func(omega, a, b, c) result ( delta )

        implicit none

        real(dp) :: omega, a, b, c
        real(dp) :: delta

        delta = min( a + b*omega, c )

    end function

    function green_func(omega, omega_IF, delta, type) result (green_val)

        use constants_util

        implicit none

        real(dp) :: omega, omega_IF, delta

        real(dp) :: d_omega_m, d_omega_p

        complex(dp) :: green_val

        real(dp) :: green_val_r, green_val_c

        character(len=*) :: type

        ! green_val = ( ( omega - omega_IF + ii*delta )**(-1) - ( omega + omega_IF - ii*delta )**(-1) )

        d_omega_m = omega - omega_IF
        d_omega_p = omega + omega_IF

        green_val_r =  ( d_omega_m / ( d_omega_m**2 + delta**2 ) ) - ( d_omega_p / ( d_omega_p**2 + delta**2 ) )

        select case ( trim(adjustl(type)) ) 

            case ( 'lorentz' )

                green_val_c = delta / ( d_omega_m**2 + delta**2 )

            case ( 'gauss' )

                green_val_c = sqrt(pi)*delta**(-1)*exp( - ( d_omega_m / delta )**2 )

        end select

        green_val = green_val_r - ii*green_val_c

    end function

end module
