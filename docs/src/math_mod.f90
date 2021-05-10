module math_mod
    !! A collection of useful math functions

    use prec
    use special_functions
    use constants

    implicit none

contains

    function integrate_power_law(b, x1, x2, x_s) result(integral)
        !!
        !! = int_{x1}^{x2} (x/x_s)^b dx = x1 * int_{1}^{x2/x1} (x1/x_s)^b y^b dy
        !! = (b + 1)^(-1) * x1 * exp( b log(x1/x_s) ) * ( exp( (b + 1) log(x2/x1) ) - 1 )
        !!
        !! The reason this function is non-trivial is because evaluating x**(b + 1)
        !! when b is large can be problematic.

        implicit none

        real(dp) :: b
        real(dp) :: x1, x2, x_s

        real(dp) :: integral

        real(dp) :: eps

        eps = 1.0e-8_dp

        if ( abs(b + 1.0_dp) < eps ) then

            integral = x_s*log(x2/x1)

        else 

            integral = ( b + 1.0_dp )**(-1)*x1*exp( b*log(x1/x_s) )*&
                        ( exp( (b + 1.0_dp)*log(x2/x1) ) - 1.0_dp )

        end if

    end function

    function power_law_fit(log_x_pts, log_y_pts) result(fit_params)
        !! Finds the best fit parameters for y = y1 * (x/x1) ^ b, given log_y and log_x
        !!
        !! fit_params = [a, b]

        implicit none

        real(dp) :: log_x_pts(2)
        real(dp) :: log_y_pts(2)

        real(dp) :: fit_params
        real(dp) :: b

        if ( log_x_pts(2) > log_x_pts(1) ) then

            fit_params = (log_y_pts(2) - log_y_pts(1))/(log_x_pts(2) - log_x_pts(1))

        else

            print*, 'ERROR'

        end if

    end function

    function generate_uniform_points_on_sphere(n_theta, n_phi) result(angular_mesh)
        !! Generates an (n_thata*n_phi, 2) list of points which are uniformly 
        !! distributed on the sphere

        implicit none

        integer :: n_theta, n_phi

        real(dp) :: angular_mesh(n_theta*n_phi, 2)

        integer :: t, p, i
        real(dp) :: theta, phi

        i = 1
        do t = 1, n_theta
            do p = 1, n_phi

                theta = acos(&
                2.0_dp*((t - 1.0_dp)/max(1.0_dp, n_theta - 1.0_dp)) - 1.0_dp&
                )

                phi = 2.0_dp*pi*(p - 1.0_dp)/(max(n_phi - 1.0_dp, 1.0_dp))

                angular_mesh(i, :) = [theta, phi]

                i = i + 1

            end do
        end do

    end function

    function Q_func(x, a, b, Q_max) result(Q)
        !! Places x in an appropriate bin
        !!
        !! Q = min(Q_max, 1 + floor((x - a)/b))
        implicit none

        integer :: Q, Q_max

        real(dp) :: x, a, b

        Q = min( Q_max, 1 + floor((x - a)/b) )

    end function

    function factorial(n) result(fact)
        implicit none
        integer :: n, fact, i

        fact = 1
        if ( n .gt. 1 ) then
            do i = 1, n
                fact = fact*i
            end do
        end if
    end function

    recursive function sph_harmonic(l, m, theta, phi) result(y_lm)
        !! Spherical harmonic function with phase convention identical
        !! to Mathematica
        implicit none

        integer :: l, m

        real(dp) :: theta, phi, c_theta
        real(dp) :: legendre_pol
        real(dp) :: norm

        complex(dp) :: y_lm

        if ( m < 0 ) then

            y_lm = (-1)**m*conjg(sph_harmonic(l, -m, theta, phi)) 

        else

            c_theta = cos(theta)

            norm = sqrt((2*l + 1)*factorial(l - m)/(4.0_dp*pi*factorial(l + m)))

            call lpmv(1.0_dp*l, m, c_theta, legendre_pol)
            y_lm = norm*exp(ii*m*phi)*legendre_pol

        end if

    end function

    function get_phi(n_hat) result(phi)
        !! returns the theta value of a unit direction vector, n_hat
        implicit none

        real(dp) :: n_hat(3)
        real(dp) :: phi

        phi = atan2(n_hat(2), n_hat(1))

        if ( phi .lt. 0.0_dp ) then
            phi = phi + 2.0_dp*pi
        end if

    end function

    function get_theta(n_hat) result(theta)
        !! returns the theta value of a unit direction vector, n_hat
        implicit none

        real(dp) :: n_hat(3)
        real(dp) :: theta

        theta = acos(n_hat(3))

    end function


end module

