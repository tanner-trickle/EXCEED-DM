module math_mod
    !! Collection of math functions.

    use prec
    use special_functions
    use constants

    implicit none

contains

    function Q_func(x, a, b, Q_max) result(Q)
        !! \( Q = \text{min} \left( Q_\text{max}, 1 + \lfloor \frac{x - a}{b} \rfloor \right) \)
        !!
        !! Note : \( 1 \leq Q \leq Q_\text{max} \).
        implicit none

        integer :: Q, Q_max

        real(dp) :: x, a, b

        Q = max(min( Q_max, 1 + floor((x - a)/b) ), 1)

    end function

    function uniform_list(N, x1, x2) result(li)
        !! Returns a list of `N` uniformly distributed numbers between `x1` and `x2`,
        !! including `x1` and `x2`.
        !!
        !! Note : `N` \( \geq  1 \)

        integer :: N
        real(dp) :: x1
        real(dp) :: x2

        real(dp) :: li(N)

        integer :: i

        do i  = 1, N
            li(i) = x1 + (x2 - x1)*(i - 1)/max(N - 1, 1)
        end do 

    end function

    subroutine calc_eigvals_33(mat, eigval)
        !! Compute the eigenvalues of a 3x3 complex matrix.
        implicit none

        complex(dp) :: mat(3, 3)

        complex(dp) :: eig_vec(3, 3)
        complex(dp) :: eigval(3)

        complex(dp) :: lwork(6), rwork(6)
        complex(dp) :: dummy_1(1, 1)
        complex(dp) :: dummy_2(1, 1)

        integer :: info

        call zgeev('N', 'N', 3, mat, 3, eigval, dummy_1, 1, dummy_2, 1, lwork, 6, rwork, info) 

    end subroutine

    subroutine calc_eig_system_33(mat, eig_val, eig_vec)
        !! Compute the eigenvalues and (right) eigenvectors of a complex 3x3 matrix.

        implicit none

        complex(dp) :: mat(3, 3)

        complex(dp) :: eig_vec(3, 3)
        complex(dp) :: eig_val(3)

        complex(dp) :: lwork(6), rwork(6)
        complex(dp) :: dummy(1, 1)

        integer :: info

        call zgeev('N', 'V', 3, mat, 3, eig_val, dummy, 1, eig_vec, 3, lwork, 6, rwork, info) 

    end subroutine

    function pauli_spin_matrix(i) result(mat)
        !! Returns the \(i\)th Pauli spin matrix, \( \sigma^i\). Identity matrix is defined as the
        !! \( 0 \)th component.

        implicit none

        integer :: i

        complex(dp) :: mat(2, 2)

        mat = (0.0_dp, 0.0_dp)

        select case ( i )

            case ( 0 )

                mat(1, 1) = (1.0_dp, 0.0_dp)
                mat(1, 2) = (0.0_dp, 0.0_dp)
                mat(2, 1) = (0.0_dp, 0.0_dp)
                mat(2, 2) = (1.0_dp, 0.0_dp)

            case ( 1 )

                mat(1, 1) = (0.0_dp, 0.0_dp)
                mat(1, 2) = (1.0_dp, 0.0_dp)
                mat(2, 1) = (1.0_dp, 0.0_dp)
                mat(2, 2) = (0.0_dp, 0.0_dp)

            case ( 2 )

                mat(1, 1) = (0.0_dp, 0.0_dp)
                mat(1, 2) = -ii
                mat(2, 1) = ii
                mat(2, 2) = (0.0_dp, 0.0_dp)
                
            case ( 3 )

                mat(1, 1) = (1.0_dp, 0.0_dp)
                mat(1, 2) = (0.0_dp, 0.0_dp)
                mat(2, 1) = (0.0_dp, 0.0_dp)
                mat(2, 2) = (-1.0_dp, 0.0_dp)

            case default

                mat(1, 1) = (1.0_dp, 0.0_dp)
                mat(1, 2) = (0.0_dp, 0.0_dp)
                mat(2, 1) = (0.0_dp, 0.0_dp)
                mat(2, 2) = (1.0_dp, 0.0_dp)

        end select

    end function

    function integrate_power_law(b, x1, x2, x_s) result(integral)
        !! Integral of a power law function. 
        
        !! $$\begin{align*}
        !!  & = \int_{x_1}^{x_2} \left( \frac{x}{x_s} \right)^b dx \\
        !!  & = x_1 \left( \frac{x_1}{x_s} \right)^b \int_{1}^{x_2/x_1} y^b dy \\
        !!  & = \frac{x_1}{b + 1} \exp( b \log(x_1/x_s) ) ( \exp( (b + 1) \log(x_2/x_1) ) - 1 )
        !! \end{align*}$$
        !!
        !! The complicated implementation is due to problems when \( b \) is large.

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
        !! Finds \( \beta \) such that \( y_2 = y_1 \left( \frac{x_2}{x_1} \right)^\beta \).

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
        !! Generates an (`n_theta` \( \times \) `n_phi`) list of \( \theta, \phi \) coordinates
        !! which are uniformly distributed on the sphere.

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

    function factorial(n) result(fact)
        !! \( n! \) function.
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
        !! Spherical harmonic function, \( Y^m_l \) with phase convention identical to Mathematica.
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
        !! Returns the \( \phi \) value of a unit direction vector, `n_hat`.
        implicit none

        real(dp) :: n_hat(3)
        real(dp) :: phi

        phi = atan2(n_hat(2), n_hat(1))

        if ( phi .lt. 0.0_dp ) then
            phi = phi + 2.0_dp*pi
        end if

    end function

    function get_theta(n_hat) result(theta)
        !! Returns the \( \theta \) value of a unit direction vector, `n_hat`.
        implicit none

        real(dp) :: n_hat(3)
        real(dp) :: theta

        theta = acos(n_hat(3))

    end function

end module

