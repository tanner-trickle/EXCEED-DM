module math_util

    use prec_util, only: dp

    use constants_util

    implicit none

contains

    function STO_radial(r_list, n, N0, Z) result( wf )
        ! Radial part of a Slater Type Orbital (STO) basis function.
        !
        ! Units : eV^(3/2)

        use constants_util

        implicit none

        real(dp) :: n
        real(dp) :: N0, Z
        real(dp) :: r_list(:)

        real(dp) :: wf(size(r_list))

        wf = a0**(-1.5_dp)*N0*(r_list/a0)**(n - 1.0_dp)*exp(-Z*r_list/a0)

    end function

    function pauli_spin_matrix(i) result(mat)
        ! Returns the \(i\)th Pauli spin matrix, \( \sigma^i\). Identity matrix is defined as the
        ! \( 0 \)th component.

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

    function reduced_mass(m1, m2) result ( mu )

        implicit none

        real(dp), intent(in) :: m1
        real(dp), intent(in) :: m2

        real(dp) :: mu

        mu = m1*m2/(m1 + m2)

    end function

    integer elemental function clamp_int(x, min_int, max_int) result( i )

        implicit none

        integer, intent(in) :: x, min_int, max_int

        i = min( max( x, min_int ), max_int )

    end function

    subroutine compute_sph_harmonics(l, m, theta_list, phi_list, sph_harm_list)
        ! Spherical harmonic function, \( Y^m_l \) with phase convention identical to Mathematica.

        use special_functions

        integer :: l, m

        real(dp) :: theta_list(:)
        real(dp) :: phi_list(:)
        complex(dp) :: sph_harm_list(:)
        real(dp) :: ylm_list(size(sph_harm_list))

        real(dp) :: norm

        real(dp) :: cos_theta_list(size(theta_list))

        integer :: i

        norm = sqrt((2*l + 1)*factorial(l - abs(m))/(4.0_dp*pi*factorial(l + abs(m))))

        if ( m < 0 ) then
            norm = (-1)**m*norm
        end if

        cos_theta_list = cos(theta_list)

        do i = 1, size(cos_theta_list)

            call lpmv(1.0_dp*l, abs(m), cos_theta_list(i), ylm_list(i))

        end do

        sph_harm_list = norm*exp(ii*m*phi_list)*ylm_list

    end subroutine

    subroutine cartesian_to_polar_list(x_vec_list, r_list, theta_list, phi_list)
        ! Takes a list of x vectors (N_x, 3) and returns lists of , \( r, \theta, \phi \) 

        use constants_util

        implicit none

        real(dp) :: x_vec_list(:, :)
        real(dp) :: r_list(:)
        real(dp) :: theta_list(:)
        real(dp) :: phi_list(:)

        real(dp) :: x_hat_list(size(x_vec_list, 1), size(x_vec_list, 2))

        integer :: d

        r_list = norm2(x_vec_list, 2) + 1.0e-8_dp

        do d = 1, 3
            x_hat_list(:, d) = x_vec_list(:, d)/r_list
        end do

        theta_list = acos(x_hat_list(:, 3))
        phi_list = mod(atan2( x_hat_list(:, 2), x_hat_list(:, 1) ) + 2.0_dp*pi, 2.0_dp*pi)

    end subroutine

    function cross_product(vec1, vec2) result( cross )
        ! Computes the cross product of two vectors

        implicit none

        real(dp), intent(in) :: vec1(3)
        real(dp), intent(in) :: vec2(3)
        real(dp) :: cross(3)

        cross(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
        cross(2) = -( vec1(1)*vec2(3) - vec1(3)*vec2(1) )
        cross(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)

    end function

    function calc_eigvals_33(mat) result(eigs)
        ! Computes the eigenvalues of a complex, 3x3 matrix.

        complex(dp), intent(in) :: mat(3, 3)
        complex(dp) :: eigs(3)

        complex(dp) :: lwork(6), rwork(6)
        complex(dp) :: dummy_1(1, 1)
        complex(dp) :: dummy_2(1, 1)

        integer :: info

        call zgeev('N', 'N', 3, mat, 3, eigs, dummy_1, 1, dummy_2, 1, lwork, 6, rwork, info) 

    end function

    function factorial(n) result(fact)
        ! \( n! \) function.
        implicit none
        integer :: n, fact, i

        fact = 1
        if ( n > 1 ) then
            do i = 1, n
                fact = fact*i
            end do
        end if
    end function

end module
