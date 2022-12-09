module absorption_rate_calculator

    use prec_util, only: dp

    use exdm_inputs_type

    use PiIF_calculator_type

    implicit none

contains

    subroutine absorption_rate_compute(exdm_inputs, spin_dof, PiIF_calculator, absorption_rate)

        implicit none

        type(exdm_inputs_t) :: exdm_inputs
        integer :: spin_dof
        type(PiIF_calculator_t) :: PiIF_calculator

        real(dp) :: absorption_rate(:, :)

        if ( exdm_inputs%control%verbose ) then
            print*, 'Computing absorption rate...'
            print*
        end if

        select case ( trim(adjustl(exdm_inputs%dm_model%particle_type)) )

            case ( 'vector' )

                call absorption_rate_compute_vector(exdm_inputs, & 
                    PiIF_calculator%Pi_1_1_mat, absorption_rate)

            case ( 'scalar' )

                call absorption_rate_compute_scalar(exdm_inputs, &
                    PiIF_calculator%Pi_v2_v2, absorption_rate)

            case ( 'ps' )

                if ( spin_dof == 1 ) then

                    call absorption_rate_compute_ps_SI(exdm_inputs, &
                        PiIF_calculator%Pi_v_v, absorption_rate)

                else 

                    call absorption_rate_compute_ps_SD(exdm_inputs, &
                        PiIF_calculator%Pi_vds_vds, absorption_rate)

                end if

            case ( 'mdm' )

                call absorption_rate_compute_mdm(exdm_inputs, &
                    PiIF_calculator%Pi_v_v, absorption_rate)

            case ( 'edm' )

                call absorption_rate_compute_edm(exdm_inputs, &
                    PiIF_calculator%Pi_v2_v2, absorption_rate)

        end select

        if ( exdm_inputs%control%verbose ) then
            print*, 'Done computing absorption rate!'
            print*
        end if

    end subroutine

    subroutine absorption_rate_compute_ps_SI(exdm_inputs, Pi_v_v, absorption_rate)

        use math_util
        use constants_util

        implicit none

        type(exdm_inputs_t) :: exdm_inputs
        complex(dp) :: Pi_v_v(:, :, :, :)

        real(dp) :: absorption_rate(:, :)

        integer :: w, m, j

        real(dp) :: gam

        complex(dp) :: pi_mat(3, 3)

        real(dp) :: omega

        do w = 1, size(absorption_rate, 2)
            do m = 1, size(absorption_rate, 1)

                omega = exdm_inputs%dm_model%mX(m)

                gam = 0.0_dp

                do j = 1, 3

                    gam = gam + (-1.0_dp/omega)*&
                        (omega**2/(4.0_dp*m_elec**2))*&
                        aimag(Pi_v_v(j, j, m, w))

                end do

                absorption_rate(m, w) = (exdm_inputs%dm_model%rho_X/exdm_inputs%material%rho_T)*&
                    exdm_inputs%dm_model%mX(m)**(-1)*&
                    exdm_inputs%experiment%M*exdm_inputs%experiment%T*&
                    gam

            end do
        end do

    end subroutine

    subroutine absorption_rate_compute_ps_SD(exdm_inputs, Pi_vds_vds, absorption_rate)

        use math_util
        use constants_util

        implicit none

        type(exdm_inputs_t) :: exdm_inputs
        complex(dp) :: Pi_vds_vds(:, :)

        real(dp) :: absorption_rate(:, :)

        integer :: w, m

        real(dp) :: gam

        do w = 1, size(absorption_rate, 2)
            do m = 1, size(absorption_rate, 1)

                gam = -(exdm_inputs%dm_model%mX(m))*(4.0_dp*m_elec**2)**(-1)*&
                    aimag(Pi_vds_vds(m, w))

                absorption_rate(m, w) = (exdm_inputs%dm_model%rho_X/exdm_inputs%material%rho_T)*&
                    exdm_inputs%dm_model%mX(m)**(-1)*&
                    exdm_inputs%experiment%M*exdm_inputs%experiment%T*&
                    gam

            end do
        end do

    end subroutine

    subroutine absorption_rate_compute_scalar(exdm_inputs, Pi_v2_v2, absorption_rate)

        use math_util
        use constants_util

        implicit none

        type(exdm_inputs_t) :: exdm_inputs
        complex(dp) :: Pi_v2_v2(:, :)

        real(dp) :: absorption_rate(:, :)

        integer :: w, m

        real(dp) :: gam

        do w = 1, size(absorption_rate, 2)
            do m = 1, size(absorption_rate, 1)

                gam = -(exdm_inputs%dm_model%mX(m))**(-1)*(1.0_dp/4.0_dp)*&
                    aimag(Pi_v2_v2(m, w))

                absorption_rate(m, w) = (exdm_inputs%dm_model%rho_X/exdm_inputs%material%rho_T)*&
                    exdm_inputs%dm_model%mX(m)**(-1)*&
                    exdm_inputs%experiment%M*exdm_inputs%experiment%T*&
                    gam

            end do
        end do

    end subroutine

    subroutine absorption_rate_compute_edm(exdm_inputs, Pi_v2_v2, absorption_rate)

        use math_util
        use constants_util

        implicit none

        type(exdm_inputs_t) :: exdm_inputs
        complex(dp) :: Pi_v2_v2(:, :)

        real(dp) :: absorption_rate(:, :)

        integer :: w, m

        real(dp) :: gam

        absorption_rate = (exdm_inputs%dm_model%rho_X/exdm_inputs%material%rho_T)*&
                    exdm_inputs%experiment%M*exdm_inputs%experiment%T*&
                    aimag( Pi_v2_v2 )

    end subroutine

    subroutine absorption_rate_compute_mdm(exdm_inputs, Pi_v_v, absorption_rate)

        use math_util
        use constants_util

        implicit none

        type(exdm_inputs_t) :: exdm_inputs
        complex(dp) :: Pi_v_v(:, :, :, :)

        real(dp) :: absorption_rate(:, :)

        integer :: w, m, j

        real(dp) :: gam

        complex(dp) :: pi_mat(3, 3)

        real(dp) :: omega

        do w = 1, size(absorption_rate, 2)
            do m = 1, size(absorption_rate, 1)

                omega = exdm_inputs%dm_model%mX(m)

                gam = 0.0_dp

                do j = 1, 3

                    gam = gam + (-2.0_dp/3.0_dp)*aimag(Pi_v_v(j, j, m, w))

                end do

                absorption_rate(m, w) = (exdm_inputs%dm_model%rho_X/exdm_inputs%material%rho_T)*&
                    exdm_inputs%experiment%M*exdm_inputs%experiment%T*&
                    gam

            end do
        end do

    end subroutine

    subroutine absorption_rate_compute_vector(exdm_inputs, Pi_1_1_mat, absorption_rate)

        use math_util
        use constants_util

        implicit none

        type(exdm_inputs_t) :: exdm_inputs
        complex(dp) :: Pi_1_1_mat(:, :, :, :)

        real(dp) :: absorption_rate(:, :)

        integer :: w, m, j

        real(dp) :: gam

        complex(dp) :: pi_eigvals(3)
        complex(dp) :: pi_mat(3, 3)

        do w = 1, size(absorption_rate, 2)
            do m = 1, size(absorption_rate, 1)

                pi_mat = e_EM**2*(exdm_inputs%dm_model%mX(m)/m_elec)**2*Pi_1_1_mat(:, :, m, w)

                pi_eigvals = calc_eigvals_33(pi_mat)

                gam = 0.0_dp

                do j = 1, 3
                    gam = gam + (-1.0_dp)*e_EM**(-2)*(3.0_dp*exdm_inputs%dm_model%mX(m))**(-1)*&
                        aimag( exdm_inputs%dm_model%mX(m)**2*pi_eigvals(j) / (exdm_inputs%dm_model%mX(m)**2 - pi_eigvals(j) ) )
                end do

                absorption_rate(m, w) = (exdm_inputs%dm_model%rho_X/exdm_inputs%material%rho_T)*&
                    exdm_inputs%dm_model%mX(m)**(-1)*&
                    exdm_inputs%experiment%M*exdm_inputs%experiment%T*&
                    gam

            end do
        end do

    end subroutine

end module
