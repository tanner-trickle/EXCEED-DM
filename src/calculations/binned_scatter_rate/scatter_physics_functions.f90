module scatter_physics_functions

    use prec_util, only: dp

    implicit none

contains

    function fermi_factor(Zeff, Ef) result( ff )

        use constants_util

        implicit none

        real(dp), intent(in) :: Zeff
        real(dp), intent(in) :: Ef

        real(dp) :: ff

        real(dp) :: fermi_val

        fermi_val = 2.0_dp*pi*Zeff*&
            (alpha_EM*m_elec/sqrt(2.0_dp*m_elec*Ef))

        ff = fermi_val*( 1.0_dp - exp(-fermi_val) )**(-1)

    end function

    function compute_v_minus(q_vE_list, half_q2_mag_list, qm1_mag_list, &
            mX_inv, omega, v_esc) result(v_m_list)

        implicit none

        real(dp), intent(in) :: q_vE_list(:)
        real(dp), intent(in) :: half_q2_mag_list(:)
        real(dp), intent(in) :: qm1_mag_list(:)
        real(dp), intent(in) :: mX_inv
        real(dp), intent(in) :: omega
        real(dp), intent(in) :: v_esc

        real(dp) :: v_m_list(size(half_q2_mag_list))

        v_m_list = qm1_mag_list*abs(q_vE_list + half_q2_mag_list*mX_inv + omega)

    end function

end module
