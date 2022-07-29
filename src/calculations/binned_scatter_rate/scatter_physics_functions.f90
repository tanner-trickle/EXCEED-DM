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

    function compute_v_minus(q_vec_list, q_mag_list, mX, vE_vec, omega) result(v_m_list)

        implicit none

        real(dp), intent(in) :: q_vec_list(:, :)
        real(dp), intent(in) :: q_mag_list(:)
        real(dp), intent(in) :: mX
        real(dp), intent(in) :: vE_vec(:)
        real(dp), intent(in) :: omega

        real(dp) :: v_m_list(size(q_mag_list))

        v_m_list = q_mag_list**(-1)*&
            abs(&
                matmul(q_vec_list, vE_vec) + 0.5_dp*q_mag_list**2/mX + omega &
                )

    end function

end module
