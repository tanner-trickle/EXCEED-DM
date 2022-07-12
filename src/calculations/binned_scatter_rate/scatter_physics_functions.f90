module scatter_physics_functions

    use prec_util, only: dp

    implicit none

contains

    subroutine fermi_factor_extra_FF(extra_FF, Zeff_list, Ef_list)

        use constants_util

        implicit none

        real(dp) :: extra_FF(:, :)
        real(dp) :: Zeff_list(:)
        real(dp) :: Ef_list(:)

        real(dp) :: fermi_val
        real(dp) :: fermi_factor

        integer :: i, f

        do f = 1, size(Ef_list)
            do i = 1, size(Zeff_list)

                fermi_val = 2.0_dp*pi*Zeff_list(i)*&
                    (alpha_EM*m_elec/sqrt(2.0_dp*m_elec*Ef_list(f)))

                fermi_factor = fermi_val*( 1.0_dp - exp(-fermi_val) )**(-1)

                extra_FF(i, f) = fermi_factor

            end do
        end do

    end subroutine

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
