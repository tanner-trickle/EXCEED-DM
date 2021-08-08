module units
    !! Unit conversion factors.
    !!
    !! To convert a variable from \( x \) units to \( y \) units: var_y = x_to_y * var_x.

    use prec

    implicit none

    real(dp) :: inv_Ang_to_eV = 1973.37_dp
    real(dp) :: inv_eV_to_cm = 1.97327e-5_dp
    real(dp) :: Ang_to_inv_eV = 5.068e-4_dp 

    real(dp) :: g_to_eV = 5.61e32_dp
    real(dp) :: kg_to_eV = 5.61e35_dp
    real(dp) :: inv_cm_to_eV = 1.97327e-5_dp
    real(dp) :: eV_to_inv_cm = 5.06773e4_dp 

    real(dp) :: km_per_sec_to_none = 3.33563e-6_dp

    real(dp) :: yr_to_inv_eV = 4.7912e22_dp 

end module
