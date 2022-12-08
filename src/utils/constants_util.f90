module constants_util
    ! Defines useful constants and conversion factors.
    !
    ! .. note:: No external dependancies.

    use iso_fortran_env

    implicit none

    ! Conversion factors
    real(real64), parameter :: inv_Ang_to_eV      = 1973.37_real64
    real(real64), parameter :: inv_eV_to_cm       = 1.97327e-5_real64
    real(real64), parameter :: Ang_to_inv_eV      = 5.068e-4_real64
    real(real64), parameter :: g_to_eV            = 5.61e32_real64
    real(real64), parameter :: kg_to_eV           = 5.61e35_real64
    real(real64), parameter :: inv_cm_to_eV       = 1.97327e-5_real64
    real(real64), parameter :: eV_to_inv_cm       = 5.06773e4_real64
    real(real64), parameter :: km_per_sec_to_none = 3.33563e-6_real64
    real(real64), parameter :: yr_to_inv_eV       = 4.7912e22_real64
    real(real64), parameter :: inv_AMU_to_inv_eV  = 1.074e-9_real64

    ! Constants
    real(real64), parameter :: m_elec = 511.0e3_real64 ! Electron mass, :math:`m_e`, eV

    real(real64), parameter :: alpha_EM = 1.0_real64/137.0_real64 ! Fine structure constant, :math:`\alpha`

    real(real64), parameter :: a0 = 2.681336e-4_real64 ! Bohr radius, :math:`a_0`, :math:`\text{eV}^{-1}`

    real(real64), parameter :: pi = 4.0_real64*atan(1.0_real64) ! :math:`\pi`

    complex(real64), parameter :: ii = (0.0_real64, 1.0_real64) ! Imaginary unit

    real(real64), parameter :: e_EM = sqrt(4*pi*alpha_EM) ! Unit of electric charge, :math:`e`

end module
