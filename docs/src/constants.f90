module constants
    !! Collection of useful constants
    use prec

    implicit none

    real(dp), parameter :: m_elec = 511.0e3_dp
        !! Electron mass
        !!
        !! Units : eV

    real(dp), parameter :: alpha_EM = 1.0_dp/137.0_dp
        !! Fine structure constant

    real(dp), parameter :: a0 = 2.681336e-4_dp
       !! Bohr radius
       !!
       !! Units : eV^(-1)

    real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)

    complex(dp), parameter :: ii = (0.0_dp, 1.0_dp)
        !! Imaginary unit

    real(dp), parameter :: e_EM = sqrt(4*pi*alpha_EM)
        !! Unit of electric charge

end module
