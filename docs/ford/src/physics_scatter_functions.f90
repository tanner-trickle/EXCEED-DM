module physics_scatter_functions
    !! Physics functions needed for DM-electron scattering rate calcuations.

    use prec
    use constants
    use units

    use dm_model_type
    
    implicit none

contains

   function v_minus(q_vec, mX, vE_vec, omega) result(v_m)
       !! \( v_- = \frac{1}{q} 
       !! | \mathbf{q} \cdot \mathbf{v}_E + \frac{q^2}{2 m_\chi} + \omega | \)
       !!
       !! Note : assumes \( v_- < v_\text{esc} \).
       !!
       !! Units : None

       implicit none

       real(dp) :: mX, omega
       real(dp) :: q_mag

       real(dp) :: v_m

       real(dp) :: q_vec(3)
       real(dp) :: vE_vec(3)

       q_mag = norm2(q_vec)
       v_m = (1/q_mag)*abs(dot_product(q_vec, vE_vec) + 0.5_dp*q_mag**2/mX + omega)

   end function

   function g_func(q, v_m, dm_model) result (g_fun)
       !! Kinematic, \( g(v_-) \) function.
       !!
       !! \( g(\mathbf{q}, \omega) = 2 \pi \int d^3\mathbf{v} 
       !! f_\chi(\mathbf{v}) \delta(\omega - \omega_\mathbf{q}) \)
       !!
       !! \( g(v_-) = \frac{2 \pi^2 v_0^2}{N_0 q} \left( \exp{(-v_-^2/v_0^2)} - \exp{(-v_\text{esc}^2/v_0^2)} \right) \)
       !!
       !! Units : \( \text{eV}^{-1} \)
       implicit none

       real(dp) :: q, v_m
       type(dm_model_t) :: dm_model
       real(dp) :: g_fun

       g_fun = (dm_model%g_func_c1/q)*(exp(-(v_m/dm_model%v0)**2) - dm_model%g_func_c2)

   end function

   function red_mass(m1, m2) result(mu)
       !! \( \mu = \frac{m_1 m_2}{m_1 + m_2} \)
       implicit none

       real(dp) :: m1, m2
       real(dp) :: mu 

       mu = m1*m2/(m1 + m2)

   end function

   function F_med_sq_func(q_mag, power) result (F_med_sq_val)
       !! Mediator form factor squared, 
       !! \( \mathcal{F}_\text{med}^2 = \left( \frac{\alpha m_e}{q} \right)^{2 \beta} \).
       !!
       !! Units : None
       implicit none

       real(dp) :: q_mag
       real(dp) :: F_med_sq_val
       real(dp) :: power

       f_med_sq_val = (alpha_EM*m_elec/q_mag)**(2*power)

   end function

end module
