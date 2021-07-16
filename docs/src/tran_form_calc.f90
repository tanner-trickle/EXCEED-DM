module tran_form_calc
    !! Compute the transition form factors

    use prec
    use constants

    use DFT_parameters

    implicit none

    !! Compute different transition form factors for spin-independent
    !! and spin-dependent wave functions 
    interface calc_tran_form_1
        module procedure calc_tran_form_1_no_spin
        module procedure calc_tran_form_1_spin
    end interface

    interface calc_tran_form_v
        module procedure calc_tran_form_v_no_spin
        module procedure calc_tran_form_v_spin
    end interface

    interface calc_tran_form_v2
        module procedure calc_tran_form_v2_no_spin
        module procedure calc_tran_form_v2_spin
    end interface
    
contains

    !! TODO
    subroutine calc_tran_form_1_no_spin(tran_form, verbose)
        !! Computes the spin independent scalar transition form factor:
        !!
        !! $$\begin{align}
        !!      f_{i, i', \mathbf{k}} = \sum_\mathbf{G} \widetilde{u}_{i' \mathbf{k} \mathbf{G}}^* \widetilde{u}_{i \mathbf{k}
        !! \mathbf{G}} \nonumber
        !! \end{align}$$ 
        !!
        !! for a given i, i'
        implicit none




    end subroutine

    subroutine calc_tran_form_1_spin(tran_form, verbose)
        !! Computes the spin dependent scalar transition form factor:
        !!
        !! $$\begin{align}
        !!      f_{i, i', \mathbf{k}}^{ss'} = \sum_\mathbf{G} {\widetilde{u}_{i' \mathbf{k} \mathbf{G}}^{s'}}^* \widetilde{u}_{i \mathbf{k}
        !! \mathbf{G}}^s \nonumber
        !! \end{align}$$ 
        !!
        !! for a given i, i'

    end subroutine

    subroutine calc_tran_form_v_no_spin()

    end subroutine

    subroutine calc_tran_form_v_spin()

    end subroutine

    subroutine calc_tran_form_v2_no_spin()

    end subroutine

    subroutine calc_tran_form_v2_spin()

    end subroutine



    subroutine calc_tran_form(tran_form, tran_form_vec, omega_k, &
        val_id, cond_id, wfc_FT_i, wfc_FT_f, verbose)

        implicit none

        complex(dp) :: tran_form(n_k)
        complex(dp) :: tran_form_vec(3, n_k)

        real(dp) :: omega_k(n_k)

        integer :: val_id, cond_id

        complex(dp) :: wfc_FT_i(n_k, n_in_G)
        complex(dp) :: wfc_FT_f(n_k, n_in_G)

        logical, optional :: verbose

        integer :: k, g

        real(dp) :: p_vec(3)
        real(dp) :: p_sq

        tran_form = (0.0_dp, 0.0_dp)
        tran_form_vec = (0.0_dp, 0.0_dp)
        omega_k = 0.0_dp

        do k = 1, n_k

            omega_k(k) = energy_bands(k, cond_id) - energy_bands(k, val_id)

            do g = 1, n_in_G

                p_vec = k_grid_xyz(k, :) + in_G_grid_xyz(g, :)

                p_sq = dot_product(p_vec, p_vec)

                tran_form(k) = tran_form(k) +&
                    (p_sq/m_elec**2)*conjg(wfc_FT_f(k, g))*wfc_FT_i(k, g)

                tran_form_vec(:, k) = tran_form_vec(:, k) + &
                    (p_vec/m_elec)*conjg(wfc_FT_f(k, g))*wfc_FT_i(k, g)

            end do

        end do

    end subroutine

end module
