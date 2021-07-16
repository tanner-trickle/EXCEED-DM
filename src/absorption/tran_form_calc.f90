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
    subroutine calc_tran_form_1_no_spin(tran_form, val_id, cond_id, wfc_FT_i, wfc_FT_f, verbose)
        !! Computes the spin independent scalar transition form factor:
        !!
        !! $$\begin{align}
        !!      f_{i, i', \mathbf{k}} = \sum_\mathbf{G} \widetilde{u}_{i' \mathbf{k} \mathbf{G}}^* \widetilde{u}_{i \mathbf{k}
        !! \mathbf{G}} = \delta_{i i'} \nonumber
        !! \end{align}$$ 
        !!
        !! for a given i, i'
        !!
        !! Have this function here to have a nice interface.
        implicit none

        logical, optional :: verbose

        complex(dp) :: tran_form(:)

        integer :: val_id, cond_id

        complex(dp) :: wfc_FT_i(:, :)
        complex(dp) :: wfc_FT_f(:, :)

        integer :: k, g

        if ( val_id /= cond_id ) then

            tran_form = (0.0_dp, 0.0_dp)

        end if

    end subroutine

    subroutine calc_tran_form_1_spin(tran_form, val_id, cond_id, wfc_FT_i, wfc_FT_f, verbose)
        !! Computes the spin dependent scalar transition form factor:
        !!
        !! $$\begin{align}
        !!      f_{i, i', \mathbf{k}}^{ss'} = \sum_\mathbf{G} {\widetilde{u}_{i' \mathbf{k} \mathbf{G}}^{s'}}^* \widetilde{u}_{i \mathbf{k}
        !! \mathbf{G}}^s \nonumber
        !! \end{align}$$ 
        !!
        !! for a given i, i'.
        !!
        !! Note:
        !! $$\begin{align}
        !!      \sum_{s} f_{i, i', \mathbf{k}}^{s s} = \delta_{i i'}
        !! \end{align}$$ 

        implicit none

        logical, optional :: verbose

        complex(dp) :: tran_form(:, :, :)

        integer :: val_id, cond_id

        complex(dp) :: wfc_FT_i(:, :, :)
        complex(dp) :: wfc_FT_f(:, :, :)

        integer :: k, s, sp

        tran_form = (0.0_dp, 0.0_dp)

        if ( val_id /= cond_id ) then

            do k = 1, n_k
                do s = 1, 2
                    do sp = 1, 2

                        tran_form(k, s, sp) = tran_form(k, s, sp) + sum(conjg(wfc_FT_f(k, :, sp))*wfc_FT_i(k, :, s))

                    end do
                end do
            end do

        else

            tran_form = (0.0_dp, 0.0_dp)

        end if

    end subroutine

    subroutine calc_tran_form_v_no_spin(tran_form, val_id, cond_id, wfc_FT_i, wfc_FT_f, verbose)

        implicit none

        logical, optional :: verbose

        complex(dp) :: tran_form(:, :)

        integer :: val_id, cond_id

        complex(dp) :: wfc_FT_i(:, :)
        complex(dp) :: wfc_FT_f(:, :)

        integer :: k, g

        real(dp) :: p_vec(3)

        do k = 1, n_k

            do g = 1, n_in_G

                p_vec = k_grid_xyz(k, :) + in_G_grid_xyz(g, :)

                tran_form(:, k) = tran_form(:, k) + &
                    (p_vec/m_elec)*conjg(wfc_FT_f(k, g))*wfc_FT_i(k, g)

            end do

        end do

    end subroutine

    subroutine calc_tran_form_v_spin(tran_form, val_id, cond_id, wfc_FT_i, wfc_FT_f, verbose)

        implicit none

        logical, optional :: verbose

        complex(dp) :: tran_form(:, :, :, :)

        integer :: val_id, cond_id

        complex(dp) :: wfc_FT_i(:, :, :)
        complex(dp) :: wfc_FT_f(:, :, :)

        integer :: k, g, s, sp

        real(dp) :: p_vec(3)

        tran_form = (0.0_dp, 0.0_dp)

        do k = 1, n_k

            do s = 1, 2
                do sp = 1, 2

                    do g = 1, n_in_G

                        p_vec = k_grid_xyz(k, :) + in_G_grid_xyz(g, :)

                        tran_form(:, k, s, sp) = tran_form(:, k, s, sp) + &
                            (p_vec/m_elec)*conjg(wfc_FT_f(k, g, sp))*wfc_FT_i(k, g, s)

                    end do

                end do
            end do

        end do

    end subroutine

    subroutine calc_tran_form_v2_no_spin(tran_form, val_id, cond_id, wfc_FT_i, wfc_FT_f, verbose)

        implicit none

        logical, optional :: verbose

        complex(dp) :: tran_form(:)

        integer :: val_id, cond_id

        complex(dp) :: wfc_FT_i(:, :)
        complex(dp) :: wfc_FT_f(:, :)

        integer :: k, g

        real(dp) :: p_vec(3)
        real(dp) :: p_sq

        tran_form = (0.0_dp, 0.0_dp)

        do k = 1, n_k

            do g = 1, n_in_G

                p_vec = k_grid_xyz(k, :) + in_G_grid_xyz(g, :)

                p_sq = norm2(p_vec)**2

                tran_form(k) = tran_form(k) + &
                    (p_sq/m_elec**2)*conjg(wfc_FT_f(k, g))*wfc_FT_i(k, g)

            end do

        end do

    end subroutine

    subroutine calc_tran_form_v2_spin(tran_form, val_id, cond_id, wfc_FT_i, wfc_FT_f, verbose)

        implicit none

        logical, optional :: verbose

        complex(dp) :: tran_form(:, :, :)

        integer :: val_id, cond_id

        complex(dp) :: wfc_FT_i(:, :, :)
        complex(dp) :: wfc_FT_f(:, :, :)

        integer :: k, g, s, sp

        real(dp) :: p_vec(3)
        real(dp) :: p_sq

        tran_form = (0.0_dp, 0.0_dp)

        do k = 1, n_k

            do s = 1, 2
                do sp = 1, 2

                    do g = 1, n_in_G

                        p_vec = k_grid_xyz(k, :) + in_G_grid_xyz(g, :)

                        p_sq = norm2(p_vec)**2

                        tran_form(k, s, sp) = tran_form(k, s, sp) + &
                            (p_sq/m_elec**2)*conjg(wfc_FT_f(k, g, sp))*wfc_FT_i(k, g, s)

                    end do

                end do
            end do

        end do

    end subroutine

end module
