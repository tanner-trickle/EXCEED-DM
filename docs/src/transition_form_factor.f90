module transition_form_factor
    ! See the documentation for LaTeX'ed math eqs.
    !!
    !! Define general scattering operators such that the scattering rate for
    !! any DM model can be computed. This will generalize the transition form factor (TFF) defined as
    !!
    !! $$\begin{align}
    !!    T_{if}(\mathcal{O}_1, \mathcal{O}_2) \equiv \langle f | e^{i \mathbf{q} \cdot \mathbf{x}} \mathcal{O}_1 | i \rangle
    !! \cdot \langle i | e^{-i \mathbf{q} \cdot \mathbf{x}} \mathcal{O}_2^* | f \rangle 
    !! \end{align}$$
    !! 
    !! All spin independent results come from \( \mathcal{O}_1 = \mathcal{O}_2 = 1 \). Note that this quantity is different 
    !! than the mediator form factor and screening factor, both of which can be set independently of this function.
    !!
    !! Each \( \mathcal{O} \) that the user can specify will have a unique index.
    !!
    !! A catalog will be kept of all the operators and their corresponding index.
    !!
    !! The calculation can become more difficult technically depending on the operator so a hierarchy will be kept in
    !! order to simplify the calculation when possible. 
    !!
    !! TODO: For now we will assume that only a single T_if needs to be computed. Future work could improve on this by
    !! setting up wrapper functions which compute multiple T's and sum them. This will be necessary when the scattering operator
    !! has more than one term, e.g. \( \mathcal{O} = \frac{k}{m_e} + \mathbf{S}_e \).
    !!
    !! Catalog :
    !! <ul>
    !!     <li> 
    !!         1 - \( \mathcal{O} = 1 \) [ vc (SI/SD), cc (SI) ]
    !!     </li>
    !!     <li>
    !!         2 - \( \mathcal{O} = \mathbf{S}_e \) [vc (SD)]
    !!     </li>
    !! </ul>
    !!
    !! Notes : Bracketed quantities indicate what transition types are currently supported, <b>s</b> indicates its only
    !! supported for spin dependent wave functions. 

    use prec
    use constants
    use math_mod

    implicit none

    interface calc_tff_pw_pw
        module procedure calc_tff_pw_pw_no_spin
        module procedure calc_tff_pw_pw_spin
    end interface

contains

    subroutine calc_tff_pw_pw_spin(tff_id, TFF, wfc_i, wfc_f, n_FFT_grid, FFT_plan, verbose)
        !! Compute the transition form factor between two spin dependent wave functions defined on a uniform grid in the unit cell.
        !!
        !! TFF = FFT( \( u_{f, k_f}^* \mathcal{O}_1 u_{i, k_i} \)  ) \( \times \) FFT( \( u_{f, k_f}^* \mathcal{O}_2 u_{i, k_i} \) )
        !! \( ^* \)
        !!
        !! Dim : [ n_FFT_grid ]
        !!
        !! Units : None

        implicit none

        integer :: tff_id(2)

        logical, optional :: verbose

        integer :: n_FFT_grid(3)

        integer :: FFT_plan(8)

        complex(dp) :: wfc_i(2, n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3))
        complex(dp) :: wfc_f(2, n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3))

        complex(dp) :: Tx_s_terms(2, n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3))
            !! scalar overlaps in position space
        complex(dp) :: Tx_v_terms(3, 2, n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3))
            !! vector overlaps in position space 

        complex(dp) :: TFF_s_terms(2, n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3))
            !! scalar
        complex(dp) :: TFF_v_terms(3, 2, n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3))
            !! vector

        real(dp) :: TFF(n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3))

        complex(dp) :: op(2, 2)

        integer :: i, g1, g2, g3

        Tx_s_terms = (0.0_dp, 0.0_dp)
        Tx_v_terms = (0.0_dp, 0.0_dp)

        TFF_s_terms = (0.0_dp, 0.0_dp)
        TFF_v_terms = (0.0_dp, 0.0_dp)

        if ( tff_id(1) == tff_id(2) ) then

            ! only have to compute one element since the other is the complex conjugate
            
            select case ( tff_id(1) )

                case ( 1 )

                    ! 1

                    do g3 = 1, n_FFT_grid(3)
                        do g2 = 1, n_FFT_grid(2)
                            do g1 = 1, n_FFT_grid(1)

                                Tx_s_terms(1, g1, g2, g3) = &
                                        conjg(wfc_f(1, g1, g2, g3))*wfc_i(1, g1, g2, g3) + &
                                        conjg(wfc_f(2, g1, g2, g3))*wfc_i(2, g1, g2, g3)

                            end do
                        end do
                    end do

                    call dfftw_execute_dft(FFT_plan, Tx_s_terms(1, :, :, :) , TFF_s_terms(1, :, :, :)) 

                    ! TFF_s_terms(2, :, :, :) = conjg(TFF_s_terms(1, :, :, :))

                    TFF = abs(TFF_s_terms(1, :, :, :))**2

                case ( 2 )

                    ! pauli matrix operator, \( \sigma^i \)

                    do g3 = 1, n_FFT_grid(3)
                        do g2 = 1, n_FFT_grid(2)
                            do g1 = 1, n_FFT_grid(1)

                                Tx_v_terms(1, 1, g1, g2, g3) = &
                                    conjg(wfc_f(1, g1, g2, g3))*wfc_i(2, g1, g2, g3) + &
                                    conjg(wfc_f(2, g1, g2, g3))*wfc_i(1, g1, g2, g3)

                                Tx_v_terms(2, 1, g1, g2, g3) = &
                                    -ii*conjg(wfc_f(1, g1, g2, g3))*wfc_i(2, g1, g2, g3) + &
                                    ii*conjg(wfc_f(2, g1, g2, g3))*wfc_i(1, g1, g2, g3)

                                Tx_v_terms(3, 1, g1, g2, g3) = &
                                    conjg(wfc_f(1, g1, g2, g3))*wfc_i(1, g1, g2, g3) - &
                                    conjg(wfc_f(2, g1, g2, g3))*wfc_i(2, g1, g2, g3)

                            end do
                        end do
                    end do

                    do i = 1, 3

                        call dfftw_execute_dft(FFT_plan, Tx_v_terms(i, 1, :, :, :) , TFF_v_terms(i, 1, :, :, :)) 

                        ! TFF_v_terms(i, 2, :, :, :) = conjg(TFF_v_terms(i, 1, :, :, :))

                        TFF = TFF + abs(TFF_v_terms(i, 1, :, :, :))**2

                    end do

                case default

                    do g3 = 1, n_FFT_grid(3)
                        do g2 = 1, n_FFT_grid(2)
                            do g1 = 1, n_FFT_grid(1)

                                Tx_s_terms(1, g1, g2, g3) = &
                                        conjg(wfc_f(1, g1, g2, g3))*wfc_i(1, g1, g2, g3) + &
                                        conjg(wfc_f(2, g1, g2, g3))*wfc_i(2, g1, g2, g3)

                            end do
                        end do
                    end do

                    call dfftw_execute_dft(FFT_plan, Tx_s_terms(1, :, :, :) , TFF_s_terms(1, :, :, :)) 

                    ! TFF_s_terms(2, :, :, :) = conjg(TFF_s_terms(1, :, :, :))

                    TFF = abs(TFF_s_terms(1, :, :, :))**2

            end select

        end if

    end subroutine

    subroutine calc_tff_pw_pw_no_spin(tff_id, TFF, wfc_i, wfc_f, n_FFT_grid, FFT_plan, verbose)
        !! Compute the transition form factor between two spin independent wave functions defined on a uniform grid in the unit cell.
        !!
        !! TFF = FFT( \( u_{f, k_f}^* \mathcal{O}_1 u_{i, k_i} \)  ) \( \times \) FFT( \( u_{f, k_f}^* \mathcal{O}_2 u_{i, k_i} \) )
        !! \( ^* \)
        !!
        !! Dim : [ n_FFT_grid ]
        !!
        !! Units : None

        implicit none

        logical, optional :: verbose

        integer :: tff_id(2)

        integer :: n_FFT_grid(3)

        integer :: FFT_plan(8)

        complex(dp) :: wfc_i(n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3))
        complex(dp) :: wfc_f(n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3))

        complex(dp) :: Tx_s_terms(2, n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3))
            !! scalar overlap in position space
        complex(dp) :: Tx_v_terms(3, 2, n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3))
            !! vector overlap in position space

        complex(dp) :: TFF_s_terms(2, n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3))
            !! scalar
        complex(dp) :: TFF_v_terms(3, 2, n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3))
            !! vector

        real(dp) :: TFF(n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3))

        Tx_s_terms = (0.0_dp, 0.0_dp)
        Tx_v_terms = (0.0_dp, 0.0_dp)

        TFF_s_terms = (0.0_dp, 0.0_dp)
        TFF_v_terms = (0.0_dp, 0.0_dp)

        if ( tff_id(1) == tff_id(2) ) then

            ! only have to compute one element since the other is the complex conjugate
            
            select case ( tff_id(1) )

                case ( 1 )

                    ! 1

                    Tx_s_terms(1, :, :, :) = conjg(wfc_f)*wfc_i

                    call dfftw_execute_dft(FFT_plan, Tx_s_terms(1, :, :, :) , TFF_s_terms(1, :, :, :)) 

                    ! TFF_s_terms(2, :, :, :)  = conjg(TFF_s_terms(1, :, :, :))

                    TFF = abs(TFF_s_terms(1, :, :, :))**2

                case default

                    Tx_s_terms(1, :, :, :) = conjg(wfc_f)*wfc_i

                    call dfftw_execute_dft(FFT_plan, Tx_s_terms(1, :, :, :) , TFF_s_terms(1, :, :, :)) 

                    ! TFF_s_terms(2, :, :, :)  = conjg(TFF_s_terms(1, :, :, :))

                    TFF = abs(TFF_s_terms(1, :, :, :))**2

            end select

        end if

    end subroutine

end module
