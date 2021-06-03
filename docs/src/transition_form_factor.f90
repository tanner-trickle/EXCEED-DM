module transition_form_factor
    !! --------------------------------
    !!     Transition Form Factor (TFF)
    !! --------------------------------
    !!
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
    !!         1 - \( \mathcal{O} = 1 \) [ vc ]
    !!     </li>
    !!     <li> 
    !!         2 - \( \mathcal{O} = \mathbf{S}_e \) [ vc, <b>s</b> ]
    !!     </li>
    !! </ul>
    !!
    !! Notes : Bracketed quantities indicate what transition types are currently supported, <b>s</b> indicates its only
    !! supported for spin dependent wave functions. 

    use prec
    use constants

    implicit none

    integer :: tff_id(2) = [1, 1]

    NAMELIST /tff/ tff_id

    interface calc_tff_vc
        module procedure calc_tff_vc_no_spin
        module procedure calc_tff_vc_spin
    end interface

contains

    subroutine print_tff_input(verbose)
        implicit none

        logical, optional :: verbose

        if ( verbose ) then
            print*, '----------------------------------------'
            print*, '    ----------------------'
            print*, '    Transition Form Factor'
            print*, '    ----------------------'
            print*
            print*, '        ID  : ', tff_id
            print*
        end if

    end subroutine

    subroutine load_tff_input(filename, verbose)

        implicit none

        character(len=*) :: filename

        logical, optional :: verbose

        logical :: file_exists

        integer :: error

        if ( verbose ) then

            print*, 'Loading Transition Form Factor parameters...'
            print*

        end if

        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            open(100, file = trim(filename), iostat = error)
            read(100, nml=tff, iostat=error)
            close(100)

            call print_tff_input(verbose = verbose)

            if ( verbose ) then

                print*, '----------------------------------------'
                print*

            end if

        else

            if ( verbose ) then

                print*, '!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*
                print*, '   Input file for tff parameters : ', trim(filename), ' does NOT exist.'
                print*
                print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*

            end if

            stop

        end if

    end subroutine

    subroutine calc_tff_vc_spin(TFF, wfc_i, wfc_f, n_FFT_grid, FFT_plan, verbose)
        !! Compute the tff for the valence to conduction calculation
        !!
        !! TFF = FFT( \( u_{f, k_f}^* \mathcal{O}_1 u_{i, k_i} \)  ) \( \times \) FFT( \( u_{f, k_f}^* \mathcal{O}_2 u_{i, k_i} \) )
        !! \( ^* \)
        !!
        !! Dim : [ n_FFT_grid ]
        !!
        !! Units : None

        use DFT_parameters
        use math_mod

        implicit none

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

                    op = pauli_spin_matrix(0)

                    do g3 = 1, n_FFT_grid(3)
                        do g2 = 1, n_FFT_grid(2)
                            do g1 = 1, n_FFT_grid(1)

                                Tx_s_terms(1, g1, g2, g3) = dot_product(& 
                                        conjg(wfc_f(:, g1, g2, g3)),&
                                        matmul( op , wfc_i(:, g1, g2, g3) )&
                                        )

                            end do
                        end do
                    end do

                    call dfftw_execute_dft(FFT_plan, Tx_s_terms(1, :, :, :) , TFF_s_terms(1, :, :, :)) 

                    TFF_s_terms(2, :, :, :) = conjg(TFF_s_terms(1, :, :, :))

                    TFF = TFF_s_terms(1, :, :, :)*TFF_s_terms(2, :, :, :)

                case ( 2 )

                    ! pauli matrix operator

                    do i = 1, 3

                        do g3 = 1, n_FFT_grid(3)
                            do g2 = 1, n_FFT_grid(2)
                                do g1 = 1, n_FFT_grid(1)

                                    op = pauli_spin_matrix(i)

                                    Tx_v_terms(i, 1, g1, g2, g3) = dot_product(& 
                                            conjg(wfc_f(:, g1, g2, g3)),&
                                            matmul( op , wfc_i(:, g1, g2, g3) )&
                                            )

                                end do
                            end do
                        end do

                        call dfftw_execute_dft(FFT_plan, Tx_v_terms(i, 1, :, :, :) , TFF_v_terms(i, 1, :, :, :)) 

                        TFF_v_terms(i, 2, :, :, :) = conjg(TFF_v_terms(i, 1, :, :, :))

                        TFF = TFF + TFF_v_terms(i, 1, :, :, :)*TFF_v_terms(i, 2, :, :, :)

                    end do

                case default

                    op = pauli_spin_matrix(0)

                    do g3 = 1, n_FFT_grid(3)
                        do g2 = 1, n_FFT_grid(2)
                            do g1 = 1, n_FFT_grid(1)

                                Tx_s_terms(1, g1, g2, g3) = dot_product(& 
                                        conjg(wfc_f(:, g1, g2, g3)),&
                                        matmul( op , wfc_i(:, g1, g2, g3) )&
                                        )

                            end do
                        end do
                    end do

                    call dfftw_execute_dft(FFT_plan, Tx_s_terms(1, :, :, :) , TFF_s_terms(1, :, :, :)) 

                    TFF_s_terms(2, :, :, :) = conjg(TFF_s_terms(1, :, :, :))

                    TFF = TFF_s_terms(1, :, :, :)*TFF_s_terms(2, :, :, :)

            end select

        end if

    end subroutine

    subroutine calc_tff_vc_no_spin(TFF, wfc_i, wfc_f, n_FFT_grid, FFT_plan, verbose)
        !! Compute the tff for the valence to conduction calculation
        !!
        !! TFF = FFT( \( u_{f, k_f}^* \mathcal{O}_1 u_{i, k_i} \)  ) \( \times \) FFT( \( u_{f, k_f}^* \mathcal{O}_2 u_{i, k_i} \) )
        !! \( ^* \)
        !!
        !! Dim : [ n_FFT_grid ]
        !!
        !! Units : None

        use DFT_parameters

        implicit none

        logical, optional :: verbose

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

                    TFF_s_terms(2, :, :, :)  = conjg(TFF_s_terms(1, :, :, :))

                    TFF = TFF_s_terms(1, :, :, :)*TFF_s_terms(2, :, :, :)

                case default

                    Tx_s_terms(1, :, :, :) = conjg(wfc_f)*wfc_i

                    call dfftw_execute_dft(FFT_plan, Tx_s_terms(1, :, :, :) , TFF_s_terms(1, :, :, :)) 

                    TFF_s_terms(2, :, :, :)  = conjg(TFF_s_terms(1, :, :, :))

                    TFF = TFF_s_terms(1, :, :, :)*TFF_s_terms(2, :, :, :)

            end select

        end if

    end subroutine

end module
