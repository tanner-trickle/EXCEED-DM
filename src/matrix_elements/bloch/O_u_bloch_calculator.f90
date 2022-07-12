module O_u_bloch_calculator

    use prec_util, only: dp

    use elec_state_bloch_type

    implicit none

contains

    subroutine compute_O_u_1(u, O_u_1, init_state)

        implicit none

        complex(dp) :: u(:, :, :, :)
        complex(dp) :: O_u_1(:, :, :, :)
        class(elec_state_bloch_t) :: init_state

        O_u_1 = u

    end subroutine

    subroutine compute_O_u_v(u, O_u_v, G_grid, init_state)

        use constants_util

        use FFT_util

        implicit none

        complex(dp) :: u(:, :, :, :)
        complex(dp) :: O_u_v(:, :, :, :, :)
        real(dp) :: G_grid(:, :, :, :)
        class(elec_state_bloch_t) :: init_state

        complex(dp), allocatable :: u_FT(:, :, :, :)

        integer :: d, s

        real(dp), allocatable :: v_vec_array(:, :, :, :)

        allocate(u_FT(size(u, 1), size(u, 2), size(u, 3), size(u, 4)), source = ( 0.0_dp, 0.0_dp ))
        allocate(v_vec_array(size(u, 1), size(u, 2), size(u, 3), 3), source = 0.0_dp)

        v_vec_array = G_grid/m_elec

        ! compute FFT
        do s = 1, size(u, 4)
            call dfftw_execute_dft(init_state%FFT_plans(1, :), u(:, :, :, s), u_FT(:, :, :, s)) 
        end do

        ! Normalize
        u_FT = (1.0_dp*size(u, 1)*size(u, 2)*size(u, 3))**(-1)*u_FT

        do d = 1, 3

            v_vec_array(:, :, :, d) = v_vec_array(:, :, :, d) + init_state%k_vec(d)/m_elec

            do s = 1, size(u, 4)

                call dfftw_execute_dft(init_state%FFT_plans(2, :),&
                    v_vec_array(:, :, :, d)*u_FT(:, :, :, s),& 
                    O_u_v(:, :, :, s, d)) 

            end do
        end do

    end subroutine

    subroutine compute_O_u_v2(u, O_u_v2, G_grid, init_state)

        use constants_util
        use FFT_util

        implicit none

        complex(dp) :: u(:, :, :, :)
        complex(dp) :: O_u_v2(:, :, :, :)
        real(dp) :: G_grid(:, :, :, :)
        class(elec_state_bloch_t) :: init_state

        complex(dp), allocatable :: u_FT(:, :, :, :)

        integer :: d, s

        integer :: g1, g2, g3

        real(dp), allocatable :: v2_vec_array(:, :, :)

        allocate(u_FT(size(u, 1), size(u, 2), size(u, 3), size(u, 4)), source = ( 0.0_dp, 0.0_dp ))
        allocate(v2_vec_array(size(u, 1), size(u, 2), size(u, 3)), source = 0.0_dp)

        do g3 = 1, init_state%n_x_grid(3)
            do g2 = 1, init_state%n_x_grid(2)
                do g1 = 1, init_state%n_x_grid(1)

                    v2_vec_array(g1, g2, g3) = (norm2( init_state%k_vec + G_grid(g1, g2, g3, :) )/m_elec)**2

                end do
            end do
        end do

        ! compute FFT
        do s = 1, size(u, 4)
            call dfftw_execute_dft(init_state%FFT_plans(1, :), u(:, :, :, s), u_FT(:, :, :, s)) 
        end do

        ! Normalize
        u_FT = (1.0_dp*size(u, 1)*size(u, 2)*size(u, 3))**(-1)*u_FT

        do s = 1, size(u, 4)

            call dfftw_execute_dft(init_state%FFT_plans(2, :),&
                v2_vec_array*u_FT(:, :, :, s),& 
                O_u_v2(:, :, :, s)) 

        end do

    end subroutine

    subroutine compute_O_u_s(u, O_u_s, init_state)

        use math_util

        implicit none

        complex(dp) :: u(:, :, :, :)
        complex(dp) :: O_u_s(:, :, :, :, :)
        class(elec_state_bloch_t) :: init_state

        integer :: d, s

        integer :: g1, g2, g3

        real(dp) :: sigma_d(2, 2)

        if ( size(u, 4) == 1 ) then

            O_u_s = ( 0.0_dp, 0.0_dp )

        else

            do d = 1, 3

                sigma_d = pauli_spin_matrix(d)        

                do g3 = 1, size(u, 3)
                    do g2 = 1, size(u, 2)
                        do g1 = 1, size(u, 1)

                            O_u_s(g1, g2, g3, :, d) = &
                                matmul(sigma_d, u(g1, g2, g3, :))

                        end do
                    end do
                end do
            end do

        end if

    end subroutine

    subroutine compute_O_u_vds(u, O_u_vds, G_grid, init_state)

        use math_util
        use FFT_util
        use constants_util

        implicit none

        complex(dp) :: u(:, :, :, :)
        complex(dp) :: O_u_vds(:, :, :, :)
        class(elec_state_bloch_t) :: init_state
        real(dp) :: G_grid(:, :, :, :)

        complex(dp), allocatable :: u_FT(:, :, :, :)
        complex(dp), allocatable :: O_u_FT(:, :, :, :)

        integer :: d, s

        integer :: g1, g2, g3

        complex(dp), allocatable :: vds_array(:, :, :, :, :)

        if ( size(u, 4) == 1 ) then

            O_u_vds = ( 0.0_dp, 0.0_dp )

        else

            allocate(u_FT(size(u, 1), size(u, 2), size(u, 3), size(u, 4)), source = ( 0.0_dp, 0.0_dp ))
            allocate(O_u_FT(size(u, 1), size(u, 2), size(u, 3), size(u, 4)), source = ( 0.0_dp, 0.0_dp ))
            allocate(vds_array(size(u, 1), size(u, 2), size(u, 3), 2, 2), source = ( 0.0_dp, 0.0_dp ))

            do g3 = 1, size(u, 3)
                do g2 = 1, size(u, 2)
                    do g1 = 1, size(u, 1)

                        do d = 1, 3
                            vds_array(g1, g2, g3, :, :) = vds_array(g1, g2, g3, :, :) + &
                                m_elec**(-1)*(init_state%k_vec(d) + G_grid(g1, g2, g3, d))*pauli_spin_matrix(d)
                        end do

                    end do
                end do
            end do

            ! compute FFT
            do s = 1, size(u, 4)
                call dfftw_execute_dft(init_state%FFT_plans(1, :), u(:, :, :, s), u_FT(:, :, :, s)) 
            end do

            ! Normalize
            u_FT = (1.0_dp*size(u, 1)*size(u, 2)*size(u, 3))**(-1)*u_FT

            do g3 = 1, size(u, 3)
                do g2 = 1, size(u, 2)
                    do g1 = 1, size(u, 1)

                        O_u_FT(g1, g2, g3, :) = matmul(vds_array(g1, g2, g3, :, :), u_FT(g1, g2, g3, :))

                    end do
                end do
            end do

            do s = 1, size(u, 4)

                call dfftw_execute_dft(init_state%FFT_plans(2, :),&
                    O_u_FT(:, :, :, s), &
                    O_u_vds(:, :, :, s)) 

            end do

        end if

    end subroutine

    subroutine compute_O_u_vxs(u, O_u_vxs, G_grid, init_state)

        use math_util
        use FFT_util
        use constants_util

        implicit none

        complex(dp) :: u(:, :, :, :)
        complex(dp) :: O_u_vxs(:, :, :, :, :)
        class(elec_state_bloch_t) :: init_state
        real(dp) :: G_grid(:, :, :, :)

        real(dp) :: G_vec(3)

        complex(dp), allocatable :: u_FT(:, :, :, :)
        complex(dp), allocatable :: O_u_FT(:, :, :, :, :)

        integer :: d, s
        integer :: g1, g2, g3

        complex(dp), allocatable :: vxs_array(:, :, :, :, :, :)

        if ( size(u, 4) == 1 ) then

            O_u_vxs = ( 0.0_dp, 0.0_dp )

        else

            allocate(u_FT(size(u, 1), size(u, 2), size(u, 3), size(u, 4)), source = ( 0.0_dp, 0.0_dp ))
            allocate(O_u_FT(size(u, 1), size(u, 2), size(u, 3), size(u, 4), 3), source = ( 0.0_dp, 0.0_dp ))
            allocate(vxs_array(size(u, 1), size(u, 2), size(u, 3), 3, 2, 2), source = ( 0.0_dp, 0.0_dp ))

            do g3 = 1, size(u, 3)
                do g2 = 1, size(u, 2)
                    do g1 = 1, size(u, 1)

                        G_vec = G_grid(g1, g2, g3, :)

                        vxs_array(g1, g2, g3, 1, :, :) = m_elec**(-1)*(&
                            (init_state%k_vec(2) + G_vec(2))*pauli_spin_matrix(3)&
                            - (init_state%k_vec(3) + G_vec(3))*pauli_spin_matrix(2))

                        vxs_array(g1, g2, g3, 2, :, :) = - m_elec**(-1)*(&
                             (init_state%k_vec(1) + G_vec(1))*pauli_spin_matrix(3)&
                            - (init_state%k_vec(3) + G_vec(3))*pauli_spin_matrix(1))

                        vxs_array(g1, g2, g3, 3, :, :) = m_elec**(-1)*(&
                             (init_state%k_vec(1) + G_vec(1))*pauli_spin_matrix(2)&
                            - (init_state%k_vec(2) + G_vec(2))*pauli_spin_matrix(1))

                    end do
                end do
            end do

            ! compute FFT
            do s = 1, size(u, 4)
                call dfftw_execute_dft(init_state%FFT_plans(1, :), u(:, :, :, s), u_FT(:, :, :, s)) 
            end do

            ! Normalize
            u_FT = (1.0_dp*size(u, 1)*size(u, 2)*size(u, 3))**(-1)*u_FT

            do g3 = 1, size(u, 3)
                do g2 = 1, size(u, 2)
                    do g1 = 1, size(u, 1)

                        do d = 1, 3

                            O_u_FT(g1, g2, g3, :, d) = matmul(vxs_array(g1, g2, g3, d, :, :), u_FT(g1, g2, g3, :))

                        end do

                    end do
                end do
            end do

            do s = 1, size(u, 4)
                do d = 1, 3

                    call dfftw_execute_dft(init_state%FFT_plans(2, :),&
                        O_u_FT(:, :, :, s, d), &
                        O_u_vxs(:, :, :, s, d)) 

                end do
            end do

        end if

    end subroutine

end module
