module TIF_calculator_bloch_type

    use prec_util, only: dp

    use exdm_inputs_type

    use elec_state_bloch_type

    implicit none

    type :: TIF_calculator_bloch_t

        complex(dp), allocatable :: u_i(:, :, :, :)
        complex(dp), allocatable :: u_f(:, :, :, :)

        integer :: n_FFT_grid(3)
        integer :: backward_FFT_plan(8)

        integer :: n_q

        real(dp), allocatable :: G_vec_list(:, :)

        real(dp), allocatable :: init_G_grid(:, :, :, :)

        complex(dp), allocatable :: O_u_1(:, :, :, :)
        complex(dp), allocatable :: O_u_v(:, :, :, :, :)
        complex(dp), allocatable :: O_u_v2(:, :, :, :)
        complex(dp), allocatable :: O_u_s(:, :, :, :, :)
        complex(dp), allocatable :: O_u_vds(:, :, :, :)
        complex(dp), allocatable :: O_u_vxs(:, :, :, :, :)
        complex(dp), allocatable :: O_u_vivj(:, :, :, :, :, :)

        ! padded to size of FFT variables
        integer :: n_pad_grid(3)
        complex(dp), allocatable :: u_f_pad(:, :, :, :)
        complex(dp), allocatable :: O_u_1_pad(:, :, :, :)
        complex(dp), allocatable :: O_u_v_pad(:, :, :, :, :)
        complex(dp), allocatable :: O_u_v2_pad(:, :, :, :)
        complex(dp), allocatable :: O_u_s_pad(:, :, :, :, :)
        complex(dp), allocatable :: O_u_vds_pad(:, :, :, :)
        complex(dp), allocatable :: O_u_vxs_pad(:, :, :, :, :)
        complex(dp), allocatable :: O_u_vivj_pad(:, :, :, :, :, :)

        contains

            procedure :: setup => TIF_calculator_bloch_type_setup
            procedure :: init_initialize => TIF_calculator_bloch_type_init_initialize
            procedure :: fin_initialize => TIF_calculator_bloch_type_fin_initialize
            procedure :: compute_TIF_1 => TIF_calculator_bloch_type_compute_TIF_1
            procedure :: compute_TIF_v => TIF_calculator_bloch_type_compute_TIF_v
            procedure :: compute_TIF_v2 => TIF_calculator_bloch_type_compute_TIF_v2
            procedure :: compute_TIF_s => TIF_calculator_bloch_type_compute_TIF_s
            procedure :: compute_TIF_vds => TIF_calculator_bloch_type_compute_TIF_vds
            procedure :: compute_TIF_vxs => TIF_calculator_bloch_type_compute_TIF_vxs
            procedure :: compute_TIF_vivj => TIF_calculator_bloch_type_compute_TIF_vivj

    end type

contains

    subroutine TIF_calculator_bloch_type_setup(self, exdm_inputs, TIF_mask, init_state, fin_state, q0_limit)

        use FFT_util

        implicit none

        class(TIF_calculator_bloch_t) :: self

        type(exdm_inputs_t) :: exdm_inputs
        logical :: TIF_mask(:)

        class(elec_state_bloch_t) :: init_state
        class(elec_state_bloch_t) :: fin_state

        logical, optional :: q0_limit

        integer :: n_q, n_FFT

        integer :: i

        ! initialize u_i, u_f
        allocate(self%u_i(init_state%n_x_grid(1), &
                     init_state%n_x_grid(2), &
                     init_state%n_x_grid(3), &
                     init_state%spin_dof), source = ( 0.0_dp, 0.0_dp ) )

        allocate(self%u_f(fin_state%n_x_grid(1), &
                     fin_state%n_x_grid(2), &
                     fin_state%n_x_grid(3), &
                     fin_state%spin_dof), source = ( 0.0_dp, 0.0_dp ) )

        allocate(self%init_G_grid(init_state%n_x_grid(1), &
                     init_state%n_x_grid(2), &
                     init_state%n_x_grid(3), &
                     3), source = 0.0_dp )

        ! compute init G grid
        call create_FFT_G_vec_array(self%init_G_grid, init_state%n_x_grid, &
            exdm_inputs%material%k_red_to_xyz)

        ! only q -> 0 in q0_limit
        if ( .not. q0_limit ) then

            ! set the FFT grid size, doubled to avoid wrapping problems
            self%n_FFT_grid = init_state%n_x_grid + fin_state%n_x_grid - 1
            self%n_pad_grid = self%n_FFT_grid

            self%n_q = product(self%n_FFT_grid)

            ! create G_vec_list corresponding to the FFT grid
            allocate(self%G_vec_list(self%n_q, 3), source = 0.0_dp)
            call create_FFT_G_vec_list(self%G_vec_list, &
                                   self%n_FFT_grid, &
                                   exdm_inputs%material%k_red_to_xyz)


        else 

            do i = 1, 3
                self%n_pad_grid(i) = max( size(self%u_i, i), size(self%u_f, i) )
            end do

            self%n_q = 1

        end if

        call set_fft_plan_backward_3d(self%n_pad_grid, self%backward_FFT_plan)

        ! allocate arrays which are specified by the mask
        allocate(self%u_f_pad(self%n_pad_grid(1), self%n_pad_grid(2), self%n_pad_grid(3), fin_state%spin_dof), &
                    source = ( 0.0_dp, 0.0_dp ))

        if ( TIF_mask(1) ) then
            allocate(self%O_u_1(init_state%n_x_grid(1), &
                                init_state%n_x_grid(2), &
                                init_state%n_x_grid(3), &
                                init_state%spin_dof), source = ( 0.0_dp, 0.0_dp ))

            allocate(self%O_u_1_pad(self%n_pad_grid(1), &
                                self%n_pad_grid(2), &
                                self%n_pad_grid(3), &
                                init_state%spin_dof), source = ( 0.0_dp, 0.0_dp ))

        end if
        if ( TIF_mask(2) ) then
            allocate(self%O_u_v(init_state%n_x_grid(1), &
                                init_state%n_x_grid(2), &
                                init_state%n_x_grid(3), &
                                init_state%spin_dof, 3), source = ( 0.0_dp, 0.0_dp ))

            allocate(self%O_u_v_pad(self%n_pad_grid(1), &
                                self%n_pad_grid(2), &
                                self%n_pad_grid(3), &
                                init_state%spin_dof, 3), source = ( 0.0_dp, 0.0_dp ))

        end if
        if ( TIF_mask(3) ) then
            allocate(self%O_u_v2(init_state%n_x_grid(1), &
                                init_state%n_x_grid(2), &
                                init_state%n_x_grid(3), &
                                init_state%spin_dof), source = ( 0.0_dp, 0.0_dp ))

            allocate(self%O_u_v2_pad(self%n_pad_grid(1), &
                                self%n_pad_grid(2), &
                                self%n_pad_grid(3), &
                                init_state%spin_dof), source = ( 0.0_dp, 0.0_dp ))

        end if
        if ( TIF_mask(4) ) then
            allocate(self%O_u_s(init_state%n_x_grid(1), &
                                init_state%n_x_grid(2), &
                                init_state%n_x_grid(3), &
                                init_state%spin_dof, 3), source = ( 0.0_dp, 0.0_dp ))

            allocate(self%O_u_s_pad(self%n_pad_grid(1), &
                                self%n_pad_grid(2), &
                                self%n_pad_grid(3), &
                                init_state%spin_dof, 3), source = ( 0.0_dp, 0.0_dp ))

        end if
        if ( TIF_mask(5) ) then
            allocate(self%O_u_vds(init_state%n_x_grid(1), &
                                init_state%n_x_grid(2), &
                                init_state%n_x_grid(3), &
                                init_state%spin_dof), source = ( 0.0_dp, 0.0_dp ))

            allocate(self%O_u_vds_pad(self%n_pad_grid(1), &
                                self%n_pad_grid(2), &
                                self%n_pad_grid(3), &
                                init_state%spin_dof), source = ( 0.0_dp, 0.0_dp ))
        end if
        if ( TIF_mask(6) ) then
            allocate(self%O_u_vxs(init_state%n_x_grid(1), &
                                init_state%n_x_grid(2), &
                                init_state%n_x_grid(3), &
                                init_state%spin_dof, 3), source = ( 0.0_dp, 0.0_dp ))

            allocate(self%O_u_vxs_pad(self%n_pad_grid(1), &
                                self%n_pad_grid(2), &
                                self%n_pad_grid(3), &
                                init_state%spin_dof, 3), source = ( 0.0_dp, 0.0_dp ))
        end if
        if ( TIF_mask(7) ) then
            allocate(self%O_u_vivj(init_state%n_x_grid(1), &
                                init_state%n_x_grid(2), &
                                init_state%n_x_grid(3), &
                                init_state%spin_dof, 3, 3), source = ( 0.0_dp, 0.0_dp ))

            allocate(self%O_u_vivj_pad(self%n_pad_grid(1), &
                                self%n_pad_grid(2), &
                                self%n_pad_grid(3), &
                                init_state%spin_dof, 3, 3), source = ( 0.0_dp, 0.0_dp ))
        end if

    end subroutine

    subroutine TIF_calculator_bloch_type_init_initialize(self, init_state, TIF_mask)

        use O_u_bloch_calculator

        use FFT_util

        use timer_util

        implicit none

        class(TIF_calculator_bloch_t) :: self
        class(elec_state_bloch_t) :: init_state

        logical :: TIF_mask(:)

        integer :: s, d, i, j

        type(timer_t) :: timer

        ! Compute u_i
        call init_state%compute_u(self%u_i)

        ! Compute and pad O_u_i's
        if ( TIF_mask(1) ) then

            call compute_O_u_1(self%u_i, self%O_u_1, init_state)

            do s = 1, init_state%spin_dof
                call zero_pad_FFT_matrix(self%O_u_1(:, :, :, s), self%O_u_1_pad(:, :, :, s), &
                                            init_state%FFT_plans(1, :), self%backward_FFT_plan)
            end do

        end if
        if ( TIF_mask(2) ) then
            call compute_O_u_v(self%u_i, self%O_u_v, self%init_G_grid, init_state)

            do d = 1, 3
                do s = 1, init_state%spin_dof
                    call zero_pad_FFT_matrix(self%O_u_v(:, :, :, s, d), self%O_u_v_pad(:, :, :, s, d), &
                                            init_state%FFT_plans(1, :), self%backward_FFT_plan)
                end do
            end do

        end if
        if ( TIF_mask(3) ) then
            call compute_O_u_v2(self%u_i, self%O_u_v2, self%init_G_grid, init_state)

            do s = 1, init_state%spin_dof
                call zero_pad_FFT_matrix(self%O_u_v2(:, :, :, s), self%O_u_v2_pad(:, :, :, s), &
                                            init_state%FFT_plans(1, :), self%backward_FFT_plan)
            end do

        end if
        if ( TIF_mask(4) ) then
            call compute_O_u_s(self%u_i, self%O_u_s, init_state)

            do d = 1, 3
                do s = 1, init_state%spin_dof
                    call zero_pad_FFT_matrix(self%O_u_s(:, :, :, s, d), self%O_u_s_pad(:, :, :, s, d), &
                                            init_state%FFT_plans(1, :), self%backward_FFT_plan)
                end do
            end do

        end if
        if ( TIF_mask(5) ) then
            call compute_O_u_vds(self%u_i, self%O_u_vds, self%init_G_grid, init_state)

            do s = 1, init_state%spin_dof
                call zero_pad_FFT_matrix(self%O_u_vds(:, :, :, s), self%O_u_vds_pad(:, :, :, s), &
                                            init_state%FFT_plans(1, :), self%backward_FFT_plan)
            end do

        end if
        if ( TIF_mask(6) ) then
            call compute_O_u_vxs(self%u_i, self%O_u_vxs, self%init_G_grid, init_state)

            do d = 1, 3
                do s = 1, init_state%spin_dof
                    call zero_pad_FFT_matrix(self%O_u_vxs(:, :, :, s, d), self%O_u_vxs_pad(:, :, :, s, d), &
                                            init_state%FFT_plans(1, :), self%backward_FFT_plan)
                end do
            end do

        end if
        if ( TIF_mask(7) ) then
            call compute_O_u_vivj(self%u_i, self%O_u_vivj, self%init_G_grid, init_state)

            do j = 1, 3
                do i = j, 3 

                    do s = 1, init_state%spin_dof

                        call zero_pad_FFT_matrix(self%O_u_vivj(:, :, :, s, i, j), self%O_u_vivj_pad(:, :, :, s, i, j), &
                                                init_state%FFT_plans(1, :), self%backward_FFT_plan)

                    end do

                end do
            end do

            ! symmetrize
            do j = 1, 3 
                do i = j, 3

                    if ( i /= j ) then

                        self%O_u_vivj_pad(:, :, :, :, j, i) = self%O_u_vivj_pad(:, :, :, :, i, j)

                    end if

                end do
            end do

        end if

    end subroutine

    subroutine TIF_calculator_bloch_type_fin_initialize(self, fin_state)

        use timer_util

        use FFT_util

        implicit none

        class(TIF_calculator_bloch_t) :: self
        class(elec_state_bloch_t) :: fin_state

        integer :: s

        type(timer_t) :: timer

        ! Compute u_f
        call fin_state%compute_u(self%u_f)

        do s = 1, fin_state%spin_dof
            call zero_pad_FFT_matrix(self%u_f(:, :, :, s), self%u_f_pad(:, :, :, s), &
                                        fin_state%FFT_plans(1, :), self%backward_FFT_plan)
        end do

    end subroutine

    subroutine compute_TIF_helper_scalar(self, u_f_pad, O_u_pad, TIF)

        use FFT_util
        use timer_util

        implicit none

        class(TIF_calculator_bloch_t) :: self

        complex(dp) :: u_f_pad(:, :, :, :)
        complex(dp) :: O_u_pad(:, :, :, :)

        complex(dp) :: TIF(:)

        complex(dp), allocatable :: TIF_FT_pad(:, :, :)
        complex(dp), allocatable :: TIF_pad(:, :, :)

        type(timer_t) :: timer

        integer :: s

        allocate(TIF_FT_pad(self%n_FFT_grid(1), self%n_FFT_grid(2), self%n_FFT_grid(3)), source = ( 0.0_dp, 0.0_dp ))
        allocate(TIF_pad(self%n_FFT_grid(1), self%n_FFT_grid(2), self%n_FFT_grid(3)), source = ( 0.0_dp, 0.0_dp ))

        TIF_FT_pad = sum(conjg(u_f_pad)*O_u_pad, 4)

        ! KEY THAT THIS IS BACKWARDS
        call dfftw_execute_dft(self%backward_FFT_plan, TIF_FT_pad, TIF_pad) 

        ! normalize
        TIF_pad = (1.0_dp*product(self%n_pad_grid))**(-1)*TIF_pad

        ! pack in to list
        TIF = pack(TIF_pad, .TRUE.)

    end subroutine

    subroutine compute_TIF_helper_vector(self, u_f_pad, O_u_pad, init_state, fin_state, TIF)

        use FFT_util

        implicit none

        class(TIF_calculator_bloch_t) :: self
        class(elec_state_bloch_t) :: init_state, fin_state

        complex(dp) :: u_f_pad(:, :, :, :)
        complex(dp) :: O_u_pad(:, :, :, :, :)

        complex(dp) :: TIF(:, :)

        complex(dp), allocatable :: TIF_FT_pad(:, :, :, :)
        complex(dp), allocatable :: TIF_pad(:, :, :, :)

        integer :: d, s

        allocate(TIF_FT_pad(self%n_FFT_grid(1), self%n_FFT_grid(2), self%n_FFT_grid(3), 3), source = ( 0.0_dp, 0.0_dp ))
        allocate(TIF_pad(self%n_FFT_grid(1), self%n_FFT_grid(2), self%n_FFT_grid(3), 3), source = ( 0.0_dp, 0.0_dp ))

        do d = 1, 3

            do s = 1, min( size(u_f_pad, 4), size(O_u_pad, 4) )

                TIF_FT_pad(:, :, :, d) = TIF_FT_pad(:, :, :, d) + conjg(u_f_pad(:, :, :, s))*O_u_pad(:, :, :, s, d)

            end do

            ! KEY THAT THIS IS BACKWARDS
            call dfftw_execute_dft(self%backward_FFT_plan, TIF_FT_pad(:, :, :, d), TIF_pad(:, :, :, d)) 

            ! normalize
            TIF_pad(:, :, :, d) = (1.0_dp*product(self%n_pad_grid))**(-1)*TIF_pad(:, :, :, d)

            ! pack in to list
            TIF(:, d) = pack(TIF_pad(:, :, :, d), .TRUE.)

        end do

    end subroutine

    subroutine compute_TIF_helper_scalar_q0_limit(self, u_f_pad, O_u_pad, TIF)

        implicit none

        class(TIF_calculator_bloch_t) :: self

        complex(dp) :: u_f_pad(:, :, :, :)
        complex(dp) :: O_u_pad(:, :, :, :)

        complex(dp) :: TIF(:)

        TIF = (1.0_dp*product(self%n_pad_grid))**(-1)*sum( conjg(u_f_pad)*O_u_pad )

    end subroutine

    subroutine compute_TIF_helper_vector_q0_limit(self, u_f_pad, O_u_pad, init_state, fin_state, TIF)

        use elec_config_bloch_single_PW_type

        implicit none

        class(TIF_calculator_bloch_t) :: self
        class(elec_state_bloch_t) :: init_state, fin_state

        complex(dp) :: u_f_pad(:, :, :, :)
        complex(dp) :: O_u_pad(:, :, :, :, :)

        complex(dp) :: TIF(:, :)

        integer :: d

        do d = 1, 3
            TIF(1, d) = (1.0_dp*product(self%n_pad_grid))**(-1)*sum( conjg(u_f_pad)*O_u_pad(:, :, :, :, d) )
        end do

    end subroutine

    subroutine compute_TIF_helper_tensor_q0_limit(self, u_f_pad, O_u_pad, init_state, fin_state, TIF)

        ! use elec_config_bloch_single_PW_type

        implicit none

        class(TIF_calculator_bloch_t) :: self
        class(elec_state_bloch_t) :: init_state, fin_state

        complex(dp) :: u_f_pad(:, :, :, :)
        complex(dp) :: O_u_pad(:, :, :, :, :, :)

        complex(dp) :: TIF(:, :, :)

        integer :: i, j

        do j = 1, 3
            do i = 1, 3
                TIF(1, i, j) = (1.0_dp*product(self%n_pad_grid))**(-1)*sum( conjg(u_f_pad)*O_u_pad(:, :, :, :, i, j) )
            end do 
        end do

    end subroutine

    subroutine TIF_calculator_bloch_type_compute_TIF_1(self, TIF_1, init_state, fin_state, q0_limit)

        implicit none

        class(TIF_calculator_bloch_t) :: self
        complex(dp) :: TIF_1(:)
        class(elec_state_bloch_t) :: init_state, fin_state
        logical, optional :: q0_limit

        if ( q0_limit ) then

            ! For | I > != | F > this is zero by definition.
            ! TIF_1 = ( 0.0_dp, 0.0_dp )
            call compute_TIF_helper_scalar_q0_limit(self, self%u_f_pad, self%O_u_1_pad, TIF_1)

        else 

            call compute_TIF_helper_scalar(self, self%u_f_pad, self%O_u_1_pad, TIF_1)

        end if

    end subroutine

    subroutine TIF_calculator_bloch_type_compute_TIF_v(self, TIF_v, init_state, fin_state, q0_limit)

        implicit none

        class(TIF_calculator_bloch_t) :: self
        complex(dp) :: TIF_v(:, :)
        class(elec_state_bloch_t) :: init_state, fin_state
        logical, optional :: q0_limit

        if ( q0_limit ) then

            call compute_TIF_helper_vector_q0_limit(self, self%u_f_pad, self%O_u_v_pad, init_state, fin_state, TIF_v)

        else 

            call compute_TIF_helper_vector(self, self%u_f_pad, self%O_u_v_pad, init_state, fin_state, TIF_v)

        end if

    end subroutine

    subroutine TIF_calculator_bloch_type_compute_TIF_vivj(self, TIF_vivj, init_state, fin_state, q0_limit)

        implicit none

        class(TIF_calculator_bloch_t) :: self
        complex(dp) :: TIF_vivj(:, :, :)
        class(elec_state_bloch_t) :: init_state, fin_state
        logical, optional :: q0_limit

        if ( q0_limit ) then

            call compute_TIF_helper_tensor_q0_limit(self, self%u_f_pad, self%O_u_vivj_pad, init_state, fin_state, TIF_vivj)

        else 

            print*, 'ERROR: TIF_vivj for elec_state_bloch_t NOT implemented.'

            ! call compute_TIF_helper_vector(self, self%u_f_pad, self%O_u_v_pad, init_state, fin_state, TIF_v)

        end if

    end subroutine

    subroutine TIF_calculator_bloch_type_compute_TIF_v2(self, TIF_v2, init_state, fin_state, q0_limit)

        implicit none

        class(TIF_calculator_bloch_t) :: self
        complex(dp) :: TIF_v2(:)
        class(elec_state_bloch_t) :: init_state, fin_state
        logical, optional :: q0_limit

        if ( q0_limit ) then

            call compute_TIF_helper_scalar_q0_limit(self, self%u_f_pad, self%O_u_v2_pad, TIF_v2)

        else 

            call compute_TIF_helper_scalar(self, self%u_f_pad, self%O_u_v2_pad, TIF_v2)

        end if

    end subroutine

    subroutine TIF_calculator_bloch_type_compute_TIF_s(self, TIF_s, init_state, fin_state, q0_limit)

        implicit none

        class(TIF_calculator_bloch_t) :: self
        complex(dp) :: TIF_s(:, :)
        class(elec_state_bloch_t) :: init_state, fin_state
        logical, optional :: q0_limit

        if ( q0_limit ) then

            call compute_TIF_helper_vector_q0_limit(self, self%u_f_pad, self%O_u_s_pad, init_state, fin_state, TIF_s)

        else 

            call compute_TIF_helper_vector(self, self%u_f_pad, self%O_u_s_pad, init_state, fin_state, TIF_s)

        end if

    end subroutine

    subroutine TIF_calculator_bloch_type_compute_TIF_vds(self, TIF_vds, init_state, fin_state, q0_limit)

        implicit none

        class(TIF_calculator_bloch_t) :: self
        complex(dp) :: TIF_vds(:)
        class(elec_state_bloch_t) :: init_state, fin_state
        logical, optional :: q0_limit

        if ( q0_limit ) then

            call compute_TIF_helper_scalar_q0_limit(self, self%u_f_pad, self%O_u_vds_pad, TIF_vds)

        else 

            call compute_TIF_helper_scalar(self, self%u_f_pad, self%O_u_vds_pad, TIF_vds)

        end if

    end subroutine

    subroutine TIF_calculator_bloch_type_compute_TIF_vxs(self, TIF_vxs, init_state, fin_state, q0_limit)

        implicit none

        class(TIF_calculator_bloch_t) :: self
        complex(dp) :: TIF_vxs(:, :)
        class(elec_state_bloch_t) :: init_state, fin_state
        logical, optional :: q0_limit

        if ( q0_limit ) then

            call compute_TIF_helper_vector_q0_limit(self, self%u_f_pad, self%O_u_vxs_pad, init_state, fin_state, TIF_vxs)

        else 

            call compute_TIF_helper_vector(self, self%u_f_pad, self%O_u_vxs_pad, init_state, fin_state, TIF_vxs)

        end if

    end subroutine

end module
