module PiIF_calculator_type

    use prec_util, only: dp
    use constants_util

    use abs_physics_functions

    implicit none

    type :: PiIF_calculator_t

        integer :: n = 4
        logical :: mask(4)

        complex(dp), allocatable :: Pi_1_1_mat(:, :, :, :) ! Dim : [3, 3, n_mX, n_widths]
        complex(dp), allocatable :: Pi_v_v(:, :, :, :)
        complex(dp), allocatable :: Pi_v2_v2(:, :) ! Dim : [n_mX, n_widths]
        complex(dp), allocatable :: Pi_vds_vds(:, :)

        contains

            procedure :: setup => setup_PiIF_calculator
            procedure :: compute_all => compute_all_PiIF_calculator

    end type

contains

    subroutine PiIF_mask_to_TIF_mask(Pi_mask, TIF_mask)

        implicit none

        logical :: Pi_mask(:)
        logical :: TIF_mask(:)

        TIF_mask = .FALSE.

        if ( Pi_mask(1) ) then
            TIF_mask(2) = .TRUE.
        end if

        if ( Pi_mask(2) ) then
            TIF_mask(2) = .TRUE.
        end if

        if ( Pi_mask(3) ) then
            TIF_mask(3) = .TRUE.
        end if

        if ( Pi_mask(4) ) then
            TIF_mask(5) = .TRUE.
        end if

    end subroutine

    subroutine particle_type_to_PiIF_mask(particle_type, spin_dof, mask, n)

        implicit none

        character(len=*) :: particle_type
        integer :: spin_dof
        integer :: n
        logical :: mask(n)

        mask = .FALSE.

        select case( trim(adjustl(particle_type)) )

        case( 'vector' )

            mask(1) = .TRUE.

        case( 'ps' )

            if ( spin_dof == 1 ) then

                mask(1) = .TRUE.

            else

                mask(4) = .TRUE.

            end if

        case( 'scalar' )

            mask(3) = .TRUE.

        end select

    end subroutine

    subroutine compute_all_PiIF_calculator(self, TIF_calculator, omega_IF, mX, widths,&
            jac_k, pc_vol, spin_dof, smear_type)

        use TIF_calculator_type

        implicit none

        class(PiIF_calculator_t) :: self
        class(TIF_calculator_t) :: TIF_calculator

        real(dp) :: omega_IF
        real(dp) :: jac_k
        real(dp) :: pc_vol
        integer :: spin_dof
        
        real(dp) :: mX(:)
        real(dp) :: widths(:, :)

        character(len=*) :: smear_type

        if ( self%mask(1) ) then
            call compute_Pi_1_1_mat(self%Pi_1_1_mat, TIF_calculator%TIF_v, &
                omega_IF, mX, widths, &
                jac_k, pc_vol, spin_dof, smear_type)
        end if
        if ( self%mask(2) ) then
            call compute_Pi_v_v(self%Pi_v_v, TIF_calculator%TIF_v, &
                omega_IF, mX, widths, &
                jac_k, pc_vol, spin_dof, smear_type)
        end if
        if ( self%mask(3) ) then
            call compute_Pi_v2_v2(self%Pi_v2_v2, TIF_calculator%TIF_v2, &
                omega_IF, mX, widths, &
                jac_k, pc_vol, spin_dof, smear_type)
        end if
        if ( self%mask(4) ) then
            call compute_Pi_vds_vds(self%Pi_vds_vds, TIF_calculator%TIF_vds, &
                omega_IF, mX, widths, &
                jac_k, pc_vol, spin_dof, smear_type)
        end if

    end subroutine

    subroutine setup_PiIF_calculator(self, n_mX, n_widths)

        implicit none

        class(PiIF_calculator_t) :: self

        integer :: n_mX, n_widths

        if ( self%mask(1) ) then
            allocate(self%Pi_1_1_mat(3, 3, n_mX, n_widths), &
                     source = (0.0_dp, 0.0_dp))
        end if
        if ( self%mask(2) ) then
            allocate(self%Pi_v_v(3, 3, n_mX, n_widths), &
                     source = (0.0_dp, 0.0_dp))
        end if
        if ( self%mask(3) ) then
            allocate(self%Pi_v2_v2(n_mX, n_widths), &
                     source = (0.0_dp, 0.0_dp))
        end if
        if ( self%mask(4) ) then
            allocate(self%Pi_vds_vds(n_mX, n_widths), &
                     source = (0.0_dp, 0.0_dp))
        end if

    end subroutine

    subroutine reduce_PiIF(root_process, PiIF, sum_Pi)

        use mpi

        implicit none

        integer :: root_process

        class(PiIF_calculator_t) :: PiIF
        class(PiIF_calculator_t) :: sum_Pi

        integer :: i, err

        if ( PiIF%mask(1) ) then
            call MPI_Reduce(PiIF%Pi_1_1_mat, &
                            sum_Pi%Pi_1_1_mat, &
                            size(PiIF%Pi_1_1_mat), &
                            MPI_DOUBLE_COMPLEX, &
                            MPI_SUM, &
                            root_process, &
                            MPI_COMM_WORLD, &
                            err)
        end if

        if ( PiIF%mask(2) ) then
            call MPI_Reduce(PiIF%Pi_v_v, &
                            sum_Pi%Pi_v_v, &
                            size(PiIF%Pi_v_v), &
                            MPI_DOUBLE_COMPLEX, &
                            MPI_SUM, &
                            root_process, &
                            MPI_COMM_WORLD, &
                            err)
        end if

        if ( PiIF%mask(3) ) then
            call MPI_Reduce(PiIF%Pi_v2_v2, &
                            sum_Pi%Pi_v2_v2, &
                            size(PiIF%Pi_v2_v2), &
                            MPI_DOUBLE_COMPLEX, &
                            MPI_SUM, &
                            root_process, &
                            MPI_COMM_WORLD, &
                            err)
        end if

        if ( PiIF%mask(4) ) then
            call MPI_Reduce(PiIF%Pi_vds_vds, &
                            sum_Pi%Pi_vds_vds, &
                            size(PiIF%Pi_vds_vds), &
                            MPI_DOUBLE_COMPLEX, &
                            MPI_SUM, &
                            root_process, &
                            MPI_COMM_WORLD, &
                            err)
        end if

    end subroutine

    subroutine compute_Pi_1_1_mat(Pi_1_1, TIF_v, &
            omega_IF, mX, widths, &
            jac_k, pc_vol, spin_dof, smear_type)

        implicit none

        complex(dp) :: Pi_1_1(:, :, :, :)
        complex(dp) :: TIF_v(:, :)
        real(dp) :: omega_IF
        real(dp) :: mX(:)
        real(dp) :: widths(:, :)
        real(dp) :: jac_k
        real(dp) :: pc_vol
        integer :: spin_dof

        integer :: i, j, m

        real(dp) :: delta

        complex(dp) :: outer_TIF_v(3, 3)

        character(len=*) :: smear_type

        do j = 1, 3
            do i = 1, 3

                outer_TIF_v(i, j) = conjg(TIF_v(1, i))*TIF_v(1, j)

            end do
        end do

        do i = 1, size(widths, 1)
            do m = 1, size(mX)

                delta = width_func(mX(m), widths(i, 1), widths(i, 2), widths(i, 3))
                
                Pi_1_1(:, :, m, i) = Pi_1_1(:, :, m, i) + &
                    (2.0_dp/spin_dof)*jac_k*pc_vol**(-1)*&
                    (m_elec/omega_IF)**2*green_func(mX(m), omega_IF, delta, smear_type)*&
                    outer_TIF_v

            end do
        end do

    end subroutine

    subroutine compute_Pi_v_v(Pi_v_v, TIF_v, &
            omega_IF, mX, widths, &
            jac_k, pc_vol, spin_dof, smear_type)

        implicit none

        complex(dp) :: Pi_v_v(:, :, :, :)
        complex(dp) :: TIF_v(:, :)
        real(dp) :: omega_IF
        real(dp) :: mX(:)
        real(dp) :: widths(:, :)
        real(dp) :: jac_k
        real(dp) :: pc_vol
        integer :: spin_dof

        integer :: i, j, m

        real(dp) :: delta

        complex(dp) :: outer_TIF_v(3, 3)

        character(len=*) :: smear_type

        do j = 1, 3
            do i = 1, 3

                outer_TIF_v(i, j) = conjg(TIF_v(1, i))*TIF_v(1, j)

            end do
        end do

        do i = 1, size(widths, 1)
            do m = 1, size(mX)

                delta = width_func(mX(m), widths(i, 1), widths(i, 2), widths(i, 3))
                
                Pi_v_v(:, :, m, i) = Pi_v_v(:, :, m, i) + &
                    (2.0_dp/spin_dof)*jac_k*pc_vol**(-1)*&
                    Pi_scaling_func(mX(m), omega_IF)*green_func(mX(m), omega_IF, delta, smear_type)*&
                    outer_TIF_v

            end do
        end do

    end subroutine

    subroutine compute_Pi_v2_v2(Pi_v2_v2, TIF_v2, &
            omega_IF, mX, widths, &
            jac_k, pc_vol, spin_dof, smear_type)

        implicit none

        complex(dp) :: Pi_v2_v2(:, :)
        complex(dp) :: TIF_v2(:)
        real(dp) :: omega_IF
        real(dp) :: mX(:)
        real(dp) :: widths(:, :)
        real(dp) :: jac_k
        real(dp) :: pc_vol
        integer :: spin_dof

        integer :: i, m

        real(dp) :: delta

        character(len=*) :: smear_type

        do i = 1, size(widths, 1)
            do m = 1, size(mX)

                delta = width_func(mX(m), widths(i, 1), widths(i, 2), widths(i, 3))
                
                Pi_v2_v2(m, i) = Pi_v2_v2(m, i) + &
                    (0.25_dp)*(2.0_dp/spin_dof)*jac_k*pc_vol**(-1)*&
                    Pi_scaling_func(mX(m), omega_IF)*green_func(mX(m), omega_IF, delta, smear_type)*&
                    conjg(TIF_v2(1))*TIF_v2(1)

            end do
        end do

    end subroutine

    subroutine compute_Pi_vds_vds(Pi_vds_vds, TIF_vds, &
            omega_IF, mX, widths, &
            jac_k, pc_vol, spin_dof, smear_type)

        implicit none

        complex(dp) :: Pi_vds_vds(:, :)
        complex(dp) :: TIF_vds(:)
        real(dp) :: omega_IF
        real(dp) :: mX(:)
        real(dp) :: widths(:, :)
        real(dp) :: jac_k
        real(dp) :: pc_vol
        integer :: spin_dof

        integer :: i, m

        real(dp) :: delta

        character(len=*) :: smear_type

        do i = 1, size(widths, 1)
            do m = 1, size(mX)

                delta = width_func(mX(m), widths(i, 1), widths(i, 2), widths(i, 3))
                
                Pi_vds_vds(m, i) = Pi_vds_vds(m, i) + &
                    (0.25_dp)*(2.0_dp/spin_dof)*jac_k*pc_vol**(-1)*&
                    Pi_scaling_func(mX(m), omega_IF)*green_func(mX(m), omega_IF, delta, smear_type)*&
                    conjg(TIF_vds(1))*TIF_vds(1)

            end do
        end do

    end subroutine

end module
