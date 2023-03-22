module TIF_calculator_atomic_type

    use prec_util, only: dp

    use exdm_inputs_type

    use elec_state_atomic_type

    implicit none

    type :: TIF_calculator_atomic_t

        real(dp), allocatable :: R_F(:)
        real(dp), allocatable :: R_I(:)
        real(dp), allocatable :: dR_I(:)
        real(dp), allocatable :: d2R_I(:)

        real(dp), allocatable :: sph_x_list(:, :)
        real(dp), allocatable :: jac_r_list(:)

        integer :: lp_min, lp_max
        integer :: mp_min, mp_max

        integer :: n_q

        contains

            procedure :: setup => TIF_calculator_atomic_type_setup
            procedure :: init_initialize => TIF_calculator_atomic_type_init_initialize
            procedure :: fin_initialize => TIF_calculator_atomic_type_fin_initialize
            procedure :: compute_TIF_1 => TIF_calculator_atomic_type_compute_TIF_1
            procedure :: compute_TIF_v => TIF_calculator_atomic_type_compute_TIF_v
            procedure :: compute_TIF_v2 => TIF_calculator_atomic_type_compute_TIF_v2
            procedure :: compute_TIF_s => TIF_calculator_atomic_type_compute_TIF_s
            procedure :: compute_TIF_vds => TIF_calculator_atomic_type_compute_TIF_vds
            procedure :: compute_TIF_vxs => TIF_calculator_atomic_type_compute_TIF_vxs
            procedure :: compute_TIF_vivj => TIF_calculator_atomic_type_compute_TIF_vivj

    end type

contains

    subroutine TIF_calculator_atomic_type_setup(self, exdm_inputs, TIF_mask, init_state, fin_state, q0_limit)

        implicit none

        class(TIF_calculator_atomic_t) :: self

        type(exdm_inputs_t) :: exdm_inputs
        logical :: TIF_mask(:)

        class(elec_state_atomic_t) :: init_state
        class(elec_state_atomic_t) :: fin_state

        logical, optional :: q0_limit

        integer :: n_r

        integer :: r

        ! initialize 
        n_r = exdm_inputs%numerics_TIF_calculator_atomic%n_r 

        ! since we're focusing on absorption right now...
        self%n_q = 1

        ! create integration r grid
        allocate(self%sph_x_list(n_r, 3), source = 0.0_dp)
        allocate(self%jac_r_list(n_r), source = 1.0_dp)

        call exdm_inputs%numerics_TIF_calculator_atomic%generate_integration_grid(&
            self%sph_x_list(:, 1), self%jac_r_list)

        ! allocate wave functions to be computed
        allocate(self%R_F(n_r), source = 0.0_dp)

        allocate(self%R_I(n_r), source = 0.0_dp)
        allocate(self%dR_I(n_r), source = 0.0_dp)
        allocate(self%d2R_I(n_r), source = 0.0_dp)

    end subroutine

    subroutine TIF_calculator_atomic_type_init_initialize(self, init_state, TIF_mask)

        use elec_state_atomic_STO_basis_type

        implicit none

        class(TIF_calculator_atomic_t) :: self
        class(elec_state_atomic_t) :: init_state

        logical :: TIF_mask(:)

        init_state%sph_x_list = self%sph_x_list

        ! TODO: add spin dependent scenario
        if ( TIF_mask(1) ) then
            self%lp_min = init_state%l
            self%lp_max = init_state%l
            self%mp_min = init_state%m
            self%mp_min = init_state%m
        end if
        if ( TIF_mask(2) ) then
            self%lp_min = init_state%l - 1
            self%lp_max = init_state%l + 1
            self%mp_min = init_state%m - 1
            self%mp_max = init_state%m + 1
        end if
        if ( TIF_mask(7) ) then
            self%lp_min = init_state%l - 2
            self%lp_max = init_state%l + 2
            self%mp_min = init_state%m - 2
            self%mp_max = init_state%m + 2
        end if

        ! R_I
        call init_state%compute_radial_wf(self%R_I)

        ! for now, compute all derivatives, to optimize only compute those needed based on TIF_mask
        select type ( init_state )
        class is ( elec_state_atomic_STO_basis_t)

            ! dR_I
            call init_state%compute_d_radial_wf(self%dR_I)

            ! d2R_I
            call init_state%compute_d2_radial_wf(self%d2R_I)

        end select

    end subroutine

    subroutine TIF_calculator_atomic_type_fin_initialize(self, fin_state)

        implicit none

        class(TIF_calculator_atomic_t) :: self
        class(elec_state_atomic_t) :: fin_state

        fin_state%sph_x_list = self%sph_x_list

    end subroutine

    subroutine TIF_calculator_atomic_type_compute_TIF_1(self, TIF_1, init_state, fin_state, q0_limit)

        implicit none

        class(TIF_calculator_atomic_t) :: self
        complex(dp) :: TIF_1(:)
        class(elec_state_atomic_t) :: init_state, fin_state
        logical, optional :: q0_limit

        if ( q0_limit ) then 

            ! state orthonormalization
            TIF_1 = (0.0_dp, 0.0_dp)

        else 

            print*, 'ERROR: TIF_1, in q/= 0 limit, not implemented for atomic states, setting TIF_1 = 0.'
            print*
            TIF_1 = (0.0_dp, 0.0_dp)

        end if

    end subroutine

    subroutine TIF_calculator_atomic_type_compute_TIF_v(self, TIF_v, init_state, fin_state, q0_limit)

        use constants_util
        use derivative_VSH_coefficients

        implicit none

        class(TIF_calculator_atomic_t) :: self
        complex(dp) :: TIF_v(:, :)
        class(elec_state_atomic_t) :: init_state, fin_state
        logical, optional :: q0_limit

        integer :: i
        integer :: k, q, l, lp, m, mp

        complex(dp) :: c_i
        complex(dp) :: d_i

        real(dp) :: I01, I10

        if ( q0_limit ) then 

            l = init_state%l
            m = init_state%m

            lp = fin_state%l
            mp = fin_state%m

            k = lp - l
            q = mp - m

            if ( ( abs(k) <= 1 ) .and. ( abs(q) <= 1 ) ) then

                I01 = sum( self%jac_r_list * self%sph_x_list(:, 1)**2 * self%R_F * self%dR_I )
                I10 = sum( self%jac_r_list * self%sph_x_list(:, 1) * self%R_F * self%R_I )

                do i = 1, 3

                    c_i = derivative_VSH_c(i, k, q, l, m)
                    d_i = derivative_VSH_d(i, k, q, l, m)

                    TIF_v(1, i) = -(ii/m_elec)*(c_i*I01 + d_i*I10)

                end do

            else 

                TIF_v = (0.0_dp, 0.0_dp)

            end if

        else 

            print*, 'ERROR: TIF_v, in q/=0 limit, not implemented for atomic states, setting TIF_v = 0.'
            print*
            TIF_v = (0.0_dp, 0.0_dp)

        end if

    end subroutine

    subroutine TIF_calculator_atomic_type_compute_TIF_vivj(self, TIF_vivj, init_state, fin_state, q0_limit)

        use constants_util
        use derivative_VSH_coefficients

        implicit none

        class(TIF_calculator_atomic_t) :: self
        complex(dp) :: TIF_vivj(:, :, :)
        complex(dp) :: tr_TIF_vivj
        complex(dp) :: TIF_v2
        class(elec_state_atomic_t) :: init_state, fin_state
        logical, optional :: q0_limit

        real(dp) :: I02, I11, I20

        complex(dp) :: A_mat(3, 3), B_mat(3, 3), C_mat(3, 3)

        integer :: i, j, k, q

        complex(dp) :: c_i, c_j, d_i, d_j

        integer :: l, lp, m, mp

        if ( q0_limit ) then 

            l = init_state%l
            m = init_state%m

            lp = fin_state%l
            mp = fin_state%m

            I02 = sum( self%jac_r_list * self%sph_x_list(:, 1)**2 * self%R_F * self%d2R_I )
            I11 = sum( self%jac_r_list * self%sph_x_list(:, 1) * self%R_F * self%dR_I )
            I20 = sum( self%jac_r_list * self%R_F * self%R_I )

            A_mat = (0.0_dp, 0.0_dp)
            B_mat = (0.0_dp, 0.0_dp)
            C_mat = (0.0_dp, 0.0_dp)

            do i = 1, 3
                do j = 1, 3

                    do k = -1, 1
                        do q = -1, 1

                            if ( ( l + k >= 0 ) .and. &
                                 ( abs(m + q) <= l + k ) .and. & 
                                 ( abs(lp - l - k) <= 1 ) .and. &
                                 ( abs(mp - m - q) <= 1 ) ) then

                                c_i = derivative_VSH_c(i, lp - l - k, mp - m - q, l + k, m + q)
                                d_i = derivative_VSH_d(i, lp - l - k, mp - m - q, l + k, m + q)

                                c_j = derivative_VSH_c(j, k, q, l, m)
                                d_j = derivative_VSH_d(j, k, q, l, m)

                                A_mat(i, j) = A_mat(i, j) + c_i*c_j 

                                B_mat(i, j) = B_mat(i, j) + d_i*c_j + d_j*c_i 

                                C_mat(i, j) = C_mat(i, j) + d_j*(d_i - c_i)

                            end if

                        end do
                    end do

                    TIF_vivj(1, i, j) = - (m_elec)**(-2)*(A_mat(i, j)*I02 + B_mat(i, j)*I11 + C_mat(i, j)*I20 )

                end do 
            end do

        else 

            print*, 'ERROR: TIF_vivj, in q/=0 limit, not implemented for atomic states, setting TIF_vivj = 0.'
            print*
            TIF_vivj = (0.0_dp, 0.0_dp)

        end if

    end subroutine

    subroutine TIF_calculator_atomic_type_compute_TIF_v2(self, TIF_v2, init_state, fin_state, q0_limit)

        use constants_util

        implicit none

        class(TIF_calculator_atomic_t) :: self
        complex(dp) :: TIF_v2(:)
        class(elec_state_atomic_t) :: init_state, fin_state
        logical, optional :: q0_limit

        if ( q0_limit ) then 

            TIF_v2 = -(m_elec)**(-2)*(&
                sum( self%jac_r_list * self%sph_x_list(:, 1)**2 * self%R_F * self%d2R_I ) &
              + 2*sum( self%jac_r_list * self%sph_x_list(:, 1) * self%R_F * self%dR_I ) &
              - (init_state%l*(init_state%l + 1))*sum( self%jac_r_list * self%R_F * self%R_I ) &
                )

        else

            print*, 'ERROR: TIF_v2, in q/=0 limit, not implemented for atomic states, setting TIF_v = 0.'
            print*
            TIF_v2 = (0.0_dp, 0.0_dp)

        end if

    end subroutine

    subroutine TIF_calculator_atomic_type_compute_TIF_s(self, TIF_s, init_state, fin_state, q0_limit)

        implicit none

        class(TIF_calculator_atomic_t) :: self
        complex(dp) :: TIF_s(:, :)
        class(elec_state_atomic_t) :: init_state, fin_state
        logical, optional :: q0_limit

        if ( q0_limit ) then 

            print*, 'ERROR: TIF_s, in q=0 limit, not implemented for atomic states, setting TIF_v = 0.'
            TIF_s = (0.0_dp, 0.0_dp)

        else 

            print*, 'ERROR: TIF_s, in q/=0 limit, not implemented for atomic states, setting TIF_v = 0.'
            TIF_s = (0.0_dp, 0.0_dp)

        end if

    end subroutine

    subroutine TIF_calculator_atomic_type_compute_TIF_vds(self, TIF_vds, init_state, fin_state, q0_limit)

        implicit none

        class(TIF_calculator_atomic_t) :: self
        complex(dp) :: TIF_vds(:)
        class(elec_state_atomic_t) :: init_state, fin_state
        logical, optional :: q0_limit

        if ( q0_limit ) then 

            print*, 'ERROR: TIF_vds, in q=0 limit, not implemented for atomic states, setting TIF_v = 0.'
            TIF_vds = (0.0_dp, 0.0_dp)

        else 

            print*, 'ERROR: TIF_vds, in q/=0 limit, not implemented for atomic states, setting TIF_v = 0.'
            TIF_vds = (0.0_dp, 0.0_dp)

        end if

    end subroutine

    subroutine TIF_calculator_atomic_type_compute_TIF_vxs(self, TIF_vxs, init_state, fin_state, q0_limit)

        implicit none

        class(TIF_calculator_atomic_t) :: self
        complex(dp) :: TIF_vxs(:, :)
        class(elec_state_atomic_t) :: init_state, fin_state
        logical, optional :: q0_limit

        if ( q0_limit ) then 

            print*, 'ERROR: TIF_vxs, in q=0 limit, not implemented for atomic states, setting TIF_v = 0.'
            TIF_vxs = (0.0_dp, 0.0_dp)

        else 

            print*, 'ERROR: TIF_vxs, in q/=0 limit, not implemented for atomic states, setting TIF_v = 0.'
            TIF_vxs = (0.0_dp, 0.0_dp)

        end if

    end subroutine

end module
