! module r2_O_R_atomic_calculator
!
!     use prec_util, only: dp
!
!     use elec_state_atomic_type 
!
!     implicit none
!
! contains
!
!     subroutine compute_r2_O_R_1(r2_O_R_1, init_state)
!
!         implicit none
!
!         complex(dp) :: r2_O_R_1(:)
!         real(dp) :: R_wf(size(r2_O_R_1))
!         class(elec_state_atomic_t) :: init_state
!
!         call init_state%compute_radial_wf(R_wf)
!
!         r2_O_R_1 = init_state%sph_x_list(:, 1)**2*R_wf
!
!     end subroutine
!
!     subroutine compute_r2_O_R_v2(r2_O_R_v2, init_state)
!
!         use constants_util
!
!         use math_util
!
!         use elec_state_atomic_STO_basis_type
!
!         implicit none
!
!         complex(dp) :: r2_O_R_v2(:)
!         class(elec_state_atomic_t) :: init_state
!         real(dp) :: R_wf_j(size(r2_O_R_v2))
!
!         integer :: j
!
!         integer :: n, l
!
!         real(dp) :: coeff, zeta
!
!         l = init_state%l
!
!         ! make sure that the initial state is in the STO basis
!         select type( init_state )
!         class is ( elec_state_atomic_STO_basis_t )
!
!             r2_O_r_v2 = ( 0.0_dp, 0.0_dp )
!
!             do j = 1, size(init_state%Cnl_list)
!
!                 n = init_state%nl_list(j)
!                 coeff = init_state%Cnl_list(j)
!                 zeta = init_state%Zl_list(j)/a0
!
!                 R_wf_j = STO_radial(init_state%sph_x_list(:, 1), &
!                     init_state%nl_list(j), &
!                     init_state%norml_list(j), &
!                     init_state%Zl_list(j))
!
!                 r2_O_R_v2 = r2_O_R_v2 + &
!                     (m_elec)**(-2)*coeff*R_wf_j*(&
!                         ( n*(n - 1) - l*(l + 1) ) &
!                          - 2*n*zeta*init_state%sph_x_list(:, 1) &
!                          + zeta**2*init_state%sph_x_list(:, 1)**2 )
!
!             end do
!             
!         end select
!
!     end subroutine
!
!     subroutine compute_r2_O_R_v(r2_O_R_v, init_state)
!
!         implicit none
!
!         complex(dp) :: r2_O_R_v(:, :)
!         class(elec_state_atomic_t) :: init_state
!
!         print*, 'ERROR: v of STO wave functions not implemented yet.'
!
!     end subroutine
!
!     subroutine compute_r2_O_R_s(r2_O_R_s, init_state)
!
!         implicit none
!
!         complex(dp) :: r2_O_R_s(:, :)
!         class(elec_state_atomic_t) :: init_state
!
!         r2_O_R_s = ( 0.0_dp, 0.0_dp )
!
!     end subroutine
!
!     subroutine compute_r2_O_R_vxs(r2_O_R_vxs, init_state)
!
!         implicit none
!
!         complex(dp) :: r2_O_R_vxs(:, :)
!         class(elec_state_atomic_t) :: init_state
!
!         r2_O_R_vxs = ( 0.0_dp, 0.0_dp )
!
!     end subroutine
!
!     subroutine compute_r2_O_R_vds(r2_O_R_vds, init_state)
!
!         implicit none
!
!         complex(dp) :: r2_O_R_vds(:)
!         class(elec_state_atomic_t) :: init_state
!
!         r2_O_R_vds = ( 0.0_dp, 0.0_dp )
!
!     end subroutine
!
! end module
