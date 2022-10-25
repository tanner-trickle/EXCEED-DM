module elec_state_bloch_PW_basis_type

    use prec_util, only: dp

    use hdf5_utils

    use elec_state_bloch_type

    implicit none

    type, extends(elec_state_bloch_t) :: elec_state_bloch_PW_basis_t

        integer(HID_T) :: hdf5_file_id

        character(len=512) :: hdf5_data_filename 
        character(len=512) :: hdf5_G_list_red_path
        character(len=512) :: hdf5_u_FT_r_path ! Path in hdf5 file to real component of :math:`\widetilde{u}`
        character(len=512) :: hdf5_u_FT_c_path ! Path in hdf5 file to complex component of :math:`\widetilde{u}`

        contains

            procedure :: compute_u => elec_state_bloch_PW_basis_type_compute_u

    end type

contains

    subroutine elec_state_bloch_PW_basis_type_compute_u(self, u)

        use hdf5_utils

        use constants_util
        use timer_util
        use FFT_util

        implicit none

        class(elec_state_bloch_PW_basis_t) :: self
        complex(dp) :: u(:, :, :, :)

        integer(HID_T) :: file_id

        complex(dp), allocatable:: u_wf_FT(:, :)
        real(dp), allocatable:: u_wf_FT_r(:, :)
        real(dp), allocatable:: u_wf_FT_c(:, :)
        complex(dp), allocatable :: u_wf_FT_expand(:, :, :, :)

        integer :: g, s, i
        integer :: FFT_idx(3)

        integer :: dims(6)

        integer, allocatable :: G_list_red(:, :)

        type(timer_t) :: timer

        u = (0.0_dp, 0.0_dp)

        call timer%start()

        ! load G_list_red
        call hdf_get_dims(self%hdf5_file_id, self%hdf5_G_list_red_path, dims)
        allocate(G_list_red(dims(1), dims(2)), source = 0)
        call hdf_read_dataset(self%hdf5_file_id, self%hdf5_G_list_red_path, G_list_red)

        allocate(u_wf_FT(dims(1), self%spin_dof), source = ( 0.0_dp, 0.0_dp ))
        allocate(u_wf_FT_r(dims(1), self%spin_dof), source = 0.0_dp)
        allocate(u_wf_FT_c(dims(1), self%spin_dof), source = 0.0_dp)

        call timer%end()

        ! print*, 'time to load G_list_red = ', timer%pretty_dt_str()

        call timer%start()

        ! load u_wf_FT from hdf5 file
        call hdf_read_dataset(self%hdf5_file_id, self%hdf5_u_FT_r_path, u_wf_FT_r)
        call hdf_read_dataset(self%hdf5_file_id, self%hdf5_u_FT_c_path, u_wf_FT_c)

        call timer%end()

        ! print*, 'time to load wf coeff = ', timer%pretty_dt_str()

        u_wf_FT = u_wf_FT_r + ii*u_wf_FT_c

        ! expand to the FFT size
        allocate(u_wf_FT_expand(self%n_x_grid(1),&
                                self%n_x_grid(2),&
                                self%n_x_grid(3),&
                                self%spin_dof), source = ( 0.0_dp, 0.0_dp ) )

        do g = 1, size(G_list_red, 1)
            FFT_idx = G_red_to_FFT_idx(G_list_red(g, :), self%n_x_grid)
            u_wf_FT_expand(FFT_idx(1), FFT_idx(2), FFT_idx(3), :) = u_wf_FT(g, :)
        end do

        ! compute FFT's of each spin component
        do s = 1, self%spin_dof
            call dfftw_execute_dft(self%FFT_plans(2, :), u_wf_FT_expand(:, :, :, s), u(:, :, :, s))
        end do

    end subroutine

end module
