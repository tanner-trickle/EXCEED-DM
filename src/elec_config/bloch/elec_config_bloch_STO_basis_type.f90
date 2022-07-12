module elec_config_bloch_STO_basis_type

    use prec_util, only: dp

    use elec_state_bloch_STO_basis_type

    implicit none

    type :: elec_config_bloch_STO_basis_t

        type(elec_state_bloch_STO_basis_t), allocatable :: states(:)

        real(dp), allocatable :: jac_list(:)
        real(dp), allocatable :: Zeff_list(:)
        integer, allocatable :: k_id_list(:)

        contains

            procedure :: initialize => elec_config_bloch_STO_basis_type_initialize

    end type

contains

    subroutine elec_config_bloch_STO_basis_type_initialize(self, exdm_inputs, type)

        use hdf5_utils
        use exdm_inputs_type

        use FFT_util

        implicit none

        class(elec_config_bloch_STO_basis_t) :: self
        type(exdm_inputs_t) :: exdm_inputs
        character(len=*) :: type

        integer(HID_t) :: file_id

        integer :: n_states
        integer :: dims(6)

        integer :: rank

        integer :: i, n, r, r1, r2, r3

        character(len=512) :: data_path

        real(dp), allocatable :: energy_list(:)
        integer, allocatable :: i_list(:)
        real(dp), allocatable :: k_vec_red_list(:, :)
        integer, allocatable :: nlm_list(:, :)
        integer, allocatable :: nj_list(:)
        real(dp), allocatable :: coeff_list(:, :, :) 
        real(dp), allocatable :: eq_pos_red_list(:, :)

        integer :: n_x_grid(3)
        integer :: n_r_vec_grid(3)
        integer :: FFT_plans(2, 8)
        integer :: spin_dof

        call hdf_open_file(file_id, exdm_inputs%elec_config_input%filename, &
            status='OLD', action='READ')

        ! make sure the dataset exists
        data_path = "elec_states/"//trim(adjustl(type))//"/bloch/STO_basis/"
        if ( hdf_exists(file_id, trim(adjustl(data_path))) ) then

            ! state info
            call hdf_get_dims(file_id, trim(adjustl(data_path))//"state_info/energy_list", dims)
            n_states = dims(1)

            allocate(self%states(n_states))
            allocate(self%jac_list(n_states), source = 1.0_dp)
            allocate(self%Zeff_list(n_states), source = 1.0_dp)
            allocate(self%k_id_list(n_states), source = 1)

            allocate(energy_list(n_states), source = 0.0_dp)
            call hdf_read_dataset(file_id, trim(adjustl(data_path))//"state_info/energy_list", energy_list)

            allocate(i_list(n_states), source = 1)
            call hdf_read_dataset(file_id, trim(adjustl(data_path))//"state_info/i_list", i_list)

            allocate(k_vec_red_list(n_states, 3), source = 0.0_dp)
            call hdf_read_dataset(file_id, trim(adjustl(data_path))//"state_info/k_vec_red_list", k_vec_red_list)

            allocate(nlm_list(n_states, 3), source = 0)
            call hdf_read_dataset(file_id, trim(adjustl(data_path))//"state_info/nlm_list", nlm_list)

            allocate(nj_list(n_states), source = 0)
            call hdf_read_dataset(file_id, trim(adjustl(data_path))//"state_info/nj_list", nj_list)

            allocate(coeff_list(n_states, maxval(nj_list), 4), source = 0.0_dp)
            call hdf_read_dataset(file_id, trim(adjustl(data_path))//"state_info/coeff_list", coeff_list)

            allocate(eq_pos_red_list(n_states, 3), source = 0.0_dp)
            call hdf_read_dataset(file_id, trim(adjustl(data_path))//"state_info/eq_pos_red_list", eq_pos_red_list)

            call hdf_read_dataset(file_id, trim(adjustl(data_path))//"state_info/jac_list", self%jac_list)
            call hdf_read_dataset(file_id, trim(adjustl(data_path))//"state_info/Zeff_list", self%Zeff_list)
            call hdf_read_dataset(file_id, trim(adjustl(data_path))//"state_info/k_id_list", self%k_id_list)

            ! config
            call hdf_read_dataset(file_id, trim(adjustl(data_path))//"config/n_x_grid", n_x_grid)
            call hdf_read_dataset(file_id, trim(adjustl(data_path))//"config/n_r_vec_grid", n_r_vec_grid)

            call create_FFT_plan_pair(n_x_grid, FFT_plans)

            ! Assumption for now...
            spin_dof = 1

            ! initialize all states
            do n = 1, n_states

                ! elec state
                self%states(n)%energy = energy_list(n)
                self%states(n)%i = i_list(n)
                self%states(n)%spin_dof = spin_dof

                ! bloch state
                self%states(n)%k_vec_red = k_vec_red_list(n, :)
                self%states(n)%k_vec = matmul( exdm_inputs%material%k_red_to_xyz, k_vec_red_list(n, :) )
                self%states(n)%n_x_grid = n_x_grid
                self%states(n)%FFT_plans = FFT_plans

                ! PW basis
                self%states(n)%n = nlm_list(n, 1)
                self%states(n)%l = nlm_list(n, 2)
                self%states(n)%m = nlm_list(n, 3)

                allocate(self%states(n)%nl_list(nj_list(n)), source = 0.0_dp)
                self%states(n)%nl_list = coeff_list(n, :nj_list(n), 1)

                allocate(self%states(n)%Zl_list(nj_list(n)), source = 0.0_dp)
                self%states(n)%Zl_list = coeff_list(n, :nj_list(n), 2)

                allocate(self%states(n)%norml_list(nj_list(n)), source = 0.0_dp)
                self%states(n)%norml_list = coeff_list(n, :nj_list(n), 3)

                allocate(self%states(n)%Cnl_list(nj_list(n)), source = 0.0_dp)
                self%states(n)%Cnl_list = coeff_list(n, :nj_list(n), 4)

                self%states(n)%eq_pos_red = eq_pos_red_list(n, :)
                self%states(n)%eq_pos = matmul( exdm_inputs%material%red_to_xyz, eq_pos_red_list(n, :) )

                self%states(n)%n_r_vec_grid = n_r_vec_grid

                self%states(n)%red_to_xyz = exdm_inputs%material%red_to_xyz
                self%states(n)%pc_vol = exdm_inputs%material%pc_vol

                allocate(self%states(n)%r_vec_grid(&
                    product(n_r_vec_grid), 3), &
                    source = 0.0_dp)

                r = 0
                do r1 = 1, n_r_vec_grid(1)
                    do r2 = 1, n_r_vec_grid(2)
                        do r3 = 1, n_r_vec_grid(3)

                            r = r + 1

                            self%states(n)%r_vec_grid(r, :) = matmul( self%states(n)%red_to_xyz, &
                                                        [r1 - 1, r2 - 1, r3 - 1] - n_r_vec_grid/2 )

                        end do
                    end do
                end do

            end do

        else

            allocate(self%states(0))
            allocate(self%jac_list(0))
            allocate(self%Zeff_list(0))
            allocate(self%k_id_list(0))

        end if

        call hdf_close_file(file_id)

    end subroutine

end module
