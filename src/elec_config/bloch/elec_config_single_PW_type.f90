module elec_config_bloch_single_PW_type

    use prec_util, only: dp

    use elec_state_bloch_single_PW_type

    implicit none

    type :: elec_config_bloch_single_PW_t

        type(elec_state_bloch_single_PW_t), allocatable :: states(:)

        real(dp), allocatable :: jac_list(:)
        real(dp), allocatable :: Zeff_list(:)
        integer, allocatable :: k_id_list(:)

        contains

            procedure :: initialize => elec_config_bloch_single_PW_type_initialize

    end type

contains

    subroutine elec_config_bloch_single_PW_type_initialize(self, exdm_inputs, type)

        use hdf5_utils
        use exdm_inputs_type

        use FFT_util

        implicit none

        class(elec_config_bloch_single_PW_t) :: self
        type(exdm_inputs_t) :: exdm_inputs
        character(len=*) :: type

        integer(HID_t) :: file_id

        integer :: n_x_grid(3)
        integer :: FFT_plans(2, 8)
        integer :: spin_dof

        integer :: n, n_states

        integer :: dims(6)

        character(len=512) :: data_path

        real(dp), allocatable :: p_vec_list(:, :)
        real(dp), allocatable :: energy_list(:)
        integer, allocatable :: i_list(:)

        call hdf_open_file(file_id, exdm_inputs%elec_config_input%filename, &
            status='OLD', action='READ')

        ! make sure the dataset exists
        data_path = "elec_states/"//trim(adjustl(type))//"/bloch/single_PW/"
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

            allocate(p_vec_list(n_states, 3), source = 0.0_dp)
            call hdf_read_dataset(file_id, trim(adjustl(data_path))//"state_info/p_vec_list", p_vec_list)

            call hdf_read_dataset(file_id, trim(adjustl(data_path))//"state_info/jac_list", self%jac_list)
            call hdf_read_dataset(file_id, trim(adjustl(data_path))//"state_info/Zeff_list", self%Zeff_list)
            call hdf_read_dataset(file_id, trim(adjustl(data_path))//"state_info/k_id_list", self%k_id_list)

            ! config
            call hdf_read_dataset(file_id, trim(adjustl(data_path))//"config/n_x_grid", n_x_grid)

            call create_FFT_plan_pair(n_x_grid, FFT_plans)

            ! initialize all states
            do n = 1, n_states

                ! elec state
                self%states(n)%energy = energy_list(n)
                self%states(n)%i = i_list(n)
                self%states(n)%spin_dof = 1

                ! single PW
                self%states(n)%p_vec = p_vec_list(n, :)

                call compute_kG_from_p(p_vec_list(n, :), self%states(n)%k_vec, self%states(n)%G_vec, &
                    exdm_inputs%material%k_red_to_xyz, exdm_inputs%material%k_xyz_to_red)

                self%states(n)%G_vec_red = matmul( exdm_inputs%material%k_xyz_to_red, self%states(n)%G_vec)

                ! bloch state
                self%states(n)%k_vec_red = matmul( exdm_inputs%material%k_xyz_to_red, self%states(n)%k_vec)
                self%states(n)%n_x_grid = n_x_grid
                self%states(n)%FFT_plans = FFT_plans

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
