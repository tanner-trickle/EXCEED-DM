module elec_config_bloch_PW_basis_type

    use prec_util, only: dp

    use elec_state_bloch_PW_basis_type

    implicit none

    type :: elec_config_bloch_PW_basis_t

        type(elec_state_bloch_PW_basis_t), allocatable :: states(:)

        real(dp), allocatable :: jac_list(:)
        real(dp), allocatable :: Zeff_list(:)
        integer, allocatable :: k_id_list(:)

        integer, allocatable :: G_list_red(:, :)
        real(dp), allocatable :: G_list(:, :)

        contains

            procedure :: initialize => elec_config_bloch_PW_basis_type_initialize

    end type

contains

    subroutine elec_config_bloch_PW_basis_type_initialize(self, exdm_inputs, type)

        use hdf5_utils
        use exdm_inputs_type

        use FFT_util

        implicit none

        class(elec_config_bloch_PW_basis_t) :: self
        type(exdm_inputs_t) :: exdm_inputs
        character(len=*) :: type

        integer(HID_t) :: file_id

        integer :: n_states
        integer :: dims(6)

        integer :: rank

        integer :: i, n

        character(len=512) :: data_path

        real(dp), allocatable :: energy_list(:)
        integer, allocatable :: i_list(:)
        real(dp), allocatable :: k_vec_red_list(:, :)

        integer :: n_x_grid(3)
        integer :: FFT_plans(2, 8)
        integer :: spin_dof

        character(len=512) :: n_str

        call hdf_open_file(file_id, exdm_inputs%elec_config_input%filename, &
            status='OLD', action='READ')

        ! make sure the dataset exists
        data_path = "elec_states/"//trim(adjustl(type))//"/bloch/PW_basis/"
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

            call hdf_read_dataset(file_id, trim(adjustl(data_path))//"state_info/jac_list", self%jac_list)
            call hdf_read_dataset(file_id, trim(adjustl(data_path))//"state_info/Zeff_list", self%Zeff_list)
            call hdf_read_dataset(file_id, trim(adjustl(data_path))//"state_info/k_id_list", self%k_id_list)

            ! config
            call hdf_get_dims(file_id, trim(adjustl(data_path))//"config/G_list_red", dims)
            allocate(self%G_list_red(dims(1), dims(2)), source = 0)
            allocate(self%G_list(dims(1), dims(2)), source = 0.0_dp)
            call hdf_read_dataset(file_id, trim(adjustl(data_path))//"config/G_list_red", self%G_list_red)

            self%G_list = transpose( matmul( exdm_inputs%material%k_red_to_xyz, transpose(self%G_list_red) ) )

            do i = 1, 3
                n_x_grid(i) = maxval(self%G_list_red(:, i)) - minval(self%G_list_red(:, i)) + 1
            end do

            call create_FFT_plan_pair(n_x_grid, FFT_plans)

            call hdf_read_dataset(file_id, trim(adjustl(data_path))//"config/G_list_red", self%G_list_red)

            call hdf_get_dims(file_id, trim(adjustl(data_path))//"state_info/u_FT_r/n_1", dims)
            spin_dof = dims(2)

            ! initialize all states
            do n = 1, n_states

                write(n_str, *) n

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
                self%states(n)%hdf5_data_filename = exdm_inputs%elec_config_input%filename
                self%states(n)%hdf5_G_list_red_path = trim(adjustl(data_path))//"config/G_list_red"
                self%states(n)%hdf5_u_FT_r_path = trim(adjustl(data_path))//"state_info/u_FT_r/n_"//trim(adjustl(n_str))
                self%states(n)%hdf5_u_FT_c_path = trim(adjustl(data_path))//"state_info/u_FT_c/n_"//trim(adjustl(n_str))

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
