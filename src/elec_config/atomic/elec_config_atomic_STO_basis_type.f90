module elec_config_atomic_STO_basis_type

    use prec_util, only: dp

    use hdf5_utils 

    use elec_state_atomic_STO_basis_type

    implicit none

    type :: elec_config_atomic_STO_basis_t

        type(elec_state_atomic_STO_basis_t), allocatable :: states(:)

        real(dp), allocatable :: jac_list(:)
        real(dp), allocatable :: Zeff_list(:)

        contains

            procedure :: initialize => elec_config_atomic_STO_basis_type_initialize

    end type

contains

    subroutine elec_config_atomic_STO_basis_type_initialize(self, exdm_inputs, type, &
            hdf5_file_id)

        use exdm_inputs_type

        implicit none

        class(elec_config_atomic_STO_basis_t) :: self
        type(exdm_inputs_t) :: exdm_inputs
        character(len=*) :: type

        integer(HID_t) :: hdf5_file_id

        character(len=512) :: data_path
        integer :: dims(6)
        integer :: rank
        integer :: n_states

        real(dp), allocatable :: energy_list(:)
        integer, allocatable :: i_list(:)
        integer, allocatable :: nlm_list(:, :)
        integer, allocatable :: nj_list(:)
        real(dp), allocatable :: coeff_list(:, :, :) 

        integer :: spin_dof

        integer :: n

        ! make sure the dataset exists
        data_path = "elec_states/"//trim(adjustl(type))//"/atomic/STO_basis/"
        if ( hdf_exists(hdf5_file_id, trim(adjustl(data_path))) ) then

            ! load information from HDF5 file

            ! find the number of states in this approximation
            call hdf_get_dims(hdf5_file_id, trim(adjustl(data_path))//"state_info/energy_list", dims)
            n_states = dims(1)

            allocate(self%states(n_states))
            allocate(self%jac_list(n_states), source = 1.0_dp)
            allocate(self%Zeff_list(n_states), source = 1.0_dp)

            allocate(energy_list(n_states), source = 0.0_dp)
            call hdf_read_dataset(hdf5_file_id, trim(adjustl(data_path))//"state_info/energy_list", energy_list)

            allocate(i_list(n_states), source = 1)
            call hdf_read_dataset(hdf5_file_id, trim(adjustl(data_path))//"state_info/i_list", i_list)

            allocate(nlm_list(n_states, 3), source = 0)
            call hdf_read_dataset(hdf5_file_id, trim(adjustl(data_path))//"state_info/nlm_list", nlm_list)

            allocate(nj_list(n_states), source = 0)
            call hdf_read_dataset(hdf5_file_id, trim(adjustl(data_path))//"state_info/nj_list", nj_list)

            allocate(coeff_list(n_states, maxval(nj_list), 4), source = 0.0_dp)
            call hdf_read_dataset(hdf5_file_id, trim(adjustl(data_path))//"state_info/coeff_list", coeff_list)

            call hdf_read_dataset(hdf5_file_id, trim(adjustl(data_path))//"state_info/jac_list", self%jac_list)
            call hdf_read_dataset(hdf5_file_id, trim(adjustl(data_path))//"state_info/Zeff_list", self%Zeff_list)

            ! assumption for now...
            spin_dof = 1

            ! set variables within each state

            do n = 1, n_states

                ! elec state
                self%states(n)%energy = energy_list(n)
                self%states(n)%i = i_list(n)
                self%states(n)%spin_dof = spin_dof

                ! atomic state
                self%states(n)%l = nlm_list(n, 2)
                self%states(n)%m = nlm_list(n, 3)

                ! generate grid points

                ! allocate(self%states(n)%sph_x_list(n_r*n_theta*n_phi, 3), source = 0.0_dp)

                ! STO basis
                self%states(n)%n = nlm_list(n, 1) 

                allocate(self%states(n)%nl_list(nj_list(n)), source = 0.0_dp)
                self%states(n)%nl_list = coeff_list(n, :nj_list(n), 1)

                allocate(self%states(n)%Zl_list(nj_list(n)), source = 0.0_dp)
                self%states(n)%Zl_list = coeff_list(n, :nj_list(n), 2)

                allocate(self%states(n)%norml_list(nj_list(n)), source = 0.0_dp)
                self%states(n)%norml_list = coeff_list(n, :nj_list(n), 3)

                allocate(self%states(n)%Cnl_list(nj_list(n)), source = 0.0_dp)
                self%states(n)%Cnl_list = coeff_list(n, :nj_list(n), 4)

            end do

        else

            allocate(self%states(0))
            allocate(self%jac_list(0))
            allocate(self%Zeff_list(0))

        end if

    end subroutine

end module
