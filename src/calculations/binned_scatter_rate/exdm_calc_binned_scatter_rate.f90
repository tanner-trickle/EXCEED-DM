module exdm_calc_binned_scatter_rate

    use prec_util, only: dp
    use timer_util

    use exdm_inputs_type
    use exdm_elec_config_type

    use elec_state_type

    implicit none

contains

    subroutine exdm_binned_scatter_rate(n_proc, proc_id, root_proc_id, exdm_inputs, exdm_elec_config)

        use mpi

        use hdf5_utils

        use scatter_physics_functions

        use binned_scatter_rate_calculator

        implicit none

        integer :: proc_id, root_proc_id
        integer :: n_proc

        type(exdm_inputs_t) :: exdm_inputs
        type(exdm_elec_config_t) :: exdm_elec_config

        integer :: mpi_err

        type(timer_t) :: timer_compute
        integer(HID_T) :: file_id

        real(dp), allocatable :: total_binned_scatter_rate(:, :, :, :, :, :) ! [n_q, n_E, n_mX, n_models, n_vE, n_init_groups]
        real(dp), allocatable :: binned_scatter_rate(:, :, :, :, :, :)

        class(elec_state_t), allocatable :: init_states(:), fin_states(:)

        real(dp), allocatable :: extra_FF(:, :)
        integer :: n_init

        ! initialize 
        if ( exdm_inputs%control%verbose ) then
            print*, 'Initializing binned scatter rate calculation...'
            print*
        end if

        if ( proc_id == root_proc_id ) then
            call timer_compute%start()
        end if

        if ( proc_id == root_proc_id ) then

            allocate(total_binned_scatter_rate(&
                        exdm_inputs%numerics_binned_scatter_rate%n_q_bins, &
                        exdm_inputs%numerics_binned_scatter_rate%n_E_bins, &
                        size(exdm_inputs%dm_model%mX), &
                        size(exdm_inputs%dm_model%med_FF), &
                        size(exdm_inputs%astroph_model%v_e_list, 1), &
                        exdm_elec_config%n_init_groups), &
                source = 0.0_dp)

        end if

        allocate(binned_scatter_rate(&
                    exdm_inputs%numerics_binned_scatter_rate%n_q_bins, &
                    exdm_inputs%numerics_binned_scatter_rate%n_E_bins, &
                    size(exdm_inputs%dm_model%mX), &
                    size(exdm_inputs%dm_model%med_FF), &
                    size(exdm_inputs%astroph_model%v_e_list, 1), &
                    exdm_elec_config%n_init_groups), &
            source = 0.0_dp)

        ! check memory
        if ( 16.0_dp*product(shape(total_binned_scatter_rate)) > 1.0e9_dp ) then

            if ( exdm_inputs%control%verbose ) then

                print*, 'WARNING: total_binned_scatter_rate size greater than 1 GB.'
                print*

            end if

        end if

        if ( exdm_inputs%control%verbose ) then
            print*, 'Done initializing binned scatter rate calculation!'
            print*
        end if


        if ( exdm_inputs%control%verbose ) then
            print*, 'Computing binned scatter rate...'
            print*
        end if

        ! simplify this?

        ! bloch PW -> bloch PW
        allocate(extra_FF(size(exdm_elec_config%init_bloch_PW_basis_config%states), &
                          size(exdm_elec_config%fin_bloch_PW_basis_config%states)), &
                source = 1.0_dp)

        call binned_scatter_rate_compute(n_proc, proc_id, root_proc_id, &
                                            exdm_inputs, &
                                            exdm_elec_config%init_bloch_PW_basis_config%states, &
                                            exdm_elec_config%fin_bloch_PW_basis_config%states, &
                                            exdm_elec_config%init_bloch_PW_basis_config%jac_list, &
                                            exdm_elec_config%fin_bloch_PW_basis_config%jac_list, &
                                            extra_FF, &
                                            binned_scatter_rate)
        deallocate(extra_FF)

        ! bloch PW -> bloch STO
        allocate(extra_FF(size(exdm_elec_config%init_bloch_PW_basis_config%states), &
                          size(exdm_elec_config%fin_bloch_STO_basis_config%states)), &
                source = 1.0_dp)

        call binned_scatter_rate_compute(n_proc, proc_id, root_proc_id, &
                                            exdm_inputs, &
                                            exdm_elec_config%init_bloch_PW_basis_config%states, &
                                            exdm_elec_config%fin_bloch_STO_basis_config%states, &
                                            exdm_elec_config%init_bloch_PW_basis_config%jac_list, &
                                            exdm_elec_config%fin_bloch_STO_basis_config%jac_list, &
                                            extra_FF, &
                                            binned_scatter_rate)
        deallocate(extra_FF)

        ! bloch PW -> single PW
        allocate(extra_FF(size(exdm_elec_config%init_bloch_PW_basis_config%states), &
                          size(exdm_elec_config%fin_bloch_single_PW_config%states)), &
                source = 1.0_dp)

        if ( product(shape(extra_FF)) > 0 ) then

            call fermi_factor_extra_FF(extra_FF, &
                exdm_elec_config%init_bloch_PW_basis_config%Zeff_list, &
                exdm_elec_config%fin_bloch_single_PW_config%states(:)%energy)

        end if

        call binned_scatter_rate_compute(n_proc, proc_id, root_proc_id, &
                                            exdm_inputs, &
                                            exdm_elec_config%init_bloch_PW_basis_config%states, &
                                            exdm_elec_config%fin_bloch_single_PW_config%states, &
                                            exdm_elec_config%init_bloch_PW_basis_config%jac_list, &
                                            exdm_elec_config%fin_bloch_single_PW_config%jac_list, &
                                            extra_FF, &
                                            binned_scatter_rate)
        deallocate(extra_FF)

        ! bloch STO -> bloch PW
        allocate(extra_FF(size(exdm_elec_config%init_bloch_STO_basis_config%states), &
                          size(exdm_elec_config%fin_bloch_PW_basis_config%states)), &
                source = 1.0_dp)

        call binned_scatter_rate_compute(n_proc, proc_id, root_proc_id, &
                                            exdm_inputs, &
                                            exdm_elec_config%init_bloch_STO_basis_config%states, &
                                            exdm_elec_config%fin_bloch_PW_basis_config%states, &
                                            exdm_elec_config%init_bloch_STO_basis_config%jac_list, &
                                            exdm_elec_config%fin_bloch_PW_basis_config%jac_list, &
                                            extra_FF, &
                                            binned_scatter_rate)
        deallocate(extra_FF)

        ! bloch STO -> bloch STO
        allocate(extra_FF(size(exdm_elec_config%init_bloch_STO_basis_config%states), &
                          size(exdm_elec_config%fin_bloch_STO_basis_config%states)), &
                source = 1.0_dp)

        call binned_scatter_rate_compute(n_proc, proc_id, root_proc_id, &
                                            exdm_inputs, &
                                            exdm_elec_config%init_bloch_STO_basis_config%states, &
                                            exdm_elec_config%fin_bloch_STO_basis_config%states, &
                                            exdm_elec_config%init_bloch_STO_basis_config%jac_list, &
                                            exdm_elec_config%fin_bloch_STO_basis_config%jac_list, &
                                            extra_FF, &
                                            binned_scatter_rate)
        deallocate(extra_FF)

        ! bloch STO -> single PW
        allocate(extra_FF(size(exdm_elec_config%init_bloch_STO_basis_config%states), &
                          size(exdm_elec_config%fin_bloch_single_PW_config%states)), &
                source = 1.0_dp)

        if ( product(shape(extra_FF)) > 0 ) then

            call fermi_factor_extra_FF(extra_FF, &
                exdm_elec_config%init_bloch_STO_basis_config%Zeff_list, &
                exdm_elec_config%fin_bloch_single_PW_config%states(:)%energy)

        end if

        call binned_scatter_rate_compute(n_proc, proc_id, root_proc_id, &
                                            exdm_inputs, &
                                            exdm_elec_config%init_bloch_STO_basis_config%states, &
                                            exdm_elec_config%fin_bloch_single_PW_config%states, &
                                            exdm_elec_config%init_bloch_STO_basis_config%jac_list, &
                                            exdm_elec_config%fin_bloch_single_PW_config%jac_list, &
                                            extra_FF, &
                                            binned_scatter_rate)
        deallocate(extra_FF)

        ! single PW -> bloch PW
        allocate(extra_FF(size(exdm_elec_config%init_bloch_single_PW_config%states), &
                          size(exdm_elec_config%fin_bloch_PW_basis_config%states)), &
                source = 1.0_dp)

        call binned_scatter_rate_compute(n_proc, proc_id, root_proc_id, &
                                            exdm_inputs, &
                                            exdm_elec_config%init_bloch_single_PW_config%states, &
                                            exdm_elec_config%fin_bloch_PW_basis_config%states, &
                                            exdm_elec_config%init_bloch_single_PW_config%jac_list, &
                                            exdm_elec_config%fin_bloch_PW_basis_config%jac_list, &
                                            extra_FF, &
                                            binned_scatter_rate)
        deallocate(extra_FF)

        ! single PW -> bloch STO
        allocate(extra_FF(size(exdm_elec_config%init_bloch_single_PW_config%states), &
                          size(exdm_elec_config%fin_bloch_STO_basis_config%states)), &
                source = 1.0_dp)

        call binned_scatter_rate_compute(n_proc, proc_id, root_proc_id, &
                                            exdm_inputs, &
                                            exdm_elec_config%init_bloch_single_PW_config%states, &
                                            exdm_elec_config%fin_bloch_STO_basis_config%states, &
                                            exdm_elec_config%init_bloch_single_PW_config%jac_list, &
                                            exdm_elec_config%fin_bloch_STO_basis_config%jac_list, &
                                            extra_FF, &
                                            binned_scatter_rate)
        deallocate(extra_FF)

        ! single PW -> single PW
        allocate(extra_FF(size(exdm_elec_config%init_bloch_single_PW_config%states), &
                          size(exdm_elec_config%fin_bloch_single_PW_config%states)), &
                source = 1.0_dp)

        if ( product(shape(extra_FF)) > 0 ) then
            call fermi_factor_extra_FF(extra_FF, &
                exdm_elec_config%init_bloch_single_PW_config%Zeff_list, &
                exdm_elec_config%fin_bloch_single_PW_config%states(:)%energy)
        end if

        call binned_scatter_rate_compute(n_proc, proc_id, root_proc_id, &
                                            exdm_inputs, &
                                            exdm_elec_config%init_bloch_single_PW_config%states, &
                                            exdm_elec_config%fin_bloch_single_PW_config%states, &
                                            exdm_elec_config%init_bloch_single_PW_config%jac_list, &
                                            exdm_elec_config%fin_bloch_single_PW_config%jac_list, &
                                            extra_FF, &
                                            binned_scatter_rate)
        deallocate(extra_FF)

        if ( exdm_inputs%control%verbose ) then
            print*, 'Done computing binned scatter rate!'
            print*
        end if

        call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
        
        ! reduce output
        if ( exdm_inputs%control%verbose ) then
            print*, 'Reducing binned scatter rate...'
            print*
        end if

        call MPI_Reduce(binned_scatter_rate, &
                        total_binned_scatter_rate, &
                        size(binned_scatter_rate), &
                        MPI_DOUBLE, &
                        MPI_SUM, &
                        root_proc_id, &
                        MPI_COMM_WORLD, &
                        mpi_err)

        if ( exdm_inputs%control%verbose ) then
            print*, 'Done reducing binned scatter rate!'
            print*
        end if

        if ( proc_id == root_proc_id ) then

            ! save data
            if ( exdm_inputs%control%verbose ) then
                print*, 'Saving binned scatter rate information...'
                print*
            end if

            call timer_compute%end()

            call hdf_open_file(file_id, exdm_inputs%control%out_filename, &
                status='OLD', action='WRITE')

            ! save timing information
            call hdf_write_dataset(file_id, 'timing/dt_compute', timer_compute%dt)

            call hdf_close_file(file_id)

            call save_binned_scatter_rate(exdm_inputs%control%out_filename, total_binned_scatter_rate)

            if ( exdm_inputs%control%verbose ) then
                print*, 'Done saving binned scatter rate information!'
                print*
            end if

        end if

    end subroutine

    subroutine save_binned_scatter_rate(filename, binned_scatter_rate)

        use hdf5_utils

        implicit none

        character(len=*) :: filename
        real(dp) :: binned_scatter_rate(:, :, :, :, :, :)
            ! [n_q, n_E, n_mX, n_models, n_vE, n_init_groups]

        integer(HID_t) :: file_id

        integer :: m, n, v, i

        character(len=512) :: data_path
        character(len=512) :: base_data_path
        character(len=512) :: n_str, m_str, v_str, i_str

        call hdf_open_file(file_id, filename, status='OLD', action='WRITE')

        call hdf_create_group(file_id, 'binned_scatter_rate')

        base_data_path = 'binned_scatter_rate'

        do n = 1, size(binned_scatter_rate, 4)
            do v = 1, size(binned_scatter_rate, 5)
                do m = 1, size(binned_scatter_rate, 3)

                    data_path = base_data_path

                    if ( size(binned_scatter_rate, 4) > 1 ) then

                        write(n_str, *) n
                        data_path = trim(adjustl(data_path))//'/model_'//trim(adjustl(n_str))

                        if ( .not. hdf_exists(file_id, trim(adjustl(data_path))) ) then
                            call hdf_create_group(file_id, trim(adjustl(data_path)))
                        end if

                    end if

                    if ( size(binned_scatter_rate, 5) > 1 ) then

                        write(v_str, *) v
                        data_path = trim(adjustl(data_path))//'/v_e_'//trim(adjustl(v_str))

                        if ( .not. hdf_exists(file_id, trim(adjustl(data_path))) ) then
                            call hdf_create_group(file_id, trim(adjustl(data_path)))
                        end if

                    end if

                    if ( size(binned_scatter_rate, 3) > 1 ) then

                        write(m_str, *) m
                        data_path = trim(adjustl(data_path))//'/mass_'//trim(adjustl(m_str))

                        if ( .not. hdf_exists(file_id, trim(adjustl(data_path))) ) then
                            call hdf_create_group(file_id, trim(adjustl(data_path)))
                        end if

                    end if

                    call hdf_write_dataset(file_id, &
                        trim(adjustl(data_path))//'/total_binned_scatter_rate', &
                        sum(binned_scatter_rate(:, :, m, n, v, :), 3))

                    do i = 1, size(binned_scatter_rate, 6)

                        if ( size(binned_scatter_rate, 6) > 1 ) then

                            write(i_str, *) i
                            ! data_path = trim(adjustl(data_path))//'/i_'//trim(adjustl(i_str))

                            if ( .not. hdf_exists(file_id, trim(adjustl(data_path))//'/i_'//trim(adjustl(i_str))) ) then
                                call hdf_create_group(file_id, trim(adjustl(data_path))//'/i_'//trim(adjustl(i_str)))
                            end if

                            call hdf_write_dataset(file_id, &
                                trim(adjustl(data_path))//'/i_'//trim(adjustl(i_str))//'/binned_scatter_rate', &
                                binned_scatter_rate(:, :, m, n, v, i))

                        end if

                    end do

                end do

            end do

        end do

        call hdf_close_file(file_id)

    end subroutine

end module
