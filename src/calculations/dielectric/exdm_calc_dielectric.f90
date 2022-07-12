module exdm_calc_dielectric

    use prec_util, only: dp
    use timer_util

    use exdm_inputs_type
    use exdm_elec_config_type

    use elec_state_type

    implicit none

contains

    subroutine exdm_dielectric(n_proc, proc_id, root_proc_id, exdm_inputs, exdm_elec_config)

        use mpi

        use hdf5_utils

        use scatter_physics_functions

        use dielectric_calculator

        implicit none

        type(timer_t) :: timer_compute

        integer :: proc_id, root_proc_id
        integer :: n_proc

        type(exdm_inputs_t) :: exdm_inputs
        type(exdm_elec_config_t) :: exdm_elec_config

        integer(HID_t) :: file_id

        integer :: mpi_err

        complex(dp), allocatable :: total_dielectric(:, :, :, :, :) ! Dim : [ n_q_mag, n_q_theta, n_q_phi, n_E, n_widths]
        complex(dp), allocatable :: dielectric(:, :, :, :, :) ! Dim : [ n_q_mag, n_q_theta, n_q_phi, n_E, n_widths]

        integer :: q, q_t, q_p

        real(dp), allocatable :: extra_FF(:, :)

        ! initialize 
        if ( exdm_inputs%control%verbose ) then
            print*, 'Initializing dielectric calculation...'
            print*
        end if

        if ( proc_id == root_proc_id ) then
            call timer_compute%start()
        end if

        if ( proc_id == root_proc_id ) then

            allocate(total_dielectric(&
                        exdm_inputs%numerics_dielectric%n_q_bins, &
                        exdm_inputs%numerics_dielectric%n_q_theta, &
                        exdm_inputs%numerics_dielectric%n_q_phi, &
                        exdm_inputs%numerics_dielectric%n_E_bins, &
                        size(exdm_inputs%numerics_dielectric%widths, 1) &
                        ), &
                source = ( 0.0_dp, 0.0_dp ))

        end if

        allocate(dielectric(&
                        exdm_inputs%numerics_dielectric%n_q_bins, &
                        exdm_inputs%numerics_dielectric%n_q_theta, &
                        exdm_inputs%numerics_dielectric%n_q_phi, &
                        exdm_inputs%numerics_dielectric%n_E_bins, &
                        size(exdm_inputs%numerics_dielectric%widths, 1) &
                    ), &
            source = ( 0.0_dp, 0.0_dp ))

        ! check memory
        if ( 16.0_dp*product(shape(total_dielectric)) > 1.0e9_dp ) then

            if ( exdm_inputs%control%verbose ) then

                print*, 'WARNING: total_dielectric size greater than 1 GB.'
                print*

            end if

        end if

        if ( exdm_inputs%control%verbose ) then
            print*, 'Done initializing dielectric calculation!'
            print*
        end if

        if ( exdm_inputs%control%verbose ) then
            print*, 'Computing dielectric...'
            print*
        end if

        ! bloch PW -> bloch PW
        allocate(extra_FF(size(exdm_elec_config%init_bloch_PW_basis_config%states), &
                          size(exdm_elec_config%fin_bloch_PW_basis_config%states)), &
                source = 1.0_dp)

        call dielectric_compute(n_proc, proc_id, root_proc_id, exdm_inputs, &
                                    exdm_elec_config%init_bloch_PW_basis_config%states, &
                                    exdm_elec_config%fin_bloch_PW_basis_config%states, &
                                    exdm_elec_config%init_bloch_PW_basis_config%jac_list, &
                                    exdm_elec_config%fin_bloch_PW_basis_config%jac_list, &
                                    extra_FF, &
                                    dielectric)

        deallocate(extra_FF)

        ! bloch PW -> bloch STO
        allocate(extra_FF(size(exdm_elec_config%init_bloch_PW_basis_config%states), &
                          size(exdm_elec_config%fin_bloch_STO_basis_config%states)), &
                source = 1.0_dp)

        call dielectric_compute(n_proc, proc_id, root_proc_id, &
                                            exdm_inputs, &
                                            exdm_elec_config%init_bloch_PW_basis_config%states, &
                                            exdm_elec_config%fin_bloch_STO_basis_config%states, &
                                            exdm_elec_config%init_bloch_PW_basis_config%jac_list, &
                                            exdm_elec_config%fin_bloch_STO_basis_config%jac_list, &
                                            extra_FF, &
                                            dielectric)

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

        call dielectric_compute(n_proc, proc_id, root_proc_id, &
                                            exdm_inputs, &
                                            exdm_elec_config%init_bloch_PW_basis_config%states, &
                                            exdm_elec_config%fin_bloch_single_PW_config%states, &
                                            exdm_elec_config%init_bloch_PW_basis_config%jac_list, &
                                            exdm_elec_config%fin_bloch_single_PW_config%jac_list, &
                                            extra_FF, &
                                            dielectric)
        deallocate(extra_FF)

        ! bloch STO -> bloch PW
        allocate(extra_FF(size(exdm_elec_config%init_bloch_STO_basis_config%states), &
                          size(exdm_elec_config%fin_bloch_PW_basis_config%states)), &
                source = 1.0_dp)

        call dielectric_compute(n_proc, proc_id, root_proc_id, &
                                            exdm_inputs, &
                                            exdm_elec_config%init_bloch_STO_basis_config%states, &
                                            exdm_elec_config%fin_bloch_PW_basis_config%states, &
                                            exdm_elec_config%init_bloch_STO_basis_config%jac_list, &
                                            exdm_elec_config%fin_bloch_PW_basis_config%jac_list, &
                                            extra_FF, &
                                            dielectric)
        deallocate(extra_FF)

        ! bloch STO -> bloch STO
        allocate(extra_FF(size(exdm_elec_config%init_bloch_STO_basis_config%states), &
                          size(exdm_elec_config%fin_bloch_STO_basis_config%states)), &
                source = 1.0_dp)

        call dielectric_compute(n_proc, proc_id, root_proc_id, &
                                            exdm_inputs, &
                                            exdm_elec_config%init_bloch_STO_basis_config%states, &
                                            exdm_elec_config%fin_bloch_STO_basis_config%states, &
                                            exdm_elec_config%init_bloch_STO_basis_config%jac_list, &
                                            exdm_elec_config%fin_bloch_STO_basis_config%jac_list, &
                                            extra_FF, &
                                            dielectric)
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

        call dielectric_compute(n_proc, proc_id, root_proc_id, &
                                            exdm_inputs, &
                                            exdm_elec_config%init_bloch_STO_basis_config%states, &
                                            exdm_elec_config%fin_bloch_single_PW_config%states, &
                                            exdm_elec_config%init_bloch_STO_basis_config%jac_list, &
                                            exdm_elec_config%fin_bloch_single_PW_config%jac_list, &
                                            extra_FF, &
                                            dielectric)
        deallocate(extra_FF)

        ! single PW -> bloch PW
        allocate(extra_FF(size(exdm_elec_config%init_bloch_single_PW_config%states), &
                          size(exdm_elec_config%fin_bloch_PW_basis_config%states)), &
                source = 1.0_dp)

        call dielectric_compute(n_proc, proc_id, root_proc_id, &
                                            exdm_inputs, &
                                            exdm_elec_config%init_bloch_single_PW_config%states, &
                                            exdm_elec_config%fin_bloch_PW_basis_config%states, &
                                            exdm_elec_config%init_bloch_single_PW_config%jac_list, &
                                            exdm_elec_config%fin_bloch_PW_basis_config%jac_list, &
                                            extra_FF, &
                                            dielectric)
        deallocate(extra_FF)

        ! single PW -> bloch STO
        allocate(extra_FF(size(exdm_elec_config%init_bloch_single_PW_config%states), &
                          size(exdm_elec_config%fin_bloch_STO_basis_config%states)), &
                source = 1.0_dp)

        call dielectric_compute(n_proc, proc_id, root_proc_id, &
                                            exdm_inputs, &
                                            exdm_elec_config%init_bloch_single_PW_config%states, &
                                            exdm_elec_config%fin_bloch_STO_basis_config%states, &
                                            exdm_elec_config%init_bloch_single_PW_config%jac_list, &
                                            exdm_elec_config%fin_bloch_STO_basis_config%jac_list, &
                                            extra_FF, &
                                            dielectric)
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

        call dielectric_compute(n_proc, proc_id, root_proc_id, &
                                            exdm_inputs, &
                                            exdm_elec_config%init_bloch_single_PW_config%states, &
                                            exdm_elec_config%fin_bloch_single_PW_config%states, &
                                            exdm_elec_config%init_bloch_single_PW_config%jac_list, &
                                            exdm_elec_config%fin_bloch_single_PW_config%jac_list, &
                                            extra_FF, &
                                            dielectric)
        deallocate(extra_FF)

        if ( exdm_inputs%control%verbose ) then
            print*, 'Done computing dielectric!'
            print*
        end if

        call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
        
        ! reduce output
        if ( exdm_inputs%control%verbose ) then
            print*, 'Reducing dielectric...'
            print*
        end if

        call MPI_Reduce(dielectric, &
                        total_dielectric, &
                        size(dielectric), &
                        MPI_DOUBLE_COMPLEX, &
                        MPI_SUM, &
                        root_proc_id, &
                        MPI_COMM_WORLD, &
                        mpi_err)

        if ( exdm_inputs%control%verbose ) then
            print*, 'Done reducing dielectric!'
            print*
        end if

        if ( proc_id == root_proc_id ) then

            ! save data
            if ( exdm_inputs%control%verbose ) then
                print*, 'Saving dielectric information...'
                print*
            end if

            ! factor of 1 in the dielectric formula
            total_dielectric = (1.0_dp, 0.0_dp) + total_dielectric

            call timer_compute%end()

            call hdf_open_file(file_id, exdm_inputs%control%out_filename, &
                status='OLD', action='WRITE')

            ! save timing information
            call hdf_write_dataset(file_id, 'timing/dt_compute', timer_compute%dt)

            call hdf_close_file(file_id)

            call save_dielectric(exdm_inputs%control%out_filename, total_dielectric)

            if ( exdm_inputs%control%verbose ) then
                print*, 'Done saving dielectric information!'
                print*
            end if

        end if

    end subroutine

    subroutine save_dielectric(filename, dielectric)

        use hdf5_utils

        implicit none

        character(len=*) :: filename
        complex(dp) :: dielectric(:, :, :, :, :)
            ! [n_q_mag, n_q_theta, n_q_phi, n_E, n_widths]

        integer(HID_t) :: file_id

        integer :: w

        character(len=512) :: data_path
        character(len=512) :: w_str

        call hdf_open_file(file_id, filename, status='OLD', action='WRITE')

        call hdf_create_group(file_id, 'dielectric')

        do w = 1, size(dielectric, 5)

            write(w_str, *) w

            data_path = 'dielectric/width_'//trim(adjustl(w_str))

            call hdf_create_group(file_id, trim(adjustl(data_path)))

            call hdf_write_dataset(file_id, &
                trim(adjustl(data_path))//'/dielectric_r', &
                real(dielectric(:, :, :, :, w)))

            call hdf_write_dataset(file_id, &
                trim(adjustl(data_path))//'/dielectric_c', &
                aimag(dielectric(:, :, :, :, w)))

        end do

        call hdf_close_file(file_id)

    end subroutine

end module
