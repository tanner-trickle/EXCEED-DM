module exdm_calc_absorption_rate

    use prec_util, only: dp

    use timer_util

    use exdm_inputs_type 
    use exdm_elec_config_type 

    use elec_state_type

    implicit none

contains

    subroutine exdm_absorption_rate(n_proc, proc_id, root_proc_id, exdm_inputs, exdm_elec_config)

        use hdf5_utils

        use PiIF_calculator_type

        use self_energies_calculator

        use absorption_rate_calculator

        implicit none

        integer :: proc_id, root_proc_id
        integer :: n_proc

        type(exdm_inputs_t) :: exdm_inputs
        type(exdm_elec_config_t) :: exdm_elec_config

        integer :: mpi_err

        type(timer_t) :: timer_compute
        integer(HID_T) :: file_id

        real(dp), allocatable :: absorption_rate(:, :) ! [ n_mX, n_widths ]

        type(PiIF_calculator_t) :: total_PiIF_calculator
        type(PiIF_calculator_t) :: PiIF_calculator

        integer :: spin_dof

        ! initialize 
        if ( exdm_inputs%control%verbose ) then
            print*, 'Initializing absorption rate calculation...'
            print*
        end if

        if ( proc_id == root_proc_id ) then
            call timer_compute%start()
        end if

        if ( proc_id == root_proc_id ) then

            allocate(absorption_rate(size(exdm_inputs%dm_model%mX), &
                                 size(exdm_inputs%numerics_absorption_rate%widths, 1)), &
                    source = 0.0_dp)

        end if

        if ( ( size(exdm_elec_config%init_bloch_PW_basis_config%states) >= 1 ) &
            .and. ( size(exdm_elec_config%fin_bloch_PW_basis_config%states) >= 1 ) ) then

            spin_dof = min( exdm_elec_config%init_bloch_PW_basis_config%states(1)%spin_dof, &
                        exdm_elec_config%fin_bloch_PW_basis_config%states(1)%spin_dof )

        else

            spin_dof = 1

        end if

        ! which self energies to compute
        call particle_type_to_PiIF_mask(exdm_inputs%dm_model%particle_type, &
            spin_dof, &
            PiIF_calculator%mask, PiIF_calculator%n)

        ! allocate relevant Pi variables on each processor
        call PiIF_calculator%setup(size(exdm_inputs%dm_model%mX), &
            size(exdm_inputs%numerics_absorption_rate%widths, 1))

        ! allocate summed Pi's on main processor
        if ( proc_id == root_proc_id ) then

            total_PiIF_calculator%mask = PiIF_calculator%mask

            call total_PiIF_calculator%setup(size(exdm_inputs%dm_model%mX), &
                size(exdm_inputs%numerics_absorption_rate%widths, 1))

        end if

        if ( exdm_inputs%control%verbose ) then
            print*, 'Done initializing absorption rate calculation!'
            print*
        end if

        if ( exdm_inputs%control%verbose ) then
            print*, 'Computing self energies...'
            print*
        end if

        ! bloch PW -> bloch PW
        call self_energies_compute(n_proc, proc_id, root_proc_id, &
            exdm_inputs, &
            exdm_elec_config%init_bloch_PW_basis_config%states, &
            exdm_elec_config%fin_bloch_PW_basis_config%states, &
            exdm_elec_config%init_bloch_PW_basis_config%k_id_list, &
            exdm_elec_config%fin_bloch_PW_basis_config%k_id_list, &
            exdm_elec_config%init_bloch_PW_basis_config%jac_list, &
            PiIF_calculator)

        ! bloch PW -> bloch STO
        call self_energies_compute(n_proc, proc_id, root_proc_id, &
            exdm_inputs, &
            exdm_elec_config%init_bloch_PW_basis_config%states, &
            exdm_elec_config%fin_bloch_STO_basis_config%states, &
            exdm_elec_config%init_bloch_PW_basis_config%k_id_list, &
            exdm_elec_config%fin_bloch_STO_basis_config%k_id_list, &
            exdm_elec_config%init_bloch_PW_basis_config%jac_list, &
            PiIF_calculator)

        ! bloch PW -> single PW
        call self_energies_compute(n_proc, proc_id, root_proc_id, &
            exdm_inputs, &
            exdm_elec_config%init_bloch_PW_basis_config%states, &
            exdm_elec_config%fin_bloch_single_PW_config%states, &
            exdm_elec_config%init_bloch_PW_basis_config%k_id_list, &
            exdm_elec_config%fin_bloch_single_PW_config%k_id_list, &
            exdm_elec_config%init_bloch_PW_basis_config%jac_list, &
            PiIF_calculator)

        ! bloch STO -> bloch PW
        call self_energies_compute(n_proc, proc_id, root_proc_id, &
            exdm_inputs, &
            exdm_elec_config%init_bloch_STO_basis_config%states, &
            exdm_elec_config%fin_bloch_PW_basis_config%states, &
            exdm_elec_config%init_bloch_STO_basis_config%k_id_list, &
            exdm_elec_config%fin_bloch_PW_basis_config%k_id_list, &
            exdm_elec_config%init_bloch_STO_basis_config%jac_list, &
            PiIF_calculator)

        ! bloch STO -> bloch STO
        call self_energies_compute(n_proc, proc_id, root_proc_id, &
            exdm_inputs, &
            exdm_elec_config%init_bloch_STO_basis_config%states, &
            exdm_elec_config%fin_bloch_STO_basis_config%states, &
            exdm_elec_config%init_bloch_STO_basis_config%k_id_list, &
            exdm_elec_config%fin_bloch_STO_basis_config%k_id_list, &
            exdm_elec_config%init_bloch_STO_basis_config%jac_list, &
            PiIF_calculator)

        ! bloch STO -> single PW
        call self_energies_compute(n_proc, proc_id, root_proc_id, &
            exdm_inputs, &
            exdm_elec_config%init_bloch_STO_basis_config%states, &
            exdm_elec_config%fin_bloch_single_PW_config%states, &
            exdm_elec_config%init_bloch_STO_basis_config%k_id_list, &
            exdm_elec_config%fin_bloch_single_PW_config%k_id_list, &
            exdm_elec_config%init_bloch_STO_basis_config%jac_list, &
            PiIF_calculator)

        ! single PW -> bloch PW
        call self_energies_compute(n_proc, proc_id, root_proc_id, &
            exdm_inputs, &
            exdm_elec_config%init_bloch_single_PW_config%states, &
            exdm_elec_config%fin_bloch_PW_basis_config%states, &
            exdm_elec_config%init_bloch_single_PW_config%k_id_list, &
            exdm_elec_config%fin_bloch_PW_basis_config%k_id_list, &
            exdm_elec_config%init_bloch_single_PW_config%jac_list, &
            PiIF_calculator)

        ! single PW -> bloch STO
        call self_energies_compute(n_proc, proc_id, root_proc_id, &
            exdm_inputs, &
            exdm_elec_config%init_bloch_single_PW_config%states, &
            exdm_elec_config%fin_bloch_STO_basis_config%states, &
            exdm_elec_config%init_bloch_single_PW_config%k_id_list, &
            exdm_elec_config%fin_bloch_STO_basis_config%k_id_list, &
            exdm_elec_config%init_bloch_single_PW_config%jac_list, &
            PiIF_calculator)

        ! single PW -> single PW
        call self_energies_compute(n_proc, proc_id, root_proc_id, &
            exdm_inputs, &
            exdm_elec_config%init_bloch_single_PW_config%states, &
            exdm_elec_config%fin_bloch_single_PW_config%states, &
            exdm_elec_config%init_bloch_single_PW_config%k_id_list, &
            exdm_elec_config%fin_bloch_single_PW_config%k_id_list, &
            exdm_elec_config%init_bloch_single_PW_config%jac_list, &
            PiIF_calculator)

        if ( exdm_inputs%control%verbose ) then
            print*, 'Reducing self energies...'
            print*
        end if

        call reduce_PiIF(root_proc_id, PiIF_calculator, total_PiIF_calculator)

        if ( exdm_inputs%control%verbose ) then
            print*, 'Done reducing self energies!'
            print*
            print*, 'Done computing self energies!'
            print*
        end if

        if ( proc_id == root_proc_id ) then
            
            call absorption_rate_compute(exdm_inputs, spin_dof, total_PiIF_calculator, absorption_rate)

            ! save data
            if ( exdm_inputs%control%verbose ) then
                print*, 'Saving absorption rate information...'
                print*
            end if

            call timer_compute%end()

            call hdf_open_file(file_id, exdm_inputs%control%out_filename, &
                status='OLD', action='WRITE')

            ! save timing information
            call hdf_write_dataset(file_id, 'timing/dt_compute', timer_compute%dt)

            call hdf_close_file(file_id)

            call save_absorption_rate(exdm_inputs%control%out_filename, absorption_rate)

            if ( exdm_inputs%control%verbose ) then
                print*, 'Done saving absorption rate information!'
                print*
            end if

        end if

    end subroutine

    subroutine save_absorption_rate(filename, absorption_rate)

        use hdf5_utils

        implicit none

        character(len=*) :: filename

        real(dp) :: absorption_rate(:, :) ! [ n_mX, n_widths ]

        integer :: m, w

        integer(HID_t) :: file_id

        character(len=512) :: data_path
        character(len=512) :: base_data_path
        character(len=512) :: m_str, w_str

        call hdf_open_file(file_id, filename, status='OLD', action='WRITE')

        call hdf_create_group(file_id, 'absorption_rate')

        base_data_path = 'absorption_rate'

        do w = 1, size(absorption_rate, 2)
            do m = 1, size(absorption_rate, 1)

                data_path = base_data_path

                if ( size(absorption_rate, 2) > 1 ) then

                    write(w_str, *) w
                    data_path = trim(adjustl(data_path))//'/width_'//trim(adjustl(w_str))

                    if ( .not. hdf_exists(file_id, trim(adjustl(data_path))) ) then
                        call hdf_create_group(file_id, trim(adjustl(data_path)))
                    end if

                end if

                if ( size(absorption_rate, 1) > 1 ) then

                    write(m_str, *) m
                    data_path = trim(adjustl(data_path))//'/mass_'//trim(adjustl(m_str))

                    if ( .not. hdf_exists(file_id, trim(adjustl(data_path))) ) then
                        call hdf_create_group(file_id, trim(adjustl(data_path)))
                    end if

                end if

                call hdf_write_dataset(file_id, &
                    trim(adjustl(data_path))//'/absorption_rate', &
                    absorption_rate(m, w))

            end do

        end do

        call hdf_close_file(file_id)

    end subroutine

end module
