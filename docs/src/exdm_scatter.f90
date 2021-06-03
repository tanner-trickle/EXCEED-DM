module exdm_scatter
    !! Perform DM-electron scattering rate calculations

    use mpi
    use hdf5
    use h5lt

    use prec

    use control_input
    use particle_physics_scatter
    use numerics_input
    use in_med_scr
    use transition

    use transition_form_factor

    use exdm_scatter_vc
    use exdm_scatter_cc
    use exdm_scatter_vf
    use exdm_scatter_cf

    implicit none

    real(dp), allocatable :: binned_rate_t_scatter(:, :, :, :, :, :)
        !! Dim : [n_q_bins + 1, n_E_bins + 1, n_mX, n_FDM, n_time, n_tran_per_proc]
        !!
        !! Binned rate allocated on each processor, for each job to compute for
        !!
        !! Units : eV^2

    real(dp), allocatable :: binned_rate_if_scatter(:, :, :, :, :, :, :)
        !! Dim : [n_q_bins + 1, n_E_bins + 1, n_mX, n_FDM, n_time, n_init, n_fin]
        !!
        !! Binned rate between all transitions 
        !!
        !! Warning : memory intensive
        !!
        !! Units : eV^2

    real(dp), allocatable :: binned_rate_i_scatter(:, :, :, :, :, :)
        !! Dim : [n_q_bins + 1, n_E_bins + 1, n_mX, n_FDM, n_time, n_init]
        !!
        !! Binned rate, summed over final states
        !!
        !! Units : eV^2

    real(dp), allocatable :: binned_rate_scatter(:, :, :, :, :)
        !! Dim : [n_q_bins + 1, n_E_bins + 1, n_mX, n_FDM, n_time]
        !!
        !! Binned rate summed over all transitions
        !!
        !! Units : eV^2

    real(dp), allocatable :: rate_scatter(:, :, :)
        !! Dim : [n_mX, n_FDM, n_time]
        !!
        !! Total rate summed over all transitions
        !!
        !! Units : eV^2

    character(len=64) :: calc_modes_list(4) = ['vc', 'vf', 'cc', 'cf']

contains

    subroutine run_dme_scatter(proc_id, root_process, &
            out_filename, nml_filename, DFT_input_filename, sto_wf_filename, &
            core_elec_config_filename, n_proc, save_binned_rate_if, verbose)
        !! Scattering rate calculation
        implicit none

        integer :: proc_id
        integer :: root_process

        integer :: n_proc

        character(len=*) :: out_filename
        character(len=*) :: nml_filename
        character(len=*) :: DFT_input_filename
        character(len=*) :: sto_wf_filename
        character(len=*) :: core_elec_config_filename

        logical :: save_binned_rate_if

        logical, optional :: verbose

        integer :: status(MPI_STATUS_SIZE)
        integer :: tag = 0
        integer :: err

        integer :: t, i

        if ( verbose ) then
            print*, 'Starting scattering rate calculation...'
            print*
        end if

        call load_particle_physics_scatter(nml_filename, verbose = verbose)
        call load_numerics(nml_filename, verbose = verbose)
        call load_in_med_scr(nml_filename, verbose = verbose)
        call load_tff_input(nml_filename, verbose = verbose)

        call set_job_table(n_proc, n_init, n_fin, verbose=verbose)

        ! allocate arrays
        if ( proc_id == root_process ) then

            if ( save_binned_rate_if ) then

                allocate(binned_rate_if_scatter(n_q_bins + 1, n_E_bins + 1, n_mX, n_FDM, n_time, n_init, n_fin))
                binned_rate_if_scatter = 0.0_dp

                if ( (n_q_bins + 1)*(n_E_bins + 1)*n_mX*n_FDM*n_time*n_init*n_fin*8.0_dp > 1.0e10_dp ) then

                    if ( verbose ) then

                        print*, '~~~ WARNING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
                        print*
                        print*, '    Attempting to store more than 10 GB of data in binned_rate_if. If this is unwanted behavior'
                        print*, '    set save_binned_rate_if = .FALSE. in the control namelist.'
                        print*
                        print*, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
                        print*

                    end if

                end if

            end if

            allocate(binned_rate_i_scatter(n_q_bins + 1, n_E_bins + 1, n_mX, n_FDM, n_time, n_init))
            binned_rate_i_scatter = 0.0_dp

            allocate(binned_rate_scatter(n_q_bins + 1, n_E_bins + 1, n_mX, n_FDM, n_time))
            binned_rate_scatter = 0.0_dp

            allocate(rate_scatter(n_mX, n_FDM, n_time))
            rate_scatter = 0.0_dp

        end if

        allocate(binned_rate_t_scatter(n_q_bins + 1, n_E_bins + 1, n_mX, n_FDM, n_time, n_tran_per_proc))
        binned_rate_t_scatter = 0.0_dp

        ! specific calc modes to compute binned_rate_t_scatter
        ! the subroutines here have one purpose and that is to compute 
        ! the binned rate for a given transition id
        if ( trim(calc_mode) .eq. 'vc' ) then

            call run_dme_scatter_vc(binned_rate_t_scatter, n_tran_per_proc, DFT_input_filename, &
                out_filename, proc_id, root_process, verbose = verbose)

        else if ( trim(calc_mode) .eq. 'cc' ) then

            call run_dme_scatter_cc(binned_rate_t_scatter, n_tran_per_proc, DFT_input_filename, &
                sto_wf_filename, core_elec_config_filename, out_filename, proc_id, root_process, verbose = verbose)

        else if ( trim(calc_mode) .eq. 'vf' ) then

            call run_dme_scatter_vf(binned_rate_t_scatter, n_tran_per_proc, DFT_input_filename, &
                nml_filename, out_filename, proc_id, root_process, verbose = verbose)

        else if ( trim(calc_mode) .eq. 'cf' ) then

            call run_dme_scatter_cf(binned_rate_t_scatter, n_tran_per_proc, sto_wf_filename, &
                core_elec_config_filename, nml_filename, out_filename, proc_id, root_process, verbose = verbose)

        else

            if ( verbose ) then

                print*, '!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*
                print*, '    Calculation mode : ', trim(calc_mode), ' is not a valid option. Options are : '
                print*
                
                do i = 1, size(calc_modes_list)

                    print*, '        ', trim(calc_modes_list(i))
                    print* 

                end do

                print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*

            end if

            stop

        end if

        ! sends data to the main node
        if ( proc_id .ne. root_process ) then

            call MPI_SEND(binned_rate_t_scatter, &
               size(binned_rate_t_scatter), MPI_DOUBLE, root_process, tag, MPI_COMM_WORLD, err)

        end if

        if ( proc_id .eq. root_process ) then

            ! add main processors contribution
            call update_rates(binned_rate_t_scatter, proc_id)

            do i = 1, n_proc
                if ( (i - 1) .ne. root_process ) then

                    call MPI_RECV(binned_rate_t_scatter, &
                       size(binned_rate_t_scatter), MPI_DOUBLE, i - 1, MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                    ! add other processors contributions
                    call update_rates(binned_rate_t_scatter, i - 1)

                end if
            end do

        end if

        ! save data
        if ( proc_id .eq. root_process ) then
            call save_particle_physics_scatter(out_filename, verbose = verbose)
            call save_numerics(out_filename, verbose = verbose)
            call save_in_med_scr(out_filename, verbose = verbose)
            call save_rates(out_filename, verbose = verbose)
        end if

    end subroutine

    subroutine update_rates(binned_rate_t, proc_id)
        !! updates all the rate arrays with the output 
        !! of a single processors results
        implicit none

        real(dp) :: binned_rate_t(n_q_bins + 1, n_E_bins + 1, n_mX, n_FDM, n_time, n_tran_per_proc)

        integer :: i, f, m, t

        integer :: proc_id

        integer :: tran_id

        do i = 1, n_tran_per_proc

            tran_id = job_table(proc_id + 1, i)

            if ( tran_id .ne. 0 ) then

                if ( save_binned_rate_if ) then
                    binned_rate_if_scatter(:, :, :, :, :, &
                        tran_to_init_fin_id(tran_id, 1), &
                        tran_to_init_fin_id(tran_id, 2) ) = &
                    binned_rate_if_scatter(:, :, :, :, :, &
                        tran_to_init_fin_id(tran_id, 1), &
                        tran_to_init_fin_id(tran_id, 2) ) + binned_rate_t(:, :, :, :, :, i)
                end if

                binned_rate_i_scatter(:, :, :, :, :, &
                    tran_to_init_fin_id(tran_id, 1)) = &
                binned_rate_i_scatter(:, :, :, :, :, &
                    tran_to_init_fin_id(tran_id, 1)) + binned_rate_t(:, :, :, :, :, i)
                
                binned_rate_scatter = binned_rate_scatter + binned_rate_t(:, :, :, :, :, i)

                do t = 1, n_time
                    do f = 1, n_FDM
                        do m = 1, n_mX

                            rate_scatter(m, f, t) = rate_scatter(m, f, t) + &
                                sum(binned_rate_t(:, :, m, f, t, i))

                        end do
                    end do
                end do

            end if

        end do

    end subroutine

    ! save data functions
    subroutine save_rates(filename, verbose)

        implicit none

        character(len=*) :: filename

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id

        logical, optional :: verbose
        logical :: file_exists

        integer(HSIZE_T) :: dims1(1) = [1]
        integer(HSIZE_T) :: dims2(2)

        integer :: error

        integer :: m, f, t
        integer :: i, fin
        character(len=64) :: m_str, t_str, f_str
        character(len=64) :: i_str, fin_str

        if ( verbose ) then

            print*, 'Saving scattering rates...'
            print*

        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'rates', group_id, error)

            dims2 = [n_q_bins + 1, n_E_bins + 1]

            do t = 1, n_time
                write(t_str, *) t
                t_str = trim(adjustl(t_str))

                call h5gcreate_f(file_id,&
                    'rates'//&
                    '/t_'//trim(t_str),&
                    group_id, error)

                do f = 1, n_FDM 
                    write(f_str, *) f
                    f_str = trim(adjustl(f_str))

                    call h5gcreate_f(file_id,&
                        'rates'//&
                        '/t_'//trim(t_str)//&
                        '/f_'//trim(f_str),&
                        group_id, error)

                    do m = 1, n_mX
                        write(m_str, *) m
                        m_str = trim(adjustl(m_str))

                        call h5gcreate_f(file_id,& 
                            'rates'//&
                            '/t_'//trim(t_str)//&
                            '/f_'//trim(f_str)//&
                            '/m_'//trim(m_str),&
                            group_id, error)

                        call h5ltmake_dataset_double_f(file_id,&
                            'rates'//&
                            '/t_'//trim(t_str)//&
                            '/f_'//trim(f_str)//&
                            '/m_'//trim(m_str)//&
                            '/total', &
                            size(dims1), dims1,&
                            rate_scatter(m, f, t), error)

                        call h5ltmake_dataset_double_f(file_id,&
                            'rates'//&
                            '/t_'//trim(t_str)//&
                            '/f_'//trim(f_str)//&
                            '/m_'//trim(m_str)//&
                            '/total_binned', &
                            size(dims2), dims2,&
                            binned_rate_scatter(:, :, m, f, t), error)

                        do i = 1, n_init

                            write(i_str, *) i
                            i_str = trim(adjustl(i_str))

                            call h5gcreate_f(file_id,& 
                                'rates'//&
                                '/t_'//trim(t_str)//&
                                '/f_'//trim(f_str)//&
                                '/m_'//trim(m_str)//&
                                '/init_'//trim(i_str),&
                                group_id, error)

                            call h5ltmake_dataset_double_f(file_id,&
                                'rates'//&
                                '/t_'//trim(t_str)//&
                                '/f_'//trim(f_str)//&
                                '/m_'//trim(m_str)//&
                                '/init_'//trim(i_str)//&
                                '/binned_i', &
                                size(dims2), dims2,&
                                binned_rate_i_scatter(:, :, m, f, t, i), error)

                            if ( save_binned_rate_if ) then

                                do fin = 1, n_fin

                                    write(fin_str, *) fin
                                    fin_str = trim(adjustl(fin_str))

                                    call h5gcreate_f(file_id,& 
                                        'rates'//&
                                        '/t_'//trim(t_str)//&
                                        '/f_'//trim(f_str)//&
                                        '/m_'//trim(m_str)//&
                                        '/init_'//trim(i_str)//&
                                        '/fin_'//trim(fin_str),&
                                        group_id, error)

                                    call h5ltmake_dataset_double_f(file_id,&
                                        'rates'//&
                                        '/t_'//trim(t_str)//&
                                        '/f_'//trim(f_str)//&
                                        '/m_'//trim(m_str)//&
                                        '/init_'//trim(i_str)//&
                                        '/fin_'//trim(fin_str)//&
                                        '/binned_if',&
                                        size(dims2), dims2,&
                                        binned_rate_if_scatter(:, :, m, f, t, i, fin), error)

                                end do

                            end if

                        end do
                    end do
                end do
            end do

            call h5fclose_f(file_id, error)
            call h5close_f(error)

        else

            if ( verbose ) then

                print*, '!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*
                print*, '    Output file : ', trim(filename), ' does NOT exist.'
                print*
                print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*

            end if

            stop

        end if

    end subroutine

end module
