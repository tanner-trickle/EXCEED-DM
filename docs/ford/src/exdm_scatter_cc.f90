module exdm_scatter_cc
    !! Compute the core \( \rightarrow \) conduction DM-electron scattering rate.
    
    use mpi

    use prec
    use info_messages

    use io_input
    use control_input

    use binned_scatter_rate_type
    use dm_model_type
    use expt_type
    use material_type
    use PW_dataset_type
    use core_electron_type
    use in_med_scr_type

    use FFT_util
    use MPI_util

    use numerics_scatter_cc
    use calc_exdm_scatter_cc

    implicit none

contains

    subroutine exdm_scatter_cc_set_n_init(io_files, n_init)
        !! Sets the number of initial states to record the rate for.
        !! Specific to the core \( \rightarrow \) conduction DM-electron
        !! scattering rate calculation.
        implicit none

        type(io_files_t) :: io_files
        integer :: n_init

        type(PW_dataset_t) :: PW_dataset
        type(numerics_scatter_cc_t) :: numerics
        type(core_electron_t) :: core_electron

        call PW_dataset%load(io_files%PW_data_filename, verbose = .FALSE.)

        call core_electron%load(io_files%core_elec_config_filename, io_files%sto_data_filename, &
            verbose = .FALSE.)

        call numerics%load(io_files%nml_input_filename, PW_dataset%n_cond, verbose = .FALSE.)

        call numerics%create_core_id_list(core_electron)

        n_init = size(numerics%core_id_list)

    end subroutine

    subroutine run_exdm_scatter_cc(n_init, n_proc, proc_id, root_process, &
            binned_rate_init, main_control, io_files, target_mat, dm_model, &
            bins, expt, in_med_scr, verbose)
        !! Compute the core \( \rightarrow \) conduction DM-electron scattering rate.
        implicit none

        integer :: n_init
        integer :: n_proc
        integer :: proc_id, root_process
        type(binned_scatter_rate_t) :: binned_rate_init(n_init)
        type(control_t) :: main_control
        type(io_files_t) :: io_files
        type(material_t) :: target_mat
        type(dm_model_t) :: dm_model
        type(bins_scatter_t) :: bins
        type(expt_t) :: expt
        type(in_med_scr_t) :: in_med_scr
        logical, optional :: verbose

        type(PW_dataset_t) :: PW_dataset
        type(FFT_grid_t) :: fft_grid
        type(parallel_manager_t) :: ikf_manager
        type(numerics_scatter_cc_t) :: numerics

        type(core_electron_t) :: core_electron

        type(binned_scatter_rate_t), allocatable :: binned_rate_job(:)
            !! Dim : [ n_tran_per_proc ]
            !!
            !! Binned rate per cross section, per kg-year for each transition
            !!
            !! Units : cm^(-2)

        integer, allocatable :: job_id_to_ikf(:, :)

        integer :: init_id

        integer :: n_FFT_grid(3)

        integer :: n_jobs_per_proc
        integer :: i, j, f

        integer :: job_id
        integer :: cond_id, kf

        complex(dp), allocatable :: wfc_fkf(:, :, :)
        complex(dp), allocatable :: wfc_i(:, :, :)

        if ( verbose ) then
            print*, 'Starting c -> c scattering rate calculation...'
            print*
        end if

        call PW_dataset%load(io_files%PW_data_filename, verbose = verbose)
        call PW_dataset%do_scissor_correction(target_mat%band_gap, verbose = verbose)

        call core_electron%load(&
            io_files%core_elec_config_filename, &
            io_files%STO_data_filename, verbose = verbose)

        call numerics%load(io_files%nml_input_filename,&
           PW_dataset%n_cond, verbose = verbose)

        ! at least doubled to avoid wrapping problems
        do i = 1, 3
            n_FFT_grid(i) = max(&
                    2*(&
                        maxval(PW_dataset%G_grid_red(:, i)) &
                        - minval(PW_dataset%G_grid_red(:, i))&
                    ) + 1, &
                    numerics%n_FFT_grid(i)&
                )
        end do

        call FFT_grid%init(n_FFT_grid, PW_dataset%k_red_to_xyz, 'b', verbose = verbose)
        call FFT_grid%print(verbose = verbose)

        ! create job table, parallelizing over {i, kf}
        call ikf_manager%init(PW_dataset%n_k*n_init, verbose = verbose)

        allocate(job_id_to_ikf(ikf_manager%n_jobs, 4))
        job_id_to_ikf = 0

        call numerics%create_core_id_list(core_electron)
        call numerics%create_k_id_list(PW_dataset%n_k)

        call ikf_manager%create_job_to_2d_ID_table(&
            numerics%core_id_list, &
            numerics%k_id_list, &
            job_id_to_ikf, verbose = verbose)

        ! allocate the binned rate arrays
        allocate(binned_rate_job(ikf_manager%n_jobs_per_proc))
        do i = 1, ikf_manager%n_jobs_per_proc
            call binned_rate_job(i)%init(bins, dm_model, expt)
        end do

        ! time calculation

        if ( ( proc_id == root_process ) .and. ( main_control%timer ) ) then

            call time_exdm_scatter_cc_calc(FFT_grid, core_electron, PW_dataset, target_mat, &
                        bins, dm_model, expt, in_med_scr, numerics, ikf_manager%n_jobs_per_proc, &
                        core_electron%n_state, PW_dataset%n_val + 1, 1, verbose = verbose)

        end if

        ! allocate wave functions
        allocate(wfc_i(FFT_grid%n_grid(1), FFT_grid%n_grid(2), FFT_grid%n_grid(3)))
        allocate(wfc_fkf(FFT_grid%n_grid(1), FFT_grid%n_grid(2), FFT_grid%n_grid(3)))

        ! do calculation
        if ( verbose ) then
            print*, 'Calculating transition rates...'
            print*
        end if

        do j = 1, ikf_manager%n_jobs_per_proc
            
            job_id = ikf_manager%job_table(proc_id + 1, j)

            if ( job_id /= 0 ) then

                init_id = job_id_to_ikf(job_id, 1)
                kf = job_id_to_ikf(job_id, 2)

                call core_electron%bloch_wf_on_grid(FFT_grid%n_grid, &
                    wfc_i, init_id, target_mat%pc_vol, PW_dataset%red_to_xyz, &
                    shift = .TRUE., verbose = .FALSE.)

                do f = 1, numerics%n_cond_max

                    cond_id = PW_dataset%n_val + f
                    ! load final wave function
                    call PW_dataset%load_wfc_ik_expanded_no_spin(cond_id, kf, FFT_grid, wfc_fkf)

                    ! compute rate
                    call exdm_scatter_cc_calc(binned_rate_job(j), &
                        FFT_grid, core_electron, PW_dataset, target_mat, &
                        bins, dm_model, expt, in_med_scr, &
                        wfc_i, wfc_fkf, &
                        init_id, cond_id, kf, verbose = .FALSE.)

                end do

            end if

        end do

        if ( verbose ) then
            print*, 'Done calculating transition rates!'
            print*
        end if

        call ikf_manager%comm_scatter_binned_rate_job_init(proc_id, root_process, job_id_to_ikf, &
            binned_rate_job, binned_rate_init, verbose)

        if ( proc_id == root_process ) then
            call PW_dataset%save(io_files%out_filename, verbose = verbose)
            call core_electron%save(io_files%out_filename, verbose = verbose)
            call fft_grid%save(io_files%out_filename, verbose = verbose)
            call numerics%save(io_files%out_filename, verbose = verbose)
        end if

    end subroutine

    subroutine time_exdm_scatter_cc_calc(FFT_grid, core_electron, PW_dataset, target_mat, &
                        bins, dm_model, expt, in_med_scr, numerics, n_jobs_per_proc, &
                        init_id, cond_id, kf, verbose)
        !! Clocks the core \( \rightarrow \) conduction DM-electron scattering rate calculation
        !! by running a smaller version of the program.
        use timing 
        use mpi

        implicit none

        type(FFT_grid_t) :: FFT_grid
        type(core_electron_t) :: core_electron
        type(PW_dataset_t) :: PW_dataset
        type(material_t) :: target_mat
        type(bins_scatter_t) :: bins
        type(dm_model_t) :: dm_model
        type(expt_t) :: expt
        type(in_med_scr_t) :: in_med_scr
        type(numerics_scatter_cc_t) :: numerics
        integer :: n_jobs_per_proc
        integer :: init_id, cond_id, k, kf
        logical, optional :: verbose

        type(binned_scatter_rate_t) :: b_rate

        complex(dp), allocatable :: wfc_i(:, :, :)
        complex(dp), allocatable :: wfc_fkf(:, :, :)

        if ( verbose ) then
            print*, 'Timing c -> c calculation...'
            print*
        end if

        call b_rate%init(bins, dm_model, expt)

        allocate(wfc_i(FFT_grid%n_grid(1), FFT_grid%n_grid(2), FFT_grid%n_grid(3)))
        allocate(wfc_fkf(FFT_grid%n_grid(1), FFT_grid%n_grid(2), FFT_grid%n_grid(3)))

        time(3) = MPI_Wtime()

        call core_electron%bloch_wf_on_grid(FFT_grid%n_grid, &
            wfc_i, init_id, target_mat%pc_vol, PW_dataset%red_to_xyz, &
            shift = .TRUE., verbose = .FALSE.)

        time(4) = MPI_Wtime()

        call PW_dataset%load_wfc_ik_expanded_no_spin(cond_id, kf, FFT_grid, wfc_fkf)

        call exdm_scatter_cc_calc(b_rate,&
            FFT_grid, core_electron, PW_dataset, target_mat, &
            bins, dm_model, expt, in_med_scr, &
            wfc_i, wfc_fkf, &
            init_id, cond_id, kf, verbose = .FALSE.)

        time(5) = MPI_Wtime()

        if ( verbose ) then
            call print_section_seperator()
            print*, '    -------------'
            print*, '    Timing (TEST)'
            print*, '    -------------'
            print*
            print*, '        (TEST) Compute core WF :'
            print*, '            ', trim(pretty_time_format(time(4) - time(3)))
            print*
            print*, '        (Total) Compute core WF :'
            print*, '            ', trim(pretty_time_format(&
                n_jobs_per_proc*(time(4) - time(3))&
                ))
            print*
            print*, '        (TEST) Rate : '
            print*, '            ', trim(pretty_time_format(time(5) - time(4)))
            print*
            print*, '        (Total) Rate : '
            print*, '            ', trim(pretty_time_format(&
                n_jobs_per_proc*numerics%n_cond_max*(time(5) - time(4))&
                ))
            print*
            print*, '        Expected run time for whole calculation :'
            print*, '            ', trim(pretty_time_format(&
                n_jobs_per_proc*(time(4) - time(3)) &
                + n_jobs_per_proc*numerics%n_cond_max*(time(5) - time(4))&
                ))
            print*
            call print_section_seperator()
            print*

        end if

    end subroutine

end module
