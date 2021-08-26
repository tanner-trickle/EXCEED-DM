module exdm_scatter_cf
    !! Compute the core \( \rightarrow \) free DM-electron scattering rate.

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

    use numerics_scatter_cf
    use calc_exdm_scatter_cf

    implicit none

contains

    subroutine exdm_scatter_cf_set_n_init(io_files, n_init)
        !! Sets the number of initial states to record the rate for.
        !! Specific to the core \( \rightarrow \) free DM-electron
        !! scattering rate calculation.
        implicit none

        type(io_files_t) :: io_files
        integer :: n_init

        type(numerics_scatter_cf_t) :: numerics
        type(core_electron_t) :: core_electron

        call core_electron%load(io_files%core_elec_config_filename, io_files%sto_data_filename, &
            verbose = .FALSE.)

        call numerics%load(io_files%nml_input_filename, core_electron, &
            verbose = .FALSE.)

        call numerics%create_core_id_list(core_electron)

        n_init = size(numerics%core_id_list)

    end subroutine

    subroutine run_exdm_scatter_cf(n_init, n_proc, proc_id, root_process, &
        binned_rate_init, main_control, io_files, target_mat, dm_model, &
        bins, expt, in_med_scr, verbose)
        !! Compute the core \( \rightarrow \) free DM-electron scattering rate.
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

        type(core_electron_t) :: core_electron
        type(parallel_manager_t) :: iw_manager
        type(numerics_scatter_cf_t) :: numerics

        integer, allocatable :: job_id_to_iw(:, :)

        integer :: t, f, i, k, j
        integer :: job_id
        integer :: init_id, w_id, core_id

        real(dp) :: log_omega_min, log_omega_max

        if ( verbose ) then

            print*, 'Starting c -> f scattering rate calculation...'
            print*

        end if

        call core_electron%load(&
            io_files%core_elec_config_filename, &
            io_files%STO_data_filename, verbose = verbose)

        call numerics%load(io_files%nml_input_filename, core_electron, &
            verbose = verbose)
        
        ! parallelize over {i, omega} variables
        call iw_manager%init(n_init*numerics%n_omega, verbose = verbose)

        allocate(job_id_to_iw(iw_manager%n_jobs, 4))
        job_id_to_iw = 0

        call numerics%create_core_id_list(core_electron)
        call numerics%create_w_id_list()

        log_omega_min = log10(numerics%Ef_min)
        log_omega_max = log10(0.5_dp*dm_model%vX_max**2*maxval(dm_model%mX))

        numerics%omega_list = 10.0_dp**uniform_list(numerics%n_omega, log_omega_min, log_omega_max)

        call iw_manager%create_job_to_2d_ID_table(&
            numerics%core_id_list, &
            numerics%w_id_list, &
            job_id_to_iw, verbose = verbose)

        if ( log_omega_max > log_omega_min ) then

            if ( verbose ) then
                print*, 'Calculating transition rates...'
                print*
            end if

            do j = 1, iw_manager%n_jobs_per_proc
                
                job_id = iw_manager%job_table(proc_id + 1, j)

                if ( job_id /= 0 ) then

                    ! count states down from highest valence band
                    core_id = job_id_to_iw(job_id, 1)
                    w_id = job_id_to_iw(job_id, 2)
                    init_id = job_id_to_iw(job_id, 3)

                    ! compute rate
                    call exdm_scatter_cf_calc(binned_rate_init(init_id), &
                        core_electron, target_mat, &
                        bins, dm_model, expt, in_med_scr, numerics, &
                        core_id, w_id, verbose = .FALSE.)

                end if

            end do

            if ( verbose ) then
                print*, 'Done calculating transition rates!'
                print*
            end if

        else

            call print_warning_message(&
                'No mass in mX is large enough to warrant the c -> f calculation '//&
                'with the specified omega_min. Skipping c -> f calculation.', &
                verbose = verbose)

        end if

        if ( proc_id == root_process ) then
            call core_electron%save(io_files%out_filename, verbose = verbose)
            call numerics%save(io_files%out_filename, verbose = verbose)
        end if

    end subroutine

    ! subroutine time_exdm_scatter_cf_calc(log_omega_table, &
    !         ki_angular_mesh, kf_angular_mesh, verbose)
    !     !! Times the v -> f scattering rate calculation
    !     use timing 
    !     use mpi

    !     implicit none

    !     real(dp) :: log_omega_table(n_fin, 2)
    !     real(dp) :: log_omegas(2)

    !     real(dp) :: ki_angular_mesh(:, :)
    !     real(dp) :: kf_angular_mesh(:, :)

    !     integer :: tran_id

    !     logical, optional :: verbose
    !     real(dp) :: b_rate(n_q_bins + 1, n_E_bins + 1, n_mX, n_FDM, n_time)

    !     integer :: init_id, fin_id

    !     if ( verbose ) then

    !         print*, 'Timing c -> f calculation...'
    !         print*

    !     end if

    !     ! one of the jobs which takes longer
    !     tran_id = job_table(12, n_tran_per_proc)

    !     init_id = tran_to_init_fin_id(tran_id, 1)
    !     fin_id = tran_to_init_fin_id(tran_id, 2)

    !     log_omegas = log_omega_table(fin_id, :)

    !     time(3) = MPI_Wtime()

    !     call exdm_scatter_cf_calc(b_rate,& 
    !         init_id, log_omegas, n_ki, ki_angular_mesh, kf_angular_mesh, verbose = verbose)

    !     time(4) = MPI_Wtime()

    !     if ( verbose ) then

    !         print*, '----------------------------------------'
    !         print*, '    -------------'
    !         print*, '    Timing (TEST)'
    !         print*, '    -------------'
    !         print*
    !         print*, '        (TEST) Run time : '
    !         print*, '            ', trim(pretty_time_format(time(4) - time(3)))
    !         print*
    !         print*, '        Expected run time for whole calculation :'
    !         print*, '            ', trim(pretty_time_format(&
    !             n_tran_per_proc*(time(4) - time(3))&
    !             ))
    !         print*
    !         print*, '----------------------------------------'
    !         print*

    !     end if

    ! end subroutine

end module
