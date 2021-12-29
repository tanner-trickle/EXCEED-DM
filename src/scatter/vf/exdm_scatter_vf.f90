module exdm_scatter_vf
    !! Compute the valence \( \rightarrow \) free DM-electron scattering rate.

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
    use in_med_scr_type

    use FFT_util
    use MPI_util

    use numerics_scatter_vf
    use calc_exdm_scatter_vf

    implicit none

contains

    subroutine exdm_scatter_vf_set_n_init(io_files, n_init)
        !! Sets the number of initial states to record the rate for.
        !! Specific to the valence \( \rightarrow \) free DM-electron
        !! scattering rate calculation.
        implicit none

        type(io_files_t) :: io_files
        integer :: n_init

        type(numerics_scatter_vf_t) :: numerics
        type(PW_dataset_t) :: PW_dataset

        call PW_dataset%load(io_files%PW_data_filename, verbose = .FALSE.)

        call numerics%load(io_files%nml_input_filename, &
            PW_dataset, verbose = .FALSE.)

        call numerics%create_val_id_list(PW_dataset%n_val)

        n_init = size(numerics%val_id_list)

    end subroutine

    subroutine run_exdm_scatter_vf(n_init, n_proc, proc_id, root_process, &
        binned_rate_init, main_control, io_files, target_mat, dm_model, &
        bins, expt, in_med_scr, verbose)
        !! Compute the valence \( \rightarrow \) free DM-electron scattering rate.
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
        type(parallel_manager_t) :: ik_manager
        type(numerics_scatter_vf_t) :: numerics

        integer, allocatable :: job_id_to_ik(:, :)

        integer :: t, f, i, k, j
        integer :: job_id
        integer :: val_id, fin_id
        integer :: init_id

        complex(dp), allocatable :: wfc_FT_ik(:)

        real(dp) :: log_omega_min, log_omega_max

        if ( verbose ) then
            print*, 'Starting v -> f scattering rate calculation...'
            print*
        end if

        call PW_dataset%load(io_files%PW_data_filename, verbose = verbose)
        call PW_dataset%do_scissor_correction(target_mat%band_gap, verbose = verbose)

        call numerics%load(io_files%nml_input_filename, PW_dataset, verbose = verbose)

        ! create job table, parallelizing over {i, k}
        call ik_manager%init(PW_dataset%n_k*numerics%n_val_max, verbose = verbose)

        allocate(job_id_to_ik(ik_manager%n_jobs, 4))
        job_id_to_ik = 0

        call numerics%create_val_id_list(PW_dataset%n_val)
        call numerics%create_k_id_list(PW_dataset%n_k)

        call ik_manager%create_job_to_2d_ID_table(&
            numerics%val_id_list, &
            numerics%k_id_list, &
            job_id_to_ik, verbose = verbose)

        allocate(wfc_FT_ik(PW_dataset%n_G))

        log_omega_min = log10(PW_dataset%Ef_max)
        log_omega_max = log10(0.5_dp*dm_model%vX_max**2*maxval(dm_model%mX))

        if ( log_omega_max > log_omega_min ) then

            numerics%omega_list = 10.0_dp**uniform_list(numerics%n_omega, log_omega_min, log_omega_max)

            if ( verbose ) then
                print*, 'Calculating transition rates...'
                print*
            end if

            do j = 1, ik_manager%n_jobs_per_proc
                
                job_id = ik_manager%job_table(proc_id + 1, j)

                if ( job_id /= 0 ) then

                    ! count states down from highest valence band
                    val_id = job_id_to_ik(job_id, 1)
                    k = job_id_to_ik(job_id, 2)
                    init_id = job_id_to_ik(job_id, 3)

                    ! load initial wave function
                    call PW_dataset%load_wfc_FT_ik_no_spin(val_id, k, wfc_FT_ik)

                    ! compute rate
                    call exdm_scatter_vf_calc(binned_rate_init(init_id), &
                        PW_dataset, target_mat, &
                        bins, dm_model, expt, in_med_scr, &
                        numerics, &
                        wfc_FT_ik, val_id, k, verbose = .FALSE.)

                end if

            end do

            if ( verbose ) then
                print*, 'Done calculating transition rates!'
                print*
            end if

        else

            call print_warning_message(&
                'No mass in mX is large enough to warrant the v -> f calculation'//&
                ' with the specified Ef_max. Skipping v -> f calculation.', &
                verbose = verbose)

        end if

        if ( proc_id == root_process ) then
            call PW_dataset%save(io_files%out_filename, verbose = verbose)
            call numerics%save(io_files%out_filename, verbose = verbose)
        end if

    end subroutine

end module
