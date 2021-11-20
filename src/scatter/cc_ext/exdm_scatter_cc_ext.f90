module exdm_scatter_cc_ext
    !! Compute the core \( \rightarrow \) conduction DM-electron scattering rate by computing \( \frac{dR}{dq} \).

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

    use numerics_scatter_cc_ext
    ! use calc_exdm_scatter_cc_ext

contains

    subroutine exdm_scatter_cc_ext_set_n_init(io_files, n_init)
        !! Sets the number of initial states to record the rate for.
        !! Specific to the core \( \rightarrow \) conduction (extended) DM-electron
        !! scattering rate calculation.
        implicit none

        type(io_files_t) :: io_files
        integer :: n_init

        type(PW_dataset_t) :: PW_dataset
        type(numerics_scatter_cc_ext_t) :: numerics
        type(core_electron_t) :: core_electron

        call PW_dataset%load(io_files%PW_data_filename, verbose = .FALSE.)

        call core_electron%load(io_files%core_elec_config_filename, io_files%sto_data_filename, &
            verbose = .FALSE.)

        call numerics%load(io_files%nml_input_filename, PW_dataset%n_cond, verbose = .FALSE.)

        call numerics%create_core_id_list(core_electron)

        n_init = size(numerics%core_id_list)

    end subroutine

    subroutine run_exdm_scatter_cc_ext(n_init, n_proc, proc_id, root_process, &
            binned_rate_init, main_control, io_files, target_mat, dm_model, &
            bins, expt, in_med_scr, verbose)
        !! Compute the core \( \rightarrow \) conduction DM-electron scattering rate by computing \( \frac{dR}{dq} \).

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
        type(parallel_manager_t) :: ikf_manager
        type(numerics_scatter_cc_ext_t) :: numerics

        type(core_electron_t) :: core_electron

        integer, allocatable :: job_id_to_ikf(:, :)

        integer :: init_id, core_id

        integer :: n_jobs_per_proc
        integer :: i, j, f

        integer :: job_id
        integer :: cond_id, kf

        complex(dp), allocatable :: wfc_FT_fkf(:)

        if ( verbose ) then
            print*, 'Starting c -> c (extended) scattering rate calculation...'
            print*
        end if

        call PW_dataset%load(io_files%PW_data_filename, verbose = verbose)
        call PW_dataset%do_scissor_correction(target_mat%band_gap, verbose = verbose)

        call core_electron%load(&
            io_files%core_elec_config_filename, &
            io_files%STO_data_filename, verbose = verbose)

        call numerics%load(io_files%nml_input_filename,&
           PW_dataset%n_cond, verbose = verbose)

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

        allocate(wfc_FT_fkf(PW_dataset%n_G))

        if ( verbose ) then
            print*, 'Calculating transition rates...'
            print*
        end if

        do j = 1, ikf_manager%n_jobs_per_proc
            
            job_id = ikf_manager%job_table(proc_id + 1, j)

            if ( job_id /= 0 ) then

                core_id = job_id_to_ikf(job_id, 1)
                kf = job_id_to_ikf(job_id, 2)
                init_id = job_id_to_ikf(job_id, 3)

                do f = 1, numerics%n_cond_max

                    cond_id = PW_dataset%n_val + f

                    ! load final wave function
                    call PW_dataset%load_wfc_FT_ik_no_spin(cond_id, kf, wfc_FT_fkf)

                    ! compute rate
                    ! call exdm_scatter_cc_ext_calc(binned_rate_init(init_id), &
                    !     core_electron, PW_dataset, target_mat, &
                    !     bins, dm_model, expt, in_med_scr, &
                    !     wfc_FT_fkf, &
                    !     core_id, cond_id, kf, verbose = .FALSE.)

                end do

            end if

        end do

        if ( verbose ) then
            print*, 'Done calculating transition rates!'
            print*
        end if

        if ( proc_id == root_process ) then
            call PW_dataset%save(io_files%out_filename, verbose = verbose)
            call core_electron%save(io_files%out_filename, verbose = verbose)
            call numerics%save(io_files%out_filename, verbose = verbose)
        end if

    end subroutine

end module
