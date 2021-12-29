module exdm_scatter_vc
    !! Compute the valence \( \rightarrow \) conduction DM-electron scattering rate.

    use mpi

    use prec

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

    use numerics_scatter_vc
    use calc_exdm_scatter_vc

    implicit none

contains

    subroutine exdm_scatter_vc_set_n_init(io_files, n_init)
        !! Sets the number of initial states to record the rate for.
        !! Specific to the valence \( \rightarrow \) conduction DM-electron
        !! scattering rate calculation.
        implicit none

        type(io_files_t) :: io_files
        integer :: n_init

        type(numerics_scatter_vc_t) :: numerics
        type(PW_dataset_t) :: PW_dataset

        call PW_dataset%load(io_files%PW_data_filename, verbose = .FALSE.)

        call numerics%load(io_files%nml_input_filename, &
            PW_dataset%n_val, PW_dataset%n_cond, verbose = .FALSE.)

        call numerics%create_val_id_list(PW_dataset%n_val)

        n_init = size(numerics%val_id_list)

    end subroutine

    subroutine run_exdm_scatter_vc(n_init, n_proc, proc_id, root_process, &
        binned_rate_init, main_control, io_files, target_mat, dm_model, &
        bins, expt, in_med_scr, verbose)
        !! Compute the valence \( \rightarrow \) conduction DM-electron scattering rate.
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
        type(parallel_manager_t) :: ik_manager
        type(numerics_scatter_vc_t) :: numerics

        integer, allocatable :: job_id_to_ik(:, :)

        integer :: i, j, f

        integer :: job_id
        integer :: val_id, cond_id, init_id
        integer :: k, kf

        integer :: n_FFT_grid(3)

        complex(dp), allocatable :: wfc_ik(:, :, :)
        complex(dp), allocatable :: wfc_iks(:, :, :, :)

        complex(dp), allocatable :: wfc_fkf(:, :, :)
        complex(dp), allocatable :: wfc_fkfs(:, :, :, :)

        if ( verbose ) then
            print*, 'Starting v -> c scattering rate calculation...'
            print*
        end if

        call PW_dataset%load(io_files%PW_data_filename, verbose = verbose)
        call PW_dataset%do_scissor_correction(target_mat%band_gap, verbose = verbose)

        call numerics%load(io_files%nml_input_filename, &
            PW_dataset%n_val, PW_dataset%n_cond, verbose = verbose)

        ! doubled to avoid wrapping problems
        do i = 1, 3
            n_FFT_grid(i) = 2*(&
                maxval(PW_dataset%G_grid_red(:, i)) - minval(PW_dataset%G_grid_red(:, i))&
                ) + 1
        end do

        call FFT_grid%init(n_FFT_grid, PW_dataset%k_red_to_xyz, 'b', verbose = verbose)
        call FFT_grid%print(verbose = verbose)

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

        ! time calculation
        ! if ( ( proc_id == root_process ) .and. ( main_control%timer ) ) then

        !     call time_exdm_scatter_vc_calc(FFT_grid, PW_dataset, target_mat, &
        !                         bins, dm_model, expt, in_med_scr, &
        !                         numerics, ik_manager%n_jobs_per_proc, &
        !                         PW_dataset%n_val, PW_dataset%n_val + 1, 1, 1, verbose = verbose)

        ! end if

        ! allocate wave functions
        if ( PW_dataset%include_spin ) then

            allocate(wfc_iks(2, FFT_grid%n_grid(1), FFT_grid%n_grid(2), FFT_grid%n_grid(3)))
            allocate(wfc_fkfs(2, FFT_grid%n_grid(1), FFT_grid%n_grid(2), FFT_grid%n_grid(3)))

        else

            allocate(wfc_ik(FFT_grid%n_grid(1), FFT_grid%n_grid(2), FFT_grid%n_grid(3)))
            allocate(wfc_fkf(FFT_grid%n_grid(1), FFT_grid%n_grid(2), FFT_grid%n_grid(3)))

        end if

        ! do calculation
        if ( verbose ) then
            print*, 'Calculating transition rates...'
            print*
        end if

        do j = 1, ik_manager%n_jobs_per_proc
            
            job_id = ik_manager%job_table(proc_id + 1, j)

            if ( job_id /= 0 ) then

                ! count states down from highest valence band
                val_id = job_id_to_ik(job_id, 1)
                init_id = job_id_to_ik(job_id, 3)
                k = job_id_to_ik(job_id, 2)

                if ( PW_dataset%include_spin ) then

                    ! load initial wave function
                    call PW_dataset%load_wfc_ik_expanded_spin(val_id, k, FFT_grid, wfc_iks)

                    do f = 1, numerics%n_cond_max

                        cond_id = PW_dataset%n_val + f

                        do kf = 1, PW_dataset%n_k

                            ! load final wave function
                            call PW_dataset%load_wfc_ik_expanded_spin(cond_id, kf, FFT_grid, wfc_fkfs)

                            ! compute rate
                            call exdm_scatter_vc_calc(binned_rate_init(init_id), &
                                FFT_grid, PW_dataset, target_mat, &
                                bins, dm_model, expt, in_med_scr, &
                                wfc_iks, wfc_fkfs, &
                                val_id, cond_id, k, kf, verbose = .FALSE.)

                        end do
                    end do

                else

                    ! load initial wave function
                    call PW_dataset%load_wfc_ik_expanded_no_spin(val_id, k, FFT_grid, wfc_ik)

                    do f = 1, numerics%n_cond_max

                        cond_id = PW_dataset%n_val + f

                        do kf = 1, PW_dataset%n_k

                            ! load final wave function
                            call PW_dataset%load_wfc_ik_expanded_no_spin(cond_id, kf, FFT_grid, wfc_fkf)

                            ! compute rate
                            call exdm_scatter_vc_calc(binned_rate_init(init_id), &
                                FFT_grid, PW_dataset, target_mat, &
                                bins, dm_model, expt, in_med_scr, &
                                wfc_ik, wfc_fkf, &
                                val_id, cond_id, k, kf, verbose = .FALSE.)

                        end do
                    end do

                end if

            end if

        end do

        if ( verbose ) then
            print*, 'Done calculating transition rates!'
            print*
        end if

        if ( proc_id == root_process ) then
            call PW_dataset%save(io_files%out_filename, verbose = verbose)
            call fft_grid%save(io_files%out_filename, verbose = verbose)
            call numerics%save(io_files%out_filename, verbose = verbose)
        end if

    end subroutine

    ! subroutine time_exdm_scatter_vc_calc(FFT_grid, PW_dataset, target_mat, &
    !         bins, dm_model, expt, in_med_scr, &
    !         numerics, n_jobs_per_proc, val_id, cond_id, k, kf, verbose)
    !     !! Clocks the valence \( \rightarrow \) conduction DM-electron scattering rate calculation
    !     !! by running a smaller version of the program.
    !     use timing 
    !     use mpi

    !     implicit none

    !     type(FFT_grid_t) :: FFT_grid
    !     type(PW_dataset_t) :: PW_dataset
    !     type(material_t) :: target_mat
    !     type(bins_scatter_t) :: bins
    !     type(dm_model_t) :: dm_model
    !     type(expt_t) :: expt
    !     type(in_med_scr_t) :: in_med_scr
    !     type(numerics_scatter_vc_t) :: numerics
    !     integer :: n_jobs_per_proc
    !     integer :: val_id, cond_id, k, kf
    !     logical, optional :: verbose

    !     complex(dp), allocatable :: wfc_ik(:, :, :)
    !     complex(dp), allocatable :: wfc_iks(:, :, :, :)

    !     complex(dp), allocatable :: wfc_fkf(:, :, :)
    !     complex(dp), allocatable :: wfc_fkfs(:, :, :, :)

    !     type(binned_scatter_rate_t) :: b_rate

    !     if ( verbose ) then
    !         print*, 'Timing v -> c calculation...'
    !         print*
    !     end if

    !     if ( PW_dataset%include_spin ) then

    !         allocate(wfc_iks(2, FFT_grid%n_grid(1), FFT_grid%n_grid(2), FFT_grid%n_grid(3)))
    !         allocate(wfc_fkfs(2, FFT_grid%n_grid(1), FFT_grid%n_grid(2), FFT_grid%n_grid(3)))

    !         call PW_dataset%load_wfc_ik_expanded_spin(val_id, k, FFT_grid, wfc_iks)

    !     else

    !         allocate(wfc_ik(FFT_grid%n_grid(1), FFT_grid%n_grid(2), FFT_grid%n_grid(3)))
    !         allocate(wfc_fkf(FFT_grid%n_grid(1), FFT_grid%n_grid(2), FFT_grid%n_grid(3)))

    !         call PW_dataset%load_wfc_ik_expanded_no_spin(val_id, k, FFT_grid, wfc_ik)

    !     end if

    !     call b_rate%init(bins, dm_model, expt)

    !     time(3) = MPI_Wtime()

    !     if ( PW_dataset%include_spin ) then

    !         call PW_dataset%load_wfc_ik_expanded_spin(cond_id, kf, FFT_grid, wfc_fkfs)

    !         call exdm_scatter_vc_calc(b_rate, &
    !             FFT_grid, PW_dataset, target_mat, &
    !             bins, dm_model, expt, in_med_scr, &
    !             wfc_iks, wfc_fkfs, &
    !             val_id, cond_id, k, kf, verbose = .FALSE.)

    !     else

    !         call PW_dataset%load_wfc_ik_expanded_no_spin(cond_id, kf, FFT_grid, wfc_fkf)

    !         call exdm_scatter_vc_calc(b_rate, &
    !             FFT_grid, PW_dataset, target_mat, &
    !             bins, dm_model, expt, in_med_scr, &
    !             wfc_ik, wfc_fkf, &
    !             val_id, cond_id, k, kf, verbose = .FALSE.)

    !     end if

    !     time(4) = MPI_Wtime()

    !     if ( verbose ) then
    !         call print_section_seperator()
    !         print*, '    -------------'
    !         print*, '    Timing (TEST)'
    !         print*, '    -------------'
    !         print*
    !         print*, '        (TEST) Run time: '
    !         print*, '            ', trim(pretty_time_format(time(4) - time(3)))
    !         print*
    !         print*, '        Expected run time for whole calculation : '
    !         print*, '            ', trim(pretty_time_format(&
    !             n_jobs_per_proc*PW_dataset%n_k*numerics%n_cond_max*(time(4) - time(3))&
    !             ))
    !         print*
    !         call print_section_seperator()
    !         print*
    !     end if

    ! end subroutine

end module
