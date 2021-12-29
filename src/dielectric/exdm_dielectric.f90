module exdm_dielectric
    !! Compute the dielectric function, \( \varepsilon(\mathbf{q}, \omega) \).
    !!
    !! Note : currently only includes valence \( \rightarrow \) conduction contributions.
    !!
    !! Equations:
    !!
    !!$$\begin{align*}
    !!      \varepsilon(\mathbf{q}, \omega) & = 1 - \frac{e^2}{q^2} \Pi_{11}(\mathbf{q}, \omega) \\
    !!      \Pi_{11} & = \frac{1}{V} \sum_{II'} G(\omega, \Delta \omega_{II'})
    !!      | \langle I' | e^{i \mathbf{q} \cdot \mathbf{x}} | I \rangle |^2 \\
    !!      G(\omega, \Delta \omega) & = \frac{1}{\omega - \Delta \omega + i \delta} 
    !!      - \frac{1}{\omega + \Delta \omega - i \delta }
    !! \end{align*}$$

    use mpi

    use prec

    use io_input
    use control_input

    use bins_dielectric_type
    use dm_model_type
    use expt_type
    use material_type
    use PW_dataset_type
    use width_parameters_type

    use FFT_util
    use MPI_util

    use numerics_dielectric
    use calc_dielectric_vc

    implicit none

contains

    subroutine run_exdm_dielectric(proc_id, root_process, n_proc, &
            io_files, main_control, target_mat, verbose)
        !! Compute the dielectric function, \( \varepsilon(\mathbf{q}, \omega) \).
        !!
        !! Note : currently only includes valence \( \rightarrow \) conduction contributions.
        !!
        !! Equations:
        !!
        !!$$\begin{align*}
        !!      \varepsilon(\mathbf{q}, \omega) & = 1 - \frac{e^2}{q^2} \Pi_{11}(\mathbf{q}, \omega) \\
        !!      \Pi_{11} & = \frac{1}{V} \sum_{II'} G(\omega, \Delta \omega_{II'})
        !!      | \langle I' | e^{i \mathbf{q} \cdot \mathbf{x}} | I \rangle |^2 \\
        !!      G(\omega, \Delta \omega) & = \frac{1}{\omega - \Delta \omega + i \delta} 
        !!      - \frac{1}{\omega + \Delta \omega - i \delta }
        !! \end{align*}$$

        implicit none

        integer :: proc_id, root_process, n_proc
        type(io_files_t) :: io_files
        type(control_t) :: main_control
        type(material_t) :: target_mat
        logical, optional :: verbose

        logical :: file_exists

        type(PW_dataset_t) :: PW_dataset
        type(bins_dielectric_t) :: bins
        type(numerics_dielectric_t) :: numerics
        type(FFT_grid_t) :: FFT_grid
        type(parallel_manager_t) :: ik_manager
        type(width_parameters_t) :: widths

        integer :: i
        integer :: n_FFT_grid(3)

        integer, allocatable :: job_id_to_ik(:, :)

        complex(dp), allocatable :: wfc_ik(:, :, :)
        complex(dp), allocatable :: wfc_iks(:, :, :, :)

        complex(dp), allocatable :: wfc_fkf(:, :, :)
        complex(dp), allocatable :: wfc_fkfs(:, :, :, :)

        integer :: job_id, j
        integer :: val_id, cond_id, f, k, kf

        complex(dp), allocatable :: dielec(:, :, :, :)
            !! Dim : [ bins%n_E, bins%n_q, bins%n_q_theta, bins%n_q_phi ]
            !!
            !!$$\begin{align*}
            !!      \eps(\mathbf{q}, \omega) & = 1 - \frac{e^2}{q^2} \Pi_{11}(\mathbf{q}, \omega)
            !!      \Pi_{11} & = \frac{1}{V} \sum_{II'} G(\omega, \Delta \omega_{II'})
            !!      | \langle I' | e^{i \mathbf{q} \cdot \mathbf{x} | I \rangle |^2 \\
            !!      G(\omega, \Delta \omega) & = \frac{1}{\omega - \Delta \omega + i \delta} 
            !!      - \frac{1}{\omega + \Delta \omega - i \delta }
            !! \end{align*}
            !!
            !! Units : None

        character(len=512) :: dielectric_output_filename

        if ( trim(main_control%process) == 'dielectric' ) then

            dielectric_output_filename = io_files%out_filename

        else

            dielectric_output_filename = io_files%dielectric_filename

        end if

        call PW_dataset%load(io_files%PW_data_filename, verbose = verbose)
        call PW_dataset%do_scissor_correction(target_mat%band_gap, verbose = verbose)

        call numerics%load(io_files%nml_input_filename, &
            PW_dataset%n_val, PW_dataset%n_cond, verbose = verbose)

        if ( verbose ) then
            print*, 'Starting dielectric calculation...'
            print*
        end if

        call bins%load(io_files%nml_input_filename, verbose = verbose)

        call widths%load(io_files%nml_input_filename, verbose = verbose)

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

        call check_dielectric_memory(ik_manager%n_jobs_per_proc*&
            bins%n_E*bins%n_q*bins%n_q_theta*bins%n_q_phi, verbose = verbose)

        allocate(job_id_to_ik(ik_manager%n_jobs, 4))
        job_id_to_ik = 0

        call numerics%create_val_id_list(PW_dataset%n_val)
        call numerics%create_k_id_list(PW_dataset%n_k)

        call ik_manager%create_job_to_2d_ID_table(&
            numerics%val_id_list, &
            numerics%k_id_list, &
            job_id_to_ik, verbose = verbose)

        allocate(dielec(bins%n_E, bins%n_q, bins%n_q_theta, bins%n_q_phi))
        ! factor of 1 in dielectric formula
        dielec = (1.0_dp, 0.0_dp)

        call numerics%define_q_grid(bins%n_q*bins%q_width, &
            PW_dataset, FFT_grid)

        call numerics%compute_n_q_bin(bins, PW_dataset, verbose = verbose)

        ! allocate wave functions
        if ( PW_dataset%include_spin ) then

            allocate(wfc_iks(2, FFT_grid%n_grid(1), FFT_grid%n_grid(2), FFT_grid%n_grid(3)))
            allocate(wfc_fkfs(2, FFT_grid%n_grid(1), FFT_grid%n_grid(2), FFT_grid%n_grid(3)))

        else

            allocate(wfc_ik(FFT_grid%n_grid(1), FFT_grid%n_grid(2), FFT_grid%n_grid(3)))
            allocate(wfc_fkf(FFT_grid%n_grid(1), FFT_grid%n_grid(2), FFT_grid%n_grid(3)))

        end if

        if ( verbose ) then
            print*, 'Calculating the dielectric...'
            print*
        end if

        do j = 1, ik_manager%n_jobs_per_proc
            
            job_id = ik_manager%job_table(proc_id + 1, j)

            if ( job_id /= 0 ) then

                ! count states down from highest valence band
                val_id = job_id_to_ik(job_id, 1)
                k = job_id_to_ik(job_id, 2)

                if ( PW_dataset%include_spin ) then

                    ! load initial wave function
                    call PW_dataset%load_wfc_ik_expanded_spin(val_id, k, FFT_grid, wfc_iks)

                    do f = 1, numerics%n_cond_max

                        cond_id = PW_dataset%n_val + f

                        do kf = 1, PW_dataset%n_k

                            ! load final wave function
                            call PW_dataset%load_wfc_ik_expanded_spin(cond_id, kf, FFT_grid, wfc_fkfs)

                            ! compute dielectric
                            call dielectric_calc_vc(dielec, &
                                FFT_grid, PW_dataset, target_mat, bins, widths, numerics, &
                                wfc_iks, wfc_fkfs, val_id, cond_id, k, kf, verbose = .FALSE.)

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

                            ! compute dielectric
                            call dielectric_calc_vc(dielec, &
                                FFT_grid, PW_dataset, target_mat, bins, widths, numerics, &
                                wfc_ik, wfc_fkf, val_id, cond_id, k, kf, verbose = .FALSE.)

                        end do
                    end do

                end if

            end if

        end do

        if ( verbose ) then
            print*, 'Done calculating the dielectric!'
            print*
        end if

        call comm_reduce_dielectric(proc_id, root_process, dielec, verbose = verbose)

        if ( proc_id == root_process ) then

            call numerics%save(dielectric_output_filename, verbose = verbose)
            call bins%save(dielectric_output_filename, verbose = verbose)
            call save_dielectric(dielectric_output_filename, dielec, verbose = verbose)

        end if

    end subroutine

    subroutine save_dielectric(filename, dielectric, verbose)
        !! Saves the dielectric.

        implicit none

        character(len=*) :: filename
        complex(dp) :: dielectric(:, :, :, :)
        logical, optional :: verbose

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id

        logical :: file_exists

        integer(HSIZE_T) :: dims1(1) = [1]
        ! integer(HSIZE_T) :: dims2(2)
        integer(HSIZE_T) :: dims4(4)

        integer :: error

        if ( verbose ) then

            print*, 'Saving dielectric...'
            print*

        end if

        ! open the file to write
        call h5open_f(error)
        call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

        call h5gcreate_f(file_id, 'dielectric', group_id, error)

        dims4 = [size(dielectric, 1), &
                 size(dielectric, 2), &
                 size(dielectric, 3), &
                 size(dielectric, 4)]

        call h5ltmake_dataset_double_f(file_id, 'dielectric/dielectric_r', &
            size(dims4), dims4,&
            real(dielectric), error)

        call h5ltmake_dataset_double_f(file_id, 'dielectric/dielectric_c', &
            size(dims4), dims4,&
            aimag(dielectric), error)

        call h5fclose_f(file_id, error)
        call h5close_f(error)

    end subroutine

end module
