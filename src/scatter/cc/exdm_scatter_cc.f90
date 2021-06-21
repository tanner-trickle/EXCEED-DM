module exdm_scatter_cc
    !! Compute the scattering rate from core to conduction states

    use prec

    use control_input
    use numerics_input
    use material_input
    use DFT_parameters
    use core_electrons
    use FFT_util
    use transition

    use calc_exdm_scatter_cc

    implicit none

contains

    subroutine run_dme_scatter_cc(binned_rate_t, n_tran_per_proc, DFT_input_filename, &
            sto_wf_filename, core_elec_config_filename, out_filename,&
            proc_id, root_process, verbose)

        implicit none

        real(dp) :: binned_rate_t(:, :, :, :, :, :)

        logical, optional :: verbose

        integer :: proc_id, root_process

        integer :: n_tran_per_proc
        integer :: t

        integer :: tran_id
        integer :: init_id, cond_id

        character(len=*) :: DFT_input_filename
        character(len=*) :: out_filename
        character(len=*) :: sto_wf_filename
        character(len=*) :: core_elec_config_filename

        integer :: n_FFT_grid(3)
        integer :: n_FFT

        integer :: wfc_fft_plan(8)
        integer :: Tif_fft_plan(8)

        complex(dp), allocatable :: wfc_i(:, :, :)
        complex(dp), allocatable :: wfc_FT_f(:, :)

        if ( verbose ) then

            print*, 'Starting c -> c scattering rate calculation...'
            print*

        end if

        ! calculation setup
        call load_DFT_parameters(trim(DFT_input_filename), verbose=verbose)
        call do_scissor_correction(band_gap, verbose = verbose)

        call load_core_elec_config(trim(core_elec_config_filename), verbose=verbose)
        call load_core_sto_data(trim(sto_wf_filename), verbose=verbose)

        ! doubled to avoid wrapping problems
        n_FFT_grid(1) = max(&
            2*(maxval(in_G_grid_red(:, 1)) - minval(in_G_grid_red(:, 1))) + 1, &
            n_FFT_grid_input(1))
        n_FFT_grid(2) = max(&
            2*(maxval(in_G_grid_red(:, 2)) - minval(in_G_grid_red(:, 2))) + 1, &
            n_FFT_grid_input(2))
        n_FFT_grid(3) = max(&
            2*(maxval(in_G_grid_red(:, 3)) - minval(in_G_grid_red(:, 3))) + 1, &
            n_FFT_grid_input(3))

        n_FFT = n_FFT_grid(1)*n_FFT_grid(2)*n_FFT_grid(3)

        if ( verbose ) then

            print*, '----------------------------------------'
            print*, '    ---'
            print*, '    FFT'
            print*, '    ---'
            print*
            print*, '        Dimensions of FFT grid : '
            print*, '        ', n_FFT_grid
            print*
            print*, '        Number of G points in FFT grid = ', n_FFT
            print*
            print*, '----------------------------------------'
            print*

        end if

        call set_fft_plan_backward_3d(n_FFT_grid, wfc_fft_plan)
        call set_fft_plan_backward_3d(n_FFT_grid, Tif_fft_plan)

        call set_sym_FFT_G_grid_xyz(n_FFT_grid, k_red_to_xyz, verbose = verbose)

        ! time calculation

        if ( ( proc_id == root_process ) .and. ( timer ) ) then

            call time_exdm_scatter_cc_calc(DFT_input_filename, 1, wfc_fft_plan, &
                Tif_fft_plan, n_FFT_grid, verbose = verbose)

        end if

        allocate(wfc_FT_f(n_k, n_in_G))
        allocate(wfc_i(n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3)))
        wfc_i = (0.0_dp, 0.0_dp)

        ! do calculation

        if ( verbose ) then

            print*, 'Calculating transition rates...'
            print*

        end if

        do t = 1, n_tran_per_proc

            tran_id = job_table(proc_id + 1, t)

            if ( tran_id .ne. 0 ) then

                init_id = tran_to_init_fin_id(tran_id, 1)
                cond_id = tran_to_init_fin_id(tran_id, 2) + n_val

                call calc_core_sto_wf_grid(wfc_i, init_id, n_FFT_grid, red_to_xyz, &
                    shift = .TRUE., verbose = verbose)

                wfc_i = sqrt(pc_vol)*wfc_i

                call get_in_wfc_FT(DFT_input_filename, cond_id, wfc_FT_f)

                call dme_scatter_cc_calc(binned_rate_t(:, :, :, :, :, t),& 
                    wfc_i, wfc_FT_f, init_id, cond_id, n_FFT_grid, n_k, &
                    wfc_fft_plan, Tif_FFT_plan, verbose = verbose)

            end if

        end do

        if ( verbose ) then
            print*, 'Done calculating transition rates!'
            print*
        end if

        if ( proc_id == root_process ) then

            call save_DFT_parameters(out_filename, verbose = verbose)
            call save_core_electrons(out_filename, verbose = verbose)

        end if

    end subroutine

    subroutine time_exdm_scatter_cc_calc(DFT_input_filename, tran_id, wfc_fft_plan, &
            Tif_fft_plan, n_FFT_grid, verbose)
        !! Times the c -> c scattering rate calculation
        use timing 
        use mpi

        implicit none

        character(len=*) :: DFT_input_filename

        integer :: wfc_fft_plan(8)
        integer :: Tif_fft_plan(8)
        integer :: tran_id

        integer :: n_FFT_grid(3)

        logical, optional :: verbose
        real(dp) :: b_rate(n_q_bins + 1, n_E_bins + 1, n_mX, n_FDM, n_time)

        complex(dp), allocatable :: wfc_i(:, :, :)
        complex(dp), allocatable :: wfc_FT_f(:, :)

        integer :: init_id, cond_id

        if ( verbose ) then

            print*, 'Timing c -> c calculation...'
            print*

        end if

        allocate(wfc_FT_f(n_k, n_in_G))
        allocate(wfc_i(n_FFT_grid(1), n_FFT_grid(2), n_FFT_grid(3)))

        init_id = tran_to_init_fin_id(tran_id, 1)
        cond_id = tran_to_init_fin_id(tran_id, 2) + n_val

        time(3) = MPI_Wtime()

        call calc_core_sto_wf_grid(wfc_i, init_id, n_FFT_grid, red_to_xyz, &
            shift = .TRUE., verbose = verbose)
        wfc_i = sqrt(pc_vol)*wfc_i

        call get_in_wfc_FT(DFT_input_filename, cond_id, wfc_FT_f)

        time(4) = MPI_Wtime()

        call dme_scatter_cc_calc(b_rate,& 
            wfc_i, wfc_FT_f, init_id, cond_id, n_FFT_grid, 1, &
            wfc_fft_plan, Tif_FFT_plan, verbose = verbose)

        time(5) = MPI_Wtime()

        if ( verbose ) then

            print*, '----------------------------------------'
            print*, '    -------------'
            print*, '    Timing (TEST)'
            print*, '    -------------'
            print*
            print*, '        (TEST) Core WF :'
            print*, '            ', trim(pretty_time_format(time(4) - time(3)))
            print*
            print*, '        (Total) Core WF :'
            print*, '            ', trim(pretty_time_format(&
                n_tran_per_proc*(time(4) - time(3))&
                ))
            print*
            print*, '        (TEST) Rate : '
            print*, '            ', trim(pretty_time_format(time(5) - time(4)))
            print*
            print*, '        (Total) Rate : '
            print*, '            ', trim(pretty_time_format(&
                n_tran_per_proc*n_k*(time(5) - time(4))&
                ))
            print*
            print*, '        Expected run time for whole calculation :'
            print*, '            ', trim(pretty_time_format(&
                n_tran_per_proc*(time(4) - time(3)) + n_tran_per_proc*n_k*(time(5) - time(4))&
                ))
            print*
            print*, '----------------------------------------'
            print*

        end if

    end subroutine

end module
