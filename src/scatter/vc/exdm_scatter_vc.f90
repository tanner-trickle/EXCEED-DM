module exdm_scatter_vc
    !! Compute the scattering rate from valence to conduction states

    use prec

    use numerics_input
    use control_input
    use material_input
    use DFT_parameters

    use FFT_util
    use transition

    use calc_exdm_scatter_vc

    implicit none

contains

    subroutine run_dme_scatter_vc(binned_rate_t, n_tran_per_proc, DFT_input_filename, &
            out_filename, proc_id, root_process, verbose)
        implicit none

        real(dp) :: binned_rate_t(:, :, :, :, :, :)

        logical, optional :: verbose

        integer :: proc_id, root_process

        integer :: n_tran_per_proc
        integer :: t

        integer :: tran_id
        integer :: val_id, cond_id

        character(len=*) :: DFT_input_filename
        character(len=*) :: out_filename

        integer :: n_FFT_grid(3)
        integer :: n_FFT

        integer :: wfc_fft_plan(8)
        integer :: Tif_fft_plan(8)

        complex(dp), allocatable :: wfc_FT_i(:, :)
        complex(dp), allocatable :: wfc_FT_f(:, :)

        if ( verbose ) then

            print*, 'Starting v -> c scattering rate calculation...'
            print*

        end if

        ! calculation setup
        call load_DFT_parameters(DFT_input_filename, verbose = verbose)
        call do_scissor_correction(band_gap, verbose = verbose)

        ! doubled to avoid wrapping problems
        n_FFT_grid(1) = 2*(maxval(in_G_grid_red(:, 1)) - minval(in_G_grid_red(:, 1))) + 1
        n_FFT_grid(2) = 2*(maxval(in_G_grid_red(:, 2)) - minval(in_G_grid_red(:, 2))) + 1
        n_FFT_grid(3) = 2*(maxval(in_G_grid_red(:, 3)) - minval(in_G_grid_red(:, 3))) + 1

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

            call time_exdm_scatter_vc_calc(DFT_input_filename, 1, wfc_fft_plan, &
                Tif_fft_plan, n_FFT_grid, verbose = verbose)

        end if

        allocate(wfc_FT_i(n_k, n_in_G))
        allocate(wfc_FT_f(n_k, n_in_G))

        ! do calculation

        if ( verbose ) then

            print*, 'Calculating transition rates...'
            print*

        end if

        do t = 1, n_tran_per_proc

            tran_id = job_table(proc_id + 1, t)

            if ( tran_id .ne. 0 ) then

                val_id = tran_to_init_fin_id(tran_id, 1)
                cond_id = tran_to_init_fin_id(tran_id, 2) + n_val

                call get_in_wfc_FT(DFT_input_filename, val_id, wfc_FT_i)
                call get_in_wfc_FT(DFT_input_filename, cond_id, wfc_FT_f)

                call dme_scatter_vc_calc(binned_rate_t(:, :, :, :, :, t),& 
                    wfc_FT_i, wfc_FT_f, val_id, cond_id, n_FFT_grid, n_k, &
                    wfc_fft_plan, Tif_FFT_plan, verbose = verbose)
            end if

        end do

        if ( verbose ) then

            print*, 'Done calculating transition rates!'
            print*

        end if

        if ( proc_id == root_process ) then

            call save_DFT_parameters(out_filename, verbose = verbose)

        end if

    end subroutine

    subroutine time_exdm_scatter_vc_calc(DFT_input_filename, tran_id, wfc_fft_plan, &
            Tif_fft_plan, n_FFT_grid, verbose)
        !! Times the v -> c scattering rate calculation
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

        complex(dp), allocatable :: wfc_FT_i(:, :)
        complex(dp), allocatable :: wfc_FT_f(:, :)

        integer :: val_id, cond_id

        if ( verbose ) then

            print*, 'Timing v -> c calculation...'
            print*

        end if

        allocate(wfc_FT_i(n_k, n_in_G))
        allocate(wfc_FT_f(n_k, n_in_G))

        val_id = tran_to_init_fin_id(tran_id, 1)
        cond_id = tran_to_init_fin_id(tran_id, 2) + n_val

        call get_in_wfc_FT(DFT_input_filename, val_id, wfc_FT_i)
        call get_in_wfc_FT(DFT_input_filename, cond_id, wfc_FT_f)

        time(3) = MPI_Wtime()

        call dme_scatter_vc_calc(b_rate,& 
            wfc_FT_i, wfc_FT_f, val_id, cond_id, n_FFT_grid, 1, &
            wfc_fft_plan, Tif_FFT_plan, verbose = verbose)

        time(4) = MPI_Wtime()

        if ( verbose ) then

            print*, '----------------------------------------'
            print*, '    -------------'
            print*, '    Timing (TEST)'
            print*, '    -------------'
            print*
            print*, '        (TEST) Run time: '
            print*, '            ', trim(pretty_time_format(time(4) - time(3)))
            print*
            print*, '        Expected run time for whole calculation : '
            print*, '            ', trim(pretty_time_format(&
                n_tran_per_proc*n_k**2*(time(4) - time(3))&
                ))
            print*
            print*, '----------------------------------------'
            print*

        end if

    end subroutine

end module
