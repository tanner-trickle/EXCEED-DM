module dme_scatter_vc
    !! Compute the scattering rate from valence to conduction states

    use prec

    use numerics_input
    use material_input
    use DFT_parameters

    use FFT_util
    use transition

    use calc_dme_scatter_vc

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

        ! calculation setup
        call load_DFT_parameters(DFT_input_filename, verbose = verbose)
        call do_scissor_correction(band_gap, verbose = verbose)

        ! doubled to avoid wrapping problems
        n_FFT_grid(1) = 2*(maxval(in_G_grid_red(:, 1)) - minval(in_G_grid_red(:, 1))) + 1
        n_FFT_grid(2) = 2*(maxval(in_G_grid_red(:, 2)) - minval(in_G_grid_red(:, 2))) + 1
        n_FFT_grid(3) = 2*(maxval(in_G_grid_red(:, 3)) - minval(in_G_grid_red(:, 3))) + 1

        n_FFT = n_FFT_grid(1)*n_FFT_grid(2)*n_FFT_grid(3)

        if ( verbose ) then

            print*, 'Dimensions of FFT grid : '
            print*, n_FFT_grid
            print*
            print*, 'Number of G points in FFT grid = ', n_FFT
            print*
            print*, '----------'
            print*

        end if

        call set_fft_plan_backward_3d(n_FFT_grid, wfc_fft_plan)
        call set_fft_plan_backward_3d(n_FFT_grid, Tif_fft_plan)

        call set_sym_FFT_G_grid_xyz(n_FFT_grid, k_red_to_xyz, verbose = verbose)

        ! time calculation

        allocate(wfc_FT_i(n_k, n_in_G))
        allocate(wfc_FT_f(n_k, n_in_G))

        ! do calculation

        if ( verbose ) then

            print*, 'Calculating rate...'
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
            print*, '----------'
            print*
        end if

        if ( proc_id == root_process ) then

            call save_DFT_parameters(out_filename, verbose = verbose)

        end if

    end subroutine

end module
