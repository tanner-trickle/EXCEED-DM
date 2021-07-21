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

        complex(dp), allocatable :: wfc_FT_i_s(:, :, :)
        complex(dp), allocatable :: wfc_FT_f_s(:, :, :)

        real(dp) :: q_max_FFT

        if ( verbose ) then

            print*, 'Starting v -> c scattering rate calculation...'
            print*

        end if

        ! calculation setup
        !! only load these variables if they haven't been loaded before
        !! to compute the dielectric
        if ( trim(screen_type) /= 'numeric' ) then
            call load_DFT_parameters(DFT_input_filename, verbose = verbose)
            call do_scissor_correction(band_gap, verbose = verbose)
        else 
            if ( load_dielectric_from_file ) then
                call load_DFT_parameters(DFT_input_filename, verbose = verbose)
                call do_scissor_correction(band_gap, verbose = verbose)
            end if
        end if

        ! doubled to avoid wrapping problems
        n_FFT_grid(1) = 2*(maxval(in_G_grid_red(:, 1)) - minval(in_G_grid_red(:, 1))) + 1
        n_FFT_grid(2) = 2*(maxval(in_G_grid_red(:, 2)) - minval(in_G_grid_red(:, 2))) + 1
        n_FFT_grid(3) = 2*(maxval(in_G_grid_red(:, 3)) - minval(in_G_grid_red(:, 3))) + 1

        n_FFT = n_FFT_grid(1)*n_FFT_grid(2)*n_FFT_grid(3)

        call find_q_max_FFT(q_max_FFT, k_red_to_xyz, n_FFT_grid, verbose = .FALSE.)

        ! call calc_eig_system_33((1.0_dp, 0.0_dp)*k_red_to_xyz, k_red_eig_vals, k_red_eig_vecs)
        ! q_s_FFT = minval(abs(k_red_eig_vals))/2.0_dp

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
            print*, '        q_max_FFT = ', q_max_FFT/1.0e3_dp, ' keV'
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
                Tif_fft_plan, n_FFT_grid, q_max_FFT, verbose = verbose)

        end if

        ! do calculation

        if ( verbose ) then

            print*, 'Calculating transition rates...'
            print*

        end if

        if ( include_spin ) then

            allocate(wfc_FT_i_s(n_k, n_in_G, 2))
            allocate(wfc_FT_f_s(n_k, n_in_G, 2))

        else
            
            allocate(wfc_FT_i(n_k, n_in_G))
            allocate(wfc_FT_f(n_k, n_in_G))

        end if

        do t = 1, n_tran_per_proc

            tran_id = job_table(proc_id + 1, t)

            if ( tran_id .ne. 0 ) then

                val_id = tran_to_init_fin_id(tran_id, 1)
                cond_id = tran_to_init_fin_id(tran_id, 2) + n_val

                if ( include_spin ) then

                    call get_in_wfc_FT(DFT_input_filename, val_id, wfc_FT_i_s)
                    call get_in_wfc_FT(DFT_input_filename, cond_id, wfc_FT_f_s)

                    call dme_scatter_vc_calc(binned_rate_t(:, :, :, :, :, t),& 
                        wfc_FT_i_s, wfc_FT_f_s, val_id, cond_id, n_FFT_grid, n_k, &
                        wfc_fft_plan, Tif_FFT_plan, q_max_FFT, verbose = verbose)

                else

                    call get_in_wfc_FT(DFT_input_filename, val_id, wfc_FT_i)
                    call get_in_wfc_FT(DFT_input_filename, cond_id, wfc_FT_f)

                    call dme_scatter_vc_calc(binned_rate_t(:, :, :, :, :, t),& 
                        wfc_FT_i, wfc_FT_f, val_id, cond_id, n_FFT_grid, n_k, &
                        wfc_fft_plan, Tif_FFT_plan, q_max_FFT, verbose = verbose)

                end if

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
            Tif_fft_plan, n_FFT_grid, q_max_FFT, verbose)
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

        complex(dp), allocatable :: wfc_FT_i_s(:, :, :)
        complex(dp), allocatable :: wfc_FT_f_s(:, :, :)

        integer :: val_id, cond_id

        real(dp) :: q_max_FFT

        if ( verbose ) then

            print*, 'Timing v -> c calculation...'
            print*

        end if

        if ( include_spin ) then

            allocate(wfc_FT_i_s(n_k, n_in_G, 2))
            allocate(wfc_FT_f_s(n_k, n_in_G, 2))

        else
            
            allocate(wfc_FT_i(n_k, n_in_G))
            allocate(wfc_FT_f(n_k, n_in_G))

        end if

        val_id = tran_to_init_fin_id(tran_id, 1)
        cond_id = tran_to_init_fin_id(tran_id, 2) + n_val

        if ( include_spin ) then

            call get_in_wfc_FT(DFT_input_filename, val_id, wfc_FT_i_s)
            call get_in_wfc_FT(DFT_input_filename, cond_id, wfc_FT_f_s)

            time(3) = MPI_Wtime()
            
            call dme_scatter_vc_calc(b_rate,& 
                wfc_FT_i_s, wfc_FT_f_s, val_id, cond_id, n_FFT_grid, 1, &
                wfc_fft_plan, Tif_FFT_plan, q_max_FFT, verbose = verbose)

            time(4) = MPI_Wtime()

        else

            call get_in_wfc_FT(DFT_input_filename, val_id, wfc_FT_i)
            call get_in_wfc_FT(DFT_input_filename, cond_id, wfc_FT_f)

            time(3) = MPI_Wtime()

            call dme_scatter_vc_calc(b_rate,& 
                wfc_FT_i, wfc_FT_f, val_id, cond_id, n_FFT_grid, 1, &
                wfc_fft_plan, Tif_FFT_plan, q_max_FFT, verbose = verbose)

            time(4) = MPI_Wtime()

        end if

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
