program binned_wfc
    !! Computes the binned wave function coefficients
    use hdf5
    use h5lt

    use prec
    use timing
    use units
    use DFT_parameters
    use core_electrons

    use FFT_util

    implicit none

    real(dp) :: q_bin_width
    integer :: n_q_bins

    integer :: q, k, g, i, n, g1, g2, g3

    character(len=512) :: DFT_input_filename = ''
    character(len=512) :: core_elec_config_filename = ''
    character(len=512) :: sto_wf_filename = ''
    character(len=512) :: out_filename = './binned_wfcs.hdf5'

    character(len=512) :: q_bin_width_in
    character(len=512) :: n_q_bins_in
    character(len=512) :: pc_vol_A_in

    real(dp) :: pc_vol_A
    real(dp) :: pc_vol

    real(dp) :: q_vec(3)
    real(dp) :: q_mag
    integer :: q_bin_num

    integer, allocatable :: n_q(:)

    real(dp), allocatable :: dft_u_sq(:, :)
    real(dp), allocatable :: core_u_sq(:, :)

    complex(dp), allocatable :: wfc_FT(:, :)

    logical :: verbose = .TRUE.

    integer :: n_grid(3)

    integer :: core_wfc_fft_plan(8)

    complex(dp), allocatable :: core_wfc(:, :, :)
    complex(dp), allocatable :: core_wfc_FT(:, :, :)

    integer(HSIZE_T) :: dims2(2)
    integer(HID_T) :: file_id

    integer :: error

    integer :: n_grid_tot

    real(dp) :: time_start, time_end

    call get_command_argument(1, q_bin_width_in)
    call get_command_argument(2, n_q_bins_in)
    call get_command_argument(3, pc_vol_A_in)
    call get_command_argument(4, DFT_input_filename)
    call get_command_argument(5, core_elec_config_filename)
    call get_command_argument(6, sto_wf_filename)
    call get_command_argument(7, out_filename)

    read(q_bin_width_in, *) q_bin_width
    read(n_q_bins_in, *) n_q_bins
    read(pc_vol_A_in, *) pc_vol_A

    pc_vol = Ang_to_inv_eV**3*pc_vol_A

    print*
    print*, '- binned wfc -'
    print*

    call cpu_time(time_start) 

    ! load data...
    call load_DFT_parameters(trim(DFT_input_filename), verbose=verbose)
    call load_core_elec_config(trim(core_elec_config_filename), verbose=verbose)
    call load_core_sto_data(trim(sto_wf_filename), verbose=verbose)

    allocate(dft_u_sq(n_bands, n_q_bins + 1))
    allocate(core_u_sq(n_core_states, n_q_bins + 1))
    allocate(wfc_FT(n_k, n_in_G))
    allocate(n_q(n_q_bins + 1))

    dft_u_sq = 0.0_dp
    core_u_sq = 0.0_dp
    wfc_FT = (0.0_dp, 0.0_dp)
    n_q = 0

    ! compute on the same q grid

    print*, 'Running dft piece...'

    ! dft
    do i = 1, n_bands

        ! get the wave function coefficients
        call get_in_wfc_FT(DFT_input_filename, i, wfc_FT)
        
        do g = 1, n_in_G 
            do k = 1, n_k

                q_vec = k_grid_xyz(k, :) + in_G_grid_xyz(g, :)
                q_mag = norm2(q_vec)

                q_bin_num = min(1 + floor(q_mag/q_bin_width), n_q_bins + 1)
                
                if ( i == 1 ) then 
                    n_q(q_bin_num) = n_q(q_bin_num) + 1
                end if

                dft_u_sq(i, q_bin_num) = dft_u_sq(i, q_bin_num) + abs(wfc_FT(k, g))**2

            end do

        end do

    end do

    do q = 1, n_q_bins + 1
        
        dft_u_sq(:, q) = dft_u_sq(:, q)/max(1.0_dp, 1.0_dp*n_q(q))

    end do

    !!!

    n_q = 0

    n_grid = [10, 10, 10]

    n_grid_tot = n_grid(1)*n_grid(2)*n_grid(3)

    allocate(core_wfc_FT(n_grid(1), n_grid(2), n_grid(3)))
    allocate(core_wfc(n_grid(1), n_grid(2), n_grid(3)))

    core_wfc = (0.0_dp, 0.0_dp)
    core_wfc_FT = (0.0_dp, 0.0_dp)

    ! plan the fft
    call set_fft_plan_forward_3d(n_grid, core_wfc_fft_plan)
    call set_sym_FFT_G_grid_xyz(n_grid, k_red_to_xyz, verbose)

    print*, 'Running core piece...'

    ! core
    do n = 1, n_core_states

        call calc_core_sto_wf_grid(core_wfc, n, n_grid, red_to_xyz, shift=.TRUE.)

        do k = 1, n_k

            ! compute fft
            call dfftw_execute_dft(core_wfc_fft_plan, (1.0_dp/n_grid_tot)*sqrt(pc_vol)*core_wfc, core_wfc_FT) 

            do g3 = 1, n_grid(3)
                do g2 = 1, n_grid(2)
                    do g1 = 1, n_grid(1)

                        q_vec = k_grid_xyz(k, :) + sym_FFT_G_grid_xyz(g1, g2, g3, :)
                        q_mag = norm2(q_vec)

                        q_bin_num = min(1 + floor(q_mag/q_bin_width), n_q_bins + 1)
                        
                        if ( n == 1 ) then
                            n_q(q_bin_num) = n_q(q_bin_num) + 1
                        end if

                        core_u_sq(n, q_bin_num) = core_u_sq(n, q_bin_num) + abs(core_wfc_FT(g1, g2, g3))**2

                    end do
                end do
            end do

        end do
    end do

    do q = 1, n_q_bins + 1
        
        core_u_sq(:, q) = core_u_sq(:, q)/max(1.0_dp, 1.0_dp*n_q(q))

    end do

    ! save data

    print*, 'Done computing. Saving data...'

    call h5open_f(error)
    call h5fcreate_f(out_filename, H5F_ACC_TRUNC_F, file_id, error)

    dims2 = [n_bands, n_q_bins + 1]
    call h5ltmake_dataset_double_f(file_id, 'dft_binned_wfc', size(dims2), dims2,&
        dft_u_sq, error)
    dims2 = [n_core_states, n_q_bins + 1]
    call h5ltmake_dataset_double_f(file_id, 'core_binned_wfc', size(dims2), dims2,&
        core_u_sq, error)

    call h5fclose_f(file_id, error)
    call h5close_f(error)

    call save_DFT_parameters(out_filename, verbose = verbose)
    call save_core_electrons(out_filename, verbose = verbose)

    call cpu_time(time_end)

    print*, 'Run time : ', trim(pretty_time_format(time_end - time_start))
    print*

    print*, '----------'
    print*

end program
