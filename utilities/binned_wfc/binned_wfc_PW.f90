program binned_wfc_PW
    !! Computes the binned wave function coefficients for data in a plane wave basis.

    use hdf5
    use h5lt

    use prec
    use timing
    use math_mod

    use PW_dataset_type
    use units

    implicit none

    real(dp) :: q_bin_width
    integer :: n_q_bins

    integer :: q, k, g, i, n, g1, g2, g3

    character(len=512) :: PW_dataset_filename = ''
    character(len=512) :: out_filename = './binned_wfcs.hdf5'

    character(len=512) :: q_bin_width_in
    character(len=512) :: n_q_bins_in

    real(dp) :: pc_vol_A
    real(dp) :: pc_vol

    real(dp) :: q_vec(3)
    real(dp) :: q_mag
    integer :: q_bin_num

    integer, allocatable :: n_q(:)

    real(dp), allocatable :: binned_wfc_FT_sq(:, :)

    complex(dp), allocatable :: wfc_ik_FT_no_spin(:)
    complex(dp), allocatable :: wfc_ik_FT_spin(:, :)

    logical :: verbose = .TRUE.

    integer :: n_grid(3)

    integer(HSIZE_T) :: dims2(2)
    integer(HID_T) :: file_id

    integer :: error

    type(PW_dataset_t) :: PW_dataset

    real(dp) :: time_start, time_end

    call get_command_argument(1, q_bin_width_in)
    call get_command_argument(2, n_q_bins_in)
    call get_command_argument(3, PW_dataset_filename)
    call get_command_argument(4, out_filename)

    read(q_bin_width_in, *) q_bin_width
    read(n_q_bins_in, *) n_q_bins

    print*
    print*, '- binned wfc -'
    print*

    call cpu_time(time_start) 

    call PW_dataset%load(trim(PW_dataset_filename), verbose = verbose)

    allocate(binned_wfc_FT_sq(PW_dataset%n_cond + PW_dataset%n_val, n_q_bins))
    allocate(n_q(n_q_bins))
    binned_wfc_FT_sq = 0.0_dp
    n_q = 0

    print*, 'Calculating...'

    if ( .not. PW_dataset%include_spin ) then

        allocate(wfc_ik_FT_no_spin(PW_dataset%n_G))
        wfc_ik_FT_no_spin = (0.0_dp, 0.0_dp)

        do i = 1, PW_dataset%n_val + PW_dataset%n_cond
            do k = 1, PW_dataset%n_k

                call PW_dataset%load_wfc_FT_ik_no_spin(i, k, wfc_ik_FT_no_spin)

                do g = 1, PW_dataset%n_G

                    q_vec = PW_dataset%k_grid_xyz(k, :) + PW_dataset%G_grid_xyz(g, :)
                    q_mag = norm2(q_vec)

                    q_bin_num = Q_func(q_mag, 0.0_dp, q_bin_width, n_q_bins)
                    
                    if ( i == 1 ) then 
                        n_q(q_bin_num) = n_q(q_bin_num) + 1
                    end if

                    binned_wfc_FT_sq(i, q_bin_num) = binned_wfc_FT_sq(i, q_bin_num) + abs(wfc_ik_FT_no_spin(g))**2

                end do

            end do
        end do

    else 

        allocate(wfc_ik_FT_spin(PW_dataset%n_G, 2))
        wfc_ik_FT_spin = (0.0_dp, 0.0_dp)

        do i = 1, PW_dataset%n_val + PW_dataset%n_cond
            do k = 1, PW_dataset%n_k

                call PW_dataset%load_wfc_FT_ik_spin(i, k, wfc_ik_FT_spin)

                do g = 1, PW_dataset%n_G

                    q_vec = PW_dataset%k_grid_xyz(k, :) + PW_dataset%G_grid_xyz(g, :)
                    q_mag = norm2(q_vec)

                    q_bin_num = Q_func(q_mag, 0.0_dp, q_bin_width, n_q_bins)
                    
                    if ( i == 1 ) then 
                        n_q(q_bin_num) = n_q(q_bin_num) + 1
                    end if

                    binned_wfc_FT_sq(i, q_bin_num) = binned_wfc_FT_sq(i, q_bin_num) + &
                        conjg(wfc_ik_FT_spin(g, 1))*wfc_ik_FT_spin(g, 1) + &
                        conjg(wfc_ik_FT_spin(g, 2))*wfc_ik_FT_spin(g, 2)

                end do

            end do
        end do

    end if

    do q = 1, n_q_bins
        
        binned_wfc_FT_sq(:, q) = binned_wfc_FT_sq(:, q)/max(1.0_dp, 1.0_dp*n_q(q))

    end do

    ! save data

    print*, 'Done computing. Saving data...'

    call h5open_f(error)
    call h5fcreate_f(out_filename, H5F_ACC_TRUNC_F, file_id, error)

    dims2 = [size(binned_wfc_FT_sq, 1), size(binned_wfc_FT_sq, 2)]
    call h5ltmake_dataset_double_f(file_id, 'binned_wfc_FT_sq', size(dims2), dims2,&
        binned_wfc_FT_sq, error)

    call h5fclose_f(file_id, error)
    call h5close_f(error)

    call PW_dataset%save(out_filename, verbose = verbose)

    call cpu_time(time_end)

    print*, 'Run time : ', trim(pretty_time_format(time_end - time_start))
    print*

end program
