module exdm_scatter
    !! Compute the DM-electron scattering rate.

    use mpi
    use hdf5
    use h5lt

    use prec

    use control_input
    use io_input
    use MPI_util

    use material_type
    use dm_model_type
    use expt_type
    use binned_scatter_rate_type
    use bins_scatter_type
    use in_med_scr_type

    use exdm_scatter_vc
    use exdm_scatter_cc
    use exdm_scatter_vf
    use exdm_scatter_cf

    implicit none

contains

    subroutine run_exdm_scatter(proc_id, root_process, n_proc, &
            io_files, main_control, target_mat, expt, dm_model, &
            verbose)
        !! Compute the DM-electron scattering rate.
        implicit none

        integer          :: proc_id
        integer          :: root_process
        integer          :: n_proc
        type(io_files_t) :: io_files
        type(control_t)  :: main_control
        type(material_t) :: target_mat
        type(dm_model_t) :: dm_model
        type(expt_t)     :: expt
        type(bins_scatter_t)     :: bins

        type(in_med_scr_t) :: in_med_scr

        logical, optional :: verbose

        integer :: status(MPI_STATUS_SIZE)
        integer :: tag = 0
        integer :: err
        integer :: ierr

        integer :: t, i

        type(binned_scatter_rate_t), allocatable :: binned_rate_init(:)
            !! Dim : [ n_init ]
            !!
            !! Binned rate per cross section, per kg-year for each initial state
            !!
            !! Units : cm^(-2)
        integer :: n_init
            !! Number of initial states

        character(len=64) :: calc_modes_list(4) = ['vc', 'vf', 'cc', 'cf']

        if ( verbose ) then
            print*, 'Starting scattering rate calculation...'
            print*
        end if

        ! load the screening parameters
        call in_med_scr%load(proc_id, root_process, n_proc, io_files, main_control, &
            target_mat, verbose = verbose)

        ! load the binning parameters
        call bins%load(io_files%nml_input_filename, verbose = verbose)

        ! set the n_init parameter
        if ( trim(main_control%calc_mode) == 'vc' ) then
            call exdm_scatter_vc_set_n_init(io_files, n_init)
        else if ( trim(main_control%calc_mode) == 'cc' ) then
            call exdm_scatter_cc_set_n_init(io_files, n_init)
        else if ( trim(main_control%calc_mode) == 'vf' ) then
            call exdm_scatter_vf_set_n_init(io_files, n_init)
        else if ( trim(main_control%calc_mode) == 'cf' ) then
            call exdm_scatter_cf_set_n_init(io_files, n_init)
        end if

        allocate(binned_rate_init(n_init))
        do i = 1, n_init
            call binned_rate_init(i)%init(bins, dm_model, expt)
        end do

        if ( trim(main_control%calc_mode) == 'vc' ) then

            call run_exdm_scatter_vc(n_init, n_proc, proc_id, root_process, &
                binned_rate_init, main_control, io_files, target_mat, dm_model, &
                bins, expt, in_med_scr, verbose = verbose)

        else if ( trim(main_control%calc_mode) == 'cc' ) then

            call run_exdm_scatter_cc(n_init, n_proc, proc_id, root_process, &
                binned_rate_init, main_control, io_files, target_mat, dm_model, &
                bins, expt, in_med_scr, verbose = verbose)

        else if ( trim(main_control%calc_mode) == 'vf' ) then

            call run_exdm_scatter_vf(n_init, n_proc, proc_id, root_process, &
                binned_rate_init, main_control, io_files, target_mat, dm_model, &
                bins, expt, in_med_scr, verbose = verbose)

        else if ( trim(main_control%calc_mode) == 'cf' ) then

            call run_exdm_scatter_cf(n_init, n_proc, proc_id, root_process, &
                binned_rate_init, main_control, io_files, target_mat, dm_model, &
                bins, expt, in_med_scr, verbose = verbose)

        else

            call print_error_message('Calculation mode : '//trim(main_control%calc_mode)//' is not a valid option.', &
                verbose = verbose)
            stop

        end if

        ! communicate computed data
        call comm_reduce_binned_rate_init(proc_id, root_process, binned_rate_init, &
            verbose)

        ! save data
        if ( proc_id == root_process ) then

            call bins%save(io_files%out_filename, verbose = verbose)
            call in_med_scr%save(io_files%out_filename, verbose = verbose)
            call save_scatter_rates(io_files%out_filename, &
                bins, expt, dm_model, binned_rate_init, verbose = verbose)

        end if

    end subroutine

    subroutine save_scatter_rates(filename, bins, expt, dm_model, &
            binned_rate_init, verbose)
        !! Saves the binned and total scattering rates.

        implicit none

        character(len=*) :: filename
        type(bins_scatter_t) :: bins
        type(expt_t) :: expt
        type(dm_model_t) :: dm_model
        type(binned_scatter_rate_t) :: binned_rate_init(:)

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id

        integer :: n_init

        type(binned_scatter_rate_t) :: total_binned_rate
            !! Total binned scattering rate.
            !!
            !! Units : cm^(-2)

        real(dp) :: rate(dm_model%n_mX, dm_model%n_med_FF, expt%n_time)

        logical, optional :: verbose
        logical :: file_exists

        integer(HSIZE_T) :: dims1(1) = [1]
        integer(HSIZE_T) :: dims2(2)

        integer :: error

        integer :: m, f, t
        integer :: i, fin

        n_init = size(binned_rate_init)

        call total_binned_rate%init(bins, dm_model, expt)
        ! compute total binned rate
        do i = 1, n_init
            total_binned_rate%binned_rate = total_binned_rate%binned_rate + &
                binned_rate_init(i)%binned_rate
        end do

        rate = total_binned_rate%compute_rate()

        if ( verbose ) then
            print*, 'Saving scattering rates...'
            print*
        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'scatter_rates', group_id, error)

            dims2 = [bins%n_q, bins%n_E]

            do t = 1, expt%n_time

                call h5gcreate_f(file_id,&
                    'scatter_rates'//&
                    '/t_'//trim(adjustl(int_to_str(t))),&
                    group_id, error)

                do f = 1, dm_model%n_med_FF 

                    call h5gcreate_f(file_id,&
                        'scatter_rates'//&
                        '/t_'//trim(adjustl(int_to_str(t)))//&
                        '/f_'//trim(adjustl(int_to_str(f))),&
                        group_id, error)

                    do m = 1, dm_model%n_mX

                        call h5gcreate_f(file_id,& 
                            'scatter_rates'//&
                            '/t_'//trim(adjustl(int_to_str(t)))//&
                            '/f_'//trim(adjustl(int_to_str(f)))//&
                            '/m_'//trim(adjustl(int_to_str(m))),&
                            group_id, error)

                        call h5ltmake_dataset_double_f(file_id,&
                            'scatter_rates'//&
                            '/t_'//trim(adjustl(int_to_str(t)))//&
                            '/f_'//trim(adjustl(int_to_str(f)))//&
                            '/m_'//trim(adjustl(int_to_str(m)))//&
                            '/total', &
                            size(dims1), dims1,&
                            rate(m, f, t), error)

                        call h5ltmake_dataset_double_f(file_id,&
                            'scatter_rates'//&
                            '/t_'//trim(adjustl(int_to_str(t)))//&
                            '/f_'//trim(adjustl(int_to_str(f)))//&
                            '/m_'//trim(adjustl(int_to_str(m)))//&
                            '/total_binned', &
                            size(dims2), dims2,&
                            total_binned_rate%binned_rate(:, :, m, f, t), error)

                        do i = 1, n_init

                            call h5gcreate_f(file_id,& 
                                'scatter_rates'//&
                                '/t_'//trim(adjustl(int_to_str(t)))//&
                                '/f_'//trim(adjustl(int_to_str(f)))//&
                                '/m_'//trim(adjustl(int_to_str(m)))//&
                                '/init_'//trim(adjustl(int_to_str(i))),&
                                group_id, error)

                            call h5ltmake_dataset_double_f(file_id,&
                                'scatter_rates'//&
                                '/t_'//trim(adjustl(int_to_str(t)))//&
                                '/f_'//trim(adjustl(int_to_str(f)))//&
                                '/m_'//trim(adjustl(int_to_str(m)))//&
                                '/init_'//trim(adjustl(int_to_str(i)))//&
                                '/binned_i', &
                                size(dims2), dims2,&
                                binned_rate_init(i)%binned_rate(:, :, m, f, t), error)

                        end do
                    end do
                end do
            end do

            call h5fclose_f(file_id, error)
            call h5close_f(error)

        else

            call print_error_message('Output file : '//trim(filename)//' does NOT exist.', &
                verbose = verbose)
            stop

        end if

    end subroutine

end module
