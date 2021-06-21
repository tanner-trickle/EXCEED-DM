module dme_absorption
    !! Perform DM absorption rate calculations
    !! 
    !! Note : currently specific to vc transitions

    use mpi
    use hdf5
    use h5lt

    use prec
    use control_input
    use numerics_input
    use material_input
    use transition
    use DFT_parameters
    use particle_physics_abs

    use tran_form_calc
    use calc_PI

    use rate_calc_vector
    use rate_calc_ps
    use rate_calc_scalar

    implicit none

    !! building blocks of all absorption calculations
    complex(dp), allocatable :: tran_form(:, :, :)
        !! Dim : [n_init, n_fin, n_k]
        !!
        !! Scalar transition form factor for absorption calculations
        !!
        !! = <f| (p/m_e)^2 |i>
        !!
        !! Units : None

    complex(dp), allocatable :: tran_form_vec(:, :, :, :)
        !! Dim : [3, n_init, n_fin, n_k]
        !!
        !! Vector transition form factor for absorption calculations
        !! 
        !! = <f| (p_vec/m_e)  |i>
        !!
        !! Units : None

    real(dp), allocatable :: omega_iipk(:, :, :)
        !! Dim : [n_init, n_fin, n_k]
        !!
        !! Energy difference of the transition from i -> f at k
        !!
        !! Units : None

    complex(dp), allocatable :: tran_form_t(:, :)
        !! Dim : [n_tran_per_proc, n_k]
        !!
        !! Scalar transition form factor for absorption calculations
        !!
        !! Units : None

    complex(dp), allocatable :: tran_form_vec_t(:, :, :)
        !! Dim : [3, n_tran_per_proc, n_k]
        !!
        !! Vector transition form factor for absorption calculations
        !!
        !! Units : None

    real(dp), allocatable :: omega_iipk_t(:, :)
        !! Dim : [n_tran_per_proc, n_k]
        !!
        !! Energy difference of the transition from i -> f at k
        !!
        !! Units : None

    real(dp), allocatable :: abs_rate(:, :, :)
        !! Dim : [n_omega, n_widths, n_times] 
        !!
        !! Absorption rate per kg-year
        !!
        !! Units : None

contains

    subroutine run_dme_absorption(proc_id, root_process, &
            out_filename, nml_filename, DFT_input_filename, &
            n_proc, verbose)
        !! Absorption rate calculation

        implicit none

        integer :: proc_id
        integer :: root_process

        integer :: n_proc

        character(len=*) :: out_filename
        character(len=*) :: nml_filename
        character(len=*) :: DFT_input_filename

        logical, optional :: verbose

        integer :: status(MPI_STATUS_SIZE)
        integer :: tag = 0
        integer :: err

        integer :: t, tran_id, val_id, cond_id, i

        complex(dp), allocatable :: wfc_FT_i(:, :)
        complex(dp), allocatable :: wfc_FT_f(:, :)

        call load_numerics(nml_filename, verbose = verbose)
        call load_DFT_parameters(DFT_input_filename, verbose = verbose)
        call do_scissor_correction(band_gap, verbose = verbose)

        call load_absorption_input(nml_filename, verbose = verbose)

        call set_job_table(n_proc, n_init, n_fin, verbose=verbose)

        allocate(tran_form_t(n_tran_per_proc, n_k))
        tran_form_t = (0.0_dp, 0.0_dp)

        allocate(tran_form_vec_t(3, n_tran_per_proc, n_k))
        tran_form_vec_t = (0.0_dp, 0.0_dp)

        allocate(omega_iipk_t(n_tran_per_proc, n_k))
        omega_iipk_t = 0.0_dp

        ! these arrays are small enough for each processor to hold their
        ! own copy...
        allocate(tran_form(n_init, n_fin, n_k))
        tran_form = (0.0_dp, 0.0_dp)

        allocate(tran_form_vec(3, n_init, n_fin, n_k))
        tran_form_vec = (0.0_dp, 0.0_dp)

        allocate(omega_iipk(n_init, n_fin, n_k))
        omega_iipk = 0.0_dp

        allocate(wfc_FT_i(n_k, n_in_G))
        allocate(wfc_FT_f(n_k, n_in_G))

        ! compute tran_form data
        do t = 1, n_tran_per_proc

            tran_id = job_table(proc_id + 1, t)

            if ( tran_id .ne. 0 ) then

                val_id = tran_to_init_fin_id(tran_id, 1)
                cond_id = tran_to_init_fin_id(tran_id, 2) + n_val

                call get_in_wfc_FT(DFT_input_filename, val_id, wfc_FT_i)
                call get_in_wfc_FT(DFT_input_filename, cond_id, wfc_FT_f)

                call calc_tran_form(tran_form_t(t, :), tran_form_vec_t(:, t, :), &
                    omega_iipk_t(t, :), val_id, cond_id, wfc_FT_i, wfc_FT_f, verbose = verbose)

            end if

        end do

        ! update the transition form factors
        if ( proc_id .ne. root_process ) then

            call MPI_SEND(tran_form_t, &
               size(tran_form_t), MPI_DOUBLE_COMPLEX, root_process, tag, MPI_COMM_WORLD, err)

            call MPI_SEND(tran_form_vec_t, &
               size(tran_form_vec_t), MPI_DOUBLE_COMPLEX, root_process, tag, MPI_COMM_WORLD, err)

            call MPI_SEND(omega_iipk_t, &
               size(omega_iipk_t), MPI_DOUBLE, root_process, tag, MPI_COMM_WORLD, err)

        end if

        if ( proc_id == root_process ) then

            ! add main processors contribution
            call update_tran_form(tran_form_t, tran_form_vec_t, omega_iipk_t, proc_id, verbose = verbose)

            do i = 1, n_proc
                if ( (i - 1) .ne. root_process ) then

                    call MPI_RECV(tran_form_t, &
                       size(tran_form_t), MPI_DOUBLE_COMPLEX, i - 1, MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                    call MPI_RECV(tran_form_vec_t, &
                       size(tran_form_vec_t), MPI_DOUBLE_COMPLEX, i - 1, MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                    call MPI_RECV(omega_iipk_t, &
                       size(omega_iipk_t), MPI_DOUBLE, i - 1, MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                    ! add other processors contributions
                    call update_tran_form(tran_form_t, tran_form_vec_t, omega_iipk_t, i - 1, verbose = verbose)

                end if
            end do

        end if

        ! save the transition form factors
        if ( proc_id == root_process ) then

            call save_numerics(out_filename, verbose = verbose) 
            call save_DFT_parameters(out_filename, verbose = verbose)
            call save_absorption_input(out_filename, verbose = verbose)

            call save_tran_form(out_filename, verbose = verbose)

        end if

        ! send the full summed dataset back out to each processor
        if ( proc_id == root_process ) then

            do i = 1, n_proc
                if ( (i - 1) .ne. root_process ) then

                    call MPI_SEND(tran_form, &
                       size(tran_form), MPI_DOUBLE_COMPLEX, i - 1, tag, MPI_COMM_WORLD, err)

                    call MPI_SEND(tran_form_vec, &
                       size(tran_form_vec), MPI_DOUBLE_COMPLEX, i - 1, tag, MPI_COMM_WORLD, err)

                    call MPI_SEND(omega_iipk, &
                       size(omega_iipk), MPI_DOUBLE, i - 1, tag, MPI_COMM_WORLD, err)

                end if
            end do

        end if

        ! recv the full dataset from the main processor
        if ( proc_id /= root_process ) then

            call MPI_RECV(tran_form, &
               size(tran_form), MPI_DOUBLE_COMPLEX, root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

            call MPI_RECV(tran_form_vec, &
               size(tran_form_vec), MPI_DOUBLE_COMPLEX, root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

            call MPI_RECV(omega_iipk, &
               size(omega_iipk), MPI_DOUBLE, root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

        end if

        ! compute Pi's
        call compute_self_energies(nml_filename, tran_form, &
            tran_form_vec, omega_iipk, n_init, n_fin, n_k, n_proc, &
            proc_id, root_process, verbose = verbose)

        if ( proc_id == root_process ) then

            call save_self_energies(out_filename, verbose = verbose)

        end if

        if ( proc_id == root_process ) then

            call load_particle_physics_abs(nml_filename, verbose = verbose)

            ! compute rate for different models
            allocate(abs_rate(n_omega, n_widths, n_time))
            abs_rate = 0.0_dp

            if ( trim(calc_mode) == 'vector' ) then

                call calc_rate_vector(pi_11_mat, abs_rate, verbose = verbose)

            else if ( trim(calc_mode) == 'ps' ) then

                call calc_rate_ps(pi_11_mat, abs_rate, verbose = verbose)

            else if ( trim(calc_mode) == 'scalar' ) then

                call calc_rate_scalar(pi_11_mat, pi_1v, pi_v1, pi_vv, abs_rate, &
                    verbose = verbose)

            end if

            call save_abs_rate(out_filename, verbose = verbose)

        end if

    end subroutine

    subroutine save_abs_rate(filename, verbose)

        implicit none
        character(len=*) :: filename

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id

        logical, optional :: verbose
        logical :: file_exists

        integer(HSIZE_T) :: dims1(1) = [1]
        integer(HSIZE_T) :: dims2(2)

        integer :: error

        integer :: p, t
        integer :: i, fin
        character(len=64) :: t_str, p_str

        if ( verbose ) then

            print*, '    Saving absorption rate...'
            print*

        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'abs_rate', group_id, error)

            do t = 1, n_time

                write(t_str, *) t
                t_str = trim(adjustl(t_str))

                call h5gcreate_f(file_id,& 
                    'abs_rate'//&
                    '/t_'//trim(t_str),&
                    group_id, error)
                
                do p = 1, n_widths

                    write(p_str, *) p
                    p_str = trim(adjustl(p_str))

                    dims1 = [n_omega]

                    call h5ltmake_dataset_double_f(file_id,&
                        'abs_rate'//&
                        '/t_'//trim(t_str)//&
                        '/width_'//trim(p_str),&
                        size(dims1), dims1,&
                        abs_rate(:, p, t), error)

                end do 

            end do

            call h5fclose_f(file_id, error)
            call h5close_f(error)

            if ( verbose ) then
                print*, '----------'
                print*
            end if

        else

            if ( verbose ) then

                print*, '!! ERROR !!'
                print*
                print*, '   Output file : ', trim(filename), ' does NOT exist.'
                print*
                print*, '!!!!!!!!!!!'
                print*

            end if

            stop

        end if

    end subroutine

    subroutine save_tran_form(filename, verbose)

        implicit none

        character(len=*) :: filename

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id

        logical, optional :: verbose
        logical :: file_exists

        integer(HSIZE_T) :: dims1(1) = [1]
        integer(HSIZE_T) :: dims2(2)

        integer :: error

        integer :: m, f, t
        integer :: i, fin
        character(len=64) :: i_str, fin_str

        if ( verbose ) then

            print*, 'Saving absorption transition form factors...'
            print*

        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'tran_form_data', group_id, error)

            do i = 1, n_init

                write(i_str, *) i
                i_str = trim(adjustl(i_str))

                call h5gcreate_f(file_id,& 
                    'tran_form_data'//&
                    '/init_'//trim(i_str),&
                    group_id, error)

                do fin = 1, n_fin

                    write(fin_str, *) fin
                    fin_str = trim(adjustl(fin_str))

                    call h5gcreate_f(file_id,& 
                        'tran_form_data'//&
                        '/init_'//trim(i_str)//&
                        '/fin_'//trim(fin_str),&
                        group_id, error)

                    dims1 = [n_k]

                    call h5ltmake_dataset_double_f(file_id,&
                        'tran_form_data'//&
                        '/init_'//trim(i_str)//&
                        '/fin_'//trim(fin_str)//&
                        '/tran_form_r', &
                        size(dims1), dims1,&
                        real(tran_form(i, fin, :)), error)

                    call h5ltmake_dataset_double_f(file_id,&
                        'tran_form_data'//&
                        '/init_'//trim(i_str)//&
                        '/fin_'//trim(fin_str)//&
                        '/tran_form_c', &
                        size(dims1), dims1,&
                        aimag(tran_form(i, fin, :)), error)

                    call h5ltmake_dataset_double_f(file_id,&
                        'tran_form_data'//&
                        '/init_'//trim(i_str)//&
                        '/fin_'//trim(fin_str)//&
                        '/omega_iip', &
                        size(dims1), dims1,&
                        omega_iipk(i, fin, :), error)

                    dims2 = [3, n_k]

                    call h5ltmake_dataset_double_f(file_id,&
                        'tran_form_data'//&
                        '/init_'//trim(i_str)//&
                        '/fin_'//trim(fin_str)//&
                        '/tran_form_vec_r', &
                        size(dims2), dims2,&
                        real(tran_form_vec(:, i, fin, :)), error)

                    call h5ltmake_dataset_double_f(file_id,&
                        'tran_form_data'//&
                        '/init_'//trim(i_str)//&
                        '/fin_'//trim(fin_str)//&
                        '/tran_form_vec_c', &
                        size(dims2), dims2,&
                        aimag(tran_form_vec(:, i, fin, :)), error)

                end do
            end do


            call h5fclose_f(file_id, error)
            call h5close_f(error)

            if ( verbose ) then
                print*, '----------'
                print*
            end if

        else

            if ( verbose ) then

                print*, '!! ERROR !!'
                print*
                print*, '   Output file : ', trim(filename), ' does NOT exist.'
                print*
                print*, '!!!!!!!!!!!'
                print*

            end if

            stop

        end if

    end subroutine

    subroutine update_tran_form(tran_form_t, tran_form_vec_t, omega_iipk_t, proc_id, verbose)

        implicit none

        complex(dp) :: tran_form_t(n_tran_per_proc, n_k)
        complex(dp) :: tran_form_vec_t(3, n_tran_per_proc, n_k)

        real(dp) :: omega_iipk_t(n_tran_per_proc, n_k)

        integer :: t, tran_id, init_id, fin_id, proc_id

        logical, optional :: verbose

        do t = 1, n_tran_per_proc

            tran_id = job_table(proc_id + 1, t)

            init_id = tran_to_init_fin_id(tran_id, 1)
            fin_id = tran_to_init_fin_id(tran_id, 2)

            if ( tran_id .ne. 0 ) then

                tran_form(init_id, fin_id, :) = tran_form_t(t, :)
                tran_form_vec(:, init_id, fin_id, :) = tran_form_vec_t(:, t, :)

                omega_iipk(init_id, fin_id, :) = omega_iipk_t(t, :)

            end if

        end do

    end subroutine

end module
