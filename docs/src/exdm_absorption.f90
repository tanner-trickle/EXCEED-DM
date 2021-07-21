module exdm_absorption
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
    use rate_calc_scalar_LO

    implicit none

    !! building blocks of all absorption calculations
    complex(dp), allocatable :: tran_form_1_no_spin(:, :, :)
        !! Dim : [n_init, n_fin, n_k]
        !!
        !! Units : None

    complex(dp), allocatable :: tran_form_1_spin(:, :, :, :, :)
        !! Dim : [n_init, n_fin, n_k, s, sp]
        !!
        !! Units : None

    complex(dp), allocatable :: tran_form_v_no_spin(:, :, :, :)
        !! Dim : [3, n_init, n_fin, n_k]
        !!
        !! Units : None

    complex(dp), allocatable :: tran_form_v_spin(:, :, :, :, :, :)
        !! Dim : [3, n_init, n_fin, n_k, s, sp]
        !!
        !! Units : None

    complex(dp), allocatable :: tran_form_v2_no_spin(:, :, :)
        !! Dim : [n_init, n_fin, n_k]
        !!
        !! Units : None

    complex(dp), allocatable :: tran_form_v2_spin(:, :, :, :, :)
        !! Dim : [n_init, n_fin, n_k, s, sp]
        !!
        !! Units : None

    real(dp), allocatable :: omega_iipk(:, :, :)
        !! Dim : [n_init, n_fin, n_k]
        !!
        !! Energy difference of the transition from i -> f at k
        !!
        !! Units : None

    !! processor specific
    complex(dp), allocatable :: tran_form_1_no_spin_t(:, :)
        !! Dim : [n_tran_per_proc, n_k]
        !!
        !! Units : None

    complex(dp), allocatable :: tran_form_1_spin_t(:, :, :, :)
        !! Dim : [n_tran_per_proc, n_k, s, sp]
        !!
        !! Units : None

    complex(dp), allocatable :: tran_form_v_no_spin_t(:, :, :)
        !! Dim : [3, n_tran_per_proc, n_k]
        !!
        !! Units : None

    complex(dp), allocatable :: tran_form_v_spin_t(:, :, :, :, :)
        !! Dim : [3, n_tran_per_proc, n_k, s, sp]
        !!
        !! Units : None

    complex(dp), allocatable :: tran_form_v2_no_spin_t(:, :)
        !! Dim : [n_tran_per_proc, n_k]
        !!
        !! Units : None

    complex(dp), allocatable :: tran_form_v2_spin_t(:, :, :, :)
        !! Dim : [n_tran_per_proc, n_k, s, sp]
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

    real(dp), allocatable :: abs_rate_buff(:, :, :)
        !! Dim : [n_omega, n_widths, n_times] 
        !!
        !! Absorption rate per kg-year
        !!
        !! Units : None

    integer :: n_v_per_proc

    real(dp), allocatable :: v_mesh(:, :)
        !! Dim : [ n_v_mag*n_v_theta*n_v_phi, 3 ]
        !!
        !! Mesh of velocity vectors
        !!
        !! Units : None

    real(dp) :: v_vec(3)

    real(dp) :: abs_v_max

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
        integer :: v_id

        complex(dp), allocatable :: wfc_FT_i(:, :)
        complex(dp), allocatable :: wfc_FT_f(:, :)

        complex(dp), allocatable :: wfc_FT_i_spin(:, :, :)
        complex(dp), allocatable :: wfc_FT_f_spin(:, :, :)

        call load_numerics(nml_filename, verbose = verbose)
        call load_DFT_parameters(DFT_input_filename, verbose = verbose)
        call do_scissor_correction(band_gap, verbose = verbose)

        call load_absorption_input(nml_filename, verbose = verbose)

        call set_job_table(n_proc, n_init, n_fin, verbose=verbose)

        !! allocate/initialize variables

        if ( proc_id == root_process ) then

            allocate(tran_form_1_no_spin(n_init, n_fin, n_k))
            tran_form_1_no_spin  = (0.0_dp, 0.0_dp)

            allocate(tran_form_1_spin(n_init, n_fin, n_k, 2, 2))
            tran_form_1_spin     = (0.0_dp, 0.0_dp)

            allocate(tran_form_v_no_spin(3, n_init, n_fin, n_k))
            tran_form_v_no_spin  = (0.0_dp, 0.0_dp)

            allocate(tran_form_v_spin(3, n_init, n_fin, n_k, 2, 2))
            tran_form_v_spin     = (0.0_dp, 0.0_dp)

            allocate(tran_form_v2_no_spin(n_init, n_fin, n_k))
            tran_form_v2_no_spin = (0.0_dp, 0.0_dp)

            allocate(tran_form_v2_spin(n_init, n_fin, n_k, 2, 2))
            tran_form_v2_spin    = (0.0_dp, 0.0_dp)

            allocate(omega_iipk(n_init, n_fin, n_k))
            omega_iipk           = 0.0_dp

        end if

        !! processor specific

        allocate(tran_form_1_no_spin_t(n_tran_per_proc, n_k))
        tran_form_1_no_spin_t  = (0.0_dp, 0.0_dp)

        allocate(tran_form_1_spin_t(n_tran_per_proc, n_k, 2, 2))
        tran_form_1_spin_t     = (0.0_dp, 0.0_dp)

        allocate(tran_form_v_no_spin_t(3, n_tran_per_proc, n_k))
        tran_form_v_no_spin_t  = (0.0_dp, 0.0_dp)

        allocate(tran_form_v_spin_t(3, n_tran_per_proc, n_k, 2, 2))
        tran_form_v_spin_t     = (0.0_dp, 0.0_dp)

        allocate(tran_form_v2_no_spin_t(n_tran_per_proc, n_k))
        tran_form_v2_no_spin_t = (0.0_dp, 0.0_dp)

        allocate(tran_form_v2_spin_t(n_tran_per_proc, n_k, 2, 2))
        tran_form_v2_spin_t    = (0.0_dp, 0.0_dp)

        allocate(omega_iipk_t(n_tran_per_proc, n_k))
        omega_iipk_t           = 0.0_dp

        !! allow for spin dependent wave functions
        if ( include_spin ) then

            allocate(wfc_FT_i_spin(n_k, n_in_G, 2))
            allocate(wfc_FT_f_spin(n_k, n_in_G, 2))

        else

            allocate(wfc_FT_i(n_k, n_in_G))
            allocate(wfc_FT_f(n_k, n_in_G))

        end if

        ! compute transition form factors
        do t = 1, n_tran_per_proc

            tran_id = job_table(proc_id + 1, t)

            if ( tran_id /= 0 ) then

                val_id = tran_to_init_fin_id(tran_id, 1)
                cond_id = tran_to_init_fin_id(tran_id, 2) + n_val

                omega_iipk_t(t, :) = energy_bands(:, cond_id) - energy_bands(:, val_id)

                if ( include_spin ) then

                    call get_in_wfc_FT(DFT_input_filename, val_id, wfc_FT_i_spin)
                    call get_in_wfc_FT(DFT_input_filename, cond_id, wfc_FT_f_spin)

                    !! compute transition form factors
                    call calc_tran_form_1(tran_form_1_spin_t(t, :, :, :), &
                        val_id, cond_id, wfc_FT_i_spin, wfc_FT_f_spin, verbose = verbose)

                    call calc_tran_form_v(tran_form_v_spin_t(:, t, :, :, :), &
                        val_id, cond_id, wfc_FT_i_spin, wfc_FT_f_spin, verbose = verbose)

                    call calc_tran_form_v2(tran_form_v2_spin_t(t, :, :, :), &
                        val_id, cond_id, wfc_FT_i_spin, wfc_FT_f_spin, verbose = verbose)

                else

                    call get_in_wfc_FT(DFT_input_filename, val_id, wfc_FT_i)
                    call get_in_wfc_FT(DFT_input_filename, cond_id, wfc_FT_f)

                    !! compute transition form factors
                    call calc_tran_form_1(tran_form_1_no_spin_t(t, :), &
                        val_id, cond_id, wfc_FT_i, wfc_FT_f, verbose = verbose)

                    call calc_tran_form_v(tran_form_v_no_spin_t(:, t, :), &
                        val_id, cond_id, wfc_FT_i, wfc_FT_f, verbose = verbose)

                    call calc_tran_form_v2(tran_form_v2_no_spin_t(t, :), &
                        val_id, cond_id, wfc_FT_i, wfc_FT_f, verbose = verbose)

                end if

            end if

        end do

        ! update the transition form factors
        if ( proc_id /= root_process ) then

            if ( include_spin ) then

                call MPI_SEND(tran_form_1_spin_t, &
                   size(tran_form_1_spin_t), MPI_DOUBLE_COMPLEX, &
                   root_process, tag, MPI_COMM_WORLD, err)

                call MPI_SEND(tran_form_v_spin_t, &
                   size(tran_form_v_spin_t), MPI_DOUBLE_COMPLEX, &
                   root_process, tag, MPI_COMM_WORLD, err)

                call MPI_SEND(tran_form_v2_spin_t, &
                   size(tran_form_v2_spin_t), MPI_DOUBLE_COMPLEX, &
                   root_process, tag, MPI_COMM_WORLD, err)

            else

                call MPI_SEND(tran_form_1_no_spin_t, &
                   size(tran_form_1_no_spin_t), MPI_DOUBLE_COMPLEX, &
                   root_process, tag, MPI_COMM_WORLD, err)

                call MPI_SEND(tran_form_v_no_spin_t, &
                   size(tran_form_v_no_spin_t), MPI_DOUBLE_COMPLEX, &
                   root_process, tag, MPI_COMM_WORLD, err)

                call MPI_SEND(tran_form_v2_no_spin_t, &
                   size(tran_form_v2_no_spin_t), MPI_DOUBLE_COMPLEX, &
                   root_process, tag, MPI_COMM_WORLD, err)

            end if

            call MPI_SEND(omega_iipk_t, &
               size(omega_iipk_t), MPI_DOUBLE, root_process, tag, MPI_COMM_WORLD, err)

        end if

        if ( proc_id == root_process ) then

            ! main processors contribution
            call update_tran_form(tran_form_1_no_spin_t, &
                                  tran_form_v_no_spin_t, &
                                  tran_form_v2_no_spin_t, &
                                  tran_form_1_spin_t, &
                                  tran_form_v_spin_t, &
                                  tran_form_v2_spin_t, &
                                  omega_iipk_t, proc_id, verbose = verbose)

            do i = 1, n_proc
                if ( (i - 1) .ne. root_process ) then

                    if ( include_spin ) then

                        call MPI_RECV(tran_form_1_spin_t, &
                           size(tran_form_1_spin_t), MPI_DOUBLE_COMPLEX, i - 1, MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                        call MPI_RECV(tran_form_v_spin_t, &
                           size(tran_form_v_spin_t), MPI_DOUBLE_COMPLEX, i - 1, MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                        call MPI_RECV(tran_form_v2_spin_t, &
                           size(tran_form_v2_spin_t), MPI_DOUBLE_COMPLEX, i - 1, MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                    else

                        call MPI_RECV(tran_form_1_no_spin_t, &
                           size(tran_form_1_no_spin_t), MPI_DOUBLE_COMPLEX, i - 1, MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                        call MPI_RECV(tran_form_v_no_spin_t, &
                           size(tran_form_v_no_spin_t), MPI_DOUBLE_COMPLEX, i - 1, MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                        call MPI_RECV(tran_form_v2_no_spin_t, &
                           size(tran_form_v2_no_spin_t), MPI_DOUBLE_COMPLEX, i - 1, MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                    end if

                    call MPI_RECV(omega_iipk_t, &
                       size(omega_iipk_t), MPI_DOUBLE, i - 1, MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                    ! other processors contributions
                    call update_tran_form(tran_form_1_no_spin_t, &
                                          tran_form_v_no_spin_t, &
                                          tran_form_v2_no_spin_t, &
                                          tran_form_1_spin_t, &
                                          tran_form_v_spin_t, &
                                          tran_form_v2_spin_t, &
                                          omega_iipk_t, i - 1, verbose = verbose)

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

        ! compute Pi's
        call compute_self_energies_setup(proc_id, root_process, verbose = verbose)

        if ( include_spin ) then

            call compute_self_energies(nml_filename, &
                tran_form_1_spin_t, &
                tran_form_v_spin_t, &
                tran_form_v2_spin_t, &
                omega_iipk_t, n_init, n_fin, n_k, n_proc, &
                proc_id, root_process, verbose = verbose)

        else

            call compute_self_energies(nml_filename, &
                tran_form_1_no_spin_t, &
                tran_form_v_no_spin_t, &
                tran_form_v2_no_spin_t, &
                omega_iipk_t, n_init, n_fin, n_k, n_proc, &
                proc_id, root_process, verbose = verbose)

        end if

        call comm_self_energies(proc_id, root_process, n_proc, verbose = verbose)

        if ( proc_id == root_process ) then

            call save_self_energies(out_filename, verbose = verbose)

        end if

        call bcast_self_energies(proc_id, root_process, verbose = verbose)

        !! compute absorption rates in parallel across v points
        call load_particle_physics_abs(nml_filename, verbose = verbose)

        !! reset job_table to parallelize over v points
        !! TODO: create a new job_table when transition is rewritten...
        call set_job_table(n_proc, 1, n_v_mag*n_v_theta*n_v_phi, verbose=verbose)

        !! number of v points each processor has to compute for.
        n_v_per_proc = n_tran_per_proc

        call create_v_mesh(n_v_mag, n_v_theta, n_v_phi)

        allocate(abs_rate(n_omega, n_widths, n_time))
        abs_rate = 0.0_dp

        !! compute abs_rate for each v point

        do t = 1, n_v_per_proc

            v_id = job_table(proc_id + 1, t)

            if ( v_id /= 0 ) then

                v_vec = v_mesh(v_id, :)

                if ( trim(calc_mode) == 'vector' ) then

                    call calc_rate_vector(pi_vi_vj, v_vec, &
                        abs_v_max, abs_rate, verbose = verbose)

                else if ( trim(calc_mode) == 'ps' ) then

                    call calc_rate_ps(pi_vi_vj, v_vec, &
                        abs_v_max, abs_rate, verbose = verbose)

                else if ( trim(calc_mode) == 'scalar' ) then

                    call calc_rate_scalar(pi_1_1_mat, pi_v2_v2, v_vec, abs_v_max, abs_rate, &
                        verbose = verbose)

                else if ( trim(calc_mode) == 'scalar_LO' ) then

                    call calc_rate_scalar_LO(pi_v2_v2, v_vec, abs_v_max, abs_rate, &
                        verbose = verbose)

                end if

            end if

        end do

        !! send back to the main processor

        if ( proc_id /= root_process ) then

            call MPI_SEND(abs_rate, &
                size(abs_rate), MPI_DOUBLE, root_process, &
                tag, MPI_COMM_WORLD, err)

        end if

        if ( proc_id == root_process ) then

            allocate(abs_rate_buff(n_omega, n_widths, n_time))

            do i = 1, n_proc

                abs_rate_buff = 0.0_dp

                if ( i - 1 /= root_process ) then

                    call MPI_RECV(abs_rate_buff, &
                        size(abs_rate_buff), MPI_DOUBLE, i - 1, &
                        MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                    abs_rate = abs_rate + abs_rate_buff

                end if

            end do

        end if

        if ( proc_id == root_process ) then

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
        integer(HSIZE_T) :: dims3(3)
        integer(HSIZE_T) :: dims4(4)

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
                        '/omega_iip', &
                        size(dims1), dims1,&
                        omega_iipk(i, fin, :), error)

                    !! These are zero

                    ! call h5ltmake_dataset_double_f(file_id,&
                    !     'tran_form_data'//&
                    !     '/init_'//trim(i_str)//&
                    !     '/fin_'//trim(fin_str)//&
                    !     '/tran_form_1_r', &
                    !     size(dims1), dims1,&
                    !     real(tran_form_1_no_spin(i, fin, :)), error)

                    ! call h5ltmake_dataset_double_f(file_id,&
                    !     'tran_form_data'//&
                    !     '/init_'//trim(i_str)//&
                    !     '/fin_'//trim(fin_str)//&
                    !     '/tran_form_1_c', &
                    !     size(dims1), dims1,&
                    !     aimag(tran_form_1_no_spin(i, fin, :)), error)

                    call h5ltmake_dataset_double_f(file_id,&
                        'tran_form_data'//&
                        '/init_'//trim(i_str)//&
                        '/fin_'//trim(fin_str)//&
                        '/tran_form_v2_r', &
                        size(dims1), dims1,&
                        real(tran_form_v2_no_spin(i, fin, :)), error)

                    call h5ltmake_dataset_double_f(file_id,&
                        'tran_form_data'//&
                        '/init_'//trim(i_str)//&
                        '/fin_'//trim(fin_str)//&
                        '/tran_form_v2_c', &
                        size(dims1), dims1,&
                        aimag(tran_form_v2_no_spin(i, fin, :)), error)

                    dims2 = [3, n_k]

                    call h5ltmake_dataset_double_f(file_id,&
                        'tran_form_data'//&
                        '/init_'//trim(i_str)//&
                        '/fin_'//trim(fin_str)//&
                        '/tran_form_v_r', &
                        size(dims2), dims2,&
                        real(tran_form_v_no_spin(:, i, fin, :)), error)

                    call h5ltmake_dataset_double_f(file_id,&
                        'tran_form_data'//&
                        '/init_'//trim(i_str)//&
                        '/fin_'//trim(fin_str)//&
                        '/tran_form_v_c', &
                        size(dims2), dims2,&
                        aimag(tran_form_v_no_spin(:, i, fin, :)), error)

                    if ( include_spin ) then

                        dims3 = [n_k, 2, 2]

                        call h5ltmake_dataset_double_f(file_id,&
                            'tran_form_data'//&
                            '/init_'//trim(i_str)//&
                            '/fin_'//trim(fin_str)//&
                            '/tran_form_1_spin_r', &
                            size(dims1), dims1,&
                            real(tran_form_1_spin(i, fin, :, :, :)), error)

                        call h5ltmake_dataset_double_f(file_id,&
                            'tran_form_data'//&
                            '/init_'//trim(i_str)//&
                            '/fin_'//trim(fin_str)//&
                            '/tran_form_1_spin_c', &
                            size(dims1), dims1,&
                            aimag(tran_form_1_spin(i, fin, :, :, :)), error)

                        call h5ltmake_dataset_double_f(file_id,&
                            'tran_form_data'//&
                            '/init_'//trim(i_str)//&
                            '/fin_'//trim(fin_str)//&
                            '/tran_form_v2_spin_r', &
                            size(dims1), dims1,&
                            real(tran_form_v2_spin(i, fin, :, :, :)), error)

                        call h5ltmake_dataset_double_f(file_id,&
                            'tran_form_data'//&
                            '/init_'//trim(i_str)//&
                            '/fin_'//trim(fin_str)//&
                            '/tran_form_v2_spin_c', &
                            size(dims1), dims1,&
                            aimag(tran_form_v2_spin(i, fin, :, :, :)), error)

                        dims4 = [3, n_k, 2, 2]

                        call h5ltmake_dataset_double_f(file_id,&
                            'tran_form_data'//&
                            '/init_'//trim(i_str)//&
                            '/fin_'//trim(fin_str)//&
                            '/tran_form_v_spin_r', &
                            size(dims2), dims2,&
                            real(tran_form_v_spin(:, i, fin, :, :, :)), error)

                        call h5ltmake_dataset_double_f(file_id,&
                            'tran_form_data'//&
                            '/init_'//trim(i_str)//&
                            '/fin_'//trim(fin_str)//&
                            '/tran_form_v_spin_c', &
                            size(dims2), dims2,&
                            aimag(tran_form_v_spin(:, i, fin, :, :, :)), error)

                    end if

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

    subroutine update_tran_form(tran_form_1_no_spin_t, &
                                tran_form_v_no_spin_t, &
                                tran_form_v2_no_spin_t, &
                                tran_form_1_spin_t, &
                                tran_form_v_spin_t, &
                                tran_form_v2_spin_t, &
                                omega_iipk_t, proc_id, verbose)

        implicit none

        complex(dp) :: tran_form_1_no_spin_t(n_tran_per_proc, n_k)
        complex(dp) :: tran_form_v_no_spin_t(3, n_tran_per_proc, n_k)
        complex(dp) :: tran_form_v2_no_spin_t(n_tran_per_proc, n_k)

        complex(dp) :: tran_form_1_spin_t(n_tran_per_proc, n_k, 2, 2)
        complex(dp) :: tran_form_v_spin_t(3, n_tran_per_proc, n_k, 2, 2)
        complex(dp) :: tran_form_v2_spin_t(n_tran_per_proc, n_k, 2, 2)

        real(dp) :: omega_iipk_t(n_tran_per_proc, n_k)

        integer :: t, tran_id, init_id, fin_id, proc_id

        logical, optional :: verbose

        do t = 1, n_tran_per_proc

            tran_id = job_table(proc_id + 1, t)

            init_id = tran_to_init_fin_id(tran_id, 1)
            fin_id = tran_to_init_fin_id(tran_id, 2)

            if ( tran_id .ne. 0 ) then

                tran_form_1_no_spin(init_id, fin_id, :) = tran_form_1_no_spin_t(t, :)
                tran_form_v_no_spin(:, init_id, fin_id, :) = tran_form_v_no_spin_t(:, t, :)
                tran_form_v2_no_spin(init_id, fin_id, :) = tran_form_v2_no_spin_t(t, :)

                tran_form_1_spin(init_id, fin_id, :, :, :) = tran_form_1_spin_t(t, :, :, :)
                tran_form_v_spin(:, init_id, fin_id, :, :, :) = tran_form_v_spin_t(:, t, :, :, :)
                tran_form_v2_spin(init_id, fin_id, :, :, :) = tran_form_1_spin_t(t, :, :, :)

                omega_iipk(init_id, fin_id, :) = omega_iipk_t(t, :)

            end if

        end do

    end subroutine

    subroutine create_v_mesh(n_v_mag, n_v_theta, n_v_phi, verbose)

        implicit none

        integer :: n_v_mag, n_v_theta, n_v_phi
        integer :: n_v

        logical, optional :: verbose

        integer :: v, a

        real(dp) :: v_mag, v_theta, v_phi
        real(dp) :: v_angular_mesh(n_v_theta*n_v_phi, 2)

        integer :: v_id

        n_v = n_v_mag*n_v_theta*n_v_phi

        allocate(v_mesh(n_v, 3))
        v_mesh = 0.0_dp

        abs_v_max = vEsc + vE

        v_angular_mesh = generate_uniform_points_on_sphere(n_v_theta, n_v_phi)

        v_id = 0

        do v = 1, n_v_mag

            v_mag = abs_v_max*(v - 1.0_dp)/max(1.0_dp, n_v_mag - 1.0_dp)
            
            do a = 1, n_v_theta*n_v_phi

                v_theta = v_angular_mesh(a, 1) 
                v_phi = v_angular_mesh(a, 2) 

                v_id = v_id + 1

                v_mesh(v_id, 1) = v_mag*sin(v_theta)*cos(v_phi)
                v_mesh(v_id, 2) = v_mag*sin(v_theta)*sin(v_phi)
                v_mesh(v_id, 3) = v_mag*cos(v_theta)

            end do

        end do

    end subroutine

end module
