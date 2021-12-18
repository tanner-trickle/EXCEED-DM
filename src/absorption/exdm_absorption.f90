module exdm_absorption
    !! Computes DM absorption rates.
    !! 
    !! Note : currently specific to valence to conduction transitions

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
    use PW_dataset_type
    use width_parameters_type

    use numerics_abs
    use abs_tran_form_calc
    use abs_calc_PI

    use rate_calc_vector 
    use rate_calc_ps
    use rate_calc_scalar
    use rate_calc_scalar_LO

    implicit none

contains

    subroutine run_exdm_absorption(proc_id, root_process, n_proc, &
            io_files, main_control, target_mat, expt, dm_model, & 
            verbose)
        !! Computes DM absorption rates.
        !! 
        !! Note : currently specific to valence to conduction transitions

        implicit none

        integer :: proc_id
        integer :: root_process
        integer :: n_proc
        type(io_files_t) :: io_files
        type(control_t) :: main_control
        type(material_t) :: target_mat
        type(expt_t) :: expt
        type(dm_model_t) :: dm_model
        logical, optional :: verbose

        type(PW_dataset_t) :: PW_dataset
        type(width_parameters_t) :: widths
        type(numerics_abs_t) :: numerics

        type(parallel_manager_t) :: ik_manager
        type(parallel_manager_t) :: v_manager
        integer, allocatable :: job_id_to_ik(:, :)

        complex(dp), allocatable :: tran_form_1_no_spin_job(:, :)
            !! Dim : [n_jobs_per_proc, n_cond_max]
            !!
            !! Units : None
        complex(dp), allocatable :: tran_form_1_spin_job(:, :, :, :)
            !! Dim : [n_jobs_per_proc, n_cond_max, 2, 2]
            !!
            !! Units : None
        complex(dp), allocatable :: tran_form_v_no_spin_job(:, :, :)
            !! Dim : [3, n_jobs_per_proc, n_cond_max]
            !!
            !! Units : None
        complex(dp), allocatable :: tran_form_v_spin_job(:, :, :, :, :)
            !! Dim : [3, n_jobs_per_proc, n_cond_max, 2, 2]
            !!
            !! Units : None
        complex(dp), allocatable :: tran_form_v2_no_spin_job(:, :)
            !! Dim : [n_jobs_per_proc, n_cond_max]
            !!
            !! Units : None
        complex(dp), allocatable :: tran_form_v2_spin_job(:, :, :, :)
            !! Dim : [n_jobs_per_proc, n_cond_max, 2, 2]
            !!
            !! Units : None
        complex(dp), allocatable :: tran_form_vs_spin_job(:, :)
            !! Dim : [n_jobs_per_proc, n_cond_max]
            !!
            !! Units : None

        complex(dp), allocatable :: pi_v2_v2(:, :)
            !! Dim : [n_mX, n_widths]
            !! 
            !! \( \Pi_{\bar v^2, \bar v^2} \) - self energy with two v^2 insertions.
            !!
            !! Units : eV^2
        complex(dp), allocatable :: pi_1_1_mat(:, :, :, :)
            !! Dim : [3, 3, n_omega, n_widths] 
            !!
            !! \( \mathbf{\Pi}_{11} \) self energy, without q_vec's, i.e.
            !!
            !! $$\begin{align*}
            !!      \Pi_{11} = \frac{\mathbf{q}}{m_e} \cdot \mathbf{Pi}_{11} \cdot \frac{ \mathbf{q} }{m_e}
            !! \end{align*}$$
            !!
            !! Units : eV^2
        complex(dp), allocatable :: pi_vi_vj(:, :, :, :)
            !! Dim : [3, 3, n_omega, n_widths] 
            !!
            !! \( \Pi_{v^i, v^j} \) self energy.
            !!
            !! Units : eV^2

        complex(dp), allocatable :: pi_vs_vs(:, :)
            !! Dim : [n_mX, n_widths]
            !! 
            !! \( \Pi_{\mathbf{v} \cdot \mathbf{\sigma}, \mathbf{v} \cdot \mathbf{\sigma}} \) - self energy with two v^2 insertions.
            !!
            !! Units : eV^2

        real(dp), allocatable :: abs_rate(:, :, :)
            !! Dim : [ dm_model%n_mX, widths%n, expt%n_time ]
            !!
            !! Absorption rate, assuming the mediator-electron coupling, \( g = 1 \).
            !!
            !! Units : None

        real(dp) :: v_vec(3)
        integer :: j, v_id

        call PW_dataset%load(io_files%PW_data_filename, verbose = verbose)
        call PW_dataset%do_scissor_correction(target_mat%band_gap, verbose = verbose)

        call widths%load(io_files%nml_input_filename, verbose = verbose)

        call numerics%load(io_files%nml_input_filename, &
            PW_dataset%n_val, PW_dataset%n_cond, dm_model, verbose = verbose)

        ! Compute the relevant transition form factors

        ! parallelize over {i, k}
        call ik_manager%init(PW_dataset%n_k*numerics%n_val_max, verbose = verbose)

        allocate(job_id_to_ik(ik_manager%n_jobs, 4))
        job_id_to_ik = 0

        call numerics%create_val_id_list(PW_dataset%n_val)
        call numerics%create_k_id_list(PW_dataset%n_k)

        call ik_manager%create_job_to_2d_ID_table(&
            numerics%val_id_list, &
            numerics%k_id_list, &
            job_id_to_ik, verbose = verbose)

        ! allocate/initialize processor specific variables

        allocate(tran_form_1_no_spin_job(ik_manager%n_jobs_per_proc, numerics%n_cond_max))
        tran_form_1_no_spin_job  = (0.0_dp, 0.0_dp)
        allocate(tran_form_1_spin_job(ik_manager%n_jobs_per_proc, numerics%n_cond_max, 2, 2))
        tran_form_1_spin_job  = (0.0_dp, 0.0_dp)
        allocate(tran_form_v_no_spin_job(3, ik_manager%n_jobs_per_proc, numerics%n_cond_max))
        tran_form_v_no_spin_job  = (0.0_dp, 0.0_dp)
        allocate(tran_form_v_spin_job(3, ik_manager%n_jobs_per_proc, numerics%n_cond_max, 2, 2))
        tran_form_v_spin_job  = (0.0_dp, 0.0_dp)
        allocate(tran_form_v2_no_spin_job(ik_manager%n_jobs_per_proc, numerics%n_cond_max))
        tran_form_v2_no_spin_job  = (0.0_dp, 0.0_dp)
        allocate(tran_form_v2_spin_job(ik_manager%n_jobs_per_proc, numerics%n_cond_max, 2, 2))
        tran_form_v2_spin_job  = (0.0_dp, 0.0_dp)
        allocate(tran_form_vs_spin_job(ik_manager%n_jobs_per_proc, numerics%n_cond_max))
        tran_form_vs_spin_job  = (0.0_dp, 0.0_dp)

        ! compute the transition form factors
        call calc_abs_tran_form(proc_id, root_process, &
                                  tran_form_1_no_spin_job, tran_form_1_spin_job, &
                                  tran_form_v_no_spin_job, tran_form_v_spin_job, &
                                  tran_form_v2_no_spin_job, tran_form_v2_spin_job, &
                                  tran_form_vs_spin_job, &
                                  PW_dataset, ik_manager, job_id_to_ik, numerics, &
                                  io_files%out_filename, &
                                  verbose = verbose)

        ! compute Pi's

        ! allocate pi variables
        allocate(pi_v2_v2(dm_model%n_mX, widths%n))
        pi_v2_v2 = (0.0_dp, 0.0_dp)
        allocate(pi_1_1_mat(3, 3, dm_model%n_mX, widths%n))
        pi_1_1_mat = (0.0_dp, 0.0_dp)
        allocate(pi_vi_vj(3, 3, dm_model%n_mX, widths%n))
        pi_vi_vj = (0.0_dp, 0.0_dp)

        allocate(pi_vs_vs(dm_model%n_mX, widths%n))
        pi_vs_vs = (0.0_dp, 0.0_dp)

        if ( PW_dataset%include_spin ) then

            call calc_abs_Pi(proc_id, root_process, &
                tran_form_1_spin_job, tran_form_v_spin_job, tran_form_v2_spin_job, &
                tran_form_vs_spin_job, &
                pi_v2_v2, pi_vi_vj, pi_1_1_mat, pi_vs_vs, ik_manager, job_id_to_ik, &
                PW_dataset, widths, dm_model, target_mat, io_files%out_filename, verbose = verbose)

        else

            call calc_abs_Pi(proc_id, root_process, &
                tran_form_1_no_spin_job, tran_form_v_no_spin_job, tran_form_v2_no_spin_job, &
                tran_form_vs_spin_job, &
                pi_v2_v2, pi_vi_vj, pi_1_1_mat, pi_vs_vs, ik_manager, job_id_to_ik, &
                PW_dataset, widths, dm_model, target_mat, io_files%out_filename, verbose = verbose)

        end if

        ! compute rates

        if ( verbose ) then

            print*, 'Computing absorption rates...'
            print*

        end if

        allocate(abs_rate(dm_model%n_mX, widths%n, expt%n_time))
        abs_rate = 0.0_dp

        ! parallelize velocity integral
        call v_manager%init(size(numerics%v_mesh, 1), verbose = verbose)

        do j = 1, v_manager%n_jobs_per_proc

            v_id = v_manager%job_table(proc_id + 1, j)

            if ( v_id /= 0 ) then

                v_vec = numerics%v_mesh(v_id, :)

                if ( trim(dm_model%particle_type) == 'vector' ) then

                    call calc_rate_vector(pi_1_1_mat, v_vec, &
                        dm_model, expt, widths, target_mat, numerics, &
                        abs_rate, verbose = verbose)

                else if ( trim(dm_model%particle_type) == 'ps' ) then

                    call calc_rate_ps(pi_1_1_mat, v_vec, &
                        dm_model, expt, widths, target_mat, numerics, &
                        abs_rate, verbose = verbose)

                else if ( trim(dm_model%particle_type) == 'scalar' ) then

                    if ( trim(main_control%calc_mode) == 'LO' ) then

                        call calc_rate_scalar_LO(pi_v2_v2, v_vec, &
                            dm_model, expt, widths, target_mat, numerics, &
                            abs_rate, verbose = verbose)

                    else

                        call calc_rate_scalar(pi_1_1_mat, pi_v2_v2, v_vec, &
                            dm_model, expt, widths, target_mat, numerics, &
                            abs_rate, verbose = verbose)

                    end if

                end if

            end if

        end do

        if ( verbose ) then

            print*, 'Done computing absorption rates!'
            print*

        end if

        call ik_manager%comm_abs_rate(proc_id, root_process, abs_rate, &
            verbose = verbose)

        if ( proc_id == root_process ) then

            call PW_dataset%save(io_files%out_filename, verbose = verbose)
            call widths%save(io_files%out_filename, verbose = verbose)
            call numerics%save(io_files%out_filename, verbose = verbose)

            call save_abs_rate(io_files%out_filename, abs_rate, verbose = verbose)

        end if

    end subroutine

    subroutine save_abs_rate(filename, abs_rate, verbose)
        !! Save the absorption rate data.

        implicit none
        character(len=*) :: filename
        real(dp) :: abs_rate(:, :, :)
        logical, optional :: verbose

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id

        logical :: file_exists

        integer(HSIZE_T) :: dims1(1) = [1]
        integer(HSIZE_T) :: dims2(2)

        integer :: error

        integer :: p, t
        integer :: i, fin

        if ( verbose ) then

            print*, 'Saving absorption rates...'
            print*

        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'abs_rates', group_id, error)

            do t = 1, size(abs_rate, 3)

                call h5gcreate_f(file_id,& 
                    'abs_rates'//&
                    '/t_'//trim(adjustl(int_to_str(t))),&
                    group_id, error)
                
                do p = 1, size(abs_rate, 2)

                    dims1 = [size(abs_rate, 1)]

                    call h5ltmake_dataset_double_f(file_id,&
                        'abs_rates'//&
                        '/t_'//trim(adjustl(int_to_str(t)))//&
                        '/width_'//trim(adjustl(int_to_str(p))),&
                        size(dims1), dims1,&
                        abs_rate(:, p, t), error)

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
