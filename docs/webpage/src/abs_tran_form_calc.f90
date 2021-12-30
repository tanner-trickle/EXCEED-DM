module abs_tran_form_calc
    !! Computes the transition form factors, \( f \), for absorption rate calculation.

    use prec
    use constants

    use MPI_util

    use PW_dataset_type

    use numerics_abs

    implicit none

    ! Compute different transition form factors for spin-independent
    ! and spin-dependent wave functions 
    interface calc_tran_form_1
        module procedure calc_tran_form_1_no_spin
        module procedure calc_tran_form_1_spin
    end interface

    interface calc_tran_form_v
        module procedure calc_tran_form_v_no_spin
        module procedure calc_tran_form_v_spin
    end interface

    interface calc_tran_form_v2
        module procedure calc_tran_form_v2_no_spin
        module procedure calc_tran_form_v2_spin
    end interface
    
contains

    subroutine calc_abs_tran_form(proc_id, root_process, &
                                  tran_form_1_no_spin_job, tran_form_1_spin_job, &
                                  tran_form_v_no_spin_job, tran_form_v_spin_job, &
                                  tran_form_v2_no_spin_job, tran_form_v2_spin_job, &
                                  tran_form_vs_spin_job, &
                                  PW_dataset, ik_manager, job_id_to_ik, numerics, &
                                  filename, &
                                  verbose)
        !! Computes all of the transition form factors, 
        !! \( f_{i, i', \mathbf{k}}, \mathbf{f}_{i, i', \mathbf{k}}, \widetilde{f}_{i, i', \mathbf{k}} \) 
        !! (with extra \( s, s' \) spin indicies if the wave functions are spin dependent) for the 
        !! DM absorption calculation.

        implicit none

        integer :: proc_id
        integer :: root_process

        complex(dp) :: tran_form_1_no_spin_job(:, :) 
        complex(dp) :: tran_form_1_spin_job(:, :, :, :) 
        complex(dp) :: tran_form_v_no_spin_job(:, :, :)
        complex(dp) :: tran_form_v_spin_job(:, :, :, :, :)
        complex(dp) :: tran_form_v2_no_spin_job(:, :)
        complex(dp) :: tran_form_v2_spin_job(:, :, :, :)
        complex(dp) :: tran_form_vs_spin_job(:, :)

        type(PW_dataset_t) :: PW_dataset
        type(parallel_manager_t) :: ik_manager
        integer :: job_id_to_ik(:, :)
        type(numerics_abs_t) :: numerics
        character(len=*) :: filename
        logical, optional :: verbose

        complex(dp), allocatable :: tran_form_1_no_spin(:, :, :)
            !! Dim : [n_init, n_fin, n_k]
            !!
            !! Units : None
        complex(dp), allocatable :: tran_form_1_spin(:, :, :, :, :)
            !! Dim : [n_init, n_fin, n_k, 2, 2]
            !!
            !! Units : None
        complex(dp), allocatable :: tran_form_v_no_spin(:, :, :, :)
            !! Dim : [3, n_init, n_fin, n_k]
            !!
            !! Units : None
        complex(dp), allocatable :: tran_form_v_spin(:, :, :, :, :, :)
            !! Dim : [3, n_init, n_fin, n_k, 2, 2]
            !!
            !! Units : None
        complex(dp), allocatable :: tran_form_v2_no_spin(:, :, :)
            !! Dim : [n_init, n_fin, n_k]
            !!
            !! Units : None
        complex(dp), allocatable :: tran_form_v2_spin(:, :, :, :, :)
            !! Dim : [n_init, n_fin, n_k, 2, 2]
            !!
            !! Units : None
        complex(dp), allocatable :: tran_form_vs_spin(:, :, :)
            !! Dim : [n_init, n_fin, n_k]
            !!
            !! Units : None

        integer :: j, job_id, f
        integer :: val_id, k, cond_id

        complex(dp), allocatable :: wfc_FT_ik(:)
        complex(dp), allocatable :: wfc_FT_fk(:)

        complex(dp), allocatable :: wfc_FT_iks(:, :)
        complex(dp), allocatable :: wfc_FT_fks(:, :)

        if ( verbose ) then

            print*, 'Computing transition form factors...'
            print*

        end if

        ! allocate wave function data
        if ( PW_dataset%include_spin ) then

            allocate(wfc_FT_iks(PW_dataset%n_G, 2))
            allocate(wfc_FT_fks(PW_dataset%n_G, 2))

        else

            allocate(wfc_FT_ik(PW_dataset%n_G))
            allocate(wfc_FT_fk(PW_dataset%n_G))

        end if

        ! compute processor specific variables
        
        do j = 1, ik_manager%n_jobs_per_proc

            job_id = ik_manager%job_table(proc_id + 1, j)

            if ( job_id /= 0 ) then

                val_id = job_id_to_ik(job_id, 1)
                k = job_id_to_ik(job_id, 2)

                if ( PW_dataset%include_spin ) then

                    ! load initial wave function
                    call PW_dataset%load_wfc_FT_ik_spin(val_id, k, wfc_FT_iks)

                    do f = 1, numerics%n_cond_max

                        cond_id = PW_dataset%n_val + f

                        ! load final wave function
                        call PW_dataset%load_wfc_FT_ik_spin(cond_id, k, wfc_FT_fks)

                        call calc_tran_form_1(tran_form_1_spin_job(j, f, :, :), &
                            PW_dataset, val_id, cond_id, k, wfc_FT_iks, wfc_FT_fks)

                        call calc_tran_form_v(tran_form_v_spin_job(:, j, f, :, :), &
                            PW_dataset, val_id, cond_id, k, wfc_FT_iks, wfc_FT_fks)

                        call calc_tran_form_v2(tran_form_v2_spin_job(j, f, :, :), &
                            PW_dataset, val_id, cond_id, k, wfc_FT_iks, wfc_FT_fks)

                        call calc_tran_form_vs_spin(tran_form_vs_spin_job(j, f), &
                            PW_dataset, val_id, cond_id, k, wfc_FT_iks, wfc_FT_fks)

                    end do

                else

                    ! load initial wave function
                    call PW_dataset%load_wfc_FT_ik_no_spin(val_id, k, wfc_FT_ik)

                    do f = 1, numerics%n_cond_max

                        cond_id = PW_dataset%n_val + f

                        call PW_dataset%load_wfc_FT_ik_no_spin(cond_id, k, wfc_FT_fk)

                        call calc_tran_form_1(tran_form_1_no_spin_job(j, f), &
                            PW_dataset, val_id, cond_id, k, wfc_FT_ik, wfc_FT_fk)

                        call calc_tran_form_v(tran_form_v_no_spin_job(:, j, f), &
                            PW_dataset, val_id, cond_id, k, wfc_FT_ik, wfc_FT_fk)

                        call calc_tran_form_v2(tran_form_v2_no_spin_job(j, f), &
                            PW_dataset, val_id, cond_id, k, wfc_FT_ik, wfc_FT_fk)

                    end do

                end if

            end if

        end do

        if ( verbose ) then

            print*, 'Done computing transition form factors!'
            print*

        end if

        if ( numerics%save_tran_form ) then

            if ( proc_id == root_process ) then

                allocate(tran_form_1_no_spin(numerics%n_val_max, numerics%n_cond_max, PW_dataset%n_k))
                tran_form_1_no_spin     = (0.0_dp, 0.0_dp)
                allocate(tran_form_1_spin(numerics%n_val_max, numerics%n_cond_max, PW_dataset%n_k, 2, 2))
                tran_form_1_spin     = (0.0_dp, 0.0_dp)
                allocate(tran_form_v_no_spin(3, numerics%n_val_max, numerics%n_cond_max, PW_dataset%n_k))
                tran_form_v_no_spin  = (0.0_dp, 0.0_dp)
                allocate(tran_form_v_spin(3,numerics%n_val_max, numerics%n_cond_max, PW_dataset%n_k, 2, 2))
                tran_form_v_spin     = (0.0_dp, 0.0_dp)
                allocate(tran_form_v2_no_spin(numerics%n_val_max, numerics%n_cond_max, PW_dataset%n_k))
                tran_form_v2_no_spin = (0.0_dp, 0.0_dp)
                allocate(tran_form_v2_spin(numerics%n_val_max, numerics%n_cond_max, PW_dataset%n_k, 2, 2))
                tran_form_v2_spin    = (0.0_dp, 0.0_dp)
                allocate(tran_form_vs_spin(numerics%n_val_max, numerics%n_cond_max, PW_dataset%n_k))
                tran_form_vs_spin     = (0.0_dp, 0.0_dp)

            end if

            ! communicate the jobs to the main processor
            call ik_manager%comm_abs_tran_form(proc_id, root_process, job_id_to_ik, &
                tran_form_1_no_spin_job, tran_form_1_spin_job, &
                tran_form_v_no_spin_job, tran_form_v_spin_job, &
                tran_form_v2_no_spin_job, tran_form_v2_spin_job, &
                tran_form_vs_spin_job, &
                tran_form_1_no_spin, tran_form_1_spin, &
                tran_form_v_no_spin, tran_form_v_spin, &
                tran_form_v2_no_spin, tran_form_v2_spin, &
                tran_form_vs_spin, numerics, verbose = verbose)

            ! save data
            if ( proc_id == root_process ) then

                call save_tran_form(filename, &
                       tran_form_1_no_spin, tran_form_1_spin, &
                       tran_form_v_no_spin, tran_form_v_spin, &
                       tran_form_v2_no_spin, tran_form_v2_spin, &
                       PW_dataset%include_spin, verbose = verbose)

            end if

        end if

    end subroutine

    subroutine calc_tran_form_1_no_spin(tran_form_1, &
            PW_dataset, val_id, cond_id, k, wfc_FT_ik, wfc_FT_fk)
        !! Computes the spin independent scalar transition form factor:
        !!
        !! $$\begin{align}
        !!      f_{i, i', \mathbf{k}} = \sum_\mathbf{G} \widetilde{u}_{i' \mathbf{k} \mathbf{G}}^* \widetilde{u}_{i \mathbf{k}
        !! \mathbf{G}} = \delta_{i i'} \nonumber
        !! \end{align}$$ 
        !!
        !! for a given i, i', k.
        !!
        !! Have this function here for a consistent interface. This vanished by orthogonality of the
        !! wave functions.
        complex(dp) :: tran_form_1
        type(PW_dataset_t) :: PW_dataset
        integer :: val_id
        integer :: cond_id
        integer :: k
        complex(dp) :: wfc_FT_ik(:)
        complex(dp) :: wfc_FT_fk(:)

        if ( val_id /= cond_id ) then
            tran_form_1 = (0.0_dp, 0.0_dp)
        end if

    end subroutine

    subroutine calc_tran_form_1_spin(tran_form_1, &
            PW_dataset, val_id, cond_id, k, wfc_FT_ik, wfc_FT_fk)
        !! Computes the spin dependent scalar transition form factor:
        !!
        !! $$\begin{align}
        !!      f_{i, i', \mathbf{k}}^{ss'} = \sum_\mathbf{G} {\widetilde{u}_{i' \mathbf{k} \mathbf{G}}^{s'}}^* \widetilde{u}_{i \mathbf{k}
        !! \mathbf{G}}^s \nonumber
        !! \end{align}$$ 
        !!
        !! for a given i, i', k.
        !!
        !! Note:
        !! $$\begin{align}
        !!      \sum_{s} f_{i, i', \mathbf{k}}^{s s} = \delta_{i i'}
        !! \end{align}$$ 
        complex(dp) :: tran_form_1(:, :)
        type(PW_dataset_t) :: PW_dataset
        integer :: val_id
        integer :: cond_id
        integer :: k
        complex(dp) :: wfc_FT_ik(:, :)
        complex(dp) :: wfc_FT_fk(:, :)

        integer :: s, sp

        tran_form_1 = (0.0_dp, 0.0_dp)

        if ( val_id /= cond_id ) then

            do s = 1, 2
                do sp = 1, 2

                    tran_form_1(s, sp) = tran_form_1(s, sp) + sum(conjg(wfc_FT_fk(:, sp))*wfc_FT_ik(:, s))

                end do
            end do

        else

            tran_form_1 = (0.0_dp, 0.0_dp)

        end if

    end subroutine

    subroutine calc_tran_form_v_no_spin(tran_form_v, &
            PW_dataset, val_id, cond_id, k, wfc_FT_ik, wfc_FT_fk)
        !! Computes the spin independent vector transition form factor:
        !!
        !! $$\begin{align*}
        !!      \mathbf{f}_{i, i', \mathbf{k}} = \frac{1}{m_e} \sum_\mathbf{G} (\mathbf{k} + \mathbf{G}) \widetilde{u}_{i' \mathbf{k} \mathbf{G}}^* \widetilde{u}_{i \mathbf{k}
        !! \mathbf{G}}
        !! \end{align*}$$ 
        !!
        !! for a given i, i', k.
        complex(dp) :: tran_form_v(:)
        type(PW_dataset_t) :: PW_dataset
        integer :: val_id
        integer :: cond_id
        integer :: k
        complex(dp) :: wfc_FT_ik(:)
        complex(dp) :: wfc_FT_fk(:)

        integer :: g

        real(dp) :: p_vec(3)

        do g = 1, PW_dataset%n_G

            p_vec = PW_dataset%k_grid_xyz(k, :) + PW_dataset%G_grid_xyz(g, :)

            tran_form_v = tran_form_v + &
                (p_vec/m_elec)*conjg(wfc_FT_fk(g))*wfc_FT_ik(g)

        end do

    end subroutine

    subroutine calc_tran_form_v_spin(tran_form_v, &
            PW_dataset, val_id, cond_id, k, wfc_FT_ik, wfc_FT_fk)
        !! Computes the spin dependent vector transition form factor:
        !!
        !! $$\begin{align*}
        !!      \mathbf{f}_{i, i', \mathbf{k}, s, s'} = \frac{1}{m_e} \sum_\mathbf{G} (\mathbf{k} + \mathbf{G}) 
        !! \widetilde{u}_{i' \mathbf{k} \mathbf{G} s'}^* \widetilde{u}_{i \mathbf{k} \mathbf{G} s}
        !! \end{align*}$$ 
        !!
        !! for a given i, i', k.
        complex(dp) :: tran_form_v(:, :, :)
        type(PW_dataset_t) :: PW_dataset
        integer :: val_id
        integer :: cond_id
        integer :: k
        complex(dp) :: wfc_FT_ik(:, :)
        complex(dp) :: wfc_FT_fk(:, :)

        integer :: g, s, sp

        real(dp) :: p_vec(3)
        
        do s = 1, 2
            do sp = 1, 2

                do g = 1, PW_dataset%n_G

                    p_vec = PW_dataset%k_grid_xyz(k, :) + PW_dataset%G_grid_xyz(g, :)

                    tran_form_v(:, s, sp) = tran_form_v(:, s, sp) + &
                        (p_vec/m_elec)*conjg(wfc_FT_fk(g, sp))*wfc_FT_ik(g, s)

                end do

            end do

        end do

    end subroutine

    subroutine calc_tran_form_vs_spin(tran_form_vs, &
            PW_dataset, val_id, cond_id, k, wfc_FT_ik, wfc_FT_fk)
        !* Computes the spin dependent \( \mathbf{v} \cdot \mathbf{\sigma} \) transition form factor:
        !
        ! $$
        ! \begin{align*}
        !      f_{i, i', \mathbf{k}} = \frac{1}{m_e} \sum_\mathbf{G} \mathbf{\widetilde{u}}_{i' \mathbf{k} \mathbf{G}}^* \cdot (
        !      (\mathbf{k} + \mathbf{G}) \cdot \sigma ) \cdot \mathbf{\widetilde{u}}_{i \mathbf{k} \mathbf{G}} 
        !  \end{align*}$$ 
        !
        ! for a given i, i', k.

        complex(dp) :: tran_form_vs
        type(PW_dataset_t) :: PW_dataset
        integer :: val_id
        integer :: cond_id
        integer :: k
        complex(dp) :: wfc_FT_ik(:, :)
        complex(dp) :: wfc_FT_fk(:, :)

        integer :: g, s, sp

        real(dp) :: p_vec(3)

        complex(dp) :: p_dot_sigma(2, 2)

        do g = 1, PW_dataset%n_G

            p_vec = PW_dataset%k_grid_xyz(k, :) + PW_dataset%G_grid_xyz(g, :)

            p_dot_sigma(1, 1) = p_vec(3)
            p_dot_sigma(1, 2) = p_vec(1) - ii*p_vec(2)
            p_dot_sigma(2, 1) = p_vec(1) + ii*p_vec(2)
            p_dot_sigma(2, 2) = -p_vec(3)

            tran_form_vs = tran_form_vs + &
                (1.0_dp/m_elec)*&
                p_dot_sigma(1, 1)*conjg(wfc_FT_fk(g, 1))*wfc_FT_ik(g, 1)

            tran_form_vs = tran_form_vs + &
                (1.0_dp/m_elec)*&
                p_dot_sigma(1, 2)*conjg(wfc_FT_fk(g, 1))*wfc_FT_ik(g, 2)

            tran_form_vs = tran_form_vs + &
                (1.0_dp/m_elec)*&
                p_dot_sigma(2, 1)*conjg(wfc_FT_fk(g, 2))*wfc_FT_ik(g, 1)

            tran_form_vs = tran_form_vs + &
                (1.0_dp/m_elec)*&
                p_dot_sigma(2, 2)*conjg(wfc_FT_fk(g, 2))*wfc_FT_ik(g, 2)

        end do

    end subroutine

    subroutine calc_tran_form_v2_no_spin(tran_form_v2, &
            PW_dataset, val_id, cond_id, k, wfc_FT_ik, wfc_FT_fk)
        !! Computes the spin independent \( v^2 \) transition form factor:
        !!
        !! $$\begin{align}
        !!      \widetilde{f}_{i, i', \mathbf{k}} = \frac{1}{m_e^2} 
        !!  \sum_\mathbf{G} ( \mathbf{k} + \mathbf{G} )^2 \widetilde{u}_{i' \mathbf{k} \mathbf{G}}^* 
        !! \widetilde{u}_{i \mathbf{k} \mathbf{G}}
        !! \end{align}$$ 
        !!
        !! for a given i, i', k.
        complex(dp) :: tran_form_v2
        type(PW_dataset_t) :: PW_dataset
        integer :: val_id
        integer :: cond_id
        integer :: k
        complex(dp) :: wfc_FT_ik(:)
        complex(dp) :: wfc_FT_fk(:)

        integer :: g

        real(dp) :: p_vec(3)
        real(dp) :: p_sq

        do g = 1, PW_dataset%n_G

            p_vec = PW_dataset%k_grid_xyz(k, :) + PW_dataset%G_grid_xyz(g, :)

            p_sq = norm2(p_vec)**2

            tran_form_v2 = tran_form_v2 + &
                (p_sq/m_elec**2)*conjg(wfc_FT_fk(g))*wfc_FT_ik(g)

        end do

    end subroutine

    subroutine calc_tran_form_v2_spin(tran_form_v2, &
            PW_dataset, val_id, cond_id, k, wfc_FT_ik, wfc_FT_fk)
        !! Computes the spin dependent \( v^2 \) transition form factor:
        !!
        !! $$\begin{align}
        !!      \widetilde{f}_{i, i', \mathbf{k}, s, s'} = \frac{1}{m_e^2} 
        !!  \sum_\mathbf{G} ( \mathbf{k} + \mathbf{G} )^2 \widetilde{u}_{i' \mathbf{k} \mathbf{G} s'}^* 
        !! \widetilde{u}_{i \mathbf{k} \mathbf{G} s}
        !! \end{align}$$ 
        !!
        !! for a given i, i', k.
        complex(dp) :: tran_form_v2(:, :)
        type(PW_dataset_t) :: PW_dataset
        integer :: val_id
        integer :: cond_id
        integer :: k
        complex(dp) :: wfc_FT_ik(:, :)
        complex(dp) :: wfc_FT_fk(:, :)

        integer :: g, s, sp

        real(dp) :: p_vec(3)
        real(dp) :: p_sq

        do s = 1, 2
            do sp = 1, 2

                do g = 1, PW_dataset%n_G

                    p_vec = PW_dataset%k_grid_xyz(k, :) + PW_dataset%G_grid_xyz(g, :)

                    p_sq = norm2(p_vec)**2

                    tran_form_v2(s, sp) = tran_form_v2(s, sp) + &
                        (p_sq/m_elec**2)*conjg(wfc_FT_fk(g, sp))*wfc_FT_ik(g, s)

                end do

            end do
        end do

    end subroutine

    subroutine save_tran_form(filename, &
           tran_form_1_no_spin, tran_form_1_spin, &
           tran_form_v_no_spin, tran_form_v_spin, &
           tran_form_v2_no_spin, tran_form_v2_spin, &
           include_spin, verbose)
        !! Saves the transition form factors, \( f_{i, i', \mathbf{k}}, 
        !! \mathbf{f}_{i, i', \mathbf{k}}, \widetilde{f}_{i, i', k}\), 
        !! with \( s, s' \) indicies if the wave functions are spin dependent.

        use info_messages

        implicit none

        character(len=*) :: filename
        logical, optional :: verbose
        logical :: include_spin

        complex(dp) :: tran_form_1_no_spin(:, :, :)
            !! Dim : [n_init, n_fin, n_k]
            !!
            !! Units : None
        complex(dp) :: tran_form_1_spin(:, :, :, :, :)
            !! Dim : [n_init, n_fin, n_k, 2, 2]
            !!
            !! Units : None
        complex(dp) :: tran_form_v_no_spin(:, :, :, :)
            !! Dim : [3, n_init, n_fin, n_k]
            !!
            !! Units : None
        complex(dp) :: tran_form_v_spin(:, :, :, :, :, :)
            !! Dim : [3, n_init, n_fin, n_k, 2, 2]
            !!
            !! Units : None
        complex(dp) :: tran_form_v2_no_spin(:, :, :)
            !! Dim : [n_init, n_fin, n_k]
            !!
            !! Units : None
        complex(dp) :: tran_form_v2_spin(:, :, :, :, :)
            !! Dim : [n_init, n_fin, n_k, 2, 2]
            !!
            !! Units : None

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id

        logical :: file_exists

        integer(HSIZE_T) :: dims1(1) = [1]
        integer(HSIZE_T) :: dims2(2)
        integer(HSIZE_T) :: dims3(3)
        integer(HSIZE_T) :: dims4(4)

        integer :: error

        integer :: m, f, t
        integer :: i, fin

        if ( verbose ) then

            print*, 'Saving absorption transition form factors...'
            print*

        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'abs_tran_form', group_id, error)

            do i = 1, size(tran_form_1_no_spin, 1)

                call h5gcreate_f(file_id,& 
                    'abs_tran_form'//&
                    '/init_'//trim(adjustl(int_to_str(i))),&
                    group_id, error)

                do fin = 1, size(tran_form_1_no_spin, 2)

                    call h5gcreate_f(file_id,& 
                        'abs_tran_form'//&
                        '/init_'//trim(adjustl(int_to_str(i)))//&
                        '/fin_'//trim(adjustl(int_to_str(fin))),&
                        group_id, error)

                    if ( include_spin ) then

                        dims3 = [size(tran_form_1_no_spin, 3), 2, 2]

                        call h5ltmake_dataset_double_f(file_id,&
                            'abs_tran_form'//&
                            '/init_'//trim(adjustl(int_to_str(i)))//&
                            '/fin_'//trim(adjustl(int_to_str(fin)))//&
                            '/tran_form_1_spin_r', &
                            size(dims3), dims3,&
                            real(tran_form_1_spin(i, fin, :, :, :)), error)

                        call h5ltmake_dataset_double_f(file_id,&
                            'abs_tran_form'//&
                            '/init_'//trim(adjustl(int_to_str(i)))//&
                            '/fin_'//trim(adjustl(int_to_str(fin)))//&
                            '/tran_form_1_spin_c', &
                            size(dims3), dims3,&
                            aimag(tran_form_1_spin(i, fin, :, :, :)), error)

                        call h5ltmake_dataset_double_f(file_id,&
                            'abs_tran_form'//&
                            '/init_'//trim(adjustl(int_to_str(i)))//&
                            '/fin_'//trim(adjustl(int_to_str(fin)))//&
                            '/tran_form_v2_spin_r', &
                            size(dims3), dims3,&
                            real(tran_form_v2_spin(i, fin, :, :, :)), error)

                        call h5ltmake_dataset_double_f(file_id,&
                            'abs_tran_form'//&
                            '/init_'//trim(adjustl(int_to_str(i)))//&
                            '/fin_'//trim(adjustl(int_to_str(fin)))//&
                            '/tran_form_v2_spin_c', &
                            size(dims3), dims3,&
                            aimag(tran_form_v2_spin(i, fin, :, :, :)), error)

                        dims4 = [3, size(tran_form_1_no_spin, 3), 2, 2]

                        call h5ltmake_dataset_double_f(file_id,&
                            'abs_tran_form'//&
                            '/init_'//trim(adjustl(int_to_str(i)))//&
                            '/fin_'//trim(adjustl(int_to_str(fin)))//&
                            '/tran_form_v_spin_r', &
                            size(dims4), dims4,&
                            real(tran_form_v_spin(:, i, fin, :, :, :)), error)

                        call h5ltmake_dataset_double_f(file_id,&
                            'abs_tran_form'//&
                            '/init_'//trim(adjustl(int_to_str(i)))//&
                            '/fin_'//trim(adjustl(int_to_str(fin)))//&
                            '/tran_form_v_spin_c', &
                            size(dims4), dims4,&
                            aimag(tran_form_v_spin(:, i, fin, :, :, :)), error)

                    else

                        dims1 = [size(tran_form_1_no_spin, 3)]

                        call h5ltmake_dataset_double_f(file_id,&
                            'abs_tran_form'//&
                            '/init_'//trim(adjustl(int_to_str(i)))//&
                            '/fin_'//trim(adjustl(int_to_str(fin)))//&
                            '/tran_form_v2_r', &
                            size(dims1), dims1,&
                            real(tran_form_v2_no_spin(i, fin, :)), error)

                        call h5ltmake_dataset_double_f(file_id,&
                            'abs_tran_form'//&
                            '/init_'//trim(adjustl(int_to_str(i)))//&
                            '/fin_'//trim(adjustl(int_to_str(fin)))//&
                            '/tran_form_v2_c', &
                            size(dims1), dims1,&
                            aimag(tran_form_v2_no_spin(i, fin, :)), error)

                        dims2 = [3, size(tran_form_1_no_spin, 3)]

                        call h5ltmake_dataset_double_f(file_id,&
                            'abs_tran_form'//&
                            '/init_'//trim(adjustl(int_to_str(i)))//&
                            '/fin_'//trim(adjustl(int_to_str(fin)))//&
                            '/tran_form_v_r', &
                            size(dims2), dims2,&
                            real(tran_form_v_no_spin(:, i, fin, :)), error)

                        call h5ltmake_dataset_double_f(file_id,&
                            'abs_tran_form'//&
                            '/init_'//trim(adjustl(int_to_str(i)))//&
                            '/fin_'//trim(adjustl(int_to_str(fin)))//&
                            '/tran_form_v_c', &
                            size(dims2), dims2,&
                            aimag(tran_form_v_no_spin(:, i, fin, :)), error)

                    end if

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
