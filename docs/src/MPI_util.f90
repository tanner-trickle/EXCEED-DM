module MPI_util
    !! Utilities and wrapper functions for MPI routines. Generally, whenever data has to be communicated
    !! from one processor to another the `MPI_Send` and `MPI_Recv` calls are inside subroutines
    !! in this module.

    use mpi
    use prec

    use info_messages

    implicit none

    type parallel_manager_t
        !! Manages indexing the jobs that each processor needs to do and 
        !! communicates the results to the main/root processor.
        integer :: n_jobs
            !! Total number of jobs
        integer :: n_jobs_per_proc
            !! Number of transitions per processor
        integer, allocatable :: job_table(:, :)
            !! Dim : [ n_proc, n_jobs_per_proc ]
            !!
            !! `job_table(i, :)` gives the list of job ID's that processor
            !! i needs to compute for.

        contains

            procedure :: init => parallel_manager_init
            procedure :: create_job_to_2d_ID_table
            procedure :: comm_scatter_binned_rate_job_init
            procedure :: comm_abs_tran_form
            procedure :: comm_self_energies
            procedure :: comm_abs_rate
            procedure :: comm_dielectric

    end type

contains

    subroutine comm_dielectric(self, proc_id, root_process, &
            dielectric_job, dielectric, verbose)
        !! Communicate the dielectric pieces computed at each processor to the total dielectric
        !! on the main processor.
        class(parallel_manager_t) :: self
        integer :: proc_id, root_process
        complex(dp) :: dielectric_job(:, :, :, :, :)
        complex(dp) :: dielectric(:, :, :, :)
        logical, optional :: verbose

        integer :: n_proc

        integer :: status(MPI_STATUS_SIZE)
        integer :: tag = 0
        integer :: err

        integer :: i

        if ( verbose ) then
            print*, 'Communicating dielectric...'
            print*
        end if

        call MPI_COMM_SIZE(MPI_COMM_WORLD, n_proc, err)

        if ( proc_id /= root_process ) then

            ! send to main processor
            call MPI_SEND(dielectric_job, &
               size(dielectric_job), MPI_DOUBLE_COMPLEX, &
               root_process, tag, MPI_COMM_WORLD, err)

        end if

        if ( proc_id == root_process ) then

            ! add the main processors contribution
            dielectric = dielectric + sum(dielectric_job, 1)

            ! receive the other processors contributions
            do i = 1, n_proc

                if ( i - 1 /= root_process ) then

                    call MPI_RECV(dielectric_job, &
                       size(dielectric_job), MPI_DOUBLE_COMPLEX, &
                       i - 1, MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                    dielectric = dielectric + sum(dielectric_job, 1)

                end if

            end do

        end if

        if ( verbose ) then
            print*, 'Done communicating dielectric!'
            print*
        end if

    end subroutine


    subroutine comm_abs_rate(self, proc_id, root_process, &
            abs_rate, verbose)
        !! Communicate the absorption rate contribution from each processor to the main processor.

        class(parallel_manager_t) :: self
        integer :: proc_id, root_process
        real(dp) :: abs_rate(:, :, :)
        logical, optional :: verbose

        integer :: n_proc

        real(dp), allocatable :: abs_rate_buff(:, :, :)
            !! Dim : [ dm_model%n_mX, widths%n, expt%n_time ]
            !!
            !! Absorption rate, assuming the mediator-electron coupling, \( g = 1 \).
            !!
            !! Buffer variable, copy of `abs_rate`
            !!
            !! Units : None

        integer :: status(MPI_STATUS_SIZE)
        integer :: tag = 0
        integer :: err

        integer :: i

        if ( verbose ) then

            print*, 'Communicating absorption rates...'
            print*

        end if

        call MPI_COMM_SIZE(MPI_COMM_WORLD, n_proc, err)

        if ( proc_id /= root_process ) then

            call MPI_SEND(abs_rate, &
                size(abs_rate), MPI_DOUBLE, root_process, &
                tag, MPI_COMM_WORLD, err)

        end if

        if ( proc_id == root_process ) then

            allocate(abs_rate_buff(&
                size(abs_rate, 1), size(abs_rate, 2), size(abs_rate, 3)))

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

        if ( verbose ) then

            print*, 'Done communicating absorption rate data!'
            print*

        end if

    end subroutine

    subroutine comm_self_energies(self, proc_id, root_process, &
           pi_v2_v2_job, pi_1_1_mat_job, pi_vi_vj_job, &
           pi_v2_v2, pi_1_1_mat, pi_vi_vj, verbose)
        !! Communicate the self energies computed at each processor to the main processor.
        implicit none

        class(parallel_manager_t) :: self
        integer :: proc_id, root_process
        complex(dp) :: pi_v2_v2_job(:, :, :)
            !! Dim : [n_omega, n_widths, n_tran_per_proc]
            !! 
            !! self energy with two v^2 insertions
            !!
            !! Units : eV^2

        complex(dp):: pi_1_1_mat_job(:, :, :, :, :)
            !! Dim : [3, 3, n_omega, n_widths, n_tran_per_proc] 
            !!
            !! self energy, without q_vec's
            !!
            !! Units : eV^2

        complex(dp) :: pi_vi_vj_job(:, :, :, :, :)
            !! Dim : [3, 3, n_omega, n_widths, n_tran_per_proc] 
            !!
            !! Units : eV^2
        complex(dp) :: pi_v2_v2(:, :)
        complex(dp) :: pi_vi_vj(:, :, :, :)
        complex(dp) :: pi_1_1_mat(:, :, :, :)
        logical, optional :: verbose

        integer :: n_proc
        integer :: i, j

        integer :: status(MPI_STATUS_SIZE)
        integer :: tag = 0
        integer :: err

        call MPI_COMM_SIZE(MPI_COMM_WORLD, n_proc, err)

        if ( verbose ) then

            print*, 'Communicating self energies...'
            print*

        end if

        if ( proc_id /= root_process ) then

            call MPI_SEND(pi_v2_v2_job, &
               size(pi_v2_v2_job), MPI_DOUBLE_COMPLEX, root_process, tag, MPI_COMM_WORLD, err)

            call MPI_SEND(pi_1_1_mat_job, &
               size(pi_1_1_mat_job), MPI_DOUBLE_COMPLEX, root_process, tag, MPI_COMM_WORLD, err)

            call MPI_SEND(pi_vi_vj_job, &
               size(pi_vi_vj_job), MPI_DOUBLE_COMPLEX, root_process, tag, MPI_COMM_WORLD, err)

        end if

        if ( proc_id == root_process ) then

            ! add main processors contribution
            pi_v2_v2 = pi_v2_v2 + sum(pi_v2_v2_job, 3)
            pi_1_1_mat = pi_1_1_mat + sum(pi_1_1_mat_job, 5)
            pi_vi_vj = pi_vi_vj + sum(pi_vi_vj_job, 5)

            do i = 1, n_proc
                if ( (i - 1) /= root_process ) then

                    call MPI_RECV(pi_v2_v2_job, &
                       size(pi_v2_v2_job), MPI_DOUBLE_COMPLEX, i - 1, MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                    call MPI_RECV(pi_1_1_mat_job, &
                       size(pi_1_1_mat_job), MPI_DOUBLE_COMPLEX, i - 1, MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                    call MPI_RECV(pi_vi_vj_job, &
                       size(pi_vi_vj_job), MPI_DOUBLE_COMPLEX, i - 1, MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                    ! add other processors contributions
                    pi_v2_v2 = pi_v2_v2 + sum(pi_v2_v2_job, 3)
                    pi_1_1_mat = pi_1_1_mat + sum(pi_1_1_mat_job, 5)
                    pi_vi_vj = pi_vi_vj + sum(pi_vi_vj_job, 5)

                end if
            end do

        end if

        ! all processors now know all self energies

        call MPI_Bcast(pi_v2_v2, size(pi_v2_v2), &
            MPI_DOUBLE_COMPLEX, root_process, MPI_COMM_WORLD, err)

        call MPI_Bcast(pi_1_1_mat, size(pi_1_1_mat), &
            MPI_DOUBLE_COMPLEX, root_process, MPI_COMM_WORLD, err)

        call MPI_Bcast(pi_vi_vj, size(pi_vi_vj), &
            MPI_DOUBLE_COMPLEX, root_process, MPI_COMM_WORLD, err)

        if ( verbose ) then

            print*, 'Done communicating self energies!'
            print*

        end if

    end subroutine

    subroutine comm_abs_tran_form(self, proc_id, root_process, job_id_to_ik, &
            tran_form_1_no_spin_job, tran_form_1_spin_job, &
            tran_form_v_no_spin_job, tran_form_v_spin_job, &
            tran_form_v2_no_spin_job, tran_form_v2_spin_job, &
            tran_form_1_no_spin, tran_form_1_spin, &
            tran_form_v_no_spin, tran_form_v_spin, &
            tran_form_v2_no_spin, tran_form_v2_spin, numerics, verbose)
        !! Communicate the transition form factors computed at each processor to the main processor. 

        use numerics_abs

        class(parallel_manager_t) :: self
        integer :: proc_id, root_process
        integer :: job_id_to_ik(:, :)
        type(numerics_abs_t) :: numerics

        ! buffer variables
        complex(dp), allocatable:: tran_form_1_no_spin_job_buff(:, :)
            !! Dim : [n_jobs_per_proc, n_cond_max]
            !!
            !! Units : None
        complex(dp), allocatable :: tran_form_1_spin_job_buff(:, :, :, :)
            !! Dim : [n_jobs_per_proc, n_cond_max, 2, 2]
            !!
            !! Units : None
        complex(dp), allocatable :: tran_form_v_no_spin_job_buff(:, :, :)
            !! Dim : [3, n_jobs_per_proc, n_cond_max]
            !!
            !! Units : None
        complex(dp), allocatable :: tran_form_v_spin_job_buff(:, :, :, :, :)
            !! Dim : [3, n_jobs_per_proc, n_cond_max, 2, 2]
            !!
            !! Units : None
        complex(dp), allocatable :: tran_form_v2_no_spin_job_buff(:, :)
            !! Dim : [n_jobs_per_proc, n_cond_max]
            !!
            !! Units : None
        complex(dp), allocatable :: tran_form_v2_spin_job_buff(:, :, :, :)
            !! Dim : [n_jobs_per_proc, n_cond_max, 2, 2]
            !!
            !! Units : None

        complex(dp) :: tran_form_1_no_spin_job(:, :)
            !! Dim : [n_jobs_per_proc, n_cond_max]
            !!
            !! Units : None
        complex(dp) :: tran_form_1_spin_job(:, :, :, :)
            !! Dim : [n_jobs_per_proc, n_cond_max, 2, 2]
            !!
            !! Units : None
        complex(dp) :: tran_form_v_no_spin_job(:, :, :)
            !! Dim : [3, n_jobs_per_proc, n_cond_max]
            !!
            !! Units : None
        complex(dp) :: tran_form_v_spin_job(:, :, :, :, :)
            !! Dim : [3, n_jobs_per_proc, n_cond_max, 2, 2]
            !!
            !! Units : None
        complex(dp) :: tran_form_v2_no_spin_job(:, :)
            !! Dim : [n_jobs_per_proc, n_cond_max]
            !!
            !! Units : None
        complex(dp) :: tran_form_v2_spin_job(:, :, :, :)
            !! Dim : [n_jobs_per_proc, n_cond_max, 2, 2]
            !!
            !! Units : None

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
        logical, optional :: verbose

        integer :: n_proc

        integer :: init_id, k, job_id, j, i

        integer :: status(MPI_STATUS_SIZE)
        integer :: tag = 0
        integer :: err

        call MPI_COMM_SIZE(MPI_COMM_WORLD, n_proc, err)

        if ( verbose ) then

            print*, 'Communicating transition form factors for absorption calculation...'
            print*

        end if

        if ( proc_id == root_process ) then

            allocate(tran_form_1_no_spin_job_buff(self%n_jobs_per_proc, numerics%n_cond_max))
            tran_form_1_no_spin_job_buff  = (0.0_dp, 0.0_dp)
            allocate(tran_form_1_spin_job_buff(self%n_jobs_per_proc, numerics%n_cond_max, 2, 2))
            tran_form_1_spin_job_buff  = (0.0_dp, 0.0_dp)
            allocate(tran_form_v_no_spin_job_buff(3, self%n_jobs_per_proc, numerics%n_cond_max))
            tran_form_v_no_spin_job_buff  = (0.0_dp, 0.0_dp)
            allocate(tran_form_v_spin_job_buff(3, self%n_jobs_per_proc, numerics%n_cond_max, 2, 2))
            tran_form_v_spin_job_buff  = (0.0_dp, 0.0_dp)
            allocate(tran_form_v2_no_spin_job_buff(self%n_jobs_per_proc, numerics%n_cond_max))
            tran_form_v2_no_spin_job_buff  = (0.0_dp, 0.0_dp)
            allocate(tran_form_v2_spin_job_buff(self%n_jobs_per_proc, numerics%n_cond_max, 2, 2))
            tran_form_v2_spin_job_buff  = (0.0_dp, 0.0_dp)

        end if

        if ( proc_id /= root_process ) then

            call MPI_SEND(tran_form_1_no_spin_job, &
                size(tran_form_1_no_spin_job), MPI_DOUBLE_COMPLEX, root_process, &
                tag, MPI_COMM_WORLD, err)

            call MPI_SEND(tran_form_1_spin_job, &
                size(tran_form_1_spin_job), MPI_DOUBLE_COMPLEX, root_process, &
                tag, MPI_COMM_WORLD, err)

            call MPI_SEND(tran_form_v_no_spin_job, &
                size(tran_form_v_no_spin_job), MPI_DOUBLE_COMPLEX, root_process, &
                tag, MPI_COMM_WORLD, err)

            call MPI_SEND(tran_form_v_spin_job, &
                size(tran_form_v_spin_job), MPI_DOUBLE_COMPLEX, root_process, &
                tag, MPI_COMM_WORLD, err)

            call MPI_SEND(tran_form_v2_no_spin_job, &
                size(tran_form_v2_no_spin_job), MPI_DOUBLE_COMPLEX, root_process, &
                tag, MPI_COMM_WORLD, err)

            call MPI_SEND(tran_form_v2_spin_job, &
                size(tran_form_v2_spin_job), MPI_DOUBLE_COMPLEX, root_process, &
                tag, MPI_COMM_WORLD, err)

        end if

        ! receive data at main processor and put in total

        if ( proc_id == root_process ) then

            ! add main processors contribution
            do j = 1, self%n_jobs_per_proc

                job_id = self%job_table(proc_id + 1, j)

                if ( job_id /= 0 ) then

                    init_id = job_id_to_ik(job_id, 3)
                    k = job_id_to_ik(job_id, 4)

                    tran_form_1_no_spin(init_id, :, k) = tran_form_1_no_spin_job(j, :)
                    tran_form_1_spin(init_id, :, k, :, :) = tran_form_1_spin_job(j, :, :, :)
                    tran_form_v_no_spin(:, init_id, :, k) = tran_form_v_no_spin_job(:, j, :)
                    tran_form_v_spin(:, init_id, :, k, :, :) = tran_form_v_spin_job(:, j, :, :, :)
                    tran_form_v2_no_spin(init_id, :, k) = tran_form_v2_no_spin_job(j, :)
                    tran_form_v2_spin(init_id, :, k, :, :) = tran_form_v2_spin_job(j, :, :, :)

                end if

            end do

            ! receive other contributions to buffer

            do i = 1, n_proc

                if ( i - 1 /= root_process ) then

                    call MPI_RECV(tran_form_1_no_spin_job_buff, &
                        size(tran_form_1_no_spin_job_buff), MPI_DOUBLE_COMPLEX, i - 1, &
                        MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                    call MPI_RECV(tran_form_1_spin_job_buff, &
                        size(tran_form_1_spin_job_buff), MPI_DOUBLE_COMPLEX, i - 1, &
                        MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                    call MPI_RECV(tran_form_v_no_spin_job_buff, &
                        size(tran_form_v_no_spin_job_buff), MPI_DOUBLE_COMPLEX, i - 1, &
                        MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                    call MPI_RECV(tran_form_v_spin_job_buff, &
                        size(tran_form_v_spin_job_buff), MPI_DOUBLE_COMPLEX, i - 1, &
                        MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                    call MPI_RECV(tran_form_v2_no_spin_job_buff, &
                        size(tran_form_v2_no_spin_job_buff), MPI_DOUBLE_COMPLEX, i - 1, &
                        MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                    call MPI_RECV(tran_form_v2_spin_job_buff, &
                        size(tran_form_v2_spin_job_buff), MPI_DOUBLE_COMPLEX, i - 1, &
                        MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                    ! add other processors contributions
                    do j = 1, self%n_jobs_per_proc

                        job_id = self%job_table(i, j)

                        if ( job_id /= 0 ) then

                            init_id = job_id_to_ik(job_id, 3)
                            k = job_id_to_ik(job_id, 4)

                            tran_form_1_no_spin(init_id, :, k) = tran_form_1_no_spin_job_buff(j, :)
                            tran_form_1_spin(init_id, :, k, :, :) = tran_form_1_spin_job_buff(j, :, :, :)
                            tran_form_v_no_spin(:, init_id, :, k) = tran_form_v_no_spin_job_buff(:, j, :)
                            tran_form_v_spin(:, init_id, :, k, :, :) = tran_form_v_spin_job_buff(:, j, :, :, :)
                            tran_form_v2_no_spin(init_id, :, k) = tran_form_v2_no_spin_job_buff(j, :)
                            tran_form_v2_spin(init_id, :, k, :, :) = tran_form_v2_spin_job_buff(j, :, :, :)

                        end if

                    end do

                end if

            end do

        end if

        if ( verbose ) then

            print*, 'Done communicating transition form factors for absorption calculation!'
            print*

        end if

    end subroutine

    subroutine comm_scatter_binned_rate_job_init(self, proc_id, root_process, job_id_to_init_id, &
            binned_rate_job, binned_rate_init, verbose)
        !! Communicate the `binned_rate_job`'s computed at each processor to the `binned_rate_init`'s
        !! which hold the binned rates for each initial state.

        use binned_scatter_rate_type

        implicit none

        class(parallel_manager_t) :: self
        integer :: proc_id, root_process
        integer :: job_id_to_init_id(:, :)
        type(binned_scatter_rate_t) :: binned_rate_job(:)
        type(binned_scatter_rate_t) :: binned_rate_init(:)
        logical, optional :: verbose

        integer :: status(MPI_STATUS_SIZE)
        integer :: tag = 0
        integer :: err

        integer :: i, j, init_id, job_id

        integer :: n_proc

        call MPI_COMM_SIZE(MPI_COMM_WORLD, n_proc, err)

        if ( verbose ) then
            print*, 'Communicating scattering rate data...'
            print*
        end if

        ! send data to main processor
        if ( proc_id /= root_process ) then

            do j = 1, self%n_jobs_per_proc

                call MPI_SEND(binned_rate_job(j)%binned_rate, &
                   size(binned_rate_job(j)%binned_rate), &
                   MPI_DOUBLE, root_process, tag, MPI_COMM_WORLD, err)

           end do

        end if

        if ( proc_id == root_process ) then

            ! add main processors contribution to binned_rate_init
            do j = 1, self%n_jobs_per_proc

                job_id = self%job_table(proc_id + 1, j)

                if ( job_id /= 0 ) then

                    init_id = job_id_to_init_id(job_id, 3)

                    binned_rate_init(init_id)%binned_rate = binned_rate_init(init_id)%binned_rate + &
                        binned_rate_job(j)%binned_rate

                end if

            end do

            do i = 1, n_proc
                if ( (i - 1) /= root_process ) then

                    do j = 1, self%n_jobs_per_proc

                        call MPI_RECV(binned_rate_job(j)%binned_rate, &
                           size(binned_rate_job(j)%binned_rate), &
                           MPI_DOUBLE, i - 1, MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                        job_id = self%job_table(i, j)

                        if ( job_id /= 0 ) then

                            init_id = job_id_to_init_id(job_id, 3)

                            binned_rate_init(init_id)%binned_rate = binned_rate_init(init_id)%binned_rate + &
                                binned_rate_job(j)%binned_rate

                        end if

                    end do

                end if
            end do

        end if

        if ( verbose ) then
            print*, 'Done communicating scattering rate data!'
            print*
        end if

    end subroutine

    subroutine create_job_to_2d_ID_table(self, list1, list2, job_id_to_2d_ID_table, &
           verbose)
        !! Given two lists, {list1, list2} create a single list which has
        !! n[list1]*n[list_2] elements, each of which has 4 elements. The first two elements
        !! are a pair from {list1, list2}, the second two index the value in the first pair, i.e.
        !! 
        !! `job_id_to_2d_ID_table(1, 1) = list1(1)`
        !!
        !! `job_id_to_2d_ID_table(1, 2) = list2(7)`
        !!
        !! `job_id_to_2d_ID_table(1, 3) = 1`
        !!
        !! `job_id_to_2d_ID_table(1, 4) = 7`

        implicit none

        class(parallel_manager_t) :: self

        integer :: list1(:)
        integer :: list2(:)
        integer :: job_id_to_2d_ID_table(:, :)
            !! Dim : [n_jobs, 4]
        logical, optional :: verbose

        integer :: id 
        integer :: i1, i2

        id = 0
        do i1 = 1, size(list1)
            do i2 = 1, size(list2)

                id = id + 1

                job_id_to_2d_ID_table(id, 1) = list1(i1)
                job_id_to_2d_ID_table(id, 2) = list2(i2)
                job_id_to_2d_ID_table(id, 3) = i1
                job_id_to_2d_ID_table(id, 4) = i2

            end do
        end do 

    end subroutine

    subroutine parallel_manager_init(self, n_jobs, verbose)
        !! Initializes a `parallel_manager` instance.

        implicit none
        class(parallel_manager_t) :: self
        integer :: n_jobs
            !! Number of initial states
        logical, optional :: verbose

        integer :: n_proc
            !! Number of processors
        integer :: i, j, job_id, id, f
        integer :: err

        call MPI_COMM_SIZE(MPI_COMM_WORLD, n_proc, err)

        self%n_jobs = n_jobs

        if ( mod(self%n_jobs, n_proc) == 0 ) then 
            self%n_jobs_per_proc = self%n_jobs/n_proc
        else
            self%n_jobs_per_proc = self%n_jobs/n_proc + 1
        end if 

        allocate(self%job_table(n_proc, self%n_jobs_per_proc))

        job_id = 0
        
        do j = 1, self%n_jobs_per_proc
            do i = 1, n_proc
                
                job_id = job_id + 1

                if (job_id > self%n_jobs) then 
                    self%job_table(i, j) = 0
                else
                    self%job_table(i, j) = job_id
                end if 

            end do
        end do  
        
        if ( verbose ) then             

            call print_section_seperator()
            print*
            
            if ( mod(self%n_jobs, n_proc) == 0 ) then 
        
                print*, '    Equal processor load.'
                print*
            
            else if ( self%n_jobs_per_proc == 1 ) then 

                print*, '    Number of processors is greater than the number of jobs.'//&
                        ' Consider lowering the number of processors.'
                print*
                print*, '    Number of jobs = ', self%n_jobs
                print*

            else

                print*, '    Unequal processor load. Some processors will be given null jobs.'
                print*
                print*, '    Number of jobs = ', self%n_jobs
                print*

            end if 

            print*, '    Number of jobs per processor = ', self%n_jobs_per_proc
            print*   
            call print_section_seperator()
            print*
            
        end if

    end subroutine

end module
