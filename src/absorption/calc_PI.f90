module calc_PI
    !! Computes self-energies needed for absorption calculations

    use mpi
    use hdf5
    use h5lt

    use prec
    use constants
    use material_input
    use absorption_input
    use DFT_parameters
    use transition

    implicit none

    complex(dp), allocatable :: pi_v2_v2(:, :)
        !! Dim : [n_omega, n_widths]
        !! 
        !! self energy with two v^2 insertions
        !!
        !! Units : eV^2

    complex(dp), allocatable :: pi_v2_v2_t(:, :, :)
        !! Dim : [n_omega, n_widths, n_tran_per_proc]
        !! 
        !! self energy with two v^2 insertions
        !!
        !! Units : eV^2

    complex(dp), allocatable :: pi_1_1_mat(:, :, :, :)
        !! Dim : [3, 3, n_omega, n_widths] 
        !!
        !! self energy, without q_vec's
        !!
        !! Pi_11 = (q/m_elec) . Pi_11_mat . (q/m_elec)
        !!
        !! Units : eV^2

    complex(dp), allocatable :: pi_1_1_mat_t(:, :, :, :, :)
        !! Dim : [3, 3, n_omega, n_widths, n_tran_per_proc] 
        !!
        !! self energy, without q_vec's
        !!
        !! Units : eV^2

    complex(dp), allocatable :: pi_vi_vj(:, :, :, :)
        !! Dim : [3, 3, n_omega, n_widths] 
        !!
        !! Units : eV^2

    complex(dp), allocatable :: pi_vi_vj_t(:, :, :, :, :)
        !! Dim : [3, 3, n_omega, n_widths, n_tran_per_proc] 
        !!
        !! Units : eV^2

    interface compute_self_energies
        module procedure compute_self_energies_no_spin
        module procedure compute_self_energies_spin
    end interface

    interface calc_pi_v2_v2
        module procedure calc_pi_v2_v2_no_spin
        module procedure calc_pi_v2_v2_spin
    end interface

    interface calc_pi_1_1_mat
        module procedure calc_pi_1_1_mat_no_spin
        module procedure calc_pi_1_1_mat_spin
    end interface

    interface calc_pi_vi_vj
        module procedure calc_pi_vi_vj_no_spin
        module procedure calc_pi_vi_vj_spin
    end interface

contains

    subroutine compute_self_energies_setup(proc_id, root_process, verbose)

        implicit none

        integer :: proc_id, root_process

        logical, optional :: verbose

        if ( verbose ) then

            print*, 'Calculating self energies...'
            print*

        end if

        allocate(pi_v2_v2(n_omega, n_widths))
        pi_v2_v2 = (0.0_dp, 0.0_dp)

        allocate(pi_v2_v2_t(n_omega, n_widths, n_tran_per_proc))
        pi_v2_v2_t = (0.0_dp, 0.0_dp)

        allocate(pi_1_1_mat(3, 3, n_omega, n_widths))
        pi_1_1_mat = (0.0_dp, 0.0_dp)

        allocate(pi_1_1_mat_t(3, 3, n_omega, n_widths, n_tran_per_proc))
        pi_1_1_mat_t = (0.0_dp, 0.0_dp)

        allocate(pi_vi_vj(3, 3, n_omega, n_widths))
        pi_vi_vj = (0.0_dp, 0.0_dp)

        allocate(pi_vi_vj_t(3, 3, n_omega, n_widths, n_tran_per_proc))
        pi_vi_vj_t = (0.0_dp, 0.0_dp)

    end subroutine

    subroutine bcast_self_energies(proc_id, root_process, verbose)
        !! After the self energies have been summed together on the main
        !! processor, send the summed Pi's back out to each processor
        !! so that the rate calculation can be parallelized over v

        implicit none

        integer :: proc_id, root_process
        logical, optional :: verbose
        integer :: err

        call MPI_Bcast(pi_v2_v2, size(pi_v2_v2), &
            MPI_DOUBLE_COMPLEX, root_process, MPI_COMM_WORLD, err)

        call MPI_Bcast(pi_1_1_mat, size(pi_1_1_mat), &
            MPI_DOUBLE_COMPLEX, root_process, MPI_COMM_WORLD, err)

        call MPI_Bcast(pi_vi_vj, size(pi_vi_vj), &
            MPI_DOUBLE_COMPLEX, root_process, MPI_COMM_WORLD, err)

    end subroutine

    subroutine comm_self_energies(proc_id, root_process, n_proc, verbose)
        !! Communicate the self energies that each processor is computing

        implicit none

        integer :: proc_id, root_process, n_proc

        logical, optional :: verbose

        integer :: status(MPI_STATUS_SIZE)
        integer :: tag = 0
        integer :: err

        integer :: i

        ! sum together again 
        if ( proc_id .ne. root_process ) then

            call MPI_SEND(pi_v2_v2_t, &
               size(pi_v2_v2_t), MPI_DOUBLE_COMPLEX, root_process, tag, MPI_COMM_WORLD, err)

            call MPI_SEND(pi_1_1_mat_t, &
               size(pi_1_1_mat_t), MPI_DOUBLE_COMPLEX, root_process, tag, MPI_COMM_WORLD, err)

            call MPI_SEND(pi_vi_vj_t, &
               size(pi_vi_vj_t), MPI_DOUBLE_COMPLEX, root_process, tag, MPI_COMM_WORLD, err)

        end if

        if ( proc_id .eq. root_process ) then

            ! add main processors contribution
            pi_v2_v2 = pi_v2_v2 + sum(pi_v2_v2_t, 3)
            pi_1_1_mat = pi_1_1_mat + sum(pi_1_1_mat_t, 5)
            pi_vi_vj = pi_vi_vj + sum(pi_vi_vj_t, 5)

            do i = 1, n_proc
                if ( (i - 1) .ne. root_process ) then

                    call MPI_RECV(pi_v2_v2_t, &
                       size(pi_v2_v2_t), MPI_DOUBLE_COMPLEX, i - 1, MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                    call MPI_RECV(pi_1_1_mat_t, &
                       size(pi_1_1_mat_t), MPI_DOUBLE_COMPLEX, i - 1, MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                    call MPI_RECV(pi_vi_vj_t, &
                       size(pi_vi_vj_t), MPI_DOUBLE_COMPLEX, i - 1, MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                    ! add other processors contributions
                    pi_v2_v2 = pi_v2_v2 + sum(pi_v2_v2_t, 3)
                    pi_1_1_mat = pi_1_1_mat + sum(pi_1_1_mat_t, 5)
                    pi_vi_vj = pi_vi_vj + sum(pi_vi_vj_t, 5)

                end if
            end do

        end if

        if ( verbose ) then
            print*, '----------------------------------------'
            print*
        end if

    end subroutine

    !! fill in both of these
    subroutine compute_self_energies_no_spin(nml_filename, tran_form_1_t, &
            tran_form_v_t, tran_form_v2_t, omega_iipk_t, n_init, n_fin, n_k, n_proc, &
            proc_id, root_process, verbose)

        implicit none

        character(len=*) :: nml_filename

        complex(dp) :: tran_form_1_t(:, :)
        complex(dp) :: tran_form_v_t(:, :, :)
        complex(dp) :: tran_form_v2_t(:, :)

        real(dp) :: omega_iipk_t(:, :)

        integer :: n_init, n_fin
        integer :: n_k
        integer :: n_proc, proc_id, root_process

        logical, optional :: verbose

        integer :: t
        integer :: tran_id, init_id, fin_id

        ! compute the pi's
        do t = 1, n_tran_per_proc

            tran_id = job_table(proc_id + 1, t)

            if ( tran_id /= 0 ) then

                init_id = tran_to_init_fin_id(tran_id, 1)
                fin_id = tran_to_init_fin_id(tran_id, 2)

                call calc_pi_v2_v2(pi_v2_v2_t(:, :, t), tran_form_v2_t(t, :), omega_iipk_t(t, :))
                call calc_pi_1_1_mat(pi_1_1_mat_t(:, :, :, :, t), tran_form_v_t(:, t, :), omega_iipk_t(t, :))
                call calc_pi_vi_vj(pi_vi_vj_t(:, :, :, :, t), tran_form_v_t(:, t, :), omega_iipk_t(t, :))

            end if

        end do

    end subroutine

    subroutine compute_self_energies_spin(nml_filename, tran_form_1_t, &
            tran_form_v_t, tran_form_v2_t, omega_iipk_t, n_init, n_fin, n_k, n_proc, &
            proc_id, root_process, verbose)

        implicit none

        character(len=*) :: nml_filename

        complex(dp) :: tran_form_1_t(:, :, :, :)
        complex(dp) :: tran_form_v_t(:, :, :, :, :)
        complex(dp) :: tran_form_v2_t(:, :, :, :)

        real(dp) :: omega_iipk_t(:, :)

        integer :: n_init, n_fin
        integer :: n_k
        integer :: n_proc, proc_id, root_process

        logical, optional :: verbose

        integer :: t
        integer :: tran_id, init_id, fin_id

        ! compute the pi's
        do t = 1, n_tran_per_proc

            tran_id = job_table(proc_id + 1, t)

            if ( tran_id .ne. 0 ) then

                init_id = tran_to_init_fin_id(tran_id, 1)
                fin_id = tran_to_init_fin_id(tran_id, 2)

                call calc_pi_v2_v2(pi_v2_v2_t(:, :, t), tran_form_v2_t(t, :, :, :), omega_iipk_t(t, :))
                call calc_pi_1_1_mat(pi_1_1_mat_t(:, :, :, :, t), tran_form_v_t(:, t, :, :, :), omega_iipk_t(t, :))
                call calc_pi_vi_vj(pi_vi_vj_t(:, :, :, :, t), tran_form_v_t(:, t, :, :, :), omega_iipk_t(t, :))

            end if

        end do

    end subroutine

    subroutine calc_pi_v2_v2_no_spin(pi_v2_v2, tran_form_v2, omega_iipk)

        implicit none

        complex(dp) :: pi_v2_v2(:, :)

        complex(dp) :: tran_form_v2(:)

        real(dp) :: omega_iipk(:)

        integer :: k, w, p

        real(dp) :: omega, width, omega_k

        do p = 1, n_widths
            do w = 1, n_omega

                omega = omega_list(w)
                width = delta_list(w, p)

                do k = 1, n_k

                    omega_k = omega_iipk(k)

                    if ( abs(omega - omega_k) < sigma_gamma*width ) then

                        pi_v2_v2(w, p) = pi_v2_v2(w, p) + &
                           (pc_vol)**(-1)*k_weight(k)*(0.25_dp)*&
                           (&
                                ( omega - omega_k + ii*width )**(-1) - &
                                ( omega + omega_k - ii*width )**(-1) &
                            )*&
                           abs(tran_form_v2(k))**2*(spin_degen/2.0_dp)

                    end if

                end do

            end do
        end do

    end subroutine

    subroutine calc_pi_v2_v2_spin(pi_v2_v2, tran_form_v2, omega_iipk)

        implicit none

        complex(dp) :: pi_v2_v2(:, :)

        complex(dp) :: tran_form_v2(:, :, :)

        complex(dp) :: tf

        real(dp) :: omega_iipk(:)

        integer :: k, w, p
        integer :: s, sp

        real(dp) :: omega, width, omega_k

        do p = 1, n_widths
            do w = 1, n_omega

                omega = omega_list(w)
                width = delta_list(w, p)

                do k = 1, n_k

                    omega_k = omega_iipk(k)

                    tf = (0.0_dp, 0.0_dp)

                    do s = 1, 2
                        tf = tf + tran_form_v2(k, s, s) 
                    end do

                    if ( abs(omega - omega_k) < sigma_gamma*width ) then

                        pi_v2_v2(w, p) = pi_v2_v2(w, p) + &
                           (pc_vol)**(-1)*k_weight(k)*(0.25_dp)*&
                           (&
                                ( omega - omega_k + ii*width )**(-1) - &
                                ( omega + omega_k - ii*width )**(-1) &
                            )*&
                           abs(tf)**2*(spin_degen/2.0_dp)

                    end if

                end do

            end do

        end do

    end subroutine

    subroutine calc_pi_1_1_mat_no_spin(pi_1_1_mat, tran_form_v, omega_iipk)

        implicit none

        complex(dp) :: pi_1_1_mat(:, :, :, :)

        complex(dp) :: tran_form_v(:, :)

        real(dp) :: omega_iipk(:)

        complex(dp) :: outer_ff(3, 3, n_k)

        real(dp) :: omega_k

        integer :: i, j
        integer :: k

        integer :: p, w

        real(dp) :: width, omega

        do k = 1, n_k
            do i = 1, 3
                do j = 1, 3

                    omega_k = omega_iipk(k)

                    outer_ff(i, j, k) = (m_elec/omega_k)**2*conjg(tran_form_v(i, k))*tran_form_v(j, k)

                end do
            end do
        end do

        do p = 1, n_widths
            do w = 1, n_omega

                omega = omega_list(w)
                width = delta_list(w, p)

                do k = 1, n_k

                    omega_k = omega_iipk(k)

                    if ( abs(omega - omega_k) < sigma_gamma*width ) then

                        pi_1_1_mat(:, :, w, p) = pi_1_1_mat(:, :, w, p) + &
                           (pc_vol)**(-1)*k_weight(k)*&
                           ( &
                                ( omega - omega_k + ii*width )**(-1) - &
                                ( omega + omega_k - ii*width )**(-1) &
                            )*&
                           outer_ff(:, :, k)*(spin_degen/2.0_dp)

                    end if

                end do

            end do
        end do

    end subroutine

    subroutine calc_pi_1_1_mat_spin(pi_1_1_mat, tran_form_v, omega_iipk)

        implicit none

        complex(dp) :: pi_1_1_mat(:, :, :, :)

        complex(dp) :: tran_form_v(:, :, :, :)

        real(dp) :: omega_iipk(:)

        complex(dp) :: outer_ff(3, 3, n_k)

        real(dp) :: omega_k

        integer :: i, j
        integer :: k

        integer :: s

        integer :: p, w

        real(dp) :: width, omega

        complex(dp) :: tf_v(3)

        do k = 1, n_k

            tf_v = (0.0_dp, 0.0_dp)

            do s = 1, 2

                tf_v = tf_v + tran_form_v(:, k, s, s)

            end do

            do i = 1, 3
                do j = 1, 3

                    omega_k = omega_iipk(k)

                    outer_ff(i, j, k) = (m_elec/omega_k)**2*conjg(tf_v(i))*tf_v(j)

                end do
            end do
        end do

        do p = 1, n_widths
            do w = 1, n_omega

                omega = omega_list(w)
                width = delta_list(w, p)

                do k = 1, n_k

                    omega_k = omega_iipk(k)

                    if ( abs(omega - omega_k) < sigma_gamma*width ) then

                        pi_1_1_mat(:, :, w, p) = pi_1_1_mat(:, :, w, p) + &
                           (pc_vol)**(-1)*k_weight(k)*&
                           ( &
                                ( omega - omega_k + ii*width )**(-1) - &
                                ( omega + omega_k - ii*width )**(-1) &
                            )*&
                           outer_ff(:, :, k)*(spin_degen/2.0_dp)

                    end if

                end do

            end do
        end do

    end subroutine

    subroutine calc_pi_vi_vj_no_spin(pi_vi_vj, tran_form_v, omega_iipk)

        implicit none

        complex(dp) :: pi_vi_vj(:, :, :, :)

        complex(dp) :: tran_form_v(:, :)

        real(dp) :: omega_iipk(:)

        complex(dp) :: outer_ff(3, 3, n_k)

        real(dp) :: omega_k

        integer :: i, j
        integer :: k

        integer :: p, w

        real(dp) :: width, omega

        do k = 1, n_k
            do i = 1, 3
                do j = 1, 3

                    outer_ff(i, j, k) = conjg(tran_form_v(i, k))*tran_form_v(j, k)

                end do
            end do
        end do

        do p = 1, n_widths
            do w = 1, n_omega

                omega = omega_list(w)
                width = delta_list(w, p)

                do k = 1, n_k

                    omega_k = omega_iipk(k)

                    if ( abs(omega - omega_k) < sigma_gamma*width ) then

                        pi_vi_vj(:, :, w, p) = pi_vi_vj(:, :, w, p) + &
                           (pc_vol)**(-1)*k_weight(k)*&
                           ( &
                                ( omega - omega_k + ii*width )**(-1) - &
                                ( omega + omega_k - ii*width )**(-1) &
                            )*&
                           outer_ff(:, :, k)*(spin_degen/2.0_dp)

                    end if

                end do

            end do
        end do

    end subroutine

    subroutine calc_pi_vi_vj_spin(pi_vi_vj, tran_form_v, omega_iipk)

        implicit none

        complex(dp) :: pi_vi_vj(:, :, :, :)

        complex(dp) :: tran_form_v(:, :, :, :)

        real(dp) :: omega_iipk(:)

        complex(dp) :: outer_ff(3, 3, n_k)

        real(dp) :: omega_k

        integer :: i, j
        integer :: k

        integer :: s

        integer :: p, w

        real(dp) :: width, omega

        complex(dp) :: tf_v(3)

        do k = 1, n_k

            tf_v = (0.0_dp, 0.0_dp)

            do s = 1, 2

                tf_v = tf_v + tran_form_v(:, k, s, s)

            end do

            do i = 1, 3
                do j = 1, 3

                    outer_ff(i, j, k) = conjg(tf_v(i))*tf_v(j)

                end do
            end do
        end do

        do p = 1, n_widths
            do w = 1, n_omega

                omega = omega_list(w)
                width = delta_list(w, p)

                do k = 1, n_k

                    omega_k = omega_iipk(k)

                    if ( abs(omega - omega_k) < sigma_gamma*width ) then

                        pi_vi_vj(:, :, w, p) = pi_vi_vj(:, :, w, p) + &
                           (pc_vol)**(-1)*k_weight(k)*&
                           ( &
                                ( omega - omega_k + ii*width )**(-1) - &
                                ( omega + omega_k - ii*width )**(-1) &
                            )*&
                           outer_ff(:, :, k)*(spin_degen/2.0_dp)

                    end if

                end do

            end do
        end do

    end subroutine

    subroutine save_self_energies(filename, verbose)

        implicit none

        character(len=*) :: filename

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id

        logical, optional :: verbose
        logical :: file_exists

        integer(HSIZE_T) :: dims1(1) = [1]
        integer(HSIZE_T) :: dims2(2)
        integer(HSIZE_T) :: dims3(3)

        integer :: error

        integer :: m, f, t
        integer :: i, fin
        character(len=64) :: i_str, fin_str

        if ( verbose ) then

            print*, 'Saving self energies...'
            print*

        end if

        ! make sure the file exists
        inquire(file = trim(filename), exist = file_exists)

        if ( file_exists ) then

            call h5open_f(error)
            call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, error)

            call h5gcreate_f(file_id, 'self_energies', group_id, error)

            do i = 1, n_widths

                write(i_str, *) i
                i_str = trim(adjustl(i_str))

                call h5gcreate_f(file_id,& 
                    'self_energies'//&
                    '/width_'//trim(i_str),&
                    group_id, error)

                dims1 = [n_omega]

                call h5ltmake_dataset_double_f(file_id,&
                    'self_energies'//&
                    '/width_'//trim(i_str)//&
                    '/pi_v2_v2_r', &
                    size(dims1), dims1,&
                    real(pi_v2_v2(:, i)), error)

                call h5ltmake_dataset_double_f(file_id,&
                    'self_energies'//&
                    '/width_'//trim(i_str)//&
                    '/pi_v2_v2_c', &
                    size(dims1), dims1,&
                    aimag(pi_v2_v2(:, i)), error)

                dims3 = [3, 3, n_omega]

                call h5ltmake_dataset_double_f(file_id,&
                    'self_energies'//&
                    '/width_'//trim(i_str)//&
                    '/pi_1_1_mat_r', &
                    size(dims3), dims3,&
                    real(pi_1_1_mat(:, :, :, i)), error)

                call h5ltmake_dataset_double_f(file_id,&
                    'self_energies'//&
                    '/width_'//trim(i_str)//&
                    '/pi_1_1_mat_c', &
                    size(dims3), dims3,&
                    aimag(pi_1_1_mat(:, :, :, i)), error)

                call h5ltmake_dataset_double_f(file_id,&
                    'self_energies'//&
                    '/width_'//trim(i_str)//&
                    '/pi_vi_vj_r', &
                    size(dims3), dims3,&
                    real(pi_vi_vj(:, :, :, i)), error)

                call h5ltmake_dataset_double_f(file_id,&
                    'self_energies'//&
                    '/width_'//trim(i_str)//&
                    '/pi_vi_vj_c', &
                    size(dims3), dims3,&
                    aimag(pi_vi_vj(:, :, :, i)), error)

            end do

            call h5fclose_f(file_id, error)
            call h5close_f(error)

            if ( verbose ) then
                print*, '----------------------------------------'
                print*
            end if

        else

            if ( verbose ) then

                print*, '!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*
                print*, '   Output file : ', trim(filename), ' does NOT exist.'
                print*
                print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print*

            end if

            stop

        end if

    end subroutine

end module
