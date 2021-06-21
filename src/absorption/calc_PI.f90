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

    complex(dp), allocatable :: pi_vv(:, :)
        !! Dim : [n_omega, n_widths]
        !! 
        !! self energy with two v^2 insertions
        !!
        !! Units : eV^2

    complex(dp), allocatable :: pi_11_mat(:, :, :, :)
        !! Dim : [3, 3, n_omega, n_widths] 
        !!
        !! self energy, without q_vec's
        !!
        !! Pi_11 = (q/m_elec) . Pi_11_mat . (q/m_elec)
        !!
        !! Units : eV^2

    complex(dp), allocatable :: pi_v1(:, :, :)
        !! Dim : [3, n_omega, n_widths]
        !! 
        !! self energy Pi_{v^2, 1}
        !!
        !! Units : eV^2

    complex(dp), allocatable :: pi_1v(:, :, :)
        !! Dim : [3, n_omega, n_widths]
        !! 
        !! self energy Pi_{1, v^2} = Pi_{v^2, 1}^*(-delta)
        !!
        !! Units : eV^2

    !real(dp), allocatable :: im_pi_vv(:, :)
    !    !! Dim : [n_omega, n_widths]
    !    !! 
    !    !! imaginary part of self energy with two v^2 insertions
    !    !!
    !    !! Units : None

    !real(dp), allocatable :: im_pi_11_mat(:, :, :, :)
    !    !! Dim : [3, 3, n_omega, n_widths]
    !    !! 
    !    !! imaginary part of self energy (without q vectors)
    !    !!
    !    !! Units : eV^(-2)

    complex(dp), allocatable :: pi_vv_t(:, :, :)
        !! Dim : [n_omega, n_widths, n_tran_per_proc]
        !! 
        !! self energy with two v^2 insertions
        !!
        !! Units : eV^2

    complex(dp), allocatable :: pi_11_mat_t(:, :, :, :, :)
        !! Dim : [3, 3, n_omega, n_widths, n_tran_per_proc] 
        !!
        !! self energy, without q_vec's
        !!
        !! Units : eV^2

    complex(dp), allocatable :: pi_v1_t(:, :, :, :)
        !! Dim : [3, n_omega, n_widths, n_tran_per_proc]
        !! 
        !! self energy Pi_{v^2, 1}
        !!
        !! Units : eV^2

    complex(dp), allocatable :: pi_1v_t(:, :, :, :)
        !! Dim : [3, n_omega, n_widths, n_tran_per_proc]
        !! 
        !! self energy Pi_{1, v^2} = Pi_{v^2, 1}^*(-delta)
        !!
        !! Units : eV^2

    !real(dp), allocatable :: im_pi_vv_t(:, :, :)
    !    !! Dim : [n_omega, n_widths, n_tran_per_proc]
    !    !! 
    !    !! imaginary part of self energy with two v^2 insertions
    !    !!
    !    !! Units : None

    !real(dp), allocatable :: im_pi_11_mat_t(:, :, :, :, :)
    !    !! Dim : [3, 3, n_omega, n_widths, n_tran_per_proc]
    !    !! 
    !    !! imaginary part of self energy (without q vectors)
    !    !!
    !    !! Units : eV^(-2)

contains

    subroutine compute_self_energies(nml_filename, tran_form, &
            tran_form_vec, omega_iipk, n_init, n_fin, n_k, n_proc, &
            proc_id, root_process, verbose)

        implicit none

        character(len=*) :: nml_filename

        integer :: n_init, n_fin, n_k

        integer :: n_proc, proc_id, root_process

        integer :: t, init_id, fin_id, tran_id
        integer :: i

        complex(dp) :: tran_form(n_init, n_fin, n_k)
        complex(dp) :: tran_form_vec(3, n_init, n_fin, n_k)
        real(dp) :: omega_iipk(n_init, n_fin, n_k)

        complex(dp) :: tran_form_k(n_k)
        complex(dp) :: tran_form_vec_k(3, n_k)
        real(dp) :: omega_iipk_k(n_k)

        integer :: status(MPI_STATUS_SIZE)
        integer :: tag = 0
        integer :: err

        logical, optional :: verbose

        if ( verbose ) then

            print*, 'Calculating self energies...'
            print*

        end if

        ! initialize variables
        allocate(pi_vv_t(n_omega, n_widths, n_tran_per_proc))
        pi_vv_t = (0.0_dp, 0.0_dp)

        allocate(pi_11_mat_t(3, 3, n_omega, n_widths, n_tran_per_proc))
        pi_11_mat_t = (0.0_dp, 0.0_dp)

        allocate(pi_1v_t(3, n_omega, n_widths, n_tran_per_proc))
        pi_1v_t = (0.0_dp, 0.0_dp)

        allocate(pi_v1_t(3, n_omega, n_widths, n_tran_per_proc))
        pi_v1_t = (0.0_dp, 0.0_dp)

        ! allocate(im_pi_vv_t(n_omega, n_widths, n_tran_per_proc))
        ! im_pi_vv_t = 0.0_dp

        ! allocate(im_pi_11_mat_t(3, 3, n_omega, n_widths, n_tran_per_proc))
        ! im_pi_11_mat_t = 0.0_dp

        if ( proc_id == root_process ) then

            allocate(pi_vv(n_omega, n_widths))
            pi_vv = (0.0_dp, 0.0_dp)

            allocate(pi_11_mat(3, 3, n_omega, n_widths))
            pi_11_mat = (0.0_dp, 0.0_dp)

            allocate(pi_1v(3, n_omega, n_widths))
            pi_1v = (0.0_dp, 0.0_dp)

            allocate(pi_v1(3, n_omega, n_widths))
            pi_v1 = (0.0_dp, 0.0_dp)

            ! allocate(im_pi_vv(n_omega, n_widths))
            ! im_pi_vv = 0.0_dp

            ! allocate(im_pi_11_mat(3, 3, n_omega, n_widths))
            ! im_pi_11_mat = 0.0_dp

        end if

        ! each processor computes
        do t = 1, n_tran_per_proc

            tran_id = job_table(proc_id + 1, t)

            if ( tran_id .ne. 0 ) then

                init_id = tran_to_init_fin_id(tran_id, 1)
                fin_id = tran_to_init_fin_id(tran_id, 2)

                tran_form_k = tran_form(init_id, fin_id, :)
                tran_form_vec_k = tran_form_vec(:, init_id, fin_id, :)
                omega_iipk_k = omega_iipk(init_id, fin_id, :)

                call calc_pi_vv(pi_vv_t(:, :, t), tran_form_k, omega_iipk_k)
                call calc_pi_11_mat(pi_11_mat_t(:, :, :, :, t), tran_form_vec_k, &
                    omega_iipk_k)
                call calc_pi_1v(pi_1v_t(:, :, :, t), tran_form_k, tran_form_vec_k, &
                    omega_iipk_k)
                call calc_pi_v1(pi_v1_t(:, :, :, t), tran_form_k, tran_form_vec_k, &
                    omega_iipk_k)

                ! call calc_im_pi_vv(im_pi_vv_t(:, :, t), tran_form_k, omega_iipk_k)
                ! call calc_im_pi_11_mat(im_pi_11_mat_t(:, :, :, :, t), tran_form_vec_k, &
                !     omega_iipk_k)

            end if

        end do

        ! sum together again 
        if ( proc_id .ne. root_process ) then

            call MPI_SEND(pi_vv_t, &
               size(pi_vv_t), MPI_DOUBLE_COMPLEX, root_process, tag, MPI_COMM_WORLD, err)

            call MPI_SEND(pi_11_mat_t, &
               size(pi_11_mat_t), MPI_DOUBLE_COMPLEX, root_process, tag, MPI_COMM_WORLD, err)

            call MPI_SEND(pi_1v_t, &
               size(pi_1v_t), MPI_DOUBLE_COMPLEX, root_process, tag, MPI_COMM_WORLD, err)

            call MPI_SEND(pi_v1_t, &
               size(pi_v1_t), MPI_DOUBLE_COMPLEX, root_process, tag, MPI_COMM_WORLD, err)

            ! call MPI_SEND(im_pi_vv_t, &
            !    size(im_pi_vv_t), MPI_DOUBLE, root_process, tag, MPI_COMM_WORLD, err)

            ! call MPI_SEND(im_pi_11_mat_t, &
            !    size(im_pi_11_mat_t), MPI_DOUBLE, root_process, tag, MPI_COMM_WORLD, err)

        end if

        if ( proc_id .eq. root_process ) then

            ! add main processors contribution
            pi_vv = pi_vv + sum(pi_vv_t, 3)
            pi_11_mat = pi_11_mat + sum(pi_11_mat_t, 5)
            pi_1v = pi_1v + sum(pi_1v_t, 4)
            pi_v1 = pi_v1 + sum(pi_v1_t, 4)

            ! im_pi_vv = im_pi_vv + sum(im_pi_vv_t, 3)
            ! im_pi_11_mat = im_pi_11_mat + sum(im_pi_11_mat_t, 5)

            do i = 1, n_proc
                if ( (i - 1) .ne. root_process ) then

                    call MPI_RECV(pi_vv_t, &
                       size(pi_vv_t), MPI_DOUBLE_COMPLEX, i - 1, MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                    call MPI_RECV(pi_11_mat_t, &
                       size(pi_11_mat_t), MPI_DOUBLE_COMPLEX, i - 1, MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                    call MPI_RECV(pi_1v_t, &
                       size(pi_1v_t), MPI_DOUBLE_COMPLEX, i - 1, MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                    call MPI_RECV(pi_v1_t, &
                       size(pi_v1_t), MPI_DOUBLE_COMPLEX, i - 1, MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                    ! call MPI_RECV(im_pi_vv_t, &
                    !    size(im_pi_vv_t), MPI_DOUBLE, i - 1, MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                    ! call MPI_RECV(im_pi_11_mat_t, &
                    !    size(im_pi_11_mat_t), MPI_DOUBLE, i - 1, MPI_ANY_TAG, MPI_COMM_WORLD, status, err)

                    ! add other processors contributions
                    pi_vv = pi_vv + sum(pi_vv_t, 3)
                    pi_11_mat = pi_11_mat + sum(pi_11_mat_t, 5)
                    pi_1v = pi_1v + sum(pi_1v_t, 4)
                    pi_v1 = pi_v1 + sum(pi_v1_t, 4)

                    ! im_pi_vv = im_pi_vv + sum(im_pi_vv_t, 3)
                    ! im_pi_11_mat = im_pi_11_mat + sum(im_pi_11_mat_t, 5)

                end if
            end do

        end if

        if ( verbose ) then
            print*, '----------'
            print*
        end if

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
                    '/pi_vv_r', &
                    size(dims1), dims1,&
                    real(pi_vv(:, i)), error)

                call h5ltmake_dataset_double_f(file_id,&
                    'self_energies'//&
                    '/width_'//trim(i_str)//&
                    '/pi_vv_c', &
                    size(dims1), dims1,&
                    aimag(pi_vv(:, i)), error)

                dims3 = [3, 3, n_omega]

                call h5ltmake_dataset_double_f(file_id,&
                    'self_energies'//&
                    '/width_'//trim(i_str)//&
                    '/pi_11_mat_r', &
                    size(dims3), dims3,&
                    real(pi_11_mat(:, :, :, i)), error)

                call h5ltmake_dataset_double_f(file_id,&
                    'self_energies'//&
                    '/width_'//trim(i_str)//&
                    '/pi_11_mat_c', &
                    size(dims3), dims3,&
                    aimag(pi_11_mat(:, :, :, i)), error)

                dims2 = [3, n_omega]

                call h5ltmake_dataset_double_f(file_id,&
                    'self_energies'//&
                    '/width_'//trim(i_str)//&
                    '/pi_1v_r', &
                    size(dims2), dims2,&
                    real(pi_1v(:, :, i)), error)

                call h5ltmake_dataset_double_f(file_id,&
                    'self_energies'//&
                    '/width_'//trim(i_str)//&
                    '/pi_1v_c', &
                    size(dims2), dims2,&
                    aimag(pi_1v(:, :, i)), error)

                call h5ltmake_dataset_double_f(file_id,&
                    'self_energies'//&
                    '/width_'//trim(i_str)//&
                    '/pi_v1_r', &
                    size(dims2), dims2,&
                    real(pi_v1(:, :, i)), error)

                call h5ltmake_dataset_double_f(file_id,&
                    'self_energies'//&
                    '/width_'//trim(i_str)//&
                    '/pi_v1_c', &
                    size(dims2), dims2,&
                    aimag(pi_v1(:, :, i)), error)

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

    ! subroutine calc_im_pi_11_mat(im_pi_11_mat, tran_form_vec_t, omega_iipk_t)

    !     implicit none

    !     real(dp) :: im_pi_11_mat(3, 3, n_omega, n_widths)
    !     complex(dp) :: outer_ff(3, 3, n_k)

    !     complex(dp) :: tran_form_vec_t(3, n_k)
    !     real(dp) :: omega_iipk_t(n_k)

    !     integer :: k, w, p, i, j

    !     real(dp) :: omega, width, omega_k

    !     do k = 1, n_k
    !         do i = 1, 3
    !             do j = 1, 3

    !                 outer_ff(i, j, k) = conjg(tran_form_vec_t(i, k))*tran_form_vec_t(j, k)

    !             end do
    !         end do
    !     end do

    !     do p = 1, n_widths
    !         do w = 1, n_omega

    !             omega = omega_list(w)
    !             width = delta_list(w, p)

    !             do k = 1, n_k

    !                 omega_k = omega_iipk_t(k)

    !                 ! im_pi_11_mat(:, :, w, p) = im_pi_11_mat(:, :, w, p) + &
    !                 !    (pc_vol)**(-1)*k_weight(k)*&
    !                 !    (4.0_dp*omega_k**2*width)*&
    !                 !    ( ( omega**2 - omega_k**2 )**2 + ( 2.0_dp*omega_k*width )**2 )**(-1)*&
    !                 !    outer_ff(:, :, k)

    !                 im_pi_11_mat(:, :, w, p) = im_pi_11_mat(:, :, w, p) + &
    !                    (pc_vol)**(-1)*k_weight(k)*&
    !                    (width)*&
    !                    (&
    !                     ( ( omega - omega_k )**2 + ( width )**2 )**(-1) + &
    !                     ( ( omega + omega_k )**2 + ( width )**2 )**(-1) &
    !                    )*&
    !                    outer_ff(:, :, k)

    !             end do

    !         end do
    !     end do

    ! end subroutine

    ! subroutine calc_im_pi_vv(im_pi_vv, tran_form_t, omega_iipk_t)

    !     implicit none

    !     real(dp) :: im_pi_vv(n_omega, n_widths)

    !     complex(dp) :: tran_form_t(n_k)
    !     real(dp) :: omega_iipk_t(n_k)

    !     integer :: k, w, p

    !     real(dp) :: omega, width, omega_k

    !     do p = 1, n_widths
    !         do w = 1, n_omega

    !             omega = omega_list(w)
    !             width = delta_list(w, p)

    !             print*, 'omega = ', omega
    !             print*, 'width = ', width

    !             do k = 1, n_k

    !                 omega_k = omega_iipk_t(k)

    !                 ! im_pi_vv(w, p) = im_pi_vv(w, p) - &
    !                 !    (pc_vol)**(-1)*k_weight(k)*&
    !                 !    (4.0_dp*omega_k**2*width)*&
    !                 !    ( ( omega**2 - omega_k**2 )**2 + ( 2.0_dp*omega_k*width )**2 )**(-1)*&
    !                 !    abs(tran_form_t(k))**2 

    !                 im_pi_vv(w, p) = im_pi_vv(w, p) - &
    !                    (pc_vol)**(-1)*k_weight(k)*&
    !                    (width)*&
    !                    (&
    !                     ( ( omega - omega_k )**2 + ( width )**2 )**(-1) + &
    !                     ( ( omega + omega_k )**2 + ( width )**2 )**(-1) &
    !                    )*&
    !                    abs(tran_form_t(k))**2 

    !             end do

    !         end do
    !     end do

    ! end subroutine

    subroutine calc_pi_vv(pi_vv, tran_form_t, omega_iipk_t)

        implicit none

        complex(dp) :: pi_vv(n_omega, n_widths)

        complex(dp) :: tran_form_t(n_k)
        real(dp) :: omega_iipk_t(n_k)

        integer :: k, w, p

        real(dp) :: omega, width, omega_k

        do p = 1, n_widths
            do w = 1, n_omega

                omega = omega_list(w)
                width = delta_list(w, p)

                do k = 1, n_k

                    omega_k = omega_iipk_t(k)

                    if ( abs(omega - omega_k) < sigma_gamma*width ) then

                        pi_vv(w, p) = pi_vv(w, p) + &
                           (pc_vol)**(-1)*k_weight(k)*&
                           (&
                                ( omega - omega_k + ii*width )**(-1) - &
                                ( omega + omega_k - ii*width )**(-1) &
                            )*&
                           abs(tran_form_t(k))**2*(1.0_dp/2.0_dp)**2 

                    end if

                end do

            end do
        end do

    end subroutine

    subroutine calc_pi_11_mat(pi_11_mat, tran_form_vec_t, omega_iipk_t)

        implicit none

        complex(dp) :: pi_11_mat(3, 3, n_omega, n_widths)
        complex(dp) :: outer_ff(3, 3, n_k)

        complex(dp) :: tran_form_vec_t(3, n_k)
        real(dp) :: omega_iipk_t(n_k)

        integer :: k, w, p, i, j

        real(dp) :: omega, width, omega_k

        do k = 1, n_k
            do i = 1, 3
                do j = 1, 3

                    omega_k = omega_iipk_t(k)

                    outer_ff(i, j, k) = (m_elec/omega_k)**2*conjg(tran_form_vec_t(i, k))*tran_form_vec_t(j, k)

                end do
            end do
        end do

        do p = 1, n_widths
            do w = 1, n_omega

                omega = omega_list(w)
                width = delta_list(w, p)

                do k = 1, n_k

                    omega_k = omega_iipk_t(k)

                    if ( abs(omega - omega_k) < sigma_gamma*width ) then

                        pi_11_mat(:, :, w, p) = pi_11_mat(:, :, w, p) + &
                           (pc_vol)**(-1)*k_weight(k)*&
                           ( &
                                ( omega - omega_k + ii*width )**(-1) - &
                                ( omega + omega_k - ii*width )**(-1) &
                            )*&
                           outer_ff(:, :, k)

                    end if

                end do

            end do
        end do

    end subroutine

    subroutine calc_pi_v1(pi_v1, tran_form_t, tran_form_vec_t, omega_iipk_t)

        implicit none

        complex(dp) :: pi_v1(3, n_omega, n_widths)
        complex(dp) :: tran_form_t(n_k)
        complex(dp) :: tran_form_vec_t(3, n_k)
        real(dp) :: omega_iipk_t(n_k)

        integer :: k, w, p

        real(dp) :: omega, width, omega_k

        do p = 1, n_widths
            do w = 1, n_omega

                omega = omega_list(w)
                width = delta_list(w, p)

                do k = 1, n_k

                    omega_k = omega_iipk_t(k)

                    if ( abs(omega - omega_k) < sigma_gamma*width ) then

                        pi_v1(:, w, p) = pi_v1(:, w, p) + &
                           (pc_vol)**(-1)*k_weight(k)*&
                           (&
                                ( omega - omega_k + ii*width )**(-1)*&
                                conjg(tran_form_vec_t(:, k))*tran_form_t(k) + &
                                ( omega + omega_k - ii*width )**(-1)*&
                                conjg(tran_form_t(k))*tran_form_vec_t(:, k) &
                            )*(1.0_dp/2.0_dp)*(m_elec/omega_k)

                    end if

                end do

            end do
        end do

    end subroutine

    subroutine calc_pi_1v(pi_1v, tran_form_t, tran_form_vec_t, omega_iipk_t)

        implicit none

        complex(dp) :: pi_1v(3, n_omega, n_widths)
        complex(dp) :: tran_form_t(n_k)
        complex(dp) :: tran_form_vec_t(3, n_k)
        real(dp) :: omega_iipk_t(n_k)

        integer :: k, w, p

        real(dp) :: omega, width, omega_k

        do p = 1, n_widths
            do w = 1, n_omega

                omega = omega_list(w)
                width = delta_list(w, p)

                do k = 1, n_k

                    omega_k = omega_iipk_t(k)

                    if ( abs(omega - omega_k) < sigma_gamma*width ) then

                        pi_1v(:, w, p) = pi_1v(:, w, p) + &
                           (pc_vol)**(-1)*k_weight(k)*&
                           (&
                                ( omega - omega_k + ii*width )**(-1)*&
                                conjg(tran_form_t(k))*tran_form_vec_t(:, k) + &
                                ( omega + omega_k - ii*width )**(-1)*&
                                 conjg(tran_form_vec_t(:, k))*tran_form_t(k)&
                            )*(1.0_dp/2.0_dp)*(m_elec/omega_k)
                    end if

                end do

            end do
        end do

    end subroutine

end module
