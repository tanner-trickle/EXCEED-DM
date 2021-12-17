module abs_calc_PI
    !! Computes self-energies, \( \Pi_{\mathcal{O}_1, \mathcal{O}_2} \), for absorption rate calculations.

    use mpi
    use hdf5
    use h5lt

    use prec
    use constants

    use width_parameters_type
    use PW_dataset_type
    use dm_model_type
    use material_type

    use MPI_util

    implicit none

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

    interface calc_abs_Pi
        module procedure calc_abs_Pi_no_spin
        module procedure calc_abs_Pi_spin
    end interface

contains

    subroutine calc_abs_Pi_no_spin(proc_id, root_process, &
            tran_form_1_job, tran_form_v_job, tran_form_v2_job, &
            pi_v2_v2, pi_vi_vj, pi_1_1_mat, ik_manager, job_id_to_ik, &
            PW_dataset, widths, dm_model, target_mat, filename, verbose)
        !! Compute the self energies with spin independent wave functions.
        implicit none

        integer :: proc_id
        integer :: root_process
        complex(dp) :: tran_form_1_job(:, :)
        complex(dp) :: tran_form_v_job(:, :, :)
        complex(dp) :: tran_form_v2_job(:, :)
        complex(dp) :: pi_v2_v2(:, :)
        complex(dp) :: pi_vi_vj(:, :, :, :)
        complex(dp) :: pi_1_1_mat(:, :, :, :)
        type(parallel_manager_t) :: ik_manager
        integer :: job_id_to_ik(:, :)
        type(PW_dataset_t) :: PW_dataset
        type(width_parameters_t) :: widths
        type(dm_model_t) :: dm_model
        type(material_t) :: target_mat
        character(len=*) :: filename
        logical, optional :: verbose

        complex(dp), allocatable :: pi_v2_v2_job(:, :, :)
            !! Dim : [n_omega, n_widths, n_tran_per_proc]
            !! 
            !! self energy with two v^2 insertions
            !!
            !! Units : eV^2

        complex(dp), allocatable :: pi_1_1_mat_job(:, :, :, :, :)
            !! Dim : [3, 3, n_omega, n_widths, n_tran_per_proc] 
            !!
            !! self energy, without q_vec's
            !!
            !! Units : eV^2

        complex(dp), allocatable :: pi_vi_vj_job(:, :, :, :, :)
            !! Dim : [3, 3, n_omega, n_widths, n_tran_per_proc] 
            !!
            !! Units : eV^2

        integer :: j, job_id, val_id, k

        if ( verbose ) then

            print*, 'Computing self energies...'
            print*

        end if

        ! allocate job variables
        allocate(pi_v2_v2_job(size(pi_v2_v2, 1), size(pi_v2_v2, 2), ik_manager%n_jobs_per_proc))
        pi_v2_v2_job = (0.0_dp, 0.0_dp)
        allocate(pi_vi_vj_job(3, 3, size(pi_v2_v2, 1), size(pi_v2_v2, 2), ik_manager%n_jobs_per_proc))
        pi_vi_vj_job = (0.0_dp, 0.0_dp)
        allocate(pi_1_1_mat_job(3, 3, size(pi_v2_v2, 1), size(pi_v2_v2, 2), ik_manager%n_jobs_per_proc))
        pi_1_1_mat_job = (0.0_dp, 0.0_dp)

        do j = 1, ik_manager%n_jobs_per_proc

            job_id = ik_manager%job_table(proc_id + 1, j)

            if ( job_id /= 0 ) then

                val_id = job_id_to_ik(job_id, 1)
                k = job_id_to_ik(job_id, 2)

                call calc_pi_v2_v2(pi_v2_v2_job(:, :, j), tran_form_v2_job(j, :), &
                    dm_model, widths, PW_dataset, target_mat, val_id, k)
                call calc_pi_1_1_mat(pi_1_1_mat_job(:, :, :, :, j), tran_form_v_job(:, j, :), &
                    dm_model, widths, PW_dataset, target_mat, val_id, k)
                call calc_pi_vi_vj(pi_vi_vj_job(:, :, :, :, j), tran_form_v_job(:, j, :), &
                    dm_model, widths, PW_dataset, target_mat, val_id, k)

            end if

        end do

        if ( verbose ) then

            print*, 'Done computing self energies!'
            print*

        end if

        call ik_manager%comm_self_energies(proc_id, root_process, &
            pi_v2_v2_job, pi_1_1_mat_job, pi_vi_vj_job, &
            pi_v2_v2, pi_1_1_mat, pi_vi_vj, verbose = verbose)

        if ( proc_id == root_process ) then

            call save_self_energies(filename, pi_v2_v2, pi_1_1_mat, pi_vi_vj, verbose = verbose)

        end if

    end subroutine

    subroutine calc_abs_Pi_spin(proc_id, root_process, &
            tran_form_1_job, tran_form_v_job, tran_form_v2_job, &
            pi_v2_v2, pi_vi_vj, pi_1_1_mat, ik_manager, job_id_to_ik, &
            PW_dataset, widths, dm_model, target_mat, filename, verbose)
        !! Compute the self energies with spin dependent wave functions.
        implicit none

        integer :: proc_id
        integer :: root_process
        complex(dp) :: tran_form_1_job(:, :, :, :)
        complex(dp) :: tran_form_v_job(:, :, :, :, :)
        complex(dp) :: tran_form_v2_job(:, :, :, :)
        complex(dp) :: pi_v2_v2(:, :)
        complex(dp) :: pi_vi_vj(:, :, :, :)
        complex(dp) :: pi_1_1_mat(:, :, :, :)
        type(parallel_manager_t) :: ik_manager
        integer :: job_id_to_ik(:, :)
        type(PW_dataset_t) :: PW_dataset
        type(width_parameters_t) :: widths
        type(dm_model_t) :: dm_model
        type(material_t) :: target_mat
        character(len=*) :: filename
        logical, optional :: verbose

        complex(dp), allocatable :: pi_v2_v2_job(:, :, :)
            !! Dim : [n_omega, n_widths, n_tran_per_proc]
            !! 
            !! self energy with two v^2 insertions
            !!
            !! Units : eV^2

        complex(dp), allocatable :: pi_1_1_mat_job(:, :, :, :, :)
            !! Dim : [3, 3, n_omega, n_widths, n_tran_per_proc] 
            !!
            !! self energy, without q_vec's
            !!
            !! Units : eV^2

        complex(dp), allocatable :: pi_vi_vj_job(:, :, :, :, :)
            !! Dim : [3, 3, n_omega, n_widths, n_tran_per_proc] 
            !!
            !! Units : eV^2

        integer :: j, job_id, val_id, k

        if ( verbose ) then

            print*, 'Computing self energies...'
            print*

        end if

        ! allocate job variables
        allocate(pi_v2_v2_job(size(pi_v2_v2, 1), size(pi_v2_v2, 2), ik_manager%n_jobs_per_proc))
        pi_v2_v2_job = (0.0_dp, 0.0_dp)
        allocate(pi_vi_vj_job(3, 3, size(pi_v2_v2, 1), size(pi_v2_v2, 2), ik_manager%n_jobs_per_proc))
        pi_vi_vj_job = (0.0_dp, 0.0_dp)
        allocate(pi_1_1_mat_job(3, 3, size(pi_v2_v2, 1), size(pi_v2_v2, 2), ik_manager%n_jobs_per_proc))
        pi_1_1_mat_job = (0.0_dp, 0.0_dp)

        do j = 1, ik_manager%n_jobs_per_proc

            job_id = ik_manager%job_table(proc_id + 1, j)

            if ( job_id /= 0 ) then

                val_id = job_id_to_ik(job_id, 1)
                k = job_id_to_ik(job_id, 2)

                call calc_pi_v2_v2(pi_v2_v2_job(:, :, j), tran_form_v2_job(j, :, :, :), &
                    dm_model, widths, PW_dataset, target_mat, val_id, k)
                call calc_pi_1_1_mat(pi_1_1_mat_job(:, :, :, :, j), tran_form_v_job(:, j, :, :, :), &
                    dm_model, widths, PW_dataset, target_mat, val_id, k)
                call calc_pi_vi_vj(pi_vi_vj_job(:, :, :, :, j), tran_form_v_job(:, j, :, :, :), &
                    dm_model, widths, PW_dataset, target_mat, val_id, k)

            end if

        end do

        if ( verbose ) then

            print*, 'Done computing self energies!'
            print*

        end if

        call ik_manager%comm_self_energies(proc_id, root_process, &
            pi_v2_v2_job, pi_1_1_mat_job, pi_vi_vj_job, &
            pi_v2_v2, pi_1_1_mat, pi_vi_vj, verbose = verbose)

        if ( proc_id == root_process ) then

            call save_self_energies(filename, pi_v2_v2, pi_1_1_mat, pi_vi_vj, verbose = verbose)

        end if

    end subroutine

    subroutine calc_pi_v2_v2_no_spin(pi_v2_v2, tran_form_v2, &
            dm_model, widths, PW_dataset, target_mat, val_id, k)
        !* Compute \( \Pi_{\bar v^2, \bar v^2} \) from spin independent wave functions.
        !
        ! Note: this calculation has some subtleties due to having to take care of divergent terms. Specifically, it can be shown
        ! that \( \mathrm{Im} \Pi_{\bar{v}^2, \bar{v}^2} \sim \omega^2 \) as \( \omega \rightarrow 0 \) (for the gapped targets of
        ! interest here). See discussion in calculation of \( \Pi_{v_i, v_j} \) for more details.
        !
        ! [TODO] Check if this procedure is valid for \( \mathrm{Re} \Pi_{\bar{v}^2, \bar{v}^2} \) or 
        ! \( \mathrm{Re}\Pi_{\bar{v}^2, \bar{v}^2} + \) (constant term).

        implicit none

        complex(dp) :: pi_v2_v2(:, :)
        complex(dp) :: tran_form_v2(:)
        type(dm_model_t) :: dm_model
        type(width_parameters_t) :: widths
        type(PW_dataset_t) :: PW_dataset
        type(material_t) :: target_mat
        integer :: val_id, k

        integer :: w, p, f, cond_id

        real(dp) :: omega, width, omega_if

        do p = 1, widths%n
            do w = 1, dm_model%n_mX

                omega = dm_model%mX(w)
                width = widths%get_width(p, omega)

                do f = 1, size(tran_form_v2, 1)

                    cond_id = f + PW_dataset%n_val

                    omega_if = PW_dataset%energy_bands(k, cond_id) - PW_dataset%energy_bands(k, val_id)

                    if ( abs(omega - omega_if) < widths%sigma*width ) then

                        pi_v2_v2(w, p) = pi_v2_v2(w, p) + &
                           (target_mat%pc_vol)**(-1)*PW_dataset%k_weight(k)*(0.25_dp)*&
                           ( omega / omega_if )**2*&
                           (&
                                ( omega - omega_if + ii*width )**(-1) - &
                                ( omega + omega_if - ii*width )**(-1) &
                            )*&
                           abs(tran_form_v2(f))**2*(PW_dataset%spin_degen/2.0_dp)

                    end if

                end do

            end do
        end do

    end subroutine

    subroutine calc_pi_v2_v2_spin(pi_v2_v2, tran_form_v2, &
            dm_model, widths, PW_dataset, target_mat, val_id, k)
        !* Compute \( \Pi_{\bar v^2, \bar v^2} \) from spin dependent wave functions.
        !
        ! Note: this calculation has some subtleties due to having to take care of divergent terms. Specifically, it can be shown
        ! that \( \mathrm{Im} \Pi_{\bar{v}^2, \bar{v}^2} \sim \omega^2 \) as \( \omega \rightarrow 0 \) (for the gapped targets of
        ! interest here). See discussion in calculation of \( \Pi_{v_i, v_j} \) for more details.
        !
        ! [TODO] Check if this procedure is valid for \( \mathrm{Re} \Pi_{\bar{v}^2, \bar{v}^2} \) or 
        ! \( \mathrm{Re}\Pi_{\bar{v}^2, \bar{v}^2} + \) (constant term).

        implicit none

        complex(dp) :: pi_v2_v2(:, :)
        complex(dp) :: tran_form_v2(:, :, :)
        type(dm_model_t) :: dm_model
        type(width_parameters_t) :: widths
        type(PW_dataset_t) :: PW_dataset
        type(material_t) :: target_mat
        integer :: val_id, k

        integer :: w, p, f, cond_id, s

        real(dp) :: omega, width, omega_if

        complex(dp) :: tf

        do p = 1, widths%n
            do w = 1, dm_model%n_mX

                omega = dm_model%mX(w)
                width = widths%get_width(p, omega)

                do f = 1, size(tran_form_v2, 1)

                    tf = (0.0_dp, 0.0_dp)

                    do s = 1, 2
                        tf = tf + tran_form_v2(f, s, s) 
                    end do

                    cond_id = f + PW_dataset%n_val

                    omega_if = PW_dataset%energy_bands(k, cond_id) - PW_dataset%energy_bands(k, val_id)

                    if ( abs(omega - omega_if) < widths%sigma*width ) then

                        pi_v2_v2(w, p) = pi_v2_v2(w, p) + &
                           (target_mat%pc_vol)**(-1)*PW_dataset%k_weight(k)*(0.25_dp)*&
                           ( omega / omega_if )**2*&
                           (&
                                ( omega - omega_if + ii*width )**(-1) - &
                                ( omega + omega_if - ii*width )**(-1) &
                            )*&
                           abs(tf)**2*(PW_dataset%spin_degen/2.0_dp)

                    end if

                end do

            end do
        end do

    end subroutine

    subroutine calc_pi_1_1_mat_no_spin(pi_1_1_mat, tran_form_v, &
            dm_model, widths, PW_dataset, target_mat, val_id, k)
        !! Compute \( \mathbf{\Pi}_{11} \) from spin independent wave functions.

        implicit none

        complex(dp) :: pi_1_1_mat(:, :, :, :)
        complex(dp) :: tran_form_v(:, :)
        type(dm_model_t) :: dm_model
        type(width_parameters_t) :: widths
        type(PW_dataset_t) :: PW_dataset
        type(material_t) :: target_mat
        integer :: val_id, k

        integer :: w, p, f, cond_id, i, j

        real(dp) :: omega, width, omega_if

        complex(dp) :: outer_ff(3, 3, size(tran_form_v, 2))

        do f = 1, size(tran_form_v, 2)

            cond_id = f + PW_dataset%n_val

            do i = 1, 3
                do j = 1, 3

                    omega_if = PW_dataset%energy_bands(k, cond_id) - PW_dataset%energy_bands(k, val_id)

                    outer_ff(i, j, f) = (m_elec/omega_if)**2*conjg(tran_form_v(i, f))*tran_form_v(j, f)

                end do
            end do
        end do

        do p = 1, widths%n
            do w = 1, dm_model%n_mX

                omega = dm_model%mX(w)
                width = widths%get_width(p, omega)

                do f = 1, size(tran_form_v, 2)

                    cond_id = f + PW_dataset%n_val

                    omega_if = PW_dataset%energy_bands(k, cond_id) - PW_dataset%energy_bands(k, val_id)

                    if ( abs(omega - omega_if) < widths%sigma*width ) then

                        pi_1_1_mat(:, :, w, p) = pi_1_1_mat(:, :, w, p) + &
                           (target_mat%pc_vol)**(-1)*PW_dataset%k_weight(k)*&
                           ( &
                                ( omega - omega_if + ii*width )**(-1) - &
                                ( omega + omega_if - ii*width )**(-1) &
                            )*&
                           outer_ff(:, :, f)*(PW_dataset%spin_degen/2.0_dp)

                    end if

                end do

            end do
        end do

    end subroutine

    subroutine calc_pi_1_1_mat_spin(pi_1_1_mat, tran_form_v, &
            dm_model, widths, PW_dataset, target_mat, val_id, k)
        !! Compute \( \mathbf{\Pi}_{11} \) from spin dependent wave functions.

        implicit none

        complex(dp) :: pi_1_1_mat(:, :, :, :)
        complex(dp) :: tran_form_v(:, :, :, :)
        type(dm_model_t) :: dm_model
        type(width_parameters_t) :: widths
        type(PW_dataset_t) :: PW_dataset
        type(material_t) :: target_mat
        integer :: val_id, k

        integer :: w, p, f, cond_id, s, i, j

        real(dp) :: omega, width, omega_if

        complex(dp) :: outer_ff(3, 3, size(tran_form_v, 2))

        complex(dp) :: tf_v(3)

        do f = 1, size(tran_form_v, 2)

            tf_v = (0.0_dp, 0.0_dp)

            do s = 1, 2

                tf_v = tf_v + tran_form_v(:, f, s, s)

            end do

            cond_id = f + PW_dataset%n_val

            do i = 1, 3
                do j = 1, 3

                    omega_if = PW_dataset%energy_bands(k, cond_id) - PW_dataset%energy_bands(k, val_id)

                    outer_ff(i, j, f) = (m_elec/omega_if)**2*conjg(tf_v(i))*tf_v(j)

                end do
            end do
        end do

        do p = 1, widths%n
            do w = 1, dm_model%n_mX

                omega = dm_model%mX(w)
                width = widths%get_width(p, omega)

                do f = 1, size(tran_form_v, 2)

                    cond_id = f + PW_dataset%n_val

                    omega_if = PW_dataset%energy_bands(k, cond_id) - PW_dataset%energy_bands(k, val_id)

                    if ( abs(omega - omega_if) < widths%sigma*width ) then

                        pi_1_1_mat(:, :, w, p) = pi_1_1_mat(:, :, w, p) + &
                           (target_mat%pc_vol)**(-1)*PW_dataset%k_weight(k)*&
                           ( &
                                ( omega - omega_if + ii*width )**(-1) - &
                                ( omega + omega_if - ii*width )**(-1) &
                            )*&
                           outer_ff(:, :, f)*(PW_dataset%spin_degen/2.0_dp)

                    end if

                end do

            end do
        end do

    end subroutine

    subroutine calc_pi_vi_vj_no_spin(pi_vi_vj, tran_form_v, &
            dm_model, widths, PW_dataset, target_mat, val_id, k)
        !* Compute \( \Pi_{v^i, v^j} \) from spin independent wave functions.
        ! 
        ! Note: this calculation has some subtleties due to having to take care of divergent terms. Specifically, it can be shown
        ! that \( \mathrm{Im} \Pi_{v_i, v_j} \sim \omega^2 \) as \( \omega \rightarrow 0 \) (for the gapped targets of interest
        ! here). If one naively computes \( \Pi_{v_i, v_j} \) numerically you find that \( \Pi_{v_i, v_j} \sim \delta \) as \(
        ! \omega \rightarow 0 \). While alone this does not look problematic since \( \Pi_{v_i, v_j} \) itself is converging, the
        ! calculation of other physical quantities, such as the dielectric \( \propto \Pi_{v_i, v_j}/\omega^2 \), will diverge.
        !
        ! This is remedied by explicitly only computing the first non-divergent term, which can be shown, see 
        ! https://journals.aps.org/prb/pdf/10.1103/PhysRevB.48.11705, to be related to the naive calculation with the addition of a
        ! \( ( \omega / \omega_{ii'\mathbf{k}})^2 \) term inside the 1BZ sum. We follow that procedure here.
        !
        ! [TODO] Check if this procedure is valid for \( \mathrm{Re} \Pi_{v_i, v_j} \) or \( \mathrm{Re} \Pi_{v_i, v_j} + \) 
        ! (constant term).

        implicit none

        complex(dp) :: pi_vi_vj(:, :, :, :)
        complex(dp) :: tran_form_v(:, :)
        type(dm_model_t) :: dm_model
        type(width_parameters_t) :: widths
        type(PW_dataset_t) :: PW_dataset
        type(material_t) :: target_mat
        integer :: val_id, k

        integer :: w, p, f, cond_id, i, j

        real(dp) :: omega, width, omega_if

        complex(dp) :: outer_ff(3, 3, size(tran_form_v, 2))

        do f = 1, size(tran_form_v, 2)

            cond_id = f + PW_dataset%n_val

            do i = 1, 3
                do j = 1, 3

                    omega_if = PW_dataset%energy_bands(k, cond_id) - PW_dataset%energy_bands(k, val_id)

                    outer_ff(i, j, f) = conjg(tran_form_v(i, f))*tran_form_v(j, f)

                end do
            end do
        end do

        do p = 1, widths%n
            do w = 1, dm_model%n_mX

                omega = dm_model%mX(w)
                width = widths%get_width(p, omega)

                do f = 1, size(tran_form_v, 2)

                    cond_id = f + PW_dataset%n_val

                    omega_if = PW_dataset%energy_bands(k, cond_id) - PW_dataset%energy_bands(k, val_id)

                    if ( abs(omega - omega_if) < widths%sigma*width ) then

                        pi_vi_vj(:, :, w, p) = pi_vi_vj(:, :, w, p) + &
                           (target_mat%pc_vol)**(-1)*PW_dataset%k_weight(k)*&
                           ( omega / omega_if )**2*&
                           ( &
                                ( omega - omega_if + ii*width )**(-1) - &
                                ( omega + omega_if - ii*width )**(-1) &
                            )*&
                           outer_ff(:, :, f)*(PW_dataset%spin_degen/2.0_dp)

                    end if

                end do

            end do
        end do

    end subroutine

    subroutine calc_pi_vi_vj_spin(pi_vi_vj, tran_form_v, &
            dm_model, widths, PW_dataset, target_mat, val_id, k)
        !* Compute \( \Pi_{v^i, v^j} \) from spin dependent wave functions.
        !
        ! Note: this calculation has some subtleties due to having to take care of divergent terms. Specifically, it can be shown
        ! that \( \mathrm{Im} \Pi_{v_i, v_j} \sim \omega^2 \) as \( \omega \rightarrow 0 \) (for the gapped targets of interest
        ! here). If one naively computes \( \Pi_{v_i, v_j} \) numerically you find that \( \Pi_{v_i, v_j} \sim \delta \) as \(
        ! \omega \rightarow 0 \). While alone this does not look problematic since \( \Pi_{v_i, v_j} \) itself is converging, the
        ! calculation of other physical quantities, such as the dielectric \( \propto \Pi_{v_i, v_j}/\omega^2 \), will diverge.
        !
        ! This is remedied by explicitly only computing the first non-divergent term, which can be shown, see 
        ! https://journals.aps.org/prb/pdf/10.1103/PhysRevB.48.11705, to be related to the naive calculation with the addition of a
        ! \( ( \omega / \omega_{ii'\mathbf{k}})^2 \) term inside the 1BZ sum. We follow that procedure here.
        !
        ! [TODO] Check if this procedure is valid for \( \mathrm{Re} \Pi_{v_i, v_j} \) or \( \mathrm{Re} \Pi_{v_i, v_j} + \) 
        ! (constant term).

        implicit none

        complex(dp) :: pi_vi_vj(:, :, :, :)
        complex(dp) :: tran_form_v(:, :, :, :)
        type(dm_model_t) :: dm_model
        type(width_parameters_t) :: widths
        type(PW_dataset_t) :: PW_dataset
        type(material_t) :: target_mat
        integer :: val_id, k

        integer :: w, p, f, cond_id, s, i, j

        real(dp) :: omega, width, omega_if

        complex(dp) :: outer_ff(3, 3, size(tran_form_v, 2))

        complex(dp) :: tf_v(3)

        do f = 1, size(tran_form_v, 2)

            tf_v = (0.0_dp, 0.0_dp)

            do s = 1, 2

                tf_v = tf_v + tran_form_v(:, f, s, s)

            end do

            cond_id = f + PW_dataset%n_val

            do i = 1, 3
                do j = 1, 3

                    omega_if = PW_dataset%energy_bands(k, cond_id) - PW_dataset%energy_bands(k, val_id)

                    outer_ff(i, j, f) = conjg(tf_v(i))*tf_v(j)

                end do
            end do
        end do

        do p = 1, widths%n
            do w = 1, dm_model%n_mX

                omega = dm_model%mX(w)
                width = widths%get_width(p, omega)

                do f = 1, size(tran_form_v, 2)

                    cond_id = f + PW_dataset%n_val

                    omega_if = PW_dataset%energy_bands(k, cond_id) - PW_dataset%energy_bands(k, val_id)

                    if ( abs(omega - omega_if) < widths%sigma*width ) then

                        pi_vi_vj(:, :, w, p) = pi_vi_vj(:, :, w, p) + &
                           (target_mat%pc_vol)**(-1)*PW_dataset%k_weight(k)*&
                           ( omega / omega_if )**2*&
                           ( &
                                ( omega - omega_if + ii*width )**(-1) - &
                                ( omega + omega_if - ii*width )**(-1) &
                            )*&
                           outer_ff(:, :, f)*(PW_dataset%spin_degen/2.0_dp)

                    end if

                end do

            end do
        end do

    end subroutine

    subroutine save_self_energies(filename, &
            pi_v2_v2, pi_1_1_mat, pi_vi_vj, verbose)
        !! Saves the self energies, \( \Pi_{\mathcal{O}_1, \mathcal{O}_2} \).

        implicit none

        character(len=*) :: filename
        complex(dp) :: pi_v2_v2(:, :)
        complex(dp) :: pi_vi_vj(:, :, :, :)
        complex(dp) :: pi_1_1_mat(:, :, :, :)
        logical, optional :: verbose

        integer(HID_T) :: file_id
        integer(HID_T) :: group_id

        logical :: file_exists

        integer(HSIZE_T) :: dims1(1) = [1]
        integer(HSIZE_T) :: dims2(2)
        integer(HSIZE_T) :: dims3(3)

        integer :: error

        integer :: m, f, t
        integer :: i, fin

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

            do i = 1, size(pi_v2_v2, 2)

                call h5gcreate_f(file_id,& 
                    'self_energies'//&
                    '/width_'//trim(adjustl(int_to_str(i))),&
                    group_id, error)

                dims1 = [size(pi_v2_v2, 1)]

                call h5ltmake_dataset_double_f(file_id,&
                    'self_energies'//&
                    '/width_'//trim(adjustl(int_to_str(i)))//&
                    '/pi_v2_v2_r', &
                    size(dims1), dims1,&
                    real(pi_v2_v2(:, i)), error)

                call h5ltmake_dataset_double_f(file_id,&
                    'self_energies'//&
                    '/width_'//trim(adjustl(int_to_str(i)))//&
                    '/pi_v2_v2_c', &
                    size(dims1), dims1,&
                    aimag(pi_v2_v2(:, i)), error)

                dims3 = [3, 3, size(pi_v2_v2, 1)]

                call h5ltmake_dataset_double_f(file_id,&
                    'self_energies'//&
                    '/width_'//trim(adjustl(int_to_str(i)))//&
                    '/pi_1_1_mat_r', &
                    size(dims3), dims3,&
                    real(pi_1_1_mat(:, :, :, i)), error)

                call h5ltmake_dataset_double_f(file_id,&
                    'self_energies'//&
                    '/width_'//trim(adjustl(int_to_str(i)))//&
                    '/pi_1_1_mat_c', &
                    size(dims3), dims3,&
                    aimag(pi_1_1_mat(:, :, :, i)), error)

                call h5ltmake_dataset_double_f(file_id,&
                    'self_energies'//&
                    '/width_'//trim(adjustl(int_to_str(i)))//&
                    '/pi_vi_vj_r', &
                    size(dims3), dims3,&
                    real(pi_vi_vj(:, :, :, i)), error)

                call h5ltmake_dataset_double_f(file_id,&
                    'self_energies'//&
                    '/width_'//trim(adjustl(int_to_str(i)))//&
                    '/pi_vi_vj_c', &
                    size(dims3), dims3,&
                    aimag(pi_vi_vj(:, :, :, i)), error)

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
