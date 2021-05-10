module dme_scatter_cf
    !! Compute the scattering rate from core to free states

    use prec
    use math_mod

    use numerics_input
    use material_input
    use particle_physics_scatter
    use core_electrons
    use transition
    use Zeff_input

    use calc_dme_scatter_cf

contains

    subroutine run_dme_scatter_cf(binned_rate_t, n_tran_per_proc, &
            sto_wf_filename, core_elec_config_filename, nml_input_filename, &
            out_filename, proc_id, root_process, verbose)

        implicit none

        real(dp) :: binned_rate_t(:, :, :, :, :, :)

        integer :: n_tran_per_proc

        character(len=*) :: sto_wf_filename
        character(len=*) :: core_elec_config_filename
        character(len=*) :: nml_input_filename
        character(len=*) :: out_filename

        integer :: proc_id, root_process

        logical, optional :: verbose

        integer :: tran_id, init_id, fin_id
        integer :: t, i, f

        real(dp) :: log_omega_min, log_omega_max
        real(dp) :: log_omega
        real(dp) :: log_omegas(2)

        real(dp), allocatable :: log_omega_list(:)
        real(dp), allocatable :: log_omega_table(:, :)

        real(dp) :: ki_angular_mesh(n_ki_theta*n_ki_phi, 2)
        real(dp) :: kf_angular_mesh(n_kf_theta*n_kf_phi, 2)

        ! calculation setup
        
        call load_core_elec_config(trim(core_elec_config_filename), verbose=verbose)
        call load_core_sto_data(trim(sto_wf_filename), verbose=verbose)

        call load_Zeff_parameters(trim(nml_input_filename), verbose = verbose)

        ki_angular_mesh = generate_uniform_points_on_sphere(n_ki_theta, n_ki_phi)
        kf_angular_mesh = generate_uniform_points_on_sphere(n_kf_theta, n_kf_phi)

        ! do calculation

        ! create array of omega's that each processor will compute the transition rate
        ! for
        allocate(log_omega_table(n_fin, 2))
        allocate(log_omega_list(n_fin + 1))

        log_omega_min = log10(Ef_max)
        log_omega_max = log10(0.5_dp*v_max**2*maxval(mX))

        if ( log_omega_max > log_omega_min ) then

            do f = 1, n_fin + 1

                log_omega = log_omega_min + (log_omega_max - log_omega_min)*&
                    (f - 1.0_dp)/max(n_fin + 1.0_dp - 1.0_dp, 1.0_dp)

                log_omega_list(f) = log_omega

            end do

            do f = 1, n_fin
                do i = 0, 1

                    log_omega_table(f, i + 1) = log_omega_list(f + i)

                end do
            end do

            if ( verbose ) then

                print*, 'Calculating rate...'
                print*

            end if

            do t = 1, n_tran_per_proc

                tran_id = job_table(proc_id + 1, t)

                if ( tran_id .ne. 0 ) then

                    init_id = tran_to_init_fin_id(tran_id, 1)
                    fin_id = tran_to_init_fin_id(tran_id, 2)

                    log_omegas = log_omega_table(fin_id, :)

                    call dme_scatter_cf_calc(binned_rate_t(:, :, :, :, :, t),& 
                        init_id, log_omegas, n_ki, ki_angular_mesh, kf_angular_mesh, verbose = verbose)

                end if

            end do

        else

            if ( verbose ) then

                print*, '!!! WARNING !!!'
                print*
                print*, '   No mass in mX is large enough to warrant the c -> f calculation'
                print*, '   with the specified Ef_max. Skipping c -> f calculation.'
                print*
                print*, '!!!!!!!!!!!!!!!'
                print*

            end if

        end if

        if ( verbose ) then
            print*, '----------'
            print*
        end if

        if ( proc_id == root_process ) then

            call save_core_electrons(out_filename, verbose = verbose)

        end if

    end subroutine

end module
