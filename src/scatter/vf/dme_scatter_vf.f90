module dme_scatter_vf
    !! Compute the scattering rate from valence to free states 

    use prec
    use math_mod

    use numerics_input
    use material_input
    use DFT_parameters
    use particle_physics_scatter
    use transition
    use Zeff_input

    use calc_dme_scatter_vf

    implicit none

contains

    subroutine run_dme_scatter_vf(binned_rate_t, n_tran_per_proc, DFT_input_filename,&
            nml_input_filename, out_filename, proc_id, root_process, verbose)

        implicit none

        real(dp) :: binned_rate_t(:, :, :, :, :, :)

        logical, optional :: verbose

        integer :: proc_id, root_process

        integer :: n_tran_per_proc
        integer :: t, f, i

        integer :: tran_id
        integer :: val_id, fin_id

        character(len=*) :: DFT_input_filename
        character(len=*) :: out_filename
        character(len=*) :: nml_input_filename

        complex(dp), allocatable :: wfc_FT_i(:, :)

        real(dp) :: angular_mesh(n_kf_theta*n_kf_phi, 2)

        real(dp), allocatable :: log_omega_table(:, :)
        real(dp), allocatable :: log_omega_list(:)
        real(dp) :: log_omega_min, log_omega_max
        real(dp) :: log_omega
        real(dp) :: log_omegas(2)

        ! calculation setup
        call load_DFT_parameters(trim(DFT_input_filename), verbose = verbose)
        call do_scissor_correction(band_gap, verbose = verbose)

        call load_Zeff_parameters(trim(nml_input_filename), verbose = verbose)

        angular_mesh = generate_uniform_points_on_sphere(n_kf_theta, n_kf_phi)

        allocate(wfc_FT_i(n_k, n_in_G))

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

                    val_id = tran_to_init_fin_id(tran_id, 1)
                    fin_id = tran_to_init_fin_id(tran_id, 2)

                    call get_in_wfc_FT(DFT_input_filename, val_id, wfc_FT_i)

                    log_omegas = log_omega_table(fin_id, :)

                    call dme_scatter_vf_calc(binned_rate_t(:, :, :, :, :, t),& 
                        wfc_FT_i, val_id, log_omegas, n_k, angular_mesh, verbose = verbose)

                end if

            end do

        else

            if ( verbose ) then

                print*, '!!! WARNING !!!'
                print*
                print*, '   No mass in mX is large enough to warrant the v -> f calculation'
                print*, '   with the specified Ef_max. Skipping v -> f calculation.'
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

            call save_DFT_parameters(out_filename, verbose = verbose)

        end if

    end subroutine

end module
