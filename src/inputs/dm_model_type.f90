module dm_model_type
    ! Defines the :code:`dm_model_t` type which contains all of the inputs regarding the DM model.

    use prec_util, only: dp

    implicit none

    type :: dm_model_t

        real(dp), allocatable :: mX(:)
            ! DM masses
            !
            ! Units : eV
        real(dp), allocatable :: med_FF(:)

        real(dp) :: rho_X_GeV_per_cm3
        real(dp) :: rho_X

        character(len=512) :: particle_type
        character(len=512) :: FIF_id

        contains

            procedure :: set_defaults => dm_model_type_set_defaults
            procedure :: get_values => dm_model_type_get_values
            procedure :: save => dm_model_type_save

    end type

contains

    subroutine dm_model_type_save(self, filename)

        use hdf5_utils

        implicit none

        class(dm_model_t) :: self
        character(len=*) :: filename

        integer(HID_T) :: file_id

        call hdf_open_file(file_id, filename, status='OLD', action='WRITE')

        call hdf_create_group(file_id, 'dm_model')

        call hdf_write_dataset(file_id, 'dm_model/mX', self%mX)
        call hdf_write_dataset(file_id, 'dm_model/med_FF', self%med_FF)
        call hdf_write_dataset(file_id, 'dm_model/rho_X', self%rho_X_GeV_per_cm3)
        call hdf_write_dataset(file_id, 'dm_model/particle_type', self%particle_type)
        call hdf_write_dataset(file_id, 'dm_model/FIF_id', self%FIF_id)

        call hdf_close_file(file_id)

    end subroutine

    subroutine dm_model_type_set_defaults(self, cfg)
        !* Defines the default values inside the dm_model
        use m_config

        implicit none

        class(dm_model_t) :: self
        type(CFG_t) :: cfg

        call CFG_add(cfg,                       &
                     "dm_model%mX",             &
                     [ 0.0_dp ],              &
                     "Dark matter masses, $m_\chi$<br />"//&
                     "Units : eV<br />"//&
                     "Dim : [:]", &
                     dynamic_size = .TRUE.)

        call CFG_add(cfg,                       &
                     "dm_model%med_FF",             &
                     [ 0.0_dp ],              &
                     "Mediator form factor powers, \( \beta \)<br />"//&
                     "Formula: \( \mathcal{F}_\text{med} = \left( \frac{\alpha m_e}{q} \right)^\beta \)<br />"//&
                     "Dim : [:]", &
                     dynamic_size = .TRUE.)

        call CFG_add(cfg,                       &
                     "dm_model%rho_X_GeV_per_cm3",             &
                     0.4_dp,              &
                     "Local dark matter density (GeV/cm^3)")

        call CFG_add(cfg,                       &
                     "dm_model%particle_type",             &
                     'fermion',              &
                     "Type of incoming dark matter particle")

        call CFG_add(cfg,                       &
                     "dm_model%FIF_id",             &
                     'SI',              &
                     "Form factor ID to compute for")

        call CFG_add(cfg,                       &
                     "dm_model%mX_logspace",             &
                     [ 0.0_dp, 1.0_dp, 1.0_dp ],              &
                     "Add logarithmically spaced dark matter masses. [ N, mX_min, mX_max] (eV)")

        call CFG_add(cfg,                       &
                     "dm_model%mX_linspace",             &
                     [ 0.0_dp, 1.0_dp, 1.0_dp ],              &
                     "Add linearly spaced dark matter masses. [ N, mX_min, mX_max] (eV)")

    end subroutine

    subroutine dm_model_type_get_values(self, cfg)

        use m_config

        use constants_util

        implicit none

        class(dm_model_t) :: self
        type(CFG_t) :: cfg

        integer :: n, i

        real(dp), allocatable :: extra_mX(:)

        real(dp) :: mX_logspace(3)
        real(dp) :: mX_linspace(3)

        call CFG_get_size(cfg, "dm_model%mX", n)
        call CFG_get(cfg, "dm_model%mX_logspace", mX_logspace)
        call CFG_get(cfg, "dm_model%mX_linspace", mX_linspace)

        allocate(self%mX(n + int(mX_logspace(1)) + int(mX_linspace(1))))

        call CFG_get(cfg, "dm_model%mX", self%mX(:n))

        do i = 1, int(mX_logspace(1))
            self%mX(n + i) = 10.0_dp**(&
                log10(mX_logspace(2)) + &
                ( log10(mX_logspace(3)) - log10(mX_logspace(2)) )*(i - 1.0_dp)/&
                max(1, int(mX_logspace(1)) - 1)&
                )
        end do

        do i = 1, int(mX_linspace(1))
            self%mX(n + int(mX_logspace(1)) + i) = mX_linspace(2) + &
                ( mX_linspace(3) - mX_linspace(2) )*(i - 1.0_dp)/&
                max(1, int(mX_linspace(1)) - 1)
        end do

        ! Get rid of 0 or negative values
        self%mX = pack(self%mX, self%mX > 1.0e-8_dp)

        call CFG_get_size(cfg, "dm_model%med_FF", n)
        allocate(self%med_FF(n))
        call CFG_get(cfg, "dm_model%med_FF", self%med_FF)

        call CFG_get(cfg, "dm_model%rho_X_GeV_per_cm3", self%rho_X_GeV_per_cm3)
        self%rho_X = self%rho_X_GeV_per_cm3*inv_cm_to_eV**3*10**9

        call CFG_get(cfg, "dm_model%particle_type", self%particle_type)
        call CFG_get(cfg, "dm_model%FIF_id", self%FIF_id)

    end subroutine

end module
