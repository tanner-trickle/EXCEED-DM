module material_type

    use prec_util

    implicit none

    type :: material_t
        ! .. note: 
        !   .. math:
        !       \mathbf{b}_1 = \frac{2 \pi}{\Omega} \mathbf{a}_1 \times \mathbf{a}_2 

        character(len=512) :: name
        character(len=512) :: materials_project_ID

        real(dp) :: pc_vol_A
        real(dp) :: pc_vol

        real(dp) :: rho_T_g_per_cm3
        real(dp) :: rho_T

        real(dp) :: band_gap

        real(dp) :: a_vecs_Ang(3, 3)
        real(dp) :: a_vecs(3, 3)

        real(dp) :: b_vecs(3, 3)

        real(dp) :: k_red_to_xyz(3, 3) ! Convert momentum from reduced to XYZ coordinates
        real(dp) :: k_xyz_to_red(3, 3) ! Convert momentum from XYZ to reduced coordinates
        real(dp) :: red_to_xyz(3, 3)

        contains

            procedure :: set_defaults => material_type_set_defaults
            procedure :: get_values => material_type_get_values
            procedure :: save => material_type_save

    end type

contains

    subroutine material_type_save(self, filename)

        use hdf5_utils

        implicit none

        class(material_t) :: self
        character(len=*) :: filename

        integer(HID_T) :: file_id

        call hdf_open_file(file_id, filename, status='OLD', action='WRITE')

        call hdf_create_group(file_id, 'material')

        call hdf_write_dataset(file_id, 'material/name', self%name)
        call hdf_write_dataset(file_id, 'material/pc_vol', self%pc_vol_A)
        call hdf_write_dataset(file_id, 'material/rho_T', self%rho_T_g_per_cm3)
        call hdf_write_dataset(file_id, 'material/band_gap', self%band_gap)

        call hdf_close_file(file_id)

    end subroutine

    subroutine material_type_set_defaults(self, cfg)

        use m_config

        implicit none

        class(material_t) :: self
        type(CFG_t) :: cfg

        call CFG_add(cfg, &
                     "material%name", &
                     "", &
                     "Name of the target material")

        call CFG_add(cfg, &
                     "material%materials_project_ID", &
                     "", &
                     "Materials Project ID of the target material")

        call CFG_add(cfg, &
                     "material%rho_T_g_per_cm3", &
                     1.0_dp, &
                     "Mass density of the target material, \( \rho_T \)<br />Units: \( \mathrm{g}/\mathrm{cm}^3 \)")

        call CFG_add(cfg, &
                     "material%band_gap", &
                     0.0_dp, &
                     "Band gap of the target material, \( E_g \)<br />Units: \( \mathrm{g}/\mathrm{cm}^3 \)")

        call CFG_add(cfg, &
                     "material%a_vecs_Ang", &
                     [ 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp ], &
                     "Lattice vectors of the target material, \( \mathbf{a}_i \). Each row is a different lattice vector"//&
                     "<br />Units : \( \mathrm{\AA} \)<br />Dim : [3, 3]", &
                    dynamic_size = .TRUE.)

    end subroutine

    subroutine material_type_get_values(self, cfg)

        use m_config

        use constants_util
        use math_util

        implicit none

        class(material_t) :: self
        type(CFG_t) :: cfg

        real(dp) :: a_vecs_list(9)

        complex(dp) :: a_vec_mat(3, 3)
        complex(dp) :: eigs(3)
        real(dp) :: det

        call CFG_get(cfg, "material%name", self%name)
        call CFG_get(cfg, "material%materials_project_ID", self%materials_project_ID)

        call CFG_get(cfg, "material%rho_T_g_per_cm3", self%rho_T_g_per_cm3)
        self%rho_T = self%rho_T_g_per_cm3*inv_cm_to_eV**3*g_to_eV

        call CFG_get(cfg, "material%band_gap", self%band_gap)

        call CFG_get(cfg, "material%a_vecs_Ang", a_vecs_list)
        self%a_vecs_Ang = transpose(reshape(a_vecs_list, [ 3, 3 ]))

        self%a_vecs = self%a_vecs_Ang*Ang_to_inv_eV

        self%red_to_xyz = transpose(self%a_vecs)

        a_vec_mat = (1.0_dp, 0.0_dp)*self%a_vecs_Ang

        ! compute primitive cell volume
        eigs = calc_eigvals_33(a_vec_mat)
        det = abs(eigs(1)*eigs(2)*eigs(3))

        self%pc_vol_A = det
        self%pc_vol = self%pc_vol_A*Ang_to_inv_eV**3

        ! Compute b_vec's
        self%b_vecs(1, :) = (2*pi/self%pc_vol)*cross_product(self%a_vecs(2, :), self%a_vecs(3, :))
        self%b_vecs(2, :) = (2*pi/self%pc_vol)*cross_product(self%a_vecs(3, :), self%a_vecs(1, :))
        self%b_vecs(3, :) = (2*pi/self%pc_vol)*cross_product(self%a_vecs(1, :), self%a_vecs(2, :))

        self%k_red_to_xyz = transpose(self%b_vecs)

        self%k_xyz_to_red = transpose(self%red_to_xyz)/(2.0_dp*pi)

    end subroutine

end module
