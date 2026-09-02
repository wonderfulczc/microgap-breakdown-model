! Stage E1 benchmark-specific Afivo user hook.
!
! This file is copied temporarily over programs/standard_3d/m_user.f90 for E1
! electrostatic runs, then the pinned checkout is restored. It defines a
! non-axisymmetric triangular-foil electrode and compact diagnostics only; it
! does not modify Afivo core src/.
module m_user
#include "../../afivo/src/cpp_macros.h"
  use m_af_all
  use m_chemistry
  use m_config
  use m_field, only: field_get_E_vector
  use m_gas
  use m_lookup_table
  use m_streamer
  use m_transport_data
  use m_types
  use m_units_constants
  use m_user_methods

  implicit none
  private

  public :: user_initialize

  real(dp) :: e1_foil_thickness = 10.0e-6_dp
  real(dp) :: e1_foil_width = 40.0e-6_dp
  real(dp) :: e1_foil_length = 50.0e-6_dp
  real(dp) :: e1_tip_radius = 3.0e-6_dp
  real(dp) :: e1_edge_radius = 1.0e-6_dp
  real(dp) :: e1_gap = 70.0e-6_dp
  real(dp) :: e1_ground_z = 0.0_dp
  logical :: e1_write_field_export = .false.
  character(len=512) :: e1_field_export_prefix = "e1_triangular_foil"
  integer, save :: e1_export_counter = 0
  character(len=32) :: e2_mode = "electrostatic"
  real(dp) :: e2_seed_n0 = 1.0e16_dp
  real(dp) :: e2_seed_sigma = 3.0e-6_dp
  real(dp) :: e2_seed_center(3) = [0.0_dp, 0.0_dp, 60.0e-6_dp]
  real(dp) :: e2_head_threshold = 1.0e14_dp
  logical :: e2_write_source_export = .false.
  character(len=512) :: e2_source_export_prefix = "e2_triangular_foil"
  integer :: e2_source_export_stride = 1
  integer, save :: e2_source_export_counter = 0
  integer, save :: e2_source_call_counter = 0
  integer, save :: e2_i_neg = -1

contains

  subroutine user_initialize(cfg, tree)
    type(CFG_t), intent(inout) :: cfg
    type(af_t), intent(inout) :: tree

    call CFG_add_get(cfg, "e1%foil_thickness", e1_foil_thickness, &
         "Stage E1 foil thickness in y direction")
    call CFG_add_get(cfg, "e1%foil_width", e1_foil_width, &
         "Stage E1 full foil width in x direction")
    call CFG_add_get(cfg, "e1%foil_length", e1_foil_length, &
         "Stage E1 triangular foil length in z direction")
    call CFG_add_get(cfg, "e1%tip_radius", e1_tip_radius, &
         "Stage E1 finite in-plane tip rounding radius")
    call CFG_add_get(cfg, "e1%edge_radius", e1_edge_radius, &
         "Stage E1 finite transverse edge rounding radius")
    call CFG_add_get(cfg, "e1%gap", e1_gap, &
         "Stage E1 gap from ground plane to foil tip")
    call CFG_add_get(cfg, "e1%ground_z", e1_ground_z, &
         "Stage E1 grounded plane z location")
    call CFG_add_get(cfg, "e1%field_export_prefix", e1_field_export_prefix, &
         "Stage E1 field-cell CSV export prefix")
    call CFG_add_get(cfg, "e1%write_field_export", e1_write_field_export, &
         "Write Stage E1 compact field-cell CSV export")
    call CFG_add_get(cfg, "e2%mode", e2_mode, &
         "Stage E2 mode: electrostatic or dynamic_gaussian")
    call CFG_add_get(cfg, "e2%seed_n0", e2_seed_n0, &
         "Stage E2 Gaussian seed peak density")
    call CFG_add_get(cfg, "e2%seed_sigma", e2_seed_sigma, &
         "Stage E2 Gaussian seed sigma")
    call CFG_add_get(cfg, "e2%seed_center", e2_seed_center, &
         "Stage E2 Gaussian seed center")
    call CFG_add_get(cfg, "e2%head_threshold", e2_head_threshold, &
         "Stage E2 electron-density threshold for leading edge")
    call CFG_add_get(cfg, "e2%source_export_prefix", e2_source_export_prefix, &
         "Stage E2 source snapshot CSV export prefix")
    call CFG_add_get(cfg, "e2%source_export_stride", e2_source_export_stride, &
         "Write one Stage E2 source snapshot every this many output calls")
    call CFG_add_get(cfg, "e2%write_source_export", e2_write_source_export, &
         "Write Stage E2 rho/ne/E/J source snapshots")

    user_lsf => e1_triangular_foil_lsf
    if (trim(e2_mode) == "dynamic_gaussian") then
       user_initial_conditions => e2_initial_conditions
    end if
    user_log_variables => e1_log_variables
  end subroutine user_initialize

  integer function e2_negative_ion_index()
    if (e2_i_neg < 1) then
       e2_i_neg = species_itree(species_index("N_min"))
    end if
    e2_negative_ion_index = e2_i_neg
  end function e2_negative_ion_index

  subroutine e2_initial_conditions(box)
    type(box_t), intent(inout) :: box
    integer :: i, j, k, nc
    integer :: ijk(3)
    real(dp) :: rr(3), dens, r2

    if (e2_negative_ion_index() < 1) error stop "Stage E2 requires N_min species"
    nc = box%n_cell
    do k = 0, nc+1
       do j = 0, nc+1
          do i = 0, nc+1
             ijk = [i, j, k]
             rr = af_r_cc(box, ijk)
             if (e1_triangular_foil_lsf(rr) <= 0.0_dp) then
                dens = 0.0_dp
             else
                r2 = (rr(1) - e2_seed_center(1))**2 + &
                     (rr(2) - e2_seed_center(2))**2 + &
                     (rr(3) - e2_seed_center(3))**2
                dens = e2_seed_n0 * exp(-r2 / (2.0_dp * e2_seed_sigma**2))
             end if
             box%cc(i, j, k, i_electron) = dens
             box%cc(i, j, k, i_1pos_ion) = dens
             box%cc(i, j, k, e2_negative_ion_index()) = 0.0_dp
          end do
       end do
    end do
  end subroutine e2_initial_conditions

  real(dp) function e1_triangular_foil_lsf(rr)
    real(dp), intent(in) :: rr(NDIM)
    real(dp) :: phi_xz, phi_y, half_thickness_core

    phi_xz = e1_signed_distance_rounded_triangle(rr(1), rr(3))
    half_thickness_core = max(0.5_dp * e1_foil_thickness - e1_edge_radius, &
         0.1_dp * e1_foil_thickness)
    phi_y = abs(rr(2)) - half_thickness_core
    e1_triangular_foil_lsf = max(phi_xz, phi_y) - e1_edge_radius
  end function e1_triangular_foil_lsf

  real(dp) function e1_signed_distance_rounded_triangle(x, z)
    real(dp), intent(in) :: x, z
    real(dp) :: ax, az, bx, bz, cx, cz, sd, dist
    logical :: inside

    ax = 0.0_dp
    az = e1_ground_z + e1_gap + e1_tip_radius - e1_edge_radius
    bx = 0.5_dp * e1_foil_width - e1_edge_radius
    bz = e1_ground_z + e1_gap + e1_foil_length - e1_tip_radius
    cx = -0.5_dp * e1_foil_width + e1_edge_radius
    cz = bz

    dist = min(e1_dist_segment(x, z, ax, az, bx, bz), &
         min(e1_dist_segment(x, z, bx, bz, cx, cz), &
         e1_dist_segment(x, z, cx, cz, ax, az)))
    inside = e1_same_side_inside(x, z, ax, az, bx, bz, cx, cz)

    if (inside) then
       sd = -dist
    else
       sd = dist
    end if
    e1_signed_distance_rounded_triangle = sd - e1_tip_radius + e1_edge_radius
  end function e1_signed_distance_rounded_triangle

  real(dp) function e1_dist_segment(px, pz, ax, az, bx, bz)
    real(dp), intent(in) :: px, pz, ax, az, bx, bz
    real(dp) :: vx, vz, wx, wz, denom, t, qx, qz

    vx = bx - ax
    vz = bz - az
    wx = px - ax
    wz = pz - az
    denom = max(vx * vx + vz * vz, tiny(1.0_dp))
    t = max(0.0_dp, min(1.0_dp, (wx * vx + wz * vz) / denom))
    qx = ax + t * vx
    qz = az + t * vz
    e1_dist_segment = sqrt((px - qx)**2 + (pz - qz)**2)
  end function e1_dist_segment

  logical function e1_same_side_inside(px, pz, ax, az, bx, bz, cx, cz)
    real(dp), intent(in) :: px, pz, ax, az, bx, bz, cx, cz
    real(dp) :: c1, c2, c3

    c1 = e1_cross2(bx - ax, bz - az, px - ax, pz - az)
    c2 = e1_cross2(cx - bx, cz - bz, px - bx, pz - bz)
    c3 = e1_cross2(ax - cx, az - cz, px - cx, pz - cz)
    e1_same_side_inside = (c1 >= 0.0_dp .and. c2 >= 0.0_dp .and. &
         c3 >= 0.0_dp) .or. (c1 <= 0.0_dp .and. c2 <= 0.0_dp .and. &
         c3 <= 0.0_dp)
  end function e1_same_side_inside

  pure real(dp) function e1_cross2(ax, az, bx, bz)
    real(dp), intent(in) :: ax, az, bx, bz
    e1_cross2 = ax * bz - az * bx
  end function e1_cross2

  subroutine e1_log_variables(tree, n_vars, var_names, var_values)
    type(af_t), intent(in) :: tree
    integer, intent(out) :: n_vars
    character(len=name_len), intent(inout) :: var_names(user_max_log_vars)
    real(dp), intent(inout) :: var_values(user_max_log_vars)
    integer :: lvl, n, id, nc, i, j, k
    real(dp) :: vol, emax, eabs, lsf
    real(dp) :: vol_high_05, vol_high_08, conductor_vol
    real(dp) :: ne, np, nn, charge_number, total_charge_C
    real(dp) :: total_electrons, sum_x, sum_y, sum_z, sum_r2, sum_x2_minus_y2
    real(dp) :: sum_rho_abs, sum_rho_x2_minus_y2, min_ne, head_z
    type(af_loc_t) :: loc_field
    real(dp) :: r_field(3)
    real(dp) :: rr(3)
    integer :: ijk(3)

    call af_tree_max_cc(tree, i_electric_fld, emax, loc_field)
    if (loc_field%id > 0) then
       r_field = af_r_loc(tree, loc_field)
    else
       r_field = 0.0_dp
    end if

    vol_high_05 = 0.0_dp
    vol_high_08 = 0.0_dp
    conductor_vol = 0.0_dp
    total_charge_C = 0.0_dp
    total_electrons = 0.0_dp
    sum_x = 0.0_dp
    sum_y = 0.0_dp
    sum_z = 0.0_dp
    sum_r2 = 0.0_dp
    sum_x2_minus_y2 = 0.0_dp
    sum_rho_abs = 0.0_dp
    sum_rho_x2_minus_y2 = 0.0_dp
    min_ne = huge(1.0_dp)
    head_z = huge(1.0_dp)
    do lvl = 1, tree%highest_lvl
       vol = product(af_lvl_dr(tree, lvl))
       do n = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(n)
          nc = tree%boxes(id)%n_cell
          do k = 1, nc
             do j = 1, nc
                do i = 1, nc
                   eabs = tree%boxes(id)%cc(i, j, k, i_electric_fld)
                   lsf = tree%boxes(id)%cc(i, j, k, i_lsf)
                   ne = tree%boxes(id)%cc(i, j, k, i_electron)
                   np = tree%boxes(id)%cc(i, j, k, i_1pos_ion)
                   if (e2_negative_ion_index() > 0) then
                      nn = tree%boxes(id)%cc(i, j, k, e2_negative_ion_index())
                   else
                      nn = 0.0_dp
                   end if
                   charge_number = np - ne - nn
                   total_charge_C = total_charge_C + UC_elem_charge * charge_number * vol
                   total_electrons = total_electrons + ne * vol
                   min_ne = min(min_ne, ne)
                   if (lsf <= 0.0_dp) then
                      conductor_vol = conductor_vol + vol
                   else if (emax > 0.0_dp) then
                      if (eabs > 0.5_dp * emax) vol_high_05 = vol_high_05 + vol
                      if (eabs > 0.8_dp * emax) vol_high_08 = vol_high_08 + vol
                   end if
                   ijk = [i, j, k]
                   rr = af_r_cc(tree%boxes(id), ijk)
                   if (ne > 0.0_dp) then
                      sum_x = sum_x + rr(1) * ne * vol
                      sum_y = sum_y + rr(2) * ne * vol
                      sum_z = sum_z + rr(3) * ne * vol
                      sum_r2 = sum_r2 + (rr(1)**2 + rr(2)**2) * ne * vol
                      sum_x2_minus_y2 = sum_x2_minus_y2 + &
                           (rr(1)**2 - rr(2)**2) * ne * vol
                   end if
                   if (abs(charge_number) > 0.0_dp) then
                      sum_rho_abs = sum_rho_abs + abs(charge_number) * vol
                      sum_rho_x2_minus_y2 = sum_rho_x2_minus_y2 + &
                           abs(charge_number) * (rr(1)**2 - rr(2)**2) * vol
                   end if
                   if (ne >= e2_head_threshold) head_z = min(head_z, rr(3))
                end do
             end do
          end do
       end do
    end do

    if (e1_write_field_export) call e1_write_field_cells(tree)
    if (e2_write_source_export) then
       if (mod(e2_source_call_counter, max(e2_source_export_stride, 1)) == 0) &
            call e2_write_source_cells(tree)
       e2_source_call_counter = e2_source_call_counter + 1
    end if

    n_vars = 17
    var_names(1) = "e1_vol_Egt_0p5"
    var_names(2) = "e1_vol_Egt_0p8"
    var_names(3) = "e1_conductor_vol"
    var_names(4) = "e1_Emax_x"
    var_names(5) = "e1_Emax_y"
    var_names(6) = "e1_Emax_z"
    var_names(7) = "e1_min_dx"
    var_names(8) = "e2_x_cm"
    var_names(9) = "e2_y_cm"
    var_names(10) = "e2_r_cm"
    var_names(11) = "e2_head_z"
    var_names(12) = "e2_total_charge_C"
    var_names(13) = "e2_min_ne"
    var_names(14) = "e2_ne_moment_asym"
    var_names(15) = "e2_rho_moment_asym"
    var_names(16) = "e2_total_electrons"
    var_names(17) = "e3_z_cm"
    var_values(1) = vol_high_05
    var_values(2) = vol_high_08
    var_values(3) = conductor_vol
    var_values(4) = r_field(1)
    var_values(5) = r_field(2)
    var_values(6) = r_field(3)
    var_values(7) = af_min_dr(tree)
    if (total_electrons > 0.0_dp) then
       var_values(8) = sum_x / total_electrons
       var_values(9) = sum_y / total_electrons
       var_values(10) = sqrt((sum_x / total_electrons)**2 + &
            (sum_y / total_electrons)**2)
       var_values(14) = abs(sum_x2_minus_y2) / max(sum_r2, 1.0e-300_dp)
       var_values(17) = sum_z / total_electrons
    else
       var_values(8:10) = 0.0_dp
       var_values(14) = 0.0_dp
       var_values(17) = 0.0_dp
    end if
    if (head_z < huge(1.0_dp) / 10.0_dp) then
       var_values(11) = head_z
    else
       var_values(11) = -huge(1.0_dp) / 10.0_dp
    end if
    var_values(12) = total_charge_C
    var_values(13) = min_ne
    if (sum_rho_abs > 0.0_dp) then
       var_values(15) = abs(sum_rho_x2_minus_y2) / &
            max(sum_rho_abs * e1_gap**2, 1.0e-300_dp)
    else
       var_values(15) = 0.0_dp
    end if
    var_values(16) = total_electrons
  end subroutine e1_log_variables

  subroutine e1_write_field_cells(tree)
    type(af_t), intent(in) :: tree
    integer :: unit, lvl, n, id
    character(len=1024) :: path

    write(path, "(A,'_field_cells_',I6.6,'.csv')") &
         trim(e1_field_export_prefix), e1_export_counter
    open(newunit=unit, file=trim(path), status="replace", action="write")
    write(unit, "(A)") "x_m,y_m,z_m,phi_V,Ex_Vpm,Ey_Vpm,Ez_Vpm,Eabs_Vpm,lsf_m,level"
    do lvl = 1, tree%highest_lvl
       do n = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(n)
          call e1_write_box_cells(unit, tree%boxes(id), lvl)
       end do
    end do
    close(unit)
    e1_export_counter = e1_export_counter + 1
  end subroutine e1_write_field_cells

  subroutine e1_write_box_cells(unit, box, lvl)
    integer, intent(in) :: unit, lvl
    type(box_t), intent(in) :: box
    integer :: i, j, k, nc
    integer :: ijk(3)
    real(dp) :: rr(3)
    real(dp) :: evec(DTIMES(1:box%n_cell), NDIM)

    nc = box%n_cell
    evec = field_get_E_vector(box)
    do k = 1, nc
       do j = 1, nc
          do i = 1, nc
             ijk = [i, j, k]
             rr = af_r_cc(box, ijk)
             write(unit, "(ES25.16E3,8(',',ES25.16E3),',',I0)") rr(1), rr(2), rr(3), &
                  box%cc(i, j, k, i_phi), evec(i, j, k, 1), &
                  evec(i, j, k, 2), evec(i, j, k, 3), &
                  box%cc(i, j, k, i_electric_fld), &
                  box%cc(i, j, k, i_lsf), lvl
          end do
       end do
    end do
  end subroutine e1_write_box_cells

  subroutine e2_write_source_cells(tree)
    type(af_t), intent(in) :: tree
    integer :: unit, lvl, n, id
    character(len=1024) :: path

    write(path, "(A,'_source_',I6.6,'.csv')") &
         trim(e2_source_export_prefix), e2_source_export_counter
    open(newunit=unit, file=trim(path), status="replace", action="write")
    write(unit, "(A)") "time_s,x_m,y_m,z_m,cell_volume_m3,rho_Cpm3,ne_m3," // &
         "Ex_Vpm,Ey_Vpm,Ez_Vpm,Jx_Apm2,Jy_Apm2,Jz_Apm2," // &
         "Jrf_x_Apm2,Jrf_y_Apm2,Jrf_z_Apm2,Eabs_Vpm,lsf_m,level"
    do lvl = 1, tree%highest_lvl
       do n = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(n)
          call e2_write_box_source_cells(unit, tree%boxes(id), lvl, &
               product(af_lvl_dr(tree, lvl)))
       end do
    end do
    close(unit)
    e2_source_export_counter = e2_source_export_counter + 1
  end subroutine e2_write_source_cells

  subroutine e2_write_box_source_cells(unit, box, lvl, cell_vol)
    integer, intent(in) :: unit, lvl
    type(box_t), intent(in) :: box
    real(dp), intent(in) :: cell_vol
    integer :: i, j, k, nc
    integer :: ijk(3)
    real(dp) :: rr(3)
    real(dp) :: evec(DTIMES(1:box%n_cell), NDIM)
    real(dp) :: ne, np, nn, rho, eabs, Td, sigma
    real(dp) :: jx, jy, jz
    real(dp) :: jrf_x, jrf_y, jrf_z

    nc = box%n_cell
    evec = field_get_E_vector(box)
    do k = 1, nc
       do j = 1, nc
          do i = 1, nc
             ijk = [i, j, k]
             rr = af_r_cc(box, ijk)
             ne = box%cc(i, j, k, i_electron)
             np = box%cc(i, j, k, i_1pos_ion)
             nn = box%cc(i, j, k, e2_negative_ion_index())
             rho = UC_elem_charge * (np - ne - nn)
             eabs = box%cc(i, j, k, i_electric_fld)
             Td = SI_to_Townsend * eabs * gas_inverse_number_density
             sigma = LT_get_col(td_tbl, td_mobility, Td) * &
                  gas_inverse_number_density * max(ne, 0.0_dp) * UC_elem_charge
             jx = sigma * evec(i, j, k, 1)
             jy = sigma * evec(i, j, k, 2)
             jz = sigma * evec(i, j, k, 3)
             jrf_x = -0.5_dp * UC_elem_charge * &
                  (box%fc(i, j, k, 1, flux_elec) + box%fc(i + 1, j, k, 1, flux_elec))
             jrf_y = -0.5_dp * UC_elem_charge * &
                  (box%fc(i, j, k, 2, flux_elec) + box%fc(i, j + 1, k, 2, flux_elec))
             jrf_z = -0.5_dp * UC_elem_charge * &
                  (box%fc(i, j, k, 3, flux_elec) + box%fc(i, j, k + 1, 3, flux_elec))
             write(unit, "(18(ES25.16E3,','),I0)") global_time, rr(1), rr(2), rr(3), &
                  cell_vol, rho, ne, evec(i, j, k, 1), evec(i, j, k, 2), &
                  evec(i, j, k, 3), jx, jy, jz, jrf_x, jrf_y, jrf_z, &
                  eabs, box%cc(i, j, k, i_lsf), lvl
          end do
       end do
    end do
  end subroutine e2_write_box_source_cells

end module m_user
