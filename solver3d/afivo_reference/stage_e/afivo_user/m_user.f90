! Stage E1 benchmark-specific Afivo user hook.
!
! This file is copied temporarily over programs/standard_3d/m_user.f90 for E1
! electrostatic runs, then the pinned checkout is restored. It defines a
! non-axisymmetric triangular-foil electrode and compact diagnostics only; it
! does not modify Afivo core src/.
module m_user
#include "../../afivo/src/cpp_macros.h"
  use m_af_all
  use m_config
  use m_field, only: field_get_E_vector
  use m_streamer
  use m_types
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

    user_lsf => e1_triangular_foil_lsf
    user_log_variables => e1_log_variables
  end subroutine user_initialize

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
    type(af_loc_t) :: loc_field
    real(dp) :: r_field(3)

    call af_tree_max_cc(tree, i_electric_fld, emax, loc_field)
    if (loc_field%id > 0) then
       r_field = af_r_loc(tree, loc_field)
    else
       r_field = 0.0_dp
    end if

    vol_high_05 = 0.0_dp
    vol_high_08 = 0.0_dp
    conductor_vol = 0.0_dp
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
                   if (lsf <= 0.0_dp) then
                      conductor_vol = conductor_vol + vol
                   else if (emax > 0.0_dp) then
                      if (eabs > 0.5_dp * emax) vol_high_05 = vol_high_05 + vol
                      if (eabs > 0.8_dp * emax) vol_high_08 = vol_high_08 + vol
                   end if
                end do
             end do
          end do
       end do
    end do

    if (e1_write_field_export) call e1_write_field_cells(tree)

    n_vars = 7
    var_names(1) = "e1_vol_Egt_0p5"
    var_names(2) = "e1_vol_Egt_0p8"
    var_names(3) = "e1_conductor_vol"
    var_names(4) = "e1_Emax_x"
    var_names(5) = "e1_Emax_y"
    var_names(6) = "e1_Emax_z"
    var_names(7) = "e1_min_dx"
    var_values(1) = vol_high_05
    var_values(2) = vol_high_08
    var_values(3) = conductor_vol
    var_values(4) = r_field(1)
    var_values(5) = r_field(2)
    var_values(6) = r_field(3)
    var_values(7) = af_min_dr(tree)
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

end module m_user
