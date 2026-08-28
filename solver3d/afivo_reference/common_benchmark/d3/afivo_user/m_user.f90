! Stage D3 benchmark-specific Afivo user hook.
! This file is copied temporarily over programs/standard_3d/m_user.f90 for D3
! runs, then the pinned checkout is restored. It does not modify Afivo core src/.
module m_user
  use m_af_all
  use m_config
  use m_user_methods
  use m_types

  implicit none
  private

  public :: user_initialize

  character(len=32) :: d3_mode = "none"
  real(dp) :: d3_seed_n0 = 1.0e16_dp
  real(dp) :: d3_seed_sigma = 3.0e-6_dp
  real(dp) :: d3_seed_z0 = 65.0e-6_dp
  real(dp), parameter :: d3_ground_z = 0.0_dp
  real(dp), parameter :: d3_rod_radius = 5.0e-6_dp
  real(dp), parameter :: d3_rod_cap_z = 75.0e-6_dp
  real(dp) :: d3_homogeneous_ne = 0.0_dp
  real(dp) :: d3_homogeneous_np = 0.0_dp
  real(dp) :: d3_homogeneous_nn = 0.0_dp

contains

  subroutine user_initialize(cfg, tree)
    type(CFG_t), intent(inout) :: cfg
    type(af_t), intent(inout) :: tree

    call CFG_add_get(cfg, "d3%mode", d3_mode, &
         "Stage D3 benchmark mode: none, gaussian, homogeneous")
    call CFG_add_get(cfg, "d3%seed_n0", d3_seed_n0, &
         "Stage D3 Gaussian seed peak density")
    call CFG_add_get(cfg, "d3%seed_sigma", d3_seed_sigma, &
         "Stage D3 common Gaussian sigma")
    call CFG_add_get(cfg, "d3%seed_z0", d3_seed_z0, &
         "Stage D3 Gaussian seed z center")
    call CFG_add_get(cfg, "d3%homogeneous_ne", d3_homogeneous_ne, &
         "Stage D3 homogeneous electron density")
    call CFG_add_get(cfg, "d3%homogeneous_np", d3_homogeneous_np, &
         "Stage D3 homogeneous positive ion density")
    call CFG_add_get(cfg, "d3%homogeneous_nn", d3_homogeneous_nn, &
         "Stage D3 homogeneous negative ion density")

    if (trim(d3_mode) == "homogeneous") then
       user_initial_conditions => d3_initial_conditions
    end if
    user_log_variables => d3_log_variables
  end subroutine user_initialize

  subroutine d3_initial_conditions(box)
    use m_chemistry
    use m_streamer
    type(box_t), intent(inout) :: box
    integer :: i, j, k, nc, i_neg
    integer :: ijk(3)
    real(dp) :: rr(3), dens, r2

    nc = box%n_cell
    i_neg = species_itree(species_index("N_min"))
    if (i_neg < 1) error stop "Stage D3 requires N_min species"
    do k = 0, nc+1
       do j = 0, nc+1
          do i = 0, nc+1
             if (trim(d3_mode) == "gaussian") then
                ijk = [i, j, k]
                rr = af_r_cc(box, ijk)
                if (d3_inside_conductor(rr)) then
                   dens = 0.0_dp
                else
                   r2 = rr(1)**2 + rr(2)**2 + (rr(3) - d3_seed_z0)**2
                   dens = d3_seed_n0 * exp(-r2 / (2.0_dp * d3_seed_sigma**2))
                end if
                box%cc(i, j, k, i_electron) = dens
                box%cc(i, j, k, i_1pos_ion) = dens
                box%cc(i, j, k, i_neg) = 0.0_dp
             else if (trim(d3_mode) == "homogeneous") then
                box%cc(i, j, k, i_electron) = d3_homogeneous_ne
                box%cc(i, j, k, i_1pos_ion) = d3_homogeneous_np
                box%cc(i, j, k, i_neg) = d3_homogeneous_nn
             end if
          end do
       end do
    end do
  end subroutine d3_initial_conditions

  logical function d3_inside_conductor(rr)
    real(dp), intent(in) :: rr(3)
    real(dp) :: rho, cap_dist

    rho = sqrt(rr(1)**2 + rr(2)**2)
    cap_dist = sqrt(rho**2 + (rr(3) - d3_rod_cap_z)**2)
    d3_inside_conductor = rr(3) <= d3_ground_z .or. &
         cap_dist <= d3_rod_radius .or. &
         (rho <= d3_rod_radius .and. rr(3) >= d3_rod_cap_z)
  end function d3_inside_conductor

  subroutine d3_log_variables(tree, n_vars, var_names, var_values)
    use m_streamer
    type(af_t), intent(in) :: tree
    integer, intent(out) :: n_vars
    character(len=name_len), intent(inout) :: var_names(user_max_log_vars)
    real(dp), intent(inout) :: var_values(user_max_log_vars)
    integer :: lvl, n, id, nc, i, j, k
    integer :: ijk(3)
    real(dp) :: vol, ne, rr(3)
    real(dp) :: sum_ne, sum_x, sum_y, sum_r2, sum_x2_minus_y2

    sum_ne = 0.0_dp
    sum_x = 0.0_dp
    sum_y = 0.0_dp
    sum_r2 = 0.0_dp
    sum_x2_minus_y2 = 0.0_dp

    do lvl = 1, tree%highest_lvl
       vol = product(af_lvl_dr(tree, lvl))
       do n = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(n)
          nc = tree%boxes(id)%n_cell
          do k = 1, nc
             do j = 1, nc
                do i = 1, nc
                   ijk = [i, j, k]
                   rr = af_r_cc(tree%boxes(id), ijk)
                   ne = tree%boxes(id)%cc(i, j, k, i_electron)
                   sum_ne = sum_ne + ne * vol
                   sum_x = sum_x + rr(1) * ne * vol
                   sum_y = sum_y + rr(2) * ne * vol
                   sum_r2 = sum_r2 + (rr(1)**2 + rr(2)**2) * ne * vol
                   sum_x2_minus_y2 = sum_x2_minus_y2 + (rr(1)**2 - rr(2)**2) * ne * vol
                end do
             end do
          end do
       end do
    end do

    n_vars = 4
    var_names(1) = "d3_x_cm"
    var_names(2) = "d3_y_cm"
    var_names(3) = "d3_xy_moment_asym"
    var_names(4) = "d3_r_rms"
    if (sum_ne > 0.0_dp) then
       var_values(1) = sum_x / sum_ne
       var_values(2) = sum_y / sum_ne
       var_values(3) = abs(sum_x2_minus_y2) / max(sum_r2, 1.0e-300_dp)
       var_values(4) = sqrt(max(sum_r2 / sum_ne, 0.0_dp))
    else
       var_values(1:4) = 0.0_dp
    end if
  end subroutine d3_log_variables

end module m_user
