# Stage 2 Formula Audit

## Reopen implementation audit addendum (2026-07-22)

- OpenCharge uses the axisymmetric ring kernel with `std::comp_ellint_1(sqrt(m))`: the C++ function takes modulus `k`, while the formula constructs parameter `m=k²`. Gaussian cell quadrature is the baseline and `MPI_Allreduce` assembles global charge contributions.
- The boundary is recomputed from the supplied charge field before every potential solve. Tests cover on-axis and off-axis rings, Gaussian volume charge, threshold response, and changed-charge response.
- The Liu 2007 Robin cross coefficients are nonzero and enter the boundary right-hand side. The normal derivative uses a second-order one-sided stencil; sequential relaxed fixed-point solves stop only after measured boundary-residual convergence.
- Former unconditional Robin checks and scripted observations are invalid evidence and are excluded by the evidence audit.

This audit was completed before implementation. Kulikovsky (1995) is a seven-page scan; every page was rendered at 1.5× resolution and inspected visually. Text extraction was used only as a navigation aid for the born-digital papers.

| equation_id | reference and PDF page | original equation | program symbol and sign | SI units | limit / implementation | automated test | ambiguity |
|---|---|---|---|---|---|---|---|
| S2-EQ-001 | Kulikovsky 1995 p.1 | (2),(4) | `sg_flux`; converted from electron current convention to positive-direction particle flux `W n-D dn/dx` | m^-2 s^-1 | Bernoulli series, `expm1`, overflow branches; `cpp/src/transport.cpp` | Bernoulli, diffusion, upwind, antisymmetry | paper's `j` direction differs from project particle-flux convention |
| S2-EQ-002 | Kulikovsky 1995 pp.3-4 | (17),(21) | `h_virtual=sqrt(2 epsilon D h/abs(W_R-W_L))` | m | SG fallback for small gradient or `h_v>=h` | ISG regular/fallback | epsilon is numerical, not physical |
| S2-EQ-003 | Kulikovsky 1995 p.4 | (22),(23) | virtual positions and velocities | m; m/s | symmetric about face | ISG regular branch | none |
| S2-EQ-004 | Kulikovsky 1995 p.3 | (20) | `u=n/n_ref`; interpolate `log1p(u)` | dimensionless | explicit positive `n_ref`, scaling test | n_ref scaling | paper's `n+1` is dimensionless; interface resolves dimensional ambiguity |
| S2-EQ-005 | Kulikovsky 1995 pp.4-5 | (24),(31) | ISG-0 flux and zero-width limit | m^-2 s^-1 | zero-width evaluated separately | zero-width finite | printed formula uses electron-current notation |
| S2-EQ-006 | Kulikovsky 1995 pp.5-6 | (32)-(36), Figs.3-5 | linear-field shock and Gaussian profiles | paper nondimensional units | 100/200 nodes; CFL and epsilon sweeps | Python regression | no pixel-level reproduction claimed |
| S2-EQ-007 | standard axisymmetric calculus; Shi 2016 p.3 and Shi 2017 p.3 lineage | Poisson | `(1/r)d_r(r d_r phi)+d_zz phi=-rho/eps0` | V/m2 | exact annular cell volumes; axis flux zero | grid volume, PETSc solve | boundary details are Stage-2 choices |
| S2-EQ-008 | electrostatic azimuthal integration | ring Green function | complete elliptic integral `K(k)` | V | axis analytic limit and `k^2<1` guard | ring potential, elliptic kernel | near-singular boundary cell requires quadrature control |
| S2-EQ-009 | Bourdon 2007 PDF p.3 | (1)-(5) | Zheleznyak direct integral | m^-3 s^-1 | Python small-grid reference | kernel-fit comparison | emission source is supplied, chemistry excluded |
| S2-EQ-010 | Bourdon 2007 PDF p.8 | (14), Table 3 | `A=[.0067,.0346,.3059]`, `lambda=[.0447,.1121,.5994]` | original cm^-1 Torr^-1; matrix m^-1 | multiply by `p_O2` and 100 cm/m | coefficient conversion | explicitly not Helmholtz Table 2 |
| S2-EQ-011 | Bourdon 2007 PDF p.9 | (16)-(19) | `gamma`, `kappa^2`, two fields per group | SI | computed from radicals, six elliptic solves | SP3 constants and solve | Ségur misprints corrected by Bourdon |
| S2-EQ-012 | Bourdon 2007 PDF p.9 | (20) | `S_ph=sum(a_j psi_0j)` | m^-3 s^-1 | Table-3 group weights | finite SP3 solution | input already contains emission efficiency |
| S2-EQ-013 | Liu et al. 2007 PDF p.2 | (1),(2) | coupled outward-normal Robin conditions | m^-1 | fixed-point coupling; no-reflection/no-emission | boundary sensitivity registry | extracted glyphs cross-checked against rendered equations |
| S2-EQ-014 | Liu and Pasko 2004 PDF pp.5-7 | (6) and transport discussion | integral reference and future transport table | mixed originals; SI on ingest | PCHIP/no extrapolation decision | table range test | digitization deferred to Stage 3 |
| S2-EQ-015 | Shi 2016 PDF pp.3-4; Shi 2017 PDF pp.2-3 | model lineage | Poisson and SP3 module use | SI | scope reference only | integration interface test | not a Case-I implementation |
| S2-EQ-016 | Ihaddadene & Celestin 2015 PDF pp.2-4 | collision model context | no implemented Stage-2 equation | n/a | exclusion audit | forbidden-output validator | collision is outside Stage 2 |
| S2-EQ-017 | Shi 2019 PDF pp.2-4 | model and radiation chain | no implemented Stage-2 radiation equation | n/a | exclusion audit | forbidden-output validator | Case I and current moment excluded |

The local Luque (2017) paper was also checked as full-Maxwell/collision context; it does not authorize any Stage 2 FDTD implementation. No batch OCR formula was accepted as sole evidence.
