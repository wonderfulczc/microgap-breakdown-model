# Stage C-R2 Reaction Source Decomposition

This diagnostic exposes the electron reaction source already used by the
frozen Stage C solver. It does not change transport, chemistry,
photoionization, timestep control, boundary conditions, current definitions,
or `Gb/Rb`.

## Existing Electron Source Equation

For each gas cell, `evaluate_reactions()` assembles the Stage C reaction
model as

```text
S_e =
  S_impact
  + S_photo
  - L_attach2
  - L_attach3
  - L_recomb_ep
```

with units `m^-3 s^-1`.

The code mapping is:

- `S_impact = ionization_frequency * ne`
- `L_attach2 = attachment_two_body_frequency * ne`
- `L_attach3 = attachment_three_body_frequency * ne`
- `L_recomb_ep = beta_ep(mobility, diffusion) * ne * np`
- `S_photo = sph`

`positive-negative` ion recombination is present in the ion equations, but it
does not enter the electron equation.

## Sign Convention

The component fields named `*_source` are positive production magnitudes.
The component fields named `*_loss` are positive loss magnitudes. The net
electron reaction source is signed:

```text
net_electron_reaction_source =
  impact_source
  + photo_source
  - attachment2_loss
  - attachment3_loss
  - electron_recombination_loss
```

## Numerical Stage

The decomposition is `PRE_LIMITER`. It represents the physical reaction model
before the discrete positivity check and negative-density clamp in
`StreamerSolver::step()`. The positivity protection remains unchanged.

## Volume-Integrated Rates

When enabled, the diagnostic integrates gas-cell volumetric rates using the
axisymmetric cell volume:

```text
dV = 2 pi r dr dz
```

The output rates have units `s^-1`:

- `impact_rate_s_1`
- `photo_rate_s_1`
- `attach2_rate_s_1`
- `attach3_rate_s_1`
- `recomb_e_rate_s_1`
- `net_electron_reaction_rate_s_1`

The source closure is

```text
source_closure =
  impact_rate
  + photo_rate
  - attach2_rate
  - attach3_rate
  - recomb_e_rate
  - net_electron_reaction_rate
```

This is only a chemical-source closure. It excludes electron transport
divergence, electrode absorption, outer-boundary loss, displacement current,
and total particle-number continuity.

## SP3 Semantics

With SP3 off, `sph` is set to zero and the photoionization component is zero.
With SP3 on, `sph` is the existing Stage C SP3 solution masked to gas cells.
The decomposition reads this existing field; it does not solve another
photoionization problem.

## Output

`stage_c2_dynamic --source-decomposition` writes
`reaction_source_diagnostics.csv`. The option is off by default, so baseline
Stage C output and solver state are unchanged unless the diagnostic is
explicitly requested.

## Scope

C-R2 does not classify discharge stages, identify streamer heads, or make a
publication-level physical interpretation. Those analyses require later
diagnostic layers.
