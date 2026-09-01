# Stage F3 RF Spectral and Trusted-Bandwidth Framework

Stage F3 analyzes observer waveforms. It does not attribute spectral bands to
discharge stages and does not make VHF/UHF/GHz claims from Stage E smoke data.

## Fourier Convention

The implemented convention is

`E_tilde(f) = integral E(t) exp(-i 2*pi*f*t) dt`

The discrete implementation uses `dt * rfft(E)` and keeps the internal complex
spectrum. One-sided spectra are used for real-signal reporting. Metadata always
records `dt`, `N`, `T`, `df`, and Nyquist frequency.

## Windowing

Supported windows:

- `rectangular`
- `hann`

The window metadata records coherent gain and energy correction. Band-energy
and fluence integrations use Parseval-consistent one-sided ESD normalization.
Hann-windowed absolute energy uses the mean-square window correction.

## Field ESD and Fluence

Field ESD is

`S_E(f) = |E_tilde(f)|^2`

with one-sided weighting. Vector field power is

`|Ex_tilde|^2 + |Ey_tilde|^2 + |Ez_tilde|^2`.

For far-field radiative components,

`dF/df = |E_rad_perp_tilde(f)|^2 / Z0`, `Z0 = mu0 c`.

If the observer fails the far-field audit, spectral fluence is marked invalid
even though field ESD may still be output.

## Trusted Frequency Bounds

The trust report records:

- source validity and current provenance
- continuity and remap status
- `dt`, duration, `df`, Nyquist
- low-frequency cycle bound
- derivative trust frequency
- retarded-interpolation trust frequency
- mesh trust frequency
- sampling trust frequency
- far-field validity

The trusted upper frequency is the minimum of all active upper bounds. It is
not defined as a fixed fraction of Nyquist.

## Stage E Smoke Gate

Stage E source data remain:

`SCIENTIFIC_RF_VALID = false`

because the exported current is drift-only, the source time window is too short,
and the source snapshots are sparse. F3 can run FFT/ESD/CWT pipeline smoke on
these data, but VHF/UHF/GHz band results are not interpretable.

