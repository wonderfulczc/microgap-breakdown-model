# Stage 1.2 Analytical Lifecycle and Radiation Chain Method

## Scope and source

This stage reconstructs only Shi et al. 2019 Eqs. (4)-(6): the analytical lifecycle current moment, its derivative and spectrum, and analytical energy spectral density. It does not use Figshare data and does not reconstruct collision or isolated-streamer simulation curves.

## Lifecycle current moment

With x=t-T0 and q=alpha+beta,

I_CM(t) = I0 exp(alpha x) / [1 + exp(q x)].

The implementation evaluates `log(I0) + alpha*x - logaddexp(0,q*x)` and then exponentiates, avoiding overflow. Internal units are s, s^-1, and A m.

The analytical derivative is

dI_CM/dt = I_CM [alpha - q sigmoid(qx)].

Its time-domain unit is A m/s.

## Actual peak

Solving dI_CM/dt=0 gives

t_peak = T0 + ln(alpha/beta)/(alpha+beta),

I_peak = I0 [beta/(alpha+beta)] (alpha/beta)^[alpha/(alpha+beta)].

Therefore T0 is a transition parameter that controls the peak location, but for alpha != beta it is not the actual maximum time. The published equation is retained unchanged.

## Continuous Fourier transform

The convention is

I_tilde(omega) = integral I(t) exp(-i omega t) dt,

which yields

I_tilde(omega) = [pi I0/(alpha+beta)] exp(-i omega T0) / sin{pi(alpha-i omega)/(alpha+beta)}.

Magnitude evaluation uses |sin(a-ib)|^2 = sin^2(a)+sinh^2(b) in the log domain. Complex evaluation scales cosh/sinh by exp(-|b|), remaining finite through 100 GHz. At 100 GHz the true primary amplitude falls below binary64 representation and is explicitly reported as a legitimate underflow, not treated as a physical zero.

Independent validation uses SciPy `quad`. Low and moderate frequencies integrate the transformed dimensionless real-axis integral with cosine/sine weights. At high frequencies the contour is shifted below the real axis without crossing the nearest pole at -i pi; real and imaginary components are still independently evaluated by `quad`. This prevents exponentially small valid values from being lost to absolute quadrature noise.

## Shi 2019 Eq. (5) and dimensions

|F[dI_CM/dt]| = omega |I_tilde_CM(omega)|.

I_CM has A m, I_tilde has A m s, and multiplication by omega (rad/s, with rad dimensionless) gives A m. This equals the time Fourier transform unit of dI_CM/dt: (A m/s)*s = A m.

For frequency f in Hz, omega=2 pi f.

## FFT normalization and trust

The time window is expanded independently on both sides until I_CM(edge)/I_peak <= 1e-14, stricter than the required 1e-10. The sample interval provides at least 20 samples per target-highest-frequency period and a Nyquist frequency at least ten times that target. The grid is uniform and rounded up to a power of two; exceeding the configured point cap raises an error.

`dt * rfft(I_CM)` approximates the continuous transform. The time-origin phase factor is restored. The primary validation transforms the analytical derivative sampled on the grid; five-point and spectral derivative modes are also implemented. No mean removal and no window are applied. Zero padding is allowed only as interpolation in frequency and does not improve the physical resolution 1/(window duration).

A frequency is automatically marked untrusted if it is below five FFT bins, above one tenth of Nyquist, or if the analytical amplitude is below 1e-6 of the largest sampled validation amplitude. This prevents finite-window noise from being reported as a failed physical spectrum.

## ESD and frequency convention

The paper angular-frequency quantity is

P_omega(omega) = |omega I_tilde|^2 / (6 pi epsilon0 c^3),

and Eq. (6) uses W = 2 integral P_omega d omega. Since d omega = 2 pi df, the one-sided per-Hz ESD is

ESD_f(f) = 2 |2 pi f I_tilde|^2 / (3 epsilon0 c^3),

with J/Hz. Direct numerical integration verifies `integral ESD_f df = 2 integral P_omega d omega`.

VHF, UHF, and SHF use standard boundaries 30-300 MHz, 0.3-3 GHz, and 3-30 GHz. These energies apply only to the analytical lifecycle model and are not Shi 2019 Figure 4b collision energies.

## Data limitations

Figshare absence does not block Stage 1.2 Analytical Lifecycle and Radiation Chain because Eqs. (4)-(6) and their parameters are published. It does prevent reconstruction of Figure 3 inset, Figure 4a green collision FFT, Figure 4a purple isolated-streamer curve, Figure 2, and Figure 4b. No values were digitized or fabricated for these items.
