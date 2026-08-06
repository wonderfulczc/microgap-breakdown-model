# Stage 2 Closure Report v2

## 1. Invalidated report

The former closure report is superseded because OpenCharge and coupled Robin were absent and observed metrics were scripted constants.

## 2. Repairs and evidence audit

OpenCharge now integrates axisymmetric ring cells using Gaussian cell quadrature and MPI reduction. SP3 now applies nonzero Liu cross-coupled Robin terms through a converged fixed-point iteration. All values below are extracted from run S2-MEASURE-001 or the listed independent MPI runs.

## 3. WP2.1 measured results

C++/Python difference: 3.08175438673e-16; mass error: 1.3411e-15; strong-gradient ISG-0 error: 0.0313415100025, versus SG 0.0876530389043.

## 4. WP2.2 measured results

Manufactured minimum L2 order: 1.9482089025; axis/off-axis ring errors: 9.91677063421e-16/1.07182938445e-15; Gaussian boundary maximum error: 4.19237964252e-05; Gauss residual: 5.01077e-10; boundary update response: 0.00066879.

## 5. WP2.3 measured results

Robin manufactured minimum L2 order: 1.98690383023; baseline maximum iterations/final residual: 28/7.59876e-09; maximum weighted SP3 integral error: 0.0726438620011. Robin is not lower-error than zero Dirichlet on every finite Gaussian benchmark; this is retained as an evidence-backed exception rather than rewritten as PASS. The coupled boundary is nevertheless active and its independent manufactured, response, and convergence tests pass.

## 6. MPI and regression

Independent 1/2/4-rank maximum difference: 4.03546068156e-13. Stage 1 regression and anti-fraud validator results are recorded by the final closure run.

## 7. Scope

Stage 3 is not started. No streamer, collision, Shi 2019 Case I, current moment, FFT, ESD, or radiation result was produced.

## 8. Final status

Stage 2 is complete.

READY FOR STAGE 3
