# Stage 1 Method Summary

## 1. Objective

Stage 1 determines the target paper's model structure, identifies published parameters and equations, records missing information, and verifies the analytical chain from current moment to radiation spectrum. It does not test whether a streamer solver can reproduce the collision.

## 2. Literature chain

The local chain contains nine papers plus one unavailable data source: Shi 2019 is the reproduction target; Bourdon 2007 defines SP3 and Helmholtz photoionization parameters; Liu 2007 supplies an SP3 benchmark; Kulikovsky 1995 defines ISG-0; Ihaddadene and Celestin 2015 supplies collision context; Shi 2016 and Shi 2017 establish the model lineage; Luque 2017 supplies full-Maxwell cross-validation context; Liu and Pasko 2004 supplies earlier fluid/photoionization lineage; and Shi 2019 Figshare data is registered as unavailable.

## 3. Frozen physical model

The intended later reconstruction uses a two-dimensional axisymmetric, three-species drift-diffusion fluid model, quasistatic Poisson field, three-group SP3 photoionization, ISG-0 electron transport, current-moment integration, far-field radiation, the lifecycle model, and Fourier/ESD analysis. The fluid, Poisson, SP3, and ISG-0 portions are frozen only at literature-definition level and are not implemented.

## 4. Implemented content

Stage 1 includes parameter/equation audit, a stable Shi 2019 Eq. (4) implementation, its analytical derivative and continuous Fourier transform, sampled-FFT cross-validation, ESD and band energy, the Figure 3 main curves, and three analytical Figure 4a curves. The primary amplitude is 0.44 A·m; 0.40 A·m is a caption-conflict sensitivity variant. `T0` is not the actual peak time: the evaluated peak is 7.758250939 ns.

## 5. Not implemented

Electron transport, Poisson, SP3, streamer propagation, collision, spatial current-moment integration, Figure 2, the Figure 3 inset, Figure 4a's green collision-FFT and purple isolated-streamer curves, and Figure 4b are not implemented or reconstructed.

## 6. Conclusion

Stage 1 answers what must be reproduced, which equations and parameters apply, what is known or unknown, and whether the analytical radiation chain is correct. It does not answer whether a streamer solver can produce the paper's collision process.
