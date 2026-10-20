# Literature Contract, Reviewed 2026-10-16

The separately named Chinese attachment was not located. The full user message
is the authorized decision/protocol; no missing file is claimed read.

|ID|Original publication|Actually inspected|Bounded evidence and limits|
|---|---|---|---|
|L1|Taka & Fujiwara, 2010, Verification of Spark-Resistance Formulae for Micro-Gap ESD, E93-B(7),1801-1806, [DOI](https://doi.org/10.1587/transcom.E93.B.1801)|IEICE publisher metadata and abstract only|Human ESD and waveform-derived spark comparisons; abstract favors RW for its tested conditions. No full experimental arrays or formula constants obtained.|
|L2|Taka & Fujiwara, 2018, Validity and Applicability of Spark-resistance Formulae for Micro-gap Sparks in Parallel Disc Electrodes,138(6),316-321, [DOI](https://doi.org/10.1541/ieejfms.138.316)|J-STAGE metadata/abstract; PDF open failed|Parallel-disc comparisons distinguish500V and4000V; waveform agreement and inferred gap agreement are different criteria. No universal law selection follows.|
|L3|Fujiwara,2025, EMC Research Efforts in the IEEJ Technical Committee and State-of-the-Art Case Studies on ESD Found in its Technical Papers,20,4-18, [DOI](https://doi.org/10.1002/tee.24152)|Publisher HTML sections3.3-3.4, equations2-3 and unit definitions|Toepler/RW history formulas and constant units are legible. Constants vary by convention, pressure and inferred gap; no values transferred to this device. Review article, not the original Toepler/RW derivation.|
|L4|Schrama, Bavor, Flammer & Durfee,2026, Energy partitioning in electrostatic discharge with variable series load resistor,140,104253, [DOI](https://doi.org/10.1016/j.elstat.2026.104253)|Publisher search-delivered abstract and HTML sections1/3/5; subsequent direct open failed|Series-load energy allocation supports coupled circuit/channel review. Effective electron energy is empirical, not a measured gas-heating or RF fraction for this TENG. Equations not fully rendered, no implementation borrowed.|
|L5|Fujita, Zhou & Pommerenke,2017, Hybrid full-wave/circuit modelling of spark gaps and its experimental validation,53(7),484-486, [DOI](https://doi.org/10.1049/el.2016.4338)|Publisher HTML abstract/introduction/model text; equation1 is an image not transcribed|RW/SPICE and FDTD can be coupled. Its differently parameterized coefficient has units m^2/(V^2*s), not the K_RW below; no assumed numerical conversion or device calibration.|
|L6|Popov & Starikovskaia,2022, Relaxation of electronic excitation in nitrogen/oxygen and fuel/air mixtures: fast gas heating in plasma-assisted ignition and flame stabilization,91,100928, [DOI](https://doi.org/10.1016/j.pecs.2021.100928)|Publisher abstract/highlights only|Electronic relaxation and mixture/reduced-field regimes matter for fast heating. No device-specific partition factor obtained.|
|A1|Marszalek, Memristive Fingerprints of Electric Arcs, [arXiv1601.01612](https://arxiv.org/html/1601.01612)|Author manuscript HTML, section1 equations1-3|Hybrid Cassie-Mayr, u=i/g, external circuit explicit. Conventional pure-model forms below are specializations with offset/radiation terms removed, not a transcription of the complete hybrid equation or proof of micro-gap applicability. Original1939/1943 texts not obtained.|

Verification is per accessible section, not a claim to have reviewed every original
full article. L1's abstract describes RW as preferable across its settings, while
L2's retrospective distinguishes earlier200/2000V conclusions. This discrepancy
is retained: observable, assumed constants and fit criterion must be checked in
full originals before a physical choice. It does not block synthetic resistor cases.

## Equations, Dimensions and Domain

Current is the signed **column conduction current**, not port current including Cp.

```text
qstar = q0 + integral Icol dt                         [C]
R_T = K_T*ell/qstar                                  [ohm]
K_T: ohm*C/m = V*s/m; ell:m; q0:C
Sstar = S0 + integral Icol^2 dt                       [A^2*s]
R_RW = K_RW*ell/sqrt(Sstar)                          [ohm]
K_RW: ohm*A*sqrt(s)/m = V*sqrt(s)/m
Mayr: (1/g)*dg/dt = (u*Icol/P_M - 1)/theta_M          [1/s]
Cassie: (1/g)*dg/dt = (u^2/U_C^2 - 1)/theta_C         [1/s]
u=Icol/g; g:S; P_M:W; U_C:V; theta:s
```

Toepler requires sourced positive q0 and qstar>0; signed charge reversal may
exit its domain. Return OUT_OF_MODEL_DOMAIN at zero/negative accumulated charge;
never replace signed I with abs(I), floor q or create a recovered plasma state.
RW needs sourced S0 and constants; its nondecreasing squared-current history
cannot model resistance recovery by itself. Neither law supplies complete port
impedance or energy partition. Mayr/Cassie need a compatible arc stage, initial
conductance and loss/time constants; they are deferred, not executed or rejected.
No classic spark law is integrated in this task.

## Decision and Energy Ownership

Stop extending the10V/artificial10000K/near2atm LTE cylinder as a real70um R(t).
Preserve it as engineering evidence and a sensitivity counterexample. A collapsed
post-breakdown voltage near10V is not intrinsically impossible: the unsupported
part was replacing pre-event stored charge/topology/history by a sustained ideal
source. Do not multiply its deposited work by(800/10)^2.

Circuit/field accounting and material accounting are separate. Conduction work
does not specify electron/heavy/chemical/photon/electrode partitions. Escape of
optical energy and subsequent photon absorption/photoionization must not be
counted twice. RF boundary power and loaded receiver energy need independent
spatial propagation/receiver maps. Missing components remain UNAVAILABLE, not
residual-derived heat or radiation. C1 internal GR is not a resistor.
