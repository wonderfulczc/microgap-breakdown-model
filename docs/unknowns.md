# Unknowns

| unknown_id | description | affected_module | impact_level | why_unresolved | current_handling | resolution_method | blocking_stage | status |
|---|---|---|---|---|---|---|---|---|
| UNK-0001 | Figshare original data unavailable | streamer_collision | critical | link only; files absent | independent reconstruction | request archive; observable comparison | Stage 1 closure:no; transfer_to:Stage 5; Pointwise reproduction | open |
| UNK-0002 | original grid spacing unknown | mesh | critical | not reported | no value | grid convergence | Stage 1 closure:no; transfer_to:Stage 4 | open |
| UNK-0003 | output sampling interval unknown | fft | critical | not reported | no value | sampling sensitivity | Stage 1 closure:no; transfer_to:Stage 4; Pointwise reproduction | open |
| UNK-0004 | numerical differentiation unknown | fft | high | not reported | unfrozen | operator sensitivity | Stage 1 closure:no; transfer_to:Stage 4; Pointwise reproduction | open |
| UNK-0005 | FFT mean removal unknown | fft | high | not stated | unfrozen | with/without mean | Stage 1 closure:no; transfer_to:Pointwise reproduction | open |
| UNK-0006 | FFT window unknown | fft | high | not stated | unfrozen | window sensitivity | Stage 1 closure:no; transfer_to:Pointwise reproduction | open |
| UNK-0007 | Poisson tolerance unknown | poisson | critical | not stated | no value | residual sweep | Stage 1 closure:no; transfer_to:Stage 2; Pointwise reproduction | open |
| UNK-0008 | SP3 tolerance unknown | photoionization | critical | not stated | no value | residual sweep | Stage 1 closure:no; transfer_to:Stage 2; Pointwise reproduction | open |
| UNK-0009 | Case II and III domains unknown | domain | high | not reported | no bounds | boundary sensitivity | Stage 1 closure:no; transfer_to:Stage 5 | open |
| UNK-0010 | electron transport table unavailable | transport | critical | not published | blocked | obtain table or later independent reconstruction | Stage 1 closure:no; transfer_to:Stage 2 | open |
| UNK-0011 | diffusion current inclusion ambiguous | current_moment | critical | definition insufficient | do not freeze | drift vs drift+diffusion | Stage 1 closure:no; transfer_to:Stage 4 | open |
| UNK-0012 | ion current inclusion ambiguous | current_moment | high | not unique | do not freeze | electron vs all species | Stage 1 closure:no; transfer_to:Stage 4 | open |
| UNK-0013 | displacement current inclusion ambiguous | radiation | critical | not settled | do not add | conduction vs total current | Stage 1 closure:no; transfer_to:Stage 4 | open |
| UNK-0014 | Figure 4b integration bins unknown | esd | high | not specified | unfrozen | endpoint/bin sensitivity | Stage 1 closure:no; transfer_to:Stage 5; Pointwise reproduction | open |
| UNK-0015 | isolated-streamer control construction unknown | collision | critical | insufficient detail | unfrozen | compare control constructions | Stage 1 closure:no; transfer_to:Stage 4 | open |
| UNK-0016 | Figure 1 domain needs sensitivity validation | domain | high | graphical inference | mark inferred | expand boundaries | Stage 1 closure:no; transfer_to:Stage 4 | open |
| UNK-0017 | SP3 open boundary may differ from original | photoionization | critical | code absent | literature-level boundary definition frozen; implementation validation deferred | Liu benchmark and boundary sweeps | Stage 1 closure:no; transfer_to:Stage 2; Pointwise reproduction | open |
| UNK-0018 | I_CM0 conflict 0.44 versus 0.4 A m | lifecycle | high | paper internal conflict | primary = 0.44 A m; Figure 3 caption sensitivity variant = 0.40 A m | researcher decision | Stage 1 closure:no; transfer_to:Stage 4 sensitivity branch | open |
| UNK-0019 | ISG exponential interpolation dimensional reference density has not been frozen. | numerical_flux | critical | Kulikovsky Eq. (20) uses n+1 in normalized variables; a dimensional implementation requires an explicit reference density convention. | no numerical value assigned | freeze n_floor convention and validate nondimensional invariance before WP2.1 implementation | Stage 1 closure:no; transfer_to:Stage 2 | deferred_to_stage2 |
