# MuSiC component intervention design

Status: completed on 2026-09-25; results and technical checks are in results/music_components/RESULTS.md. Designed after inspecting the PBMC and pancreas results, then frozen with source/input hashes before intervention predictions. The original protocol bytes are retained as results/music_components/protocol.frozen.md. This is a post-hoc computational explanation experiment, not a preregistered biological mechanism test. The design below records the fixed rules; it is not an instruction to rerun completed fits. The earlier diagnostic protocol and results remain unchanged.

## Question and source

Which changes in mean donor expression profiles, cell-size estimates and between-donor variance contribute to the allocation-by-budget contrast within the tested MuSiC calculation?

Use unchanged MuSiC 1.0.0, commit f21fe67f5670d5e9fca0ad7550abaae3423eb59c, with the existing project-local R runtime. The installed music.iter function was inspected. It accepts Y, D, S and Sigma, aligns genes, scales Y by 100 when normalize=FALSE, and calls music.basic. Existing reconstruction in run_music_mechanism_v1.R establishes D[,k] = M.theta[,k] * M.S[k]. Therefore D and S must not be treated as independent primitive factors. Do not swap final weights: weights depend on evolving residuals and coefficients as well as Sigma.

## Fixed inputs

Reuse the existing mechanism selection: all 14 PBMC donors, 14 reference triples, original three reference blocks, six types, and the same eight identifier-selected targets per donor. Use budgets 60 and 300, balanced and 10:1:1 allocations, and every dominant-donor choice. The selection was already inspected for diagnostics, so it is not an untouched validation set. Do not select donors or targets by effect size.

For each triple/block, define one common gene set from the intersection of valid official supports across both budgets and all four allocations. Preserve the original full-inventory profile denominators and library-size calculation before restricting rows. Preserve official target scaling and any target-specific filtering identically across interventions. Do not renormalize profiles over the common set. Report the original-support result and the common-support result separately: gene filtering is itself an intervention.

## Eight controlled combinations

At each budget and paired balanced/unequal reference, extract three primitive components: mean relative-expression profile Theta; mean cell-size vector S; between-donor variance Sigma. Give each component either its balanced value (0) or unequal value (1), producing all eight 000–111 combinations. Rebuild D = Theta * S by columns at every combination, and pass the same selected S to the original iterative fitter. Keep Y, cell types, gene order, iteration parameters and starting procedure fixed. Residuals and weights are allowed to update normally; they are downstream calculations, not additional independently swapped factors.

000 and 111 reproduce the respective common-support balanced and unequal fits. The six hybrids are artificial computational interventions and must never be described as actual single-cell references or a new validated method. Parameter defaults remain iter.max=1000, nu=0.0001, eps=0.01, centered=FALSE and normalize=FALSE.

## Endpoints and aggregation

For every target let f(A) be MAE in percentage points when the components in subset A use unequal values and the rest use balanced values. For each of the three components j compute its Shapley contribution:

phi_j = sum over A not containing j of [|A|! (2-|A|)! / 3!] * [f(A union {j}) - f(A)].

This averages the six possible replacement orders and distributes nonlinear interactions by an explicit convention. Check sum(phi_j) = f(111)-f(000) numerically. It is a decomposition of the specified computational intervention, not a biological mediation proportion.

The single primary endpoint is mean(phi_Sigma at 300 - phi_Sigma at 60). The total common-support allocation-by-budget contrast and corresponding Theta and S contributions are required context. Secondary outputs: all eight absolute MAEs, per-type errors, donor and dominant-choice distributions, component interactions, and the difference between common-support and original-support contrasts. Retain contributions of either sign. Do not divide by a near-zero total contrast or present a percentage of mechanism explained.

Aggregate paired target differences, then dominant choices, triples and held-out donors using the existing equal-donor convention. Report individual block contrasts and SD over the three complete blocks. Three draws do not justify a precise Monte Carlo or population uncertainty claim. Shared donors, cells and targets prohibit treating individual predictions as independent biological replicates.

## Technical controls and bounded execution

Before dispatch, use the lexicographically first triple, block 0, both budgets, balanced plus the first dominant donor, and all 32 existing selected targets eligible for that triple. Verify the unchanged full-support adapter against corresponding saved predictions with maximum absolute proportion difference <=1e-10, without rerunning unrelated historical cases. Verify the basis reconstruction and 000/111 common-support endpoint identities at <=1e-10. Check explicit gene/type identities, finite nonnegative Sigma, unit-sum finite predictions, convergence status and the Shapley sum identity at <=1e-10 percentage points. Save failures and stop; do not relax tolerances to obtain a scientific answer.

The complete design has at most 14*3*2*3*32*8 = 64,512 logical weighted fits before deduplicating identical endpoints. Record actual unique fit count. Time the control batch and project the remaining cost; use at most three workers, one hour of fitting and 1 GiB of new outputs. Reuse completed controls. A limit or implementation failure yields an explicitly incomplete experiment, not selective omission of difficult cases. No new cell sampling is required.

## Interpretation and continuation

A negative primary contribution would support a role for the Sigma pathway in the negative budget contrast within these interventions. A zero or positive contribution would weaken that explanation even if profile dispersion falls. A dominant Theta or S contribution changes the explanation rather than invalidating the experiment. Strong dependence on the common gene set limits all attribution claims.

Do not launch a donor-weight correction from this result alone: changing donor weights also changes the biological mean being estimated. If a correction is later proposed, choose its tuning without held-out target truth and evaluate it on a separately fixed dataset. A mechanism null result is retained; do not replace the primary endpoint or start additional component searches to rescue the story.

Official interface: https://xuranw.github.io/MuSiC/reference/music.iter.html
Fixed source: https://github.com/xuranw/MuSiC/tree/f21fe67f5670d5e9fca0ad7550abaae3423eb59c

Implementation clarification: full-inventory denominators above means full-gene library totals within each selected reference, not profiles computed from all eligible source cells. The frozen design and original wording remain in the recorded protocol hash; no input selection or endpoint changed.
