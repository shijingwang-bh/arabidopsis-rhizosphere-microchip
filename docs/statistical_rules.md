# Statistical rules for release 2026.11.25

This file follows Methods 2.8 and Supplementary Methods SM2/SM4 in the associated manuscript. **P** denotes an unadjusted probability and **q** denotes a Benjamini–Hochberg-adjusted P value. Unadjusted model tests use P <0.05; BH-corrected analyses use q <0.05. All reported effect confidence intervals are unadjusted 95% intervals.

| Analysis | Test and correction family | Current correspondence |
|---|---|---|
| Three chip–non-chip traits and pre-inoculation lateral roots | Two-sided Welch tests; all four comparisons jointly BH-corrected | S1d and 2a; one family despite figure relocation |
| Bacterial concentration within each day | Two-sided Welch tests; two day-specific comparisons jointly BH-corrected | 3c |
| Bacterial gravity × sampling-day model | Type II F tests; model-term P values unadjusted | Results 3.3; interaction F(1,13)=2.22, P≈0.160 |
| Day-7 plant-growth model terms | Type II F tests in gravity × inoculation models; primary and HC3 sensitivity P values unadjusted | S5a |
| Day-7 average inoculation effects | Linear-model contrasts with equal weighting across gravity; four traits jointly BH-corrected | 5d/e and S5b |
| RNA genes | Fractional summed Salmon counts; filterByExpr min.count=5; TMM, ordinary voom, robust eBayes; BH across 20,449 tested genes | 2b, S2, S1 table |
| GO BP gene sets | Fixed 15–500-gene memberships; 1,705 sets; BH separately within each method; primary CAMERA correlation=0.01; directional FRY and fgsea complements | 2d, S3, S2 table |
| Global chemical model | Seven conditional terms in the full gravity × inoculation × day model; 9,999 Freedman–Lane permutations, seed 20260912; BH across seven terms | 4b and S3a table |
| LC–MS feature contrasts | Two-sided Welch–Satterthwaite tests; BH across all 3,126 features separately within each contrast; q <0.05 and absolute effect >log2(1.5) | 4c/d, 5a–c, S6/S7, S3b/S4 tables |
| Root-zone fluorescence and conditioned-medium curves | Descriptive means/SE, individual trajectories; no inferential P/q tests | 3b/d/e and S4/S5 |

The growth family contains stem-to-tip root length, mature-root diameter, leaf fresh weight and leaf area, including the two traits not plotted in main Figure 5. Its full-precision P/q pairs are retained in `reference/plot_data/growth/Growth_traits_average_inoculation_effects.csv`. The q values are approximately 0.432, 0.00408, 0.00226 and 0.431, respectively. Selection of two displayed panels does not reduce the family to two tests.

The Type II inoculation main-effect test and the equal-weighted average inoculation contrast answer different hypotheses in the unbalanced design. For leaf fresh weight, the Type II P≈0.000176 and contrast P≈0.000566 are therefore retained as separate quantities; the latter is the input to the new four-trait BH correction.

LC–MS notation is G (gravity difference among uninoculated plants), I (inoculated-minus-uninoculated within gravity), J (I under simulated microgravity minus I under 1 g), and T=J7−J3. The four within-gravity/day inoculation comparisons identify a union of 579 responsive features. Day-7 J has 137 responsive features; day-3 J has none; T has five. Candidate examples retain q from their complete 3,126-feature contrast family. No correction is recalculated within the selected compounds or displayed heatmap rows.

Broader LC–MS families (3,414 named features and 14,031 QC-passing exported features), dispersion diagnostics and six neutral leave-one-out analyses are recomputed by `analyse_LCMS_global.R`. They are sensitivity analyses, not additional filters on the primary 3,126-feature results. No omission is chosen to change the primary sample set. GO estimated-correlation CAMERA is a separate sensitivity fit; no set remains q <0.05 under that fit.
