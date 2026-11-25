# Data dictionary and experimental units

## Group labels

`stationary_1g` is the stationary 1-g condition and `clinorotation` is the simulated-microgravity condition. `uninoculated`/`inoculated` are U/I. These labels describe experimental conditions; no code converts clinorotation into a spaceflight measurement.

In `phenotype_pre.csv`, the historical `gravity` field has two meanings. For **Lateral root number**, it denotes gravity. For root length, leaf area and leaf fresh weight in the separate cultivation assessment, `stationary_1g` denotes the non-chip control and `clinorotation` denotes chip-grown plants; both were physically under stationary 1 g. This original storage convention is retained to reproduce the four Welch tests. The plots and captions use Control/Chip labels for those three traits.

In `growth_day7.csv`, the historical metric string **Root elongation (mm)** corresponds to the day-7 stem-base-to-longest-root-tip length described in the manuscript. It is not interpreted as a paired increment from an earlier root-length measurement. The label is preserved in archived tables to maintain exact row correspondence.

## Observation tables

| File | Main fields and units | Independent unit |
|---|---|---|
| `data/observations/phenotype_pre.csv` | `metric`, historical `gravity`, replicate fields, `value`; count, mm, cm² or mg according to metric | Plant in each separate comparison |
| `data/observations/growth_day7.csv` | `metric`, `gravity`, `bacteria`, replicate fields, `value`; mm, µm, mg or cm² | Individual measured plant; trait-specific sample sizes |
| `data/observations/rfp.csv` | `gravity`, `region`, `day`, chip/replicate fields, `value` in arbitrary fluorescence units | Same three tracked chips per gravity, repeatedly imaged |
| `data/observations/bacterial_counts.csv` | `gravity`, `day`, replicate fields, `cells_1e8_per_ml` | Chip; concentration is in units of 10^8 cells/mL |
| `data/observations/loading_simulation.csv` | `Inlet_rate_uL_per_min`, `Simulated_outlet_arrival_time_s` | Archived deterministic simulation setting |

Historical Figure1d/Figure3 CSV names in `reference/plot_data/` identify original source tables. Current manuscript figure numbering is listed separately in `figure_data_map.csv`.

## Conditioned-medium workbooks and trajectories

`data/conditioned_medium/day3.xlsx` and `day7.xlsx` retain the original workbooks. The extraction script reads cached cell values and preserves cell addresses. `growth_measurements_raw.csv` includes:

- `collection_day`: 3 or 7; `well`: original well ID; `role`: assay/MS_blank/LB_blank.
- `group`: medium derived from 1-g or clinorotated uninoculated chips, half-strength MS control or LB control.
- `included_main`: whether the assay well contributes to the inoculated-well mean. Day-7 E2 is uninoculated and excluded from that mean.
- `time_index`: 1–289; `time_h`: 0–72 in 0.25-h steps; `od_raw`: uncorrected OD600.
- `source_sheet`, `source_cell`: workbook location. The separate cached-calculation CSV allows comparison with the workbook's background means and corrected values.

Correction subtracts the appropriate time-specific MS or LB blank mean. The plotted bands are SE across technical wells; blank-control panels use SD. Each gravity × collection-day conditioned-medium group is one pooled preparation. Technical well count is not the number of independent biological collections.

## RNA matrices and pools

`cdna.count` and `cdna.tpm` are tab-separated tables. `Name` contains transcript IDs with terminal isoform suffixes; A1–A5 are five 1-g leaf pools and B1–B4 are four simulated-microgravity leaf pools. Each pool comprises leaf material from five chips. Transcript fractional estimated counts and TPM are summed separately after removing the terminal isoform suffix. Raw matrix values are not rounded.

`RNA_sample_metadata_pool_units.csv` records pool and tissue definitions. The author-supplied QC table and its definitions are retained. `GO_membership.csv` maps `gene_id` to `GO_ID`; `GO_dictionary.csv` supplies term names. These frozen files avoid changes caused by later annotation database releases. Gene reference tables supply existing gene descriptions and the prespecified displayed rows; they do not replace recomputed gene tests.

Gene output uses `logFC` (sµg minus 1 g on log2 scale), `P.Value` (raw), `adj.P.Val` (BH q), and other limma statistics. GO output includes primary CAMERA `FDR`, directional `FRY_q`, `estimated_rho_q` and `fgsea_NES`. The plot x-axis may show fgsea NES while the adjacent q columns refer to CAMERA/FRY; these are distinct method outputs.

## LC–MS exports

The three `.xls` exports are GB18030-encoded TSV text. They contain sample intensities and provider annotation fields; they are not binary Excel workbooks. `All_metabolites_ori.xls` contains 15,593 exported features, `metabolites_exp.xls` contains 3,414 named features, and the primary QC-passing family contains 3,126 features. The primary and named exports differ only at floating-point text precision for shared intensities; the analysis preserves the named-export values as in the manuscript workflow.

`corrected_sample_metadata_with_QC.csv` defines 41 independent culture-dish collections plus seven QC injections. N in sample IDs denotes uninoculated and P denotes inoculated. The metadata, rather than a guessed letter interpretation, determines gravity, inoculation, day and group. These samples are independent of the chip-imaging and pre-inoculation RNA pools. No RNA/LC–MS sample pairing is asserted.

Feature IDs such as `pos_5231` and `neg_3047` preserve provider polarity/feature identity. `qc_cv` is SD/mean across QC injections. `log2_effect`, `se`, `df`, `p`, `q`, `ci95_low`, `ci95_high` and `pass` describe each complete feature contrast. Annotation fields, reported mass error and low-intensity flags are retained as supplied/derived; names are annotations rather than new structural identifications.

## Images and references

Platform artwork, microscopy images and longitudinal plates are preserved under `assets/`. Image preprocessing and scale bars are retained from the archived study graphics. Full acquisition stacks and the original COMSOL model are outside this statistical package. Figure S7's submitted raster and its vector numeric redraw are explicitly distinguished in the README. Current SI CSV tables retain manuscript rounding; full precision is available in the statistical source tables.
