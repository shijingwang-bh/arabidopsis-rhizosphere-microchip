# Arabidopsis rhizosphere microchip data and code

Data, analysis scripts and figure sources associated with **An Arabidopsis rhizosphere microchip reveals enhanced early accumulation of Bacillus subtilis under simulated microgravity**.

Public repository: [shijingwang-bh/arabidopsis-rhizosphere-microchip](https://github.com/shijingwang-bh/arabidopsis-rhizosphere-microchip).

Release **2026.11.25** corresponds to the latest manuscript and Supplementary Information identified by SHA256 in `metadata/manuscript_version.json`. It includes joint Benjamini–Hochberg (BH) correction across the four day-7 average inoculation effects on plant growth. The current manuscript contains Figures 1–5, Supplementary Figures S1–S7 (S5 has separate a/b files), and Supplementary Tables S1–S5. Figure S1d contains the cultivation assessment; RNA PCA appears only in Figure 2c.

## Contents

| Directory | Contents |
|---|---|
| `data/observations/` | Plant measurements, root-zone fluorescence, channel-liquid bacterial counts and archived loading-simulation values |
| `data/conditioned_medium/` | Original OD workbooks and cell-linked extracted measurements |
| `data/rna/` | Fractional Salmon count/TPM matrices, leaf-pool metadata, frozen GO memberships and annotation/selection references |
| `data/lcms/` | Provider intensity/annotation exports and corrected sample metadata |
| `code/` | Analysis, validation, plotting and figure-assembly scripts |
| `reference/plot_data/` | Complete archived plotting inputs and statistical families used by the current figures |
| `reference/statistics/` | Complete RNA gene and gene-set results |
| `reference/figures/` | Current submission figure references |
| `reference/tables/` | Current SI table contents as CSV |
| `assets/` | Preserved platform artwork and experimental image plates |
| `docs/` | Current legends, figure/data map, statistical rules, data dictionary and upload instructions |
| `environment/` | Verified and original package versions; optional dependency installer |
| `validation/` | Checks actually completed for this release |

The statistical source tables retain every tested feature/gene/trait rather than only significant rows. Semantic directory names are used for plotting data. Some original function and CSV names contain historical figure numbers; `metadata/plot_data_aliases.json` translates these implementation keys. They are not the current manuscript figure numbers. Use `docs/figure_data_map.csv` for manuscript correspondence.

## Reproduce

Use Python 3.10 or later. The verified R version is 4.5.1. Install Python dependencies explicitly:

```sh
python -m pip install -r requirements.txt
```

The R packages are listed in `environment/verified_package_versions.csv`. Statistical modules require `car`, `dplyr`, `jsonlite`, `vegan`, `edgeR`, `limma`, `statmod`, `fgsea` and `BiocParallel`. Figures additionally require `ggplot2`, `patchwork`, `magick`, `ggrepel` and `circlize`. Cairo, Arial and the system libraries required by the R `magick` package are used for plots. The optional command below installs missing packages from the repositories configured on your system; it does not pin all transitive versions:

```sh
Rscript --vanilla environment/install_R_packages.R
```

Verify the downloaded snapshot without R:

```sh
python validate_package.py
```

Run all analyses, compare their complete results with the manuscript references, and redraw the figures:

```sh
python reproduce.py --mode all --output ../microchip_run --rscript Rscript
```

The output directory must be new/empty and outside this package. It contains `results/`, `plot_data/`, `figures/`, step logs, `run_report.json` and `reproduction_report.json`. The runner verifies dependencies and checksums before execution. It never installs packages or changes the packaged data. Use `--rscript "path/to/Rscript.exe"` on Windows if Rscript is not on PATH; `--r-library "path/to/library"` can select an existing R library.

Separate entry points are available:

```sh
python reproduce.py --mode analysis --output ../analysis_run --rscript Rscript
python reproduce.py --mode figures --output ../figure_run --rscript Rscript
```

`figures` redraws from the archived reference tables; `analysis` recomputes statistics from the provided matrices/observations. `all` synchronizes recomputed tables into an isolated plotting tree before drawing. No original workspace, external cached RDS, installed R-library copy or manuscript file is required.

## Statistical and experimental scope

`docs/statistical_rules.md` records the test families. BH-adjusted P values are denoted **q**. The four average inoculation contrasts across stem-to-tip root length, mature-root diameter, leaf fresh weight and leaf area form one family. Their q values are approximately 0.432, 0.00408, 0.00226 and 0.431, respectively. Factorial model terms for bacterial concentration and plant growth retain unadjusted P values. Confidence intervals are unadjusted 95% intervals.

RNA uses nine leaf pools (five 1-g and four simulated-microgravity pools), each containing material from five chips. Fractional transcript-level estimated counts are summed to AGI genes, followed by filtering, TMM, ordinary voom and robust empirical Bayes. No effective-length offset is introduced. The family contains 20,449 genes; eight are significant. Gene-set correction is separate for each method across 1,705 frozen GO BP sets. Estimated-correlation CAMERA results are retained as a sensitivity analysis.

LC–MS samples are 41 independent culture-dish collections plus seven QC injections. N denotes uninoculated and P denotes inoculated in provider sample IDs; corrected metadata controls all analyses. The 3,126-feature family uses QC CV ≤30%. BH correction covers all features within each contrast, with an additional absolute-effect threshold of log2(1.5). Global seven-term tests use 9,999 Freedman–Lane permutations. Broader feature-family and leave-one-out sensitivity analyses are retained in the analysis code.

Conditioned-medium OD wells are technical replicates of one pooled preparation per gravity-by-day condition. The uninoculated day-7 E2 well is retained in the individual-well display and excluded from inoculated-well summaries. No smoothing, time-point deletion or clipping of negative corrected OD values is applied.

## Reproduction boundaries

RNA reproduction starts from the supplied Salmon count/TPM matrices. FASTQ reads, `quant.sf`, effective lengths and the original upstream command are not bundled, so the package does not reconstruct sequencing or upstream quantification. The manuscript identifies BioProject PRJNA1541220. LC–MS starts from provider intensity/annotation exports, not instrument acquisition files or de novo peak calling. The `.xls` LC–MS files are GB18030 tab-separated text, not binary Excel workbooks; preserve their bytes.

Image-based platform artwork and microscopy plates are archived inputs. Their original acquisition, chip-fabrication and COMSOL model workflows are not reconstructed by the statistical scripts. Loading-simulation values are included as a table, without claiming to regenerate the simulation. Figure S7 has a new vector redraw from the same observations, confidence intervals and full-family q values because its original vector drawing script was unavailable. Its exact submitted raster appearance is preserved separately in `assets/archived_charts/` and `reference/figures/Figure_S7.pdf`. This package does not modify the manuscript or SI figures.

## Validation, citation and reuse

See `validation/VALIDATION.md` for the executed checks and their limits. `checksums_sha256.txt` covers the released files. `metadata/source_provenance.csv` records original project-relative data/asset sources without local user paths; `metadata/code_provenance.csv` records adapted analysis kernels and plotting sources. Cite the associated manuscript and this release; `CITATION.cff` contains the author metadata and repository URL. An archive DOI has not been assigned.

The existing **MIT code license** is retained in `LICENSE`. Data retain the existing separate terms in `DATA_REUSE_TERMS.txt`; the code license does not assign a new data license. See `docs/UPLOAD.md` for uploading the unpacked repository files.
