# DRD_2025_Project

**Comprehensive Analysis of DNA Methylation Using the Illumina HumanMethylation450 BeadChip Platform**

This repository contains the code, data, and report for a DNA methylation analysis project completed as part of the final examination for the course **DNA/RNA Dynamics** at the University of Bologna (2025).

---

## Project Overview

This project explores genome-wide DNA methylation patterns using the **Illumina HumanMethylation450 BeadChip** platform. The main goal is to identify differentially methylated regions (DMRs) and evaluate epigenetic modifications relevant to biological processes or disease states, using real data and a reproducible R-based pipeline.

The analysis includes:
- Data import and preprocessing of raw `.idat` files
- Quality control and normalization (e.g., preprocessNoob)
- Probe annotation and filtering
- Calculation of beta and M values
- Principal Component Analysis (PCA)
- Differential methylation analysis (t-test, multiple testing correction)
- Visualization and interpretation of results

---

## Workflow Summary

The main steps of the analysis pipeline are:

1. **Project Setup**: Run `scripts/00setup_project_structure.R` to create the required folder structure and project files.
2. **Data Import**: Load raw data and sample sheet using `minfi`.
3. **Quality Control**: Assess data quality using signal intensities, detection p-values, and control probes.
4. **Normalization**: Apply normalization (e.g., preprocessNoob) to correct technical variation.
5. **Exploratory Analysis**: Calculate and visualize beta/M values, perform PCA, and check for batch effects.
6. **Differential Methylation**: Identify differentially methylated probes between groups using statistical tests and multiple testing correction.
7. **Reporting**: All code and results are documented in `docs/FinalReport.Rmd` and rendered as HTML/PDF.

---

## Repository Structure

```bash
DRD_2025_Project/
├── data/                # Raw and processed data files
│   ├── raw/             # Raw .idat files and sample sheet
│   └── processed/       # RData files generated during analysis
├── docs/                # Project report (Rmd, HTML, PDF)
├── scripts/             # R scripts for setup and analysis
├── results/             # Output results and figures
├── lib/                 # Local R package tarballs (e.g., SummarizedExperiment)
├── .gitignore           # Git ignore rules
├── LICENSE              # License file
├── CITATION.cff         # Citation file
└── README.md            # This file
```

---

## How to Run the Analysis

1. **Clone this repository:**
   ```bash
   git clone https://github.com/kianinsilico/DRD_2025_Project.git
   cd DRD_2025_Project
   ```

2. **Set up the project structure (optional, if not already present):**
   ```r
   source("scripts/00setup_project_structure.R")
   ```

3. **Install required R packages:**
   ```r
   install.packages(c("BiocManager", "gplots", "qqman", "tinytex"))
   BiocManager::install(c(
     "minfi", "minfiData", "sva", "shinyMethyl",
     "IlluminaHumanMethylation450kmanifest",
     "IlluminaHumanMethylation450kanno.ilmn12.hg19",
     "AnnotationDbi"
   ))
   # Install SummarizedExperiment from local tarball
   install.packages("./lib/SummarizedExperiment_1.38.0.tar.gz", repos = NULL)
   ```

4. **Run the analysis:**
   - Follow the steps in `docs/FinalReport.Rmd` or run the scripts in `scripts/` as needed.
   - Output and figures will be saved in `results/` and `docs/`.

---

## Report

The full project report, including methods, code, figures, and interpretation, can be found here:
👉 [FinalReport.html](./docs/FinalReport.html)

---

## Authors

- **Keshani, K.**
- **Castro Vargas, P.**
- **Gustini, B.**
- **Aynede, A.**
- **Dobrokhotov, I.**
- **Otoijamun, F.**

As part of the course *DNA/RNA Dynamics*, University of Bologna (2025).

---

## License

[LICENSE](./LICENSE)

---

## Contact

For questions or feedback, please open an [issue](https://github.com/kianinsilico/DRD_2025_Project/issues) or contact any of the contributors.
