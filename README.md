<p align="center">
  <img src="https://img.shields.io/badge/DRD_2025_Project-DNA%20Methylation%20Pipeline-blueviolet?style=for-the-badge&logo=dna&logoColor=white" alt="DRD_2025_Project Logo" />
</p>

# DRD_2025_Project

[![Documentation](https://img.shields.io/badge/docs-online-blue.svg)](docs/FinalReport.html)
[![PDF Report](https://img.shields.io/badge/report-pdf-green.svg)](docs/FinalReport.pdf)
[![Results](https://img.shields.io/badge/results-up--to--date-brightgreen.svg)](results/)
[![Reproducible](https://img.shields.io/badge/reproducibility-yes-success.svg)](README.md)
[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/license-CC--BY--NC--SA%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-sa/4.0/)
[![Last Commit](https://img.shields.io/github/last-commit/kianinsilico/DRD_2025_Project.svg)](https://github.com/kianinsilico/DRD_2025_Project/commits/main)
[![R >= 4.0](https://img.shields.io/badge/R-%3E=4.0-blue.svg)](https://cran.r-project.org/)
[![Issues](https://img.shields.io/github/issues/kianinsilico/DRD_2025_Project.svg)](https://github.com/kianinsilico/DRD_2025_Project/issues)
[![Issues Closed](https://img.shields.io/github/issues-closed/kianinsilico/DRD_2025_Project.svg)](https://github.com/kianinsilico/DRD_2025_Project/issues?q=is%3Aissue+is%3Aclosed)

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
├── docs/                # Project report, documentation, and [README](docs/README.md)
├── scripts/             # R scripts for setup and analysis
├── results/             # Output results, figures, and [README](results/README.md)
├── lib/                 # Local R package tarballs (e.g., SummarizedExperiment)
├── .gitignore           # Git ignore rules
├── LICENSE              # License file
├── CITATION.cff         # Citation file
└── README.md            # This file
```

- See [`docs/README.md`](docs/README.md) for details on documentation and reports.
- See [`results/README.md`](results/README.md) for details on results and figures.

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
   - Follow the steps in [`docs/FinalReport.Rmd`](docs/FinalReport.Rmd) or run the scripts in `scripts/` as needed.
   - Output and figures will be saved in [`results/`](results/) and [`docs/`](docs/).

---

## Report & Results

- 📄 **HTML Report:** [docs/FinalReport.html](docs/FinalReport.html)
- 📄 **PDF Report:** [docs/FinalReport.pdf](docs/FinalReport.pdf)
- 📊 **Results & Figures:** [results/](results/)

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
