# DRD_2025_Project

**Comprehensive Analysis of DNA Methylation Using the Illumina HumanMethylation450 BeadChip Platform**

This repository contains the code and report for a DNA methylation analysis project completed as part of the final examination for the course **DNA/RNA Dynamics** at the University of Bologna.

---

##  Project Overview

This project explores genome-wide DNA methylation patterns using publicly available datasets from the **Illumina HumanMethylation450 BeadChip** platform. The goal was to identify differentially methylated regions (DMRs) and evaluate epigenetic modifications that may be relevant to biological processes or disease states.

Key steps in the analysis:

- Loading and preprocessing `.idat` files  
- Quality control and normalization  
- Identification of differentially methylated sites  
- Visualization of methylation profiles and clustering  
- Biological interpretation of significant results

---

## Technologies Used

- **Language**: R  
- **Key Packages**:
  - `minfi`
  - `limma`
  - `GEOquery`
  - `IlluminaHumanMethylation450kmanifest`
  - `IlluminaHumanMethylation450kanno.ilmn12.hg19`
  - `ggplot2`, `pheatmap`, `ComplexHeatmap`

---

##  Repository Contents

```bash
 DRD_2025_Project/
├── docs/		   # Project report and results (viewable in browser)
├── scripts/               # R scripts used for preprocessing, analysis, and visualization
├── data/                  # Sample input data
└── README.md              # This file
```

---

## Report

The full project report, including methods, figures, and results interpretation, can be found here:  
👉 [FinalReport.html](./docs/FinalReport.html)

---

## How to Run

To replicate or adapt the analysis:

1. Clone this repository:
   ```bash
   git clone https://github.com/kianinsilico/DRD_2025_Project.git
   ```

2. Open R or RStudio.

3. Install required packages:
   ```r
   install.packages(c("BiocManager", "ggplot2", "pheatmap"))
   BiocManager::install(c("minfi", "limma", "GEOquery",
                          "IlluminaHumanMethylation450kmanifest",
                          "IlluminaHumanMethylation450kanno.ilmn12.hg19"))
   ```

4. Run the scripts in the `scripts/` folder in order, or follow the steps in the report.

---

## Authors

This project was developed by:

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

For questions or feedback, feel free to open an [issue](https://github.com/kianinsilico/DRD_2025_Project/issues) or contact any of the contributors.
