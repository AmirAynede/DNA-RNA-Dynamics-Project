# Create standard project directories
dir.create("data", showWarnings = FALSE)
dir.create("data/raw", showWarnings = FALSE)
dir.create("data/processed", showWarnings = FALSE)
dir.create("scripts", showWarnings = FALSE)
dir.create("scripts/analysis", showWarnings = FALSE)
dir.create("scripts/preprocessing", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
dir.create("docs", showWarnings = FALSE)

# Create a .gitignore file if it doesn't exist
if (!file.exists(".gitignore")) {
  writeLines(c(
    "# History files",
    ".Rhistory",
    ".RData",
    ".Rproj.user",
    "",
    "# Data folders",
    "data/raw/",
    "",
    "# Results",
    "results/",
    "",
    "# MacOS files",
    ".DS_Store"
  ), ".gitignore")
}

# Create a basic README.md if not already there
if (!file.exists("README.md")) {
  writeLines("# DRD_2025_Project\n\nThis project contains data analysis and documentation for the DRD 2025 research initiative.", "README.md")
}

# Create an empty R project file (if not already created from RStudio)
if (!file.exists("DRD_2025_Project.Rproj")) {
  fileConn <- file("DRD_2025_Project.Rproj")
  writeLines(c(
    "Version: 1.0",
    "",
    "RestoreWorkspace: Default",
    "SaveWorkspace: Default",
    "AlwaysSaveHistory: Default",
    "",
    "EnableCodeIndexing: Yes",
    "UseSpacesForTab: Yes",
    "NumSpacesForTab: 2",
    "Encoding: UTF-8",
    "",
    "RnwWeave: knitr",
    "LaTeX: pdfLaTeX"
  ), fileConn)
  close(fileConn)
}
