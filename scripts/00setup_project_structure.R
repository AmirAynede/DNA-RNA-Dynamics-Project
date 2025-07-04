# Create standard project directories using a vector and loop
project_dirs <- c(
  "data",
  "data/raw",
  "data/processed",
  "scripts",
  "scripts/analysis",
  "scripts/preprocessing",
  "results",
  "docs"
)
lapply(project_dirs, dir.create, showWarnings = FALSE)

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
