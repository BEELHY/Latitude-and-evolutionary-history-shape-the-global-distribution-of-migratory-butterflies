# Install HPC packages.

dir.create(Sys.getenv("R_LIBS_USER"), recursive=TRUE, showWarnings=FALSE)
.libPaths(Sys.getenv("R_LIBS_USER"))
install.packages(
  c("mgcv", "classInt", "spdep", "car", "ggplot2"),
  repos="https://cloud.r-project.org",
  lib=Sys.getenv("R_LIBS_USER")
)
