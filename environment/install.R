# R 4.4.3. The panels were rendered with these package versions; missing packages are installed at that version,
# installed ones are left as they are (a version difference is only reported).
pkgs <- c(ggplot2 = "4.0.3", ragg = "1.5.2", systemfonts = "1.3.2", patchwork = "1.3.2", data.table = "1.18.6.1",
          ggrepel = "0.9.8", dplyr = "1.2.1", tidyr = "1.3.2", ggtext = "0.2.0", cowplot = "1.2.0")
missing <- names(pkgs)[!vapply(names(pkgs), requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "https://cloud.r-project.org")
  for (p in missing) remotes::install_version(p, version = pkgs[[p]], repos = "https://cloud.r-project.org", upgrade = "never")
}
for (p in names(pkgs)) {
  v <- as.character(packageVersion(p))
  cat(sprintf("%-12s %s%s\n", p, v, if (v == pkgs[[p]]) "" else sprintf("   (rendered with %s)", pkgs[[p]])))
}
