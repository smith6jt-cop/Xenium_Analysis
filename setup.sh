#!/usr/bin/env bash
# setup.sh — provision Claude Code cloud env for the Xenium Shiny app
set -euo pipefail
export DEBIAN_FRONTEND=noninteractive

# -------------------------------------------------------------------
# 1. system deps
# -------------------------------------------------------------------
sudo apt-get update
sudo apt-get install -y --no-install-recommends \
  build-essential cmake gfortran pkg-config ca-certificates \
  curl wget gnupg lsb-release software-properties-common \
  git \
  libcurl4-openssl-dev libssl-dev libxml2-dev \
  libhdf5-dev hdf5-tools \
  libfontconfig1-dev libfreetype6-dev \
  libpng-dev libtiff5-dev libjpeg-dev \
  libcairo2-dev libxt-dev \
  libharfbuzz-dev libfribidi-dev \
  libgit2-dev \
  libglpk-dev libgmp3-dev libmpfr-dev libnlopt-dev \
  zlib1g-dev liblzma-dev libbz2-dev \
  pandoc \
  python3-venv

# -------------------------------------------------------------------
# 2. R 4.4+ (CRAN apt repo)
# -------------------------------------------------------------------
if ! command -v R >/dev/null 2>&1 || \
   awk -v ver="$(R --version 2>/dev/null | head -1 | awk '{print $3}')" \
       'BEGIN{exit !(ver < "4.4.0")}'; then
  wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
    | sudo tee /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc >/dev/null
  CODENAME=$(lsb_release -cs)
  echo "deb https://cloud.r-project.org/bin/linux/ubuntu ${CODENAME}-cran40/" \
    | sudo tee /etc/apt/sources.list.d/cran-r.list
  sudo apt-get update
  sudo apt-get install -y --no-install-recommends r-base r-base-dev
fi

# -------------------------------------------------------------------
# 3. Quarto (for the export tab's static report)
# -------------------------------------------------------------------
if ! command -v quarto >/dev/null 2>&1; then
  QV="1.6.40"
  wget -q "https://github.com/quarto-dev/quarto-cli/releases/download/v${QV}/quarto-${QV}-linux-amd64.deb" \
    -O /tmp/quarto.deb
  sudo dpkg -i /tmp/quarto.deb
  rm /tmp/quarto.deb
fi

# -------------------------------------------------------------------
# 4. Site-wide ~/.Rprofile — Posit Public Package Manager binaries.
#    PPM URL format is /cran/__linux__/<codename>/latest; HTTPUserAgent
#    must be set for binary delivery.
# -------------------------------------------------------------------
CODENAME=$(lsb_release -cs)
mkdir -p "$HOME"
cat > "$HOME/.Rprofile" <<EOF
local({
  options(
    repos = c(
      RSPM = "https://packagemanager.posit.co/cran/__linux__/${CODENAME}/latest",
      CRAN = "https://cloud.r-project.org",
      BioC = "https://bioconductor.org/packages/release/bioc"
    ),
    HTTPUserAgent = sprintf(
      "R/%s R (%s)",
      getRversion(),
      paste(getRversion(), R.version[["platform"]],
            R.version[["arch"]], R.version[["os"]])
    ),
    Ncpus = max(1L, parallel::detectCores() - 1L),
    download.file.method = "libcurl",
    timeout = 600
  )
  Sys.setenv(
    RENV_CONFIG_PAK_ENABLED = "TRUE",
    RENV_CONFIG_INSTALL_VERBOSE = "FALSE",
    LIBARROW_BINARY = "true",
    NOT_CRAN = "true"
  )
})
EOF

# -------------------------------------------------------------------
# 5. Persistent R library + renv cache
# -------------------------------------------------------------------
mkdir -p "$HOME/.R/library" "$HOME/.cache/R/renv"
echo 'R_LIBS_USER="$HOME/.R/library"' >> "$HOME/.Renviron"
echo "RENV_PATHS_CACHE=$HOME/.cache/R/renv" >> "$HOME/.Renviron"

# -------------------------------------------------------------------
# 6. Pre-install heavy packages (CRAN + Bioconductor + GitHub).
#    Use pak — it resolves system deps in one pass and uses PPM binaries.
# -------------------------------------------------------------------
R -q -e 'install.packages(c("pak","renv"))'

R -q -e 'pak::pak(c(
  "shiny","bslib","shinyWidgets","DT","plotly","ggplot2","patchwork",
  "Seurat","SeuratObject","clustree","harmony",
  "viridisLite","RColorBrewer","qs2","arrow","data.table",
  "fs","glue","future","future.apply","promises","later",
  "shinycssloaders","shinyjs","waiter","ggrastr",
  "shinyFiles","mclust","leidenAlg",
  "rmarkdown","knitr","quarto",
  "testthat","shinytest2","devtools","pkgload",
  "bioc::rhdf5","bioc::zellkonverter",
  "immunogenomics/presto"
))'

# -------------------------------------------------------------------
# 7. Sanity check
# -------------------------------------------------------------------
R -q -e 'cat("R:",as.character(getRversion()),"\n");
         pkgs <- c("Seurat","clustree","harmony","arrow","rhdf5","presto","qs2");
         stopifnot(all(sapply(pkgs, requireNamespace, quietly=TRUE)));
         cat("All required packages load OK\n")'
