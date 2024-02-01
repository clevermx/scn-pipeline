#!/usr/bin/env bash

set -euo pipefail

## lines below require GITHUB_PAT in ~/.Renviron


R -e "devtools::install_github('satijalab/seurat-wrappers', ref='master')"
