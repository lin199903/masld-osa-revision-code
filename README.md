# MASLD-OSA Revision Code Package (Code-only)

This repository contains code and configuration for revision analyses.

## Included
- `scripts/*.R`
- `run_revision_pack.R`
- `config.yml`
- `outputs/sessionInfo/environment.yml` (if available)
- `manual_data/negative_controls/README_download_and_place.md`

## Quick start
1. Open R (>=4.3 recommended).
2. Install required packages listed in scripts and/or environment file.
3. Set working directory to repo root.
4. Run: `Rscript run_revision_pack.R`

## Notes
- Public source datasets should be obtained from GEO/FinnGen as documented in manuscript and README files.
- Manual negative-control files should be placed under `manual_data/negative_controls/` following README instructions.
