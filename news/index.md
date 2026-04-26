# Changelog

## AbSolution 1.0.1

- Fixed installation from a clean R environment: added
  `Additional_repositories` to `DESCRIPTION` so that Bioconductor
  packages required by `alakazam` (specifically `GenomicAlignments`) are
  resolved automatically during
  [`install.packages()`](https://rdrr.io/r/utils/install.packages.html).

## AbSolution 1.0.0

CRAN release: 2026-04-21

- Initial CRAN release.
- AIRR-Seq parsing with automatic germline reconstruction.
- Sequence feature extraction: nucleotide/amino acid properties,
  mutation analysis, physicochemical descriptors.
- Interactive Shiny interface with PCA/UMAP projections, violin plots,
  UpSet plots.
- Clonal exploration: multiple clonotype definitions, shared clone
  detection, dominance analysis.
- Differential variable selection between user-defined groups.
- ENCORE-compatible reproducible export with Docker support.
