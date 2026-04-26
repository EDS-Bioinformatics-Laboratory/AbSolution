## 26/04/2026 - Patch release

This is a patch release fixing an installation failure reported in clean R
environments. `alakazam` depends on `GenomicAlignments` (Bioconductor), but
`install.packages()` does not propagate `Additional_repositories` transitively
from a dependency's DESCRIPTION. The fix adds Bioconductor to
`Additional_repositories` in this package's own DESCRIPTION so that the
dependency is resolved correctly on first install.

Apologies for the trouble, this was something I did not expect.

### Test environments

* local: Ubuntu 24.04.3 LTS, R 4.5.3.  0 errors | 0 warnings | 1 note
* win-builder: release and devel. 0 errors | 0 warnings | 1 note

#### R CMD check results

0 errors | 0 warnings | 1 note

##### Notes

###### Maintainer / Days since last update
Expected for a resubmission within the same week.

###### Availability using Additional_repositories specification: ? ?
`R CMD check` reports `?` for `https://bioconductor.org/packages/release/bioc`
on all platforms (win-builder, Ubuntu). This is a known behaviour: the check
tool cannot verify availability of Bioconductor repos deterministically, but
the repository is reachable and the packages are available:

    available.packages(repos = "https://bioconductor.org/packages/release/bioc")
    # "GenomicAlignments" %in% pkgs[, "Package"] -> TRUE


## 17/04/2026 - Resubmission: reviewer comments addressed

* **Title**: Shortened to under 65 characters.

* **Description references**: Added the GitHub URL in the required
  `<https://...>` format. This will be updated with a proper
  `authors (year) <doi:...>` reference once the associated manuscript
  is published.

* **TRUE/FALSE**: Replaced all instances of `T` and `F` with `TRUE`
  and `FALSE` throughout all R source files.

* **.GlobalEnv**: Replaced the `.GlobalEnv` reference in
  `R/utils_6_reproducibility.R` (line 298) with `baseenv()` to avoid
  modifying the user's global workspace. As AbSolution is a
  Shiny application, it may fall under the Shiny exception to this rule,
  but we chose `baseenv()` to depend only on explicitly passed variables.

* **Executable examples**: Added `\donttest{}` blocks to all 13 exported
  functions. Functions other than `run_app()` are tagged with
  `\keyword{internal}` and exported solely for `::` access by shinymeta
  during report rendering, They require a live Shiny reactive context
  and real AIRR-seq data to run meaningfully. List of related functions:
- `parse_AIRRSeq_file()`
- `Feature_1()`
- `Feature__dataset()`
- `merge_FBMs()`
- `filter_merged()`
- `big_PCA()`
- `calculate_clone()`
- `draw_violinplots()`
- `draw_upsetplot()`
- `draw_sharedclonesplot()`
- `draw_feature_violinplot()`
- `Ab_palette()`
- `validate_before_run()`

The package has been tested with no ERRORs or WARNINGs
(only the expected NOTE for new submission) on: my
development machine (Ubuntu 24.04, R 4.5.3),
win-builder with R-devel (Windows Server 2022, R
4.6.0 alpha), and via rhub::rhub_check() on linux
(ubuntu-latest, R-devel), macos-arm64 (macos-latest,
R-devel), windows (windows-latest, R-devel), and
atlas (Fedora 42, R-devel).

Thanks for your work once again!

## 06/04/2026 - First submission

AbSolution is a Shiny application for interactive
analysis of B-cell and T-cell receptor repertoires.
We are aware of the number of Imports — this is
because the app provides an end-to-end pipeline
within a single tool, covering sequence parsing
(Biostrings, pwalign, IRanges, alakazam), scalable
computation on large datasets (bigstatsr,
bigparallelr), interactive visualization (plotly, DT,
reactable, upsetjs), the Shiny UI and its extensions
(bs4Dash, shinyFiles, shinyWidgets, shinymeta,
shinymanager), and reproducible export with
containerization support (golem, dockerfiler,
rmarkdown). I confirm none can be removed or moved to
suggests without breaking core functionality.

The package has been tested with no ERRORs or WARNINGs
(only the expected NOTE for new submission) on: my
development machine (Ubuntu 24.04, R 4.5.3),
win-builder with R-devel (Windows Server 2022, R
4.6.0 alpha), and via rhub::rhub_check() on linux
(ubuntu-latest, R-devel), macos-arm64 (macos-latest,
R-devel), windows (windows-latest, R-devel), and
atlas (Fedora 42, R-devel).

Thanks for your work!

