## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## 12/04/2026 - Resubmission: reviewer comments addressed

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
