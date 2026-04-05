# Run the AbSolution Shiny Application Launches the interactive AbSolution Shiny application for exploring B-cell and T-cell receptor repertoire data. The application provides tools for loading AIRR-formatted datasets, inspecting sequence characteristics, visualizing clonotype distributions, and generating summary reports.

Internally, \`run_app()\` calls \[\`shiny::shinyApp()\`\] to start the
application, using the UI defined in \[`app_ui()`\] and the server logic
defined in \[`app_server()`\].

## Usage

``` r
run_app(
  verbose = FALSE,
  onStart = NULL,
  options = list(),
  enableBookmarking = "server",
  uiPattern = "/",
  ...
)
```

## Arguments

- verbose:

  Logical. If `FALSE` (default), suppresses all warnings and messages
  during app execution. Set to `TRUE` to see full console output for
  debugging.

- onStart:

  A function that will be called before the app is actually run. This is
  only needed for `shinyAppObj`, since in the `shinyAppDir` case, a
  `global.R` file can be used for this purpose.

- options:

  Named options that should be passed to the `runApp` call (these can be
  any of the following: "port", "launch.browser", "host", "quiet",
  "display.mode" and "test.mode"). You can also specify `width` and
  `height` parameters which provide a hint to the embedding environment
  about the ideal height/width for the app.

- enableBookmarking:

  Can be one of `"url"`, `"server"`, or `"disable"`. The default value,
  `NULL`, will respect the setting from any previous calls to
  [`enableBookmarking()`](https://rdrr.io/pkg/shiny/man/enableBookmarking.html).
  See
  [`enableBookmarking()`](https://rdrr.io/pkg/shiny/man/enableBookmarking.html)
  for more information on bookmarking your app.

- uiPattern:

  A regular expression that will be applied to each `GET` request to
  determine whether the `ui` should be used to handle the request. Note
  that the entire request path must match the regular expression in
  order for the match to be considered successful.

- ...:

  arguments to pass to golem_opts. See \`?golem::get_golem_options\` for
  more details.

## Details

This function is the entry point that end-users should call after
installing the package.

## Examples

``` r
if (interactive()) {
  ### Launch the full application with default settings:
  run_app()
}
```
