#' Run the AbSolution Shiny Application
#' Launches the interactive AbSolution Shiny application for exploring
#' B-cell and T-cell receptor repertoire data. The application provides
#' tools for loading AIRR-formatted datasets, inspecting sequence
#' characteristics, visualizing clonotype distributions, and generating
#' summary reports.
#'
#' Internally, `run_app()` calls [`shiny::shinyApp()`] to start the
#' application, using the UI defined in [\code{app_ui()}] and the server
#' logic defined in [\code{app_server()}].
#'
#' This function is the entry point that end-users should call after
#' installing the package.
#'
#' @param verbose Logical. If \code{FALSE} (default), suppresses all warnings
#' and messages during app execution. Set to \code{TRUE} to see full console
#' output for debugging.
#' @param ... arguments to pass to golem_opts.
#' See `?golem::get_golem_options` for more details.
#'
#' @inheritParams shiny::shinyApp
#'
#' @export
#' @importFrom shiny shinyApp
#' @importFrom golem with_golem_options
#' @examples
#' if (interactive()) {
#'   ### Launch the full application with default settings:
#'   run_app()
#' }
run_app <- function(verbose = FALSE, onStart = NULL, options = list(),
                    enableBookmarking = "server", uiPattern = "/", ...) {
  # set.seed(1234)

  with_golem_options(
    app = shinyApp(
      ui = app_ui,
      server = app_server,
      onStart = onStart,
      options = options,
      enableBookmarking = enableBookmarking,
      uiPattern = uiPattern
    ),
    golem_opts = list(verbose = verbose, ...)
  )
}

