#' buildRmdBundle_alt
#'
#' @description Modified shinymeta function to build a Rmd bundle following
#' ENCORE guidelines.
#' Original authors:  	Joe Cheng [aut], Carson Sievert ORCID iD [cre, aut], RStudio [cph]
#' Original source: https://github.com/rstudio/shinymeta/blob/main/R/report.R
#' Original license: GPL-3 (see https://cran.r-project.org/package=shinymeta)
#' @return NULL
#' @import shinymeta
#' @noRd
buildRmdBundle_alt<-function (report_template, output_zip_path, vars = list(), include_files = list(),
                              render = TRUE, render_args = list(), golem_renv=FALSE)
{
  force(report_template)
  force(vars)
  shinymeta_with_progress_obj(function(progress) {
    progress$set(value = 0)
    progress$set(message = "Generating code")
    if (is.list(vars)) {
      vars <- lapply(vars, function(x) {
        if (is.language(x)) {
          paste(formatCode(x), collapse = "\n")
        }
        else {
          x
        }
      })
    }
    progress$set(value = 0.1)
    progress$set(message = "Expanding Rmd template")
    rmd_source <- shinymeta_knit_expand_safe(report_template, vars = vars)
    rmd_filename <- shinymeta_template_rename(report_template, "Rmd")
    build_bundle_alt(rmd_source, rmd_filename, output_zip_path,
                     include_files = include_files, render = render, render_args = render_args,
                     progress = progress,  golem_renv= golem_renv)
  })
}


#' build_bundle_alt
#'
#' @description Modified shinymeta function to build a Rmd bundle following
#' ENCORE guidelines.
#' Original authors:  	Joe Cheng [aut], Carson Sievert ORCID iD [cre, aut], RStudio [cph]
#' Original source: https://github.com/rstudio/shinymeta/blob/main/R/report.R
#' Original license: GPL-3 (see https://cran.r-project.org/package=shinymeta)
#' @return NULL
#' @import shinymeta
#' @import attachment
#' @importFrom dockerfiler dock_from_desc
#' @import golem
#' @noRd
build_bundle_alt <-function (input_src, input_filename, output_zip_path, include_files = list(),
                             render = TRUE, render_args = list(), progress, golem_renv=FALSE)
{
  force(input_src)
  force(input_filename)
  force(output_zip_path)
  force(include_files)
  force(render)
  force(render_args)
  progress$set(value = 0.2)
  progress$set(message = "Adding items to zip archive")
  x <-  shinymeta_zip_archive()
  dest_filename_full <- fs::path(shinymeta_archive_basedir(x), input_filename)
  writeLines(input_src, dest_filename_full)
  shinymeta_add_items(x, !!!include_files)
  progress$set(value = 0.3)
  if (render) {
    progress$set(message = "Rendering report")
    shinymeta_render_with_args(dest_filename_full, render_args)
  }
  newdir=fs::path(shinymeta_archive_basedir(x), "Notebook")
  dir.create(newdir, showWarnings = FALSE)
  file.copy(dest_filename_full, fs::path(newdir, input_filename))
  file.remove(dest_filename_full)
  if(golem_renv) {
    progress$set(value = 0.8)
    progress$set(message = "Exporting package and renv")
    newdir=fs::path(shinymeta_archive_basedir(x), "0_SoftwareEnvironment/R")
    dir.create(newdir, showWarnings = FALSE, recursive = TRUE)
    suppressWarnings(
      suppressMessages(
        utils::capture.output(
          golem::add_dockerfile_with_renv(output_dir = newdir),
          type = "output"
        )
      )
    )
    file.rename(from=paste0(newdir, "/renv.lock.prod"),
                to=paste0(newdir, "/renv.lock"))
  }
  progress$set(value = 0.9)
  progress$set(message = "Compressing bundle")
  archive <-  shinymeta_build_archive(x, output_zip_path)
  progress$set(value = 1)
  archive
}

#' shinymeta_zip_archive
#'
#' @description Original code from package Shinymeta, internal function zip_archive.
#' Original authors:  	Joe Cheng [aut], Carson Sievert ORCID iD [cre, aut], RStudio [cph]
#' Original source: https://github.com/rstudio/shinymeta/blob/main/R/archive.R
#' Original license: GPL-3 (see https://cran.r-project.org/package=shinymeta)
#' @return structure with dir for zip file
#' @noRd
shinymeta_zip_archive <- function(temp_dir = NULL) {
  if (!is.null(temp_dir) && (!is.character(temp_dir) || length(temp_dir) != 1)) {
    stop("temp_dir must be a single-element character vector")
  }

  if (is.null(temp_dir)) {
    temp_dir <- tempfile("archive")
    fs::dir_create(temp_dir, mode = "u=rwx,go=")
  } else if (!dir.exists(temp_dir)) {
    stop("temp_dir directory does not exist")
  }

  structure(
    list(basedir = temp_dir),
    class = "pending_zip_archive"
  )
}

#' shinymeta_add_items
#'
#' @description Original code from package Shinymeta, internal function add_items.
#' Original authors:  	Joe Cheng [aut], Carson Sievert ORCID iD [cre, aut], RStudio [cph]
#' Original source: https://github.com/rstudio/shinymeta/blob/main/R/archive.R
#' Original license: GPL-3 (see https://cran.r-project.org/package=shinymeta)
#' @return structure with dir for zip file
#' @noRd
shinymeta_add_items <- function(x, ...) {
  stopifnot(inherits(x, "pending_zip_archive"))

  include_files <- rlang::dots_list(..., .homonyms = "last", .check_assign = TRUE)

  if (is.null(names(include_files))) {
    names(include_files) <- as.character(include_files)
  }

  mapply(names(include_files), include_files, FUN = function(to, from) {
    if (nchar(from) == 0) {
      from <- to
    }
    shinymeta_add_item(x, from, to)
    NULL
  })

  x
}


#' shinymeta_archive_basedir
#'
#' @description Adapted code from package Shinymeta, internal function archive_basedir
#' Original authors:  	Joe Cheng [aut], Carson Sievert ORCID iD [cre, aut], RStudio [cph]
#' Original source: https://github.com/rstudio/shinymeta/blob/main/R/archive.R
#' Original license: GPL-3 (see https://cran.r-project.org/package=shinymeta)
#' @return baseline dir for zip file
#' @noRd
shinymeta_archive_basedir <- function(x) {
  stopifnot(inherits(x, "pending_zip_archive"))

  x[["basedir"]]
}


#' shinymeta_template_rename
#'
#' @description Original code from package Shinymeta, internal function template_rename
#' Original authors:  	Joe Cheng [aut], Carson Sievert ORCID iD [cre, aut], RStudio [cph]
#' Original source: https://github.com/rstudio/shinymeta/blob/main/R/report.R
#' Original license: GPL-3 (see https://cran.r-project.org/package=shinymeta)
#' @return filename
#' @noRd
shinymeta_template_rename <- function(input_template, extension = "Rmd") {
  stopifnot(is.character(extension) && length(extension) == 1 && identical(TRUE, nzchar(extension)))

  filename <- fs::path_ext_remove(fs::path_file(input_template))
  if (tolower(fs::path_ext(filename)) == tolower(extension)) {
    filename
  } else {
    paste0(filename, ".", extension)
  }
}

#' shinymeta_build_archive
#'
#' @description Adapted code from package Shinymeta, internal function build_archive
#' Original authors:  	Joe Cheng [aut], Carson Sievert ORCID iD [cre, aut], RStudio [cph]
#' Original source: https://github.com/rstudio/shinymeta/blob/main/R/archive.R
#' Original license: GPL-3 (see https://cran.r-project.org/package=shinymeta)
#' @return NULL
#' @noRd
shinymeta_build_archive <- function(x, output_file) {

  basedir <- shinymeta_archive_basedir(x)

  olddir <- getwd()
  setwd(basedir)
  on.exit(setwd(olddir))

  utils::zip(fs::path_abs(output_file, olddir), ".",
             flags = "-rq")
  invisible(output_file)
}

#' shinymeta_with_progress_obj
#'
#' @description Original code from package Shinymeta, internal function with_progress_obj
#' Original authors:  	Joe Cheng [aut], Carson Sievert ORCID iD [cre, aut], RStudio [cph]
#' Original source: https://github.com/rstudio/shinymeta/blob/main/R/report.R
#' Original license: GPL-3 (see https://cran.r-project.org/package=shinymeta)
#' @return session progress
#' @noRd
shinymeta_with_progress_obj <- function(callback) {
  # Note that `session` may be NULL.
  session <- shiny::getDefaultReactiveDomain()
  if (!is.null(session$userData$shinymeta_last_progress)) {
    # If the last progress object we created for this session is still visible,
    # close it. This would be in the case of an error, we never auto-dismiss in
    # that case.
    suppressWarnings(session$userData$shinymeta_last_progress$close())
  }

  progress <- shinymeta_make_progress()
  # Register our newly created progress object.
  session$userData$shinymeta_last_progress <- progress

  tryCatch(shiny::captureStackTraces({
    callback(progress)
    progress$close()
    session$userData$shinymeta_last_progress <- NULL
  }), error = function(err) {
    progress$set(value = 1, message = "An error has occurred:",
                 detail = conditionMessage(err))
    stop(err)
  })
}

#' shinymeta_make_progress
#'
#' @description Original code from package Shinymeta, internal function make_progress
#' Original authors:  	Joe Cheng [aut], Carson Sievert ORCID iD [cre, aut], RStudio [cph]
#' Original source: https://github.com/rstudio/shinymeta/blob/main/R/report.R
#' Original license: GPL-3 (see https://cran.r-project.org/package=shinymeta)
#' @return session progress
#' @noRd
shinymeta_make_progress <- function(...) {
  session <- shiny::getDefaultReactiveDomain()
  if (!is.null(session)) {
    shiny::Progress$new(session = session, ...)
  } else {
    # Return a dummy progress object
    nothing <- function(...) {}
    list(
      set = nothing,
      inc = nothing,
      getMin = nothing,
      getMax = nothing,
      getValue = nothing,
      close = nothing,
      clone = nothing
    )
  }
}



#' shinymeta_knit_expand_safe
#'
#' @description Original code from package Shinymeta, internal function knit_expand_safe
#' Original authors:  	Joe Cheng [aut], Carson Sievert ORCID iD [cre, aut], RStudio [cph]
#' Original source: https://github.com/rstudio/shinymeta/blob/main/R/utils.R
#' Original license: GPL-3 (see https://cran.r-project.org/package=shinymeta)
#' @return knitr::knit_expand
#' @noRd
shinymeta_knit_expand_safe <- function(file, vars = list(), text = xfun::read_utf8(file), delim = c("{{", "}}")) {
  # The approach we take here is to detect all knitr md patterns before and
  # after expansion, and fail if anything was either added or removed. We tried
  # just testing the output of each {{expansion}} for the patterns, but, that
  # doesn't catch cases where an inline.code is started in one expansion and
  # finished in another (see test in test-report.R).

  # Code chunk delimiter regexes
  # TODO: Can we assume `md`?
  patterns <- unname(unlist(knitr::all_patterns$md))

  matches_before <- shinymeta_count_matches_by_pattern(text, patterns)

  # Create an environment that contains nothing but the variables we want to
  # make available for template expansion, plus .GlobalEnv.
  eval_envir <- list2env(vars, parent = .GlobalEnv)

  # Use a knitr hook to ensure that only the ... arguments plus stuff in the
  # global environment are available when evaluating {{/}} expressions.
  orig_eval_inline <- knitr::knit_hooks$get("evaluate.inline")
  knitr::knit_hooks$set(evaluate.inline = function(code, envir) {
    # ignore envir, as it includes the parent frame of `knit_expand` which we
    # explicitly do not want to be used for evaluation--only ... arguments to
    # knit_expand_safe should be used.
    orig_eval_inline(code, eval_envir)
  })
  on.exit(knitr::knit_hooks$set(evaluate.inline = orig_eval_inline), add = TRUE)

  res <- knitr::knit_expand(text = text, delim = delim)

  matches_after <- shinymeta_count_matches_by_pattern(xfun::split_lines(res), patterns)

  if (!identical(matches_before, matches_after)) {
    # The process of knit_expand-ing introduced new (or removed existing?) code
    # chunks
    stop("Can't build report--user input values must not contain code chunk delimiters")
  }

  res
}



#' shinymeta_render_with_args
#'
#' @description Original code from package Shinymeta, internal function render_with_args
#' Original authors:  	Joe Cheng [aut], Carson Sievert ORCID iD [cre, aut], RStudio [cph]
#' Original source: https://github.com/rstudio/shinymeta/blob/main/R/report.R
#' Original license: GPL-3 (see https://cran.r-project.org/package=shinymeta)
#' @return NULL
#' @noRd
shinymeta_render_with_args <- function(input_file, render_args = list(), switch_dirs = TRUE, fork = TRUE) {

  if (switch_dirs) {
    old_wd <- getwd()
    setwd(fs::path_dir(input_file))
    on.exit(setwd(old_wd))
  }

  if (fork) {
    callr::r(
      function(...) rmarkdown::render(...),
      # https://github.com/rstudio/rmarkdown/issues/1204#issuecomment-344884823
      args = c(list(input_file, envir = globalenv()), render_args)
    )
  } else {
    do.call(rmarkdown::render, c(list(input_file), render_args), quote = TRUE)
  }
}


#' shinymeta_count_matches_by_pattern
#'
#' @description Original code from package Shinymeta, internal function count_matches_by_pattern
#' Original authors:  	Joe Cheng [aut], Carson Sievert ORCID iD [cre, aut], RStudio [cph]
#' Original source: https://github.com/rstudio/shinymeta/blob/main/R/utils.R
#' Original license: GPL-3 (see https://cran.r-project.org/package=shinymeta)
#' @return vector of length(pattern)
#' @noRd
shinymeta_count_matches_by_pattern <- function(string, pattern) {
  vapply(pattern, function(regex) {
    matches <- stringr::str_locate_all(string, regex)
    sum(vapply(matches, nrow, integer(1)))
  }, integer(1), USE.NAMES = FALSE)
}


#' shinymeta_add_item
#'
#' @description Adapted code from package Shinymeta, internal function add_item
#' Original authors:  	Joe Cheng [aut], Carson Sievert ORCID iD [cre, aut], RStudio [cph]
#' Original source: https://github.com/rstudio/shinymeta/blob/main/R/archive.R
#' Original license: GPL-3 (see https://cran.r-project.org/package=shinymeta)
#' @return vector of length(pattern)
#' @noRd
shinymeta_add_item <- function(x, source_file, target_file) {
  stopifnot(inherits(x, "pending_zip_archive"))

  if (!is.character(source_file) || length(source_file) != 1) {
    stop("source_file must be a single-element character vector")
  }
  if (!is.character(target_file) || length(target_file) != 1) {
    stop("target_file must be a single-element character vector")
  }

  if (fs::is_absolute_path(target_file)) {
    stop("target_file must be a relative path")
  }

  full_src <- fs::path_abs(source_file)

  basedir <- shinymeta_archive_basedir(x)

  if (fs::dir_exists(full_src)) {
    full_dest <- fs::path(basedir, target_file)
    fs::dir_copy(full_src, full_dest)
  } else {
    if (grepl("[/\\]$", target_file)) {
      # If source is a file, but target is a directory name, ensure
      # that the file gets copied into the target, rather than as
      # the target. Without this line, fs::file_copy would treat
      # the target as a filename (it would strip off the slash).
      target_file <- fs::path(target_file, fs::path_file(source_file))
    }
    full_dest <- fs::path(basedir, target_file)

    if (!fs::path_dir(target_file) %in% c("", "."))
      fs::dir_create(fs::path_dir(full_dest), recurse = TRUE)
    fs::file_copy(full_src, full_dest)
  }

  x
}
