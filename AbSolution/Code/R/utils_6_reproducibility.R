#' buildRmdBundle_alt
#'
#' @description Modified shinymeta function to build a Rmd bundle following
#' ENCORE guidelines.
#'
#' @return NULL
#' @import shinymeta
#' @noRd
buildRmdBundle_alt<-function (report_template, output_zip_path, vars = list(), include_files = list(),
                              render = TRUE, render_args = list(), golem_renv=F)
{
  force(report_template)
  force(vars)
  shinymeta:::with_progress_obj(function(progress) {
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
    rmd_source <- shinymeta:::knit_expand_safe(report_template, vars = vars)
    rmd_filename <- shinymeta:::template_rename(report_template, "Rmd")
    build_bundle_alt(rmd_source, rmd_filename, output_zip_path,
                     include_files = include_files, render = render, render_args = render_args,
                     progress = progress,  golem_renv= golem_renv)
  })
}


#' build_bundle_alt
#'
#' @description Modified shinymeta function to build a Rmd bundle following
#' ENCORE guidelines.
#'
#' @return NULL
#' @import shinymeta
#' @import golem
#' @noRd
build_bundle_alt <-function (input_src, input_filename, output_zip_path, include_files = list(),
                             render = TRUE, render_args = list(), progress, golem_renv=F)
{
  force(input_src)
  force(input_filename)
  force(output_zip_path)
  force(include_files)
  force(render)
  force(render_args)
  progress$set(value = 0.2)
  progress$set(message = "Adding items to zip archive")
  x <-  shinymeta:::zip_archive()
  dest_filename_full <- fs::path(shinymeta:::archive_basedir(x), input_filename)
  writeLines(input_src, dest_filename_full)
  shinymeta:::add_items(x, !!!include_files)
  progress$set(value = 0.3)
  if (render) {
    progress$set(message = "Rendering report")
    shinymeta:::render_with_args(dest_filename_full, render_args)
  }
  newdir=fs::path(shinymeta:::archive_basedir(x), "Notebook")
  dir.create(newdir)
  file.copy(dest_filename_full, fs::path(newdir, input_filename))
  file.remove(dest_filename_full)
  if(golem_renv) {
    progress$set(value = 0.8)
    progress$set(message = "Exporting package and renv")
    newdir=fs::path(shinymeta:::archive_basedir(x), "0_SoftwareEnvironment/R")
    dir.create(newdir)
    golem::add_dockerfile_with_renv(output_dir = newdir)
  }
  progress$set(value = 0.9)
  progress$set(message = "Compressing bundle")
  archive <-  shinymeta:::build_archive(x, output_zip_path)
  progress$set(value = 1)
  archive
}



