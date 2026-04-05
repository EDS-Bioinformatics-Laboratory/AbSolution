# Building a Prod-Ready, Robust Shiny Application.
#
# README: each step of the dev files is optional, and you don't have to
# fill every dev scripts before getting started.
# 01_start.R should be filled at start.
# 02_dev.R should be used to keep track of your development during the project.
# 03_deploy.R should be used once you need to deploy your app.
#
#
########################################
#### CURRENT FILE: ON START SCRIPT #####
########################################

## Fill the DESCRIPTION ----
## Add meta data about your application
##
## /!\ Note: if you want to change the name of your app during development,
## either re-run this function, call golem::set_golem_name(), or don't forget
## to change the name in the app_sys() function in app_config.R /!\
##
golem::fill_desc(
  pkg_name = "AbSolution", # The Name of the package containing the App
  pkg_title = "B-Cell and T-Cell Variable Region Sequence Feature Analysis’", # The Title of the package containing the App
  pkg_description = paste(
    "An interactive and scalable Shiny application for analyzing",
    "sequence features of B-cell and T-cell receptor repertoires.",
    "Supports identification and visualization of variable region",
    "characteristics, clonotypes, and mutation patterns from",
    "immunoglobulin and TCR sequence data."
  ),
  authors = c(person(given = "Rodrigo", family =  "García-Valiente",
                     email = "r.garciavaliente@amsterdamumc.nl", role = c("cre",
                                                                          "aut"),
                     comment =  c(ORCID = "0000-0003-0444-5587")),
              person(given = "Charisios", family =  "Triantafyllou",
                     email = "charisios.triantafyllou@genyo.es", role = c("ctb",
                                                                          "aut"),
                     comment =  c(ORCID = "0000-0002-0283-141X")),
              person(given = "Antoine", family =  "van Kampen",
                     email = "a.h.c.vankampen@amsterdamumc.nl", role = c("aut",
                                                                          "ths"),
                     comment =  c(ORCID = "0000-0003-1025-7232"))),
  repo_url = "https://github.com/EDS-Bioinformatics-Laboratory/AbSolution", # The URL of the GitHub Repo (optional),
  pkg_version = "1.0.0" # The Version of the package containing the App
)

## Set {golem} options ----
golem::set_golem_options()

## Install the required dev dependencies ----
golem::install_dev_deps()

## Create Common Files ----
## See ?usethis for more information
usethis::use_gpl_license(version = 3, include_future = TRUE) # You can set another license here


usethis::use_readme_rmd(open = FALSE)  # NOT FOR NOW
# devtools::build_readme()

## Init Testing Infrastructure ----
## Create a template for tests
golem::use_recommended_tests() # NOT FOR NOW

## Favicon ----
# If you want to change the favicon (default is golem's one)
golem::use_favicon("extra/AbSolution.ico") # path = "path/to/ico". Can be an online file.
# golem::remove_favicon() # Uncomment to remove the default favicon

## Add helper functions ----
golem::use_utils_ui(with_test = TRUE)
golem::use_utils_server(with_test = TRUE)

# You're now set! ----

# go to dev/02_dev.R
rstudioapi::navigateToFile("dev/02_dev.R")

# EXTRA: ----
usethis::use_package("ggplot2", min_version=T)
usethis::use_pipe()
