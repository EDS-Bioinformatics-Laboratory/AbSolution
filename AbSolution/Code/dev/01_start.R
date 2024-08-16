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
  pkg_title = "AbSolution: BRC and TCR Fab sequence feature analysis", # The Title of the package containing the App
  pkg_description = "A flexible and scalable Shiny app for BRC and TCR Fab sequence feature analysis.", # The Description of the package containing the App
  author_first_name = "Rodrigo", # Your First Name
  author_last_name = "García-Valiente", # Your Last Name
  author_email = "r.garciavaliente@amsterdamumc.nl", # Your Email
  repo_url = NULL, # The URL of the GitHub Repo (optional),
  pkg_version = "0.0.0.9000" # The Version of the package containing the App
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
