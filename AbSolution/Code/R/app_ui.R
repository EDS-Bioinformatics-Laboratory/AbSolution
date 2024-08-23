#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @import DT
#' @import shinydashboard
#' @import shinydashboardPlus
#' @import shinyWidgets
#' @import dashboardthemes
#' @import plotly
#' @import shinycssloaders
#' @import shinyFiles
#' @import shinyjs
#' @import shinymanager
#' @import shinythemes
#' @import shinyWidgets
#' @import upsetjs
#' @import utils
#' @noRd
app_ui <- function(request) {
    tagList(
        golem_add_external_resources(), dashboardPage(
            dashboardHeader(
                title = logo_absolution, shinydashboard::dropdownMenu(
                  type = "messages", messageItem(from = "Rodrigo Garcia", message = "Have fun using this app!")
              )
            ),
            dashboardSidebar(
                sidebarMenu(
                  menuItem("Home", tabName = "home", icon = icon("home")),
                  menuItem(
                    "AbSolution", icon = icon("chart-bar"),
                    tabName = "pipeline"
                ),
                  menuItem(
                    "About us", icon = icon("info"),
                    tabName = "About"
                )
              )
            ),
            dashboardBody(
                useShinyjs(), tags$head(tags$link(rel = "shortcut icon", href = "www/img/favicon.ico")),
                shinyDashboardThemes(theme = "poor_mans_flatly"),
                tabItems(
                  tabItem(
                    tabName = "home", fluidRow(
                      userBox(
                        title = NULL, type = 2, src = NULL, color = "yellow", width = 12,
                        boxToolSize = "lg", closable = F, collapsible = F, HTML(
                          "<img src=\"www/img/AbSolution_logo.png\" alt=\"AB Solution\" width=\"653\" height=\"150\">"
                      ),
                        footer = HTML(
                          "<b>A tool for antibody Fab sequence feature analysis</b>
                    <br>
                    <p align='justify'>AbSolution is a fancy pipeline. <br><br>Contact: <a href='mailto:r.garciavaliente@amsterdamumc.nl'> Biolab</a>
                                  </p>
                                  "
                      )
                    )
                  )
                ),
                  tabItem(
                    tabName = "pipeline", h2("Process your data"),
                    HTML(
                      "If you have previous files, you don't have to redo all the steps. Process and upload your files where/when needed. But do provide the sample information.<br> <br>"
                  ),
                    tabsetPanel(
                      id = "menutabset", type = "tabs", tabPanel(
                        "0.Project information", wellPanel(
                          h3("Upload your sample information table"),
                          radioButtons(
                            "TCRBCR_input_file", "The repertoire to analyze is:",
                            choiceNames = list("BCR", "TCR"),
                            choiceValues = list("BCR", "TCR"),
                            selected = "BCR"
                        ),
                          radioButtons(
                            "radiobutton_input_file", "The input file type is:",
                            choiceNames = list("AIRR-Seq format"),
                            choiceValues = list("airrseq"),
                            selected = "airrseq"
                        ),
                          uiOutput("Conditional_Action_Filetype_description"),
                          uiOutput("Conditional_Action_Filetype_upload"),
                          HTML("<br>"),
                          DT::dataTableOutput("raw_sample_file_out"),
                          HTML("<br><br><h3>Select work folder</h3>"),
                          shinyDirButton(
                            "base_folder", "Select the folder where the data will be saved",
                            title = "Please select the work folder:", buttonType = "default",
                            class = NULL, icon = icon("folder", lib = "font-awesome"),
                            multiple = F
                        ),
                          HTML("<br><br><br><h3>Next step</h3>"),
                          uiOutput("Conditional_Action_Move_to_1"),
                          HTML("<br><br>"),
                          uiOutput("Conditional_Action_Move_to_Analysis")
                      )
                    ),
                      tabPanel(
                        "1.AIRR-Seq conversion", wellPanel(
                          h3("Upload your AIRR-Seq  files"),
                          HTML("Select the main folder that contains the different AIRR-Seq files.<br> <br>"),
                          shinyDirButton(
                            "preinfolder_AIRR", "Select main folder", title = "Please select the main folder:",
                            buttonType = "default", class = NULL, icon = icon("folder", lib = "font-awesome"),
                            multiple = F
                        ),
                          HTML("<br>"),
                          DT::dataTableOutput("folder_information_AIRR"),
                          HTML("<br>"),
                          h3("Set parameters for germline reconstruction"),
                          materialSwitch(
                            "Dgene_reconstruct_airr", "Reconstruct using the D gene information instead of keeping that part of the CDR3 as is.",
                            value = T
                        ),
                          materialSwitch(
                            "FWR1partial_airr", "Is there a partial FWR1 in the sequences?",
                            value = FALSE
                        ),
                          materialSwitch(
                            "FWR4partial_airr", "Is there a partial FWR4 in the sequences?",
                            value = FALSE
                        ),
                          materialSwitch(
                            "C_region_included_airr", "Is part or totatily of the C region included in the sequence?",
                            value = FALSE
                        ),
                          HTML("<br>"),
                          h3("Process AIRR-Seq files and continue to the next step"),
                          HTML("<br>"),
                          uiOutput("Conditional_Action_Preprocess_AIRR"),
                          shinyjs::hidden(
                            div(
                              id = "pb_AIRR_vis", shinyWidgets::progressBar(
                                id = "pb_AIRR", value = 0, total = 100, title = "Processing...",
                                display_pct = TRUE
                            )
                          )
                        ),
                          uiOutput("Conditional_Action_Move_to_2_AIRR")
                      )
                    ),
                      tabPanel(
                        "2.Sequence feature determination", wellPanel(
                          actionButton(
                            "Feature_determination", "Feature determination", icon("forward"),
                            style = "color: #fff; background-color: #18bc9c"
                        ),
                          div(
                            id = "pb_Feature_vis", shinyWidgets::progressBar(
                              id = "pb_Feature", value = 0, total = 100, title = "Processing...",
                              display_pct = TRUE
                          )
                        ),
                          uiOutput("Conditional_Action_Move_to_3")
                      )
                    ),
                      tabPanel(
                        "1&2. Select your bigmem files", wellPanel(
                          shinyDirButton(
                            "bigmem_folder", "Select the folder with your .rds and .bk files",
                            title = "Please select the bigmem folder:", buttonType = "default",
                            class = NULL, icon = icon("folder", lib = "font-awesome"),
                            multiple = F
                        ),
                          HTML("<br>"),
                          DT::dataTableOutput("folder_information_AIRR_steps_1_to_3"),
                          HTML("<br><br><br>"),
                          uiOutput("make_dummy"),
                          conditionalPanel(
                            "input.dummy_dataset == true", selectizeInput(
                              "make_dummy_menu", "Make a dummy version of the dataset:",
                              choices = NULL, multiple = F, options = list(plugins = list("remove_button"))
                          )
                        ),
                          conditionalPanel(
                            "input.make_dummy_menu.length > 0 & input.dummy_dataset == true",
                            numericInput(
                              "dummy_ratio", "How many sequences per original sequence?",
                              1, min = 1, max = 10
                          ),
                            actionButton(
                              "final_dummy", "Produce the dummy dataset (click and wait until the table refreshes)",
                              icon("forward"),
                              style = "color: #fff; background-color: #18bc9c"
                          )
                        ),
                          HTML("<br><br><br>"),
                          uiOutput("Conditional_Action_Move_to_Analysis_Real")
                      )
                    ),
                      tabPanel(
                        "3.Dataset exploration and variable selection", tabsetPanel(
                          tabPanel(
                            "Menu", column(
                              4, wellPanel(
                                HTML(
                                  "Explore the dataset using both plots. The dataset and variables that you will use for the feature selection process are represented in the Selection plot. <br>"
                              ),
                                h3("1. Select your sequences"),
                                selectizeInput(
                                  "use_what", "Include sequences:", choices = c("Repertoire", "Reconstructed germline"),
                                  selected = c("Repertoire"),
                                  multiple = TRUE, options = list(plugins = list("remove_button"))
                              ),
                                selectizeInput(
                                  "use_productive_or_not", "That are:", choices = c("Productive"),
                                  selected = c("Productive"),
                                  multiple = TRUE, options = list(plugins = list("remove_button"))
                              ),
                                uiOutput("Sample_selection"),
                                sliderInput(
                                  inputId = "Rmut_filter", label = "Include only sequences with a number of R mutations between",
                                  min = 0, max = 100, value = c(0, 100),
                                  step = 1
                              ),
                                uiOutput("VJ_selection"),
                                HTML("<br>"),
                                h3("2. Select your variables"),
                                HTML(
                                  "<i>NOTE: Variables with NAs/Infinite/Zero variance are automatically removed.</i><br><br>"
                              ),
                                selectizeInput(
                                  "my_regions", "Select the region(s)", choices = c(
                                    "Whole", "FWR1", "CDR1", "FWR2", "CDR2", "FWR3",
                                    "CDR3", "FWR4"
                                ),
                                  selected = c("CDR3"),
                                  multiple = TRUE, options = list(plugins = list("remove_button"))
                              ),
                                selectizeInput(
                                  "my_var_elements", "Study features from", choices = c("NT", "AA"),
                                  selected = c("NT", "AA"),
                                  multiple = TRUE, options = list(plugins = list("remove_button"))
                              ),
                                selectizeInput(
                                  "my_vars", "Select the variable(s) class(es)",
                                  choices = c(
                                    "Length", "Composition", "NGly sites", "Hot/Cold motifs",
                                    "Leveshtein distance", "Substitutions", "Insertions",
                                    "Deletions", "Translocations", "Transitions and transversions",
                                    "Replacement and silent mutations", "Mutations from X to Y",
                                    "Peptide features"
                                ),
                                  selected = c(
                                    "Length", "Composition", "NGly sites", "Hot/Cold motifs",
                                    "Peptide features"
                                ),
                                  multiple = TRUE, options = list(plugins = list("remove_button"))
                              ),
                                selectizeInput(
                                  "my_vartypes", "Select the variable type(s)", choices = c("Baseline", "Germline diff"),
                                  selected = c("Baseline"),
                                  multiple = TRUE, options = list(plugins = list("remove_button"))
                              ),
                                uiOutput("individual_variables"),
                                materialSwitch(
                                  "use_UMAP", "Also represent the data in a UMAP (increases computational time)",
                                  value = FALSE, width = NULL, status = "primary"
                              ),
                                HTML("<br>"),
                                h3("3. Define the groups to study"),
                                materialSwitch(
                                  "work_as_categories", "I want to compare between groups",
                                  value = FALSE, width = NULL, status = "primary"
                              ),
                                conditionalPanel(
                                  "input.work_as_categories == true", uiOutput("Group_selection"),
                                  uiOutput("Group_comparison"),
                                  hr(), conditionalPanel(
                                    "input.work_as_categories == true & input.group_A.length > 0 & input.group_B.length > 0",
                                    materialSwitch(
                                      "use_sharedVDJ", "Include only sequences with V(D)J genes present in all the groups",
                                      value = FALSE, width = NULL, status = "primary"
                                  )
                                ),
                                  conditionalPanel(
                                    "input.work_as_categories == true & input.use_sharedVDJ == true & input.group_A.length > 0 & input.group_B.length > 0",
                                    materialSwitch(
                                      "VDJ_normalized_per_size", "Use the same number of sequences for each VJ combination on each group",
                                      value = T, width = NULL, status = "primary"
                                  ),
                                    conditionalPanel(
                                      "input.work_as_categories == true & input.use_sharedVDJ == true & input.group_A.length > 0 & input.group_B.length > 0 & input.VDJ_normalized_per_size == true",
                                      materialSwitch(
                                        "VDJ_maximize_clones", "Maximize the number of clones for each VJ combination on each group",
                                        value = T, width = NULL, status = "primary"
                                    ),
                                      conditionalPanel(
                                        "input.work_as_categories == true & input.use_sharedVDJ == true & input.group_A.length > 0 & input.group_B.length > 0 & input.VDJ_normalized_per_size== true & input.VDJ_maximize_clones == true",
                                        selectizeInput(
                                          "my_clone_def", "Select the clone definition",
                                          choices = c(NULL),
                                          selected = c(NULL),
                                          multiple = F, options = list(plugins = list("remove_button"))
                                      ),
                                    ),
                                      materialSwitch(
                                        "VDJ_normalized_per_sample", "And apply these rules also to each individual sample",
                                        value = F, width = NULL, status = "primary"
                                    )
                                  ),
                                    uiOutput("VDJ_subsetting")
                                ),
                                  conditionalPanel(
                                    "input.work_as_categories == true & input.group_A.length > 0 & input.group_B.length > 0",
                                    materialSwitch(
                                      "use_univlog", "Include only features that are correlated to the groups (0/1)",
                                      value = FALSE, width = NULL, status = "primary"
                                  ),
                                    conditionalPanel(
                                      "input.work_as_categories == true & input.use_univlog == true & input.group_A.length > 0 & input.group_B.length > 0",
                                      selectizeInput(
                                        "pval_type", "Select the P-value to use",
                                        choices = c("RAW", "Corrected by Bonferroni", "Corrected by B-H"),
                                        selected = c("Corrected by Bonferroni"),
                                        multiple = F, options = list(plugins = list("remove_button"))
                                    ),
                                      numericInput(
                                        "pval_cutoff", "P-value cut-off:", 0.05,
                                        min = 0, max = 1
                                    ),
                                      numericInput(
                                        "estimate_cutoff", "Estimate cut-off (absolute value):",
                                        2, min = 0
                                    ),
                                      selectizeInput(
                                        "number_selected_vars", "How many selected features to use?",
                                        choices = c(
                                          "All", "Top 10 (according to estimate)",
                                          "Top 50 (according to estimate)"
                                      ),
                                        selected = c("All"),
                                        multiple = F, options = list(plugins = list("remove_button"))
                                    )
                                  )
                                )
                              ),
                                h3("4. Update your changes"),
                                actionButton(
                                  "Update_selection", "Update Selection set, for ML",
                                  icon("forward"),
                                  style = "color: #fff; background-color: #18bc9c"
                              ),
                                actionButton(
                                  "Update_exploration", "Update Exploration set, compare with Selection",
                                  icon("forward"),
                                  style = "color: #fff; background-color: #CC2936"
                              )
                            )
                          ),
                            column(
                              8, HTML("<br><br>"),
                              HTML("Use the lasso tool or box select to get the sequence(s) information.<br><br>"),
                              fluidRow(
                                column(
                                  4, HTML("<h2>Selection plot</h2>"),
                                  fluidRow(
                                    column(
                                      6, selectInput(
                                        "Selection_plot_type", "Dimensionality Reduction:",
                                        c(PCA = "PCA")
                                    ),
                                      selectizeInput(
                                        "Selection_plot_type_dim", "Select the dimensions to plot",
                                        choices = c(1:5),
                                        selected = c(1, 2),
                                        multiple = TRUE, options = list(
                                          plugins = list("remove_button"),
                                          maxItems = 2
                                      )
                                    ),
                                      uiOutput("Plot_color")
                                  )
                                )
                              ),
                                column(8, plotlyOutput("Selection_plot"))
                            ),
                              fluidRow(
                                column(
                                  4, HTML("<h2>Exploration plot</h2>"),
                                  fluidRow(
                                    column(
                                      6, selectInput(
                                        "Exploration_plot_type", "Dimensionality Reduction:",
                                        c(PCA = "PCA")
                                    ),
                                      selectizeInput(
                                        "Exploration_plot_type_dim", "Select the dimensions to plot",
                                        choices = c(1:5),
                                        selected = c(1, 2),
                                        multiple = TRUE, options = list(
                                          plugins = list("remove_button"),
                                          maxItems = 2
                                      )
                                    ),
                                      uiOutput("Plot_color_expl")
                                  )
                                )
                              ),
                                column(8, plotlyOutput("Exploration_plot"))
                            ),
                              fluidRow(
                                column(6, uiOutput("Selection_plot_error")),
                                column(6, uiOutput("Exploration_plot_error"))
                            ),
                              fluidRow(
                                column(6, infoBoxOutput("RAM_S_Memory_Box", width = 12)),
                                column(6, infoBoxOutput("RAM_E_Memory_Box", width = 12))
                            ),
                              DT::dataTableOutput("Ab_table")
                          )
                        ),
                          tabPanel(
                            "About the variables", uiOutput("about_vars"),
                            tabsetPanel(
                              tabPanel(
                                "Used variables", fluidRow(
                                  column(
                                    10, HTML(
                                      "<br>Here you can see which variables are represented on each set.<br><br>The plots can take a few minutes to load!<br><br>"
                                  ),
                                    HTML("<h2>Selection set</h2>"),
                                    HTML("<h3>Whole sequence</h3>")
                                )
                              )
                            ),
                              tabPanel(
                                "Variable correlation", HTML(
                                  "<br>Here you can see which variables in the selection set are correlated >=0.8 for their non-zero values. Normalized and raw variables are not considered as correlated in the plot.<br><br>The plots can take a few minutes to load!<br><br>"
                              ),
                                htmlOutput("correlated_variables") %>%
                                  withSpinner(type = 6, color = "#cf395c")
                            )
                          )
                        ),
                          tabPanel(
                            "About the samples", HTML(
                              "<br>Here you can see information about the samples that are present on the selection set.<br><br>Remember you can have included (reconstructed) germline sequences!"
                          ),
                            uiOutput("about_samples")
                        ),
                          tabPanel(
                            "Add covariables", HTML("<br>Here you can add covariables to the selection dataset!<br><br>"),
                            uiOutput("about_meh"),
                            tabsetPanel(
                              tabPanel(
                                "Numerical covariables (e.g. gene expression)", HTML(
                                  "Use column names to define the sequence ID (column name Sequence_ID) and the covariables (use your own nomenclature)<br>"
                              ),
                                fileInput(
                                  "numerical_covariable_file", "Upload", multiple = F,
                                  accept = c(".tab", ".csv", ".txt")
                              ),
                                uiOutput("numerical_covariables_plot")
                            ),
                              tabPanel(
                                "Categorical covariables (e.g. your own clone IDs or new groupings)",
                                HTML(
                                  "Use column names to define the sequence ID (column name Sequence_ID) and the covariables (use your own nomenclature)<br>"
                              ),
                                fileInput(
                                  "categorical_covariable_file", "Upload", multiple = F,
                                  accept = c(".tab", ".csv", ".txt")
                              )
                            )
                          )
                        )
                      )
                    ),
                      tabPanel(
                        "4.Feature exploration", wellPanel(
                          HTML("<br></br>"),
                          HTML(
                            "Here you can see several Antibody features per region and per group (e.g. see where the NGly sites are, see the number of mutations over the whole sequence, etc.)"
                        ),
                          selectizeInput(
                            "plot_feature", label = NULL, choices = c(NULL),
                            selected = NULL
                        ),
                          uiOutput("Group_selection_for_feature"),
                          conditionalPanel(
                            "input.use_what.indexOf('Reconstructed germline') <= -1",
                            materialSwitch(
                              "show_reconstructed", "Include reconstructed sequences",
                              value = FALSE, width = NULL, status = "primary"
                          )
                        ),
                          materialSwitch(
                            "hide_points", "Remove non-selected data points", value = FALSE,
                            width = NULL, status = "primary"
                        ),
                          plotlyOutput("Violin_feature_plot")
                      )
                    ),
                      tabPanel(
                        "5.Clonal exploration", wellPanel(
                          HTML("<br></br>"),
                          HTML(
                            "Here you can group the sequences per different criteria, see how they behave on each sample and select them so you can mark their subclones in 3."
                        ),
                          hr(), fluidRow(
                            class = "clonal_row", column(
                              6, tabsetPanel(
                                tabPanel(
                                  "Clonal definition", HTML("<br>"),
                                  uiOutput("clonal_group_output"),
                                  conditionalPanel(
                                    "input.clonal_group != ''", numericInput(
                                      "dominance_threshold", "Dominance threshold (%)",
                                      0.5, min = 0, max = 100
                                  ),
                                    numericInput(
                                      "filter_clonal_group", "Do not plot clones with a frequency lower than (TO BE IMPLEMENTED)",
                                      0, min = 0, max = 100
                                  )
                                )
                              ),
                                tabPanel(
                                  "Calculate new clonal definition", HTML("<br>"),
                                  selectizeInput(
                                    "new_clonal_group", "Using base definition",
                                    choices = c("VCDR3J", "Reconstructed germline"),
                                    selected = c("VCDR3J"),
                                    multiple = F, options = list(plugins = list("remove_button"))
                                ),
                                  selectizeInput(
                                    "clonal_region_group", "using region(s)  (TO BE IMPLEMENTED)",
                                    choices = c(
                                      "Whole", "FWR1", "CDR1", "FWR2", "CDR2", "FWR3",
                                      "CDR3", "FWR4"
                                  ),
                                    selected = c("CDR3"),
                                    multiple = F, options = list(plugins = list("remove_button"))
                                ),
                                  numericInput(
                                    "identity_clonal_group", "with a % identity of  (TO BE IMPLEMENTED)",
                                    100, min = 50, max = 100
                                ),
                                  selectizeInput(
                                    "clonal_level_group", "at the level of  (TO BE IMPLEMENTED)",
                                    choices = c("AA", "NT"),
                                    selected = c("AA"),
                                    multiple = F, options = list(plugins = list("remove_button"))
                                ),
                                  materialSwitch(
                                    "calculate_shared_clones", "Look for shared clones between samples",
                                    value = FALSE, width = NULL, status = "primary"
                                ),
                                  actionButton("calculate_new_clone", "Calculate clonal definition")
                              )
                            ),
                              HTML("<br><br>")
                          )
                        ),
                          plotlyOutput("Violin_plot"),
                          tabsetPanel(
                            tabPanel(
                              "Shared clones", HTML("<br><br>"),
                              upsetjs::upsetjsOutput("upset_plot"),
                              HTML("<br>"),
                              hr(), HTML(
                                "Here you can see every sample included in the selected comparison in the UpSet plot, one by one"
                            ),
                              HTML("<br>"),
                              plotlyOutput("Comparison_plot")
                          ),
                            tabPanel("Clonal diversity"),
                            tabPanel("Clonal VDJ usage"),
                            tabPanel("Subclonal diversity and networks"),
                            tabPanel("Phylogenies")
                        )
                      )
                    ),
                      tabPanel("5.Cell exploration", wellPanel(HTML("<br></br>"))),
                      tabPanel(
                        "6.ML analysis", wellPanel(
                          tabsetPanel(
                            tabPanel(
                              "Training and validating a model", fluidRow(
                                HTML("<br></br>"),
                                HTML(
                                  "Notes for me: Once the results are out, make it also possible to compare sequences with the same number of mutations and see how they are classified by the method"
                              )
                            ),
                              fluidRow(
                                column(
                                  4, textInput("ML_name", "ML run name", "Pick a name for this ML run"),
                                  uiOutput("ML_group_selection"),
                                  uiOutput("ML_method_selection"),
                                  numericInput(
                                    "ML_use_percentage_selection_dataset", "Set the percentage of sequences from the subject samples that will be used for the ML model",
                                    100, min = 80, max = 100
                                ),
                                  numericInput(
                                    "ML_perc_training", "Set the percentage of the samples used for the training set",
                                    70, min = 30, max = 100
                                ),
                                  checkboxInput(
                                    "ML_equal_set_sizes", "Use the same size for all the groups",
                                    value = FALSE, width = NULL
                                ),
                                  checkboxInput(
                                    "ML_equal_VJs", "Use only the VJ combinations present in all the groups",
                                    value = FALSE, width = NULL
                                ),
                                  checkboxInput(
                                    "ML_only_one_subclone", "Use only one subclone from each clone",
                                    value = FALSE, width = NULL
                                ),
                                  uiOutput("ML_clone_def_to_use"),
                                  uiOutput("ML_assign_sample_to_training"),
                                  numericInput("seedML", label = "Numeric seed", value = 1234),
                                  verbatimTextOutput("onlynumbers"),
                                  actionButton(
                                    "update_ML_group_parameters", "Update training and validation sets",
                                    icon("forward"),
                                    style = "color: #fff; background-color: #c90076"
                                ),
                                  actionButton(
                                    "ML_run", "Run ML!", icon("forward"),
                                    style = "color: #fff; background-color: #CC2936"
                                )
                              ),
                                column(
                                  8, HTML("<h3>Size of groups in ML training and validation sets</h3>"),
                                  fluidRow(
                                    htmlOutput("ML_sunburst_dataset_ML") %>%
                                      withSpinner(type = 6, color = "#cf395c")
                                ),
                                  fluidRow(
                                    htmlOutput("ML_sunburst_dataset_grouping") %>%
                                      withSpinner(type = 6, color = "#cf395c")
                                )
                              )
                            )
                          ),
                            tabPanel(
                              "Further analysis - Artificial datasets, extending the predictions",
                              HTML("<br></br>"),
                              HTML(
                                "Here you can create artificial datasets (e.g. sequences with the same number of mutations as one of your groups or intermediate sequences from germline to repertoire), and see how the different ML kernels perform classifying them"
                            ),
                              hr(), HTML("Artificial dataset"),
                              hr(), HTML("ML exploration"),
                              uiOutput("ML_exploration_selection"),
                              plotlyOutput("ML_exploration")
                          )
                        )
                      )
                    ),
                      tabPanel(
                        "Help", wellPanel(
                          HTML(
                            "Welcome to AbSolution. The programme runs in different steps, each on a tab. On each step, this Help section will indicate what to do.<br><br>"
                        ),
                          uiOutput("HELP_output"),
                          HTML(
                            "<br>If at some point you need to reset AbSolution, you will find the appropiate button below."
                        ),
                          img(
                            src = "www/img/AbSolution_scheme.png", height = "100%",
                            width = "100%"
                        ),
                          HTML(
                            "<br>do a good help guide here and check how they did the documentation for immuneML, write down strenghts and points of improvement, let download the data behind a figure, maybe work in reproducibility</br>"
                        ),
                          actionButton(
                            "Reset_all", "Reset AbSolution", icon("rotate-left"),
                            style = "color: #fff; background-color: #18bc9c"
                        ),
                          materialSwitch(
                            "colorblind_mode", "Use colorblind-friendly palettes",
                            value = F, width = NULL, status = "primary"
                        )
                      )
                    )
                  )
                ),
                  tabItem(
                    tabName = "About", userBox(
                      title = userDescription(
                        title = "   Bioinformatics Laboratory", subtitle = "   Amsterdam UMC – location AMC",
                        image = "www/img/BioinformaticsLaboratory.png", type = 2
                    ),
                      status = "warning", width = 12, closable = F, footer = HTML(
                        "<p align='justify'>
                                  The Bioinformatics Laboratory was initiated in 1997 to strengthen biomedical research in the Academic Medical Center (AMC, Amsterdam). The laboratory is part of the department of Clinical Epidemiology, Biostatistics and Bioinformatics (CEBB) and is headed by Antoine van Kampen. The bioinformatics group members represent a broad range of bioinformatics, statistics, and (systems) biology expertise. They engage in research and support projects, and significantly contribute to various education and training programmes.
                                  <br><br>
                                  Our research focuses on the bioinformatics analysis of single-cell (RNAseq) data, the reconstruction of gene networks, and the computational modelling of biological systems. Research is strongly focussed on immunology.
                                  </p>"
                    ),
                      shiny::actionButton(
                        inputId = "github_visit", icon = icon("github"),
                        label = "GitHub", value = "Open pop-up", onclick = "window.open('https://github.com/EDS-Bioinformatics-Laboratory/','_blank')"
                    )
                  )
                )
              )
            )
        )
    )
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
    add_resource_path("www", app_sys("app/www"))
    tags$head(
        favicon(), bundle_resources(
            path = app_sys("app/www"),
            app_title = "AbSolution"
        )
    )
}


#' Prepare AbSolution logo
#'
#' This function is internally used to prepare the AbSolution logo
#'
#' @import dashboardthemes
#' @noRd
logo_absolution <- dashboardthemes::shinyDashboardLogoDIY(
    boldText = "Ab", mainText = "Solution", textSize = 16, badgeText = "BETA", badgeTextColor = "white",
    badgeTextSize = 2, badgeBackColor = "#cf395c", badgeBorderRadius = 3
)
