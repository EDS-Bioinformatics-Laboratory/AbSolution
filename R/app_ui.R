#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @import plotly
#' @import shinycssloaders
#' @import shinyFiles
#' @import shinymanager
#' @import shinythemes
#' @import upsetjs
#' @import utils
#' @import shinymeta
#' @import fresh
#' @noRd
app_ui <- function(request) {
    shiny::tagList(
        golem_add_external_resources(),
        bs4Dash::bs4DashPage(
          freshTheme =fresh::create_theme(
            bs4dash_vars(
              body_bg = "#2D4552",
              navbar_light_color = "#09BC8A",
              navbar_light_active_color = "#f4a261",
              navbar_light_hover_color = "#e76f51"
            ),

            bs4dash_yiq(
              contrasted_threshold = 10,
              text_dark = "#f8f9fa",
              text_light = "#2D4552"
            ),
            bs4dash_layout(
              main_bg = "#f8f9fa"
            ),
            bs4dash_sidebar_light(
              bg = "#2D4552",
              color = "#f8f9fa",
              hover_color = "#09BC8A",
              submenu_bg = "#2D4552",
              submenu_color = "#f8f9fa",
              submenu_hover_color = "#09BC8A"
            ),
            bs4dash_status(
              primary = "#2D4552", danger = "#e76f51", light = "#f8f9fa",
              secondary = "#09BC8A", warning="#f4a261"
            ),
            bs4dash_color(
              gray_900 = "#2D4552", white = "#f8f9fa"
            ),
            bs4dash_button(
              default_background_color = "#09BC8A",
              default_color="#f8f9fa"
            )
          ),
          header= bs4Dash::bs4DashNavbar(
              title=bs4Dash::bs4DashBrand(title="Solution", color = "primary",
                                 image = "www/favicon.png", opacity = 1),
              leftUi = tagList(
                # bs4Dash::dropdownMenu(
                #   badgeStatus = "info",
                #   type = "tasks",
                #   bs4Dash::taskItem(
                #     inputId = "triggerAction3",
                #     text = "My progress",
                #     color = "orange",
                #     value = 10
                #   )
                # )
              ),
              controlbarIcon = icon("circle-plus"),
              rightUi = bs4Dash::dropdownMenu(
                badgeStatus = "danger",
                type = "messages",
                bs4Dash::messageItem(
                  inputId = "triggerAction1",
                  message = "Have fun, AbSolutely! Do not forget to check the Help tab on each step",
                  from = "Rodrigo Garcia-Valiente",
                  image = "www/img/Rodrigo.jpg",
                  time = "today",
                  color = "lime"
                )
              )
            ),
            sidebar =bs4Dash::bs4DashSidebar(
              skin = "light",
              bs4Dash::bs4SidebarMenu(
                  id = "sidebar",
                  bs4Dash::bs4SidebarMenuItem("Home",
                           tabName = "home",
                           icon = bs4Dash::ionicon(name ="home"),
                           selected=T),
                  bs4Dash::bs4SidebarMenuItem(
                    "AbSolution",
                    icon =  bs4Dash::ionicon(name ="analytics"),
                    tabName = "pipeline"
                ),
                bs4Dash::bs4SidebarMenuItem(
                  "Export results",
                  icon =  bs4Dash::ionicon(name ="download"),
                  tabName = "export_results"
                ),
                bs4Dash::bs4SidebarMenuItem(
                  "Manual",
                  icon =  bs4Dash::ionicon(name ="map"),
                  tabName = "manual"
                ),
                bs4Dash::bs4SidebarMenuItem(
                    "About us", icon = bs4Dash::ionicon(name ="information-circle"),
                    tabName = "About"
                )
              )
            ),
          body =bs4Dash::bs4DashBody(
                shinyjs::useShinyjs(),
                tags$head(tags$link(rel = "shortcut icon",
                                    href = "www/img/favicon.ico")),
                tags$style(HTML(".fade:not(.show) {  opacity: 1; !important}"),
                           '.popover ul {padding: 15px; text-align: justify;}',
                           '.popover-content {text-align: justify;}'),
                # shinyDashboardThemes(theme = "poor_mans_flatly"),
                bs4Dash::bs4TabItems(
                  bs4Dash::tabItem(
                    tabName = "home", fluidRow(
                      bs4Dash::userBox(
                        title = NULL, type = 2, src = NULL, color = "#FFA500", width = 12,
                        boxToolSize = "lg", closable = F, collapsible = F, uiOutput("AbLogo"),
                        footer = HTML(
                          "<strong>AbSolution - Interactive Immune Repertoire Feature Analysis</strong>
                              <br><br>
                            <p align='justify'>
                              AbSolution is a user-friendly interactive platform for exploring B-cell and T-cell immune repertoires based on their sequence features.
                              It allows researchers to examine a wide range of physicochemical, structural, and mutation-related properties across receptor regions, clonotypes, samples, and conditions, helping to reveal patterns of immune selection and diversity.
                              By combining robust statistical analysis, intuitive visualizations, and built-in reproducibility through the ENCORE framework, AbSolution makes immune repertoire analysis accessible, transparent, and fully reproducible for both computational and experimental immunologists.
                              <br><br>
                              Use the menu on the left to navigate through the application:
                              <br><br>
                              <strong>AbSolution</strong> - Launch the interactive application to explore your data. A dedicated <em>Help</em> subtab provides personalized, step-by-step guidance for your analysis, and additional contextual help is available via the help button in the top-right corner.
                              <br>
                              <strong>Export results</strong> - Download your analysis and set up your project in a fully reproducible, ENCORE-compliant project workspace, including data, code, results, and documentation.
                              <br>
                              <strong>Manual</strong> - Access detailed documentation, tutorials, and practical examples.
                              <br>
                              <strong>About us</strong> - Learn more about the team behind AbSolution and our work in computational immunology and reproducible research.
                          <br><br><strong>Contact</strong>: <a href='mailto:r.garciavaliente@amsterdamumc.nl'>Rodrigo Garcia-Valiente</a> and <a href='mailto:a.h.vankampen@amsterdamumc.nl'> Antoine van Kampen</a>
                          </p>"
                      ),
                    HTML("<br><br>"),
                    textInput("username", "User name", "User", width="30%"),
                    textInput("usermail", "Contact email", "Email", width="30%")
                    )
                  )
                ),
                bs4Dash::tabItem(
                    tabName = "pipeline", h2("Analyze and explore your data"),

                  bs4Dash::tabBox(id = "menutabset", width=12, collapsible = F,

                      # 0.Project information  ######
                      tabPanel(
                        "0.Project information",
                          wellPanel(
                            fluidRow(

                            column(9,
                            h3("Upload your sample information table"),
                            bs4Dash::tooltip(
                              radioButtons(
                                inputId="TCRBCR_input_file",
                                label="The repertoire to analyze is:",
                                choiceNames = list("BCR", "TCR"),
                                choiceValues = list("BCR", "TCR"),
                                selected = "BCR"
                              ),
                              title = "TCR and BCR repertoires should be analyzed separatedly.",
                              placement = "bottom"
                            ) ,
                            bs4Dash::tooltip(radioButtons(
                              inputId="radiobutton_input_file",
                              label= "The input file type is:",
                              choiceNames = list("AIRR-Seq format"),
                              choiceValues = list("airrseq"),
                              selected = "airrseq"
                               ),
                              title = "Data should follow the AIRR-Seq format https://docs.airr-community.org/en/latest/datarep/rearrangements.html",
                              placement = "top"
                            ),
                            uiOutput("Conditional_Action_Filetype_description"),
                            uiOutput("Conditional_Action_Filetype_upload"),

                            HTML("<br>"),
                            DT::dataTableOutput("raw_sample_file_out"),
                            HTML("<br><br>"),
                            textAreaInput("user_0_comments",
                                          "User comments",
                                          "",
                                          width = "80%"),
                            HTML("<br><br><h3>Select work folder</h3>"),

                           shinyDirButton(
                              id="base_folder",
                              label="Select the folder where the data will be saved",
                              title = "Please select the work folder:", buttonType = "default",
                              class = NULL, icon = icon("folder", lib = "font-awesome"),
                              multiple = F
                              ),
                            HTML("<br><br><br><h3>Next step</h3>"),
                           HTML(
                             "If you have previous files, you don't have to redo all the steps. Process and upload your files where/when needed. But do provide the sample information.<br> <br>"
                           ),
                            uiOutput("Conditional_Action_Move_to_1"),
                            HTML("<br><br>"),
                            uiOutput("Conditional_Action_Move_to_Analysis")
                        ),
                          column(3,
                                 HTML("<br><br>"),
                                 uiOutput("Demo_Analysis")
                          )
                        )
                      )
                      ),
                      # 1.AIRR-Seq conversion  ######
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
                          bs4Dash::tooltip(
                            shinyWidgets::materialSwitch(
                              inputId="Dgene_reconstruct_airr",
                              label="Reconstruct CDR3 using the assigned D gene sequence instead of keeping that part of the CDR3 as is in the repertoire.",
                              value = T
                            ),
                            title= paste0(
                              "D genes, located in the CDR3, are of short length",
                              " and due to NPJ additions they are",
                              " hard to idenfity. You can decide to reconstruct it",
                              " using the inferred D germline and NPJs ",
                              " or keep the fragment of the sequence as they appear",
                              " in the repertoire"
                            ),
                            placement = "right"
                          ),
                          bs4Dash::tooltip(
                            shinyWidgets::materialSwitch(
                              inputId="FWR1partial_airr",
                              label= "Is there a partial FWR1 in the sequences?",
                              value = FALSE
                            ),
                            title= paste0(
                              "Depending on the primer location, sometimes the BCR",
                              "/TCR sequence is not fully sequenced. If this ",
                              "happens, AbSolution may be biased as not all ",
                              "sequences will have the same starting point. ",
                              "Activate this option if the FWR1 is not fully sequenced",
                              ". The FWR1 will be removed from the analysis"
                            ),
                            placement = "right"
                          ),
                          bs4Dash::tooltip(
                            shinyWidgets::materialSwitch(
                              inputId="FWR4partial_airr",
                              label="Is there a partial FWR4 in the sequences?",
                              value = FALSE
                            ),
                            title= paste0(
                              "Depending on the primer location, sometimes the BCR",
                              "/TCR sequence is not fully sequenced. If this ",
                              "happens, AbSolution may be biased as not all ",
                              "sequences will have the same ending point. ",
                              "Activate this option if the FWR4 is not fully sequenced.",
                              " The FWR4 will be removed from the analysis."
                            ),
                            placement = "right"
                          ),
                          bs4Dash::tooltip(
                            shinyWidgets::materialSwitch(
                              inputId="C_region_included_airr",
                              label="Is part or totatily of the pre-FWR1 or the C region included in the sequence?",
                              value = FALSE
                            ),
                            title=paste0(
                              "Depending on the primer location, sometimes the BCR",
                              "/TCR sequence pre-FWR1 or C region (or UMIs) are included in the",
                              "sequence. If this  happens, AbSolution may be biased",
                              " due to the sequence extra composition",
                              "Activate this option if the preFWR1/C gene/UMIs are present. ",
                              "The FWR4 will be removed from the analysis."
                            ),
                            placement="right"
                          ),

                          HTML("<br> Preprocessing preview (3 sequences) - Repertoire and germline extraction"),
                          DT::dataTableOutput("DT_example_parsed_sequence"),
                          HTML("<br>"),
                          textAreaInput("user_1_comments",
                                        "User comments",
                                        "",
                                        width = "80%"),
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
                      # 2.Sequence feature determination  ######
                      tabPanel(
                        "2.Sequence feature determination", wellPanel(
                          textAreaInput("user_2_comments",
                                        "User comments",
                                        "",
                                        width = "80%"),
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
                          HTML("<br>"),
                          uiOutput("make_dummy_2"),
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
                              "final_dummy",
                              "Produce the dummy dataset (click and wait until the table refreshes)",
                              icon("forward"),
                              style = "color: #fff; background-color: #18bc9c"
                            )
                          ),
                          uiOutput("Conditional_Action_Move_to_3")
                        )
                      ),
                      #Steps 1&2. Select your FBM and associated files  ######
                      tabPanel(
                        "Steps 1&2. Select your FBM and associated files", wellPanel(
                          shinyDirButton(
                            "FBM_folder", "Select the folder with your .rds and .bk files",
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
                              "final_dummy",
                              "Produce the dummy dataset (click and wait until the table refreshes)",
                              icon("forward"),
                              style = "color: #fff; background-color: #18bc9c"
                            )
                          ),
                          HTML("<br><br><br>"),
                          uiOutput("Conditional_Action_Move_to_Analysis_Real")
                        )
                      ),
                      # 3.Dataset exploration and variable selection ######
                      tabPanel(
                        "3.Dataset exploration and variable selection",
                        bs4Dash::tabBox(id = "3menu", width=12, collapsible = F,
                          tabPanel(
                            "Menu", fluidRow(column(
                              4, wellPanel(
                                HTML(
                                  "Explore the dataset using both plots. The dataset and variables that you will use for the feature selection process are represented in the Selection plot. <br><br>"
                                ),
                                h3("1. Select your sequences"),
                                bs4Dash::tooltip(
                                  selectizeInput(
                                    inputId="use_what",
                                    label = "Include sequences:",
                                    choices = c("Repertoire", "Reconstructed germline"),
                                    selected = c("Repertoire"),
                                    multiple = TRUE, options = list(plugins = list("remove_button"))
                                  ),
                                  title="Work with the sequenced sequences (repertoire) and/or their inferred germlines.",
                                  placement="right"
                                ),
                                bs4Dash::tooltip(
                                  selectizeInput(
                                    inputId="use_productive_or_not",
                                    label = "That are:",
                                    choices = c("Productive"),
                                    selected = c("Productive"),
                                    multiple = TRUE, options = list(plugins = list("remove_button"))
                                  ),
                                  title="Use productive and/or non-productive sequences. Currently only productive sequences are supported.",
                                  placement="right"
                                ),
                                uiOutput("Sample_selection"),
                                uiOutput("Chain_selection"),
                                uiOutput("Type_selection"),
                                bs4Dash::tooltip(
                                  sliderInput(
                                    inputId = "Rmut_filter",
                                    label = "Include only sequences with a number of replacement (R) mutations between",
                                    min = 0, max = 100, value = c(0, 100),
                                    step = 1
                                  ),
                                  title="Select minimum and maximum number of R mutations. This allows to e.g. filter naive sequences or outliers.",
                                  placement="right"
                                ),

                                uiOutput("VJ_selection"),
                                HTML("<br>"),
                                h3("2. Select your variables"),
                                HTML(
                                  "<i>NOTE: Variables with NAs/Infinite/Zero variance are automatically removed.</i><br><br>"
                                ),

                                bs4Dash::tooltip(
                                  selectizeInput(
                                    inputId="my_regions",
                                    label = "Select the region(s)",
                                    choices = c(
                                      "Whole", "FWR1", "CDR1", "FWR2", "CDR2", "FWR3",
                                      "CDR3", "FWR4"
                                    ),
                                    selected = c("CDR3"),
                                    multiple = TRUE, options = list(plugins = list("remove_button"))
                                  ),
                                  title= "Select which regions (and/or the whole Fab sequence) to include in the analysis.",
                                  placement="right"
                                ),
                                bs4Dash::tooltip(
                                  selectizeInput(
                                    inputId="my_var_elements",
                                    label = "Study features from",
                                    choices = c("NT", "AA"),
                                    selected = c("NT", "AA"),
                                    multiple = TRUE, options = list(plugins = list("remove_button"))
                                  ),
                                  title="Work with features calculated according to the AA and/or NT sequences",
                                  placement="right"
                                ),
                                bs4Dash::tooltip(
                                  selectizeInput(
                                    inputId="my_vars",
                                    label = "Select the variable(s) class(es)",
                                    choices = c(
                                      "Length", "Composition", "NGly sites", "Hot/Cold motifs",
                                      "Insertions",
                                      "Deletions", "Transitions and transversions",
                                      "Replacement and silent mutations", "Mutations from X to Y",
                                      "Peptide features"
                                    ),
                                    selected = c(
                                      "Length", "Composition", "NGly sites", "Hot/Cold motifs",
                                      "Peptide features"
                                    ),
                                    multiple = TRUE, options = list(plugins = list("remove_button"))
                                  ),
                                  title="Select or unselect groups of variables to include in the analysis.",
                                  placement="right"
                                ),
                                bs4Dash::tooltip(
                                  selectizeInput(
                                    inputId="my_vartypes",
                                    label = "Select the variable type(s)",
                                    choices = c("Baseline", "Germline diff"),
                                    selected = c("Baseline"),
                                    multiple = TRUE, options = list(plugins = list("remove_button"))
                                  ),
                                  title="Include variables measured at the sequence and/or its difference between the sequence and its germline.",
                                  placement="right"
                                ),

                                uiOutput("individual_variables"),
                                bs4Dash::tooltip(
                                  shinyWidgets::materialSwitch(
                                    inputId="use_UMAP",
                                    label = "Also represent the data in a UMAP (increases computational time)",
                                    value = FALSE, width = NULL, status = "primary"
                                  ),
                                  title="The UMAP will only be calculated if the data fits the RAM memory.",
                                  placement="right"
                                ),

                                HTML("<br>"),
                                h3("3. Define the groups to study"),


                                bs4Dash::tooltip(
                                  shinyWidgets::materialSwitch(
                                    inputId="work_as_categories",
                                    label = "I want to compare between groups",
                                    value = FALSE, width = NULL, status = "primary"
                                  ),
                                  title="Do you want to work only with certain groups of interest?",
                                  placement="right"
                                ),
                                uiOutput("Group_selection"),
                                uiOutput("Group_comparison"),
                                hr(),
                                div(
                                  style = "margin-left: 20px;",
                                  shinyWidgets::materialSwitch(
                                    inputId="use_sharedVDJ",
                                    label = "Include only sequences with VJ genes present in all the groups",
                                    value = FALSE, width = NULL, status = "primary"
                                  )
                                ),
                                div(
                                  style = "margin-left: 40px;",
                                  bs4Dash::tooltip(
                                    shinyWidgets::materialSwitch(
                                      inputId="VDJ_normalized_per_size",
                                      label = "Use the same number of sequences for each VJ combination on each group",
                                      value = T, width = NULL, status = "primary"
                                    ),
                                    title="This normalizes the feature selection and avoids overrepresentation of one group",
                                    placement="right"
                                  )
                                ),
                                div(
                                  style = "margin-left: 60px;",
                                  bs4Dash::tooltip(
                                    shinyWidgets::materialSwitch(
                                      inputId="VDJ_maximize_clones",
                                      label = "Maximize the number of clones for each VJ combination on each group",
                                      value = T, width = NULL, status = "primary"
                                    ),
                                    title="This normalizes the feature selection and avoids overrepresentation of major clones",
                                    placement="right"
                                  ),
                                  uiOutput("my_clone_def")
                                ),

                                HTML("<br>"),
                                div(
                                  style = "margin-left: 60px;",
                                  shinyWidgets::materialSwitch(
                                    inputId="VDJ_normalized_per_sample",
                                    label = "And apply these rules also to each individual sample",
                                    value = F, width = NULL, status = "primary"
                                  )
                                ),
                                div(
                                  style = "margin-left: 40px;",
                                  uiOutput("VDJ_subsetting")
                                ),
                                # selectizeInput(
                                #   inputId="my_clone_def",
                                #   label = "Select the clone definition",
                                #   choices = c(NULL),
                                #   selected = c(NULL),
                                #   multiple = F, options = list(plugins = list("remove_button"))
                                # ),
                                shinyWidgets::materialSwitch(
                                  inputId="use_univlog",
                                  label = "Include only features that are correlated to the groups (0/1)" ,
                                  value = FALSE, width = NULL, status = "primary"
                                ),
                                selectizeInput(
                                  inputId="pval_type",
                                  label = "Select the P-value to use" ,
                                  choices = c("Raw p-value",
                                              "Corrected by Bonferroni method",
                                              "Corrected by Benjamini & Hochberg method",
                                              "Corrected by Benjamini & Yekutieli method"),
                                  selected = c("Corrected by Benjamini & Yekutieli method"),
                                  multiple = F, options = list(plugins = list("remove_button"))
                                ),
                                numericInput(
                                  inputId="pval_cutoff",
                                  label = "P-value cut-off (<=):", 0.05,
                                  min = 0, max = 1
                                ),
                                numericInput(
                                  inputId="estimate_cutoff",
                                  label = "Estimate cut-off (absolute value, >=):",
                                  1, min = 0
                                ),
                                selectizeInput(
                                  inputId="number_selected_vars",
                                  label = "How many selected features to use?",
                                  choices = c(
                                    "All", "Top 10 (according to estimate)",
                                    "Top 50 (according to estimate)"
                                  ),
                                  selected = c("All"),
                                  multiple = F, options = list(plugins = list("remove_button"))
                                ),
                                h3("4. Update your changes"),
                                actionButton(
                                  "Update_selection", "Update Selection set, used for all the analyses",
                                  icon("forward"),
                                  style = "color: #fff; background-color: #18bc9c"
                                ),
                                actionButton(
                                  "Update_exploration", "Update Exploration set, to compare with Selection",
                                  icon("forward"),
                                  style = "color: #fff; background-color: #CC2936"
                                )
                              )
                            ),
                            column(
                              8, HTML("<br>"),
                              HTML("Use the lasso tool or box select to get the sequence(s) information.<br><br>"),
                              fluidRow(
                                column(
                                  4, HTML("<h2>Selection plot</h2>"),
                                  fluidRow(
                                    column(
                                      6, selectInput(
                                        inputId="Selection_plot_type",
                                        label = "Dimensionality Reduction:",
                                        c(PCA = "PCA")
                                      ),
                                      selectizeInput(
                                        inputId="Selection_plot_type_dim",
                                        label = "Select the (1 to 5) dimensions to plot",
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
                                column(8, plotlyOutput("Selection_plot")%>%
                                         withSpinner(type = 6, color = "#2D4552"),
                                       uiOutput("Selection_plot_error"))
                              ),
                              fluidRow(
                                column(
                                  4, HTML("<br><br><h2>Exploration plot</h2>"),
                                  fluidRow(
                                    column(
                                      6, selectInput(
                                        inputId="Exploration_plot_type",
                                        label = "Dimensionality Reduction:",
                                        c(PCA = "PCA")
                                      ),
                                      selectizeInput(
                                        inputId="Exploration_plot_type_dim",
                                        label ="Select the (1 to 5) dimensions to plot",
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
                                column(8,
                                       HTML("<br><br>"),
                                       plotlyOutput("Exploration_plot")%>%
                                         withSpinner(type = 6, color = "#2D4552"),
                                       uiOutput("Exploration_plot_error"))
                              ),
                              fluidRow(
                                column(6, bs4Dash::infoBoxOutput("RAM_S_Memory_Box", width = 12)),
                                column(6, bs4Dash::infoBoxOutput("RAM_E_Memory_Box", width = 12))
                              ),
                              textAreaInput("user_3_comments",
                                            "User comments",
                                            "",
                                            width = "80%"),
                              DT::dataTableOutput("Ab_table")
                            ))
                          ),
                          tabPanel(
                            "About the variables", uiOutput("about_vars"),
                            bs4Dash::tabBox(id = "varmenu", width=12, collapsible = F,
                              tabPanel(
                                "Used variables", fluidRow(
                                  column(
                                    10, HTML(
                                      "<br>Here you can see which variables are represented on each set.<br><br>The plots can take a few minutes to load!<br><br>"
                                    ),
                                    htmlOutput("used_variables") %>%
                                      withSpinner(type = 6, color = "#2D4552")
                                  )
                                )
                              ),
                              tabPanel(
                                "Variable correlation", HTML(
                                  "<br>Here you can see which variables in the selection set are correlated >=0.8 for their non-zero values. Normalized and raw variables are not considered as correlated in the plot.<br><br>The plots can take a few minutes to load!<br><br>"
                                ),
                                htmlOutput("correlated_variables") %>%
                                  withSpinner(type = 6, color = "#2D4552")
                              ),
                              tabPanel(
                                "PCA loadings", HTML(
                                  "<br>Here you can see the PCA loadings in both PCAs.<br><br>"
                                ),
                                htmlOutput("PCA_loadings") %>%
                                  withSpinner(type = 6, color = "#2D4552")
                              )
                            )
                          ),
                          tabPanel(
                            "About the samples", HTML(
                              "<br>Here you can see information about the samples that are present on the selection set.<br><br>Remember you can have included (reconstructed) germline sequences!"
                            ),
                            htmlOutput("about_samples") %>%
                              withSpinner(type = 6, color = "#2D4552")
                          ),
                          tabPanel(
                            "Add covariables", HTML("<br>Here you can add covariables to the selection dataset!<br><br>"),
                            uiOutput("about_meh"),
                            bs4Dash::tabBox(id = "covarmenu", width=12, collapsible = F,
                              # tabPanel(
                              #   "Numerical covariables (e.g. gene expression)", HTML(
                              #     "Use column names to define the sequence ID (column name Sequence_ID) and the covariables (use your own nomenclature)<br>"
                              #   ),
                              #   fileInput(
                              #     "numerical_covariable_file", "Upload", multiple = F,
                              #     accept = c(".tab", ".csv", ".txt")
                              #   ),
                              #   uiOutput("numerical_covariables_plot")
                              # ),
                              tabPanel(
                                "Categorical covariables (e.g. your own clone IDs or new groupings).",
                                HTML(
                                  "Use column names to define the sequence ID (column name Sequence_ID) and the covariables (use your own nomenclature). Use .csv format.<br>"
                                ),
                                fileInput(
                                  "categorical_covariable_file", "Upload", multiple = F,
                                  accept = c(".csv")
                                ),
                                uiOutput("results_categorical")
                              )
                            )
                          )
                        )
                      ),
                      # 4.Feature exploration ######
                      tabPanel(
                        "4.Feature exploration", wellPanel(
                          HTML("<br>"),
                          HTML(
                            "Here you can inspect several features of BCR and TCR Fab regions per region and per group (e.g. N-glycosylation sites in the FWR3, total number of mutations along the sequence, etc.)."
                          ),
                          fluidRow(
                            column(6,
                                   HTML("<br>"),
                                   selectizeInput(
                                     inputId="plot_feature", label = NULL, choices = c(NULL),
                                     selected = NULL
                                   ),
                                   shinyWidgets::materialSwitch(
                                     inputId="show_reconstructed",
                                     label = "Include reconstructed sequences",
                                     value = FALSE, width = NULL, status = "primary"
                                   ),
                                   shinyWidgets::materialSwitch(
                                     inputId="compare_opposites",
                                     label = "Do a t-test between repertoire and
                                     reconstructed sequences if Cohen's d > 0.2",
                                     value = FALSE, width = NULL, status = "primary"
                                   ),
                                   shinyWidgets::materialSwitch(
                                     inputId="really_hide_points",
                                     label = "Do not plot points (ergo, no interaction)",
                                     value = TRUE,
                                     width = NULL, status = "primary"
                                   ),
                                   shinyWidgets::materialSwitch(
                                     inputId="hide_points",
                                     label = "Hide non-selected data points",
                                     value = TRUE,
                                     width = NULL, status = "primary"
                                   ),
                                   column(6,
                                          colourpicker::colourInput("primary_color", "Select primary colour", "#09BC8A"),
                                          uiOutput("Group_selection_for_feature"))

                            ),
                            column(6,
                                   textAreaInput("user_4_comments",
                                                 "User comments",
                                                 "",
                                                 width = "80%"))
                          ),

                          plotlyOutput("Violin_feature_plot")%>%
                            withSpinner(type = 6, color = "#2D4552"),

                        )
                      ),
                      # 5.Clonal exploration ######
                      tabPanel(
                        "5.Clonal exploration", wellPanel(
                          HTML("<br></br>"),
                          HTML(
                            "Here you can group the sequences per different criteria, see how they behave on each sample and select them so you can mark their subclones in 3."
                          ),
                          hr(), fluidRow(
                            class = "clonal_row", column(
                              6,
                              bs4Dash::tabBox(id = "clonemenu1", width=12, collapsible = F,
                                tabPanel(
                                  "Clonal definition", HTML("<br>"),
                                  uiOutput("clonal_group_output"),
                                  conditionalPanel(
                                    "input.clonal_group != ''",
                                    bs4Dash::tooltip(
                                      numericInput(
                                        inputId="dominance_threshold",
                                        label = "Dominance threshold (%)",
                                        0.5, min = 0, max = 100
                                      ),
                                      title="Relative frequency cut-off to determinate which clones are (non)-dominant.",
                                      placement="right"
                                    ),
                                    shinyWidgets::materialSwitch(
                                      inputId="really_hide_points_clones",
                                      label = "Do not plot points (ergo, no interaction)",
                                      value = TRUE,
                                      width = NULL, status = "primary"
                                    ),
                                    bs4Dash::tooltip(
                                      numericInput(
                                        inputId="filter_clonal_group",
                                        label = "Do not plot clones with a frequency lower than (TO BE IMPLEMENTED)",
                                        0, min = 0, max = 100
                                      ),
                                      title="Hide low-frequency clones to reduce loading times.",
                                      placement="right"
                                    )
                                  )
                                ),
                                tabPanel(
                                  "Calculate new clonal definition", HTML("<br>"),
                                  bs4Dash::tooltip(
                                    selectizeInput(
                                      inputId="new_clonal_group",
                                      label = "Using base definition",
                                      choices = c("VCDR3J"), #, "Reconstructed germline"
                                      selected = c("VCDR3J"),
                                      multiple = F, options = list(plugins = list("remove_button"))
                                    ),
                                    title="Definition to calculate clones. Either by same V-J combination + CDR3 length + similarity or by same VJ combination + same reconstructed germline.",
                                    placement="right"
                                  ),

                                  numericInput(
                                    inputId="identity_clonal_group",
                                    label ="with a % identity of ",
                                    100, min = 100, max = 100
                                  ),

                                  selectizeInput(
                                    inputId="clonal_region_group",
                                    label = "in the region(s)",
                                    choices = c(
                                      # "Whole", "FWR1", "CDR1", "FWR2", "CDR2", "FWR3", "FWR4",
                                      "CDR3"
                                    ),
                                    selected = c("CDR3"),
                                    multiple = F, options = list(plugins = list("remove_button"))
                                  ),
                                  selectizeInput(
                                    inputId="clonal_level_group",
                                    label ="at the sequence level of ",
                                    choices = c("AA", "NT"),
                                    selected = c("AA"),
                                    multiple = F, options = list(plugins = list("remove_button"))
                                  ),
                                  shinyWidgets::materialSwitch(
                                    inputId="calculate_shared_clones",
                                    label ="Contemplate also the possibility of shared clones between samples",
                                    value = FALSE, width = NULL, status = "primary"
                                  ),
                                  actionButton("calculate_new_clone", "Calculate clonal definition")
                                )
                              ),
                              HTML("<br><br>")
                            ),
                            column(6,
                                   textAreaInput("user_5_comments",
                                                 "User comments",
                                                 "",
                                                 width = "80%"))
                          ),
                          plotlyOutput("Violin_plot")%>%
                            withSpinner(type = 6, color = "#2D4552"),
                          bs4Dash::tabBox(id = "clonemenu2", width=12, collapsible = F,
                            tabPanel(
                              "Shared clones", HTML("<br><br>"),
                              upsetjs::upsetjsOutput("upset_plot")%>%
                                withSpinner(type = 6, color = "#2D4552"),
                              HTML("<br>"),
                              hr(), HTML(
                                "In the scatterplot(s) below, you can observe the shared and non-shared clones between the samples involved in the selected intersection of the above UpSet plot. Each comparison is made on a one-to-one basis between the samples."
                              ),
                              uiOutput("coloring_scatterplot_clones"),
                              HTML("<br>"),
                              plotlyOutput("Comparison_plot")%>%
                                withSpinner(type = 6, color = "#2D4552")
                            ),
                            tabPanel("Selected clones",
                                     DT::dataTableOutput("Clonal_table")
                            )
                            # tabPanel("Clonal diversity"),
                            # tabPanel("Clonal VDJ usage"),
                            # tabPanel("Subclonal diversity and networks"),
                            # tabPanel("Phylogenies")
                          )
                        )
                      ),
                      # Help ######
                      tabPanel(
                        "Help", wellPanel(

                          HTML(
                            "Welcome to AbSolution. The programme runs in different steps, each on a tab. For each step, this Help section will explain  what to do. <br><br>"
                          ),
                          # verbatimTextOutput("code"),
                          # h3("About the current step"),
                          uiOutput("HELP_output")

                          # img(
                          #   src = "www/img/AbSolution_scheme.png", height = "100%",
                          #   width = "100%"
                          # )



                        )
                      )
                  )
                  ),
                bs4Dash::tabItem(
                  tabName = "export_results", fluidRow(
                    tabPanel(
                      "Export results", wellPanel(
                        HTML("<br></br>"),
                        HTML(
                          '<p align="justify">
                                Set up a ready-to-use ENCORE-standard project <strong>workspace</strong> that organizes your data, code, results, and documentation in a single, well-structured collection-making your analysis transparent, well documented, and easy to reproduce. Place your input data in the <code>/Data</code> folder, use a dedicated subfolder within <code>/Processing</code> as your working directory, and, after completing a full analysis with AbSolution, export all results to <code>/Processing</code> to keep your project fully traceable and organized.
                          </p>'
                        ),
                        hr(), fluidRow(
                          class = "clonal_row", column(
                            6,
                            bs4Dash::tabBox(id = "clonemenu1", width=12, collapsible = F,
                                            tabPanel(
                                              "Set up a reproducible ENCORE Project",
                                              HTML("If you want to <strong>create a complete <a href='https://www.nature.com/articles/s41467-024-52446-8'>ENCORE</a> project compendium</strong> -a self-contained folder that brings together all parts of your research (data, code, results, and notes) in a reproducible structure-, you can download the latest template here to get started."),
                                              HTML("<br><br>"),
                                              fluidRow(
                                                column(4,
                                                       textInput("sFSS_name", label=NULL,"Project_name", width="80%"),
                                                ),
                                                column(8,
                                                       shinyDirButton(
                                                         id="sFSS_folder",
                                                         label="Select the folder where the ENCORE sFSS will be set-up",
                                                         title = "Please select the work folder:", buttonType = "default",
                                                         class = NULL, icon = icon("folder", lib = "font-awesome"),
                                                         multiple = F
                                                       )
                                                )
                                              ),


                                              HTML("<br>"),
                                              fluidRow(
                                                column(12, align = "center",
                                                       uiOutput("Conditional_Action_downloadsFSS")
                                                )
                                              )


                                            ),
                                            tabPanel(
                                              "Export analysis",
                                              HTML("    <p ><strong>After running a full analysis</strong> with the app, you can download a .zip file that follows the <a href='https://www.nature.com/articles/s41467-024-52446-8'>ENCORE</a> structure. This .zip file contains the following subfolders:</p>
                                                <ul>
                                                  <li><strong>/0_SoftwareEnvironment/R/</strong>: Contains the AbSolution package in the version used for the analysis, a Docker file, and the renv information to reproduce everything.</li>
                                                  <li><strong>/Data/Dataset/</strong>: Contains the dataset used for the analysis.<ul>
                                                  <li><strong>/Raw/</strong>: Contains the sample input files.</li>
                                                  <li><strong>/Meta/</strong>: Contains a file with information about the samples.</li>
                                                  <li><strong>/Processed/</strong>: Contains the files produced during the analysis with AbSolution (sequence and feature information). These files are used for analysis and data exploration.</li></ul></li>
                                                  <li><strong>/Notebook/</strong>: Contains the .Rmd file to generate the figures. This .Rmd file is generated directly from the app's code, capturing the domain logic using shinymeta. Additionally, once exported, it uses relative paths instead of absolute paths, making it easy to knit successfully. The chunks <code>parse_input</code> and <code>feature_calculation</code> are not evaluated unless the user sets <code>eval=TRUE</code>, as they are used to produce the Processed files. For time efficiency, this part is skipped when producing the .zip.</li>
                                                  <li><strong>/Results/</strong>: Contains the HTML file produced from the .Rmd file with the key plots from AbSolution. The HTML file includes all the system information and package versioning details.</li>
                                              </ul>"),
                                              textInput("analysis_name", "Name of the exported .zip", paste("Analysis", Sys.time(), sep="_"), width="30%"),
                                              shinyWidgets::materialSwitch(
                                                "include_data_report", "Include Data subfolder",
                                                value = F, width = "100%", status = "primary"
                                              ),
                                              shinyWidgets::materialSwitch(
                                                "include_docker_renv", "Include 0_SoftwareEnvironment subfolder",
                                                value = F, width = "100%", status = "primary"
                                              ),
                                              fluidRow(
                                                column(12, align = "center",
                                                       uiOutput("export_cond")%>%
                                                         withSpinner(type = 6, color = "#2D4552")
                                                )
                                              )

                                            )

                                            )
                            ),
                            HTML("<br><br>")
                          )
                        ),

                      )
                    )
                ),
                bs4Dash::tabItem(
                  tabName = "manual", fluidRow(
                    fluidRow(
                      tags$p(
                        "This is the AbSolution manual."
                      ),
                      tags$p(
                        "For personalized, step-by-step guidance, refer to the ",
                        tags$strong("Help"),
                        " subtab in the AbSolution tab."
                      ),
                      tags$p(
                        "Additionally, you can activate the ",
                        tags$strong("tooltip mode"),
                        " by clicking the ? at the top right corner to see explanations for each option when you hover over them."
                      ),
                      tags$iframe(
                        style="height:900px; width:100%;",
                        src="www/doc/AbSolution.pdf"
                      )
                    )
                  )
                ),
                bs4Dash::tabItem(
                  tabName = "About", bs4Dash::userBox(
                    title = bs4Dash::userDescription(
                      title = "   Bioinformatics Laboratory", subtitle = "  Amsterdam UMC",
                      image = "www/img/BioinformaticsLaboratory.png", type = 2
                    ),
                    status = "orange", width = 12, closable = F, footer = HTML(
                      "<p align='justify'>
                            The Bioinformatics Laboratory (Biolab) at Amsterdam UMC is a leading research group built on two key pillars: the development and application of advanced computational immunology methods, and an unwavering commitment to reproducible research. The team specializes in techniques such as transcriptomic data analysis (including single-cell RNA sequencing) and gene regulatory network modeling to unravel complex immune system dynamics. Equally central is Biolab's focus on reproducible science: our projects follow the ENCORE framework, an in-house infrastructure that supports transparent, well-structured, and reusable computational analyses. By combining innovative immunological data science with rigorous reproducibility, the lab accelerates robust insights into disease and enables others to confidently build on our work. This dual focus has also fostered close collaborations with immunology experts, including previous work with Niek de Vries and Jeroen Guikema, underscoring the translational impact of Biolab's research.
                      </p>"
                    ),
                    shiny::actionButton(
                      inputId = "github_visit", icon = icon("github"),
                      label = "GitHub", value = "Open pop-up", onclick = "window.open('https://github.com/EDS-Bioinformatics-Laboratory/','_blank')"
                    ),
                    shiny::actionButton(
                      inputId = "web_visit", icon = icon("laptop"),
                      label = "Website", value = "Open pop-up", onclick = "window.open('https://bioinformaticslaboratory.eu/','_blank')"
                    )
                  ),
                  bs4Dash::userBox(
                    title = bs4Dash::userDescription(
                      title = "   AbSolution", subtitle = "  A practical implementation to improve reproducibility and transparency of computational research",
                      image = "www/img/AbSolution.png", type = 2
                    ),
                    status = "gray-dark", width = 12, closable = F, footer = HTML(
                      "<p align='justify'>
                              AbSolution is a user-friendly interactive platform for exploring B-cell and T-cell immune repertoires based on their sequence features.
                              It allows researchers to examine a wide range of physicochemical, structural, and mutation-related properties across receptor regions, clonotypes, samples, and conditions, helping to reveal patterns of immune selection and diversity.
                              By combining robust statistical analysis, intuitive visualizations, and built-in reproducibility through the ENCORE framework, AbSolution makes immune repertoire analysis accessible, transparent, and fully reproducible for both computational and experimental immunologists.
                      </p>"
                    ),
                    shiny::actionButton(
                      inputId = "github_ENCORE_visit", icon = icon("github"),
                      label = "GitHub", value = "Open pop-up", onclick = "window.open('https://github.com/EDS-Bioinformatics-Laboratory/AbSolution,'_blank')"
                    ),
                    shiny::actionButton(
                      inputId = "web_ENCORE_visit", icon = icon("newspaper"),
                      label = "Article 1", value = "Open pop-up", onclick = "window.open('https://media0.giphy.com/media/v1.Y2lkPTc5MGI3NjExMzVhc2lmN2NkY2V3d2JoYm9wOWg0c3JtOGdwNjN2aWpsNWFxaGN4NSZlcD12MV9pbnRlcm5hbF9naWZfYnlfaWQmY3Q9Zw/YQAuKJ7wf68qBHPw6Y/giphy.gif','_blank')"
                    ),
                    shiny::actionButton(
                      inputId = "web_ENCORE_visit", icon = icon("newspaper"),
                      label = "Article 2", value = "Open pop-up", onclick = "window.open('https://media0.giphy.com/media/v1.Y2lkPTc5MGI3NjExMzVhc2lmN2NkY2V3d2JoYm9wOWg0c3JtOGdwNjN2aWpsNWFxaGN4NSZlcD12MV9pbnRlcm5hbF9naWZfYnlfaWQmY3Q9Zw/YQAuKJ7wf68qBHPw6Y/giphy.gif','_blank')"
                    )
                  ),
                  bs4Dash::userBox(
                    title = bs4Dash::userDescription(
                      title = "   ENCORE", subtitle = "  A practical implementation to improve reproducibility and transparency of computational research",
                      image = "www/img/ENCORE.png", type = 2
                    ),
                    status = "gray-dark", width = 12, closable = F, footer = HTML(
                      "<p align='justify'>
                            <small>
                              van Kampen et al., Nature Communications (2024)
                              <a href='https://www.nature.com/articles/s41467-024-52446-8' target='_blank'>https://doi.org/10.1038/s41467-024-52446-8</a>
                            </small> <br><br>
                            ENCORE (ENhancing COmputational REproducibility) is a practical, group-developed framework that elevates transparency and reproducibility in computational research by guiding scientists on how to structure, document, and share all components of a project - from data and code to results and documentation - within a standardized, self-contained project compendium. Designed to be agnostic of programming language, data type, or infrastructure, ENCORE integrates best practices into everyday workflows, leverages modern version control (e.g., GitHub), and includes tools such as pre-defined documentation templates and an HTML-based navigator to make computational research clearer, more reliable, and easier to reuse or build upon. By embedding reproducibility at the core of every project, ENCORE helps teams produce research that is not only trustworthy but also ready for open science, collaboration, and long-term impact.
                      </p>"
                    ),
                    shiny::actionButton(
                      inputId = "github_ENCORE_visit", icon = icon("github"),
                      label = "GitHub", value = "Open pop-up", onclick = "window.open('https://github.com/EDS-Bioinformatics-Laboratory/ENCORE,'_blank')"
                    ),
                    shiny::actionButton(
                      inputId = "web_ENCORE_visit", icon = icon("newspaper"),
                      label = "Article", value = "Open pop-up", onclick = "window.open('https://www.nature.com/articles/s41467-024-52446-8','_blank')"
                    )
                  )
                )
                )
          ),
          controlbar = bs4Dash::dashboardControlbar(
            id = "controlbar",
            collapsed = T,
            overlay = F,
            skin="dark",
            width=300,
            bs4Dash::controlbarMenu(
              id = "controlbarmenu",
              bs4Dash::controlbarItem(
                title = "Main parameters",
                column(
                  width = 12,
                  align = "left",
                  numericInput(inputId="seed",
                               label = "Set random seed",
                               value = 1234,
                               min=1,
                               max=9999,
                               width="100%"),
                  verbatimTextOutput("onlynumbers"),
                  shinyWidgets::materialSwitch(
                    "colorblind_mode", "Use colorblind-friendly palettes",
                    value = F, width = NULL, status = "primary"
                  ),
                  selectInput(inputId = "img_type",
                              label = "Select img type for exporting individual plots",
                              choices = c("PNG" = "png",
                                          "SVG" = "svg"),
                              selected = "png",
                              width = "100%"),
                  numericInput(inputId="pixels_width",
                               label = "Png images width (px)",
                               value = 1400,
                               min=300,
                               max=30000,
                               width="100%"),
                  verbatimTextOutput("onlynumberswidth"),
                  numericInput(inputId="pixels_height",
                               label = "Png images height (px)",
                               value = 1000,
                               min=300,
                               max=30000,
                               width="100%"),
                  verbatimTextOutput("onlynumbersheight"),
                  numericInput(inputId="labelscale",
                               label = "Svg/png width and height size multiplier",
                               value = 1,
                               min=1,
                               max=10,
                               width="100%"),
                  verbatimTextOutput("onlynumberslabelscale"),
                  # actionButton(
                  #   "Reset_all", "Reload AbSolution", icon("rotate-left"),
                  #   style = "color: #fff; background-color: #18bc9c"
                  # )
                )
              )
            )

          ),
          footer=bs4Dash::dashboardFooter(
            left = a(
              href = "https://www.linkedin.com/in/rodrigogarciavaliente/",
              target = "_blank", "@RGarciaValiente"
            ),
            right = "2026, version 1.0.0"
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
    # add_resource_path("doc", app_sys("app/doc"))
    #https://github.com/ThinkR-open/golem/issues/297
    # add_resource_path(
    #   "sbs", system.file("www", package = "shinyBS")
    # )

    options(shiny.maxRequestSize=3000*1024^2)

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
