#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @import dplyr
#' @import plotly
#' @import reactable
#' @import bigstatsr
#' @import shinymanager
#' @import sortable
#' @import tools
#' @import shinymeta
#' @noRd
app_server <- function(input, output, session) {

    volumes <- shinyFiles::getVolumes()
    verbose <- golem::get_golem_options("verbose")

    server_logic <- function()  {
      # Page 0.Project information ######

      output$Conditional_Action_Filetype_description <- renderUI(
        {
          HTML(
            "Upload the table with your sample information in .txt or .tab format, delimited by tab. It must have at least 6 fields (Filename, Sample, Patient, Group and Subgroup).<br> <br>"
          )

        }
      )

      output$Conditional_Action_Filetype_upload <- renderUI(
        {
          bs4Dash::tooltip(
            fileInput("raw_sample_file",
                      "Upload Sample information",
                      multiple = F,
                      accept = c(".tab", ".txt")),
            title = paste0("Upload the table with your sample information",
                           " in .txt or.tab format, delimited by tab. It must have at least",
                           " 6 fields (Filename, Sample, Patient, Group and Subgroup)."),
            placement = "top")
        }
      )

      sample_info_react <- reactiveValues(table = NULL, summary_status = F, test_status=F, step_3=F)


      o_rws <- metaObserve2({
        shiny::req(input$raw_sample_file$datapath)
        isolate(metaExpr({
          sample_info_react$table <- read.table(..(input$raw_sample_file$datapath), header = T, sep = "\t")
        }))
      })

      # observeEvent(
      #     input$raw_sample_file$datapath, {
      #         sample_info_react$table <- read.table(input$raw_sample_file$datapath, header = T, sep = "\t")
      #     }
      # )

      o_Mt1a <- metaObserve2({
        shiny::req(input$Move_to_1_airr)
        isolate(metaExpr({
          if(sample_info_react$test_status) {
            write.table(alakazam::Example10x,
                        file = file.path(shinyFiles::parseDirPath(volumes, input$base_folder), "Example10x_alakazam.tsv"),
                        quote=F,sep="\t",row.names =F)
            shinyWidgets::updateMaterialSwitch(session= session,
                                               inputId="C_region_included_airr",
                                               value = TRUE)
          }
        }))
      })
      # observeEvent(input$Move_to_1_airr, {
      #   if(sample_info_react$test_status) {
      #     write.table(alakazam::Example10x,
      #                 file = file.path(parseDirPath(volumes, input$base_folder), "Example10x_alakazam.tsv"),
      #                 quote=F,sep="\t",row.names =F)
      #     updateMaterialSwitch(session= session,
      #                          inputId="C_region_included_airr",
      #                          value = TRUE)
      #   }
      # })
      output$Conditional_Action_Move_to_1 <- renderUI(
        {
          if ((sample_info_react$test_status == F && is.null(input$raw_sample_file)) ||
              all(is.numeric(unlist(input$base_folder)))) {
            HTML(
              "IMPORTANT! <i>Fill first the required fields. Only after that you will be able to proceed with the analysis.<i>"
            )
          } else {
            sample_info <- sample_info_react$table

            actionButton(
              "Move_to_1_airr", "Proceed to the next step: Pre-process your data",
              icon("angle-right")
            )
          }
        }
      )

      output$Conditional_Action_Move_to_Analysis <- renderUI(
        {
          if (is.null(input$raw_sample_file) ||
              all(is.numeric(unlist(input$base_folder)))) {
          } else {
            sample_info <- sample_info_react$table
            actionButton(
              "Move_to_analysis", "Skip next steps: You already have AbSolution-feature-calculated files.",
              icon("forward"),
              style = "color: #fff; background-color: #18bc9c"
            )
          }
        }
      )

      output$Demo_Analysis <- renderUI(
        {
          if (!is.null(input$raw_sample_file) ||
              !all(is.numeric(unlist(input$base_folder)))) {
          } else {
            bs4Dash::tooltip(
              actionButton(
                "lets_demo_analysis", "Instead, load a test dataset.",
                icon("mug-saucer"),
                style = "color: #fff; background-color: #18bc9c"
              ),
              # title = paste0(
              #   "Load a ",
              #   a("small example database subset",
              #     href = "https://alakazam.readthedocs.io/en/stable/topics/Example10x/",
              #     target="_blank"),
              #   " from ",
              #   a("10x Genomics",
              #     href = "https://support.10xgenomics.com/single-cell-vdj/datasets/2.2.0/vdj_v1_hs_cd19_b",
              #     target="_blank")
              # ),
              title = "Load a small example database subset from 10x Genomics. More info in the Help tab",
              placement = "bottom"
            )
          }
        }
      )

      o_lda <- metaObserve2({
        shiny::req(input$lets_demo_analysis)
        isolate(metaExpr({
          sample_info_react$table=data.frame(Filename=c("Example10x_alakazam"),
                                             Sample="PBMC",
                                             Patient=c("HealthyDonor"),
                                             Group=c("CD19+"),
                                             Subgroup=c(""),
                                             Additional_info="10xGenomics")
          sample_info_react$test_status=T
          shinyjs::hide("raw_sample_file")
        }))
      })

      # observeEvent(
      #   input$lets_demo_analysis, {
      #     sample_info_react$table=data.frame(Filename=c("Example10x_alakazam"),
      #                                        Sample="PBMC",
      #                                        Patient=c("HealthyDonor"),
      #                                        Group=c("CD19+"),
      #                                        Subgroup=c(""),
      #                                        Additional_info="10xGenomics")
      #     sample_info_react$test_status=T
      #     hide("raw_sample_file")
      #   }
      # )

      sample_table_react <- metaReactive2(
        {
          validate(
            need(
              any(input$raw_sample_file != "") || sample_info_react$test_status ==T,
              "You need to upload your file with the sample information."
            )
          )
          metaExpr({
            if (is.null(..(input$raw_sample_file)) && ..(sample_info_react$test_status) ==F) {
              return(NULL)
            } else {

              sample_info <- ..(sample_info_react$table)
              if (all(
                c("Filename", "Sample", "Patient", "Group", "Subgroup") %in%
                colnames(sample_info)
              )) {
                return(sample_info)
              } else {
                error.table <- data.frame(
                  Error = paste(
                    "Modify and reupload again. The following column(s) are missing:",
                    paste(
                      c("Filename", "Sample", "Patient", "Group", "Subgroup")[which(
                        !(c("Filename", "Sample", "Patient", "Group", "Subgroup") %in%
                            colnames(sample_info))
                      )],
                      collapse = ", "
                    ),
                    sep = " "
                  )
                )
                return(error.table)
              }


            }
          })

        }, varname = "sample_table_react"
      )


      parameters_only_for_shinymeta_info <- metaReactive2(
        {
          validate(
            need(
              T,
              "Parameters for info report"
            )
          )
          metaExpr({
            if (shiny::isRunning()) {
              return(NULL)
            } else {

              author <- ..(input$username)
              usermail <- ..(input$usermail)
              user_0_comments <- ..(input$user_0_comments)
              user_1_comments <- ..(input$user_1_comments)
              user_2_comments <- ..(input$user_2_comments)
              user_3_comments <- ..(input$user_3_comments)
              user_4_comments <- ..(input$user_4_comments)
              user_5_comments <- ..(input$user_5_comments)

              return(list(author=author,
                          usermail=usermail,
                          user_0_comments=user_0_comments,
                          user_1_comments=user_1_comments,
                          user_2_comments=user_2_comments,
                          user_3_comments=user_3_comments,
                          user_4_comments=user_4_comments,
                          user_5_comments=user_5_comments))

            }
          })

        }, varname = "parameters_only_for_shinymeta_info"
      )


      parameters_only_for_shinymeta <- metaReactive2(
        {
          validate(
            need(
              any(input$raw_sample_file != "") || sample_info_react$test_status ==T,
              "Parameters for report"
            )
          )
          metaExpr({
            if (shiny::isRunning()) {
              return(NULL)
            } else {

              ## AAAAAAAAAAAAAAAAA #####
              folder_values <- list()
              folder_values$Featured <- ..(folder_values$Featured)
              folder_values$AIRR_parsed <- ..(folder_values$AIRR_parsed)
              Big_mem_values <- list()
              # Big_mem_values$Header=..(Big_mem_values$Header)
              # Big_mem_values$Short_DF=..(Big_mem_values$Short_DF)
              # Big_mem_values$Big_DF=..(Big_mem_values$Big_DF)
              Big_mem_values$Run=..(Big_mem_values$Run)
              # Big_mem_values$Patient_Sample=..(Big_mem_values$Patient_Sample)
              # Big_mem_values$VJs=..(Big_mem_values$VJs)
              Selection_values <- list()
              Features_values <- list()
              # reactiveValues(
              #   rows = NULL, columns = NULL, Scores = NULL, Variance_explained = NULL, UMAP = NULL,
              #   Parameters = NULL, Cell = F
              # )
              Exploration_values  <- list()
              Features_values <-  list(Current = 0, Total = 0)
              input=list()
              input$use_what <- ..(input$use_what)
              input$use_productive_or_not <- ..(input$use_productive_or_not)
              input$my_regions <- ..(input$my_regions )
              input$my_var_elements <- ..(input$my_var_elements)
              input$my_vars <- ..(input$my_vars)
              input$my_vartypes <- ..(input$my_vartypes)
              input$use_sharedVDJ <- ..(input$use_sharedVDJ)
              input$VJ_included <- ..(input$VJ_included)
              input$groups_selected <- ..(input$groups_selected)
              input$group_A <- ..(input$group_A)
              input$group_B <- ..(input$group_B)
              input$group_C <- ..(input$group_C)
              input$use_univlog <- ..(input$use_univlog)
              input$samples_selected <- ..(input$samples_selected)
              input$exclude_variables <- ..(input$exclude_variables)
              input$pval_type <- ..(input$pval_type)
              input$pval_cutoff <- ..(input$pval_cutoff)
              input$estimate_cutoff <- ..(input$estimate_cutoff)
              input$number_selected_vars <- ..(input$number_selected_vars )
              input$VJ_deselected <- ..(input$VJ_deselected)
              input$VDJ_normalized_per_size <- ..(input$VDJ_normalized_per_size)
              input$Rmut_filter <- ..(input$Rmut_filter)
              input$work_as_categories <- ..(input$work_as_categories)
              input$VDJ_maximize_clones <- ..(input$VDJ_maximize_clones)
              input$VDJ_normalized_per_sample <- ..(input$VDJ_normalized_per_sample)
              input$clonal_group <- ..(input$clonal_group)
              input$clonal_region_group <- ..(input$clonal_region_group)
              input$plot_color <- ..(input$plot_color)
              input$plot_color_expl <- ..(input$plot_color_expl)
              input$primary_color <- ..(input$primary_color)
              input$color_by <- ..(input$color_by)
              input$Exploration_plot_type <- ..(input$Exploration_plot_type )
              input$Exploration_plot_type_dim <- ..(input$Exploration_plot_type_dim)
              input$Selection_plot_type <- ..(input$Selection_plot_type )
              input$Selection_plot_type_dim <- ..(input$Selection_plot_type_dim)
              input$plot_feature <- ..(input$plot_feature)
              input$plot_color_feature <- ..(input$plot_color_feature)
              input$show_reconstructed <- ..(input$show_reconstructed)
              input$compare_opposites <- ..(input$compare_opposites)
              input$img_type <- ..(input$img_type)
              input$pixels_height <- ..(input$pixels_height)
              input$pixels_width <- ..(input$pixels_width)
              input$labelscale <- ..(input$labelscale)
              input$really_hide_points <- ..(input$really_hide_points)
              input$really_hide_points_clones <- ..(input$really_hide_points_clones)
              input$hide_points <- ..(input$hide_points)
              input$use_UMAP <- ..(input$use_UMAP)
              input$clonal_level_group <- ..(input$clonal_level_group)
              input$identity_clonal_group <- ..(input$identity_clonal_group)
              input$clones <- ..(input$clones)
              input$filter_clonal_group <- ..(input$filter_clonal_group)
              input$dominance_threshold <- ..(input$dominance_threshold)
              input$seed <- ..(input$seed)
              input$colorblind_mode <- ..(input$colorblind_mode)
              input$new_clonal_group <- ..(input$new_clonal_group)
              input$calculate_shared_clones <- ..(input$calculate_shared_clones)
              input$preinfolder_AIRR<- ..(input$preinfolder_AIRR)
              input$C_region_included_airr<- ..(input$C_region_included_airr)
              input$Dgene_reconstruct_airr<- ..(input$Dgene_reconstruct_airr)
              input$TCRBCR_input_file<- ..(input$TCRBCR_input_file)
              input$base_folder<- ..(input$base_folder)
              input$FWR1partial_airr<- ..(input$FWR1partial_airr)
              input$FWR4partial_airr<- ..(input$FWR4partial_airr)
              input$raw_sample_file <- ..(input$raw_sample_file)
              input$include_data_report <- ..(input$include_data_report)
              Big_mem_color_values <- list()
              sample_info_react <- list()
              sample_info_react$table <- ..(sample_info_react$table)
              sample_info_react$summary_status <- ..(sample_info_react$summary_status)
              sample_info_react$test_status <- ..(sample_info_react$test_status)
              sample_info_react$step_3 <- ..(sample_info_react$step_3)
              ID_selected_values <- list()
              ID_selected_values$subclones <- ..(ID_selected_values$subclones)
              ID_selected_values$clones <- ..(ID_selected_values$clones)
              ID_selected_values$intersection_samples <- ..(ID_selected_values$intersection_samples)
              ID_selected_values$clones_subclones_id <- ..(ID_selected_values$clones_subclones_id)
              if(..(input$include_data_report)) {
                folder_values$Featured <- "../Data/Dataset/Processed/2.Feature_determination"
                folder_values$AIRR_parsed <- "../Data/Dataset/Processed/1.Files_parsed"
                input$base_folder<- "../Data/Dataset/Meta"
                input$preinfolder_AIRR<- "../Data/Dataset/Raw"
              }
              # reactiveValues(
              #   rows = NULL, columns = NULL, Scores = NULL, Variance_explained = NULL, UMAP = NULL,
              #   Parameters = NULL
              # )
              return(list(input=input,
                          folder_values=folder_values,
                          Exploration_values=Exploration_values,
                          Selection_values=Selection_values,
                          Big_mem_values=Big_mem_values,
                          Big_mem_color_values=Big_mem_color_values,
                          ID_selected_values=ID_selected_values,
                          sample_info_react=sample_info_react,
                          Features_values=Features_values))

            }
          })

        }, varname = "parameters_only_for_shinymeta"
      )
      output$raw_sample_file_out <- metaRender(DT::renderDataTable,
                                               {
                                                 DT::datatable(..(sample_table_react()),
                                                 options = list(
                                                 pageLength = -1, lengthMenu = list(
                                                   c(15, -1),
                                                   c("15", "All")
                                                 ),
                                                 info = FALSE, initComplete = htmlwidgets::JS(
                                                   "function(settings, json) {", "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
                                                   "}"
                                                 )
                                               ))}
      )

      shinyDirChoose(input, "base_folder", roots = volumes(), session = session)



      # Steps 1&2. Select your FBM and associated files ######

      output$folder_information_AIRR_steps_1_to_3 <- DT::renderDataTable(
        {

          folder_information_table_AIRR_steps_1_to_3()

        }, options = list(
          pageLength = -1, lengthMenu = list(
            c(15, -1),
            c("15", "All")
          ),
          info = FALSE, initComplete = htmlwidgets::JS(
            "function(settings, json) {", "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
            "}"
          )
        )
      )

      folder_information_table_AIRR_steps_1_to_3 <- metaReactive(
        {

          if (all(
            unlist(input$FBM_folder) ==
            0
          )) {
            return(NULL)
          } else {
            sample_info <- sample_info_react$table

            sample_summary_status <- data.frame(
              Filename = paste(sample_info$Patient, sample_info$Group, sep = "."),
              Rds_and_Bk_files_found_in_folder = rep("No", nrow(sample_info))
            )

            for (dummy_row in grep("_dummy", sample_summary_status$Filename)) {
              sample_summary_status$Filename[dummy_row] <- sample_info$Filename[dummy_row]
            }
            if (!all(is.numeric(unlist(input$FBM_folder)))) {
              sample_summary_status$Rds_and_Bk_files_found_in_folder[intersect(
                which(
                  paste(sample_summary_status$Filename, ".bk", sep = "") %in%
                    list.files(
                      path = shinyFiles::parseDirPath(volumes(), input$FBM_folder),
                      full.names = F
                    )
                ),
                which(
                  paste(sample_summary_status$Filename, ".rds", sep = "") %in%
                    list.files(
                      path = shinyFiles::parseDirPath(volumes(), input$FBM_folder),
                      full.names = F
                    )
                )
              )] <- "Yes"

            }
            if (!all(is.numeric(unlist(input$FBM_folder))) &
                all(sample_summary_status$Rds_and_Bk_files_found_in_folder == "Yes")) {
              updateSelectizeInput(
                session, "make_dummy_menu", choices = sample_summary_status$Filename[which(!grepl("_dummy", sample_summary_status$Filename))]
              )

            }

            write.table(
              sample_info_react$table, file.path(
                shinyFiles::parseDirPath(volumes, input$base_folder),
                "Sample_summary.txt"
              ),
              append = F, row.names = F, col.names = T, sep = "\t", quote = F
            )

            return(sample_summary_status)
          }
        }, varname = "folder_information_table_AIRR_steps_1_to_3"
      )



      shinyDirChoose(input, "FBM_folder", roots = volumes(), session = session)

      ###BUGFIX: even with NO files, you can proceed

      o_Ff <- metaObserve2({
        shiny::req(input$FBM_folder)
        isolate(metaExpr({
          folder_values$Featured <- file.path(shinyFiles::parseDirPath(volumes(), input$FBM_folder))
          sample_info <- sample_info_react$table
          all_ok <- c()
          basenames <- basename(list.files(path = folder_values$Featured, pattern = glob2rx("*.rds")))
          prefixes <- file_path_sans_ext(basenames)
          prefixes <- unique(prefixes)[which(
            unique(prefixes) !=
              "Merged_bm"
          )]
          for (n_grouped_by in c(1:length(strsplit(prefixes[1], split = ".", fixed = T)))) {
            tmp_all_ok <- F
            for (sample_columna in c("Filename", "Sample", "Patient", "Group", "Subgroup")) {
              ### FUUUUUUUUUUUU #####
              if (all(
                sample_info[, sample_columna] %in% sapply(prefixes, function(z) strsplit(z, split = ".", fixed = T)[[1]][n_grouped_by])
              )) {
                tmp_all_ok <- T
              }
            }

            all_ok <- c(all_ok, tmp_all_ok)
          }

          #### TEMP ######
          all_ok <- T
          if (all(all_ok)) {
            can_show_button_to_step$step_2_3 <- TRUE
          } else {
            can_show_button_to_step$step_2_3 <- FALSE
          }
        }))
      })
      # observeEvent(
      #   input$FBM_folder, {
      #     folder_values$Featured <- file.path(parseDirPath(volumes(), input$FBM_folder))
      #     sample_info <- sample_info_react$table
      #     all_ok <- c()
      #     basenames <- basename(list.files(path = folder_values$Featured, pattern = glob2rx("*.rds")))
      #     prefixes <- file_path_sans_ext(basenames)
      #     prefixes <- unique(prefixes)[which(
      #       unique(prefixes) !=
      #         "Merged_bm"
      #     )]
      #     for (n_grouped_by in c(1:length(strsplit(prefixes[1], split = ".", fixed = T)))) {
      #       tmp_all_ok <- F
      #       for (sample_columna in c("Filename", "Sample", "Patient", "Group", "Subgroup")) {
      #         ### FUUUUUUUUUUUU #####
      #         if (all(
      #           sample_info[, sample_columna] %in% sapply(prefixes, function(z) strsplit(z, split = ".", fixed = T)[[1]][n_grouped_by])
      #         )) {
      #           tmp_all_ok <- T
      #         }
      #       }
      #
      #       all_ok <- c(all_ok, tmp_all_ok)
      #     }
      #
      #     #### TEMP ######
      #     all_ok <- T
      #     if (all(all_ok)) {
      #       can_show_button_to_step$step_2_3 <- TRUE
      #     } else {
      #       can_show_button_to_step$step_2_3 <- FALSE
      #     }
      #   }
      # )

      output$make_dummy <- renderUI(
        {
          if (is.null(input$raw_sample_file) ||
              all(is.numeric(unlist(input$FBM_folder)))) {
          } else {

            if (can_show_button_to_step$step_2_3 ) {
              shinyWidgets::materialSwitch("dummy_dataset", "Do you want to make a dummy dataset?", value = F)
            }

          }
        }
      )

      output$make_dummy_2 <- renderUI(
        {
          if (is.null(input$raw_sample_file) ||
              all(is.numeric(unlist(input$FBM_folder)))) {
          } else {

            if (can_show_button_to_step$step_3) {
              shinyWidgets::materialSwitch("dummy_dataset", "Do you want to make a dummy dataset?", value = F)
            }

          }
        }
      )

      o_fd <- metaObserve2({
        shiny::req(input$final_dummy)
        isolate(metaExpr({
          set.seed(input$seed)
          moment <- paste0(
            Sys.Date(), "-", data.table::hour(Sys.time()),
            ":", data.table::minute(Sys.time())
          )
          dummy_name <- paste("Dummy", input$make_dummy_menu, moment, sep = ".")

          og_FBM <- bigstatsr::big_attach(file.path(folder_values$Featured, paste0(input$make_dummy_menu, ".rds")))
          og_DT <- data.table::fread(file.path(folder_values$Featured, paste0(input$make_dummy_menu, ".info")))
          og_header <- data.table::fread(file.path(folder_values$Featured, paste0(input$make_dummy_menu, ".example")))

          list_new_DT <- list()

          if ("V_and_D_and_J" %in% colnames(og_DT)) {
            end <- grep("V_and_D_and_J", colnames(og_DT))
          } else {
            end <- grep("V_and_J", colnames(og_DT))
          }
          FWR1partial <- if (sum(og_FBM[, grep("NT_FWR1_length", colnames(og_header))]) ==
                             0) {
            T
          } else {
            F
          }
          FWR4partial <- if (sum(og_FBM[, grep("NT_FWR4_length", colnames(og_header))]) ==
                             0) {
            T
          } else {
            F
          }

          regions <- c("FWR1", "CDR1", "FWR2", "CDR2", "FWR3", "CDR3", "FWR4")
          init <- if (FWR1partial) {
            2
          } else {
            1
          }
          finit <- if (FWR4partial) {
            6
          } else {
            7
          }
          og_DT <- as.data.frame(og_DT)
          index_to_remove <- c()
          for (iteration in c(1:input$dummy_ratio)) {
            tmp_list <- og_DT[, 1:(end)]
            region_columns <- paste("NT", regions, sep = "_")
            number <- 1
            for (i_num in seq(
              from = 1, to = nrow(tmp_list),
              by = 2
            )) {
              tmp_list[i_num + 1, region_columns] <- tmp_list[i_num, region_columns]

              number <- number + 1
              if (og_DT$ORF[i_num] != og_DT$ORF[i_num + 1]) {
                index_to_remove <- c(index_to_remove, c(i_num, i_num + 1))
              } else {
                for (region in regions) {
                  intermezzo <- grep(region, regions)
                  for (mut in c("Replacement_muts_counts")) {
                    tmp_num_mut <- og_FBM[i_num + 1, which(
                      colnames(og_header) ==
                        paste("AA", region, mut, sep = "_")
                    )]
                    if (tmp_num_mut > 0) {
                      tmp_iter <- tmp_num_mut

                      used_positions <- c()
                      previous_sequences <- c()
                      while (tmp_iter > 0) {
                        sequence <- tmp_list[i_num + 1, paste("NT", region, sep = "_")]
                        tmp_tmp_position <- sample(
                          c(1:nchar(sequence))[which(
                            !(c(1:nchar(sequence)) %in%
                                used_positions)
                          )],
                          1, replace = F
                        )

                        sequence <- strsplit(sequence, split = "")[[1]]

                        sequence[[tmp_tmp_position]] <- sample(
                          c("A", "C", "G", "T")[which(
                            c("A", "C", "G", "T") !=
                              sequence[[tmp_tmp_position]]
                          )],
                          1
                        )
                        sequence <- paste(sequence, collapse = "")



                        if (intermezzo == init) {
                          whole_seq <- paste(
                            c(
                              sequence, tmp_list[i_num + 1, paste("NT", regions[(intermezzo + 1):(finit)], sep = "_")]
                            ),
                            collapse = ""
                          )
                        } else if (intermezzo == finit) {
                          whole_seq <- paste(
                            c(
                              tmp_list[i_num + 1, paste("NT", regions[init:(intermezzo - 1)], sep = "_")],
                              sequence
                            ),
                            collapse = ""
                          )
                        } else {
                          whole_seq <- paste(
                            c(
                              tmp_list[i_num + 1, paste("NT", regions[init:(intermezzo - 1)], sep = "_")],
                              sequence, tmp_list[i_num + 1, paste("NT", regions[(intermezzo + 1):(finit)], sep = "_")]
                            ),
                            collapse = ""
                          )
                        }
                        prev_whole_seq <- paste(
                          c(tmp_list[i_num + 1, paste("NT", regions, sep = "_")]),
                          collapse = ""
                        )

                        if (og_DT$ORF[i_num] != 1) {
                          whole_seq <- paste(
                            strsplit(whole_seq, split = "")[[1]][og_DT$ORF[i_num]:nchar(whole_seq)],
                            collapse = ""
                          )
                          prev_whole_seq <- paste(
                            strsplit(prev_whole_seq, split = "")[[1]][og_DT$ORF[i_num]:nchar(prev_whole_seq)],
                            collapse = ""
                          )
                        }
                        whole_seq <- Biostrings::DNAString(as.character(whole_seq))
                        prev_whole_seq <- Biostrings::DNAString(as.character(prev_whole_seq))


                        suppressWarnings(
                          {
                            whole_seq_AA <- as.character(Biostrings::translate(whole_seq))
                          }
                        )
                        suppressWarnings(
                          {
                            prev_whole_seq_AA <- as.character(Biostrings::translate(prev_whole_seq))
                          }
                        )

                        previous_sequences <- unique(c(previous_sequences, prev_whole_seq_AA))
                        if (!(whole_seq_AA %in% previous_sequences)) {
                          seq_diff <- T
                        } else {
                          seq_diff <- F
                        }

                        if (!grepl(pattern = "*", whole_seq_AA, fixed = T, perl = F)) {
                          productive <- T
                        } else {
                          productive <- F
                        }

                        if (mut == "Replacement_muts_counts" && seq_diff && productive) {
                          tmp_iter <- tmp_iter - 1
                          used_positions <- c(used_positions, tmp_tmp_position)
                          tmp_list[i_num + 1, paste("NT", region, sep = "_")] <- sequence
                        }


                      }
                    }

                  }
                }
              }


            }

            tmp_list$ID <- paste(tmp_list$ID, "Dummy", iteration, sep = "_")
            tmp_list$Patient <- paste("Dummy", tmp_list$Patient, sep = "_")
            tmp_list$Sample <- paste("Dummy", tmp_list$Sample, sep = "_")
            tmp_list$Group <- paste("Dummy", tmp_list$Group, sep = "_")
            tmp_list$Subgroup <- paste("Dummy", tmp_list$Subgroup, iteration, sep = "_")

            tmp_list <- tmp_list[c(1:nrow(tmp_list))[which(
              !c(1:nrow(tmp_list)) %in%
                index_to_remove
            )],
            ]
            list_new_DT[[length(list_new_DT) +
                           1]] <- tmp_list
          }
          new_DT <- data.table::rbindlist(list_new_DT)

          # new_FBM



          if(verbose) {
            AbSolution::Feature__dataset(
              path_base = folder_values$Featured, DF_to_parse = new_DT, name_DF_to_parse = dummy_name,
              FWR1partial = FWR1partial, FWR4partial = FWR4partial, standard = F
            )
          } else {
            suppressMessages(
              AbSolution::Feature__dataset(
                path_base = folder_values$Featured, DF_to_parse = new_DT, name_DF_to_parse = dummy_name,
                FWR1partial = FWR1partial, FWR4partial = FWR4partial, standard = F
              )
            )
          }


          sample_info <- sample_info_react$table


          first <- T
          for (coincidence in which(
            paste(sample_info$Patient, sample_info$Group, sep = ".") ==
            input$make_dummy_menu
          )) {
            tmp_dummy_sample_info <- unname(
              t(
                as.data.frame(
                  c(
                    dummy_name, paste(
                      sample_info[coincidence, 2:ncol(sample_info)],
                      "dummy", sep = "_"
                    )
                  )
                )
              )
            )
            colnames(tmp_dummy_sample_info) <- colnames(sample_info)
            if (first) {

              dummy_sample_info <- tmp_dummy_sample_info
            } else {
              dummy_sample_info <- rbind(dummy_sample_info, tmp_dummy_sample_info)
            }
            first <- F
          }


          sample_info <- rbind(sample_info, dummy_sample_info)
          sample_info_react$table <- sample_info
        }))
      })


      output$Conditional_Action_Move_to_Analysis_Real <- renderUI(
        {
          if (is.null(input$raw_sample_file) ||
              all(is.numeric(unlist(input$FBM_folder)))) {
          } else {

            if (can_show_button_to_step$step_2_3) {
              actionButton(
                "Move_to_analysis_real", "Continue the analysis", icon("forward"),
                style = "color: #fff; background-color: #18bc9c"
              )
            } else {
              HTML(
                "<p style=\"color:red\">Check your folder! <br> <br>

             Some files are missing/mislabelled, check your folder with your sample info table</p>"
              )
            }
          }
        }
      )

      # Page 1.AIRR-Seq conversion ######

      shinyDirChoose(input, "preinfolder_AIRR", roots = volumes(), session = session)

      dirname_AIRR <- metaReactive(
        {
          shinyFiles::parseDirPath(volumes, input$preinfolder_AIRR)
        }, varname = "dirname_AIRR"
      )

      folder_information_table_AIRR <- metaReactive(
        {
          if (all(
            unlist(input$preinfolder_AIRR) ==
            0
          )) {
            return(NULL)
          } else {
            sample_info <- sample_info_react$table
            sample_summary_status <- data.frame(
              Filename = sample_info$Filename, Found_in_folder = rep("No", nrow(sample_info))
            )

            if (!all(is.numeric(unlist(input$preinfolder_AIRR)))) {
              sample_summary_status$Found_in_folder[which(
                paste(sample_summary_status$Filename, ".tsv", sep = "") %in%
                  list.files(
                    path = shinyFiles::parseDirPath(volumes(), input$preinfolder_AIRR),
                    full.names = F
                  )
              )] <- "Yes"
            }
            return(sample_summary_status)
          }
        }, varname = "folder_information_table_AIRR"
      )

      example_parsed_sequence <- metaReactive(
        {
          if (all(
            unlist(input$preinfolder_AIRR) ==
            0
          ) || is.null(folder_information_table_AIRR()) ||
          any(folder_information_table_AIRR()$Found_in_folder != "Yes")
          ) {
            return(NULL)
          } else {
            inf_table <- (sample_table_react())
            RAW_sequence_DF=data.table::fread(file.path(shinyFiles::parseDirPath(volumes(),
                                                                                 input$preinfolder_AIRR),
                                                        paste(inf_table$Filename[1],
                                                              ".tsv",
                                                              sep = "")
            ),
            nrows=20
            )


            parsed_sequence_DF= tryCatch(
              {
                AbSolution::parse_AIRRSeq_file(
                  file = inf_table$Filename[1],
                  group = inf_table$Group[1],
                  patient = inf_table$Patient[1],
                  subgroup = inf_table$Subgroup[1],
                  sample = inf_table$Sample[1],
                  input_path = paste(
                    shinyFiles::parseDirPath(volumes(), input$preinfolder_AIRR),
                    "/", sep = ""
                  ),
                  C_region_included = input$C_region_included_airr,
                  FWR1partial = input$FWR1partial_airr,
                  FWR4partial = input$FWR4partial_airr,
                  output_path = paste(folder_values$AIRR_parsed, "/", sep = ""),
                  D_gene = input$Dgene_reconstruct_airr,
                  repertoire = input$TCRBCR_input_file,
                  is_example=T
                )
              } ,
              error = function(cond) {
                data.frame(Error="These parameters are not appropiate for these datasets!")
              })

            regions= c("NT_FWR1",  "NT_CDR1", "NT_FWR2",
                       "NT_CDR2", "NT_FWR3", "NT_CDR3", "NT_FWR4")

            if(ncol(parsed_sequence_DF)!=1) {


              for (ij in c(1:min(3, length(unique(parsed_sequence_DF$ID))))) {
                tmp_RAW_sequence_DF=RAW_sequence_DF[which(RAW_sequence_DF$sequence_id==unique(parsed_sequence_DF$ID)[ij]),]


                tmp_parsed_sequence <- data.frame(
                  Filename = rep( sample_info_react$table$Filename[1], 3),
                  Sequence_name=rep(unique(parsed_sequence_DF$ID)[ij], 3),
                  Sequence_type=c("Repertoire","Repertoire","Germline"),
                  Sequence_origin=c("Input_file","Parsed_by_AbSolution","Parsed_by_AbSolution"),
                  NT_FWR1=c(toupper(tmp_RAW_sequence_DF$fwr1), parsed_sequence_DF[c(ij*2, ij*2-1),]$NT_FWR1),
                  NT_CDR1=c(toupper(tmp_RAW_sequence_DF$cdr1), parsed_sequence_DF[c(ij*2, ij*2-1),]$NT_CDR1),
                  NT_FWR2=c(toupper(tmp_RAW_sequence_DF$fwr2),parsed_sequence_DF[c(ij*2, ij*2-1),]$NT_FWR2),
                  NT_CDR2=c(toupper(tmp_RAW_sequence_DF$cdr2),parsed_sequence_DF[c(ij*2, ij*2-1),]$NT_CDR2),
                  NT_FWR3=c(toupper(tmp_RAW_sequence_DF$fwr3),parsed_sequence_DF[c(ij*2, ij*2-1),]$NT_FWR3),
                  NT_CDR3=c(toupper(tmp_RAW_sequence_DF$cdr3),parsed_sequence_DF[c(ij*2, ij*2-1),]$NT_CDR3),
                  NT_FWR4=c(toupper(tmp_RAW_sequence_DF$fwr4),parsed_sequence_DF[c(ij*2, ij*2-1),]$NT_FWR4),
                  NT_Whole=c(paste0(toupper(tmp_RAW_sequence_DF$fwr1),
                                    toupper(tmp_RAW_sequence_DF$cdr1),
                                    toupper(tmp_RAW_sequence_DF$fwr2),
                                    toupper(tmp_RAW_sequence_DF$cdr2),
                                    toupper(tmp_RAW_sequence_DF$fwr3),
                                    toupper(tmp_RAW_sequence_DF$cdr3),
                                    toupper(tmp_RAW_sequence_DF$fwr4)),
                             sapply(c(ij*2, ij*2-1),
                                    function(z) paste(parsed_sequence_DF[z,regions],
                                                      collapse="")))
                )

                if(ij==1) {
                  parsed_sequence=tmp_parsed_sequence
                } else {
                  parsed_sequence=rbind(parsed_sequence, tmp_parsed_sequence)
                }

              }

            } else {
              parsed_sequence=parsed_sequence_DF
            }




            return(parsed_sequence)
          }
        }, varname = "example_parsed_sequence"
      )

      conditional_button_preprocess_AIRR <- metaReactive(
        {
          if (all(
            unlist(input$preinfolder_AIRR) ==
            0
          )) {
            HTML(
              "IMPORTANT! <i>Fill the required information and select the folder with the AIRR-Seq files. Only after that you will be able to proceed with the analysis.<i>"
            )
          } else {
            if (!is.null(folder_information_table_AIRR()) &&
                all(folder_information_table_AIRR()$Found_in_folder == "Yes")) {

              actionButton(
                "Preprocess_AIRR", "Process AIRR datasets", icon("forward"),
                style = "color: #fff; background-color: #18bc9c"
              )

            } else {
              HTML(
                "<p style=\"color:red\">Check your folder! <br> <br>

             Something is missing/mislabelled, perhaps the filenames (without the .tsv extension) are wrong or the files are missing</p>"
              )
            }
          }
        }, varname = "conditional_button_preprocess_AIRR"
      )

      output$Conditional_Action_Preprocess_AIRR <- renderUI(
        {
          conditional_button_preprocess_AIRR()
        }
      )

      output$folder_information_AIRR <- DT::renderDataTable(
        {
          folder_information_table_AIRR()
        }, options = list(
          pageLength = -1, lengthMenu = list(
            c(15, -1),
            c("15", "All")
          ),scrollX = TRUE,
          info = FALSE, initComplete = htmlwidgets::JS(
            "function(settings, json) {", "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
            "}"
          )
        )
      )

      output$DT_example_parsed_sequence <- DT::renderDataTable(
        {
          if(is.null(example_parsed_sequence())){

          } else if(ncol(example_parsed_sequence())==1) {
            DT::datatable(example_parsed_sequence(),
                          caption = "Preview")
          } else {
            DT::datatable(example_parsed_sequence(),
                          caption = "Preview",
                          options = list(
                            scrollX = TRUE,
                            pageLength=3)) |> DT::formatStyle(
                              c("NT_FWR1",  "NT_CDR1", "NT_FWR2",
                                "NT_CDR2", "NT_FWR3", "NT_CDR3", "NT_FWR4","NT_Whole"),
                              # target = 'column',
                              fontFamily = "Courier"
                            )
          }

        }, options = list(
          info = FALSE, initComplete = htmlwidgets::JS(
            "function(settings, json) {", "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
            "}"
          )
        )
      )

      can_show_button_to_step <- reactiveValues(step_2 = FALSE, step_2_AIRR = FALSE, step_3 = FALSE, step_2_3 = FALSE)
      folder_values <- reactiveValues(AIRR_parsed = "", Featured = "")


      o_PA <- metaObserve2({
        shiny::req(input$Preprocess_AIRR)
        isolate(metaExpr( {
          if(shiny::isRunning()) {
            shinyjs::hide("Preprocess_AIRR", anim = TRUE)
            session$sendCustomMessage(type = "testmessage", message = "Preprocessing")
          }

          if (shiny::isRunning()) {
            folder_values$AIRR_parsed <- file.path(
              shinyFiles::parseDirPath(volumes, input$base_folder),
              "1.Files_parsed"
            )
          } else {
            volumes <- shinyFiles::getVolumes()
          }


          unlink(folder_values$AIRR_parsed, recursive = TRUE)
          dir.create(folder_values$AIRR_parsed)
          if (shiny::isRunning() || !(input$include_data_report))  {
            write.table(
              sample_info_react$table, file.path(
                shinyFiles::parseDirPath(volumes, input$base_folder),
                "Sample_summary.txt"
              ),
              append = F, row.names = F, col.names = T, sep = "\t", quote = F
            )
          } else {
            write.table(
              sample_info_react$table,  file.path(input$base_folder, "Sample_summary.txt"),
              append = F, row.names = F, col.names = T, sep = "\t", quote = F
            )
          }

          if(shiny::isRunning()){
            shinyjs::show("pb_AIRR_vis", anim = TRUE)
          }


          if(shiny::isRunning()) {
            inf_table <- sample_table_react()
          } else {
            inf_table <- sample_info_react$table
          }


          for (row_number_sample_table in c(1:nrow(inf_table))) {


            AbSolution::parse_AIRRSeq_file(
              file = inf_table$Filename[row_number_sample_table], group = inf_table$Group[row_number_sample_table],
              patient = inf_table$Patient[row_number_sample_table], subgroup = inf_table$Subgroup[row_number_sample_table],
              sample = inf_table$Sample[row_number_sample_table], input_path = if(shiny::isRunning()|| !(input$include_data_report)) {
                paste(
                  shinyFiles::parseDirPath(volumes(), input$preinfolder_AIRR),
                  "/", sep = ""
                )
              } else {
                paste(input$preinfolder_AIRR, "/", sep = "")
              },
              C_region_included = input$C_region_included_airr, FWR1partial = input$FWR1partial_airr,
              FWR4partial = input$FWR4partial_airr, output_path = paste(folder_values$AIRR_parsed, "/", sep = ""),
              D_gene = input$Dgene_reconstruct_airr, repertoire = input$TCRBCR_input_file
            )

            if(shiny::isRunning()){
              shinyWidgets::updateProgressBar(
                session = session, id = "pb_AIRR", value = 100 * (row_number_sample_table/nrow(inf_table)),
                total = 100, title = paste("Process", trunc(100 * (row_number_sample_table/nrow(inf_table))/10))
              )
              Sys.sleep(0.1)
            }

          }

          if(shiny::isRunning()) {
            shinyjs::hide("pb_AIRR_vis", anim = TRUE)
            can_show_button_to_step$step_2_AIRR <- TRUE
          }

        }))
      })


      output$Conditional_Action_Move_to_2_AIRR <- renderUI(
        {
          if (can_show_button_to_step$step_2_AIRR) {
            actionButton("Move_to_2", "Next step", icon("angle-right"))
          }

        }
      )

      # 2.Sequence feature determination ####
      shinyjs::hide("pb_Feature_vis")

      Features_values <- reactiveValues(Current = 0, Total = 0)

      o_Featdet <- metaObserve2({
        shiny::req(input$Feature_determination)
        isolate(metaExpr({
          if(shiny::isRunning()) {
            shinyjs::hide("Feature_determination")
            session$sendCustomMessage(type = "testmessage", message = "Calculating features, be patient!")
          }

          if(shiny::isRunning()) {
            folder_values$Featured <- file.path(
              shinyFiles::parseDirPath(volumes(), input$base_folder),
              "2.Feature_determination"
            )
          }


          unlink(folder_values$Featured, recursive = TRUE)
          dir.create(folder_values$Featured)
          if(shiny::isRunning()) {
            shinyjs::show("pb_Feature_vis")
          }



          Features_values$Total <- AbSolution::Feature_1(
            if(shiny::isRunning()) {
              shinyFiles::parseDirPath(volumes, input$base_folder)
            } else {
              dirname(folder_values$Featured)
            },
            grouping_by = c("Patient", "Group")
          )



          List_dfs <- split(
            data.table::as.data.table(
              data.table::fread(
                if(shiny::isRunning()) {
                  file.path(
                    paste(
                      shinyFiles::parseDirPath(volumes, input$base_folder),
                      "/2.Feature_determination", sep = ""
                    ),
                    "IMGT_parsed_index_extended.txt"
                  )
                } else {
                  file.path(folder_values$Featured, "IMGT_parsed_index_extended.txt")
                },
                header = T, sep = "\t", quote = FALSE
              )
            ),
            by = c("Patient", "Group")
          )


          for (i in c(1:length(List_dfs))) {

            Features_values$Current <- Features_values$Current + 1
            if(verbose) {
              AbSolution::Feature__dataset(
                path_base = if(shiny::isRunning()) {shinyFiles::parseDirPath(volumes,
                                                                             input$base_folder)
                } else {dirname(folder_values$Featured)},
                DF_to_parse = List_dfs[[i]], name_DF_to_parse = names(List_dfs)[i],
                FWR1partial = input$FWR1partial_airr, FWR4partial = input$FWR4partial_airr
              )
            } else {
              suppressMessages(
                AbSolution::Feature__dataset(
                  path_base = if(shiny::isRunning()) {shinyFiles::parseDirPath(volumes,
                                                                               input$base_folder)
                  } else {dirname(folder_values$Featured)},
                  DF_to_parse = List_dfs[[i]], name_DF_to_parse = names(List_dfs)[i],
                  FWR1partial = input$FWR1partial_airr, FWR4partial = input$FWR4partial_airr
                )
              )

            }

          }
          rm(List_dfs)
          if(shiny::isRunning()) {
            shinyjs::hide("pb_Feature_vis")
            can_show_button_to_step$step_3 <- TRUE
          }

        }))
      })


      output$Conditional_Action_Move_to_3 <- renderUI(
        {
          if (can_show_button_to_step$step_3) {
            actionButton(
              "Move_to_analysis_real", "Proceed to next step: Visualization and exploration (be patient)",
              icon("angle-right")
            )
          }

        }
      )


      o_FeatCurr <- metaObserve2({
        shiny::req(Features_values$Current)
        isolate(metaExpr({
          if (Features_values$Current == 1) {
            shinyWidgets::updateProgressBar(
              session = session, id = "pb_Feature", value = 10, total = 100,
              title = paste("Now it will take some time, be patient!")
            )
          } else {
            shinyWidgets::updateProgressBar(
              session = session, id = "pb_Feature", value = 10 + 90 * ((Features_values$Current -
                                                                          1)/Features_values$Total), total = 100, title = paste("Now it will take some time, be patient!")
            )
          }

        }))
      })



      # 3.Dataset exploration and variable selection ####

      Big_mem_values <- reactiveValues(
        Header = NULL, Short_DF = NULL, Big_DF = NULL, Run = 0, Patient_Sample = NULL,
        VJs = NULL, categorical_newfields=c(),
        categorical_missing_ids=c()
      )
      Big_mem_color_values <- reactiveValues(V = NULL, D = NULL, J = NULL, VJ = NULL, VDJ = NULL,
                                             start_calculating=F, list_values=list())
      Selection_values <- reactiveValues(
        rows = NULL, columns = NULL, Scores = NULL, Variance_explained = NULL, UMAP = NULL,
        Parameters = NULL, Cell = F, Variance = NULL
      )
      Exploration_values <- reactiveValues(
        rows = NULL, columns = NULL, Scores = NULL, Variance_explained = NULL, UMAP = NULL,
        Parameters = NULL, Variance = NULL
      )

      o_MtAr <- metaObserve2({
        shiny::req(input$Move_to_analysis_real)
        isolate(metaExpr({


          info <- AbSolution::merge_FBMs(folder_values$Featured)

          Big_mem_values$Header <- info[[1]]
          Big_mem_values$Short_DF <- info[[2]]
          Big_mem_values$Short_DF$Patient_Sample <- paste(
            Big_mem_values$Short_DF$Patient, Big_mem_values$Short_DF$Sample,
            sep = "__"
          )
          Big_mem_values$Short_DF$Text_ID <- paste(
            Big_mem_values$Short_DF$ID, Big_mem_values$Short_DF$Sequence_type,
            sep = "_&_"
          )
          Big_mem_values$Big_DF <- info[[3]]
          if(shiny::isRunning()) {
            updateSliderInput(
              inputId = "Rmut_filter", min = min(
                Big_mem_values$Big_DF[, which(Big_mem_values$Header == "AA_Whole_Replacement_muts_counts")]
              ),
              max = max(
                Big_mem_values$Big_DF[, which(Big_mem_values$Header == "AA_Whole_Replacement_muts_counts")]
              )
            )
          }


          Big_mem_values$Patient_Sample <- unique(Big_mem_values$Short_DF$Patient_Sample)
          Big_mem_values$VJs <- sort(unique(Big_mem_values$Short_DF$V_and_J))
          counts_VJs <- table(Big_mem_values$VJs)/2
          names(Big_mem_values$VJs) <- sapply(
            Big_mem_values$VJs, function(z) paste(
              z, counts_VJs[which(
                names(counts_VJs) ==
                  z
              )],
              sep = " - Counts:"
            )
          )
          rm(info)


          if (any(grepl(",", Big_mem_values$Short_DF$Best_V, fixed = T))) {
            Big_mem_values$Short_DF <- Big_mem_values$Short_DF |>
              mutate(Best_V = gsub(",.*", "", .data$Best_V)) |>
              mutate(Best_J = gsub(",.*", "", .data$Best_J)) |>
              mutate(Best_D = gsub(",.*", "", .data$Best_D)) |>
              mutate(V_and_J = paste0(.data$Best_V, "_", .data$Best_J)) |>
              mutate(V_and_D_and_J = paste0(.data$Best_V, "_", .data$Best_D, "_", .data$Best_J))
          }
          Big_mem_color_values$list_values <- list()
          Big_mem_color_values$list_values[["V"]] <- list()
          Big_mem_color_values$list_values[["D"]] <- list()
          Big_mem_color_values$list_values[["J"]] <- list()

          for (gene in unique(Big_mem_values$Short_DF$Best_V)) {
            Big_mem_color_values$list_values[["V"]][[strsplit(
              strsplit(
                strsplit(gene, split = "*", fixed = T)[[1]][1],
                split = "/"
              )[[1]][1],
              split = "-"
            )[[1]][1]]] <- sort(
              unique(
                c(
                  Big_mem_color_values$list_values[["V"]][[strsplit(
                    strsplit(
                      strsplit(gene, split = "*", fixed = T)[[1]][1],
                      split = "/"
                    )[[1]][1],
                    split = "-"
                  )[[1]][1]]],
                  gene
                )
              )
            )
          }


          for (gene in sort(unique(Big_mem_values$Short_DF$Best_J))) {
            Big_mem_color_values$list_values[["J"]][[strsplit(
              strsplit(
                strsplit(gene, split = "*", fixed = T)[[1]][1],
                split = "/"
              )[[1]][1],
              split = "-"
            )[[1]][1]]] <- sort(
              unique(
                c(
                  Big_mem_color_values$list_values[["J"]][[strsplit(
                    strsplit(
                      strsplit(gene, split = "*", fixed = T)[[1]][1],
                      split = "/"
                    )[[1]][1],
                    split = "-"
                  )[[1]][1]]],
                  gene
                )
              )
            )
          }

          if ("Best_D" %in% colnames(Big_mem_values$Short_DF)) {
            for (gene in sort(unique(Big_mem_values$Short_DF$Best_D))) {
              Big_mem_color_values$list_values[["D"]][[strsplit(
                strsplit(
                  strsplit(gene, split = "*", fixed = T)[[1]][1],
                  split = "/"
                )[[1]][1],
                split = "-"
              )[[1]][1]]] <- sort(
                unique(
                  c(
                    Big_mem_color_values$list_values[["D"]][[strsplit(
                      strsplit(
                        strsplit(gene, split = "*", fixed = T)[[1]][1],
                        split = "/"
                      )[[1]][1],
                      split = "-"
                    )[[1]][1]]],
                    gene
                  )
                )
              )
            }
          }
          # for (gene in sort(unique(test$Best_V))) {
          # list_values[['V']][[strsplit(strsplit(gene, split='/')[[1]][1],
          # split='-')[[1]][1]]]=
          # sort(unique(c(list_values[['V']][[strsplit(strsplit(gene,
          # split='/')[[1]][1], split='-')[[1]][1]]], gene))) } for (gene in
          # sort(unique(test$Best_D))) {
          # list_values[['D']][[strsplit(strsplit(gene, split='/')[[1]][1],
          # split='-')[[1]][1]]]=
          # sort(unique(c(list_values[['D']][[strsplit(strsplit(gene,
          # split='/')[[1]][1], split='-')[[1]][1]]], gene))) } for (gene in
          # sort(unique(test$Best_J))) {
          # list_values[['J']][[strsplit(strsplit(gene, split='/')[[1]][1],
          # split='-')[[1]][1]]]=
          # sort(unique(c(list_values[['J']][[strsplit(strsplit(gene,
          # split='/')[[1]][1], split='-')[[1]][1]]], gene))) }

          Big_mem_color_values$start_calculating = T
          # test$Chain[1] seecol(usecol(pal = AbSolution:::Ab_palette(list_values,
          # vect_genes_comb=sort(unique(test$V_and_J)), type_values=c('VJ')),
          # n = 'all')) seecol(usecol(pal = AbSolution:::Ab_palette(list_values[['V']],
          # vect_genes_comb=NA, type_values=c('V')), n = 'all'))
          # seecol(usecol(pal = AbSolution:::Ab_palette(list_values[['D']],
          # vect_genes_comb=NA, type_values=c('D')), n = 'all'))
          # seecol(usecol(pal = AbSolution:::Ab_palette(list_values[['J']],
          # vect_genes_comb=NA, type_values=c('J')), n = 'all'))
          # seecol(usecol(pal = AbSolution:::Ab_palette(list_values[c('V','D','J')],
          # vect_genes_comb=sort(unique(test$V_and_D_and_J)),
          # type_values=c('VDJ')), n = 'all'))


          tmp_rows_cols <- AbSolution::filter_merged(
            FBM=Big_mem_values$Big_DF,
            merged_df=Big_mem_values$Short_DF,
            merged_header=Big_mem_values$Header,
            use_rgermline="Reconstructed germline" %in% input$use_what,
            use_repertoire="Repertoire" %in% input$use_what,
            use_productive="Productive" %in% input$use_productive_or_not,
            use_nonproductive= "Non-productive" %in% input$use_productive_or_not,
            my_regions= input$my_regions,
            my_var_elements=input$my_var_elements,
            my_vars=input$my_vars,
            my_vartypes=input$my_vartypes,
            use_sharedVDJ=input$use_sharedVDJ,
            V_J_to_use=input$VJ_included,
            groups=input$groups_selected,
            group_A=input$group_A,
            group_B=input$group_B,
            group_C=input$group_C,
            univlog=input$use_univlog,
            samples_to_keep=input$samples_selected,
            variables_to_remove=input$exclude_variables,
            pval_type=input$pval_type,
            pval_cutoff=input$pval_cutoff,
            estimate_cutoff=input$estimate_cutoff,
            number_selected_vars=input$number_selected_vars,
            VJ_deselected=input$VJ_deselected,
            VDJ_normalized_per_size=input$VDJ_normalized_per_size,
            R_mut_threshold_min=input$Rmut_filter[1],
            R_mut_threshold_max=input$Rmut_filter[2],
            to_compare_groups=input$work_as_categories,
            VDJ_maximize_clones=input$VDJ_maximize_clones,
            VDJ_normalized_per_sample=input$VDJ_normalized_per_sample,
            my_clone_def=input$clonal_group,
            seed=input$seed,
            chains=input$chains_selected,
            igsubtypes=input$subtype_selected
          )

          Exploration_values$rows <- tmp_rows_cols$ROWS
          Exploration_values$columns <- tmp_rows_cols$COLUMNS
          Selection_values$rows <- tmp_rows_cols$ROWS
          Selection_values$columns <- tmp_rows_cols$COLUMNS


          rm(tmp_rows_cols)
          ## PCAS base



        }))
      })

      recalculateVDJcoloring <- metaReactive(
        {
          list(Big_mem_color_values$start_calculating,
               input$seed,
               input$colorblind_mode)
        }, varname = "recalculateVDJcoloring"
      )


      o_calcColor <- metaObserve2({
        shiny::req(recalculateVDJcoloring())
        isolate(metaExpr({
          if (Big_mem_color_values$start_calculating) {
            suppressWarnings(
              {
                Big_mem_color_values$V <- AbSolution::Ab_palette(Big_mem_color_values$list_values[["V"]],
                                                                 vect_genes_comb = NA,
                                                                 type_values = c("V"), seed=input$seed,
                                                                 colorblind = input$colorblind_mode)

                Big_mem_color_values$J <- AbSolution::Ab_palette(Big_mem_color_values$list_values[["J"]],
                                                                 vect_genes_comb = NA,
                                                                 type_values = c("J"), seed=input$seed,
                                                                 colorblind = input$colorblind_mode)

                Big_mem_color_values$VJ <- AbSolution::Ab_palette(
                  Big_mem_color_values$list_values[c("V", "J")],
                  vect_genes_comb = sort(unique(Big_mem_values$Short_DF$V_and_J)),
                  type_values = c("VJ"), seed=input$seed,
                  colorblind = input$colorblind_mode
                )

                if ("Best_D" %in% colnames(Big_mem_values$Short_DF)) {
                  Big_mem_color_values$D <- AbSolution::Ab_palette(Big_mem_color_values$list_values[["D"]], vect_genes_comb = NA, type_values = c("D"), seed=input$seed,
                                                                   colorblind = input$colorblind_mode)
                  Big_mem_color_values$VDJ <- AbSolution::Ab_palette(
                    Big_mem_color_values$list_values[c("V", "D", "J")],
                    vect_genes_comb = sort(unique(Big_mem_values$Short_DF$V_and_D_and_J)),
                    type_values = c("VDJ"), seed=input$seed,
                    colorblind = input$colorblind_mode
                  )
                }
              }
            )
          }

        }))
      })


      toListenSelection <- metaReactive(
        {
          list(Selection_values$rows, Selection_values$columns, Selection_values$Parameters)
        }, varname = "toListenSelection"
      )


      o_tLS <- metaObserve2({
        shiny::req(toListenSelection())
        isolate(metaExpr({
          if (!is.null(Selection_values$rows) &&
              length(Selection_values$rows) >
              0 && !is.null(Selection_values$columns) &&
              length(Selection_values$columns) >
              0) {


            if (is.null(input$plot_color)) {
              selection <- "Group"
            } else {
              selection <- ..(input$plot_color)
            }


            tmp_PCA <- AbSolution::big_PCA(
              FBM = Big_mem_values$Big_DF, rows = Selection_values$rows, columns = Selection_values$columns
            )
            if(!is.null(tmp_PCA[[1]] )) {

              tmp_PCA[[1]] <- as.data.frame(tmp_PCA[[1]])
              colnames(tmp_PCA[[1]]) <- paste(
                "Dim_", c(1:5),
                sep = ""
              )

              tmp_PCA[[1]]$Color <- as.factor(unlist(Big_mem_values$Short_DF[Selection_values$rows, selection, with = FALSE ]))

              tmp_PCA[[1]]$Text_ID <- Big_mem_values$Short_DF[Selection_values$rows,
                                                              get("Text_ID")]
              tmp_PCA[[1]]$Seq_type <- factor(
                unlist(Big_mem_values$Short_DF[Selection_values$rows, get("Sequence_type")]),
                levels = c("Reconstructed_germline", "Repertoire")
              )
            }
            Selection_values$Scores <- tmp_PCA[[1]]

            Selection_values$Variance_explained <- tmp_PCA[[2]]
            Selection_values$Variance <- tmp_PCA[[3]]
            Big_mem_values$Run <- Big_mem_values$Run + 1

            if (input$use_UMAP == T && as.numeric(
              unname(
                strsplit(
                  as.character(utils::capture.output(
                    print(benchmarkme::get_ram(), unit_system = "iec"))
                  ),
                  split = " "
                )[[1]][1]
              )
            ) >
            round(
              length(Selection_values$rows) *
              length(Selection_values$columns) *
              8/2^{
                20
              }/1024, 2
            )) {

              tmp_umap <- as.data.frame(
                umap::umap(Big_mem_values$Big_DF[Selection_values$rows, Selection_values$columns])$layout
              )

              colnames(tmp_umap) <- paste(
                "Dim_", c(1:2),
                sep = ""
              )
              tmp_umap$Color <- as.factor(unlist(Big_mem_values$Short_DF[Selection_values$rows, selection, with = FALSE ]))
              tmp_umap$Text_ID <- paste(
                unlist(Big_mem_values$Short_DF[Selection_values$rows, get("ID")]),
                unlist(Big_mem_values$Short_DF[Selection_values$rows, get("Sequence_type")]),
                sep = "_&_"
              )
              tmp_umap$Seq_type <- factor(
                unlist(Big_mem_values$Short_DF[Selection_values$rows, get("Sequence_type")]),
                levels = c("Reconstructed_germline", "Repertoire")
              )
              dim(tmp_umap)
              Selection_values$UMAP <- tmp_umap
              if(shiny::isRunning()) {
                updateSelectInput(
                  session, "Selection_plot_type", choices = c(PCA = "PCA", UMAP = "UMAP"),
                  selected = "PCA"
                )
              }


            } else {
              if(shiny::isRunning()) {
                updateSelectInput(session, "Selection_plot_type", choices = c(PCA = "PCA"))
              }

            }
          } else {
            Selection_values$Scores <- NULL
            Selection_values$Variance_explained <- NULL
            Selection_values$UMAP <- NULL
          }

        }))
      })



      toListenSelection_plot_type <- metaReactive(
        {
          list(input$Selection_plot_type)
        }, varname = "toListenSelection_plot_type"
      )

      o_tLSpt <- metaObserve2({
        shiny::req(toListenSelection_plot_type())
        isolate(metaExpr({
          if (!is.null(input$Selection_plot_type)) {
            shiny::req(input$Selection_plot_type)
            if (input$Selection_plot_type == "PCA") {
              updateSelectizeInput(
                session, "Selection_plot_type_dim", choices = c(1:5),
                selected = c(1, 2)
              )
            } else if (input$Selection_plot_type == "UMAP") {
              updateSelectizeInput(
                session, "Selection_plot_type_dim", choices = c(1:2),
                selected = c(1, 2)
              )
            }
          }

        }))
      })


      toListenExploration_plot_type <- metaReactive(
        {
          list(input$Exploration_plot_type)
        }, varname = "toListenExploration_plot_type"
      )

      o_tLExpt <- metaObserve2({
        shiny::req(toListenExploration_plot_type())
        isolate(metaExpr( {
          if (!is.null(input$Exploration_plot_type)) {
            shiny::req(input$Exploration_plot_type)
            if (input$Exploration_plot_type == "PCA") {
              updateSelectizeInput(
                session, "Exploration_plot_type_dim", choices = c(1:5),
                selected = c(1, 2)
              )
            } else if (input$Exploration_plot_type == "UMAP") {
              updateSelectizeInput(
                session, "Exploration_plot_type_dim", choices = c(1:2),
                selected = c(1, 2)
              )
            }
          }
        }))
      })



      toListenExploration <- metaReactive(
        {
          list(
            Exploration_values$rows, Exploration_values$columns, Exploration_values$Parameters
          )
        }, varname = "toListenExploration"
      )

      toListenMenuOptions <-  metaReactive(
        {
          list(
            input$work_as_categories, input$group_A, input$group_B,
            input$use_sharedVDJ, input$VDJ_normalized_per_size,
            input$VDJ_maximize_clones, input$use_univlog
          )
        }, varname = "toListenMenuOptions"
      )
      observeEvent(
        toListenMenuOptions(), {
          shinyjs::hide("VJ_included", anim = TRUE)
          if (input$work_as_categories) {
            shinyjs::show("Group_selection", anim = TRUE)
            shinyjs::show("Group_comparison", anim = TRUE)
            if(length(input$group_A) >0 && length(input$group_B)>0) {
              shinyjs::show("use_sharedVDJ", anim = TRUE)
              if (input$use_sharedVDJ) {
                shinyjs::show("VJ_included", anim = TRUE)
                shinyjs::show("VDJ_normalized_per_size", anim = TRUE)
                if(input$VDJ_normalized_per_size) {
                  shinyjs::show("VDJ_maximize_clones", anim = TRUE)
                  shinyjs::show("VDJ_normalized_per_sample", anim = TRUE)
                  if(input$VDJ_maximize_clones){
                    shinyjs::show("my_clone_def", anim = TRUE)
                  } else {
                    shinyjs::hide("my_clone_def", anim = TRUE)
                  }
                } else {
                  shinyjs::hide("VDJ_maximize_clones", anim = TRUE)
                  shinyWidgets::updateMaterialSwitch(session, inputId="VDJ_maximize_clones",
                                                     value = FALSE)
                  shinyjs::hide("VDJ_normalized_per_sample", anim = TRUE)
                }
              } else {
                shinyjs::hide("VDJ_normalized_per_size", anim = TRUE)
                shinyjs::hide("VJ_included", anim = TRUE)
                shinyjs::hide("VDJ_normalized_per_sample", anim = TRUE)
                shinyjs::hide("VDJ_maximize_clones", anim = TRUE)
                shinyWidgets::updateMaterialSwitch(session, inputId="VDJ_maximize_clones",
                                                   value = FALSE)
              }

              shinyjs::show("use_univlog", anim = TRUE)

              if(input$use_univlog){
                shinyjs::show("pval_cutoff", anim = TRUE)
                shinyjs::show("pval_type", anim = TRUE)
                shinyjs::show("estimate_cutoff", anim = TRUE)
                shinyjs::show("number_selected_vars", anim = TRUE)
              } else {
                shinyjs::hide("pval_cutoff", anim = TRUE)
                shinyjs::hide("pval_type", anim = TRUE)
                shinyjs::hide("estimate_cutoff", anim = TRUE)
                shinyjs::hide("number_selected_vars", anim = TRUE)
              }
            } else {
              shinyjs::hide("use_sharedVDJ", anim = TRUE)
              shinyjs::hide("VJ_included", anim = TRUE)
              shinyjs::hide("use_univlog", anim = TRUE)
              shinyWidgets::updateMaterialSwitch(session, inputId="use_univlog",
                                                 value = FALSE)

            }
          } else {
            shinyjs::hide("Group_selection", anim = TRUE)
            shinyjs::hide("Group_comparison", anim = TRUE)
            shinyjs::hide("VDJ_normalized_per_size", anim = TRUE)
            shinyjs::hide("VDJ_normalized_per_sample", anim = TRUE)
            shinyjs::hide("VJ_included", anim = TRUE)
            shinyjs::hide("use_sharedVDJ", anim = TRUE)
            shinyjs::hide("use_univlog", anim = TRUE)
            shinyjs::hide("pval_type", anim = TRUE)
            shinyjs::hide("pval_cutoff", anim = TRUE)
            shinyjs::hide("estimate_cutoff", anim = TRUE)
            shinyjs::hide("number_selected_vars", anim = TRUE)
            shinyjs::hide("my_clone_def", anim = TRUE)
            shinyjs::hide("VDJ_maximize_clones", anim = TRUE)
            shinyWidgets::updateMaterialSwitch(session, inputId="use_univlog",
                                               value = FALSE)
            shinyWidgets::updateMaterialSwitch(session, inputId="use_sharedVDJ",
                                               value = FALSE)
            shinyWidgets::updateMaterialSwitch(session, inputId="VDJ_maximize_clones",
                                               value = FALSE)
            shinyWidgets::updateMaterialSwitch(session, inputId="VDJ_normalized_per_size",
                                               value = FALSE)
            # sortable:::update_rank_list("group_A", labels=NULL)
            # sortable:::update_rank_list("group_B", labels = NULL)
            # sortable:::update_rank_list("group_C", labels = NULL)
          }
        }
      )

      o_tLEx <- metaObserve2({
        shiny::req(toListenExploration())
        isolate(metaExpr({
          if (!is.null(Exploration_values$rows) &&
              length(Exploration_values$rows) >
              0 && !is.null(Exploration_values$columns) &&
              length(Exploration_values$columns) >
              0) {


            if (is.null(input$plot_color_expl)) {
              selection <- "Group"
            } else {
              selection <- input$plot_color_expl
            }

            tmp_PCAex <- AbSolution::big_PCA(
              FBM = Big_mem_values$Big_DF, rows = Exploration_values$rows, columns = Exploration_values$columns
            )

            if(!is.null(tmp_PCAex[[1]])) {

              tmp_PCAex[[1]] <- as.data.frame(tmp_PCAex[[1]])
              colnames(tmp_PCAex[[1]]) <- paste(
                "Dim_", c(1:5),
                sep = ""
              )
              tmp_PCAex[[1]]$Color <- as.factor(unlist(Big_mem_values$Short_DF[Exploration_values$rows, selection, with = FALSE ]))
              tmp_PCAex[[1]]$Text_ID <- Big_mem_values$Short_DF[Exploration_values$rows,
                                                                get("Text_ID")]
              tmp_PCAex[[1]]$Seq_type <- factor(
                unlist(Big_mem_values$Short_DF[Exploration_values$rows, get("Sequence_type")]),
                levels = c("Reconstructed_germline", "Repertoire")
              )
            }

            Exploration_values$Scores <- tmp_PCAex[[1]]

            Exploration_values$Variance_explained <- tmp_PCAex[[2]]

            Exploration_values$Variance <- tmp_PCAex[[3]]
            if (input$use_UMAP && as.numeric(
              unname(
                strsplit(
                  as.character(utils::capture.output(
                    print(benchmarkme::get_ram(), unit_system = "iec"))
                  ),
                  split = " "
                )[[1]][1]
              )
            ) >
            round(
              length(Exploration_values$rows) *
              length(Exploration_values$columns) *
              8/2^{
                20
              }/1024, 2
            )) {

              tmp_umap <- as.data.frame(
                umap::umap(Big_mem_values$Big_DF[Exploration_values$rows, Exploration_values$columns])$layout
              )

              colnames(tmp_umap) <- paste(
                "Dim_", c(1:2),
                sep = ""
              )
              tmp_umap$Color <- as.factor(unlist(Big_mem_values$Short_DF[Exploration_values$rows, selection, with = FALSE ]))
              tmp_umap$Text_ID <- paste(
                unlist(Big_mem_values$Short_DF[Exploration_values$rows, get("ID")]),
                unlist(Big_mem_values$Short_DF[Exploration_values$rows, get("Sequence_type")]),
                sep = "_&_"
              )
              tmp_umap$Seq_type <- factor(
                unlist(Big_mem_values$Short_DF[Exploration_values$rows, get("Sequence_type")]),
                levels = c("Reconstructed_germline", "Repertoire")
              )

              Exploration_values$UMAP <- tmp_umap
              if(shiny::isRunning()){
                updateSelectInput(
                  session, "Exploration_plot_type", choices = c(PCA = "PCA", UMAP = "UMAP"),
                  selected = "PCA"
                )
              }


              rm(tmp_umap)
            } else {

              if(shiny::isRunning()) {
                updateSelectInput(
                  session, "Exploration_plot_type", choices = c(PCA = "PCA"),
                  selected = "PCA"
                )
              }
            }



            rm(tmp_PCAex)
            Big_mem_values$Run <- Big_mem_values$Run + 1


            # }

          } else {
            Exploration_values$Scores <- NULL
            Exploration_values$Variance_explained <- NULL
            Exploration_values$UMAP <- NULL
            Exploration_values$Variance <- NULL

          }

        }))
      })



      ID_selected_values <- reactiveValues(subclones = NULL, clones = NULL, clones_subclones_id=NULL, intersection_samples = NULL)

      output$Selection_plot <- metaRender(renderPlotly,
                                          {
                                            if (!is.null(Selection_values$rows) &&
                                                length(Selection_values$rows) >  0 &&
                                                !is.null(Selection_values$columns) &&
                                                length(Selection_values$columns) > 0) {




                                              if (input$plot_color == "Best_V") {
                                                Sel_colors <- Big_mem_color_values$V
                                              } else if (input$plot_color == "Best_J") {
                                                Sel_colors <- Big_mem_color_values$J
                                              } else if (input$plot_color == "Best_D") {

                                                Sel_colors <- Big_mem_color_values$D

                                              } else if (input$plot_color == "V_and_D_and_J") {

                                                Sel_colors <- Big_mem_color_values$VDJ

                                              } else if (input$plot_color == "V_and_J") {

                                                Sel_colors <- Big_mem_color_values$VJ
                                              } else {


                                                if(is.null(Selection_values$Scores)) {
                                                  Sel_colors <- "black"
                                                } else {
                                                  Sel_colors <- AbSolution::Ab_palette(
                                                    list_values = unique(Selection_values$Scores$Color),
                                                    vect_genes_comb = NA, type_values = "cualitative", colorblind = input$colorblind_mode,
                                                    seed=input$seed
                                                  )
                                                }

                                                names(Sel_colors) <- sort(unique(Selection_values$Scores$Color))
                                              }


                                              Sel_colors_BORDER <- c(Sel_colors, "#2D2926", "#00979f", "#FFA500")
                                              names(Sel_colors_BORDER) <- c(
                                                names(Sel_colors),
                                                "User selected", "Counterpart", "Clones"
                                              )
                                              Sel_colors_BORDER_width <- c(
                                                rep(1, length(Sel_colors)),
                                                4, 4, 4
                                              )
                                              names(Sel_colors_BORDER_width) <- c(
                                                names(Sel_colors),
                                                "User selected", "Counterpart", "Clones"
                                              )
                                              Sel_colors_BORDER_symbol <- c("circle", "star-diamond")
                                              names(Sel_colors_BORDER_symbol) <- c("Repertoire", "Reconstructed_germline")
                                              shiny::req(input$Selection_plot_type)
                                              if (input$Selection_plot_type == "PCA") {

                                                test <- Selection_values$Scores

                                              } else if (input$Selection_plot_type == "UMAP") {
                                                test <- Selection_values$UMAP

                                              }

                                              test$Colors <- Sel_colors[match(test$Color, names(Sel_colors))]
                                              test$Color_border <- Sel_colors_BORDER[match(test$Selected, names(Sel_colors_BORDER))]
                                              test$Width_border <- Sel_colors_BORDER_width[match(test$Selected, names(Sel_colors_BORDER_width))]
                                              test$Symbol <- Sel_colors_BORDER_symbol[match(test$Seq_type, names(Sel_colors_BORDER_symbol))]

                                              if (length(input$Selection_plot_type_dim) <
                                                  2) {
                                                dim1 <- 1
                                                dim2 <- 2
                                              } else {
                                                dim1 <- as.numeric(input$Selection_plot_type_dim[1])
                                                dim2 <- as.numeric(input$Selection_plot_type_dim[2])

                                              }
                                              tmpDim_1 <- test[, which(
                                                colnames(test) ==
                                                  paste0("Dim_", dim1)
                                              )]
                                              tmpDim_2 <- test[, which(
                                                colnames(test) ==
                                                  paste0("Dim_", dim2)
                                              )]
                                              test$Dim_1 <- tmpDim_1
                                              test$Dim_2 <- tmpDim_2

                                              fig <- plotly::plot_ly(data = test, type = "scattergl", mode = "markers",
                                                                     colors = Sel_colors) |>
                                                plotly::config(
                                                  toImageButtonOptions = list(
                                                    format = input$img_type,
                                                    filename = paste("PCA_Selection_plot",Sys.time(),sep="_"),
                                                    width = input$pixels_width,
                                                    height = input$pixels_height,
                                                    scale=input$labelscale
                                                  )
                                                )|>
                                                plotly::layout( paper_bgcolor = "transparent",
                                                                xaxis = list(showline = TRUE,
                                                                             zeroline = FALSE),
                                                                yaxis = list(showline = TRUE,
                                                                             zeroline = FALSE))


                                              test <- test[(order(test$Color)),
                                              ]

                                              tmp_test_ns <- test[intersect(
                                                intersect(
                                                  which(test$Selected != "Clones"),
                                                  which(test$Selected != "User selected")
                                                ),
                                                which(test$Selected != "Counterpart")
                                              ),
                                              ]

                                              fig <- fig |>
                                                plotly::add_trace(
                                                  data = tmp_test_ns, x = ~Dim_1, y = ~Dim_2, opacity = 0.7, text = ~Text_ID,
                                                  key = ~Text_ID, color = ~Color, colors = Sel_colors, type = "scattergl",
                                                  mode = "markers", fill = ~"", hovertemplate = paste("<b>%{text}</b>"),
                                                  marker = list(sizemode = "diameter"),
                                                  name = tmp_test_ns$Color, legendgroup = tmp_test_ns$Color
                                                )


                                              tmp_test <- test[which(test$Selected == "Counterpart"),
                                              ]

                                              if (nrow(tmp_test) >
                                                  0) {
                                                fig <- fig |>
                                                  plotly::add_trace(
                                                    data = tmp_test, x = ~Dim_1, y = ~Dim_2, opacity = 0.7, text = ~Text_ID,
                                                    key = ~Text_ID, color = ~Color, colors = Sel_colors, symbol = ~Symbol,
                                                    type = "scattergl", mode = "markers", fill = ~"", hovertemplate = paste("<b>%{text}</b>"),
                                                    marker = list(line = list(color = "#00979f", width = tmp_test$Width_border[1])),
                                                    name = tmp_test$Color, showlegend = FALSE, legendgroup = tmp_test$Color
                                                  )
                                              }

                                              tmp_test <- test[which(test$Selected == "Clones"),
                                              ]

                                              if (nrow(tmp_test) >
                                                  0) {
                                                not_in_leg <- which(
                                                  unique(tmp_test$Color) %!in%
                                                    unique(tmp_test_ns$Color)
                                                )
                                                if (length(not_in_leg) >
                                                    0) {
                                                  tmp_tmp_test <- tmp_test[which(tmp_test$Color %in% not_in_leg),
                                                  ]

                                                  fig <- fig |>
                                                    plotly::add_trace(
                                                      data = tmp_tmp_test, x = ~Dim_1, y = ~Dim_2, opacity = 0.9, text = ~Text_ID,
                                                      key = ~Text_ID, color = ~Color, symbol = ~Symbol, type = "scattergl",
                                                      mode = "markers", fill = ~"", hovertemplate = paste("<b>%{text}</b>"),
                                                      marker = list(line = list(color = "#FFA500", width = tmp_tmp_test$Width_border[1])),
                                                      name = tmp_tmp_test$Color, showlegend = TRUE, legendgroup = tmp_tmp_test$Color
                                                    )

                                                }
                                                tmp_tmp_test <- tmp_test[which(tmp_test$Color %!in% not_in_leg),
                                                ]
                                                if (nrow(tmp_tmp_test) >
                                                    0) {
                                                  fig <- fig |>
                                                    plotly::add_trace(
                                                      data = tmp_tmp_test, x = ~Dim_1, y = ~Dim_2, opacity = 0.9, text = ~Text_ID,
                                                      key = ~Text_ID, color = ~Color, symbol = ~Symbol, type = "scattergl",
                                                      mode = "markers", fill = ~"", hovertemplate = paste("<b>%{text}</b>"),
                                                      marker = list(line = list(color = "#FFA500", width = tmp_tmp_test$Width_border[1])),
                                                      name = tmp_tmp_test$Color, showlegend = FALSE, legendgroup = tmp_tmp_test$Color
                                                    )
                                                }

                                              }

                                              tmp_test <- test[which(test$Selected == "User selected"),
                                              ]

                                              if (nrow(tmp_test) >
                                                  0) {
                                                not_in_leg <- which(
                                                  unique(tmp_test$Color) %!in%
                                                    unique(tmp_test_ns$Color)
                                                )
                                                if (length(not_in_leg) >
                                                    0) {
                                                  tmp_tmp_test <- tmp_test[which(tmp_test$Color %in% not_in_leg),
                                                  ]

                                                  fig <- fig |>
                                                    plotly::add_trace(
                                                      data = tmp_tmp_test, x = ~Dim_1, y = ~Dim_2, opacity = 0.9, text = ~Text_ID,
                                                      key = ~Text_ID, color = ~Color, symbol = ~Symbol, type = "scattergl",
                                                      mode = "markers", fill = ~"", hovertemplate = paste("<b>%{text}</b>"),
                                                      marker = list(line = list(color = "#2D2926", width = tmp_tmp_test$Width_border[1])),
                                                      name = tmp_tmp_test$Color, showlegend = TRUE, legendgroup = tmp_tmp_test$Color
                                                    )
                                                }
                                                tmp_tmp_test <- tmp_test[which(tmp_test$Color %!in% not_in_leg),
                                                ]
                                                if (nrow(tmp_tmp_test) >
                                                    0) {
                                                  fig <- fig |>
                                                    plotly::add_trace(
                                                      data = tmp_tmp_test, x = ~Dim_1, y = ~Dim_2, opacity = 0.9, text = ~Text_ID,
                                                      key = ~Text_ID, color = ~Color, symbol = ~Symbol, type = "scattergl",
                                                      mode = "markers", fill = ~"", hovertemplate = paste("<b>%{text}</b>"),
                                                      marker = list(line = list(color = "#2D2926", width = tmp_tmp_test$Width_border[1])),
                                                      name = tmp_tmp_test$Color, showlegend = FALSE, legendgroup = tmp_tmp_test$Color
                                                    )
                                                }

                                              }


                                              fig <- fig |>
                                                plotly::layout(
                                                  title = "Selection plot", plot_bgcolor = "#F8F9FA", legend = list(orientation = "v", y = 0),
                                                  showlegend = T, xaxis = list(
                                                    title = if (input$Selection_plot_type == "PCA") {
                                                      paste(
                                                        "Dim ", dim1, " (", 100 * Selection_values$Variance_explained[dim1],
                                                        "%)", sep = ""
                                                      )
                                                    } else if (input$Selection_plot_type == "UMAP") {
                                                      "Dim  1"
                                                    }
                                                  ),
                                                  yaxis = list(
                                                    title = if (input$Selection_plot_type == "PCA") {
                                                      paste(
                                                        "Dim ", dim2, " (", 100 * Selection_values$Variance_explained[dim2],
                                                        "%)", sep = ""
                                                      )
                                                    } else if (input$Selection_plot_type == "UMAP") {
                                                      "Dim 2"
                                                    }, dragmode = "lasso"
                                                  )
                                                ) |>
                                                plotly::config(
                                                  displaylogo = FALSE, modeBarButtonsToRemove = c("zoomIn2d", "zoomOut2d", "autoScale2d", "hoverCompareCartesian")
                                                )

                                              if(shiny::isRunning()){
                                                ID_selected_values$subclones <- ..(event_data("plotly_selected"))
                                              }



                                              fig |>
                                                plotly::toWebGL()
                                            } else {
                                              ggplotly(ggplot2::ggplot(data.frame()) + ggplot2::geom_point() +
                                                         ggplot2::xlim(0, 10) + ggplot2::ylim(0, 10) +
                                                         ggplot2::annotate("text", x=5, y=5.6, size=30,
                                                                           col="#00979f", label=":)" ) +
                                                         ggplot2::annotate("text", x=5, y=3, size=8, col="#2D2926",
                                                                           label="No significative variables \n Try other parameters!")
                                              )
                                            }

                                          }
      )


      output$Exploration_plot <- metaRender(renderPlotly,
                                            {
                                              if (!is.null(Exploration_values$rows) &&
                                                  length(Exploration_values$rows) >
                                                  0 && !is.null(Exploration_values$columns) &&
                                                  length(Exploration_values$columns) >
                                                  0) {



                                                if (input$plot_color_expl == "Best_V") {

                                                  Ex_colors <- Big_mem_color_values$V
                                                } else if (input$plot_color_expl == "Best_J") {
                                                  Ex_colors <- Big_mem_color_values$J
                                                } else if (input$plot_color_expl == "Best_D") {

                                                  Ex_colors <- Big_mem_color_values$D

                                                } else if (input$plot_color_expl == "V_and_D_and_J") {

                                                  Ex_colors <- Big_mem_color_values$VDJ

                                                } else if (input$plot_color_expl == "V_and_J") {

                                                  Ex_colors <- Big_mem_color_values$VJ
                                                } else {

                                                  # Sel_colors=unlist(branded_colors[c(1:length(unique(Selection_values$Scores$Color)))])

                                                  if(is.null(Exploration_values$Scores)) {
                                                    Ex_colors <- "black"
                                                  } else {
                                                    Ex_colors <- AbSolution::Ab_palette(
                                                      list_values = unique(Exploration_values$Scores$Color),
                                                      vect_genes_comb = NA, type_values = "cualitative", colorblind = input$colorblind_mode,
                                                      seed=input$seed
                                                    )
                                                  }

                                                  names(Ex_colors) <- sort(unique(Exploration_values$Scores$Color))
                                                }


                                                Ex_colors_BORDER <- c(Ex_colors, "#2D2926", "#00979f", "#FFA500")
                                                names(Ex_colors_BORDER) <- c(
                                                  names(Ex_colors),
                                                  "User selected", "Counterpart", "Clones"
                                                )
                                                Ex_colors_BORDER_width <- c(
                                                  rep(1, length(Ex_colors)),
                                                  4, 4, 4
                                                )
                                                names(Ex_colors_BORDER_width) <- c(
                                                  names(Ex_colors),
                                                  "User selected", "Counterpart", "Clones"
                                                )
                                                Ex_colors_BORDER_symbol <- c("circle", "star-diamond")
                                                names(Ex_colors_BORDER_symbol) <- c("Repertoire", "Reconstructed_germline")
                                                test_ex <- Exploration_values$Scores

                                                shiny::req(input$Exploration_plot_type)
                                                if (input$Exploration_plot_type == "PCA") {

                                                  test_ex <- Exploration_values$Scores
                                                } else if (input$Exploration_plot_type == "UMAP") {
                                                  test_ex <- Exploration_values$UMAP
                                                }

                                                test_ex$Colors <- Ex_colors[match(test_ex$Color, names(Ex_colors))]
                                                test_ex$Color_border <- Ex_colors_BORDER[match(test_ex$Selected, names(Ex_colors_BORDER))]
                                                test_ex$Width_border <- Ex_colors_BORDER_width[match(test_ex$Selected, names(Ex_colors_BORDER_width))]
                                                test_ex$Symbol <- Ex_colors_BORDER_symbol[match(test_ex$Seq_type, names(Ex_colors_BORDER_symbol))]


                                                if (length(input$Exploration_plot_type_dim) <
                                                    2) {
                                                  dim1 <- 1
                                                  dim2 <- 2
                                                } else {
                                                  dim1 <- as.numeric(input$Exploration_plot_type_dim[1])
                                                  dim2 <- as.numeric(input$Exploration_plot_type_dim[2])
                                                }
                                                tmpDim_1 <- test_ex[, which(
                                                  colnames(test_ex) ==
                                                    paste0("Dim_", dim1)
                                                )]
                                                tmpDim_2 <- test_ex[, which(
                                                  colnames(test_ex) ==
                                                    paste0("Dim_", dim2)
                                                )]
                                                test_ex$Dim_1 <- tmpDim_1
                                                test_ex$Dim_2 <- tmpDim_2

                                                fig <- plotly::plot_ly(data = test_ex, type = "scattergl", mode = "markers", colors = Ex_colors) |>
                                                  plotly::config(
                                                    toImageButtonOptions = list(
                                                      format = input$img_type,
                                                      filename = paste("PCA_Exploration_plot",Sys.time(),sep="_"),
                                                      width = input$pixels_width,
                                                      height = input$pixels_height,
                                                      scale=input$labelscale
                                                    )
                                                  )|>
                                                  plotly::layout( paper_bgcolor = "transparent",
                                                                  xaxis = list(showline = TRUE,
                                                                               zeroline = FALSE),
                                                                  yaxis = list(showline = TRUE,
                                                                               zeroline = FALSE))

                                                ## NEW
                                                test_ex <- test_ex[(order(test_ex$Color)),
                                                ]

                                                tmp_test_ns <- test_ex[intersect(
                                                  intersect(
                                                    which(test_ex$Selected != "Clones"),
                                                    which(test_ex$Selected != "User selected")
                                                  ),
                                                  which(test_ex$Selected != "Counterpart")
                                                ),
                                                ]

                                                fig <- fig |>
                                                  plotly::add_trace(
                                                    data = tmp_test_ns, x = ~Dim_1, y = ~Dim_2, opacity = 0.7, text = ~Text_ID,
                                                    key = ~Text_ID, color = ~Color, colors = Ex_colors, type = "scattergl",
                                                    mode = "markers", fill = ~"", hovertemplate = paste("<b>%{text}</b>"),
                                                    marker = list(sizemode = "diameter"),
                                                    name = tmp_test_ns$Color, legendgroup = tmp_test_ns$Color
                                                  )


                                                tmp_test <- test_ex[which(test_ex$Selected == "Counterpart"),
                                                ]

                                                if (nrow(tmp_test) >
                                                    0) {
                                                  fig <- fig |>
                                                    plotly::add_trace(
                                                      data = tmp_test, x = ~Dim_1, y = ~Dim_2, opacity = 0.7, text = ~Text_ID,
                                                      key = ~Text_ID, color = ~Color, colors = Ex_colors, symbol = ~Symbol,
                                                      type = "scattergl", mode = "markers", fill = ~"", hovertemplate = paste("<b>%{text}</b>"),
                                                      marker = list(line = list(color = "#00979f", width = tmp_test$Width_border[1])),
                                                      name = tmp_test$Color, showlegend = FALSE, legendgroup = tmp_test$Color
                                                    )
                                                }

                                                tmp_test <- test_ex[which(test_ex$Selected == "Clones"),
                                                ]

                                                if (nrow(tmp_test) >
                                                    0) {
                                                  not_in_leg <- which(
                                                    unique(tmp_test$Color) %!in%
                                                      unique(tmp_test_ns$Color)
                                                  )
                                                  if (length(not_in_leg) >
                                                      0) {
                                                    tmp_tmp_test <- tmp_test[which(tmp_test$Color %in% not_in_leg),
                                                    ]

                                                    fig <- fig |>
                                                      plotly::add_trace(
                                                        data = tmp_tmp_test, x = ~Dim_1, y = ~Dim_2, opacity = 0.9, text = ~Text_ID,
                                                        key = ~Text_ID, color = ~Color, symbol = ~Symbol, type = "scattergl",
                                                        mode = "markers", fill = ~"", hovertemplate = paste("<b>%{text}</b>"),
                                                        marker = list(line = list(color = "#FFA500", width = tmp_tmp_test$Width_border[1])),
                                                        name = tmp_tmp_test$Color, showlegend = TRUE, legendgroup = tmp_tmp_test$Color
                                                      )
                                                  }
                                                  tmp_tmp_test <- tmp_test[which(tmp_test$Color %!in% not_in_leg),
                                                  ]
                                                  if (nrow(tmp_tmp_test) >
                                                      0) {
                                                    fig <- fig |>
                                                      plotly::add_trace(
                                                        data = tmp_tmp_test, x = ~Dim_1, y = ~Dim_2, opacity = 0.9, text = ~Text_ID,
                                                        key = ~Text_ID, color = ~Color, symbol = ~Symbol, type = "scattergl",
                                                        mode = "markers", fill = ~"", hovertemplate = paste("<b>%{text}</b>"),
                                                        marker = list(line = list(color = "#FFA500", width = tmp_tmp_test$Width_border[1])),
                                                        name = tmp_tmp_test$Color, showlegend = FALSE, legendgroup = tmp_tmp_test$Color
                                                      )
                                                  }

                                                }

                                                tmp_test <- test_ex[which(test_ex$Selected == "User selected"),
                                                ]

                                                if (nrow(tmp_test) >
                                                    0) {
                                                  not_in_leg <- which(
                                                    unique(tmp_test$Color) %!in%
                                                      unique(tmp_test_ns$Color)
                                                  )
                                                  if (length(not_in_leg) >
                                                      0) {
                                                    tmp_tmp_test <- tmp_test[which(tmp_test$Color %in% not_in_leg),
                                                    ]

                                                    fig <- fig |>
                                                      plotly::add_trace(
                                                        data = tmp_tmp_test, x = ~Dim_1, y = ~Dim_2, opacity = 0.9, text = ~Text_ID,
                                                        key = ~Text_ID, color = ~Color, symbol = ~Symbol, type = "scattergl",
                                                        mode = "markers", fill = ~"", hovertemplate = paste("<b>%{text}</b>"),
                                                        marker = list(line = list(color = "#2D2926", width = tmp_tmp_test$Width_border[1])),
                                                        name = tmp_tmp_test$Color, showlegend = TRUE, legendgroup = tmp_tmp_test$Color
                                                      )
                                                  }
                                                  tmp_tmp_test <- tmp_test[which(tmp_test$Color %!in% not_in_leg),
                                                  ]
                                                  if (nrow(tmp_tmp_test) >
                                                      0) {
                                                    fig <- fig |>
                                                      plotly::add_trace(
                                                        data = tmp_tmp_test, x = ~Dim_1, y = ~Dim_2, opacity = 0.9, text = ~Text_ID,
                                                        key = ~Text_ID, color = ~Color, symbol = ~Symbol, type = "scattergl",
                                                        mode = "markers", fill = ~"", hovertemplate = paste("<b>%{text}</b>"),
                                                        marker = list(line = list(color = "#2D2926", width = tmp_tmp_test$Width_border[1])),
                                                        name = tmp_tmp_test$Color, showlegend = FALSE, legendgroup = tmp_tmp_test$Color
                                                      )
                                                  }

                                                }

                                                fig <- fig |>
                                                  plotly::layout(
                                                    title = "Exploration plot", plot_bgcolor = "#F8F9FA", legend = list(orientation = "v", y = -0),
                                                    showlegend = T, xaxis = list(
                                                      title = if (input$Exploration_plot_type == "PCA") {
                                                        paste(
                                                          "Dim ", dim1, " (", 100 * Exploration_values$Variance_explained[dim1],
                                                          "%)", sep = ""
                                                        )
                                                      } else if (input$Exploration_plot_type == "UMAP") {
                                                        "Dim  1"
                                                      }
                                                    ),
                                                    yaxis = list(
                                                      title = if (input$Exploration_plot_type == "PCA") {
                                                        paste(
                                                          "Dim ", dim2, " (", 100 * Exploration_values$Variance_explained[dim2],
                                                          "%)", sep = ""
                                                        )
                                                      } else if (input$Exploration_plot_type == "UMAP") {
                                                        "Dim  2"
                                                      }, dragmode = "lasso"
                                                    )
                                                  ) |>
                                                  plotly::config(
                                                    displaylogo = FALSE, modeBarButtonsToRemove = c("zoomIn2d", "zoomOut2d", "autoScale2d","hoverCompareCartesian")
                                                  )
                                                if(shiny::isRunning()){
                                                  ID_selected_values$subclones <- ..(event_data("plotly_selected"))
                                                }
                                                # print("Finidi")

                                                fig |>
                                                  plotly::toWebGL()

                                              } else {
                                                ggplotly(ggplot2::ggplot(data.frame()) + ggplot2::geom_point() +
                                                           ggplot2::xlim(0, 10) + ggplot2::ylim(0, 10) +
                                                           ggplot2::annotate("text", x=5, y=5.6, size=30,
                                                                             col="#00979f", label=":)" ) +
                                                           ggplot2::annotate("text", x=5, y=3, size=8, col="#2D2926",
                                                                             label="No significative variables \n Try other parameters!")
                                                )
                                              }

                                            }
      )
      output$Selection_plot_error <- renderUI(
        {
          if (!is.null(Selection_values$rows) &&
              length(Selection_values$rows) >
              0 && !is.null(Selection_values$columns) &&
              length(Selection_values$columns) >
              0) {

          } else {
            HTML(
              "<p style=\"color:red\">ERROR! <br>
              Your selection lead to no rows and/or columns with this analysis. <br>
              Change some parameters!
             </p>"
            )
          }

        }
      )
      output$Exploration_plot_error <- renderUI(
        {
          if (!is.null(Exploration_values$rows) &&
              length(Exploration_values$rows) >
              0 && !is.null(Exploration_values$columns) &&
              length(Exploration_values$columns) >
              0) {

          } else {
            HTML(
              "<p style=\"color:red\">ERROR! <br>
              Select any rows and/or columns<br>
              Or perhaps your selection lead to an error with this analysis
             </p>"
            )
          }

        }
      )

      output$Plot_color <- metaRender(renderUI,
                                      {
                                        if (is.null(Selection_values$rows) ||
                                            length(Selection_values$rows) ==
                                            0 || is.null(Selection_values$columns) ||
                                            length(Selection_values$columns) ==
                                            0) {
                                          selectInput(
                                            "plot_color", "Color by", choices = c("Sample", "Patient", "Group", "Subgroup", "Chain"),
                                            selected = "Group"
                                          )
                                        } else {
                                          sample_info <- ..(sample_info_react$table)
                                          choices_color <- c(colnames(sample_info))
                                          names(choices_color) <- c(colnames(sample_info))
                                          if ("Chain" %in% colnames(Big_mem_values$Short_DF) &&
                                              "Chain" %!in% choices_color) {
                                            tmp_names_color <- names(choices_color)
                                            choices_color <- c(choices_color, "Chain")
                                            names(choices_color) <- c(tmp_names_color, "Ig Chain")
                                          }
                                          tmp_names_choices_color <- names(choices_color)
                                          choices_color <- c(choices_color, "Best_V", "Best_J", "V_and_J")
                                          names(choices_color) <- c(
                                            tmp_names_choices_color, "Main V gene isoform", "Main J gene isoform",
                                            "Main V and J gene isoforms"
                                          )
                                          choices_color <- choices_color[which(choices_color != "Additional_info")]
                                          choices_color <- choices_color[which(choices_color != "Folder")]
                                          choices_color <- choices_color[which(choices_color != "Filename")]
                                          # choices_color=choices_color[which(choices_color !=
                                          # 'V_and_J')]



                                          if ("D_region" %in% colnames(Big_mem_values$Short_DF)) {
                                            tmp_names_color <- names(choices_color)
                                            choices_color <- c(choices_color, "Best_D")
                                            names(choices_color) <- c(tmp_names_color, "Main D gene isoform")
                                          }
                                          if ("V_and_D_and_J" %in% colnames(Big_mem_values$Short_DF)) {
                                            tmp_names_color <- names(choices_color)
                                            choices_color <- c(choices_color, "V_and_D_and_J")
                                            names(choices_color) <- c(tmp_names_color, "Main V, D and J gene isoforms")
                                          }
                                          if ("C_region" %in% colnames(Big_mem_values$Short_DF)) {
                                            tmp_names_color <- names(choices_color)
                                            choices_color <- c(choices_color, "C_region")
                                            names(choices_color) <- c(tmp_names_color, "Constant region")
                                          }
                                          if ("Clone_ID" %in% colnames(Big_mem_values$Short_DF)) {
                                            tmp_names_color <- names(choices_color)
                                            choices_color <- c(choices_color, "Clone_ID")
                                            names(choices_color) <- c(tmp_names_color, "Clone ID")
                                          }
                                          other_clones=which(grepl("Clone",colnames(Big_mem_values$Short_DF)))

                                          for(index_other_clones in other_clones) {
                                            if(colnames(Big_mem_values$Short_DF)[index_other_clones] != "Clone_ID"){
                                              tmp_names_color <- names(choices_color)
                                              choices_color <- c(choices_color, colnames(Big_mem_values$Short_DF)[index_other_clones])
                                              names(choices_color) <- c(tmp_names_color, colnames(Big_mem_values$Short_DF)[index_other_clones])
                                            }
                                          }


                                          if ("Dominance" %in% colnames(Big_mem_values$Short_DF)) {
                                            tmp_names_color <- names(choices_color)
                                            choices_color <- c(choices_color, "Dominance")
                                            names(choices_color) <- c(tmp_names_color, "Dominance")
                                          }

                                          if (length(Big_mem_values$categorical_newfields)>0) {
                                            tmp_names_color <- names(choices_color)
                                            choices_color <- c(choices_color, Big_mem_values$categorical_newfields)
                                            names(choices_color) <- c(tmp_names_color, Big_mem_values$categorical_newfields)
                                          }
                                          selectInput("plot_color", "Color by", choices = choices_color, selected = "Group")
                                        }

                                      }
      )


      output$Plot_color_expl <- metaRender(renderUI,
                                           {
                                             if (is.null(Exploration_values$rows) ||
                                                 length(Exploration_values$rows) ==
                                                 0 || is.null(Exploration_values$columns) ||
                                                 length(Exploration_values$columns) ==
                                                 0) {
                                               selectInput(
                                                 "plot_color_expl", "Color by", choices = c("Sample", "Patient", "Group", "Subgroup", "Chain"),
                                                 selected = "Group"
                                               )
                                             } else {
                                               sample_info <- sample_info_react$table

                                               choices_color <- c(colnames(sample_info))
                                               names(choices_color) <- c(colnames(sample_info))
                                               if ("Chain" %in% colnames(Big_mem_values$Short_DF) &&
                                                   "Chain" %!in% choices_color) {
                                                 tmp_names_color <- names(choices_color)
                                                 choices_color <- c(choices_color, "Chain")
                                                 names(choices_color) <- c(tmp_names_color, "Ig Chain")
                                               }
                                               tmp_names_choices_color <- names(choices_color)
                                               choices_color <- c(choices_color, "Best_V", "Best_J", "V_and_J")
                                               names(choices_color) <- c(
                                                 tmp_names_choices_color, "Main V gene isoform", "Main J gene isoform",
                                                 "Main V and J gene isoforms"
                                               )
                                               choices_color <- choices_color[which(choices_color != "Additional_info")]
                                               choices_color <- choices_color[which(choices_color != "Folder")]
                                               choices_color <- choices_color[which(choices_color != "Filename")]
                                               # choices_color=choices_color[which(choices_color !=
                                               # 'V_and_J')]


                                               if ("D_region" %in% colnames(Big_mem_values$Short_DF)) {
                                                 tmp_names_color <- names(choices_color)
                                                 choices_color <- c(choices_color, "Best_D")
                                                 names(choices_color) <- c(tmp_names_color, "Main D gene isoform")
                                               }
                                               if ("V_and_D_and_J" %in% colnames(Big_mem_values$Short_DF)) {
                                                 tmp_names_color <- names(choices_color)
                                                 choices_color <- c(choices_color, "V_and_D_and_J")
                                                 names(choices_color) <- c(tmp_names_color, "Main V, D and J gene isoforms")
                                               }
                                               if ("C_region" %in% colnames(Big_mem_values$Short_DF)) {
                                                 tmp_names_color <- names(choices_color)
                                                 choices_color <- c(choices_color, "C_region")
                                                 names(choices_color) <- c(tmp_names_color, "Constant region")
                                               }
                                               if ("Clone_ID" %in% colnames(Big_mem_values$Short_DF)) {
                                                 tmp_names_color <- names(choices_color)
                                                 choices_color <- c(choices_color, "Clone_ID")
                                                 names(choices_color) <- c(tmp_names_color, "Clone ID")
                                               }

                                               other_clones=which(grepl("Clone",colnames(Big_mem_values$Short_DF)))

                                               for(index_other_clones in other_clones) {
                                                 if(colnames(Big_mem_values$Short_DF)[index_other_clones] != "Clone_ID"){
                                                   tmp_names_color <- names(choices_color)
                                                   choices_color <- c(choices_color, colnames(Big_mem_values$Short_DF)[index_other_clones])
                                                   names(choices_color) <- c(tmp_names_color, colnames(Big_mem_values$Short_DF)[index_other_clones])
                                                 }
                                               }

                                               if ("Dominance" %in% colnames(Big_mem_values$Short_DF)) {
                                                 tmp_names_color <- names(choices_color)
                                                 choices_color <- c(choices_color, "Dominance")
                                                 names(choices_color) <- c(tmp_names_color, "Dominance")
                                               }


                                               if (length(Big_mem_values$categorical_newfields)>0) {
                                                 tmp_names_color <- names(choices_color)
                                                 choices_color <- c(choices_color, Big_mem_values$categorical_newfields)
                                                 names(choices_color) <- c(tmp_names_color, Big_mem_values$categorical_newfields)
                                               }

                                               # if ('Cell_ID' %in% colnames(Big_mem_values$Short_DF)) {
                                               # tmp_names_color=names(choices_color)
                                               # choices_color=c(choices_color, 'Cell_ID')
                                               # names(choices_color)=c(tmp_names_color, 'Cell ID') }



                                               # choices_color=choices_color[match(choices_color,c('Sample',\t'Patient',\t'Group',\t'Subgroup',\t'Chain',
                                               # 'Best_V', 'Best_J', 'V_and_J',
                                               # 'D_region','C_region','Clone_ID','Cell_ID'))]
                                               # tst[match(c('D','AD','A','C')[which(c('D','AD','A','C')
                                               # %in%tst)],tst)]
                                               selectInput("plot_color_expl", "Color by", choices = choices_color, selected = "Group")
                                             }

                                           }
      )

      toListenColorsSelection <- metaReactive(
        {
          list(
            input$plot_color, input$plot_color_expl, ID_selected_values$subclones,
            ID_selected_values$clones, Exploration_values$rows, Exploration_values$columns,
            Selection_values$rows, Selection_values$columns, Big_mem_values$Run
          )
        }, varname = "toListenColorsSelection"
      )



      o_tLCSel <- metaObserve2({
        shiny::req(toListenColorsSelection())
        isolate(metaExpr({
          if (!is.null(input$plot_color_expl) &&
              !is.null(input$plot_color)) {

            tmp_sp <- (unlist(Big_mem_values$Short_DF[Exploration_values$rows, get(input$plot_color_expl)]))


            tmp_sp[which(is.na(tmp_sp))] <- "NotAvailable"
            tmp_sp[which((tmp_sp) == "  ")] <- "NotAvailable"
            tmp_sp[which((tmp_sp) == " ")] <- "NotAvailable"
            tmp_sp[which((tmp_sp) == "")] <- "NotAvailable"
            tmp_sp[which((tmp_sp) == "NotAvailable")] <- "NotAvailable"

            tmp_sp <- unname(tmp_sp)

            if (!is.null(Exploration_values$Scores)) {
              Exploration_values$Scores$Color <- tmp_sp
            }

            if (!is.null(Exploration_values$UMAP) && as.numeric(
              unname(
                strsplit(
                  as.character(utils::capture.output(
                    print(benchmarkme::get_ram(), unit_system = "iec"))
                  ),
                  split = " "
                )[[1]][1]
              )
            ) >
            round(
              length(Exploration_values$rows) *
              length(Exploration_values$columns) *
              8/2^{
                20
              }/1024, 2
            )) {
              Exploration_values$UMAP$Color <- tmp_sp
            }


            old_sp <- NULL
            all_sp <- NULL
            if (!is.null(ID_selected_values$subclones) && !is.null(Exploration_values$Scores)) {
              old_sp <- tmp_sp

              tmp_sp[Exploration_values$Scores$Text_ID %in% ID_selected_values$subclones$key] <- "User selected"



              counterpart_start <- paste(
                sapply(
                  ID_selected_values$subclones$key, function(z) strsplit(z, split = "_&_")[[1]][1]
                ),
                "_&_", sep = ""
              )
              tmp_txt_id <- paste(
                sapply(
                  Exploration_values$Scores$Text_ID, function(z) strsplit(z, split = "_&_")[[1]][1]
                ),
                "_&_", sep = ""
              )
              tmp_sp[intersect(
                which(tmp_txt_id %in% counterpart_start),
                which(tmp_sp != "User selected")
              )] <- "Counterpart"

            }

            tmp_sp <- unname(tmp_sp)

            if (!is.null(Exploration_values$Scores)) {
              Exploration_values$Scores$Selected <- tmp_sp
            }


            if (!is.null(Exploration_values$UMAP) && as.numeric(
              unname(
                strsplit(
                  as.character(utils::capture.output(
                    print(benchmarkme::get_ram(), unit_system = "iec"))
                  ),
                  split = " "
                )[[1]][1]
              )
            ) >
            round(
              length(Exploration_values$rows) *
              length(Exploration_values$columns) *
              8/2^{
                20
              }/1024, 2
            )) {
              Exploration_values$UMAP$Selected <- tmp_sp
            }


            tmp_sp <- (unlist(Big_mem_values$Short_DF[Selection_values$rows, get(input$plot_color)]))

            # if(input$plot_color %in% c('Best_V', 'Best_J')) {
            # tmp_sp=sapply(tmp_sp, function(z) strsplit(z,split='*',
            # fixed=T)[[1]][1]) } else if (input$plot_color == 'V_and_J') {
            # tmp_sp=sapply(tmp_sp, function(z) paste(strsplit(z,split='*',
            # fixed=T)[[1]][1], strsplit(strsplit(z,split='_',
            # fixed=T)[[1]][2], split='*', fixed=T)[[1]][1], sep=' & ')) }

            tmp_sp[which(is.na(tmp_sp))] <- "Not specified"
            tmp_sp[which((tmp_sp) == "  ")] <- "Not specified"
            tmp_sp[which((tmp_sp) == " ")] <- "Not specified"
            tmp_sp[which((tmp_sp) == "")] <- "Not specified"
            tmp_sp[which((tmp_sp) == "NotAvailable")] <- "Not specified"


            if (!is.null(Selection_values$Scores)) {

              Selection_values$Scores$Color <- tmp_sp

            }

            if (!is.null(Selection_values$UMAP) && as.numeric(
              unname(
                strsplit(
                  as.character(utils::capture.output(
                    print(benchmarkme::get_ram(), unit_system = "iec"))
                  ),
                  split = " "
                )[[1]][1]
              )
            ) >
            round(
              length(Selection_values$rows) *
              length(Selection_values$columns) *
              8/2^{
                20
              }/1024, 2
            )) {
              Selection_values$UMAP$Color <- tmp_sp
            }

            if (!is.null(ID_selected_values$subclones) && !is.null(Selection_values$Scores)) {
              tmp_sp[Selection_values$Scores$Text_ID %in% ID_selected_values$subclones$key] <- "User selected"

              all_sp <- rep("Non-selected", nrow(Big_mem_values$Short_DF))
              all_sp[Big_mem_values$Short_DF$Text_ID %in% ID_selected_values$subclones$key] <- "User selected"

              counterpart_start <- paste(
                sapply(
                  ID_selected_values$subclones$key, function(z) strsplit(z, split = "_&_")[[1]][1]
                ),
                "_&_", sep = ""
              )
              tmp_txt_id <- paste(
                sapply(
                  Selection_values$Scores$Text_ID, function(z) strsplit(z, split = "_&_")[[1]][1]
                ),
                "_&_", sep = ""
              )
              tmp_sp[intersect(
                which(tmp_txt_id %in% counterpart_start),
                which(tmp_sp != "User selected")
              )] <- "Counterpart"

            }

            all_sp <- unname(all_sp)

            if (!is.null(all_sp)) {
              Big_mem_values$Short_DF$Selected <- all_sp
            } else if ("Selected" %in% colnames(Big_mem_values$Short_DF)) {
              Big_mem_values$Short_DF[, "Selected" := NULL]
            }

            if (!is.null(ID_selected_values$clones)) {

              if (input$clonal_group == "Clone_ID") {
                clonepart_start <- (Big_mem_values$Short_DF)[intersect(
                  Selection_values$rows, which(
                    (Big_mem_values$Short_DF)[, get("Clone_ID")] %in%
                      ID_selected_values$clones$key
                  )
                ),
                get("ID")]

              } else {

                clonepart_start <- (Big_mem_values$Short_DF)[intersect(
                  Selection_values$rows, which(
                    (Big_mem_values$Short_DF)[, get(
                      input$clonal_group
                    )] %in%
                      ID_selected_values$clones$key
                  )
                ),
                get("ID")]
              }
              tmp_txt_id <- sapply(
                Selection_values$Scores$Text_ID, function(z) strsplit(z, split = "_&_")[[1]][1]
              )


              tmp_sp[which(tmp_txt_id %in% clonepart_start)] <- "Clones"
              ID_selected_values$clones_subclones_id=clonepart_start


            }

            all_sp_clones=rep("Non-selected", nrow(Big_mem_values$Short_DF))
            if (!is.null(ID_selected_values$clones)) {
              all_sp_clones[intersect(
                Selection_values$rows, which(
                  (Big_mem_values$Short_DF)[, get(
                    input$clonal_group
                  )] %in%
                    ID_selected_values$clones$key
                )
              )] <- "User selected"

            }

            all_sp_clones <- unname(all_sp_clones)
            Big_mem_values$Short_DF$Selected_clones <- all_sp_clones

            if (!is.null(Selection_values$Scores)) {
              Selection_values$Scores$Selected <- tmp_sp
            }

            if (  !is.null(Selection_values$UMAP) && as.numeric(
              unname(
                strsplit(
                  as.character(utils::capture.output(
                    print(benchmarkme::get_ram(), unit_system = "iec"))
                  ),
                  split = " "
                )[[1]][1]
              )
            ) >
            round(
              length(Selection_values$rows) *
              length(Selection_values$columns) *
              8/2^{
                20
              }/1024, 2
            )) {
              Selection_values$UMAP$Selected <- tmp_sp

            }

          }

        }))
      })


      observeEvent(input$categorical_covariable_file,{
        shiny::req(input$categorical_covariable_file)
        file <- read.csv(input$categorical_covariable_file$datapath)
        validate(need(any("Sequence_ID" %in% colnames(file)), "Please include a Sequence_ID colname"))
        Big_mem_values$categorical_newfields=c(Big_mem_values$categorical_newfields, colnames(file)[which(colnames(file)!= "Sequence_ID")] )
        Big_mem_values$Short_DF= merge(Big_mem_values$Short_DF, file, by.x = "ID",by.y = "Sequence_ID", all.x= TRUE, sort=FALSE)
        Big_mem_values$categorical_missing_ids=file$Sequence_ID[which(!(file$Sequence_ID %in% Big_mem_values$Short_DF$ID))]
        Big_mem_values$Short_DF[is.na(Big_mem_values$Short_DF)] ="Not_specified"
      })

      output$results_categorical <- renderUI(
        {

          if(length(Big_mem_values$categorical_newfields)>0) {
            if (length(Big_mem_values$categorical_missing_ids) == 0) {
              HTML("<p style=\"color:green\"> <b>All is OK!</b> <br>
              All Sequence_IDs were matched in the AbSolution dataset.</p>"
              )
            } else {
              HTML(
                paste("<p style=\"color:orange\"> <b>WARNING!</b> <br>
              These Sequence_IDs were not found in the AbSolution dataset",
                      paste(Big_mem_values$categorical_missing_ids,
                            collapse="<br>"),
                      "</p>", paste="")


              )
            }
          }


        }
      )

      observeEvent(
        ID_selected_values$subclones, {
          if (!is.null(ID_selected_values$subclones)) {
            # ID_selected_values$clones=NULL

          }
        }
      )
      observeEvent(
        ID_selected_values$clones, {
          if (!is.null(ID_selected_values$clones)) {
            # ID_selected_values$subclones=NULL
          }
        }
      )

      merged_data <- eventReactive(
        ID_selected_values$subclones, {

          if (!is.null(ID_selected_values$subclones) &&
              !is.null(ID_selected_values$clones)) {
            keep_cols <- colnames(Big_mem_values$Short_DF)
            keys_all <- paste(Big_mem_values$Short_DF$ID,
                              Big_mem_values$Short_DF$Sequence_type,
                              sep = "_&_")
            idx_sub <- which(keys_all %in% ID_selected_values$subclones$key)
            Big_mem_values$Short_DF[unique(
              c(idx_sub,
                Selection_values$rows[which(Selection_values$Scores$Selected == "Clones")]
              )
            ), keep_cols, with = FALSE ]

          } else if (!is.null(ID_selected_values$subclones)) {
            keep_cols <- colnames(Big_mem_values$Short_DF)
            keys_all <- paste(Big_mem_values$Short_DF$ID,
                              Big_mem_values$Short_DF$Sequence_type,
                              sep = "_&_")
            idx_sub <- which(keys_all %in% ID_selected_values$subclones$key)
            Big_mem_values$Short_DF[idx_sub, keep_cols, with = FALSE ]
          } else if (!is.null(ID_selected_values$clones)) {
            keep_cols <- colnames(Big_mem_values$Short_DF)

            Big_mem_values$Short_DF[Selection_values$rows[which(Selection_values$Scores$Selected == "Clones")], keep_cols, with = FALSE ]
          }


        }
      )  #eventReactive


      output$Ab_table <- metaRender(DT::renderDataTable,
                                    {
                                      rendered_table <- merged_data()
                                      DT::datatable(
                                        rendered_table, caption = "Selected sequences",
                                        extensions = "Buttons",
                                        options = list(scrollX = TRUE,
                                                       dom = "Bfrtip",
                                                       buttons = c("copy", "csv")
                                        )
                                      )
                                    }, options = list(
                                      dom = "Bfrtip", list(
                                        list(
                                          extend = "csv",
                                          text = "csv",
                                          exportOptions = list(modifier = list(page = "all"))
                                        ),
                                        "copy"
                                      )
                                    ), server=F
      )




      o_UpdateSel <- metaObserve2({
        shiny::req(input$Update_selection)
        isolate(metaExpr({
          # tmp_rows_cols <- AbSolution:::filter_merged(
          #   Big_mem_values$Big_DF, Big_mem_values$Short_DF, Big_mem_values$Header,
          #   "Reconstructed germline" %in% input$use_what, "Repertoire" %in% input$use_what,
          #   "Productive" %in% input$use_productive_or_not, "Non-productive" %in%
          #     input$use_productive_or_not, input$my_regions, input$my_var_elements,
          #   input$my_vars, input$my_vartypes, input$use_sharedVDJ, input$VJ_included,
          #   input$groups_selected, input$group_A, input$group_B, input$group_C,
          #   input$use_univlog, input$samples_selected, input$exclude_variables,
          #   input$pval_type, input$pval_cutoff, input$estimate_cutoff, input$number_selected_vars,
          #   input$VJ_deselected, input$VDJ_normalized_per_size, input$Rmut_filter[1],
          #   input$Rmut_filter[2], input$work_as_categories, input$VDJ_maximize_clones,
          #   input$VDJ_normalized_per_sample, input$my_clone_def, input$seed,
          #   input$chains_selected, input$subtype_selected
          # )

          ok <- AbSolution::validate_before_run(work_as_categories=input$work_as_categories,
                                                use_sharedVDJ=input$use_sharedVDJ,
                                                VJ_included=input$VJ_included,
                                                VDJ_maximize_clones=input$VDJ_maximize_clones,
                                                clonal_group=input$clonal_group)
          shiny::req(ok)

          tmp_rows_cols <- AbSolution::filter_merged(
            FBM=Big_mem_values$Big_DF,
            merged_df=Big_mem_values$Short_DF,
            merged_header=Big_mem_values$Header,
            use_rgermline="Reconstructed germline" %in% input$use_what,
            use_repertoire="Repertoire" %in% input$use_what,
            use_productive="Productive" %in% input$use_productive_or_not,
            use_nonproductive= "Non-productive" %in% input$use_productive_or_not,
            my_regions= input$my_regions,
            my_var_elements=input$my_var_elements,
            my_vars=input$my_vars,
            my_vartypes=input$my_vartypes,
            use_sharedVDJ=input$use_sharedVDJ,
            V_J_to_use=input$VJ_included,
            groups=input$groups_selected,
            group_A=input$group_A,
            group_B=input$group_B,
            group_C=input$group_C,
            univlog=input$use_univlog,
            samples_to_keep=input$samples_selected,
            variables_to_remove=input$exclude_variables,
            pval_type=input$pval_type,
            pval_cutoff=input$pval_cutoff,
            estimate_cutoff=input$estimate_cutoff,
            number_selected_vars=input$number_selected_vars,
            VJ_deselected=input$VJ_deselected,
            VDJ_normalized_per_size=input$VDJ_normalized_per_size,
            R_mut_threshold_min=input$Rmut_filter[1],
            R_mut_threshold_max=input$Rmut_filter[2],
            to_compare_groups=input$work_as_categories,
            VDJ_maximize_clones=input$VDJ_maximize_clones,
            VDJ_normalized_per_sample=input$VDJ_normalized_per_sample,
            my_clone_def=input$clonal_group,
            seed=input$seed,
            chains=input$chains_selected,
            igsubtypes=input$subtype_selected
          )
          Selection_values$rows <- tmp_rows_cols$ROWS
          Selection_values$columns <- tmp_rows_cols$COLUMNS
          Selection_values$Parameters <- list(input$use_UMAP)
        }))
      })

      # observeEvent(
      #     input$Update_selection, {
      #         tmp_rows_cols <- AbSolution:::filter_merged(
      #             Big_mem_values$Big_DF, Big_mem_values$Short_DF, Big_mem_values$Header,
      #             "Reconstructed germline" %in% input$use_what, "Repertoire" %in% input$use_what,
      #             "Productive" %in% input$use_productive_or_not, "Non-productive" %in%
      #               input$use_productive_or_not, input$my_regions, input$my_var_elements,
      #             input$my_vars, input$my_vartypes, input$use_sharedVDJ, input$VJ_included,
      #             input$groups_selected, input$group_A, input$group_B, input$group_C,
      #             input$use_univlog, input$samples_selected, input$exclude_variables,
      #             input$pval_type, input$pval_cutoff, input$estimate_cutoff, input$number_selected_vars,
      #             input$VJ_deselected, input$VDJ_normalized_per_size, input$Rmut_filter[1],
      #             input$Rmut_filter[2], input$work_as_categories, input$VDJ_maximize_clones,
      #             input$VDJ_normalized_per_sample, input$my_clone_def
      #         )
      #         Selection_values$rows <- tmp_rows_cols$ROWS
      #         Selection_values$columns <- tmp_rows_cols$COLUMNS
      #         Selection_values$Parameters <- list(input$use_UMAP)
      #     }
      # )

      o_UpdateExpl <- metaObserve2({
        shiny::req(input$Update_exploration)
        isolate(metaExpr({

          ok <- AbSolution::validate_before_run(work_as_categories=input$work_as_categories,
                                                use_sharedVDJ=input$use_sharedVDJ,
                                                VJ_included=input$VJ_included,
                                                VDJ_maximize_clones=input$VDJ_maximize_clones,
                                                clonal_group=input$clonal_group)
          shiny::req(ok)

          tmp_rows_cols <- AbSolution::filter_merged(
            FBM=Big_mem_values$Big_DF,
            merged_df=Big_mem_values$Short_DF,
            merged_header=Big_mem_values$Header,
            use_rgermline="Reconstructed germline" %in% input$use_what,
            use_repertoire="Repertoire" %in% input$use_what,
            use_productive="Productive" %in% input$use_productive_or_not,
            use_nonproductive= "Non-productive" %in% input$use_productive_or_not,
            my_regions= input$my_regions,
            my_var_elements=input$my_var_elements,
            my_vars=input$my_vars,
            my_vartypes=input$my_vartypes,
            use_sharedVDJ=input$use_sharedVDJ,
            V_J_to_use=input$VJ_included,
            groups=input$groups_selected,
            group_A=input$group_A,
            group_B=input$group_B,
            group_C=input$group_C,
            univlog=input$use_univlog,
            samples_to_keep=input$samples_selected,
            variables_to_remove=input$exclude_variables,
            pval_type=input$pval_type,
            pval_cutoff=input$pval_cutoff,
            estimate_cutoff=input$estimate_cutoff,
            number_selected_vars=input$number_selected_vars,
            VJ_deselected=input$VJ_deselected,
            VDJ_normalized_per_size=input$VDJ_normalized_per_size,
            R_mut_threshold_min=input$Rmut_filter[1],
            R_mut_threshold_max=input$Rmut_filter[2],
            to_compare_groups=input$work_as_categories,
            VDJ_maximize_clones=input$VDJ_maximize_clones,
            VDJ_normalized_per_sample=input$VDJ_normalized_per_sample,
            my_clone_def=input$clonal_group,
            seed=input$seed,
            chains=input$chains_selected,
            igsubtypes=input$subtype_selected
          )
          Exploration_values$rows <- tmp_rows_cols$ROWS
          Exploration_values$columns <- tmp_rows_cols$COLUMNS
          Exploration_values$Parameters <- list(input$use_UMAP)
        }))
      })

      # observeEvent(
      #     input$Update_exploration, {
      #         tmp_rows_cols <- AbSolution:::filter_merged(
      #             Big_mem_values$Big_DF, Big_mem_values$Short_DF, Big_mem_values$Header,
      #             "Reconstructed germline" %in% input$use_what, "Repertoire" %in% input$use_what,
      #             "Productive" %in% input$use_productive_or_not, "Non-productive" %in%
      #               input$use_productive_or_not, input$my_regions, input$my_var_elements,
      #             input$my_vars, input$my_vartypes, input$use_sharedVDJ, input$VJ_included,
      #             input$groups_selected, input$group_A, input$group_B, input$group_C,
      #             input$use_univlog, input$samples_selected, input$exclude_variables,
      #             input$pval_type, input$pval_cutoff, input$estimate_cutoff, input$number_selected_vars,
      #             input$VJ_deselected, input$VDJ_normalized_per_size, input$Rmut_filter[1],
      #             input$Rmut_filter[2], input$work_as_categories, input$VDJ_maximize_clones,
      #             input$VDJ_normalized_per_sample, input$my_clone_def
      #         )
      #         Exploration_values$rows <- tmp_rows_cols$ROWS
      #         Exploration_values$columns <- tmp_rows_cols$COLUMNS
      #         Exploration_values$Parameters <- list(input$use_UMAP)
      #     }
      # )






      output$Group_selection <- renderUI(
        {
          choices <- c("Group","Subgroup", "Patient","Sample","Patient_Sample")
          if (!is.null(Big_mem_values$Short_DF)) {
            if ("Dominance" %in% colnames(Big_mem_values$Short_DF)) {
              choices <- c(choices, "Dominance")
            }
            if ("Selected" %in% colnames(Big_mem_values$Short_DF)) {
              choices <- c(choices, "Selected")
            }

            if ("Selected_clones" %in% colnames(Big_mem_values$Short_DF)) {
              choices <- c(choices, "Selected_clones")
            }
            if ("C_region" %in% colnames(Big_mem_values$Short_DF)) {
              choices <- c(choices, "C_region")
            }
            if ("Chain" %in% colnames(Big_mem_values$Short_DF)) {
              choices <- c(choices, "Chain")
            }
            if (length(Big_mem_values$categorical_newfields)>0) {
              choices <- c(choices, Big_mem_values$categorical_newfields)
            }
          }
          selectInput("groups_selected", "Work only with", choices = choices, selected = "Group")
        }
      )
      # observeEvent(input$Move_to_analysis_real, {
      # info=merge_bigmems(folder_values$Featured)
      # Big_mem_values$Header=info[[1]] Big_mem_values$Short_DF=info[[2]]
      # Big_mem_values$Big_DF=info[[3]] rm(info)
      # Exploration_values$rows=c(1:nrow(Big_mem_values$Big_DF))
      # Exploration_values$columns=c(1:ncol(Big_mem_values$Big_DF))
      # Selection_values$rows=c(1:nrow(Big_mem_values$Big_DF))
      # Selection_values$columns=c(1:ncol(Big_mem_values$Big_DF)) })


      # in server.R create reactiveVal
      Sample_selection_current_selection <- reactiveVal(NULL)



      o_ignore_Init_samples_Rmut <- metaObserve2({
        shiny::req(list(input$samples_selected, input$Rmut_filter))
        isolate(metaExpr({

          if(sample_info_react$step_3) {
            filter_muts <- intersect(
              which(
                Big_mem_values$Big_DF[, which(Big_mem_values$Header == "AA_Whole_Replacement_muts_counts")] >=
                  input$Rmut_filter[1]
              ),
              which(
                Big_mem_values$Big_DF[, which(Big_mem_values$Header == "AA_Whole_Replacement_muts_counts")] <=
                  input$Rmut_filter[2]
              )
            )
            Big_mem_values$VJs <- sort(
              unique(
                Big_mem_values$Short_DF$V_and_J[intersect(
                  intersect(
                    which(Big_mem_values$Short_DF$Patient_Sample %in% input$samples_selected),
                    which(Big_mem_values$Short_DF$Sequence_type %in% input$use_what)
                  ),
                  filter_muts
                )]
              )
            )
            counts_VJs <- table(
              Big_mem_values$Short_DF$V_and_J[intersect(
                intersect(
                  which(Big_mem_values$Short_DF$Patient_Sample %in% input$samples_selected),
                  which(Big_mem_values$Short_DF$Sequence_type %in% input$use_what)
                ),
                filter_muts
              )]
            )

            # counts_VJs=table(Big_mem_values$Short_DF$V_and_J[(intersect(which(Big_mem_values$Short_DF$Patient_Sample
            # %in% input$samples_selected),
            # which(Big_mem_values$Short_DF$Sequence_type %in% input$use_what))
            # )])
            names(Big_mem_values$VJs) <- sapply(
              Big_mem_values$VJs, function(z) paste(
                z, counts_VJs[which(
                  names(counts_VJs) ==
                    z
                )],
                sep = " - Counts:"
              )
            )
            Sample_selection_current_selection(input$samples_selected)
          }

        }))
      })

      output$Sample_selection <- renderUI(
        {
          choices <- c("RANDOMLETTERS")
          if (!is.null(Big_mem_values$Patient_Sample)) {
            choices <- Big_mem_values$Patient_Sample
          }
          if (!is.null(Sample_selection_current_selection()) &&
              Sample_selection_current_selection()[1] != "") {
            tmp_selection <- Sample_selection_current_selection()
          } else {
            tmp_selection <- choices
          }
          if (length(tmp_selection) ==
              0) {
            tmp_selection <- choices[1]
          }
          selectizeInput(
            "samples_selected", "Select the sample(s)", choices = choices, selected = tmp_selection,
            multiple = TRUE, options = list(plugins = list("remove_button"))
          )
        }
      )


      filterselectionValues <- reactiveValues(chain_choices=NULL)
      output$Chain_selection <- renderUI(
        {
          choices <- c("Unknown")
          if (!is.null(Big_mem_values$Short_DF)) {
            if("Chain" %in% colnames(Big_mem_values$Short_DF)) {
              choices <- unique(Big_mem_values$Short_DF$Chain[which(Big_mem_values$Short_DF$Patient_Sample %in% input$samples_selected)])
            }


          }
          filterselectionValues$chain_choices=choices
          selectizeInput(
            "chains_selected", "Select the chain(s)", choices = choices, selected = choices,
            multiple = TRUE, options = list(plugins = list("remove_button"))
          )

        }
      )

      observeEvent(input$chains_selected, {

        if(any(is.na(input$chains_selected)) || any(is.null(input$chains_selected)) || length(input$chains_selected) == 0 || identical(input$chains_selected, integer(0))){

          updateSelectizeInput(session,"chains_selected",choices=filterselectionValues$chain_choices,selected=filterselectionValues$chain_choices[1])

        }
      })
      output$Type_selection <- renderUI(
        {
          choices <- c("Unknown")
          if (!is.null(Big_mem_values$Short_DF)) {
            if("C_region" %in% colnames(Big_mem_values$Short_DF)) {
              if("Unknown" %in% input$chains_selected) {
                choices <- unique(Big_mem_values$Short_DF$C_region[which(Big_mem_values$Short_DF$Patient_Sample %in% input$samples_selected)])
              } else {
                choices <-unique( Big_mem_values$Short_DF$C_region[intersect(which(Big_mem_values$Short_DF$Patient_Sample %in% input$samples_selected),
                                                                             which(Big_mem_values$Short_DF$Chain %in% input$chains_selected))])
              }

            }

          }

          selectizeInput(
            "subtype_selected", "Select the Ig isotype(s)", choices = choices, selected = choices,
            multiple = TRUE, options = list(plugins = list("remove_button"))
          )
        }
      )

      # in server.R create reactiveVal
      VJ_selection_current_selection <- reactiveVal(NULL)

      # now store your current selection in the reactive value
      o_VJ_deselected <- metaObserve2({
        shiny::req(input$VJ_deselected)
        isolate(metaExpr({
          VJ_selection_current_selection(input$VJ_deselected)
        }))
      })

      # observeEvent(
      #     input$VJ_deselected, {
      #         VJ_selection_current_selection(input$VJ_deselected)
      #     }
      # )



      output$VJ_selection <- renderUI(
        {
          choices <- c("RANDOMLETTERS")
          if (!is.null(Big_mem_values$VJs)) {
            choices <- Big_mem_values$VJs
          }
          if (!is.null(VJ_selection_current_selection()) &&
              VJ_selection_current_selection()[1] != "") {
            tmp_selection <- VJ_selection_current_selection()
            tmp_selection <- tmp_selection[which(tmp_selection %in% Big_mem_values$VJs)]
          } else {
            tmp_selection <- NULL
          }

          # selectizeInput( 'VJ_selected', 'Include combination(s)', choices
          # = choices, selected = tmp_selection, multiple = TRUE, options =
          # list(plugins= list('remove_button')))
          if (length(tmp_selection) ==
              length(choices)) {
            tmp_selection <- choices[2:length(choices)]
          }
          shinyWidgets::pickerInput(
            inputId = "VJ_deselected", label = "Exclude the following VJ combinations:",
            choices = choices, selected = tmp_selection, options = shinyWidgets::pickerOptions(
              actionsBox = TRUE, liveSearch = TRUE, maxOptions = (length(choices)) -
                1, size = 10, selectedTextFormat = "count > 3", selectAllText = "Select All but 1"
            ),
            multiple = TRUE
          )


        }
      )

      # o_ShortDF <- metaObserve2({
      #   req(Big_mem_values$Short_DF)
      #   isolate(metaExpr({
      #     if(sample_info_react$step_3) {
      #       if (!is.null(Big_mem_values$Short_DF) &&
      #           any(
      #             startsWith(
      #               colnames(Big_mem_values$Short_DF),
      #               "Clone_"
      #             )
      #           )) {
      #         updateSelectizeInput(
      #           session, "my_clone_def", choices = colnames(Big_mem_values$Short_DF)[which(
      #             startsWith(
      #               colnames(Big_mem_values$Short_DF),
      #               "Clone_"
      #             )
      #           )]
      #         )
      #       }
      #     }
      #
      #
      #
      #   }))
      # })

      output$my_clone_def <- renderUI({
        if (is.null(input$clonal_group) || input$clonal_group == "") {
          tags$em("No clonal definition is currently selected")
        } else {
          tagList(
            tags$b("Currently using clonal definition: "),
            tags$i(input$clonal_group)
          )
        }
      })



      output$Group_comparison <- renderUI(
        {

          labels <- sort(
            unique(
              Big_mem_values$Short_DF[[input$groups_selected]][which(Big_mem_values$Short_DF$Patient_Sample %in% (input$samples_selected))]
            )
          )
          bucket_list(
            header = "Include at least one element on each group", orientation = "horizontal",
            group_name = "group_classification", add_rank_list(text = "Group A", labels = labels, input_id = "group_A"),
            add_rank_list(text = "Group B", labels = NULL, input_id = "group_B"),
            add_rank_list(text = "Ignore", labels = NULL, input_id = "group_C")
          )

        }
      )


      o_SamSelWhatProd <- metaObserve2({
        shiny::req(list(input$Sample_selection, input$use_what, input$use_productive_or_not,
                 input$subtype_selected, input$chains_selected))
        isolate(metaExpr({
          {

            if(sample_info_react$step_3) {

              rows <- c()
              if ("Repertoire" %in% input$use_what) {
                rows <- c(rows, which(Big_mem_values$Short_DF$Sequence_type == "Repertoire"))
              }
              if ("Reconstructed germline" %in% input$use_what) {
                rows <- c(rows, which(Big_mem_values$Short_DF$Sequence_type != "Repertoire"))
              }

              rows <- intersect(
                rows, which(Big_mem_values$Short_DF$Patient_Sample %in% (input$samples_selected))
              )
              if("Chain" %in% colnames(Big_mem_values$Short_DF))  {
                rows <- intersect(rows,
                                  which(Big_mem_values$Short_DF$Chain %in% (input$chains_selected))
                )
              }
              if("C_region" %in% colnames(Big_mem_values$Short_DF))  {
                rows <- intersect(rows,
                                  which(Big_mem_values$Short_DF$C_region %in% (input$subtype_selected))
                )
              }

              vals <- Big_mem_values$Big_DF[rows, which(Big_mem_values$Header == "AA_Whole_Replacement_muts_counts")]
              if (length(vals) > 0 && any(is.finite(vals))) {
                updateSliderInput(inputId = "Rmut_filter", min = min(vals, na.rm = TRUE), max = max(vals, na.rm = TRUE))
              }
            }


          }
        }))
      })


      shared_VJ_list <- reactiveVal(NULL)

      o_longlist <- metaObserve2({
        shiny::req(list(
          input$samples_selected, input$use_what, input$groups_selected, input$VJ_deselected,
          input$use_sharedVDJ, input$VDJ_normalized_per_size, input$clonal_group,
          input$Rmut_filter, input$VDJ_normalized_per_sample, input$VDJ_maximize_clones
        ))
        isolate(metaExpr({
          {
            if(sample_info_react$step_3) {
              rows <- c()

              if (input$use_sharedVDJ) {
                if ("Repertoire" %in% input$use_what) {
                  rows <- c(rows, which(Big_mem_values$Short_DF$Sequence_type == "Repertoire"))
                }
                if ("Reconstructed germline" %in% input$use_what) {
                  rows <- c(rows, which(Big_mem_values$Short_DF$Sequence_type == "Reconstructed_germline"))
                }
                rows <- sort(rows)
                intersect_mut <- intersect(
                  which(
                    Big_mem_values$Big_DF[, which(Big_mem_values$Header == "AA_Whole_Replacement_muts_counts")] >=
                      input$Rmut_filter[1]
                  ),
                  which(
                    Big_mem_values$Big_DF[, which(Big_mem_values$Header == "AA_Whole_Replacement_muts_counts")] <=
                      input$Rmut_filter[2]
                  )
                )

                intersect_groupA_samples <- intersect(
                  which(Big_mem_values$Short_DF[[input$groups_selected]] %in% input$group_A),
                  which(Big_mem_values$Short_DF$Patient_Sample %in% (input$samples_selected))
                )
                intersect_all_A <- intersect(
                  intersect(intersect_groupA_samples, rows),
                  intersect_mut
                )
                int_group_A <- Big_mem_values$Short_DF$V_and_J[intersect_all_A]
                table_int_group_A <- table(int_group_A)
                intersect_groupB_samples <- intersect(
                  which(Big_mem_values$Short_DF[[input$groups_selected]] %in% input$group_B),
                  which(Big_mem_values$Short_DF$Patient_Sample %in% (input$samples_selected))
                )
                intersect_all_B <- intersect(
                  intersect(intersect_groupB_samples, rows),
                  intersect_mut
                )
                int_group_B <- Big_mem_values$Short_DF$V_and_J[intersect_all_B]
                table_int_group_B <- table(int_group_B)

                labels <- sort(unique(intersect(int_group_A, int_group_B)))
                # labels=intersect(labels, input$VJ_deselected)
                labels <- labels[which(labels %!in% input$VJ_deselected)]
                samples_A <- (unique(Big_mem_values$Short_DF$Patient_Sample[intersect_groupA_samples]))
                samples_B <- (unique(Big_mem_values$Short_DF$Patient_Sample[intersect_groupB_samples]))
                if (input$VDJ_normalized_per_size && input$VDJ_normalized_per_sample) {

                  min_sizes_A <- sapply(
                    labels, function(z) min(
                      sapply(
                        samples_A, function(y) length(
                          which(
                            Big_mem_values$Short_DF$V_and_J[intersect(intersect_all_A, which(Big_mem_values$Short_DF$Patient_Sample %in% y))] ==
                              z
                          )
                        )
                      )
                    )
                  )
                  min_sizes_B <- sapply(
                    labels, function(z) min(
                      sapply(
                        samples_B, function(y) length(
                          which(
                            Big_mem_values$Short_DF$V_and_J[intersect(intersect_all_B, which(Big_mem_values$Short_DF$Patient_Sample %in% y))] ==
                              z
                          )
                        )
                      )
                    )
                  )
                  min_total_sizes <- sapply(
                    1:length(labels),
                    function(z) min(
                      min_sizes_A[z] * length(samples_A),
                      min_sizes_B[z] * length(samples_B)
                    )
                  )
                  min_total_sizes_per_sample_groupA <- floor(min_total_sizes/length(samples_A))
                  min_total_sizes_per_sample_groupB <- floor(min_total_sizes/length(samples_B))

                  labels <- labels[which(min_total_sizes != 0)]
                  min_total_sizes_per_sample_groupA <- min_total_sizes_per_sample_groupA[which(min_total_sizes != 0)]
                  min_total_sizes_per_sample_groupB <- min_total_sizes_per_sample_groupB[which(min_total_sizes != 0)]
                  min_total_sizes <- min_total_sizes[which(min_total_sizes != 0)]
                  min_sizes_B <- min_sizes_B[which(min_total_sizes != 0)]
                  min_sizes_A <- min_sizes_A[which(min_total_sizes != 0)]
                  if (!is.null(input$clonal_group) &&
                      input$VDJ_maximize_clones) {

                    names(labels) <- sapply(
                      c(1:length(labels)),
                      function(z) paste(
                        labels[z], paste(
                          min_total_sizes[z], paste(
                            " (", paste(
                              sapply(
                                c(1:length(samples_A)),
                                function(y) min(
                                  min_total_sizes_per_sample_groupA[z], length(
                                    unique(
                                      Big_mem_values$Short_DF[[input$clonal_group]][intersect(
                                        intersect(
                                          intersect_all_A, which(
                                            Big_mem_values$Short_DF$Patient_Sample %in%
                                              samples_A[y]
                                          )
                                        ),
                                        which(Big_mem_values$Short_DF$V_and_J == labels[z])
                                      )]
                                    )
                                  )
                                )
                              ),
                              collapse = ", "
                            ),
                            " clones on each sample)", sep = ""
                          ),
                          " & group B: ", min_total_sizes[z], paste(
                            " (", paste(
                              sapply(
                                c(1:length(samples_B)),
                                function(y) min(
                                  min_total_sizes_per_sample_groupB[z], length(
                                    unique(
                                      Big_mem_values$Short_DF[[input$clonal_group]][intersect(
                                        intersect(
                                          intersect_all_B, which(
                                            Big_mem_values$Short_DF$Patient_Sample %in%
                                              samples_B[y]
                                          )
                                        ),
                                        which(Big_mem_values$Short_DF$V_and_J == labels[z])
                                      )]
                                    )
                                  )
                                )
                              ),
                              collapse = ", "
                            ),
                            " clones on each sample)", sep = ""
                          ),
                          sep = ""
                        ),
                        sep = " - Sequence counts: group A: "
                      )
                    )

                  } else {
                    names(labels) <- paste(
                      labels, " - Counts: group A - ", min_total_sizes, " (", min_total_sizes_per_sample_groupA,
                      " per sample) and ", min_total_sizes - min_total_sizes_per_sample_groupA *
                        length(samples_A),
                      " random & group B - ", min_total_sizes, " (", min_total_sizes_per_sample_groupB,
                      " per sample) and ", min_total_sizes - min_total_sizes_per_sample_groupB *
                        length(samples_B),
                      " random", sep = ""
                    )
                  }

                } else if (input$VDJ_normalized_per_size) {

                  if (!is.null(input$clonal_group) &&
                      input$VDJ_maximize_clones) {

                    names(labels) <- sapply(
                      labels, function(z) paste(
                        z, paste(
                          min(
                            table_int_group_A[which(
                              names(table_int_group_A) ==
                                z
                            )],
                            table_int_group_B[which(
                              names(table_int_group_B) ==
                                z
                            )]
                          ),
                          paste(
                            " (", paste(
                              sapply(
                                samples_A, function(y) min(
                                  min(
                                    table_int_group_A[which(
                                      names(table_int_group_A) ==
                                        z
                                    )],
                                    table_int_group_B[which(
                                      names(table_int_group_B) ==
                                        z
                                    )]
                                  ),
                                  length(
                                    unique(
                                      Big_mem_values$Short_DF[[input$clonal_group]][intersect(
                                        intersect(
                                          intersect_all_A, which(
                                            Big_mem_values$Short_DF$Patient_Sample %in%
                                              y
                                          )
                                        ),
                                        which(Big_mem_values$Short_DF$V_and_J == z)
                                      )]
                                    )
                                  )
                                )
                              ),
                              collapse = ", "
                            ),
                            " clones on each sample)", sep = ""
                          ),
                          " & group B: ", min(
                            table_int_group_A[which(
                              names(table_int_group_A) ==
                                z
                            )],
                            table_int_group_B[which(
                              names(table_int_group_B) ==
                                z
                            )]
                          ),
                          paste(
                            " (", paste(
                              sapply(
                                samples_B, function(y) min(
                                  min(
                                    table_int_group_A[which(
                                      names(table_int_group_A) ==
                                        z
                                    )],
                                    table_int_group_B[which(
                                      names(table_int_group_B) ==
                                        z
                                    )]
                                  ),
                                  length(
                                    unique(
                                      Big_mem_values$Short_DF[[input$clonal_group]][intersect(
                                        intersect(
                                          intersect_all_B, which(
                                            Big_mem_values$Short_DF$Patient_Sample %in%
                                              y
                                          )
                                        ),
                                        which(Big_mem_values$Short_DF$V_and_J == z)
                                      )]
                                    )
                                  )
                                )
                              ),
                              collapse = ", "
                            ),
                            " clones on each sample)", sep = ""
                          ),
                          sep = ""
                        ),
                        sep = " - Sequence counts: group A: "
                      )
                    )

                  } else {
                    names(labels) <- sapply(
                      labels, function(z) paste(
                        z, paste(
                          min(
                            table_int_group_A[which(
                              names(table_int_group_A) ==
                                z
                            )],
                            table_int_group_B[which(
                              names(table_int_group_B) ==
                                z
                            )]
                          ),
                          min(
                            table_int_group_A[which(
                              names(table_int_group_A) ==
                                z
                            )],
                            table_int_group_B[which(
                              names(table_int_group_B) ==
                                z
                            )]
                          ),
                          sep = " & group B - "
                        ),
                        sep = " - Counts: group A -"
                      )
                    )

                  }



                } else {

                  if (!is.null(input$clonal_group) &&
                      input$VDJ_maximize_clones) {
                    names(labels) <- sapply(
                      labels, function(z) paste(
                        z, paste(
                          table_int_group_A[which(
                            names(table_int_group_A) ==
                              z
                          )],
                          paste(
                            " (", paste(
                              sapply(
                                samples_A, function(y) length(
                                  unique(
                                    Big_mem_values$Short_DF[[input$clonal_group]][intersect(
                                      intersect(
                                        intersect_all_A, which(
                                          Big_mem_values$Short_DF$Patient_Sample %in%
                                            y
                                        )
                                      ),
                                      which(Big_mem_values$Short_DF$V_and_J == z)
                                    )]
                                  )
                                )
                              ),
                              collapse = ", "
                            ),
                            " clones on each sample)", sep = ""
                          ),
                          " & group B - ", table_int_group_B[which(
                            names(table_int_group_B) ==
                              z
                          )],
                          paste(
                            " (", paste(
                              sapply(
                                samples_B, function(y) length(
                                  unique(
                                    Big_mem_values$Short_DF[[input$clonal_group]][intersect(
                                      intersect(
                                        intersect_all_B, which(
                                          Big_mem_values$Short_DF$Patient_Sample %in%
                                            y
                                        )
                                      ),
                                      which(Big_mem_values$Short_DF$V_and_J == z)
                                    )]
                                  )
                                )
                              ),
                              collapse = ", "
                            ),
                            " clones on each sample)", sep = ""
                          ),
                          sep = ""
                        ),
                        sep = " - Sequence counts: group A -"
                      )
                    )

                  } else {
                    names(labels) <- sapply(
                      labels, function(z) paste(
                        z, paste(
                          table_int_group_A[which(
                            names(table_int_group_A) ==
                              z
                          )],
                          table_int_group_B[which(
                            names(table_int_group_B) ==
                              z
                          )],
                          sep = " & group B - "
                        ),
                        sep = " - Counts: group A -"
                      )
                    )

                  }
                }
                shared_VJ_list(labels)
              } else {
                shared_VJ_list(NULL)
              }

            }
          }

        }))
      })


      output$VDJ_subsetting <- renderUI(
        {
          labels <- shared_VJ_list()
          shiny::req(!is.null(labels), length(labels) > 0)

          shinyWidgets::pickerInput(
            inputId = "VJ_included", label = "To include only the following VJ combinations:",
            choices = labels, options = shinyWidgets::pickerOptions(
              actionsBox = TRUE, liveSearch = TRUE, size = 10, selectedTextFormat = "count > 3"
            ),
            multiple = TRUE
          )


        }
      )




      individual_variables_current_list <- reactiveVal(NULL)

      # now store your current selection in the reactive value

      o_list_parseseq <- metaObserve2({
        shiny::req(list(
          input$samples_selected, input$my_vars, input$use_what, input$use_productive_or_not,
          input$my_regions, input$my_var_elements, input$my_vartypes
        ))
        isolate(metaExpr({

          if(sample_info_react$step_3){
            choices <- list()

            if (!is.null(input$my_vars)) {
              columns <- c(1:length(Big_mem_values$Header))
              columns <- columns[which(
                columns %in% which(
                  (sapply(Big_mem_values$Header, function(x) strsplit(x, split = "_")[[1]][2])) %in%
                    input$my_regions
                )
              )]
              columns <- columns[which(
                columns %in% which(
                  grepl(
                    paste0(
                      paste("^", input$my_var_elements, sep = ""),
                      collapse = "|"
                    ),
                    Big_mem_values$Header
                  )
                )
              )]
              if ("Germline diff" %!in% input$my_vartypes) {
                columns <- columns[which(columns %!in% which(endsWith(Big_mem_values$Header, "_Diff")))]
              }
              if ("Baseline" %!in% input$my_vartypes) {
                columns <- columns[which(columns %in% which(endsWith(Big_mem_values$Header, "_Diff")))]
              }

              merged_header_splitted_at_3 <- sapply(Big_mem_values$Header, function(x) strsplit(x, split = "_")[[1]][3])
              merged_header_splitted_at_4 <- sapply(Big_mem_values$Header, function(x) strsplit(x, split = "_")[[1]][4])
              merged_header_splitted_at_3_and_4 <- sapply(
                Big_mem_values$Header, function(x) paste(
                  strsplit(x, split = "_")[[1]][3],
                  strsplit(x, split = "_")[[1]][4],
                  sep = "_"
                )
              )


              if ("Length" %in% input$my_vars) {
                tmp_index_var <- c(
                  columns[which(columns %in% which(endsWith(Big_mem_values$Header, "length")))],
                  columns[which(columns %in% which(endsWith(Big_mem_values$Header, "length_Diff")))]
                )
                choices[["Length"]] <- as.list(sort(Big_mem_values$Header[tmp_index_var]))
              }

              if ("Composition" %in% input$my_vars) {
                nucleotides <- c("A", "G", "T", "C")
                peptides <- Peptides::aaList()

                tmp_index_var <- c(
                  columns[which(
                    columns %in% which(
                      merged_header_splitted_at_3_and_4 %in% unique(
                        c(
                          paste(
                            c(nucleotides, peptides),
                            "count", sep = "_"
                          ),
                          paste(
                            c(nucleotides, peptides),
                            "norm", sep = "_"
                          )
                        )
                      )
                    )
                  )],
                  columns[which(columns %in% which(merged_header_splitted_at_4 %in% c("counts")))]
                )
                choices[["Composition"]] <- as.list(sort(Big_mem_values$Header[tmp_index_var]))
              }

              if ("Hot/Cold motifs" %in% input$my_vars) {
                tmp_index_var <- c(
                  columns[which(
                    columns %in% which(merged_header_splitted_at_3 %in% c("hot", "cold", "potential"))
                  )]
                )
                choices[["Hot/Cold motifs"]] <- as.list(sort(Big_mem_values$Header[tmp_index_var]))
              }

              if ("Substitutions" %in% input$my_vars) {
                tmp_index_var <- c(
                  columns[which(
                    columns %in% which(
                      merged_header_splitted_at_3_and_4 %in% c("Sub", "Sub_prc", "SIDT_sum", "SID_sum")
                    )
                  )],
                  columns[which(columns %in% which(merged_header_splitted_at_3 %in% c("Sub")))]
                )
                choices[["Substitutions"]] <- as.list(sort(Big_mem_values$Header[tmp_index_var]))
              }

              if ("Insertions" %in% input$my_vars) {
                tmp_index_var <- c(
                  columns[which(
                    columns %in% which(
                      merged_header_splitted_at_3_and_4 %in% c("Ins", "Ins_prc", "SIDT_sum", "SID_sum")
                    )
                  )],
                  columns[which(columns %in% which(merged_header_splitted_at_3 %in% c("Ins")))]
                )
                choices[["Insertions"]] <- as.list(sort(Big_mem_values$Header[tmp_index_var]))
              }
              if ("Deletions" %in% input$my_vars) {
                tmp_index_var <- c(
                  columns[which(
                    columns %in% which(
                      merged_header_splitted_at_3_and_4 %in% c("Del", "Del_prc", "SIDT_sum", "SID_sum")
                    )
                  )],
                  columns[which(columns %in% which(merged_header_splitted_at_3 %in% c("Del")))]
                )
                choices[["Deletions"]] <- as.list(sort(Big_mem_values$Header[tmp_index_var]))
              }
              if ("Translocations" %in% input$my_vars) {

                tmp_index_var <- c(
                  columns[which(
                    columns %in% which(merged_header_splitted_at_3_and_4 %in% c("Trasl", "Trasl_prc", "SIDT_sum"))
                  )],
                  columns[which(columns %in% which(merged_header_splitted_at_3 %in% c("Trasl")))]
                )
                choices[["Translocations"]] <- as.list(sort(Big_mem_values$Header[tmp_index_var]))
              }

              if ("Leveshtein distance" %in% input$my_vars) {
                tmp_index_var <- c(columns[which(columns %in% which(merged_header_splitted_at_3 %in% c("lv")))])
                choices[["Leveshtein distance"]] <- as.list(sort(Big_mem_values$Header[tmp_index_var]))
              }

              if ("Transitions and transversions" %in% input$my_vars) {
                tmp_index_var <- c(
                  columns[which(
                    columns %in% which(merged_header_splitted_at_3 %in% c("Transitions", "Transversions"))
                  )],
                  columns[which(
                    columns %in% which(merged_header_splitted_at_3_and_4 %in% c("Ratio_Transitions-Transversions"))
                  )]
                )
                choices[["Transitions and transversions"]] <- as.list(sort(Big_mem_values$Header[tmp_index_var]))
              }
              if ("Replacement and silent mutations" %in% input$my_vars) {
                tmp_index_var <- c(
                  columns[which(
                    columns %in% which(merged_header_splitted_at_3 %in% c("Replacement", "Silent"))
                  )],
                  columns[which(
                    columns %in% which(merged_header_splitted_at_3_and_4 %in% c("Ratio_Silent-Replacement"))
                  )]
                )
                choices[["Replacement and silent mutations"]] <- as.list(sort(Big_mem_values$Header[tmp_index_var]))
              }

              if ("Mutations from X to Y" %in% input$my_vars) {
                # columns=columns[which(columns %!in%
                # which(sapply(merged_header, function(x) paste(strsplit(x,
                # split='_')[[1]][3], strsplit(x, split='_')[[1]][4],
                # sep='_')) %in% c('Ratio_Silent') ))]
                tmp_index_var <- c(columns[which(columns %in% which(merged_header_splitted_at_4 %in% c("to")))])
                choices[["Mutations from X to Y"]] <- as.list(sort(Big_mem_values$Header[tmp_index_var]))
              }

              if ("NGly sites" %in% input$my_vars) {
                tmp_index_var <- c(columns[which(columns %in% which(merged_header_splitted_at_3 %in% c("NGly")))])
                choices[["NGly sites"]] <- as.list(sort(Big_mem_values$Header[tmp_index_var]))
              }
              if ("Peptide features" %in% input$my_vars) {
                tmp_index_var <- c(
                  columns[intersect(which(
                    columns %in% which(
                      merged_header_splitted_at_3 %in% c(
                        "Peptides", "alkzm", "Small", "Tiny", "Aliphatic", "Charged",
                        "Polar", "Basic", "NonPolar", "Aromatic", "Acidic"
                      )
                    )
                  ), which(columns %!in% which(merged_header_splitted_at_4 %in% c("to"))))]
                )
                choices[["Peptide features"]] <- as.list(sort(Big_mem_values$Header[tmp_index_var]))
              }


            }

            individual_variables_current_list(choices)
          }

        }))
      })


      output$individual_variables <- renderUI(
        {
          choices <- c("")
          if (!is.null(Big_mem_values$Header) ) {
            choices <- individual_variables_current_list()

          }

          shinyWidgets::virtualSelectInput(
            inputId = "exclude_variables", label = "Remove individual variables:",
            choices = choices, showValueAsTags = TRUE, search = TRUE, multiple = TRUE
          )
        }
      )


      output$AbLogo <- renderUI({
        if(!is.null(input$dark_mode) && input$dark_mode){
          HTML(
            "<img src=\"www/img/AbSolution_logo_darkmode.png\" alt=\"AB Solution\" width=\"653\" height=\"150\">"
          )
        } else {
          HTML(
            "<img src=\"www/img/AbSolution_logo.png\" alt=\"AB Solution\" width=\"653\" height=\"150\">"
          )
        }
      })


      output$RAM_S_Memory_Box <-  bs4Dash::renderInfoBox(
        {
          bs4Dash::infoBox(
            title = "Selection dataset", value = tags$p(
              if (as.numeric(
                unname(
                  strsplit(
                    as.character(utils::capture.output(
                      print(benchmarkme::get_ram(), unit_system = "iec"))
                    ),
                    split = " "
                  )[[1]][1]
                )
              ) >
              round(
                length(Selection_values$rows) *
                length(Selection_values$columns) *
                8/2^{
                  20
                }/1024, 2
              )) {
                "Fits in RAM. UMAP will be available"
              } else {
                "Does not fit in RAM"
              }, style = "font-size: 50%;"
            ),
            icon = icon("info-sign", lib = "glyphicon"),
            color = if (as.numeric(
              unname(
                strsplit(
                  as.character(utils::capture.output(
                    print(benchmarkme::get_ram(), unit_system = "iec"))
                  ),
                  split = " "
                )[[1]][1]
              )
            ) >
            round(
              length(Selection_values$rows) *
              length(Selection_values$columns) *
              8/2^{
                20
              }/1024, 2
            )) {
              "teal"
            } else {
              "red"
            }
          )
        }
      )
      output$RAM_E_Memory_Box <-  bs4Dash::renderInfoBox(
        {
          bs4Dash::infoBox(
            title = "Exploration dataset", value = tags$p(
              if (as.numeric(
                unname(
                  strsplit(
                    as.character(utils::capture.output(
                      print(benchmarkme::get_ram(), unit_system = "iec"))
                    ),
                    split = " "
                  )[[1]][1]
                )
              ) >
              round(
                length(Exploration_values$rows) *
                length(Exploration_values$columns) *
                8/2^{
                  20
                }/1024, 2
              )) {
                "Fits in RAM. UMAP will be available."
              } else {
                "Does not fit in RAM"
              }, style = "font-size: 50%;"
            ),
            icon = icon("info-sign", lib = "glyphicon"),
            color = if (as.numeric(
              unname(
                strsplit(
                  as.character(utils::capture.output(
                    print(benchmarkme::get_ram(), unit_system = "iec"))
                  ),
                  split = " "
                )[[1]][1]
              )
            ) >
            round(
              length(Exploration_values$rows) *
              length(Exploration_values$columns) *
              8/2^{
                20
              }/1024, 2
            )) {
              "teal"
            } else {
              "red"
            }
          )
        }
      )



      output$sunburst_selection_Whole <- renderUI(
        {
          if (!is.null(Selection_values$rows) &&
              length(Selection_values$rows) >
              0 && !is.null(Selection_values$columns) &&
              length(Selection_values$columns) >
              0) {
            tmp_sunburst <- unname(unlist(Big_mem_values$Header[Selection_values$columns]))
            # regions=c('Whole','FWR1','CDR1','FWR2','CDR2','FWR3','CDR3','FWR4')
            region <- "Whole"
            sunburst_df <- show_selected_features(tmp_sunburst, region = region)
            # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel=
            # 'Variables analyzed per each region in the selection set') )

            if (dim(sunburst_df)[1] ==
                0) {
              HTML(
                paste(
                  c("There is no", region, "sequence variables analyzed in the selection set"),
                  collapse = " "
                )
              )
            } else {
              sunburstR::sund2b(
                sunburst_df, showLabels = F, rootLabel = paste(
                  c(region, "sequence variables analyzed in the selection set"),
                  collapse = " "
                )
              )
            }

          }


        }
      )
      output$sunburst_selection_FWR1 <- renderUI(
        {
          if (!is.null(Selection_values$rows) &&
              length(Selection_values$rows) >
              0 && !is.null(Selection_values$columns) &&
              length(Selection_values$columns) >
              0) {
            tmp_sunburst <- unname(unlist(Big_mem_values$Header[Selection_values$columns]))
            # regions=c('Whole','FWR1','CDR1','FWR2','CDR2','FWR3','CDR3','FWR4')
            region <- "FWR1"
            sunburst_df <- show_selected_features(tmp_sunburst, region = region)
            # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel=
            # 'Variables analyzed per each region in the selection set') )
            if (dim(sunburst_df)[1] ==
                0) {
              HTML(
                paste(
                  c("There is no", region, "sequence variables analyzed in the selection set"),
                  collapse = " "
                )
              )
            } else {
              sunburstR::sund2b(
                sunburst_df, showLabels = F, rootLabel = paste(
                  c(region, "sequence variables analyzed in the selection set"),
                  collapse = " "
                )
              )
            }
          }
        }
      )
      output$sunburst_selection_CDR1 <- renderUI(
        {
          if (!is.null(Selection_values$rows) &&
              length(Selection_values$rows) >
              0 && !is.null(Selection_values$columns) &&
              length(Selection_values$columns) >
              0) {
            tmp_sunburst <- unname(unlist(Big_mem_values$Header[Selection_values$columns]))
            # regions=c('Whole','FWR1','CDR1','FWR2','CDR2','FWR3','CDR3','FWR4')
            region <- "CDR1"
            sunburst_df <- show_selected_features(tmp_sunburst, region = region)
            # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel=
            # 'Variables analyzed per each region in the selection set') )
            if (dim(sunburst_df)[1] ==
                0) {
              HTML(
                paste(
                  c("There is no", region, "sequence variables analyzed in the selection set"),
                  collapse = " "
                )
              )
            } else {
              sunburstR::sund2b(
                sunburst_df, showLabels = F, rootLabel = paste(
                  c(region, "sequence variables analyzed in the selection set"),
                  collapse = " "
                )
              )
            }
          }
        }
      )
      output$sunburst_selection_FWR2 <- renderUI(
        {
          if (!is.null(Selection_values$rows) &&
              length(Selection_values$rows) >
              0 && !is.null(Selection_values$columns) &&
              length(Selection_values$columns) >
              0) {
            tmp_sunburst <- unname(unlist(Big_mem_values$Header[Selection_values$columns]))
            # regions=c('Whole','FWR1','CDR1','FWR2','CDR2','FWR3','CDR3','FWR4')
            region <- "FWR2"
            sunburst_df <- show_selected_features(tmp_sunburst, region = region)
            # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel=
            # 'Variables analyzed per each region in the selection set') )
            if (dim(sunburst_df)[1] ==
                0) {
              HTML(
                paste(
                  c("There is no", region, "sequence variables analyzed in the selection set"),
                  collapse = " "
                )
              )
            } else {
              sunburstR::sund2b(
                sunburst_df, showLabels = F, rootLabel = paste(
                  c(region, "sequence variables analyzed in the selection set"),
                  collapse = " "
                )
              )
            }
          }
        }
      )
      output$sunburst_selection_CDR2 <- renderUI(
        {
          if (!is.null(Selection_values$rows) &&
              length(Selection_values$rows) >
              0 && !is.null(Selection_values$columns) &&
              length(Selection_values$columns) >
              0) {
            tmp_sunburst <- unname(unlist(Big_mem_values$Header[Selection_values$columns]))
            # regions=c('Whole','FWR1','CDR1','FWR2','CDR2','FWR3','CDR3','FWR4')
            region <- "CDR2"
            sunburst_df <- show_selected_features(tmp_sunburst, region = region)
            # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel=
            # 'Variables analyzed per each region in the selection set') )
            if (dim(sunburst_df)[1] ==
                0) {
              HTML(
                paste(
                  c("There is no", region, "sequence variables analyzed in the selection set"),
                  collapse = " "
                )
              )
            } else {
              sunburstR::sund2b(
                sunburst_df, showLabels = F, rootLabel = paste(
                  c(region, "sequence variables analyzed in the selection set"),
                  collapse = " "
                )
              )
            }
          }
        }
      )
      output$sunburst_selection_FWR3 <- renderUI(
        {
          if (!is.null(Selection_values$rows) &&
              length(Selection_values$rows) >
              0 && !is.null(Selection_values$columns) &&
              length(Selection_values$columns) >
              0) {
            tmp_sunburst <- unname(unlist(Big_mem_values$Header[Selection_values$columns]))
            # regions=c('Whole','FWR1','CDR1','FWR2','CDR2','FWR3','CDR3','FWR4')
            region <- "FWR3"
            sunburst_df <- show_selected_features(tmp_sunburst, region = region)
            # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel=
            # 'Variables analyzed per each region in the selection set') )
            if (dim(sunburst_df)[1] ==
                0) {
              HTML(
                paste(
                  c("There is no", region, "sequence variables analyzed in the selection set"),
                  collapse = " "
                )
              )
            } else {
              sunburstR::sund2b(
                sunburst_df, showLabels = F, rootLabel = paste(
                  c(region, "sequence variables analyzed in the selection set"),
                  collapse = " "
                )
              )
            }
          }
        }
      )
      output$sunburst_selection_CDR3 <- renderUI(
        {
          if (!is.null(Selection_values$rows) &&
              length(Selection_values$rows) >
              0 && !is.null(Selection_values$columns) &&
              length(Selection_values$columns) >
              0) {
            tmp_sunburst <- unname(unlist(Big_mem_values$Header[Selection_values$columns]))
            # regions=c('Whole','FWR1','CDR1','FWR2','CDR2','FWR3','CDR3','FWR4')
            region <- "CDR3"
            sunburst_df <- show_selected_features(tmp_sunburst, region = region)
            # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel=
            # 'Variables analyzed per each region in the selection set') )
            if (dim(sunburst_df)[1] ==
                0) {
              HTML(
                paste(
                  c("There is no", region, "sequence variables analyzed in the selection set"),
                  collapse = " "
                )
              )
            } else {
              sunburstR::sund2b(
                sunburst_df, showLabels = F, rootLabel = paste(
                  c(region, "sequence variables analyzed in the selection set"),
                  collapse = " "
                )
              )
            }
          }
        }
      )
      output$sunburst_selection_FWR4 <- renderUI(
        {
          if (!is.null(Selection_values$rows) &&
              length(Selection_values$rows) >
              0 && !is.null(Selection_values$columns) &&
              length(Selection_values$columns) >
              0) {
            tmp_sunburst <- unname(unlist(Big_mem_values$Header[Selection_values$columns]))
            # regions=c('Whole','FWR1','CDR1','FWR2','CDR2','FWR3','CDR3','FWR4')
            region <- "FWR4"
            sunburst_df <- show_selected_features(tmp_sunburst, region = region)
            # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel=
            # 'Variables analyzed per each region in the selection set') )
            if (dim(sunburst_df)[1] ==
                0) {
              HTML(
                paste(
                  c("There is no", region, "sequence variables analyzed in the selection set"),
                  collapse = " "
                )
              )
            } else {
              sunburstR::sund2b(
                sunburst_df, showLabels = F, rootLabel = paste(
                  c(region, "sequence variables analyzed in the selection set"),
                  collapse = " "
                )
              )
            }
          }
        }
      )


      output$sunburst_exploration_Whole <- renderUI(
        {
          if (!is.null(Exploration_values$rows) &&
              length(Exploration_values$rows) >
              0 && !is.null(Exploration_values$columns) &&
              length(Exploration_values$columns) >
              0) {
            tmp_sunburst <- unname(unlist(Big_mem_values$Header[Exploration_values$columns]))
            # regions=c('Whole','FWR1','CDR1','FWR2','CDR2','FWR3','CDR3','FWR4')
            region <- "Whole"
            sunburst_df <- show_selected_features(tmp_sunburst, region = region)
            # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel=
            # 'Variables analyzed per each region in the exploration set')
            # )

            if (dim(sunburst_df)[1] ==
                0) {
              HTML(
                paste(
                  c("There is no", region, "sequence variables analyzed in the exploration set"),
                  collapse = " "
                )
              )
            } else {
              sunburstR::sund2b(
                sunburst_df, showLabels = F, rootLabel = paste(
                  c(region, "sequence variables analyzed in the exploration set"),
                  collapse = " "
                )
              )
            }

          }


        }
      )
      output$sunburst_exploration_FWR1 <- renderUI(
        {
          if (!is.null(Exploration_values$rows) &&
              length(Exploration_values$rows) >
              0 && !is.null(Exploration_values$columns) &&
              length(Exploration_values$columns) >
              0) {
            tmp_sunburst <- unname(unlist(Big_mem_values$Header[Exploration_values$columns]))
            # regions=c('Whole','FWR1','CDR1','FWR2','CDR2','FWR3','CDR3','FWR4')
            region <- "FWR1"
            sunburst_df <- show_selected_features(tmp_sunburst, region = region)
            # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel=
            # 'Variables analyzed per each region in the exploration set')
            # )
            if (dim(sunburst_df)[1] ==
                0) {
              HTML(
                paste(
                  c("There is no", region, "sequence variables analyzed in the exploration set"),
                  collapse = " "
                )
              )
            } else {
              sunburstR::sund2b(
                sunburst_df, showLabels = F, rootLabel = paste(
                  c(region, "sequence variables analyzed in the exploration set"),
                  collapse = " "
                )
              )
            }
          }
        }
      )
      output$sunburst_exploration_CDR1 <- renderUI(
        {
          if (!is.null(Exploration_values$rows) &&
              length(Exploration_values$rows) >
              0 && !is.null(Exploration_values$columns) &&
              length(Exploration_values$columns) >
              0) {
            tmp_sunburst <- unname(unlist(Big_mem_values$Header[Exploration_values$columns]))
            # regions=c('Whole','FWR1','CDR1','FWR2','CDR2','FWR3','CDR3','FWR4')
            region <- "CDR1"
            sunburst_df <- show_selected_features(tmp_sunburst, region = region)
            # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel=
            # 'Variables analyzed per each region in the exploration set')
            # )
            if (dim(sunburst_df)[1] ==
                0) {
              HTML(
                paste(
                  c("There is no", region, "sequence variables analyzed in the exploration set"),
                  collapse = " "
                )
              )
            } else {
              sunburstR::sund2b(
                sunburst_df, showLabels = F, rootLabel = paste(
                  c(region, "sequence variables analyzed in the exploration set"),
                  collapse = " "
                )
              )
            }
          }
        }
      )
      output$sunburst_exploration_FWR2 <- renderUI(
        {
          if (!is.null(Exploration_values$rows) &&
              length(Exploration_values$rows) >
              0 && !is.null(Exploration_values$columns) &&
              length(Exploration_values$columns) >
              0) {
            tmp_sunburst <- unname(unlist(Big_mem_values$Header[Exploration_values$columns]))
            # regions=c('Whole','FWR1','CDR1','FWR2','CDR2','FWR3','CDR3','FWR4')
            region <- "FWR2"
            sunburst_df <- show_selected_features(tmp_sunburst, region = region)
            # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel=
            # 'Variables analyzed per each region in the exploration set')
            # )
            if (dim(sunburst_df)[1] ==
                0) {
              HTML(
                paste(
                  c("There is no", region, "sequence variables analyzed in the exploration set"),
                  collapse = " "
                )
              )
            } else {
              sunburstR::sund2b(
                sunburst_df, showLabels = F, rootLabel = paste(
                  c(region, "sequence variables analyzed in the exploration set"),
                  collapse = " "
                )
              )
            }
          }
        }
      )
      output$sunburst_exploration_CDR2 <- renderUI(
        {
          if (!is.null(Exploration_values$rows) &&
              length(Exploration_values$rows) >
              0 && !is.null(Exploration_values$columns) &&
              length(Exploration_values$columns) >
              0) {
            tmp_sunburst <- unname(unlist(Big_mem_values$Header[Exploration_values$columns]))
            # regions=c('Whole','FWR1','CDR1','FWR2','CDR2','FWR3','CDR3','FWR4')
            region <- "CDR2"
            sunburst_df <- show_selected_features(tmp_sunburst, region = region)
            # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel=
            # 'Variables analyzed per each region in the exploration set')
            # )
            if (dim(sunburst_df)[1] ==
                0) {
              HTML(
                paste(
                  c("There is no", region, "sequence variables analyzed in the exploration set"),
                  collapse = " "
                )
              )
            } else {
              sunburstR::sund2b(
                sunburst_df, showLabels = F, rootLabel = paste(
                  c(region, "sequence variables analyzed in the exploration set"),
                  collapse = " "
                )
              )
            }
          }
        }
      )
      output$sunburst_exploration_FWR3 <- renderUI(
        {
          if (!is.null(Exploration_values$rows) &&
              length(Exploration_values$rows) >
              0 && !is.null(Exploration_values$columns) &&
              length(Exploration_values$columns) >
              0) {
            tmp_sunburst <- unname(unlist(Big_mem_values$Header[Exploration_values$columns]))
            # regions=c('Whole','FWR1','CDR1','FWR2','CDR2','FWR3','CDR3','FWR4')
            region <- "FWR3"
            sunburst_df <- show_selected_features(tmp_sunburst, region = region)
            # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel=
            # 'Variables analyzed per each region in the exploration set')
            # )
            if (dim(sunburst_df)[1] ==
                0) {
              HTML(
                paste(
                  c("There is no", region, "sequence variables analyzed in the exploration set"),
                  collapse = " "
                )
              )
            } else {
              sunburstR::sund2b(
                sunburst_df, showLabels = F, rootLabel = paste(
                  c(region, "sequence variables analyzed in the exploration set"),
                  collapse = " "
                )
              )
            }
          }
        }
      )
      output$sunburst_exploration_CDR3 <- renderUI(
        {
          if (!is.null(Exploration_values$rows) &&
              length(Exploration_values$rows) >
              0 && !is.null(Exploration_values$columns) &&
              length(Exploration_values$columns) >
              0) {
            tmp_sunburst <- unname(unlist(Big_mem_values$Header[Exploration_values$columns]))
            # regions=c('Whole','FWR1','CDR1','FWR2','CDR2','FWR3','CDR3','FWR4')
            region <- "CDR3"
            sunburst_df <- show_selected_features(tmp_sunburst, region = region)
            # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel=
            # 'Variables analyzed per each region in the exploration set')
            # )
            if (dim(sunburst_df)[1] ==
                0) {
              HTML(
                paste(
                  c("There is no", region, "sequence variables analyzed in the exploration set"),
                  collapse = " "
                )
              )
            } else {
              sunburstR::sund2b(
                sunburst_df, showLabels = F, rootLabel = paste(
                  c(region, "sequence variables analyzed in the exploration set"),
                  collapse = " "
                )
              )
            }
          }
        }
      )
      output$sunburst_exploration_FWR4 <- renderUI(
        {
          if (!is.null(Exploration_values$rows) &&
              length(Exploration_values$rows) >
              0 && !is.null(Exploration_values$columns) &&
              length(Exploration_values$columns) >
              0) {
            tmp_sunburst <- unname(unlist(Big_mem_values$Header[Exploration_values$columns]))
            # regions=c('Whole','FWR1','CDR1','FWR2','CDR2','FWR3','CDR3','FWR4')
            region <- "FWR4"
            sunburst_df <- show_selected_features(tmp_sunburst, region = region)
            # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel=
            # 'Variables analyzed per each region in the exploration set')
            # )
            if (dim(sunburst_df)[1] ==
                0) {
              HTML(
                paste(
                  c("There is no", region, "sequence variables analyzed in the exploration set"),
                  collapse = " "
                )
              )
            } else {
              sunburstR::sund2b(
                sunburst_df, showLabels = F, rootLabel = paste(
                  c(region, "sequence variables analyzed in the exploration set"),
                  collapse = " "
                )
              )
            }
          }
        }
      )




      output$used_variables <- renderUI(
        {

          sel=data.frame(
            Region_of_the_variable=as.character(),
            Variable=as.character()
          )
          exp=data.frame(
            Region_of_the_variable=as.character(),
            Variable=as.character()
          )
          if (!is.null(Selection_values$rows) &&
              length(Selection_values$rows) >
              0 && !is.null(Selection_values$columns) &&
              length(Selection_values$columns) >
              0) {


            sel=data.frame(
              Region_of_the_variable=strsplit(
                Big_mem_values$Header[Selection_values$columns],
                split = "_"
              )[[1]][2],
              Variable=Big_mem_values$Header[Selection_values$columns]
            )

          }

          if (!is.null(Exploration_values$rows) &&
              length(Exploration_values$rows) >
              0 && !is.null(Exploration_values$columns) &&
              length(Exploration_values$columns) >
              0) {

            exp=data.frame(
              Region_of_the_variable=strsplit(
                Big_mem_values$Header[Exploration_values$columns],
                split = "_"
              )[[1]][2],
              Variable=Big_mem_values$Header[Exploration_values$columns]
            )

          }


          if (nrow(sel) > 0) sel$present_in_sel <- "Yes" else sel <- NULL
          if (nrow(exp) > 0) exp$present_in_exp <- "Yes" else exp <- NULL

          if (!is.null(sel) && !is.null(exp)) {
            data <- merge(sel, exp,
                          by = c("Region_of_the_variable", "Variable"),
                          all = T)
          } else if (!is.null(sel)) {
            data <- sel
            data$present_in_exp <- "No"
          } else if (!is.null(exp)) {
            data <- exp
            data$present_in_sel <- "No"
          } else {
            data <- data.frame(
              Region_of_the_variable = character(),
              Variable = character(),
              present_in_sel = logical(),
              present_in_exp = logical(),
              stringsAsFactors = F
            )
          }

          if (nrow(data) > 0) {
            data$present_in_sel[is.na(data$present_in_sel)] <- "No"
            data$present_in_exp[is.na(data$present_in_exp)] <- "No"
          }

          if(nrow(data)>0) {
            tagList(
              downloadButton("download_vardata_table",
                             "Download variable data for both sets",
                             onclick = "Reactable.downloadDataCSV('vardata_table', 'vardata_table.csv'); return false;"
              ),
              reactable(
                data, groupBy = c("Region_of_the_variable"),
                bordered = TRUE, elementId = "vardata_table", searchable = TRUE
              )
            )
          } else {
            HTML("There is an issue with at least one set.")
          }
        }
      )
      output$correlated_variables <- renderUI(
        {
          if (!is.null(Selection_values$rows) &&
              length(Selection_values$rows) >
              0 && !is.null(Selection_values$columns) &&
              length(Selection_values$columns) >
              0) {

            correlated_BGM <- bigstatsr::big_cor(
              Big_mem_values$Big_DF, ind.col = Selection_values$columns, ind.row = Selection_values$rows
            )


            index_correlated <- c()
            pairs_correlated <- list()
            react_table <- data.frame(
              Region_of_the_variable = as.character(), Variable = as.character(),
              Region_of_the_correlated_variable = as.character(), Correlated_variable = as.character(),
              Pearson_Correlation = as.numeric(), Raw_P_value = as.numeric()
            )
            for (column in c(1:ncol(correlated_BGM))) {
              if (any(
                correlated_BGM[c(1:ncol(correlated_BGM))[which(
                  c(1:ncol(correlated_BGM)) !=
                  column
                )],
                column] >= 0.8
              )) {

                tmp_real_corr <- c()
                tmp_pval_corr <- c()
                for (column_corr in c(1:ncol(correlated_BGM))[which(
                  correlated_BGM[c(1:ncol(correlated_BGM)),
                                 column] >= 0.8
                )]) {


                  if (column != column_corr) {
                    non_zero_rows <- unique(
                      c(
                        which(
                          Big_mem_values$Big_DF[Selection_values$rows, Selection_values$columns[column]] !=
                            0
                        ),
                        which(
                          Big_mem_values$Big_DF[Selection_values$rows, Selection_values$columns[column_corr]] !=
                            0
                        )
                      )
                    )


                    tmp_corr <- stats::cor(
                      Big_mem_values$Big_DF[Selection_values$rows[non_zero_rows],
                                            Selection_values$columns[column]], Big_mem_values$Big_DF[Selection_values$rows[non_zero_rows],
                                                                                                     Selection_values$columns[column_corr]]
                    )


                    if (!is.na(tmp_corr) &&
                        tmp_corr >= 0.8) {

                      if (gsub("_count", "", Big_mem_values$Header[Selection_values$columns[column]]) !=
                          gsub("_norm", "", Big_mem_values$Header[Selection_values$columns[column_corr]])) {
                        if (gsub("_norm", "", Big_mem_values$Header[Selection_values$columns[column]]) !=
                            gsub("_count", "", Big_mem_values$Header[Selection_values$columns[column_corr]])) {
                          pval_corr <- stats::cor.test(
                            Big_mem_values$Big_DF[Selection_values$rows, Selection_values$columns[column]],
                            Big_mem_values$Big_DF[Selection_values$rows, Selection_values$columns[column_corr]]
                          )$p.value

                          if (pval_corr <= 0.05) {
                            tmp_real_corr <- c(tmp_real_corr, column_corr)
                            tmp_pval_corr <- c(tmp_pval_corr, pval_corr)
                            react_table <- rbind(
                              react_table, data.frame(
                                Region_of_the_variable = strsplit(
                                  Big_mem_values$Header[Selection_values$columns[column]],
                                  split = "_"
                                )[[1]][2],
                                Variable = Big_mem_values$Header[Selection_values$columns[column]],
                                Region_of_the_correlated_variable = strsplit(
                                  Big_mem_values$Header[Selection_values$columns[column_corr]],
                                  split = "_"
                                )[[1]][2],
                                Correlated_variable = Big_mem_values$Header[Selection_values$columns[column_corr]],
                                Pearson_Correlation = correlated_BGM[column_corr,
                                                                     column], Raw_P_value = pval_corr
                              )
                            )

                          }

                        }
                      }

                    }

                  }

                }

                if (length(tmp_real_corr) >
                    0) {
                  pairs_correlated[[length(pairs_correlated) +
                                      1]] <- tmp_real_corr
                  names(pairs_correlated) <- c(
                    names(pairs_correlated)[1:(length(pairs_correlated) -
                                                 1)], column
                  )
                  index_correlated <- unique(c(index_correlated, column))

                }

              }
            }


            if (length(pairs_correlated) ==
                0) {
              HTML(
                "There is no correlation equal or greater than 0.8 between the sequence variables analyzed in the selection set."
              )
            } else {
              tagList(
                downloadButton("download_react_table",
                               "Download Data",
                               onclick = "Reactable.downloadDataCSV('Correlation-table', 'selection_set_data_correlated.csv'); return false;"
                ),
                reactable(
                  react_table, groupBy = c("Region_of_the_variable", "Variable"),
                  columns = list(Region_of_the_correlated_variable = colDef(aggregate = "unique")),
                  bordered = TRUE, elementId = "Correlation-table", searchable = TRUE,
                )
              )

            }
          }
        }
      )
      output$PCA_loadings <- renderUI(
        {
          sel_print_table=F
          if (!is.null(Selection_values$rows) &&
              length(Selection_values$rows) >
              0 && !is.null(Selection_values$columns) &&
              length(Selection_values$columns) >
              0) {

            sel_react_table=data.frame(
              Set=rep("Selection", length(Selection_values$columns)),
              Region_of_the_variable=strsplit(
                Big_mem_values$Header[Selection_values$columns],
                split = "_"
              )[[1]][2],
              Variable=Big_mem_values$Header[Selection_values$columns],
              PCA1_perc_load=round(Selection_values$Variance[,1], 5),
              PCA2_perc_load=round(Selection_values$Variance[,2], 5),
              PCA3_perc_load=round(Selection_values$Variance[,3], 5),
              PCA4_perc_load=round(Selection_values$Variance[,4], 5),
              PCA5_perc_load=round(Selection_values$Variance[,5], 5)
            )

            sel_print_table=T
          }


          exp_print_table=T
          if (!is.null(Exploration_values$rows) &&
              length(Exploration_values$rows) >
              0 && !is.null(Exploration_values$columns) &&
              length(Exploration_values$columns) >
              0) {



            exp_react_table=data.frame(
              Set=rep("Exploration", length(Exploration_values$columns)),
              Region_of_the_variable=strsplit(
                Big_mem_values$Header[Exploration_values$columns],
                split = "_"
              )[[1]][2],
              Variable=Big_mem_values$Header[Exploration_values$columns],
              PCA1_perc_load=round(Exploration_values$Variance[,1], 5),
              PCA2_perc_load=round(Exploration_values$Variance[,2], 5),
              PCA3_perc_load=round(Exploration_values$Variance[,3], 5),
              PCA4_perc_load=round(Exploration_values$Variance[,4], 5),
              PCA5_perc_load=round(Exploration_values$Variance[,5], 5)
            )


            exp_print_table=T
          }

          if(sel_print_table && exp_print_table) {
            tagList(
              downloadButton("download_PCA_Sel_table",
                             "Download PCA loadings for the Selection Set",
                             onclick = "Reactable.downloadDataCSV('Sel_PCA_table', 'selection_set_PCA_loadings.csv'); return false;"
              ),
              reactable(
                sel_react_table, groupBy = c("Region_of_the_variable"),
                bordered = TRUE, elementId = "Sel_PCA_table", searchable = TRUE
              ),
              HTML("<br>"),
              downloadButton("download_PCA_Exp_table",
                             "Download PCA loadings for the Exploration Set",
                             onclick = "Reactable.downloadDataCSV('Exp_PCA_table', 'exploration_set_PCA_loadings.csv'); return false;"
              ),
              reactable(
                exp_react_table, groupBy = c("Region_of_the_variable"),
                bordered = TRUE, elementId = "Exp_PCA_table", searchable = TRUE
              )
            )
          } else {
            HTML("There are no PCAs currently being plotted")
          }
        }
      )

      output$about_samples <- renderUI(
        {

          if (!is.null(Selection_values$rows) &&
              length(Selection_values$rows) >
              0 && !is.null(Selection_values$columns) &&
              length(Selection_values$columns) >
              0) {



            sel_table <- as.data.frame(table(Big_mem_values$Short_DF$Patient_Sample[Selection_values$rows],
                                             Big_mem_values$Short_DF$V_and_J[Selection_values$rows]))
            names(sel_table) <- c("Group", "V_J_genes", "Number_of_sequences")

            total_group <- stats::ave(sel_table$Number_of_sequences, sel_table$Group, FUN = sum)

            sel_table$Percentages_of_total_sequences_in_group <- 100 * sel_table$Number_of_sequences / total_group

          }

          if (!is.null(Exploration_values$rows) &&
              length(Exploration_values$rows) >
              0 && !is.null(Exploration_values$columns) &&
              length(Exploration_values$columns) >
              0) {


            exp_table <- as.data.frame(table(Big_mem_values$Short_DF$Patient_Sample[Exploration_values$rows],
                                             Big_mem_values$Short_DF$V_and_J[Exploration_values$rows]))
            names(exp_table) <- c("Group", "V_J_genes", "Number_of_sequences")

            total_group <- stats::ave(exp_table$Number_of_sequences, exp_table$Group, FUN = sum)

            exp_table$Percentages_of_total_sequences_in_group <- 100 * exp_table$Number_of_sequences / total_group

          }

          if(nrow(sel_table)>0 && nrow(exp_table)>0) {
            tagList(
              downloadButton("download_Sel_gene_table",
                             "Download gene usage for the Selection Set",
                             onclick = "Reactable.downloadDataCSV('Sel_gene_table', 'selection_set_PCA_loadings.csv'); return false;"
              ),
              reactable(
                sel_table, groupBy = c("Group"),
                bordered = TRUE, elementId = "Sel_gene_table", searchable = TRUE
              ),
              HTML("<br>"),
              downloadButton("download_Exp_gene_table",
                             "Download gene usage for the Exploration Set",
                             onclick = "Reactable.downloadDataCSV('Exp_gene_table', 'exploration_set_PCA_loadings.csv'); return false;"
              ),
              reactable(
                exp_table, groupBy = c("Group"),
                bordered = TRUE, elementId = "Exp_gene_table", searchable = TRUE
              )
            )
          } else {
            HTML("There is an issue with one of the sets")
          }
        }
      )
      output$numerical_covariables_plot <- renderUI(
        {
          if (!is.null(Selection_values$rows) &&
              length(Selection_values$rows) >
              0 && !is.null(Selection_values$columns) &&
              length(Selection_values$columns) >
              0 && input$numerical_covariable_file != "") {

            plot_heatmap <- F

            if ("Sequence_ID" %in% colnames(input$numerical_covariable_file)) {


              plot_heatmap <- T
            } else if ("Cell_ID" %in% colnames(input$numerical_covariable_file)) {




              plot_heatmap <- T
            }


            if (plot_heatmap) {
              HTML("")

            }


          }
        }
      )

      # 4.Feature exploration ####

      observeEvent(
        input$use_what, {
          if ("Reconstructed germline" %!in% input$use_what) {
            shinyjs::show("show_reconstructed", anim = TRUE)
          } else {
            shinyjs::hide("show_reconstructed", anim = TRUE)
          }
        }
      )
      o_columns <- metaObserve2({
        shiny::req(Selection_values$columns)
        isolate(metaExpr({
          tmp_choices <- Big_mem_values$Header[Selection_values$columns]
          # updateSelectInput(session, inputId = 'plot_feature',
          # choices=tmp_choices, selected=tmp_choices[1])
          updateSelectizeInput(
            session, inputId = "plot_feature", choices = tmp_choices, selected = tmp_choices[1],
            server = T, options = list(maxOptions = length(tmp_choices))
          )
        }))
      })



      output$Group_selection_for_feature <- renderUI(
        {
          if (is.null(Selection_values$rows) &&
              length(Selection_values$rows) >
              0 && !is.null(Selection_values$columns) &&
              length(Selection_values$columns) >
              0) {
            selectInput(
              "plot_color_feature", "Color data points by", choices = c("Sample", "Patient", "Group", "Subgroup", "Chain"),
              selected = "Group"
            )
          } else {
            sample_info <- sample_info_react$table
            choices_color <- c(colnames(sample_info))
            names(choices_color) <- c(colnames(sample_info))
            if ("Chain" %in% colnames(Big_mem_values$Short_DF) &&
                "Chain" %!in% choices_color) {
              tmp_names_color <- names(choices_color)
              choices_color <- c(choices_color, "Chain")
              names(choices_color) <- c(tmp_names_color, "Ig Chain")
            }
            tmp_names_choices_color <- names(choices_color)
            choices_color <- c(choices_color, "Best_V", "Best_J", "V_and_J")
            names(choices_color) <- c(
              tmp_names_choices_color, "Main V gene isoform", "Main J gene isoform",
              "Main V and J gene isoforms"
            )
            choices_color <- choices_color[which(choices_color != "Additional_info")]
            choices_color <- choices_color[which(choices_color != "Folder")]
            choices_color <- choices_color[which(choices_color != "Filename")]
            # choices_color=choices_color[which(choices_color !=
            # 'V_and_J')]



            if ("D_region" %in% colnames(Big_mem_values$Short_DF)) {
              tmp_names_color <- names(choices_color)
              choices_color <- c(choices_color, "Best_D")
              names(choices_color) <- c(tmp_names_color, "Main D gene isoform")
            }
            if ("V_and_D_and_J" %in% colnames(Big_mem_values$Short_DF)) {
              tmp_names_color <- names(choices_color)
              choices_color <- c(choices_color, "V_and_D_and_J")
              names(choices_color) <- c(tmp_names_color, "Main V, D and J gene isoforms")
            }
            if ("C_region" %in% colnames(Big_mem_values$Short_DF)) {
              tmp_names_color <- names(choices_color)
              choices_color <- c(choices_color, "C_region")
              names(choices_color) <- c(tmp_names_color, "Constant region")
            }
            if ("Clone_ID" %in% colnames(Big_mem_values$Short_DF)) {
              tmp_names_color <- names(choices_color)
              choices_color <- c(choices_color, "Clone_ID")
              names(choices_color) <- c(tmp_names_color, "Clone ID")
            }

            if ("Dominance" %in% colnames(Big_mem_values$Short_DF)) {
              tmp_names_color <- names(choices_color)
              choices_color <- c(choices_color, "Dominance")
              names(choices_color) <- c(tmp_names_color, "Dominance")
            }


            selectInput("plot_color_feature", "Color data points by", choices = choices_color, selected = "Group")
          }

        }
      )

      output$Violin_feature_plot <- metaRender(renderPlotly,
                                               {
                                                 if (!is.null(Selection_values$rows) &&
                                                     length(Selection_values$rows) >
                                                     0 && !is.null(Selection_values$columns) &&
                                                     length(Selection_values$columns) >
                                                     0) {
                                                   fig_violin <- AbSolution::draw_feature_violinplot(
                                                     values = Big_mem_values$Big_DF[, which(Big_mem_values$Header == input$plot_feature)],
                                                     name_values = input$plot_feature, sequence_info_df = Big_mem_values$Short_DF,
                                                     group_info = input$groups_selected, additional_group_info = input$plot_color_feature,
                                                     show_reconstructed = input$show_reconstructed,
                                                     compare_opposites=input$compare_opposites, selected_rows = Selection_values$rows,
                                                     selected_subclones = if(length(ID_selected_values$subclones$key)>0){sapply(strsplit(ID_selected_values$subclones$key, split="_&_"), `[`, 1)}else{NULL} ,
                                                     selected_clones = ID_selected_values$clones_subclones_id,
                                                     hide_dots = input$hide_points, seed=input$seed,
                                                     really_hide_dots=input$really_hide_points,
                                                     width = input$pixels_width, height = input$pixels_height,
                                                     img_type= input$img_type, scale=input$labelscale,
                                                     primary_color=input$primary_color
                                                   )
                                                   ID_selected_values$subclones <- ..(event_data("plotly_selected"))

                                                   fig_violin
                                                 }

                                               }
      )





      # 5.Clonal exploration ####


      # in server.R create reactiveVal
      clonal_group_current_selection <- reactiveVal(NULL)

      # now store your current selection in the reactive value
      o_clonal_group <- metaObserve2({
        shiny::req(input$clonal_group)
        isolate(metaExpr({
          clonal_group_current_selection(input$clonal_group)
        }))
      })


      output$clonal_group_output <- renderUI(
        {
          choices <- c("")
          tmp_selection <- ""
          if (!is.null(Exploration_values$rows) &&
              length(Exploration_values$rows) >
              0 && !is.null(Exploration_values$columns) &&
              length(Exploration_values$columns) >
              0) {
            choices <- c(
              colnames(Big_mem_values$Short_DF)[grepl("Clone_", colnames(Big_mem_values$Short_DF))]
            )
            if (!is.null(clonal_group_current_selection()) &&
                clonal_group_current_selection()[1] != "") {
              tmp_selection <- clonal_group_current_selection()
            } else {
              tmp_selection <- choices[1]
            }
          }

          selectInput(
            "clonal_group", "Use this clonal definition:", choices = choices,
            selected = tmp_selection
          )
        }
      )

      new_dominance_calculation <- metaReactive(
        {
          list(input$clonal_group, input$dominance_threshold)
        }, varname = "new_dominance_calculation"
      )


      o_new_dominance_calculation <- metaObserve2({
        shiny::req(new_dominance_calculation())
        isolate(metaExpr({
          if (!is.null(Selection_values$rows) &&
              length(Selection_values$rows) >
              0 && !is.null(Selection_values$columns) &&
              length(Selection_values$columns) >
              0) {
            if (input$clonal_group != "") {
              Big_mem_values$Short_DF[, "Dominance"] <- ""

              for (pat_sample in unique(Big_mem_values$Short_DF$Patient_Sample)) {
                index_tmp <- which(Big_mem_values$Short_DF[, "Patient_Sample"] == pat_sample)
                dom_table <- table(Big_mem_values$Short_DF[[input$clonal_group]][index_tmp])/length(Big_mem_values$Short_DF[[input$clonal_group]][index_tmp])
                dom_table <- sapply(
                  dom_table, function(z) if (z >=
                                             input$dominance_threshold/100) {
                    "Dominant"
                  } else {
                    "Non-dominant"
                  }
                )

                Big_mem_values$Short_DF[index_tmp, "Dominance"] <- dom_table[match(Big_mem_values$Short_DF[[input$clonal_group]][index_tmp], names(dom_table))]
              }

            }
          }
        }))
      })


      o_calculate_new_clone <- metaObserve2({
        shiny::req(input$calculate_new_clone)
        isolate(metaExpr({
          if (!is.null(Selection_values$rows) &&
              length(Selection_values$rows) >
              0 && !is.null(Selection_values$columns) &&
              length(Selection_values$columns) >
              0) {
            if (!(paste(
              "Clone", input$clonal_level_group, input$clonal_region_group, input$new_clonal_group,  input$identity_clonal_group, if ( input$calculate_shared_clones) {
                "shared"
              }else {"non-shared"}, "simil", sep = "_"
            ) %in%
            colnames(Big_mem_values$Short_DF))) {

              Big_mem_values$Short_DF <- AbSolution::calculate_clone(
                seq_df = Big_mem_values$Short_DF, clonotype = input$new_clonal_group,
                AA_or_NT = input$clonal_level_group, region = input$clonal_region_group,
                percentage = input$identity_clonal_group, calculate_shared_clones = input$calculate_shared_clones
              )

            }
          }
        }))
      })

      output$Violin_plot <- metaRender(renderPlotly,
                                       {
                                         if (!is.null(Selection_values$rows) &&
                                             length(Selection_values$rows) >
                                             0 && !is.null(Selection_values$columns) &&
                                             length(Selection_values$columns) >
                                             0 && !is.null(input$clonal_group) &&
                                             input$clonal_group != "") {


                                           # TOFIX ####
                                           fig_violin <-  AbSolution::draw_violinplots(
                                             seq_df = Big_mem_values$Short_DF, group = "Patient_Sample", selected_rows = Selection_values$rows,
                                             clonotype = input$clonal_group, AA_or_NT = input$clonal_level_group,
                                             region = input$clonal_region_group, percentage = input$identity_clonal_group,
                                             freq_filter = input$filter_clonal_group, Selected_clones = ID_selected_values$clones,
                                             dominance_threshold = input$dominance_threshold, seed=input$seed,
                                             really_hide_dots=input$really_hide_points_clones,
                                             width = input$pixels_width, height = input$pixels_height,
                                             img_type= input$img_type, scale=input$labelscale
                                           )
                                           ID_selected_values$clones <- ..(event_data("plotly_selected"))
                                           fig_violin
                                         }

                                       }
      )

      output$upset_plot <- metaRender(upsetjs::renderUpsetjs,
                                      {
                                        if (!is.null(Selection_values$rows) &&
                                            length(Selection_values$rows) >
                                            0 && !is.null(Selection_values$columns) &&
                                            length(Selection_values$columns) >
                                            0 && !is.null(input$clonal_group) &&
                                            input$clonal_group != "") {

                                          fig_upset <- AbSolution::draw_upsetplot(
                                            seq_df = Big_mem_values$Short_DF, group = "Patient_Sample", selected_rows = Selection_values$rows,
                                            clonotype = input$clonal_group, AA_or_NT = input$clonal_level_group,
                                            region = input$clonal_region_group, percentage = input$identity_clonal_group,
                                            freq_filter = input$filter_clonal_group, Selected_clones = ID_selected_values$clones
                                          )
                                          fig_upset |>
                                            upsetjs::interactiveChart("click", events_nonce = TRUE) |>
                                            upsetjs::generateDistinctIntersections()
                                        }

                                      }
      )


      o_intersection_samples <- metaObserve2({
        shiny::req(ID_selected_values$intersection_samples)
        isolate(metaExpr({
          upsetjs::upsetjsProxy("upset_plot", session) |>
            upsetjs::setSelection(ID_selected_values$intersection_samples)
        }))
      })


      o_upset_plot_click <- metaObserve2({
        shiny::req(input$upset_plot_click)
        isolate(metaExpr({
          if (isTRUE(input$upset_plot_click[["isSelected"]])) {
            ID_selected_values$intersection_samples <- ""
          } else {
            ID_selected_values$intersection_samples <- input$upset_plot_click

          }
        }))
      })

      output$coloring_scatterplot_clones <- renderUI(
        {
          choices <- c("Clonal sharing", "VJ usage", Big_mem_values$categorical_newfields)
          tmp_selection <- choices[1]

          selectInput(
            "color_scatterplot_clones", "Color:", choices = choices,
            selected = tmp_selection
          )
        }
      )

      output$Comparison_plot <- metaRender(renderPlotly,
                                           {
                                             if (!is.null(ID_selected_values$intersection_samples$setNames) &&
                                                 !is.null(ID_selected_values$intersection_samples$elems) &&
                                                 !is.null(Selection_values$rows) &&
                                                 length(Selection_values$rows) >
                                                 0 && !is.null(Selection_values$columns) &&
                                                 length(Selection_values$columns) >
                                                 0 && !is.null(input$clonal_group) &&
                                                 input$clonal_group != "") {

                                               if (length(unlist(ID_selected_values$intersection_samples$setNames)) >
                                                   1) {
                                                 fig_sharedclones <- AbSolution::draw_sharedclonesplot(
                                                   seq_df = Big_mem_values$Short_DF, sets = ID_selected_values$intersection_samples$setNames,
                                                   group = "Patient_Sample", selected_rows = Selection_values$rows,
                                                   clonotype = input$clonal_group, AA_or_NT = input$clonal_level_group,
                                                   region = input$clonal_region_group, percentage = input$identity_clonal_group,
                                                   freq_filter = input$filter_clonal_group, Selected_clones = ID_selected_values$clones,
                                                   dominance_threshold = input$dominance_threshold,
                                                   color_by=input$color_scatterplot_clones, colorblind = input$colorblind_mode,
                                                   Big_mem_color_values=Big_mem_color_values, seed=input$seed,
                                                   width = input$pixels_width, height = input$pixels_height,
                                                   img_type= input$img_type, scale=input$labelscale
                                                 )
                                                 ID_selected_values$clones <- event_data("plotly_selected")
                                                 fig_sharedclones |>
                                                   plotly::toWebGL()
                                               }

                                             }

                                           }
      )

      clonal_merged_data <- eventReactive(
        ID_selected_values$clones, {

          dt_sel <- NULL
          # print("ID_selected_values$clones")
          # print(ID_selected_values$clones)
          # print("-------------")
          # print(table(Selection_values$Scores$Selected))
          # print("-------------")
          if (!is.null(ID_selected_values$clones) && length(ID_selected_values$clones) >0) {
            keep_cols <- colnames(Big_mem_values$Short_DF)
            dt_sel <- Big_mem_values$Short_DF[Selection_values$rows[which(Selection_values$Scores$Selected == "Clones")], keep_cols, with = FALSE ]
          }
          # print(dt_sel)
          # print("-------------")

          if (is.null(dt_sel) || nrow(dt_sel) == 0 ) {
            dt_sel <- data.frame(Message= c("No user-selected clones"))
          } else{
            clone_col <- input$clonal_group
            seq_col   <- "ID"

            dt_sel <- dt_sel[
              , .(Sequences = paste(unique(get(seq_col)), collapse = ";")),
              by = .(Clone = get(clone_col))
            ][order(Clone)]

          }
          dt_sel
        }
      )  #eventReactive

      output$Clonal_table <- metaRender(DT::renderDataTable,
                                        {
                                          rendered_table <- clonal_merged_data()

                                          DT::datatable(
                                            rendered_table, caption = "User-selected clones",
                                            extensions = "Buttons",
                                            options = list(scrollX = TRUE,
                                                           dom = "Bfrtip",
                                                           buttons = c("copy", "csv")
                                            )
                                          )
                                        }, options = list(
                                          dom = "Bfrtip", list(
                                            list(
                                              extend = "csv",
                                              text = "csv",
                                              exportOptions = list(modifier = list(page = "all"))
                                            ),
                                            "copy"
                                          )
                                        ), server=F
      )

      # Help ####
      shinyDirChoose(input, "sFSS_folder", roots = volumes(), session = session)
      observeEvent(input$downloadsFSS, {

        sFFS_name=gsub(" ","",gsub("/","_",input$sFSS_name))
        sFFS_name=if(!is.null(sFFS_name) && sFFS_name !="") {sFFS_name} else {"sFFS"}
        sFFS_name=paste0(sFFS_name, ".zip")

        if (all(
          unlist(input$sFSS_folder) ==
          0
        )) {

        } else {
          sFFS_name=file.path(shinyFiles::parseDirPath(volumes(), input$sFSS_folder), sFFS_name)
          tryCatch({
            progress <- Progress$new(session, min = 0, max = 1)
            on.exit(progress$close())

            progress$set(message = "Downloading ENCORE project structure.", value = 0.1)
            download.file(url = "https://github.com/EDS-Bioinformatics-Laboratory/ENCORE/archive/refs/heads/main.zip",
                          destfile = sFFS_name,
                          quiet = !verbose)
            progress$set(message = "Unzipping ENCORE project structure.", value = 0.7)
            utils::unzip(zipfile = sFFS_name, exdir =shinyFiles::parseDirPath(volumes(), input$sFSS_folder))
            progress$set(message = "Cleaning...", value = 0.9)
            file.remove(sFFS_name)
            Sys.sleep(1)
            progress$set(message = "Done!", value = 1)
          }, error = function(e) {
            shiny::showNotification("File couldn't be downloaded. Please check your internet connection.", type = "error")
          })


        }


      })

      # getData <- metaReactive(
      #   {
      #   d <- switch(
      #     input$data,
      #     fly = fly,
      #     happy = happy,
      #     custom = if (length(input$data_file)) read.csv(input$data_file$datapath)
      #   )
      #   if (is.null(d)) return(NULL)
      #   d <- mutate_all(d, as.character)
      #   if (isTRUE(input$na.rm)) na.omit(d) else d
      # }, varname = "getData")


      output$export_cond <- renderUI(
        {
          wait_for_export()
        }
      )

      wait_for_export <- metaReactive(
        {
          if (!is.null(input$clonal_group) && input$clonal_group != "") {
            ec <- newExpansionContext()


            setupRmarkdown2 <- expandChain(
              parameters_only_for_shinymeta_info(),
              rlang::parse_expr("author=parameters_only_for_shinymeta_info$author"),
              rlang::parse_expr("usermail=parameters_only_for_shinymeta_info$usermail"),
              rlang::parse_expr("user_0_comments=parameters_only_for_shinymeta_info$user_0_comments"),
              rlang::parse_expr("user_1_comments=parameters_only_for_shinymeta_info$user_1_comments"),
              rlang::parse_expr("user_2_comments=parameters_only_for_shinymeta_info$user_2_comments"),
              rlang::parse_expr("user_3_comments=parameters_only_for_shinymeta_info$user_3_comments"),
              rlang::parse_expr("user_4_comments=parameters_only_for_shinymeta_info$user_4_comments"),
              rlang::parse_expr("user_5_comments=parameters_only_for_shinymeta_info$user_5_comments"),

              # rlang::parse_expr("author=2"),
              # rlang::parse_expr("usermail=2"),
              # rlang::parse_expr("user_0_comments=0"),
              # rlang::parse_expr("user_1_comments=1"),
              # rlang::parse_expr("user_2_comments=2"),
              # rlang::parse_expr("user_3_comments=3"),
              # rlang::parse_expr("user_4_comments=4"),
              # rlang::parse_expr("user_5_comments=5"),
              .expansionContext = ec
            )

            loading <- expandChain(
              quote({
                library(AbSolution)
                # library(shiny)
                # library(dplyr)
                # library(plotly)
              }),
              .expansionContext = ec
            )

            parameter_setup <- expandChain(
              output$raw_sample_file_out(),
              parameters_only_for_shinymeta(),
              rlang::parse_expr("input=parameters_only_for_shinymeta$input"),
              rlang::parse_expr("folder_values=parameters_only_for_shinymeta$folder_values"),
              rlang::parse_expr("Exploration_values=parameters_only_for_shinymeta$Exploration_values"),
              rlang::parse_expr("Selection_values=parameters_only_for_shinymeta$Selection_values"),
              rlang::parse_expr("Big_mem_values=parameters_only_for_shinymeta$Big_mem_values"),
              rlang::parse_expr("Big_mem_color_values=parameters_only_for_shinymeta$Big_mem_color_values"),
              rlang::parse_expr("ID_selected_values=parameters_only_for_shinymeta$ID_selected_values"),
              rlang::parse_expr("sample_info_react <- parameters_only_for_shinymeta$sample_info_react"),
              rlang::parse_expr("Features_values <- parameters_only_for_shinymeta$Features_values"),
              .expansionContext = ec
            )

            parse_input <- expandChain(
              # o_rws(),
              o_PA(),
              .expansionContext = ec
            )

            feature_calculation <- expandChain(
              o_Featdet(),
              .expansionContext = ec
            )

            file_load <- expandChain(
              o_MtAr(),
              o_calcColor(),
              .expansionContext = ec
            )

            clone <- expandChain(
              o_calculate_new_clone(),
              output$Violin_plot(),
              output$upset_plot(),
              output$Comparison_plot(),
              .expansionContext = ec
            )

            PCA_ex <- expandChain(
              o_tLEx(),
              o_tLS(),
              o_tLCSel(),
              output$Selection_plot(),
              .expansionContext = ec
            )

            feature <- expandChain(
              output$Violin_feature_plot(),
              .expansionContext = ec
            )

            shiny::req(length(setupRmarkdown2)>0 )

            shiny::req(length(loading)>0 )

            shiny::req(length(parameter_setup)>0 )

            shiny::req(length(parse_input)>0 )

            shiny::req(length(feature_calculation)>0 )

            shiny::req(length(file_load)>0 )

            shiny::req(length(clone)>0 )

            shiny::req(length(PCA_ex)>0 )

            shiny::req(length(feature)>0 )

            HELP_values$stage <- 4
            downloadButton("download_script", "Export report in ENCORE format")
          } else {

            actionButton("download_fake", "Export not yet possible", class = "btn-warning",
                         style="color: #fff; background-color: #e63946; border-color: #c1121f")

          }
        }, varname = "wait_for_export"
      )

      output$download_script <- downloadHandler(
        paste0(input$analysis_name, ".zip"),
        content = function(out) {
          # saveRDS(getData(), "data.rds")
          # on.exit(unlink("data.rds"), add = TRUE)

          ec <- newExpansionContext()
          # ec$substituteMetaReactive(getData, function() {
          #   metaExpr({readRDS("data.rds")})
          # })


          setupRmarkdown2 <- expandChain(
            parameters_only_for_shinymeta_info(),
            rlang::parse_expr("author=parameters_only_for_shinymeta_info$author"),
            rlang::parse_expr("usermail=parameters_only_for_shinymeta_info$usermail"),
            rlang::parse_expr("user_0_comments=parameters_only_for_shinymeta_info$user_0_comments"),
            rlang::parse_expr("user_1_comments=parameters_only_for_shinymeta_info$user_1_comments"),
            rlang::parse_expr("user_2_comments=parameters_only_for_shinymeta_info$user_2_comments"),
            rlang::parse_expr("user_3_comments=parameters_only_for_shinymeta_info$user_3_comments"),
            rlang::parse_expr("user_4_comments=parameters_only_for_shinymeta_info$user_4_comments"),
            rlang::parse_expr("user_5_comments=parameters_only_for_shinymeta_info$user_5_comments"),

            .expansionContext = ec
          )

          loading <- expandChain(
            quote({
              library(AbSolution)
              # library(shiny)
              # library(dplyr)
              # library(plotly)
            }),
            .expansionContext = ec
          )

          parameter_setup <- expandChain(
            output$raw_sample_file_out(),
            parameters_only_for_shinymeta(),
            rlang::parse_expr("input=parameters_only_for_shinymeta$input"),
            rlang::parse_expr("folder_values=parameters_only_for_shinymeta$folder_values"),
            rlang::parse_expr("Exploration_values=parameters_only_for_shinymeta$Exploration_values"),
            rlang::parse_expr("Selection_values=parameters_only_for_shinymeta$Selection_values"),
            rlang::parse_expr("Big_mem_values=parameters_only_for_shinymeta$Big_mem_values"),
            rlang::parse_expr("Big_mem_color_values=parameters_only_for_shinymeta$Big_mem_color_values"),
            rlang::parse_expr("ID_selected_values=parameters_only_for_shinymeta$ID_selected_values"),
            rlang::parse_expr("sample_info_react <- parameters_only_for_shinymeta$sample_info_react"),
            rlang::parse_expr("Features_values <- parameters_only_for_shinymeta$Features_values"),
            .expansionContext = ec
          )

          parse_input <- expandChain(
            # o_rws(),
            o_PA(),
            .expansionContext = ec
          )

          feature_calculation <- expandChain(
            o_Featdet(),
            .expansionContext = ec
          )

          file_load <- expandChain(
            o_MtAr(),
            o_calcColor(),
            .expansionContext = ec
          )

          clone <- expandChain(
            o_calculate_new_clone(),
            output$Violin_plot(),
            output$upset_plot(),
            output$Comparison_plot(),
            .expansionContext = ec
          )

          PCA_ex <- expandChain(
            o_tLEx(),
            o_tLS(),
            o_tLCSel(),
            output$Selection_plot(),
            .expansionContext = ec
          )

          feature <- expandChain(
            output$Violin_feature_plot(),
            .expansionContext = ec
          )


          buildRmdBundle_alt(
            "inst/app/www/Report_template.Rmd", out,
            vars = list(
              setupRmarkdown2=setupRmarkdown2,
              loading=loading,
              parse_input=parse_input,
              feature_calculation=feature_calculation,
              parameter_setup = parameter_setup,
              file_load=file_load,
              clone=clone,
              PCA_ex=PCA_ex,
              feature=feature
            ),
            golem_renv=input$include_docker_renv,
            render_args = list(output_format = "all",
                               output_dir = "Results",
                               intermediates_dir = "Notebook",
                               knit_root_dir = "Notebook"),
            include_files = if (input$include_data_report) {
              if(is.null(input$preinfolder_AIRR)) {
                stats::setNames(c(folder_values$Featured,
                                  file.path(
                                    shinyFiles::parseDirPath(volumes, input$base_folder),
                                    "Sample_summary.txt"
                                  )),
                                c("Data/Dataset/Processed/2.Feature_determination",
                                  "Data/Dataset/Meta/Sample_summary.txt"))
              } else {
                if (sample_info_react$test_status) {
                  tsv_files=file.path(shinyFiles::parseDirPath(volumes, input$base_folder), "Example10x_alakazam.tsv")
                } else {
                  tsv_files=file.path(shinyFiles::parseDirPath(volumes, input$preinfolder_AIRR),  paste(sample_table_react()[,1], ".tsv", sep=""))
                }
                stats::setNames(c(folder_values$Featured,
                                  as.vector(tsv_files), file.path(
                                    shinyFiles::parseDirPath(volumes, input$base_folder),
                                    "Sample_summary.txt"
                                  )),
                                c("Data/Dataset/Processed/2.Feature_determination",
                                  sapply(tsv_files, function(z)
                                    paste("Data/Dataset/Raw",basename(z),sep="/")),
                                  "Data/Dataset/Meta/Sample_summary.txt"))
              }
            } else {
              list()
            }
            # include_files = c("data.rds")
          )
        }
      )


      numbers <- metaReactive(
        {
          validate(
            need(
              is.numeric(input$seed),
              "Please input a number"
            )
          )
        }, varname = "numbers"
      )
      output$onlynumbers <- renderPrint(
        {
          numbers()
        }
      )

      numberswidth <- metaReactive(
        {
          validate(
            need(
              is.numeric(input$pixels_width),
              "Please input a number"
            )
          )
        }, varname = "numberswidth"
      )
      output$onlynumberswidth <- renderPrint(
        {
          numberswidth()
        }
      )

      numbersheight <- metaReactive(
        {
          validate(
            need(
              is.numeric(input$pixels_height),
              "Please input a number"
            )
          )
        }, varname = "numbersheight"
      )
      output$onlynumbersheight <- renderPrint(
        {
          numbersheight()
        }
      )

      numberslabelscale <- metaReactive(
        {
          validate(
            need(
              is.numeric(input$labelscale),
              "Please input a number"
            )
          )
        }, varname = "numberslabelscale"
      )
      output$onlynumberslabelscale <- renderPrint(
        {
          numberslabelscale()
        }
      )


      output$Conditional_Action_downloadsFSS <- renderUI(
        {
          if (!all(is.numeric(unlist(input$sFSS_folder)))) {
            actionButton("downloadsFSS", "Download sFSS template",
                         icon = icon("cloud-arrow-down", lib = "font-awesome"))
          }

        }
      )

      HELP_values <- reactiveValues(stage = 0)
      output$HELP_output <- renderUI(
        {
          bs4Dash::timelineBlock(
            width = 12,
            reversed = T,
            bs4Dash::timelineEnd(shiny::icon("hourglass-start", style = "color: rgb(61, 153, 112)"), color = "teal"),
            ## 0
            if (HELP_values$stage == 0) { bs4Dash::timelineLabel("You are here", color = "pink")},
            bs4Dash::timelineItem(
              elevation =  if (HELP_values$stage == 0) { 2} else {NULL},
              title = if (HELP_values$stage == 0) {HTML( "<b>Step 0 - Project information</b>")} else { "Step 0 - Project information"},
              icon =  if (HELP_values$stage == 0) { icon("location-dot")} else {icon("ellipsis")} ,
              color = if (HELP_values$stage == 0) { "pink"} else {"gray"},
              time =  if (HELP_values$stage == 0) { "Now"} else {"Before"},
              # footer = "Here is the footer",
              if (HELP_values$stage == 0) {                HTML(
                "In this tab we will define what type of project, samples and workspace will be used. You need to indicate:<br>
                      <ul>
                      <li>If your data is <b>TCR or BCR</b> repertoire-based.</li>
                      <li>Tell us about your data. <b>Include the required descriptive table</b>.</li>
                      <li>Specify <b>the folder where all the related-files will be created</b>.</li>
                      <li>You can also <b>try AbSolution with a <a href='https://alakazam.readthedocs.io/en/stable/topics/Example10x/'>test dataset</a> </b>.</li>
                     </ul>")
              } else {NULL}
            ),

            if (HELP_values$stage == 1) { bs4Dash::timelineLabel("You are here", color = "pink")},

            if(HELP_values$stage != 2.5) {
              ## 1

              bs4Dash::timelineItem(
                elevation =  if (HELP_values$stage == 1) { 2} else {NULL},
                title = if (HELP_values$stage == 1) {HTML( "<b>Step 1 - AIRR-Seq conversion</b>")} else { "Step 1 - AIRR-Seq conversion"},
                icon =  if (HELP_values$stage == 1) { icon("location-dot")} else {icon("ellipsis")} ,
                color = if (HELP_values$stage == 1) { "pink"} else {"gray"},
                time =  if (HELP_values$stage == 1) { "Now"} else if(HELP_values$stage < 1) {"Later"} else {"Before"},
                # footer = "Here is the footer",
                if (HELP_values$stage == 1) {                HTML(
                  "In this tab we will reconstruct the repertoire germlines based on your input. AbSolution works with the full Fab sequences or with sequences with partially-sequenced FWR1/4. You need to indicate:<br>
                   <ul>
                      <li>The <b>folder where your data is located</b>. All files must be in the same folder. Filenames should match the filenames you indicated in the previous step. If they are found in the folder, it will say 'Yes' in Found_in_folder, otherwise 'No' and you need to either a) change the folder, b) move the file to the folder, c) change the input table or d) change the filename.</li>
                      <li>Additional information to better <b>reconstruct the germlines</b>. Do you want to keep the CDR3 D-gene segment as is or use the D gene information to reconstruct it? Are the FWR1 or the FWR4 partially sequenced? Is there a pre or post-Fab sequence to be removed?</li>
                  </ul> <br>
                  Please confirm with the pre-visualization that the sequences are correctly parsed and reconstructed. A lowercase in the repertoire sequence means an insertion. A lowercase in the reconstructed germline means a deletion in the repertoire.")
                } else {NULL}
              )},
            if (HELP_values$stage == 2) { bs4Dash::timelineLabel("You are here", color = "pink")},
            if (HELP_values$stage == 2.5) { bs4Dash::timelineLabel("You are here", color = "pink")},
            if(HELP_values$stage != 2.5) {
              # 2
              bs4Dash::timelineItem(
                elevation =  if (HELP_values$stage == 2) { 2} else {NULL},
                title = if (HELP_values$stage == 2) {HTML( "<b>Step 2 - Sequence feature determination</b>")} else { "Step 2 - Sequence feature determination"},
                icon =  if (HELP_values$stage == 2) { icon("location-dot")} else {icon("ellipsis")} ,
                color = if (HELP_values$stage == 2) { "pink"} else {"gray"},
                time =  if (HELP_values$stage == 2) { "Now"} else if(HELP_values$stage < 2) {"Later"} else {"Before"},
                # footer = "Here is the footer",
                if (HELP_values$stage == 2) {                HTML(
                  "In this tab we will calculate the features (physicochemical, composition, etc.) at the NT and AA level for the whole sequence (repertoire and reconstructed germline) and its individual regions.")
                } else {NULL}
              )
            } else {
              ## 2.5
              bs4Dash::timelineItem(
                elevation =  if (HELP_values$stage == 2.5) { 2} else {NULL},
                title = if (HELP_values$stage == 2.5) {HTML( "<b>Steps 1&2. Select your FBM and associated files</b>")} else { "Steps 1&2. Select your FBM and associated files"},
                icon =  if (HELP_values$stage == 2.5) { icon("location-dot")} else {icon("ellipsis")} ,
                color = if (HELP_values$stage == 2.5) { "pink"} else {"gray"},
                time =  if (HELP_values$stage == 2.5) { "Now"} else {"Later"},
                # footer = "Here is the footer",
                if (HELP_values$stage == 2.5) {                HTML(
                  "In this tab you have to indicate the <b>folder where your pre-calculated FBM files are located</b>. All files must be in the same folder. Filenames should match the filenames you indicated in the previous step. If they are found in the folder, it will say 'Yes' in Found_in_folder, otherwise 'No' and you need to either a) change the folder, b) move the file to the folder, c) change the input table or d) change the filename. You can generate dummy datasets with the same number of mutations to use as a negative reference in the next steps.>
                  </ul> ")
                } else {NULL}
              )
            },

            # 3
            if (HELP_values$stage == 3) { bs4Dash::timelineLabel("You are here, the following steps are interconnected!", color = "pink")},
            bs4Dash::timelineItem(
              elevation =  if (HELP_values$stage >= 3) { 2} else {NULL},
              title = if (HELP_values$stage >= 3) {HTML( "<b>Step 3 - Dataset exploration and variable selection</b>")} else { "Step 3 - Dataset exploration and variable selection"},
              icon =  if (HELP_values$stage >= 3) { icon("location-dot")} else {icon("ellipsis")} ,
              color = if (HELP_values$stage >= 3) { "pink"} else {"gray"},
              time =  if (HELP_values$stage >= 3) { "Now"} else {"Later"},
              # footer = "Here is the footer",
              if (HELP_values$stage >= 3) {                HTML(
                "In this tab we will define which sequences and variables will be used for the analysis. Variables with NAs/Infinite/Zero variance are automatically removed.samples and workspace will be used. You need to indicate:<br>
                      <ol>
                      <li>Select and filter the sequences.</li>
                      <li>Select and filter the variables.</li>
                      <li>Do you want to compare between groups of sequences? If so, which sequences and what type of comparison do you want to do?</li>
                      <li>Do the calculations for the Selection set (which will be shown in the 4.Feature exploration and 5.Clonal exploration tabs) or the Exploration set (which you can use to try combinations and compare with the Selection set)?</li>
                     </ol>")
              } else {NULL}
            ),
            bs4Dash::timelineItem(
              elevation =  if (HELP_values$stage >= 3) { 2} else {NULL},
              title = if (HELP_values$stage >= 3) {HTML( "<b>Step 4 - Feature exploration</b>")} else { "Step 4 - Feature exploration"},
              icon =  if (HELP_values$stage >= 3) { icon("location-dot")} else {icon("ellipsis")} ,
              color = if (HELP_values$stage >= 3) { "pink"} else {"gray"},
              time =  if (HELP_values$stage >= 3) { "Now"} else {"Later"},
              # footer = "Here is the footer",
              if (HELP_values$stage >= 3) {                HTML(
                "In this tab you can plot the filtered variables as a violin plot and see how they behave according to the selected group in field 4. of 3.Dataset exploration and variable selection tab. You can add the values as dots and select those of interest and include the reconstructed germline sequences (if they are not already included) to compare how these values have changed since the germline sequence.")
              } else {NULL}
            ),
            bs4Dash::timelineItem(
              elevation =  if (HELP_values$stage >= 3) { 2} else {NULL},
              title = if (HELP_values$stage >= 3) {HTML( "<b>Step 5 - Clonal exploration</b>")} else { "Step 5 - Clonal exploration"},
              icon =  if (HELP_values$stage >= 3) { icon("location-dot")} else {icon("ellipsis")} ,
              color = if (HELP_values$stage >= 3) { "pink"} else {"gray"},
              time =  if (HELP_values$stage >= 3) { "Now"} else {"Later"},
              # footer = "Here is the footer",
              if (HELP_values$stage >= 3) {                HTML(
                "<p>In this tab, we will define which clonal definition will be used. You can:
                    <ul>
                      <li>Use a pre-existing clone definition (Clone_ID) or calculate it de novo.</li>
                      <li>Select a dominance threshold based on relative abundance.</li>
                      <li>See shared clones between samples.</li>
                      <li>Manually select clones of interest.</li>
                    </ul></p>
                    <p>If you have a clonal definition assigned and the images in the clonal and feature tabs have loaded, you will be able to export results in the Export Results tab.</p>")
              } else {NULL}
            ),
            if (HELP_values$stage == 4) { bs4Dash::timelineLabel("You can now export your work (it may take some time to allow it)", color = "pink")},
            bs4Dash::timelineStart(shiny::icon("hourglass-end", style = "color: rgb(61, 153, 112)"), color = "teal")
          )
        }
      )


      # 8.Background functions #### Disable menuitem when the app loads

      hideTab(inputId = "menutabset", target = "1.AIRR-Seq conversion")
      hideTab(inputId = "menutabset", target = "2.Sequence feature determination")
      hideTab(inputId = "menutabset", target = "Steps 1&2. Select your FBM and associated files")
      hideTab(
        inputId = "menutabset", target = "3.Dataset exploration and variable selection"
      )
      hideTab(inputId = "menutabset", target = "4.Feature exploration")
      hideTab(inputId = "menutabset", target = "5.Clonal exploration")
      hideTab(inputId = "menutabset", target = "5.Cell exploration")
      hideTab(inputId = "menutabset", target = "6.ML analysis")


      observeEvent(
        input$Move_to_1_airr, {
          HELP_values$stage <- 1
          hideTab(inputId = "menutabset", target = "0.Project information")
          showTab(inputId = "menutabset", target = "1.AIRR-Seq conversion", select = T)
        }
      )

      observeEvent(
        input$Move_to_2, {
          HELP_values$stage <- 2
          hideTab(inputId = "menutabset", target = "1.AIRR-Seq conversion")
          showTab(
            inputId = "menutabset", target = "2.Sequence feature determination",
            select = T
          )
          can_show_button_to_step$step_2 <- FALSE
          can_show_button_to_step$step_2_AIRR <- FALSE
        }
      )

      observeEvent(
        input$Move_to_analysis, {
          HELP_values$stage <- 2.5
          hideTab(inputId = "menutabset", target = "0.Project information")
          hideTab(inputId = "menutabset", target = "2.Sequence feature determination")
          showTab(
            inputId = "menutabset", target = "Steps 1&2. Select your FBM and associated files",
            select = T
          )
        }
      )

      observeEvent(
        input$Move_to_analysis_real, {
          HELP_values$stage <- 3
          hideTab(inputId = "menutabset", target = "0.Project information")
          hideTab(inputId = "menutabset", target = "2.Sequence feature determination")
          hideTab(inputId = "menutabset", target = "Steps 1&2. Select your FBM and associated files")
          showTab(
            inputId = "menutabset", target = "3.Dataset exploration and variable selection",
            select = T
          )
          showTab(inputId = "menutabset", target = "4.Feature exploration")
          showTab(inputId = "menutabset", target = "5.Clonal exploration")
          if (Selection_values$Cell) {
            showTab(inputId = "menutabset", target = "5.Cell exploration")
          }

          showTab(inputId = "menutabset", target = "6.ML analysis")
          if (!all(is.numeric(unlist(input$FBM_folder)))) {
            folder_values$Featured <- shinyFiles::parseDirPath(volumes(), input$FBM_folder)
          }
          sample_info_react$step_3=T
        }
      )

      # observeEvent(
      #     input$Reset_all, {
      #         session$reload()
      #         rm(list = names(input))
      #         HELP_values$stage <- 0
      #          showTab(inputId = "menutabset", target = "0.Project information", select = T)
      #         hideTab(inputId = "menutabset", target = "1.AIRR-Seq conversion")
      #         hideTab(inputId = "menutabset", target = "2.Sequence feature determination")
      #         hideTab(
      #             inputId = "menutabset", target = "3.Dataset exploration and variable selection"
      #         )
      #         hideTab(inputId = "menutabset", target = "4.Feature exploration")
      #         hideTab(inputId = "menutabset", target = "5.Clonal exploration")
      #         hideTab(inputId = "menutabset", target = "5.Cell exploration")
      #         hideTab(inputId = "menutabset", target = "6.ML analysis")
      #         can_show_button_to_step$step_2 <- FALSE
      #         can_show_button_to_step$step_2_AIRR <- FALSE
      #         can_show_button_to_step$step_3 <- FALSE
      #         folder_values$Featured <- ""
      #         folder_values$AIRR_parsed <- ""
      #         Big_mem_values$Header <- NULL
      #         Big_mem_values$Short_DF <- NULL
      #         Big_mem_values$Big_DF <- NULL
      #         Big_mem_values$categorical_newfields <-NULL
      #         Big_mem_values$categorical_missing_ids <- NULL
      #         Exploration_values$rows <- NULL
      #         Exploration_values$columns <- NULL
      #         Exploration_values$Scores <- NULL
      #         Exploration_values$Variance_explained <- NULL
      #         Exploration_values$Variance <- NULL
      #         Exploration_values$UMAP <- NULL
      #         Exploration_values$Parameters <- NULL
      #         Selection_values$rows <- NULL
      #         Selection_values$columns <- NULL
      #         Selection_values$Scores <- NULL
      #         Selection_values$Variance_explained <- NULL
      #         Selection_values$UMAP <- NULL
      #         Selection_values$Parameters <- NULL
      #         Selection_values$Cell <- F
      #         ID_selected_values$subclones <- NULL
      #         ID_selected_values$clones <- NULL
      #         ID_selected_values$clones_subclones_id <- NULL
      #         sample_info_react$test_status <- F
      #         sample_info_react$step_3=F
      #         updateNumericInput(session,"num",value=1)
      #
      #     }
      # )

      #print(session$userData)


      session$allowReconnect(TRUE)
    }


    if (!verbose) {
      withCallingHandlers(
        server_logic(),
        warning = function(w) invokeRestart("muffleWarning"),
        message = function(m) invokeRestart("muffleMessage")
      )
    } else {
      server_logic()
    }

}
