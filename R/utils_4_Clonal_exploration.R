#' 4_Clonal_exploration
#'
#' @description **Internal function.** Not intended for direct use. Exported only for
#'    `shinymeta` report rendering via `::` access. Use [run_app()] instead.
#'
#' @return The return value, if any, from executing the utility.
#'
#' @keywords internal
#' @export
#' @examples
#' \donttest{
#'   # Internal function exported for shinymeta :: access during report rendering.
#'   # Requires a live Shiny reactive context and real AIRR-seq data.
#'   # Use run_app() as the user-facing entry point.
#' }
calculate_clone <- function(
    seq_df, clonotype, AA_or_NT = "NT", region = "CDR3", percentage = 100, calculate_shared_clones
) {

    if(clonotype == "Reconstructed_germline"){
      region="Whole"
    }

    if (!(paste(
      "Clone", AA_or_NT, region, clonotype, percentage, if (calculate_shared_clones) {
        "shared"
      }else {"non-shared"}, "simil", sep = "_"
    ) %in%
        colnames(seq_df))) {
        # print("calculate_clones")
        if (clonotype == "VCDR3J") {

            # tmp=paste( seq_df[selected_rows,'V_and_J'],
            # seq_df[selected_rows,paste(AA_or_NT,
            # 'CDR3',sep='_')],sep='_______')

            tmp <- paste(
                seq_df[, get("V_and_J")],
                seq_df[, get(paste(AA_or_NT, "CDR3", sep = "_"))],
                sep = "_______"
            )



        }

        if (clonotype == "Reconstructed_germline") {
            # tmp_selected_rows=sapply(selected_rows, function(z)
            # intersect(which(seq_df$Sequence_type=='Reconstructed_germline'),
            # which(seq_df[,'ID'] %in% seq_df[z,'ID'])) ) tmp=paste(
            # seq_df[tmp_selected_rows,'V_and_J'],
            # seq_df[tmp_selected_rows,paste(AA_or_NT,
            # 'Whole',sep='_')],sep='_______')

            tmp_selected_rows <- sapply(
                c(1:nrow(seq_df)),
                function(z) intersect(
                  which(seq_df$Sequence_type == "Reconstructed_germline"),
                  which(
                    seq_df[, get("ID")] %in%
                      seq_df[z, get("ID")]
                )
              )
            )
            tmp <- paste(
                seq_df[tmp_selected_rows, get("V_and_J")],
                seq_df[tmp_selected_rows, get(paste(AA_or_NT, region, sep = "_"))],
                sep = "_______"
            )

        }

        if (!calculate_shared_clones) {
            tmp <- paste(tmp, seq_df$Patient_Sample, sep = "__")
        }

        seq_df[, paste(
            "Clone", AA_or_NT, region, clonotype, percentage, if (calculate_shared_clones) {
                "shared"
            }else {"non-shared"}, "simil", sep = "_"
        )] <- tmp
        # print("calculate_clones_DONE")
        # print(
        #     paste(
        #         "Clone", AA_or_NT, region, clonotype, percentage, if (calculate_shared_clones) {
        #           "shared"
        #         } else {"non-shared"}, "simil", sep = "_"
        #     )
        # )

    }



    return(seq_df)
}

#' 4_Clonal_exploration
#'
#' @description **Internal function.** Not intended for direct use. Exported only for
#'    `shinymeta` report rendering via `::` access. Use [run_app()] instead.
#'
#' @return The return value, if any, from executing the utility.
#' @import plotly
#' @keywords internal
#' @export
#' @examples
#' \donttest{
#'   # Internal function exported for shinymeta :: access during report rendering.
#'   # Requires a live Shiny reactive context and real AIRR-seq data.
#'   # Use run_app() as the user-facing entry point.
#' }
draw_violinplots <- function(
    seq_df, group = "Patient_Sample", selected_rows, clonotype, AA_or_NT = "AA", region = "CDR3",
    percentage = 100, freq_filter = 0, Selected_clones = NULL, dominance_threshold, seed=1234,
    really_hide_dots=FALSE, width=1400, height=1000,img_type= "png", scale=4,
    source="clone_violinplot"
) {
    set.seed(seed)
    # devtools::install_github('karthik/wesanderson') pal <-
    # wes_palette(length(unique(seq_df[,..group]))*2, name = 'Zissou1', type =
    # 'continuous')
    fig <- plot_ly(type = "violin") %>%
      config(
        toImageButtonOptions = list(
          format = img_type,
          filename = paste("Violin_plot_clones",clonotype, Sys.time(),sep="_"),
          width = width,
          height = height,
          scale=scale
        )
      )

    i <- 0
    min_value=0
    max_value=0
    for (single_group in unique(seq_df[selected_rows, get(group)])) {

        # subset_seq_df=seq_df[intersect(selected_rows,
        # which(seq_df[,get(group)]==single_group)),]
        # subset_group=unique(seq_df[intersect(selected_rows,
        # which(seq_df[,get(group)]==single_group)), get('Group')])

        og_subset_seq_df <- seq_df[intersect(
            selected_rows, which(
                seq_df[, get(group)] ==
                  single_group
            )
        ),
            ]
        subset_seq_df <- seq_df[intersect(
            which(seq_df$Sequence_type %in% unique(seq_df$Sequence_type[selected_rows])),
            which(
                seq_df[, get(group)] ==
                  single_group
            )
        ),
            ]
        subset_group <- unique(
            seq_df[intersect(
                which(seq_df$Sequence_type %in% unique(seq_df$Sequence_type[selected_rows])),
                which(
                  seq_df[, get(group)] ==
                    single_group
              )
            ),
                get("Group")]
        )
        # print(single_group)


        for (seq_type in c("Repertoire", "Reconstructed_germline")) {
            i <- i + 1
            # print(seq_type)
            subsubset_seq_df <- subset_seq_df[which(subset_seq_df$Sequence_type == seq_type),
                ]
            og_subsubset_seq_df <- og_subset_seq_df[which(og_subset_seq_df$Sequence_type == seq_type),
                ]
            if (nrow(subsubset_seq_df) >
                0) {

                # if (seq_type == "Reconstructed_germline") {
                #   print(subsubset_seq_df)
                # }
                tmp_clones <- table(subsubset_seq_df[, get(clonotype)])

                tmp_clones=tmp_clones[order(names(tmp_clones))]


                tmp_subclones <- table(
                  paste(
                    subsubset_seq_df[, get(clonotype)],
                    subsubset_seq_df[, get("NT_Whole")],
                    sep = "_&_"
                )


              )
                tmp_subclones_num=strsplit(names(tmp_subclones), split="_&_")
                tmp_subclones_num=unlist(tmp_subclones_num)[2*(1:length(tmp_subclones))-1]
                tmp_subclones_num=table(tmp_subclones_num)
                tmp_subclones_num=tmp_subclones_num[order(names(tmp_subclones_num))]

                tmp_frequencies <- data.frame(
                  Clone_ID = names(tmp_clones),
                  Frequency = (100 * as.numeric(tmp_clones))/(sum(tmp_clones)),
                #   Number_of_unique_subclones = sapply(
                #     names(tmp_clones),
                #     function(z) length(
                #       which(
                #         startsWith(
                #           names(tmp_subclones),
                #           paste(z, "_&_", sep = "")
                #       )
                #     )
                #   )
                # ),
                Number_of_unique_subclones = as.numeric(unname(tmp_subclones_num)),
                  Number_of_sequences = as.numeric(tmp_clones)
              )
                tmp_frequencies <- tmp_frequencies[which(tmp_frequencies$Clone_ID %in% unique(og_subsubset_seq_df[, get(clonotype)])),
                  ]

                max_value=max(max_value, log10(max(tmp_frequencies$Frequency)))
                min_value=min(min_value, log10(min(tmp_frequencies$Frequency)))

                # print(Selected_clones$key)
                # if (any(tmp_frequencies$Clone_ID %in% Selected_clones$key)) {
                #   print("Selected clones:")
                #   print(
                #     tmp_frequencies$Clone_ID[(which(tmp_frequencies$Clone_ID %in% Selected_clones$key))]
                # )
                # }
                # if (any(tmp_frequencies$Clone_ID %in% Selected_clones$key)) {
                #   print(
                #     list(
                #       (which(tmp_frequencies$Clone_ID %in% Selected_clones$key) -
                #         1)
               #    )
               #  )
               #  }
                # if (any(tmp_frequencies$Clone_ID %in% Selected_clones$key)) {
                #   print("YESSSCLONES")
                # }
                # print("MMmmmMMm?")

                # FIX #### writing tmp_frequencies to a file and checking later

                # saving also everything in a RData???
                fig <- add_trace(
                  fig, x = paste(
                    single_group, " (", paste(subset_group, collapse = "&"),
                    ")", sep = ""
                ),
                  y = tmp_frequencies$Frequency[which(!(tmp_frequencies$Clone_ID %in% Selected_clones$key))],
                  hoveron = "points", text = paste(
                    tmp_frequencies$Clone_ID[which(!(tmp_frequencies$Clone_ID %in% Selected_clones$key))],
                    " Number of unique Fab regions=", tmp_frequencies$Number_of_unique_subclones[which(!(tmp_frequencies$Clone_ID %in% Selected_clones$key))],
                    " Number of sequences=", tmp_frequencies$Number_of_sequences[which(!(tmp_frequencies$Clone_ID %in% Selected_clones$key))],
                    sep = ""
                ),
                  key = tmp_frequencies$Clone_ID[which(!(tmp_frequencies$Clone_ID %in% Selected_clones$key))],
                  legendgroup = seq_type, scalegroup = single_group, selectedpoints = if (any(
                    tmp_frequencies$Clone_ID[which(!(tmp_frequencies$Clone_ID %in% Selected_clones$key))] %in%
                      Selected_clones$key
                )) {
                    ((which(
                      tmp_frequencies$Clone_ID[which(!(tmp_frequencies$Clone_ID %in% Selected_clones$key))] %in%
                        Selected_clones$key
                  ) -
                      1))
                  } else {
                    NULL
                  }, name = seq_type, side = if (seq_type == "Reconstructed_germline") {
                    "positive"
                  } else {
                    "negative"
                  }, box = list(visible = TRUE, line = list(color = "#2D2926")),
                  points = if(really_hide_dots){FALSE}else{"all"}, pointpos = if (seq_type == "Reconstructed_germline") {
                    0.5
                  } else {
                    -0.5
                  }, jitter = 0.25, scalemode = "count", meanline = list(visible = T),
                  color = I(
                    if (i%%2 != 0) {
                      "#E85D75"
                    } else {
                      "#707A6C"
                    }
                ),
                  marker = list(
                    opacity = 0.75, line = list(
                      width = 2, color = if (i%%2 != 0) {
                        "#E85D75"
                      } else {
                        "#707A6C"
                      }
                  )
                ),
                  unselected = list(
                    marker = list(
                      size = 2, opacity = 0.1, color = if (i%%2 != 0) {
                        "#E85D75"
                      } else {
                        "#707A6C"
                      }
                  )
                ),
                  selected = list(marker = list(size = 10, opacity = 1, color = "#FFA500", zorder = 2)),
                  showlegend = if (i > 2) {
                    F
                  } else {
                    T
                  }
              )

                if (length(which(tmp_frequencies$Clone_ID %in% Selected_clones$key)) !=
                  0) {

                  fig <- add_trace(
                    fig, x = paste(
                      single_group, " (", paste(subset_group, collapse = "&"),
                      ")", sep = ""
                  ),
                    y = tmp_frequencies$Frequency[which(tmp_frequencies$Clone_ID %in% Selected_clones$key)],
                    hoveron = "points", text = if (length(which(tmp_frequencies$Clone_ID %in% Selected_clones$key)) ==
                      0) {
                      NULL
                    } else {
                      paste(
                        tmp_frequencies$Clone_ID[which(tmp_frequencies$Clone_ID %in% Selected_clones$key)],
                        " Number of unique Fab regions=", tmp_frequencies$Number_of_unique_subclones[which(tmp_frequencies$Clone_ID %in% Selected_clones$key)],
                        " Number of sequences=", tmp_frequencies$Number_of_sequences[which(tmp_frequencies$Clone_ID %in% Selected_clones$key)],
                        sep = ""
                    )
                    }, key = tmp_frequencies$Clone_ID[which(tmp_frequencies$Clone_ID %in% Selected_clones$key)],
                    legendgroup = seq_type, scalegroup = single_group, selectedpoints = if (any(
                      tmp_frequencies$Clone_ID[which(tmp_frequencies$Clone_ID %in% Selected_clones$key)] %in%
                        Selected_clones$key
                  )) {
                      ((which(
                        tmp_frequencies$Clone_ID[which(tmp_frequencies$Clone_ID %in% Selected_clones$key)] %in%
                          Selected_clones$key
                    ) -
                        1))
                    } else {
                      NULL
                    }, name = seq_type, side = if (seq_type == "Reconstructed_germline") {
                      "positive"
                    } else {
                      "negative"
                    }, box = list(visible = TRUE, line = list(color = "#2D2926")),
                    points = if(really_hide_dots){FALSE}else{"all"}, pointpos = if (seq_type == "Reconstructed_germline") {
                      0.2
                    } else {
                      -0.2
                    }, jitter = 0.03, scalemode = "count", meanline = list(visible = T),
                    color = I(
                      if (i%%2 != 0) {
                        "#E85D75"
                      } else {
                        "#707A6C"
                      }
                  ),
                    marker = list(
                      opacity = 0.75, size=10,zorder = 2,
                      line = list(
                        width = 4,color = "#FFA500"
                      )
                  ),
                  #   unselected = list(
                  #     marker = list(
                  #       size = 2, opacity = 0.1, color = if (i%%2 != 0) {
                  #         "#A8577E"
                  #       } else {
                  #         "#707A6C"
                  #       }
                  #   )
                  # ),
                  #   selected = list(marker = list(size = 10, opacity = 1, color = "#000000", zorder = 2)),
                    showlegend = FALSE
                )
                }
            }


        }

    }

    fig <- layout(
        fig, title = paste(
            c(
                "Relative log10 frequency distribution of clones<br><i>using the criteria ",
                clonotype, "</i>"
            ),
            collapse = ""
        ),
        yaxis = list(
          title = "Log10 clonal frequency",
          type = "log", exponentformat = "power", showexponent = "all",
          range = c( min_value-0.1,
                    max(2, max_value+0.1))
        ),
        # violingap = 0, violingroupgap = 0, violinmode = "overlay",
        legend = list(tracegroupgap = 0),
        dragmode = "lasso", shapes = list(hline((dominance_threshold))),
        xaxis = list(showline = TRUE,
                     zeroline = FALSE),
        yaxis = list(showline = TRUE,
                     zeroline = FALSE)
    ) %>%
        config(
            displaylogo = FALSE, modeBarButtonsToRemove = c("autoScale2d", "pan2d", "hoverCompareCartesian")
        )


    fig <- suppressWarnings(fig %>%
        layout(plot_bgcolor = "rgba(0, 0, 0, 0)", paper_bgcolor = "rgba(0, 0, 0, 0)",
               xaxis = list(showline = TRUE,
                            zeroline = FALSE),
               yaxis = list(showline = TRUE,
                            zeroline = FALSE)) %>%
        toWebGL())

    return(fig)
}

#' 4_Clonal_exploration
#'
#' @description **Internal function.** Not intended for direct use. Exported only for
#'    `shinymeta` report rendering via `::` access. Use [run_app()] instead.
#'
#' @return The return value, if any, from executing the utility.
#' @import upsetjs
#' @keywords internal
#' @export
#' @examples
#' \donttest{
#'   # Internal function exported for shinymeta :: access during report rendering.
#'   # Requires a live Shiny reactive context and real AIRR-seq data.
#'   # Use run_app() as the user-facing entry point.
#' }
draw_upsetplot <- function(
    seq_df, group = "Patient", selected_rows, clonotype, AA_or_NT = "AA", region = "CDR3",
    percentage = 100, freq_filter = 0, Selected_clones = NULL
) {

    set_list <- list()

    for (single_group in unique(seq_df[, get(group)])) {
        subset_seq_df <- seq_df[intersect(
            selected_rows, which(
                seq_df[, get(group)] ==
                  single_group
            )
        ),
            ]
        subset_group <- unique(
            seq_df[intersect(
                selected_rows, which(
                  seq_df[, get(group)] ==
                    single_group
              )
            ),
                get("Group")]
        )
        # print(single_group)
        subsubset_seq_df <- subset_seq_df[which(subset_seq_df$Sequence_type == "Repertoire"),
            ]
        tmp_clones <- table(subsubset_seq_df[, get(clonotype)])
        set_list[[length(set_list) +
            1]] <- unique(names(tmp_clones))

    }

    names(set_list) <- unique(seq_df[, get(group)])



    upset_plot <- upsetjs::upsetjs() %>%
        upsetjs::fromList(set_list) %>%
        chartStyleFlags(export.buttons =T)

    return(upset_plot)
}

#' 4_Clonal_exploration
#'
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#'
#' @noRd
hline <- function(y = 0, color = "#2D2926") {
    list(
        type = "line", x0 = 0, x1 = 100, xref = "paper", y0 = y, y1 = y, line = list(color = color, width = 0.75, dash = "dot")
    )
}

#' 4_Clonal_exploration
#'
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#'
#' @noRd
vline <- function(x = 0, color = "#2D2926") {
    list(
        type = "line", x0 = x, x1 = x, xref = "paper", y0 = 0, y1 = 100, line = list(color = color, width = 0.75, dash = "dot")
    )
}

#' 4_Clonal_exploration
#'
#' @description **Internal function.** Not intended for direct use. Exported only for
#'    `shinymeta` report rendering via `::` access. Use [run_app()] instead.
#'
#' @return The return value, if any, from executing the utility.
#' @import plotly
#' @import utils
#' @keywords internal
#' @export
#' @examples
#' \donttest{
#'   # Internal function exported for shinymeta :: access during report rendering.
#'   # Requires a live Shiny reactive context and real AIRR-seq data.
#'   # Use run_app() as the user-facing entry point.
#' }
draw_sharedclonesplot <- function(
    seq_df, sets = NULL, group = "Patient", selected_rows, clonotype, AA_or_NT = "AA",
    region = "CDR3", percentage = 100, freq_filter = 0, Selected_clones = NULL, dominance_threshold,
    color_by="Clonal sharing", colorblind, Big_mem_color_values,
    seed=1234, width=1400, height=1000,img_type= "png",  scale=4
) {
    set.seed(seed)
    list_figs <- list()
    set_vector <- c()
    Shared_intersected_clones=c()

    if(color_by!= "Clonal sharing" && color_by !="VJ usage") {
      colors <- AbSolution::Ab_palette(
        list_values = unique(seq_df[,get(color_by)]),
        vect_genes_comb = NA, type_values = "cualitative", colorblind = colorblind,
        seed=seed
      )
      names(colors)=unique(seq_df[,get(color_by)])
    }

    for (group_A in unlist(sets)) {
      for (group_B in unlist(sets)) {
        if (group_A != group_B) {
          tmp_comparison <- paste(
            sort(c(group_A, group_B)),
            collapse = "    vs    "
          )
          if (!(tmp_comparison %in% set_vector)) {
            set_vector <- c(set_vector, tmp_comparison)
            group1 <- sort(c(group_A, group_B))[1]
            group2 <- sort(c(group_A, group_B))[2]
            subset1_seq_df <- seq_df[intersect(
              selected_rows, which(
                seq_df[, get(group)] ==
                  group1
              )
            ),
            ]
            subset1_seq_df <- subset1_seq_df[which(subset1_seq_df$Sequence_type == "Repertoire"),
            ]
            subset2_seq_df <- seq_df[intersect(
              selected_rows, which(
                seq_df[, get(group)] ==
                  group2
              )
            ),
            ]
            subset2_seq_df <- subset2_seq_df[which(subset2_seq_df$Sequence_type == "Repertoire"),
            ]

            if(length(Shared_intersected_clones)==0) {
              Shared_intersected_clones <- intersect(unique(subset1_seq_df[, get(clonotype)]),
                                                     unique(subset2_seq_df[, get(clonotype)]))
            } else {
              Shared_intersected_clones <- intersect(Shared_intersected_clones,
                                                     intersect(unique(subset1_seq_df[, get(clonotype)]),
                                                               unique(subset2_seq_df[, get(clonotype)])))
            }


          }
        }
      }}

    set_vector <- c()
    first_iteration=T
    for (group_A in unlist(sets)) {
        for (group_B in unlist(sets)) {
            if (group_A != group_B) {
                tmp_comparison <- paste(
                  sort(c(group_A, group_B)),
                  collapse = "    vs    "
              )
                if (!(tmp_comparison %in% set_vector)) {
                  set_vector <- c(set_vector, tmp_comparison)
                  group1 <- sort(c(group_A, group_B))[1]
                  group2 <- sort(c(group_A, group_B))[2]
                  subset1_seq_df <- seq_df[intersect(
                    selected_rows, which(
                      seq_df[, get(group)] ==
                        group1
                  )
                ),
                    ]
                  subset1_seq_df <- subset1_seq_df[which(subset1_seq_df$Sequence_type == "Repertoire"),
                    ]
                  subset2_seq_df <- seq_df[intersect(
                    selected_rows, which(
                      seq_df[, get(group)] ==
                        group2
                  )
                ),
                    ]
                  subset2_seq_df <- subset2_seq_df[which(subset2_seq_df$Sequence_type == "Repertoire"),
                    ]

                  tmp_clones1 <- table(subset1_seq_df[, get(clonotype)])

                  tmp_frequencies1 <- data.frame(
                    Clone_ID = names(tmp_clones1),
                    Frequency = (100 * as.numeric(tmp_clones1))/((nrow(subset1_seq_df))),
                    Number_of_subclones = as.numeric(tmp_clones1)
                )

                  tmp_clones2 <- table(subset2_seq_df[, get(clonotype)])

                  tmp_frequencies2 <- data.frame(
                    Clone_ID = names(tmp_clones2),
                    Frequency = (100 * as.numeric(tmp_clones2))/((nrow(subset2_seq_df))),
                    Number_of_subclones = as.numeric(tmp_clones2)
                )

                  tmp_summary_df <- merge(tmp_frequencies1, tmp_frequencies2, by = "Clone_ID", all = T)
                  tmp_summary_df[(is.na(tmp_summary_df))] <- 0



                  tmp_summary_df$Freq_X <- sapply(
                    tmp_summary_df$Frequency.x, function(z) if (z ==
                      0) {
                      5e-04
                    } else {
                      z
                    }
                )
                  tmp_summary_df$Freq_Y <- sapply(
                    tmp_summary_df$Frequency.y, function(z) if (z ==
                      0) {
                      5e-04
                    } else {
                      z
                    }
                )

                  tmp_vector <- tmp_summary_df$Freq_X
                  for (tmp_value in unique(tmp_vector[which(duplicated(tmp_vector))])) {
                    tmp_summary_df$Freq_X[which(tmp_vector == tmp_value)] <- jitter(tmp_summary_df$Freq_X[which(tmp_vector == tmp_value)])
                  }

                  tmp_vector <- tmp_summary_df$Freq_Y
                  for (tmp_value in unique(tmp_vector[which(duplicated(tmp_vector))])) {
                    tmp_summary_df$Freq_Y[which(tmp_vector == tmp_value)] <- jitter(tmp_summary_df$Freq_Y[which(tmp_vector == tmp_value)])
                  }
                  tmp_summary_df$info <- paste(
                    paste(
                      tmp_summary_df$Clone_ID, paste(tmp_summary_df$Frequency.x, tmp_summary_df$Frequency.y, sep = " , "),
                      sep = " ("
                  ),
                    ")", sep = ""
                )
                  if(color_by == "Clonal sharing") {
                    tmp_summary_df$shared <- sapply(
                      1:nrow(tmp_summary_df),
                      function(z) if (tmp_summary_df$Frequency.x[z] !=
                                      0 && tmp_summary_df$Frequency.y[z] != 0) {
                        if(tmp_summary_df$Clone_ID[z] %in% Shared_intersected_clones) {
                          "Shared - intersection"
                        } else {
                          "Shared"
                        }

                      } else {
                        "Unique"
                      }
                    )
                  } else {
                    if (color_by =="VJ usage"){
                      columna_color="V_and_J"
                    } else {
                      columna_color=color_by
                    }
                    # tmp_summary_df$shared <- sapply(
                    #   tmp_summary_df$Clone_ID,
                    #   function(z) names(which.max(table(c(subset1_seq_df[subset1_seq_df[, get(clonotype)]==z, get(columna_color)],
                    #                                       subset2_seq_df[subset2_seq_df[, get(clonotype)]==z, get(columna_color)]))))
                    # )

                    DT1 <- data.table::as.data.table(subset1_seq_df)[, .(Clone_ID = get(clonotype), color = get(columna_color))]
                    DT2 <- data.table::as.data.table(subset2_seq_df)[, .(Clone_ID = get(clonotype), color = get(columna_color))]

                    DT <- data.table::rbindlist(list(DT1, DT2), use.names = TRUE)
                    DT <- DT[!is.na(color)]

                    cnt <- DT[, .N, by = .(Clone_ID, color)]

                    best <- cnt[order(Clone_ID, -N)][, .SD[1], by = Clone_ID]
                    tmp_summary_df$shared <- best$color[match(tmp_summary_df$Clone_ID, best$Clone_ID)]
                  }

                  # tmp_summary_df$selected <- rep("#2D2926", nrow(tmp_summary_df))
                  # tmp_summary_df$selected[which(startsWith(tmp_summary_df$shared, "Shared"))] <- "#E63946"
                  # tmp_summary_df$selected[which(tmp_summary_df$Clone_ID %in% Selected_clones)] <- "#FFBC0A"
                  tmp_summary_df$selected <- rep("rgba(0,0,0,0)", nrow(tmp_summary_df))
                  # tmp_summary_df$selected[which(startsWith(tmp_summary_df$shared, "Shared"))] <- "#E63946"
                  tmp_summary_df$selected[which(tmp_summary_df$Clone_ID %in% Selected_clones)] <- "#FFBC0A"


                  # if(first_iteration) {print("Iteration first")} else {print("No first")}

                  tmp_summary_df$shared <- factor(tmp_summary_df$shared, levels = sort(unique(tmp_summary_df$shared)))


                  # tmp_summary_df$selected
                  # tmp_summary_df$selected[which(tmp_summary_df$Clone_ID %in% Selected_clones)] <- "Selected"
                  # tmp_summary_df$selected <- factor(tmp_summary_df$selected, levels = c("Unselected", "Selected"))
                  # print(head(tmp_summary_df))
                  tmp_fig <- plot_ly()

                  for (sub in unique(tmp_summary_df$shared)){

                    color_plot=if(color_by=="Clonal sharing"){
                            if (sub=="Unique") {
                              "#2D2926"
                            } else if ( sub=="Shared - intersection") {
                              "#FFBC0A"
                            } else {
                              "#E63946"
                            }
                    } else if (color_by=="VJ usage") {
                      Big_mem_color_values$VJ[which(names(Big_mem_color_values$VJ)==sub)]
                    } else {
                      colors[which(names(colors)==sub)]
                    }
                    tmp_fig <- add_trace(tmp_fig,
                                         data = tmp_summary_df[which(tmp_summary_df$shared==sub),],
                                         type="scattergl", mode="markers",
                                         x = ~Freq_X, y = ~Freq_Y, text = ~info,
                                         hoverinfo = "text",
                                         key = tmp_summary_df$Clone_ID[which(tmp_summary_df$shared==sub)],
                                         legendgroup=sub,name=sub,
                                         showlegend = if(first_iteration) {T} else {F},
                                         marker = list(size = 15,opacity = 0.75,
                                                       color=color_plot,
                                                       line = list(width = 3,
                                                                   color = ~selected)))
                  }
                  tmp_fig%>%
                    config(
                      toImageButtonOptions = list(
                        format = img_type,
                        filename = paste("Shared_plot",clonotype,Sys.time(),sep="_"),
                        width = width,
                        height = height,
                        scale=scale
                      )
                    )
                  tmp_fig <- layout(
                    tmp_fig, xaxis = list(
                      title = paste("Log10 clonal frequency", group1, sep = " - "),
                      type = "log", exponentformat = "power",
                      range = c(min(-5, log10(min(tmp_summary_df$Freq_X))-0.1),
                                max(2, log10(max(tmp_summary_df$Freq_X))+0.1))
                  ),
                    yaxis = list(
                      title = paste("Log10 clonal frequency", group2, sep = " - "),
                      type = "log", exponentformat = "power", showexponent = "all",
                      range = c(min(-5, log10(min(tmp_summary_df$Freq_Y))-0.1),
                                max(2, log10(max(tmp_summary_df$Freq_Y)+0.1)))
                  ),
                    shapes = list(
                      hline((dominance_threshold)),
                      vline((dominance_threshold))
                  ),
                  xaxis = list(showline = TRUE,
                               zeroline = FALSE),
                  yaxis = list(showline = TRUE,
                               zeroline = FALSE)
                )
                  list_figs[[length(list_figs) +
                    1]] <- tmp_fig

                  first_iteration=F
                }
            }
        }


    }

    sharedclonesplot <- subplot(
        list_figs, nrows = ceiling(length(list_figs)/3),
        margin = 0.05, titleY = TRUE, titleX = TRUE
    )%>%
      config(
        toImageButtonOptions = list(
          format = img_type,
          filename = paste("Shared_clones_plot",Sys.time(),sep="_"),
          width = width,
          height = height,
          scale=scale
        )
      ) %>%
      layout( paper_bgcolor = "transparent",
              xaxis = list(showline = TRUE,
                           zeroline = FALSE),
              yaxis = list(showline = TRUE,
                           zeroline = FALSE))%>%
      config(
        displaylogo = FALSE,
        modeBarButtonsToRemove = c("zoomIn2d", "zoomOut2d",
                                   "autoScale2d",
                                   "hoverCompareCartesian")
      )

    return(sharedclonesplot)
}
