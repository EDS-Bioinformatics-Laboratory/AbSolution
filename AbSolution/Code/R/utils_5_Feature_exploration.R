#' 5_Feature_exploration
#'
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#' @import plotly
#' @noRd
draw_feature_violinplot <- function(
    values, name_values, sequence_info_df, group_info, additional_group_info = "Patient",
    show_reconstructed, selected_rows, selected_subclones = NULL, hide_dots = FALSE
) {
    print("Feature violin")
    # devtools::install_github('karthik/wesanderson') pal <-
    # wes_palette(length(unique(seq_df[,..group]))*2, name = 'Zissou1', type =
    # 'continuous')
    fig <- plot_ly(type = "violin")

    if (is.null(group_info) ||
        group_info == "") {
        group_info <- "Group"
    }
    if (show_reconstructed || length(unique(sequence_info_df$Seq_type[selected_rows])) >
        2) {
        i <- 0
        print("Reconstr")
        if (show_reconstructed && !("Reconstructed germline" %in% unique(sequence_info_df$Seq_type[selected_rows]))) {
            selected_rows <- sort(c(selected_rows, selected_rows - 1))
        } else if (show_reconstructed && !("Repertoire" %in% unique(sequence_info_df$Seq_type[selected_rows]))) {
            selected_rows <- sort(c(selected_rows, selected_rows + 1))
        }

        for (single_group in unique(sequence_info_df[selected_rows, get(group_info)])) {
            print(single_group)
            index <- intersect(
                selected_rows, which(
                  sequence_info_df[, get(group_info)] ==
                    single_group
              )
            )
            subset_seq_df <- sequence_info_df[intersect(
                selected_rows, which(
                  sequence_info_df[, get(group_info)] ==
                    single_group
              )
            ),
                ]
            subset_group <- sequence_info_df[intersect(
                selected_rows, which(
                  sequence_info_df[, get(group_info)] ==
                    single_group
              )
            ),
                get("additional_group_info")]

            for (seq_type in c("Repertoire", "Reconstructed_germline")) {
                i <- i + 1
                print(seq_type)
                sub_index <- intersect(index, which(sequence_info_df$Sequence_type == seq_type))
                subsubset_seq_df <- subset_seq_df[which(subset_seq_df$Sequence_type == seq_type),
                  ]
                subsubset_group <- subset_group[which(subset_seq_df$Sequence_type == seq_type)]
                if (nrow(subsubset_seq_df) >
                  0) {

                  # if(seq_type =='Reconstructed_germline') {
                  # print(head(subsubset_seq_df)) }
                  tmp_summary <- data.frame(Seq_ID = subsubset_seq_df$ID, Value = values[sub_index])


                  print("MMmmmMMm")
                  print(selected_subclones)
                  if (any(tmp_summary$Seq_ID %in% selected_subclones)) {
                    print(tmp_summary$Seq_ID[(which(tmp_summary$Seq_ID %in% selected_subclones))])
                  }
                  if (any(tmp_summary$Seq_ID %in% selected_subclones)) {
                    print("YESSSCLONES")
                  }
                  index_nonselected <- which(!(tmp_summary$Seq_ID %in% selected_subclones))
                  index_selected <- which((tmp_summary$Seq_ID %in% selected_subclones))
                  print("MMmmmMMm?")

                  fig <- add_trace(
                    fig, x = paste(single_group, sep = ""),
                    y = tmp_summary$Value[index_nonselected], hoveron = "points",
                    text = paste(
                      tmp_summary$Seq_ID[index_nonselected], " Value=", tmp_summary$Value[index_nonselected],
                      sep = ""
                  ),
                    key = tmp_summary$Seq_ID[index_nonselected], legendgroup = seq_type,
                    scalegroup = single_group, selectedpoints = if (any(tmp_summary$Seq_ID[index_nonselected] %in% selected_subclones)) {
                      ((which(tmp_summary$Seq_ID[index_nonselected] %in% selected_subclones) -
                        1))
                    } else {
                      NULL
                    }, name = seq_type, side = if (seq_type == "Reconstructed_germline") {
                      "positive"
                    } else {
                      "negative"
                    }, box = list(visible = T, line = list(color = "#1E2D24")),
                    points = "all", pointpos = if (seq_type == "Reconstructed_germline") {
                      0.5
                    } else {
                      -0.5
                    }, jitter = 0.25, scalemode = "count", meanline = list(visible = T),
                    color = I(
                      if (i%%2 != 0) {
                        "#2C497F"
                      } else {
                        "#808A9F"
                      }
                  ),
                    marker = list(
                      opacity = if (hide_dots) {
                        0
                      } else {
                        0.75
                      }, line = list(
                        width = 2, color = if (i%%2 != 0) {
                          "#2C497F"
                        } else {
                          "#808A9F"
                        }
                    )
                  ),
                    unselected = list(
                      marker = list(
                        size = 2, opacity = if (hide_dots) {
                          0
                        } else {
                          0.1
                        }, color = if (i%%2 != 0) {
                          "#2C497F"
                        } else {
                          "#808A9F"
                        }
                    )
                  ),
                    selected = list(
                      marker = list(
                        size = 10, opacity = if (hide_dots) {
                          0.5
                        } else {
                          1
                        }, color = "#000000", zorder = 2
                    )
                  ),
                    showlegend = if (i > 2) {
                      F
                    } else {
                      T
                    }
                )



                  if (length(index_selected) !=
                    0) {
                    print("Selected")
                    fig <- add_trace(
                      fig, x = paste(single_group, sep = ""),
                      y = tmp_summary$Value[index_selected], hoveron = "points",
                      text = paste(
                        tmp_summary$Seq_ID[index_selected], " Value=", tmp_summary$Value[index_selected],
                        sep = ""
                    ),
                      key = tmp_summary$Seq_ID[index_selected], legendgroup = seq_type,
                      scalegroup = single_group, selectedpoints = if (any(tmp_summary$Seq_ID[index_selected] %in% selected_subclones)) {
                        ((which(tmp_summary$Seq_ID[index_selected] %in% selected_subclones) -
                          1))
                      } else {
                        NULL
                      }, name = seq_type, side = if (seq_type == "Reconstructed_germline") {
                        "positive"
                      } else {
                        "negative"
                      }, box = list(visible = T, line = list(color = "#1E2D24")),
                      points = "all", pointpos = if (seq_type == "Reconstructed_germline") {
                        0.5
                      } else {
                        -0.5
                      }, jitter = 0.25, scalemode = "count", meanline = list(visible = T),
                      color = I(
                        if (i%%2 != 0) {
                          "#2C497F"
                        } else {
                          "#808A9F"
                        }
                    ),
                      marker = list(
                        opacity = if (hide_dots) {
                          0
                        } else {
                          0.75
                        }, line = list(
                          width = 2, color = if (i%%2 != 0) {
                            "#2C497F"
                          } else {
                            "#808A9F"
                          }
                      )
                    ),
                      unselected = list(
                        marker = list(
                          size = 2, opacity = if (hide_dots) {
                            0
                          } else {
                            0.1
                          }, color = if (i%%2 != 0) {
                            "#2C497F"
                          } else {
                            "#808A9F"
                          }
                      )
                    ),
                      selected = list(
                        marker = list(
                          size = 10, opacity = if (hide_dots) {
                            0.5
                          } else {
                            1
                          }, color = "#000000", zorder = 2
                      )
                    ),
                      showlegend = F
                  )
                  }

                }

            }

        }

    } else {
        i <- 0
        print("Repertoire")
        for (single_group in unique(sequence_info_df[, get(group_info)])) {
            print(single_group)
            index <- intersect(
                selected_rows, which(
                  sequence_info_df[, get(group_info)] ==
                    single_group
              )
            )
            subset_seq_df <- sequence_info_df[intersect(
                selected_rows, which(
                  sequence_info_df[, get(group_info)] ==
                    single_group
              )
            ),
                ]
            subset_group <- sequence_info_df[intersect(
                selected_rows, which(
                  sequence_info_df[, get(group_info)] ==
                    single_group
              )
            ),
                get("additional_group_info")]

            tmp_summary <- data.frame(Seq_ID = subset_seq_df$ID, Value = values[index])

            index_nonselected <- which(!(tmp_summary$Seq_ID %in% selected_subclones))
            index_selected <- which((tmp_summary$Seq_ID %in% selected_subclones))
            fig <- add_trace(
                fig, x = paste(single_group, sep = ""),
                y = tmp_summary$Value[index_nonselected], hoveron = "points", text = paste(
                  tmp_summary$Seq_ID[index_nonselected], " Value=", tmp_summary$Value[index_nonselected],
                  sep = ""
              ),
                key = tmp_summary$Seq_ID[index_nonselected], legendgroup = single_group,
                scalegroup = single_group, selectedpoints = if (any(tmp_summary$Seq_ID[index_nonselected] %in% selected_subclones)) {
                  ((which(tmp_summary$Seq_ID[index_nonselected] %in% selected_subclones) -
                    1))
                } else {
                  NULL
                }, name = single_group, box = list(visible = T, line = list(color = "#1E2D24")),
                points = "all", pointpos = 0, jitter = 0.25, scalemode = "count",
                meanline = list(visible = T),
                color = I(
                  if (i%%2 != 0) {
                    "#2C497F"
                  } else {
                    "#808A9F"
                  }
              ),
                marker = list(
                  opacity = if (hide_dots) {
                    0
                  } else {
                    0.75
                  }, line = list(
                    width = 2, color = if (i%%2 != 0) {
                      "#2C497F"
                    } else {
                      "#808A9F"
                    }
                )
              ),
                unselected = list(
                  marker = list(
                    size = 2, opacity = if (hide_dots) {
                      0
                    } else {
                      0.1
                    }, color = if (i%%2 != 0) {
                      "#2C497F"
                    } else {
                      "#808A9F"
                    }
                )
              ),
                selected = list(
                  marker = list(
                    size = 10, opacity = if (hide_dots) {
                      0.5
                    } else {
                      1
                    }, color = "#000000", zorder = 2
                )
              ),
                showlegend = F
            )



            if (length(index_selected) !=
                0) {
                print("Selected")

                fig <- add_trace(
                  fig, x = paste(single_group, sep = ""),
                  y = tmp_summary$Value[index_selected], hoveron = "points", text = paste(
                    tmp_summary$Seq_ID[index_selected], " Value=", tmp_summary$Value[index_selected],
                    sep = ""
                ),
                  key = tmp_summary$Seq_ID[index_selected], legendgroup = single_group,
                  scalegroup = single_group, selectedpoints = if (any(tmp_summary$Seq_ID[index_selected] %in% selected_subclones)) {
                    ((which(tmp_summary$Seq_ID[index_selected] %in% selected_subclones) -
                      1))
                  } else {
                    NULL
                  }, name = single_group, box = list(visible = T, line = list(color = "#1E2D24")),
                  points = "all", pointpos = 0, jitter = 0.25, scalemode = "count",
                  meanline = list(visible = T),
                  color = I(
                    if (i%%2 != 0) {
                      "#2C497F"
                    } else {
                      "#808A9F"
                    }
                ),
                  marker = list(
                    opacity = if (hide_dots) {
                      0
                    } else {
                      0.75
                    }, line = list(
                      width = 2, color = if (i%%2 != 0) {
                        "#2C497F"
                      } else {
                        "#808A9F"
                      }
                  )
                ),
                  unselected = list(
                    marker = list(
                      size = 2, opacity = if (hide_dots) {
                        0
                      } else {
                        0.1
                      }, color = if (i%%2 != 0) {
                        "#2C497F"
                      } else {
                        "#808A9F"
                      }
                  )
                ),
                  selected = list(
                    marker = list(
                      size = 10, opacity = if (hide_dots) {
                        0.5
                      } else {
                        1
                      }, color = "#000000", zorder = 2
                  )
                ),
                  showlegend = F
              )
            }


        }

    }

    fig <- layout(
        fig, title = paste(
            c(
                "Violin plot of the <i>feature ", name_values, "</i> <br> between the groups ",
                group_info
            ),
            collapse = ""
        ),
        violingap = 0, violingroupgap = 0, violinmode = "overlay", legend = list(tracegroupgap = 0),
        dragmode = "lasso"
    ) %>%
        config(
            displaylogo = FALSE, modeBarButtonsToRemove = c("autoScale2d", "pan2d", "hoverCompareCartesian")
        )

    fig <- fig %>%
        layout(plot_bgcolor = "rgba(0, 0, 0, 0)", paper_bgcolor = "rgba(0, 0, 0, 0)") %>%
        toWebGL()

    return(fig)
}
