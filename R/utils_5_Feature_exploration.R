#' 5_Feature_exploration
#'
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#' @import plotly
#' @keywords internal
#' @export
draw_feature_violinplot <- function(
    values, name_values, sequence_info_df, group_info,
    additional_group_info = "Patient",
    show_reconstructed, compare_opposites=F, selected_rows,
    selected_subclones = NULL, selected_clones = NULL,hide_dots = FALSE,
    really_hide_dots=F, width=1400, height=1000, img_type= "png", scale=4,
    seed=1234, source="feature_violinplot", primary_color="#09BC8A"
) {
  set.seed(seed)
  # print("Feature violin")
  # devtools::install_github('karthik/wesanderson') pal <-
  # wes_palette(length(unique(seq_df[,..group]))*2, name = 'Zissou1', type =
  # 'continuous')
  fig <- plot_ly(type = "violin") %>%
    config(
      toImageButtonOptions = list(
        format = img_type,
        filename = paste("Violin_plot_features",name_values, Sys.time(),sep="_"),
        width = width,
        height = height,
        scale=scale
      )
    )%>%
    layout(
      xaxis = list(showline = TRUE,
                   zeroline = FALSE),
      yaxis = list(showline = TRUE,
                   zeroline = FALSE)
    )


  if (is.null(group_info) ||
      group_info == "") {
    group_info <- "Group"
  }
  if (show_reconstructed || length(unique(sequence_info_df$Seq_type[selected_rows])) >
      2) {
    i <- 0
    # print("Reconstr")
    if (show_reconstructed && !("Reconstructed germline" %in% unique(sequence_info_df$Seq_type[selected_rows]))) {
      selected_rows <- sort(c(selected_rows, selected_rows - 1))
    } else if (show_reconstructed && !("Repertoire" %in% unique(sequence_info_df$Seq_type[selected_rows]))) {
      selected_rows <- sort(c(selected_rows, selected_rows + 1))
    }

    if (compare_opposites) {
      meta_summary<- data.frame(Group=  as.character(),
                                Subgroup= as.character() ,
                                Seq_ID = as.character(),
                                Value = as.numeric())
    }

    for (single_group in sort(unique(sequence_info_df[selected_rows, get(group_info)]),decreasing = T)) {
      # print(single_group)
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
        # print(seq_type)
        sub_index <- intersect(index, which(sequence_info_df$Sequence_type == seq_type))
        subsubset_seq_df <- subset_seq_df[which(subset_seq_df$Sequence_type == seq_type),
        ]
        subsubset_group <- subset_group[which(subset_seq_df$Sequence_type == seq_type)]
        if (nrow(subsubset_seq_df) >
            0) {

          # if(seq_type =='Reconstructed_germline') {
          # print(head(subsubset_seq_df)) }
          tmp_summary <- data.frame(Group= single_group,
                                    Subgroup=seq_type ,
                                    Seq_ID = subsubset_seq_df$ID,
                                    Value = values[sub_index])

          if(compare_opposites) {
            meta_summary <- rbind(meta_summary, tmp_summary)
          }

          index_nonselected <- which(!(tmp_summary$Seq_ID %in% c(selected_subclones,selected_clones)))
          index_selected_sub <- which((tmp_summary$Seq_ID %in% selected_subclones))
          index_selected_clones <- which((tmp_summary$Seq_ID %in% selected_clones))
          index_selected<- sort(unique(c(index_selected_sub,index_selected_clones)))
          index_colors_sub<- sapply(index_selected, function(z) if (z %in%index_selected_clones && z %in% index_selected_sub) {"#C07100"} else if(z %in%index_selected_clones){"#FFA500"}else {"#000000"})

          fig <- add_trace(
            fig, x = paste(single_group, sep = ""),
            y = tmp_summary$Value[index_nonselected], hoveron = "points",
            text = paste(
              tmp_summary$Seq_ID[index_nonselected], " Value=", tmp_summary$Value[index_nonselected],
              sep = ""
            ),
            key = tmp_summary$Seq_ID[index_nonselected], legendgroup = seq_type,
            scalegroup = single_group, selectedpoints = if (any(tmp_summary$Seq_ID[index_nonselected] %in% c(selected_subclones, selected_clones))) {
              ((which(tmp_summary$Seq_ID[index_nonselected] %in% c(selected_subclones, selected_clones)) -
                  1))
            } else {
              NULL
            }, name = seq_type, side = if (seq_type == "Reconstructed_germline") {
              "positive"
            } else {
              "negative"
            }, box = list(visible = T, line = list(color = "#2D2926")),
            points = if(really_hide_dots){FALSE}else{"all"}, pointpos = if (seq_type == "Reconstructed_germline") {
              0.5
            } else {
              -0.5
            }, jitter = 0.25, scalemode = "count", meanline = list(visible = T),
            color = I(
              if (i%%2 != 0) {
                primary_color
              } else {
                "#707A6C"
              }
            ),
            marker = list(
              opacity = if (hide_dots) {
                0
              } else {
                0.75
              }, line = list(
                width = 2, color = if (i%%2 != 0) {
                  primary_color
                } else {
                  "#707A6C"
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
                  primary_color
                } else {
                  "#707A6C"
                }
              )
            ),
            selected = list(
              marker = list(
                size = 10, opacity = if (hide_dots) {
                  0.5
                } else {
                  1
                }, color = "#2D2926", zorder = 2
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
            fig <- add_trace(
              fig, x = paste(single_group, sep = ""),
              y = tmp_summary$Value[index_selected], hoveron = "points",
              text = paste(
                tmp_summary$Seq_ID[index_selected], " Value=", tmp_summary$Value[index_selected],
                sep = ""
              ),
              key = tmp_summary$Seq_ID[index_selected], legendgroup = seq_type,
              scalegroup = single_group, selectedpoints = if (any(tmp_summary$Seq_ID[index_selected]  %in% c(selected_subclones, selected_clones))) {
                ((which(tmp_summary$Seq_ID[index_selected] %in% c(selected_subclones, selected_clones)) -
                    1))
              } else {
                NULL
              }, name = seq_type, side = if (seq_type == "Reconstructed_germline") {
                "positive"
              } else {
                "negative"
              }, box = list(visible = T, line = list(color = "#2D2926")),
              points = if(really_hide_dots){FALSE}else{"all"}, pointpos = if (seq_type == "Reconstructed_germline") {
                0.5
              } else {
                -0.5
              }, jitter = 0.25, scalemode = "count", meanline = list(visible = T),
              color = I(
                if (i%%2 != 0) {
                  primary_color
                } else {
                  "#707A6C"
                }
              ),
              marker = list(
                opacity = if (hide_dots) {
                  0
                } else {
                  0.75
                }, line = list(
                  width = 2, color = if (i%%2 != 0) {
                    primary_color
                  } else {
                    "#707A6C"
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
                    primary_color
                  } else {
                    "#707A6C"
                  }
                )
              ),
              selected = list(
                marker = list(
                  size = 10, opacity = if (hide_dots) {
                    0.5
                  } else {
                    1
                  }, color = index_colors_sub, zorder = 2
                )
              ),
              showlegend = F
            )
          }

        }

      }

    }

  } else {
    i <- 1
    for (single_group in unique(sequence_info_df[, get(group_info)])) {
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

      index_nonselected <- which(!(tmp_summary$Seq_ID %in% c(selected_subclones,selected_clones)))
      index_selected_sub <- which((tmp_summary$Seq_ID %in% selected_subclones))
      index_selected_clones <- which((tmp_summary$Seq_ID %in% selected_clones))
      index_selected <- sort(unique(c(index_selected_sub,index_selected_clones)))
      index_colors_sub<- sapply(index_selected, function(z) if (z %in%index_selected_clones && z %in% index_selected_sub) {"#C07100"} else if(z %in%index_selected_clones){"#FFA500"}else {"#000000"})
      fig <- add_trace(
        fig, x = paste(single_group, sep = ""),
        y = tmp_summary$Value[index_nonselected], hoveron = "points", text = paste(
          tmp_summary$Seq_ID[index_nonselected], " Value=", tmp_summary$Value[index_nonselected],
          sep = ""
        ),
        key = tmp_summary$Seq_ID[index_nonselected], legendgroup = single_group,
        scalegroup = single_group, selectedpoints = if (any(tmp_summary$Seq_ID[index_nonselected] %in% c(selected_subclones, selected_clones))) {
          ((which(tmp_summary$Seq_ID[index_nonselected] %in% c(selected_subclones, selected_clones)) -
              1))
        } else {
          NULL
        }, name = single_group, box = list(visible = T, line = list(color = "#2D2926")),
        points = if(really_hide_dots){FALSE}else{"all"}, pointpos = 0, jitter = 0.25, scalemode = "count",
        meanline = list(visible = T),
        color = I(
          if (i%%2 != 0) {
            primary_color
          } else {
            "#707A6C"
          }
        ),
        marker = list(
          opacity = if (hide_dots) {
            0
          } else {
            0.75
          }, line = list(
            width = 2, color = if (i%%2 != 0) {
              primary_color
            } else {
              "#707A6C"
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
              primary_color
            } else {
              "#707A6C"
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

        fig <- add_trace(
          fig, x = paste(single_group, sep = ""),
          y = tmp_summary$Value[index_selected], hoveron = "points", text = paste(
            tmp_summary$Seq_ID[index_selected], " Value=", tmp_summary$Value[index_selected],
            sep = ""
          ),
          key = tmp_summary$Seq_ID[index_selected], legendgroup = single_group,
          scalegroup = single_group, selectedpoints = if (any(tmp_summary$Seq_ID[index_selected] %in% c(selected_subclones, selected_clones))) {
            ((which(tmp_summary$Seq_ID[index_selected] %in% c(selected_subclones, selected_clones)) -
                1))
          } else {
            NULL
          }, name = single_group, box = list(visible = T, line = list(color = "#2D2926")),
          points =if(really_hide_dots){FALSE}else{"all"}, pointpos = 0, jitter = 0.25, scalemode = "count",
          meanline = list(visible = T),
          color = I(
            if (i%%2 != 0) {
              primary_color
            } else {
              "#707A6C"
            }
          ),
          marker = list(
            opacity = if (hide_dots) {
              0
            } else {
              0.75
            }, line = list(
              width = 2, color = if (i%%2 != 0) {
                primary_color
              } else {
                "#707A6C"
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
                primary_color
              } else {
                "#707A6C"
              }
            )
          ),
          selected = list(
            marker = list(
              size = 10, opacity = if (hide_dots) {
                0.5
              } else {
                1
              }, color = index_colors_sub, zorder = 2
            )
          ),
          showlegend = F
        )
      }


    }

  }

  base_font <- 14
  font_size <- round(base_font * scale)


  fig <- layout(
    fig, title = list(font = list(size = font_size + 6),
                      text=paste(c(
                        "Violin plot of the <i>feature ",
                        name_values,
                        "</i> <br> between the groups ",
                        group_info
                      ),
                      collapse = ""
                      )),
    font = list(size = font_size),

    xaxis = list(
      tickfont = list(size = font_size)
    ),

    yaxis = list(
      tickfont = list(size = font_size)
    ),

    legend = list(
      font = list(size = font_size - 1),
      tracegroupgap = 0
    ),

    hoverlabel = list(
      font = list(size = font_size)
    ),
    margin = list(
      l = 80,
      r = 80,
      t = 90 + font_size,
      b = 40 ),
    # violingap = 0, violingroupgap = 0, violinmode = "overlay",
    dragmode = "lasso"
  ) %>%
    config(
      displaylogo = FALSE, modeBarButtonsToRemove = c("autoScale2d", "pan2d", "hoverCompareCartesian")
    )


  if ((show_reconstructed || length(unique(sequence_info_df$Seq_type[selected_rows])) >
       2) && compare_opposites)
  {

    lab_from_p <- function(p){
      if (is.na(p)) return("ns")
      if (p < 0.001) return("***")
      if (p < 0.01)  return("**")
      if (p < 0.05)  return("*")
      "ns"
    }

    cats <- (unique(meta_summary$Group))
    p_by_cat <- stats::setNames(numeric(length(cats)), cats)
    lab_by_cat <- stats::setNames(character(length(cats)), cats)
    ypos_by_cat <- stats::setNames(numeric(length(cats)), cats)

    for (cc in cats){
      dcc <- meta_summary[meta_summary$Group == cc, ]
      pa <- dcc$Value[dcc$Subgroup == "Repertoire"]
      pb <- dcc$Value[dcc$Subgroup != "Repertoire"]
      tt <- stats::t.test(pa, pb, paired=T)
      # print(cc)
      # print(tt$p.value)
      p_by_cat[cc] <- tt$p.value
      d <- (mean(pa) - mean(pb)) / stats::sd(c(pa,pb))
      lab_by_cat[cc] <- if(d<=0.2){""}else{lab_from_p(tt$p.value)}
      ypos_by_cat[cc] <- max(dcc$Value, na.rm = TRUE) + 0.4  # altura de la barra
    }
    n <- length(cats)


    centers <- stats::setNames(if (n > 1) (seq_len(n) - 0.6) / n else 0.5, cats)


    bar_frac <- 0.6
    half_width <- if (n > 1) (bar_frac / 2) * (1 / n) else 0.2  # para n=1, ancho razonable


    bump_by_cat <- stats::setNames(numeric(n), cats)
    for (cc in cats) {

      y_range_cc <- range(ypos_by_cat[[cc]], na.rm = TRUE)

      bump_by_cat[[cc]] <- 0.04 * diff(range(values, na.rm = TRUE))
    }


    shapes_list <- list()
    ann_list <- list()

    for (cc in cats){
      cx <- centers[[cc]]
      y0 <- ypos_by_cat[[cc]]
      bump <- bump_by_cat[[cc]]

      if(lab_by_cat[[cc]]!="") {
        shapes_list[[length(shapes_list) + 1]] <- list(
          type = "line",
          xref = "x domain", yref = "y",
          x0 = cx - half_width, x1 = cx + half_width,
          y0 = y0, y1 = y0,
          line = list(color = "black", width = 2)
        )


        ann_list[[length(ann_list) + 1]] <- list(
          xref = "x domain", yref = "y",
          x = cx, y = y0 + bump,
          text = lab_by_cat[[cc]],
          showarrow = FALSE,
          font = list(size = 19)
        )
      }


    }


    fig <- fig %>% layout(shapes = shapes_list, annotations = ann_list)
  }


  fig <- suppressWarnings(fig %>%
    layout(plot_bgcolor = "rgba(0, 0, 0, 0)",
           paper_bgcolor = "rgba(0, 0, 0, 0)",
           yaxis = list(  showgrid = TRUE,
                          gridcolor = "#495057",
           minor = list(dtick = 0.1, ticklen = 3)) ) %>%
    toWebGL())

  return(fig)
}

