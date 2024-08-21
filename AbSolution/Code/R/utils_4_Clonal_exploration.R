#' 4_Clonal_exploration
#'
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#'
#' @noRd
calculate_clone=function(seq_df,  clonotype, AA_or_NT="NT", region="CDR3", percentage=100, calculate_shared_clones) {

  if ( !(paste("Clone",AA_or_NT, region, clonotype,percentage,"simil", sep="_")  %in% colnames(seq_df)) ) {
    print("calculate_clones")
    if (clonotype=="VCDR3J") {

      # tmp=paste( seq_df[selected_rows,"V_and_J"], seq_df[selected_rows,paste(AA_or_NT, "CDR3",sep="_")],sep="_______")

      tmp=paste( seq_df[,get("V_and_J")], seq_df[,get(paste(AA_or_NT, "CDR3",sep="_"))],sep="_______")



    }

    if (clonotype =="Reconstructed_germline") {
      # tmp_selected_rows=sapply(selected_rows, function(z) intersect(which(seq_df$Sequence_type=="Reconstructed_germline"), which(seq_df[,"ID"] %in% seq_df[z,"ID"])) )
      #
      # tmp=paste( seq_df[tmp_selected_rows,"V_and_J"], seq_df[tmp_selected_rows,paste(AA_or_NT, "Whole",sep="_")],sep="_______")

      tmp_selected_rows=sapply(c(1:nrow(seq_df)), function(z) intersect(which(seq_df$Sequence_type=="Reconstructed_germline"), which(seq_df[,get("ID")] %in% seq_df[z,get("ID")])) )
      tmp=paste( seq_df[tmp_selected_rows,get("V_and_J")], seq_df[tmp_selected_rows,get(paste(AA_or_NT, "Whole",sep="_"))],sep="_______")

    }

    if(!calculate_shared_clones) {
      tmp=paste(tmp, seq_df$Patient_Sample, sep="__")
    }

    seq_df[, paste("Clone", AA_or_NT, region, clonotype,percentage, if(calculate_shared_clones){"shared"}, "simil", sep="_")] <- tmp
    print("calculate_clones_DONE")
    print(paste("Clone", AA_or_NT, region, clonotype,percentage, if(calculate_shared_clones){"shared"}, "simil", sep="_"))

  }



  return(seq_df)
}

#' 4_Clonal_exploration
#'
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#' @import plotly
#' @noRd
draw_violinplots=function(seq_df, group="Patient", selected_rows, clonotype, AA_or_NT="AA", region="CDR3", percentage=100, freq_filter=0,Selected_clones=NULL, dominance_threshold) {

  #devtools::install_github("karthik/wesanderson")
  # pal <- wes_palette(length(unique(seq_df[,..group]))*2, name = "Zissou1", type = "continuous")
  fig <- plot_ly(type = 'violin')

  i=0
  for (single_group in unique(seq_df[selected_rows,get(group)])) {

    # subset_seq_df=seq_df[intersect(selected_rows, which(seq_df[,get(group)]==single_group)),]
    # subset_group=unique(seq_df[intersect(selected_rows, which(seq_df[,get(group)]==single_group)), get("Group")])

    og_subset_seq_df=seq_df[intersect(selected_rows, which(seq_df[,get(group)]==single_group)),]
    subset_seq_df=seq_df[intersect(which(seq_df$Sequence_type %in% unique(seq_df$Sequence_type[selected_rows])), which(seq_df[,get(group)]==single_group)),]
    subset_group=unique(seq_df[intersect(which(seq_df$Sequence_type %in% unique(seq_df$Sequence_type[selected_rows])), which(seq_df[,get(group)]==single_group)), get("Group")])
    print(single_group)


    for (seq_type in c("Repertoire", "Reconstructed_germline")) {
      i=i+1
      print(seq_type)
      subsubset_seq_df=subset_seq_df[which(subset_seq_df$Sequence_type == seq_type),]
      og_subsubset_seq_df=og_subset_seq_df[which(og_subset_seq_df$Sequence_type == seq_type),]
      if(nrow(subsubset_seq_df)>0) {

        if(seq_type =="Reconstructed_germline") {
          print(subsubset_seq_df)
        }
        tmp_clones=table(subsubset_seq_df[, get(clonotype)])

        tmp_subclones=table(paste(subsubset_seq_df[, get(clonotype)], subsubset_seq_df[, get("NT_Whole")], sep="_"))
        tmp_frequencies=data.frame(Clone_ID=names(tmp_clones), Frequency=(100*as.numeric(tmp_clones))/(sum(tmp_clones)),
                                   Number_of_unique_subclones=sapply(names(tmp_clones), function(z) length(which(startsWith(names(tmp_subclones), paste(z, "_", sep=""))))  ),
                                   Number_of_sequences= as.numeric(tmp_clones))
        tmp_frequencies=tmp_frequencies[which(tmp_frequencies$Clone_ID %in% unique(og_subsubset_seq_df[, get(clonotype)])),]

        print("MMmmmMMm")
        print(Selected_clones$key)
        if(any(tmp_frequencies$Clone_ID %in% Selected_clones$key)){print(tmp_frequencies$Clone_ID[(which(tmp_frequencies$Clone_ID %in% Selected_clones$key))])}
        if(any(tmp_frequencies$Clone_ID %in% Selected_clones$key)){print(list((which(tmp_frequencies$Clone_ID %in% Selected_clones$key)-1)))}
        if(any(tmp_frequencies$Clone_ID %in% Selected_clones$key)){print("YESSSCLONES")}
        print("MMmmmMMm?")

        # FIX ####
        #writing tmp_frequencies to a file and checking later

        #saving also everything in a RData???
        fig <- add_trace(
          fig,
          x=paste( single_group, " (", paste(subset_group, collapse="&"), ")", sep=""),
          y = tmp_frequencies$Frequency[which(!(tmp_frequencies$Clone_ID %in% Selected_clones$key))],
          hoveron = "points",
          text= paste(tmp_frequencies$Clone_ID[which(!(tmp_frequencies$Clone_ID %in% Selected_clones$key))]," Number of subclones=", tmp_frequencies$Number_of_unique_subclones[which(!(tmp_frequencies$Clone_ID %in% Selected_clones$key))]," Number of sequences=", tmp_frequencies$Number_of_sequences[which(!(tmp_frequencies$Clone_ID %in% Selected_clones$key))], sep=""),
          key = tmp_frequencies$Clone_ID[which(!(tmp_frequencies$Clone_ID %in% Selected_clones$key))],
          legendgroup =seq_type,
          scalegroup  =  single_group,
          selectedpoints=if(any(tmp_frequencies$Clone_ID[which(!(tmp_frequencies$Clone_ID %in% Selected_clones$key))] %in% Selected_clones$key)){((which(tmp_frequencies$Clone_ID[which(!(tmp_frequencies$Clone_ID %in% Selected_clones$key))] %in% Selected_clones$key)-1))}else{NULL},
          name = seq_type,
          side = if(seq_type=="Reconstructed_germline") {'positive'} else {'negative'},
          box = list(
            visible = T,
            line=list(color="#1E2D24")
          ),
          points = 'all',
          pointpos = if(seq_type=="Reconstructed_germline") {0.5} else {-0.5},
          jitter = 0.25,
          scalemode = 'count',
          meanline = list(
            visible = T
          ),
          color =I( if (i%%2!=0) {"#2C497F"}else{"#808A9F"}),
          marker = list(
            opacity = 0.75,
            line = list(
              width = 2,
              color = if (i%%2!=0) {"#2C497F"}else{"#808A9F"}
            )
            # symbol = 'circle-open'
          ),
          unselected = list(
            marker = list(
              size = 2,
              opacity = 0.1,
              color = if (i%%2!=0) {"#2C497F"}else{"#808A9F"}
            )
          ),
          selected = list(
            marker = list(
              size = 10,
              opacity = 1,
              color = "#000000",
              zorder=2
            )
          ),
          showlegend = if(i>2){F}else{T}
        )

        if (length(which(tmp_frequencies$Clone_ID %in% Selected_clones$key))!=0) {
          fig <- add_trace(
            fig,
            x=paste( single_group, " (", paste(subset_group, collapse="&"), ")", sep=""),
            y = tmp_frequencies$Frequency[which(tmp_frequencies$Clone_ID %in% Selected_clones$key)],
            hoveron = "points",
            text= if (length(which(tmp_frequencies$Clone_ID %in% Selected_clones$key))==0) {NULL} else  {  paste(tmp_frequencies$Clone_ID[which(tmp_frequencies$Clone_ID %in% Selected_clones$key)]," Number of subclones=", tmp_frequencies$Number_of_unique_subclones[which(tmp_frequencies$Clone_ID %in% Selected_clones$key)]," Number of sequences=", tmp_frequencies$Number_of_sequences[which(tmp_frequencies$Clone_ID %in% Selected_clones$key)], sep="")},
            key = tmp_frequencies$Clone_ID[which(tmp_frequencies$Clone_ID %in% Selected_clones$key)],
            legendgroup =seq_type,
            scalegroup  =  single_group,
            selectedpoints=if(any(tmp_frequencies$Clone_ID[which(tmp_frequencies$Clone_ID %in% Selected_clones$key)] %in% Selected_clones$key)){((which(tmp_frequencies$Clone_ID[which(tmp_frequencies$Clone_ID %in% Selected_clones$key)] %in% Selected_clones$key)-1))}else{NULL},
            name = seq_type,
            side = if(seq_type=="Reconstructed_germline") {'positive'} else {'negative'},
            box = list(
              visible = T,
              line=list(color="#1E2D24")
            ),
            points = 'all',
            pointpos = if(seq_type=="Reconstructed_germline") {0.5} else {-0.5},
            jitter = 0.25,
            scalemode = 'count',
            meanline = list(
              visible = T
            ),
            color =I( if (i%%2!=0) {"#2C497F"}else{"#808A9F"}),
            marker = list(
              opacity = 0.75,
              line = list(
                width = 2,
                color = if (i%%2!=0) {"#2C497F"}else{"#808A9F"}
              )
              # symbol = 'circle-open'
            ),
            unselected = list(
              marker = list(
                size = 2,
                opacity = 0.1,
                color = if (i%%2!=0) {"#2C497F"}else{"#808A9F"}
              )
            ),
            selected = list(
              marker = list(
                size = 10,
                opacity = 1,
                color = "#000000",
                zorder=2
              )
            ),
            showlegend = F
          )
        }
      }


    }

  }

  fig <- layout(
    fig,
    title = paste(c("Relative frequency distribution of clones<br><i>using the criteria ", clonotype,"</i>"), collapse=""),
    # yaxis = list(
    #   zeroline = T,
    #   type = "log",
    #   exponentformat="power"
    # ),
    violingap = 0,
    violingroupgap = 0,
    violinmode = 'overlay',
    legend = list(
      tracegroupgap = 0
    ),
    dragmode = "lasso",
    shapes = list(hline((dominance_threshold)) )
  )%>%
    config(displaylogo = FALSE, modeBarButtonsToRemove = c("autoScale2d","pan2d", "hoverCompareCartesian"))

  fig=fig%>%
    layout(plot_bgcolor  = "rgba(0, 0, 0, 0)",
           paper_bgcolor = "rgba(0, 0, 0, 0)")  %>% toWebGL()

  return (fig)
}

#' 4_Clonal_exploration
#'
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#' @import upsetjs
#' @noRd
draw_upsetplot=function(seq_df, group="Patient", selected_rows, clonotype, AA_or_NT="AA", region="CDR3", percentage=100, freq_filter=0,Selected_clones=NULL){

  set_list=list()

  for (single_group in unique(seq_df[,get(group)])) {
    subset_seq_df=seq_df[intersect(selected_rows, which(seq_df[,get(group)]==single_group)),]
    subset_group=unique(seq_df[intersect(selected_rows, which(seq_df[,get(group)]==single_group)), get("Group")])
    print(single_group)
    subsubset_seq_df=subset_seq_df[which(subset_seq_df$Sequence_type == "Repertoire"),]
    tmp_clones=table(subsubset_seq_df[, get(clonotype)])
    set_list[[length(set_list)+1]]=unique(names(tmp_clones))

  }

  names(set_list)=unique(seq_df[,get(group)])



  upset_plot=upsetjs::upsetjs() %>% upsetjs::fromList(set_list)

  return(upset_plot)
}

#' 4_Clonal_exploration
#'
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#'
#' @noRd
hline <- function(y = 0, color = "darkblue") {
  list(
    type = "line",
    x0 = 0,
    x1 = 100,
    xref = "paper",
    y0 = y,
    y1 = y,
    line = list(color = color, width=0.75, dash="dot")
  )
}

#' 4_Clonal_exploration
#'
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#'
#' @noRd
vline <- function(x = 0, color = "darkblue") {
  list(
    type = "line",
    x0 = x,
    x1 = x,
    xref = "paper",
    y0 = 0,
    y1 = 100,
    line = list(color = color, width=0.75, dash="dot")
  )
}

#' 4_Clonal_exploration
#'
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#' @import plotly
#' @import utils
#' @noRd
draw_sharedclonesplot=function(seq_df, sets=NULL, group="Patient", selected_rows, clonotype, AA_or_NT="AA", region="CDR3", percentage=100, freq_filter=0,Selected_clones=NULL, dominance_threshold){

  list_figs=list()
  set_vector=c()

  for (group_A in unlist(sets)) {
    for (group_B in unlist(sets)) {
      if(group_A != group_B){
        tmp_comparison=paste(sort(c(group_A, group_B)), collapse="    vs    ")
        if(!(tmp_comparison %in% set_vector)){
          set_vector=c(set_vector, tmp_comparison)
          group1=sort(c(group_A, group_B))[1]
          group2=sort(c(group_A, group_B))[2]
          subset1_seq_df=seq_df[intersect(selected_rows, which(seq_df[,get(group)] == group1)),]
          subset1_seq_df=subset1_seq_df[which(subset1_seq_df$Sequence_type == "Repertoire"),]
          subset2_seq_df=seq_df[intersect(selected_rows, which(seq_df[,get(group)] == group2)),]
          subset2_seq_df=subset2_seq_df[which(subset2_seq_df$Sequence_type == "Repertoire"),]

          tmp_clones1=table(subset1_seq_df[, get(clonotype)])

          tmp_frequencies1=data.frame(Clone_ID=names(tmp_clones1),
                                      Frequency=(100*as.numeric(tmp_clones1))/((nrow(subset1_seq_df))),
                                      Number_of_subclones=as.numeric(tmp_clones1))

          tmp_clones2=table(subset2_seq_df[, get(clonotype)])

          tmp_frequencies2=data.frame(Clone_ID=names(tmp_clones2),
                                      Frequency=(100*as.numeric(tmp_clones2))/((nrow(subset2_seq_df))),
                                      Number_of_subclones=as.numeric(tmp_clones2))

          tmp_summary_df=merge(tmp_frequencies1,
                               tmp_frequencies2,
                               by='Clone_ID', all=T)
          tmp_summary_df[(is.na(tmp_summary_df))]=0



          tmp_summary_df$Freq_X=sapply(tmp_summary_df$Frequency.x, function(z) if(z==0){0.0005}else{z})
          tmp_summary_df$Freq_Y=sapply(tmp_summary_df$Frequency.y, function(z) if(z==0){0.0005}else{z})

          tmp_vector=tmp_summary_df$Freq_X
          for (tmp_value in unique(tmp_vector[which(duplicated(tmp_vector))])) {
            tmp_summary_df$Freq_X[which(tmp_vector == tmp_value)]=jitter(tmp_summary_df$Freq_X[which(tmp_vector == tmp_value)])
          }

          tmp_vector=tmp_summary_df$Freq_Y
          for (tmp_value in unique(tmp_vector[which(duplicated(tmp_vector))])) {
            tmp_summary_df$Freq_Y[which(tmp_vector == tmp_value)]=jitter(tmp_summary_df$Freq_Y[which(tmp_vector == tmp_value)])
          }
          tmp_summary_df$info=paste(paste(tmp_summary_df$Clone_ID, paste(tmp_summary_df$Frequency.x, tmp_summary_df$Frequency.y, sep=" , "), sep=" ("), ")",sep="")
          tmp_summary_df$shared=sapply(1:nrow(tmp_summary_df), function(z) if (tmp_summary_df$Frequency.x[z] !=0 && tmp_summary_df$Frequency.y[z] !=0){"Shared"} else {"Unique"})
          tmp_summary_df$shared=factor(tmp_summary_df$shared, levels=c("Shared", "Unique"))
          print(head(tmp_summary_df))
          tmp_fig=plot_ly(data = tmp_summary_df, x = ~Freq_X, y = ~Freq_Y, text=~info, hoverinfo="text", key=tmp_summary_df$Clone_ID, color =~shared, colors=c('#e63946', '#2b2d42'), marker = list(size = 10, opacity = 0.75))
          tmp_fig <- layout(tmp_fig, xaxis = list(title=paste("Log10 clonal frequency", group1, sep=" - "), type = "log", exponentformat="power", range=c(-5,2)), yaxis = list(title=paste("Log10 clonal frequency", group2, sep=" - "), type = "log", exponentformat="power",showexponent = "all", range=c(-5,2)),
                            shapes = list(hline((dominance_threshold)),vline((dominance_threshold)) ))
          list_figs[[length(list_figs)+1]]=tmp_fig
        }
      }
    }


  }

  sharedclonesplot= subplot(list_figs, nrows = ceiling(length(list_figs)/4),margin=0.05,titleY=T, titleX = TRUE)

  return(sharedclonesplot)
}
