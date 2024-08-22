#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @import Biostrings
#' @import shazam
#' @import alakazam
#' @import Biostrings
#' @import shinyjs
#' @import seqinr
#' @import dplyr
#' @import plotly
#' @import reactable
#' @import DT
#' @import shinyjqui
#' @importFrom benchmarkme get_ram
#' @import bigstatsr
#' @import data.table
#' @import Peptides
#' @import shinydashboard
#' @import shinyFiles
#' @import shinyjs
#' @import shinymanager
#' @import shinyWidgets
#' @import upsetjs
#' @import sortable
#' @import stats
#' @import sunburstR
#' @import tools
#' @import umap
#' @noRd
app_server <- function(input, output, session) {

  # ### Studying the input data
  lst1 <- reactive({
    validate(need(input$raw_files != "", "Select your folders"))

    if (is.null(input$raw_files)) {
      return(NULL)
    } else {

      # path_list <- as.list(input$raw_files$datapath)
      # tbl_list <- lapply(input$raw_files$datapath, read.table, header=TRUE, sep=";")

      df <- do.call(rbind, tbl_list)
      return(df)
    }
  })


  #######Page 0

  output$Conditional_Action_Filetype_description <- renderUI({
    HTML("Upload the table with your sample information in .txt or .tab format, delimited by tab. It must have at least 6 fields (Filename, Sample, Patient, Group and Subgroup).<br> <br>")

  })

  output$Conditional_Action_Filetype_upload <- renderUI({
    fileInput("raw_sample_file", "Upload", multiple = F, accept = c(".tab",".txt"))
  })


  sample_info_react<-reactiveValues(table=NULL, summary_status=F)
  observeEvent(input$raw_sample_file$datapath, {
    sample_info_react$table<-read.table(input$raw_sample_file$datapath, header = T,sep = "\t")
  })



  output$Conditional_Action_Move_to_1 <- renderUI({
    if (is.null(input$raw_sample_file) || all(is.numeric(unlist(input$base_folder))) ) {
      HTML("IMPORTANT! <i>Fill first the required fields. Only after that you will be able to proceed with the analysis.<i>")
    } else {
      sample_info<-sample_info_react$table
      actionButton("Move_to_1_airr", "Proceed to the next step: Pre-process your data", icon("angle-right"))
    }
  })
  output$Conditional_Action_Move_to_Analysis <- renderUI({
    if (is.null(input$raw_sample_file) || all(is.numeric(unlist(input$base_folder))) ) {
    } else {
      sample_info<-sample_info_react$table
      actionButton("Move_to_analysis", "Skip next steps: You already have AbSolution-feature-calculated files.",
                   icon("forward"),
                   style="color: #fff; background-color: #18bc9c")
    }
  })


  sample_table_react <- reactive({
    validate(need(input$raw_sample_file != "", "You need to upload your file with the sample information"))

    if (is.null(input$raw_sample_file)) {
      return(NULL)
    } else {

      sample_info<-sample_info_react$table
      if(all(c("Filename","Sample","Patient","Group","Subgroup" )%in%colnames(sample_info))) {
        return(sample_info)
      } else {
        error.table=data.frame(Error=paste("Modify and reupload again. The following column(s) are missing:",
                                           paste(c("Filename","Sample","Patient","Group","Subgroup" )[which(!(c("Filename","Sample","Patient","Group","Subgroup"  )%in%colnames(sample_info)))], collapse=", "),
                                           sep=" "))
        return(error.table)
      }


    }
  })
  output$raw_sample_file_out <- DT::renderDataTable({
    sample_table_react()
  },options = list(pageLength = -1,lengthMenu = list(c(15, -1), c("15", "All")),  info = FALSE,
                   initComplete = JS(
                     "function(settings, json) {",
                     "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
                     "}")
  )
  )
  volumes = getVolumes()
  shinyDirChoose(input, "base_folder", roots = volumes(), session = session)




  #######Page 1

  shinyDirChoose(input, "preinfolder_AIRR", roots = volumes(), session = session)

  dirname_AIRR <- reactive({parseDirPath(volumes, input$preinfolder_AIRR)})



  folder_information_table_AIRR <- reactive({

    if (all(unlist(input$preinfolder_AIRR) == 0)) {
      return(NULL)
    } else {
      sample_info<-sample_info_react$table
      sample_summary_status=data.frame(Filename=sample_info$Filename, Found_in_folder=rep("No", nrow(sample_info)))

      if(!all(is.numeric(unlist(input$preinfolder_AIRR)))) {
        sample_summary_status$Found_in_folder[which(paste(sample_summary_status$Filename, ".tsv",sep="")  %in% list.files(path = parseDirPath(volumes(),input$preinfolder_AIRR), full.names = F))] = "Yes"
      }
      return(sample_summary_status)
    }
  })


  conditional_button_preprocess_AIRR <- reactive({
    if (all(unlist(input$preinfolder_AIRR)== 0) ) {
      HTML("IMPORTANT! <i>Fill the required information and select the folder with the AIRR-Seq files. Only after that you will be able to proceed with the analysis.<i>")
    } else {
      if(!is.null(folder_information_table_AIRR()) && all(folder_information_table_AIRR()$Found_in_folder == "Yes")) {

        actionButton("Preprocess_AIRR", "Process AIRR datasets",
                     icon("forward"),
                     style="color: #fff; background-color: #18bc9c")

      } else {
        HTML('<p style="color:red">Check your folder! <br> <br>

             Something is missing/mislabelled, perhaps the filenames (without the .tsv extension) are wrong or the files are missing</p>')
      }
    }
  })

  output$Conditional_Action_Preprocess_AIRR <- renderUI({
    conditional_button_preprocess_AIRR()
  })

  output$folder_information <- DT::renderDataTable({
    folder_information_table()
  },options = list(pageLength = -1,lengthMenu = list(c(15, -1), c("15", "All")),  info = FALSE,
                   initComplete = JS(
                     "function(settings, json) {",
                     "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
                     "}")
  )
  )
  output$folder_information_AIRR <- DT::renderDataTable({
    folder_information_table_AIRR()
  },options = list(pageLength = -1,lengthMenu = list(c(15, -1), c("15", "All")),  info = FALSE,
                   initComplete = JS(
                     "function(settings, json) {",
                     "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
                     "}")
  )
  )


  output$folder_information_AIRR_steps_1_to_3 <- DT::renderDataTable({

    folder_information_table_AIRR_steps_1_to_3()

  },options = list(pageLength = -1,lengthMenu = list(c(15, -1), c("15", "All")),  info = FALSE,
                   initComplete = JS(
                     "function(settings, json) {",
                     "$(this.api().table().header()).css({'background-color': '#000', 'color': '#fff'});",
                     "}")
  )
  )



  folder_information_table_AIRR_steps_1_to_3 <- reactive({

    if (all(unlist(input$bigmem_folder) == 0)) {
      return(NULL)
    } else {
      sample_info<-sample_info_react$table
      print(sample_info)
      sample_summary_status=data.frame(Filename=paste(sample_info$Patient,sample_info$Group, sep="."), Rds_and_Bk_files_found_in_folder=rep("No", nrow(sample_info)))

      for(dummy_row in grep( "_dummy",sample_summary_status$Filename)){
        sample_summary_status$Filename[dummy_row]=sample_info$Filename[dummy_row]
      }
      if(!all(is.numeric(unlist(input$bigmem_folder)))) {
        sample_summary_status$Rds_and_Bk_files_found_in_folder[intersect(which(paste(sample_summary_status$Filename, ".bk",sep="")  %in% list.files(path = parseDirPath(volumes(),input$bigmem_folder), full.names = F)), which(paste(sample_summary_status$Filename, ".rds",sep="")  %in% list.files(path = parseDirPath(volumes(),input$bigmem_folder), full.names = F)))] = "Yes"

      }
      if (!all(is.numeric(unlist(input$bigmem_folder))) & all(sample_summary_status$Rds_and_Bk_files_found_in_folder == "Yes")) {
        updateSelectizeInput(session, "make_dummy_menu", choices = sample_summary_status$Filename[which(!grepl("_dummy",sample_summary_status$Filename))])

      }

      write.table(sample_info_react$table, file.path(parseDirPath(volumes, input$base_folder), "Sample_summary.txt"), append =F,row.names =F, col.names = T, sep="\t", quote = F)

      return(sample_summary_status)
    }
  })


  can_show_button_to_step<- reactiveValues(step_2=FALSE,step_2_AIRR=FALSE,step_3=FALSE, step_2_3=FALSE)
  folder_values<-reactiveValues(AIRR_parsed="", Featured="")


  # shinyjs::hide("pb_AIRR_vis", anim = TRUE)
  observeEvent(input$Preprocess_AIRR, {
    shinyjs::hide("Preprocess_AIRR", anim = TRUE)
    session$sendCustomMessage(type = 'testmessage',
                              message = 'Preprocessing')
    folder_values$AIRR_parsed=file.path(parseDirPath(volumes, input$base_folder), "1.Files_parsed")

    unlink(folder_values$AIRR_parsed, recursive = TRUE)
    dir.create(folder_values$AIRR_parsed)
    write.table(sample_info_react$table, file.path(parseDirPath(volumes, input$base_folder), "Sample_summary.txt"), append =F,row.names =F, col.names = T, sep="\t", quote = F)
    shinyjs::show("pb_AIRR_vis", anim = TRUE)

    inf_table=(sample_table_react())
    print("inf_table")
    print(inf_table)
    for (row_number_sample_table in c(1:nrow(inf_table))) {
      print(row_number_sample_table)

      parse_AIRRSeq_file(file=inf_table$Filename[row_number_sample_table],
                         group=inf_table$Group[row_number_sample_table],
                         patient=inf_table$Patient[row_number_sample_table],
                         subgroup=inf_table$Subgroup[row_number_sample_table],
                         sample=inf_table$Sample[row_number_sample_table],
                         input_path=paste(parseDirPath(volumes(), input$preinfolder_AIRR),"/",sep=""),
                         C_region_included=input$C_region_included_airr,  #####AARRRR
                         FWR1partial=input$FWR1partial_airr,
                         FWR4partial=input$FWR4partial_airr,
                         output_path=paste(folder_values$AIRR_parsed,"/",sep=""),
                         D_gene=input$Dgene_reconstruct_airr,
                         repertoire=input$TCRBCR_input_file)

      shinyWidgets::updateProgressBar(
        session = session,
        id = "pb_AIRR",
        value = 100*(row_number_sample_table/nrow(inf_table)), total = 100,
        title = paste("Process", trunc(100*(row_number_sample_table/nrow(inf_table))/10))
      )
      Sys.sleep(0.1)
    }
    print("Step done")
    shinyjs::hide("pb_AIRR_vis", anim = TRUE)
    can_show_button_to_step$step_2_AIRR=TRUE
  })


  output$Conditional_Action_Move_to_2 <- renderUI({
    if(can_show_button_to_step$step_2) {
      actionButton("Move_to_2", "Next step", icon("angle-right"))
    }

  })

  output$Conditional_Action_Move_to_2_AIRR <- renderUI({
    if(can_show_button_to_step$step_2_AIRR) {
      actionButton("Move_to_2", "Next step", icon("angle-right"))
    }

  })





  ##Page1&2

  shinyDirChoose(input, "bigmem_folder", roots = volumes(), session = session)

  observeEvent(input$bigmem_folder,{
    folder_values$Featured=file.path(parseDirPath(volumes(), input$bigmem_folder))
    sample_info<-sample_info_react$table
    all_ok=c()
    basenames=basename(list.files(path=folder_values$Featured, pattern=glob2rx("*.rds")))
    prefixes=file_path_sans_ext(basenames)
    prefixes=unique(prefixes)[which(unique(prefixes) != "Merged_bm")]
    for(n_grouped_by in c(1:length(strsplit(prefixes[1], split=".", fixed=T)))  ){
      tmp_all_ok=F
      for(sample_columna in c("Filename","Sample","Patient","Group","Subgroup")) {
        ### FUUUUUUUUUUUU #####
        if(all(sample_info[,sample_columna] %in% sapply(prefixes, function(z) strsplit(z, split=".", fixed=T)[[1]][n_grouped_by]))){
          tmp_all_ok=T
        }
      }

      all_ok=c(all_ok, tmp_all_ok)
    }

    #### TEMP ######
    all_ok=T
    if(all(all_ok)){
      can_show_button_to_step$step_2_3=TRUE
    } else {
      can_show_button_to_step$step_2_3=FALSE
    }
  })

  output$make_dummy <- renderUI({
    if (is.null(input$raw_sample_file) || all(is.numeric(unlist(input$bigmem_folder))) ) {
    } else {

      if(can_show_button_to_step$step_2_3 ) {
        materialSwitch("dummy_dataset", "Do you want to make a dummy dataset?", value = F)
      }

    }
  })



  observeEvent(input$final_dummy, {
    moment=paste0(Sys.Date(),"-", hour(Sys.time()), ':', minute(Sys.time()))
    dummy_name=paste("Dummy", input$make_dummy_menu, moment ,sep=".")
    print(dummy_name)
    print("T")
    og_FBM=big_attach(file.path(folder_values$Featured,paste0(input$make_dummy_menu, ".rds")))
    og_DT=fread(file.path(folder_values$Featured, paste0(input$make_dummy_menu, ".info")))
    og_header=fread(file.path(folder_values$Featured, paste0(input$make_dummy_menu, ".example")))
    # new_DT
    print("TA")
    list_new_DT=list()

    if("V_and_D_and_J" %in% colnames(og_DT)) {
      end=grep("V_and_D_and_J", colnames(og_DT))
    } else {
      end=grep("V_and_J", colnames(og_DT))
    }
    FWR1partial=if(sum(og_FBM[,grep("NT_FWR1_length", colnames(og_header))])==0){T}else{F}
    FWR4partial=if(sum(og_FBM[,grep("NT_FWR4_length", colnames(og_header))])==0){T}else{F}

    print("TAa")
    regions=c("FWR1","CDR1","FWR2","CDR2","FWR3","CDR3","FWR4")
    init=if(FWR1partial){2}else{1}
    finit=if(FWR4partial){6}else{7}
    og_DT=as.data.frame(og_DT)
    index_to_remove=c()
    for(iteration in c(1:input$dummy_ratio)){
      tmp_list=og_DT[,1:(end)]
      region_columns=paste("NT", regions, sep="_")
      number=1
      for(i_num in seq(from = 1, to = nrow(tmp_list), by =2)){
        tmp_list[i_num+1, region_columns]=tmp_list[i_num, region_columns]
        print(paste(number, nrow(tmp_list)/2, sep=" / "))
        number=number+1
        if(og_DT$ORF[i_num] != og_DT$ORF[i_num+1]) {
          index_to_remove=c(index_to_remove, c(i_num,i_num+1))
        } else {
          for (region in regions) {
            intermezzo=grep(region, regions)
            for (mut in c("Replacement_muts")) {
              tmp_num_mut=og_FBM[i_num+1, which(colnames(og_header) == paste("AA",region,mut,sep="_"))]
              if(tmp_num_mut > 0) {
                tmp_iter=tmp_num_mut

                used_positions=c()
                previous_sequences=c()
                while(tmp_iter >0 ){
                  sequence=tmp_list[i_num+1, paste("NT", region, sep="_")]
                  tmp_tmp_position=sample(c(1:nchar(sequence))[which(!(c(1:nchar(sequence)) %in% used_positions))], 1, replace=F)

                  sequence=strsplit(sequence, split="")[[1]]

                  sequence[[tmp_tmp_position]] = sample(c("A","C","G","T")[which(c("A","C","G","T") !=sequence[[tmp_tmp_position]])], 1)
                  sequence=paste(sequence, collapse="")



                  if(intermezzo==init) {
                    whole_seq=paste(c(sequence,
                                      tmp_list[i_num+1, paste("NT", regions[(intermezzo+1):(finit)], sep="_") ]),collapse="")
                  } else if (intermezzo==finit) {
                    whole_seq=paste(c(tmp_list[i_num+1, paste("NT", regions[init:(intermezzo-1)], sep="_")],
                                      sequence),collapse="")
                  } else {
                    whole_seq=paste(c(tmp_list[i_num+1, paste("NT", regions[init:(intermezzo-1)], sep="_")],
                                      sequence,
                                      tmp_list[i_num+1, paste("NT", regions[(intermezzo+1):(finit)], sep="_")]),collapse="")
                  }
                  prev_whole_seq=paste(c(tmp_list[i_num+1, paste("NT", regions, sep="_")]),collapse="")

                  if(og_DT$ORF[i_num]!=1) {
                    whole_seq=paste(strsplit(whole_seq, split="")[[1]][og_DT$ORF[i_num]:nchar(whole_seq)], collapse="")
                    prev_whole_seq=paste(strsplit(prev_whole_seq, split="")[[1]][og_DT$ORF[i_num]:nchar(prev_whole_seq)], collapse="")
                  }
                  whole_seq=Biostrings::DNAString(as.character(whole_seq))
                  prev_whole_seq=Biostrings::DNAString(as.character(prev_whole_seq))


                  suppressWarnings({
                    whole_seq_AA=as.character(Biostrings::translate(whole_seq))
                  }
                  )
                  suppressWarnings( {
                    prev_whole_seq_AA=as.character(Biostrings::translate(prev_whole_seq))
                  }
                  )

                  previous_sequences=unique(c(previous_sequences, prev_whole_seq_AA))
                  if (!(whole_seq_AA %in% previous_sequences)) {
                    seq_diff=T
                  } else {
                    seq_diff=F
                  }

                  if (!grepl(pattern = "*",whole_seq_AA, fixed = T,perl = F )) {
                    productive=T
                  } else {
                    productive=F
                  }

                  if(mut =="Replacement_muts" && seq_diff && productive) {
                    tmp_iter=tmp_iter-1
                    used_positions=c(used_positions, tmp_tmp_position)
                    tmp_list[i_num+1, paste("NT", region, sep="_")]=sequence
                  }


                }
              }

            }
          }
        }


      }

      tmp_list$ID=paste(tmp_list$ID, "Dummy", iteration, sep="_")
      tmp_list$Patient=paste("Dummy", tmp_list$Patient,  sep="_")
      tmp_list$Sample=paste("Dummy", tmp_list$Sample, sep="_")
      tmp_list$Group=paste( "Dummy", tmp_list$Group, sep="_")
      tmp_list$Subgroup=paste("Dummy", tmp_list$Subgroup, iteration,sep="_")

      tmp_list=tmp_list[c(1:nrow(tmp_list))[which(!c(1:nrow(tmp_list)) %in% index_to_remove)],]
      list_new_DT[[length(list_new_DT)+1]]=tmp_list
    }
    new_DT=rbindlist( list_new_DT )
    print("TAaaa")
    # new_FBM
    #



    Feature__dataset(path_base=folder_values$Featured,
                     DF_to_parse=new_DT,
                     name_DF_to_parse=dummy_name,
                     FWR1partial=FWR1partial,
                     FWR4partial=FWR4partial,
                     standard=F)
    print("TAaaaaa")
    sample_info<-sample_info_react$table
    print(sample_info)

    first=T
    for (coincidence in which(paste(sample_info$Patient,sample_info$Group, sep=".") == input$make_dummy_menu)) {
      tmp_dummy_sample_info=unname(t(as.data.frame(c(dummy_name,
                                                     paste(sample_info[coincidence,2:ncol(sample_info)],"dummy", sep="_")))))
      colnames(tmp_dummy_sample_info)=colnames(sample_info)
      if(first){

        dummy_sample_info=tmp_dummy_sample_info
      } else {
        dummy_sample_info=rbind(dummy_sample_info,tmp_dummy_sample_info )
      }
      first=F
    }


    sample_info=rbind(sample_info, dummy_sample_info)
    sample_info_react$table=sample_info
  })

  output$Conditional_Action_Move_to_Analysis_Real <- renderUI({
    if (is.null(input$raw_sample_file) || all(is.numeric(unlist(input$bigmem_folder))) ) {
    } else {

      if(can_show_button_to_step$step_2_3 ) {
        actionButton("Move_to_analysis_real", "Continue the analysis",
                     icon("forward"),
                     style="color: #fff; background-color: #18bc9c")
      } else {
        HTML('<p style="color:red">Check your folder! <br> <br>

             Some files are missing/mislabelled, check your folder with your sample info table</p>')
      }
    }
  })



  ##Page 2
  shinyjs::hide("pb_Feature_vis")

  Features_values<-reactiveValues(Current=0, Total=0)

  observeEvent(input$Feature_determination, {
    hide("Feature_determination")
    session$sendCustomMessage(type = 'testmessage',
                              message = 'Calculating features, be patient!')
    folder_values$Featured=file.path(parseDirPath(volumes(), input$base_folder), "2.Feature_determination")

    unlink(folder_values$Featured, recursive = TRUE)
    dir.create(folder_values$Featured)
    shinyjs::show("pb_Feature_vis")


    Features_values$Total= Feature_1(parseDirPath(volumes, input$base_folder),grouping_by=c("Patient", "Group"))
    print(Features_values$Total)
    print(file.path(paste(parseDirPath(volumes, input$base_folder),"/2.Feature_determination",sep=""), "IMGT_parsed_index_extended.txt"))
    List_dfs=split(as.data.table(fread(file.path(paste(parseDirPath(volumes, input$base_folder),"/2.Feature_determination",sep=""), "IMGT_parsed_index_extended.txt"), header=T, sep="\t", quote =FALSE)), by=c("Patient", "Group"))

    print("Feature calculation")
    for (i in c(1:length(List_dfs))) {
      print(paste(i, length(List_dfs), sep=" / "))
      Features_values$Current=Features_values$Current +1
      Feature__dataset(path_base=parseDirPath(volumes, input$base_folder),
                       DF_to_parse=List_dfs[[i]],
                       name_DF_to_parse=names(List_dfs)[i],
                       FWR1partial=input$FWR1partial_airr,
                       FWR4partial=input$FWR4partial_airr)
    }
    rm(List_dfs)
    hide("pb_Feature_vis")
    can_show_button_to_step$step_3=TRUE
  })

  output$Conditional_Action_Move_to_3 <- renderUI({
    if(can_show_button_to_step$step_3) {
      actionButton("Move_to_analysis_real", "Proceed to next step: Visualization and exploration (be patient)", icon("angle-right"))
    }

  })

  observeEvent(Features_values$Current, {
    if (Features_values$Current == 1) {
      shinyWidgets::updateProgressBar(
        session = session,
        id = "pb_Feature",
        value = 10, total = 100,
        title = paste("Now it will take some time, be patient!")
      )
    } else {
      shinyWidgets::updateProgressBar(
        session = session,
        id = "pb_Feature",
        value = 10 +90*((Features_values$Current -1)/Features_values$Total), total = 100,
        title = paste("Now it will take some time, be patient!")
      )
    }

  })

  # 3.Dataset exploration and variable selection ####
  Big_mem_values<-reactiveValues(Header=NULL, Short_DF=NULL, Big_DF=NULL, Run=0, Patient_Sample=NULL, VJs=NULL)
  Big_mem_color_values<-reactiveValues(V=NULL, D=NULL, J=NULL, VJ=NULL, VDJ=NULL)
  Selection_values<-reactiveValues(rows=NULL,columns=NULL, Scores=NULL, Variance_explained=NULL, UMAP=NULL, Parameters=NULL, Cell=F)
  Exploration_values<-reactiveValues(rows=NULL,columns=NULL, Scores=NULL, Variance_explained=NULL,UMAP=NULL, Parameters=NULL)

  observeEvent(input$Move_to_analysis_real, {

    info=merge_FBMs(folder_values$Featured)
    print("Merged")
    Big_mem_values$Header=info[[1]]
    Big_mem_values$Short_DF=info[[2]]
    Big_mem_values$Short_DF$Patient_Sample=paste(Big_mem_values$Short_DF$Patient,Big_mem_values$Short_DF$Sample, sep="__")
    Big_mem_values$Short_DF$Text_ID=paste(Big_mem_values$Short_DF$ID,Big_mem_values$Short_DF$Sequence_type, sep="_&_")
    Big_mem_values$Big_DF=info[[3]]
    updateSliderInput(inputId="Rmut_filter", min=min(Big_mem_values$Big_DF[,which(Big_mem_values$Header == "AA_Whole_Replacement_muts")]), max= max(Big_mem_values$Big_DF[,which(Big_mem_values$Header == "AA_Whole_Replacement_muts")]))

    Big_mem_values$Patient_Sample=unique(Big_mem_values$Short_DF$Patient_Sample)
    Big_mem_values$VJs=sort(unique(Big_mem_values$Short_DF$V_and_J))
    counts_VJs=table(Big_mem_values$VJs)/2
    names(Big_mem_values$VJs)= sapply(Big_mem_values$VJs, function(z) paste(z,  counts_VJs[which(names(counts_VJs)==z)], sep=" - Counts:") )
    rm(info)


    if(any(grepl(",", Big_mem_values$Short_DF$Best_V, fixed=T))) {
      Big_mem_values$Short_DF =  Big_mem_values$Short_DF%>%
        mutate(Best_V=gsub(",.*","",Best_V)) %>%
        mutate(Best_J=gsub(",.*","",Best_J)) %>%
        mutate(Best_D=gsub(",.*","",Best_D)) %>%
        mutate(V_and_J=paste0(Best_V,"_",Best_J)) %>%
        mutate(V_and_D_and_J=paste0(Best_V,"_",Best_D,"_",Best_J))
    }
    list_values=list()
    list_values[["V"]]=list()
    list_values[["D"]]=list()
    list_values[["J"]]=list()
    # test=test_orig
    # test_orig=test

    print(unique(Big_mem_values$Short_DF$Best_V))
    # test=test[startsWith(test$V_and_J, "IGH"),]
    for (gene in unique(Big_mem_values$Short_DF$Best_V)) {
      print(gene)
      list_values[["V"]][[strsplit(strsplit(strsplit(gene, split="*", fixed=T)[[1]][1], split="/")[[1]][1], split="-")[[1]][1]]]= sort(unique(c(list_values[["V"]][[strsplit(strsplit(strsplit(gene, split="*", fixed=T)[[1]][1], split="/")[[1]][1], split="-")[[1]][1]]], gene)))
    }

    print(unique(Big_mem_values$Short_DF$Best_J))
    for (gene in sort(unique(Big_mem_values$Short_DF$Best_J))) {
      print(gene)
      list_values[["J"]][[strsplit(strsplit(strsplit(gene, split="*", fixed=T)[[1]][1], split="/")[[1]][1], split="-")[[1]][1]]]= sort(unique(c(list_values[["J"]][[strsplit(strsplit(strsplit(gene, split="*", fixed=T)[[1]][1], split="/")[[1]][1], split="-")[[1]][1]]], gene)))
    }

    # for (gene in sort(unique(test$Best_V))) {
    #   list_values[["V"]][[strsplit(strsplit(gene, split="/")[[1]][1], split="-")[[1]][1]]]= sort(unique(c(list_values[["V"]][[strsplit(strsplit(gene, split="/")[[1]][1], split="-")[[1]][1]]], gene)))
    # }
    # for (gene in sort(unique(test$Best_D))) {
    #   list_values[["D"]][[strsplit(strsplit(gene, split="/")[[1]][1], split="-")[[1]][1]]]= sort(unique(c(list_values[["D"]][[strsplit(strsplit(gene, split="/")[[1]][1], split="-")[[1]][1]]], gene)))
    # }
    #
    # for (gene in sort(unique(test$Best_J))) {
    #   list_values[["J"]][[strsplit(strsplit(gene, split="/")[[1]][1], split="-")[[1]][1]]]= sort(unique(c(list_values[["J"]][[strsplit(strsplit(gene, split="/")[[1]][1], split="-")[[1]][1]]], gene)))
    # }

    ##########Test
    suppressWarnings({
      Big_mem_color_values$V=Ab_palette(list_values[["V"]], vect_genes_comb=NA, type_values=c("V"))

      Big_mem_color_values$J=Ab_palette(list_values[["J"]], vect_genes_comb=NA, type_values=c("J"))

      Big_mem_color_values$VJ=Ab_palette(list_values[c("V","J")], vect_genes_comb=sort(unique(Big_mem_values$Short_DF$V_and_J)), type_values=c("VJ"))

      if("Best_D" %in% colnames(Big_mem_values$Short_DF)) {
        for (gene in sort(unique(Big_mem_values$Short_DF$Best_D))) {
          list_values[["D"]][[strsplit(strsplit(strsplit(gene, split="*", fixed=T)[[1]][1], split="/")[[1]][1], split="-")[[1]][1]]]= sort(unique(c(list_values[["D"]][[strsplit(strsplit(strsplit(gene, split="*", fixed=T)[[1]][1], split="/")[[1]][1], split="-")[[1]][1]]], gene)))
        }
        Big_mem_color_values$D=Ab_palette(list_values[["D"]], vect_genes_comb=NA, type_values=c("D"))
        Big_mem_color_values$VDJ=Ab_palette(list_values[c("V","D","J")], vect_genes_comb=sort(unique(Big_mem_values$Short_DF$V_and_D_and_J)), type_values=c("VDJ"))
      }
    })





    # test$Chain[1]
    # seecol(usecol(pal = Ab_palette(list_values, vect_genes_comb=sort(unique(test$V_and_J)), type_values=c("VJ")), n = "all"))
    # seecol(usecol(pal = Ab_palette(list_values[["V"]], vect_genes_comb=NA, type_values=c("V")), n = "all"))
    # seecol(usecol(pal = Ab_palette(list_values[["D"]], vect_genes_comb=NA, type_values=c("D")), n = "all"))
    # seecol(usecol(pal = Ab_palette(list_values[["J"]], vect_genes_comb=NA, type_values=c("J")), n = "all"))
    # seecol(usecol(pal = Ab_palette(list_values[c("V","D","J")], vect_genes_comb=sort(unique(test$V_and_D_and_J)), type_values=c("VDJ")), n = "all"))


    tmp_rows_cols=filter_merged(Big_mem_values$Big_DF,
                                Big_mem_values$Short_DF,
                                Big_mem_values$Header,
                                "Reconstructed germline" %in%  input$use_what,
                                "Repertoire" %in%  input$use_what,
                                "Productive" %in% input$use_productive_or_not,
                                "Non-productive" %in% input$use_productive_or_not,
                                input$my_regions,
                                input$my_var_elements,
                                input$my_vars,
                                input$my_vartypes,
                                input$use_sharedVDJ,
                                input$VJ_included,
                                input$groups_selected,
                                input$group_A,
                                input$group_B,
                                input$group_C,
                                input$use_univlog,
                                input$samples_selected,
                                input$exclude_variables,
                                input$pval_type,
                                input$pval_cutoff,
                                input$estimate_cutoff,
                                input$number_selected_vars,
                                input$VJ_deselected,
                                input$VDJ_normalized_per_size,
                                input$Rmut_filter[1],
                                input$Rmut_filter[2],
                                input$work_as_categories,
                                input$VDJ_maximize_clones,
                                input$VDJ_normalized_per_sample,
                                input$my_clone_def)

    Exploration_values$rows=tmp_rows_cols$ROWS
    Exploration_values$columns=tmp_rows_cols$COLUMNS
    Selection_values$rows=tmp_rows_cols$ROWS
    Selection_values$columns=tmp_rows_cols$COLUMNS


    rm(tmp_rows_cols)
    ##PCAS base



  })
  toListenSelection <- reactive({
    list(Selection_values$rows,Selection_values$columns,Selection_values$Parameters)
  })

  observeEvent(toListenSelection(), {
    if(!is.null(Selection_values$rows) && length(Selection_values$rows)>0 && !is.null(Selection_values$columns) && length(Selection_values$columns)>0){
      print("PCA")

      if (is.null(input$plot_color)) {
        selection="Group"
      } else {
        selection=input$plot_color
      }

      print(paste("PCAAAA:", length(Selection_values$rows), sep=""))
      tmp_PCA=big_PCA(FBM=Big_mem_values$Big_DF,
                           rows=Selection_values$rows,
                           columns=Selection_values$columns)
      tmp_PCA[[1]]=as.data.frame(tmp_PCA[[1]])
      colnames(tmp_PCA[[1]])=paste("Dim_", c(1:5), sep="")
      print("PCA 2")
      tmp_PCA[[1]]$Color=as.factor(unlist(Big_mem_values$Short_DF[Selection_values$rows,..selection]))
      print("PCA 3")
      tmp_PCA[[1]]$Text_ID=Big_mem_values$Short_DF[Selection_values$rows,get("Text_ID")]
      tmp_PCA[[1]]$Seq_type=factor(unlist(Big_mem_values$Short_DF[Selection_values$rows,get("Sequence_type")]), levels=c("Reconstructed_germline", "Repertoire"))
      Selection_values$Scores=tmp_PCA[[1]]
      print(paste("PCAAAA Scores:", nrow(Selection_values$Scores), sep=""))
      Selection_values$Variance_explained=tmp_PCA[[2]]

      Big_mem_values$Run=Big_mem_values$Run+1

      if(input$use_UMAP == T && as.numeric(unname(strsplit(as.character(print(get_ram(), unit_system = "metric")), split=" ")[[1]][1])) >   round(length(Selection_values$rows)*length(Selection_values$columns)*8/2^{20}/1024, 2)) {
        print("UMAP")
        tmp_umap=as.data.frame(umap(Big_mem_values$Big_DF[Selection_values$rows, Selection_values$columns])$layout)
        print(dim(tmp_umap))
        colnames(tmp_umap)=paste("Dim_", c(1:2), sep="")
        tmp_umap$Color=as.factor(unlist(Big_mem_values$Short_DF[Selection_values$rows,..selection]))
        tmp_umap$Text_ID=paste(unlist(Big_mem_values$Short_DF[Selection_values$rows,get("ID")]),unlist(Big_mem_values$Short_DF[Selection_values$rows,get("Sequence_type")]), sep="_&_")
        tmp_umap$Seq_type=factor(unlist(Big_mem_values$Short_DF[Selection_values$rows,get("Sequence_type")]), levels=c("Reconstructed_germline", "Repertoire"))
        dim(tmp_umap)
        Selection_values$UMAP=tmp_umap
        updateSelectInput(session, "Selection_plot_type",
                          choices = c("PCA" = "PCA", "UMAP" = "UMAP"), selected = "PCA")
        print("DONEUMAP")
      } else {
        updateSelectInput(session, "Selection_plot_type",
                          choices = c("PCA" = "PCA"))
      }
    }

  })


  toListenSelection_plot_type<- reactive({
    list(input$Selection_plot_type)
  })


  observeEvent(toListenSelection_plot_type(), {
    if(!is.null(input$Selection_plot_type)) {
      if (input$Selection_plot_type =="PCA") {
        updateSelectizeInput(session, "Selection_plot_type_dim", choices = c(1:5), selected=c(1, 2))
      } else if (input$Selection_plot_type =="UMAP") {
        updateSelectizeInput(session,"Selection_plot_type_dim", choices = c(1:2), selected=c(1, 2))
      }
    }

  })

  toListenExploration_plot_type<- reactive({
    list(input$Exploration_plot_type)
  })

  observeEvent(toListenExploration_plot_type(), {
    if(!is.null(input$Exploration_plot_type)) {
      if (input$Exploration_plot_type =="PCA") {
        updateSelectizeInput(session,"Exploration_plot_type_dim", choices = c(1:5), selected=c(1, 2))
      } else if (input$Exploration_plot_type =="UMAP") {
        updateSelectizeInput(session,"Exploration_plot_type_dim", choices = c(1:2), selected=c(1, 2))
      }
    }
  })

  toListenExploration <- reactive({
    list(Exploration_values$rows,Exploration_values$columns,Exploration_values$Parameters)
  })

  observeEvent(toListenExploration(), {
    if(!is.null(Exploration_values$rows) && length(Exploration_values$rows)>0 && !is.null(Exploration_values$columns) && length(Exploration_values$columns)>0){

      # if (all(Exploration_values$rows==Selection_values$rows) && all(Exploration_values$columns==Selection_values$columns) ) {
      #   Exploration_values$Scores=tmp_PCA[[1]]
      #   Exploration_values$Variance_explained=tmp_PCA[[2]]
      #   rm(tmp_PCA)
      # } else {
      print("PCA ex")
      if (is.null(input$plot_color_expl)) {
        selection="Group"
      } else {
        selection=input$plot_color_expl
      }

      tmp_PCAex=big_PCA(FBM=Big_mem_values$Big_DF,
                             rows=Exploration_values$rows,
                             columns=Exploration_values$columns)
      tmp_PCAex[[1]]=as.data.frame(tmp_PCAex[[1]])
      colnames(tmp_PCAex[[1]])=paste("Dim_", c(1:5), sep="")
      tmp_PCAex[[1]]$Color=as.factor(unlist(Big_mem_values$Short_DF[Exploration_values$rows,..selection]))
      tmp_PCAex[[1]]$Text_ID=Big_mem_values$Short_DF[Exploration_values$rows,get("Text_ID")]
      tmp_PCAex[[1]]$Seq_type=factor(unlist(Big_mem_values$Short_DF[Exploration_values$rows,get("Sequence_type")]), levels=c("Reconstructed_germline", "Repertoire"))
      Exploration_values$Scores=tmp_PCAex[[1]]
      print("AQuit3")
      print(dim(tmp_PCAex[[1]]))
      print(dim(Exploration_values$Scores))
      print(colnames(Exploration_values$Scores))
      Exploration_values$Variance_explained=tmp_PCAex[[2]]


      if(input$use_UMAP && as.numeric(unname(strsplit(as.character(print(get_ram(), unit_system = "metric")), split=" ")[[1]][1])) >   round(length(Exploration_values$rows)*length(Exploration_values$columns)*8/2^{20}/1024, 2)) {
        print("UMAPEx")
        tmp_umap=as.data.frame(umap(Big_mem_values$Big_DF[Exploration_values$rows, Exploration_values$columns])$layout)
        print(dim(tmp_umap))
        colnames(tmp_umap)=paste("Dim_", c(1:2), sep="")
        tmp_umap$Color=as.factor(unlist(Big_mem_values$Short_DF[Exploration_values$rows,..selection]))
        tmp_umap$Text_ID=paste(unlist(Big_mem_values$Short_DF[Exploration_values$rows,get("ID")]),unlist(Big_mem_values$Short_DF[Exploration_values$rows,get("Sequence_type")]), sep="_&_")
        tmp_umap$Seq_type=factor(unlist(Big_mem_values$Short_DF[Exploration_values$rows,get("Sequence_type")]), levels=c("Reconstructed_germline", "Repertoire"))

        Exploration_values$UMAP=tmp_umap
        updateSelectInput(session, "Exploration_plot_type",
                          choices = c("PCA" = "PCA", "UMAP" = "UMAP"), selected = "PCA")
        print("DONEUMAPEx")
        rm(tmp_umap)
      } else {
        print("ELSEUMAPEx")
        updateSelectInput(session, "Exploration_plot_type",
                          choices = c("PCA" = "PCA"), selected = "PCA")
      }



      rm(tmp_PCAex)
      Big_mem_values$Run=Big_mem_values$Run+1


      # }

    }

  })



  ID_selected_values<-reactiveValues(subclones=NULL,clones=NULL, intersection_samples=NULL)

  output$Selection_plot <- renderPlotly({
    if(!is.null(Selection_values$rows) && length(Selection_values$rows)>0 && !is.null(Selection_values$columns) && length(Selection_values$columns)>0) {
      print("Selection plot")
      # Sel_colors=unlist(branded_colors[c(1:length(unique(Selection_values$Scores$Color)))])
      # names(Sel_colors)=sort(unique(Selection_values$Scores$Color))
      # Sel_colors_BORDER=c(Sel_colors, "#050609", "#ff4d6d", "#FFAD05")
      # names(Sel_colors_BORDER)=c(names(Sel_colors),"User selected", "Counterpart", "Clones")
      # Sel_colors_BORDER_width=c(rep(1, length(Sel_colors)),4,3,4)
      # names(Sel_colors_BORDER_width)=c(names(Sel_colors),"User selected", "Counterpart", "Clones")
      # Sel_colors_BORDER_symbol=c("circle", "star-diamond")
      # names(Sel_colors_BORDER_symbol)=c("Repertoire", "Reconstructed_germline")
      #
      #
      #
      # fig_sel=plot_ly(data = Selection_values$Scores ,x =  ~Dim_1, y = ~Dim_2, color = ~Color, colors=Sel_colors,  opacity = 0.7,
      #                 text=~Text_ID, key = ~Text_ID, type = 'scatter', mode = 'markers', hovertemplate = paste('<b>%{text}</b>'),
      #                 marker=list(
      #                   size = 7,
      #                   symbol = Sel_colors_BORDER_symbol[match(Selection_values$Scores$Seq_type, names(Sel_colors_BORDER_symbol))],
      #                   line = list(
      #                     # color = ~Selected,
      #                     color=unname(Sel_colors_BORDER[match(Selection_values$Scores$Selected, names(Sel_colors_BORDER))]),
      #                     # colors=Ex_colors_BORDER,
      #                     width =unname( Sel_colors_BORDER_width[match(Selection_values$Scores$Selected, names(Sel_colors_BORDER_width))])
      #                   )
      #                 )
      # )%>%
      #   layout(title = 'Selection plot', plot_bgcolor = "#e5ecf6",legend = list(orientation = 'h',y=-0.3),showlegend=T, xaxis = list(title = paste('Dim 1 (', Selection_values$Variance_explained[1], "%)", sep="") ),
      #          yaxis = list(title = paste('Dim 2 (', Selection_values$Variance_explained[2], "%)", sep=""),dragmode = "lasso")   )%>%
      #   config(displaylogo = FALSE, modeBarButtonsToRemove = c("zoomIn2d", "zoomOut2d", "autoScale2d","pan2d", "hoverCompareCartesian")) %>%
      #   style(legendgroup = NULL)
      #
      # ID_selected_values$subclones=event_data("plotly_selected")
      #
      # fig_sel

      if(input$plot_color == "Best_V"){
        print("VVV")
        Sel_colors=Big_mem_color_values$V
      } else if (input$plot_color == "Best_J") {
        Sel_colors=Big_mem_color_values$J
      } else if(input$plot_color == "Best_D") {

        Sel_colors=Big_mem_color_values$D

      } else if(input$plot_color == "V_and_D_and_J") {
        print("V_and_D_and_J")
        Sel_colors=Big_mem_color_values$VDJ

      } else if(input$plot_color == "V_and_J"){
        print("VV_and_JVV")
        Sel_colors=Big_mem_color_values$VJ
      } else {
        print("ADAD")
        print(input$plot_color)
        # Sel_colors=unlist(branded_colors[c(1:length(unique(Selection_values$Scores$Color)))])
        Sel_colors=Ab_palette(list_values=unique(Selection_values$Scores$Color),
                                     vect_genes_comb=NA,
                                     type_values="cualitative",
                                     colorblind=F)
        names(Sel_colors)=sort(unique(Selection_values$Scores$Color))
      }

      print("Sel selections")
      print(Sel_colors[1:10])
      print(length(Sel_colors))
      print(length(names(Sel_colors)))
      print(names(Sel_colors)[1:10])
      print(length(sort(unique(Selection_values$Scores$Color))))

      print("_____________")
      print(Selection_values$Scores$Selected[1:10])
      print(Selection_values$Scores$ID[1:10])
      print("_____________")
      Sel_colors_BORDER=c(Sel_colors, "#050609", "#ff4d6d", "#FFAD05")
      names(Sel_colors_BORDER)=c(names(Sel_colors),"User selected", "Counterpart", "Clones")
      Sel_colors_BORDER_width=c(rep(1, length(Sel_colors)),4,4,4)
      names(Sel_colors_BORDER_width)=c(names(Sel_colors),"User selected", "Counterpart", "Clones")
      Sel_colors_BORDER_symbol=c("circle", "star-diamond")
      names(Sel_colors_BORDER_symbol)=c("Repertoire", "Reconstructed_germline")

      print(input$Selection_plot_type)
      if(input$Selection_plot_type == "PCA") {

        test=Selection_values$Scores

      } else if(input$Selection_plot_type == "UMAP") {
        test=Selection_values$UMAP
        print(dim(test))
        print(colnames(test))
        print(head(test))
      }

      test$Colors=Sel_colors[match(test$Color, names(Sel_colors))]
      test$Color_border=Sel_colors_BORDER[match(test$Selected, names(Sel_colors_BORDER))]
      test$Width_border=Sel_colors_BORDER_width[match(test$Selected, names(Sel_colors_BORDER_width))]
      test$Symbol=Sel_colors_BORDER_symbol[match(test$Seq_type, names(Sel_colors_BORDER_symbol))]

      if(length(input$Selection_plot_type_dim)<2){
        dim1=1
        dim2=2
      }else{
        dim1=as.numeric(input$Selection_plot_type_dim[1])
        dim2=as.numeric(input$Selection_plot_type_dim[2])

      }
      tmpDim_1=test[,which(colnames(test)==paste0("Dim_", dim1))]
      tmpDim_2=test[,which(colnames(test)==paste0("Dim_", dim2))]
      test$Dim_1=tmpDim_1
      test$Dim_2=tmpDim_2

      fig <- plot_ly(data=test, type = 'scatter', mode = 'markers',
                     colors = Sel_colors)  %>%
        config(
          toImageButtonOptions = list(
            format = "png",
            filename = "PCA_Selection_plot",
            width = 1400,
            height = 1000
          )
        )


      ##NEW
      test=test[(order(test$Color)),]

      tmp_test_ns=test[intersect(intersect(which(test$Selected != "Clones"),
                                           which(test$Selected != "User selected")),
                                 which(test$Selected != "Counterpart")),]

      fig <- fig %>%
        add_trace(data=tmp_test_ns,
                  x = ~Dim_1,
                  y = ~Dim_2,
                  opacity = 0.7,
                  text =  ~ Text_ID,
                  key = ~ Text_ID,
                  color = ~Color,
                  colors = Sel_colors,
                  # size = ~Width_border,
                  type = 'scatter',
                  mode = 'markers',
                  fill = ~'',
                  hovertemplate = paste('<b>%{text}</b>'),
                  marker = list(sizemode = 'diameter'),
                  name= tmp_test_ns$Color,
                  legendgroup= tmp_test_ns$Color
        )


      tmp_test=test[which(test$Selected == "Counterpart"),]

      if(nrow(tmp_test)>0) {
        fig <- fig %>%
          add_trace(data=tmp_test,
                    x = ~Dim_1,
                    y = ~Dim_2,
                    opacity = 0.7,
                    text =  ~ Text_ID,
                    key = ~ Text_ID,
                    color = ~Color,
                    colors = Sel_colors,
                    # size = ~Width_border,
                    symbol= ~Symbol,
                    type = 'scatter',
                    mode = 'markers',
                    fill = ~'',
                    hovertemplate = paste('<b>%{text}</b>'),
                    marker =  list(
                      line=list(color="#ff4d6d",
                                width=tmp_test$Width_border[1])),
                    name= tmp_test$Color, showlegend=FALSE,
                    legendgroup=tmp_test$Color
          )
      }

      tmp_test=test[which(test$Selected == "Clones"),]

      if(nrow(tmp_test)>0) {
        not_in_leg=which(unique(tmp_test$Color) %!in% unique(tmp_test_ns$Color))
        if(length(not_in_leg)>0) {
          tmp_tmp_test=tmp_test[which(tmp_test$Color %in% not_in_leg),]

          fig <- fig %>%
            add_trace(data=tmp_tmp_test,
                      x = ~x,
                      y = ~y,
                      opacity = 0.9,
                      text =  ~ Text_ID,
                      key = ~ Text_ID,
                      color = ~Color,
                      symbol= ~Symbol,
                      # size = ~Width_border,
                      type = 'scatter',
                      mode = 'markers',
                      fill = ~'',
                      hovertemplate = paste('<b>%{text}</b>'),
                      marker = list(
                        line=list(color="black",
                                  width=tmp_tmp_test$Width_border[1])),
                      name= tmp_tmp_test$Color, showlegend=TRUE,
                      legendgroup=tmp_tmp_test$Color
            )
        }
        tmp_tmp_test=tmp_test[which(tmp_test$Color %!in% not_in_leg),]
        if(nrow(tmp_tmp_test)>0) {
          fig <- fig %>%
            add_trace(data=tmp_tmp_test,
                      x = ~x,
                      y = ~y,
                      opacity = 0.9,
                      text =  ~ Text_ID,
                      key = ~ Text_ID,
                      color = ~Color,
                      symbol= ~Symbol,
                      # size = ~Width_border,
                      type = 'scatter',
                      mode = 'markers',
                      fill = ~'',
                      hovertemplate = paste('<b>%{text}</b>'),
                      marker = list(
                        line=list(color="black",
                                  width=tmp_tmp_test$Width_border[1])),
                      name= tmp_tmp_test$Color, showlegend=FALSE,
                      legendgroup=tmp_tmp_test$Color
            )
        }

      }

      tmp_test=test[which(test$Selected == "User selected"),]

      if(nrow(tmp_test)>0) {
        not_in_leg=which(unique(tmp_test$Color) %!in% unique(tmp_test_ns$Color))
        if(length(not_in_leg)>0) {
          tmp_tmp_test=tmp_test[which(tmp_test$Color %in% not_in_leg),]

          fig <- fig %>%
            add_trace(data=tmp_tmp_test,
                      x = ~x,
                      y = ~y,
                      opacity = 0.9,
                      text =  ~ Text_ID,
                      key = ~ Text_ID,
                      color = ~Color,
                      symbol= ~Symbol,
                      # size = ~Width_border,
                      type = 'scatter',
                      mode = 'markers',
                      fill = ~'',
                      hovertemplate = paste('<b>%{text}</b>'),
                      marker = list(
                        line=list(color="black",
                                  width=tmp_tmp_test$Width_border[1])),
                      name= tmp_tmp_test$Color, showlegend=TRUE,
                      legendgroup=tmp_tmp_test$Color
            )
        }
        tmp_tmp_test=tmp_test[which(tmp_test$Color %!in% not_in_leg),]
        if(nrow(tmp_tmp_test)>0) {
          fig <- fig %>%
            add_trace(data=tmp_tmp_test,
                      x = ~x,
                      y = ~y,
                      opacity = 0.9,
                      text =  ~ Text_ID,
                      key = ~ Text_ID,
                      color = ~Color,
                      symbol= ~Symbol,
                      # size = ~Width_border,
                      type = 'scatter',
                      mode = 'markers',
                      fill = ~'',
                      hovertemplate = paste('<b>%{text}</b>'),
                      marker = list(
                        line=list(color="black",
                                  width=tmp_tmp_test$Width_border[1])),
                      name= tmp_tmp_test$Color, showlegend=FALSE,
                      legendgroup=tmp_tmp_test$Color
            )
        }

      }
      ##OLD
      # dfk_unselected <- data.frame(x = as.numeric(), y=as.numeric(),
      #                              Text_ID=as.character(), Color=as.character(),
      #                              Symbol=as.character())
      #
      # dfk_selected <- data.frame(x = as.numeric(), y=as.numeric(),
      #                            Text_ID=as.character(), Color=as.character(),
      #                            Width_border=as.numeric(),
      #                            Symbol=as.character())
      # dfk_selected_legend <- data.frame(x = as.numeric(), y=as.numeric(),
      #                            Text_ID=as.character(), Color=as.character(),
      #                            Width_border=as.numeric(),
      #                            Symbol=as.character())
      #
      # dfk_selected_counterpart <- data.frame(x = as.numeric(), y=as.numeric(),
      #                            Text_ID=as.character(), Color=as.character(),
      #                            Width_border=as.numeric(),
      #                            Symbol=as.character())
      #
      # for (Color in sort(unique(test$Color))){
      #   tmp_index=which(test$Color == Color)
      #   SELECTED=F
      #
      #
      #   if(any("User selected" %in% c(unique(test$Selected[tmp_index])))) {
      #     tmp_tmp_index=which(test$Selected[tmp_index] == "User selected")
      #     tmp_tmp_index_counterpart=which(test$Selected[tmp_index] == "Counterpart")
      #     tmp_tmp_index_no=intersect(which(test$Selected[tmp_index] != "User selected"), which(test$Selected[tmp_index] != "Counterpart"))
      #
      #     SELECTED=T
      #   }
      #
      #   if(any("Clones" %in% c(unique(test$Selected[tmp_index])))) {
      #     tmp_tmp_index=which(test$Selected[tmp_index] == "Clones")
      #     tmp_tmp_index_counterpart=which(test$Selected[tmp_index] == "Counterpart")
      #     tmp_tmp_index_no=intersect(which(test$Selected[tmp_index] != "User selected"), which(test$Selected[tmp_index] != "Counterpart"))
      #
      #     SELECTED=T
      #   }
      #   if(SELECTED) {
      #
      #     tmp_dfku <- data.frame(x = test$Dim_1[tmp_index[tmp_tmp_index_no]],
      #                           y=test$Dim_2[tmp_index[tmp_tmp_index_no]],
      #                           Text_ID=test$Text_ID[tmp_index[tmp_tmp_index_no]],
      #                           Color=test$Color[tmp_index[tmp_tmp_index_no]],
      #                           Symbol=test$Symbol[tmp_index[tmp_tmp_index_no]])
      #     dfk_unselected=rbind(dfk_unselected, tmp_dfku)
      #
      #     tmp_dfk=data.frame(x = test$Dim_1[tmp_index[tmp_tmp_index]],
      #                        y=test$Dim_2[tmp_index[tmp_tmp_index]],
      #                        Text_ID=test$Text_ID[tmp_index[tmp_tmp_index]],
      #                        Width_border=test$Width_border[tmp_index[tmp_tmp_index]],
      #                        Color=test$Color[tmp_index[tmp_tmp_index]],
      #                        Symbol=test$Symbol[tmp_index[tmp_tmp_index]])
      #
      #     if(nrow(tmp_dfk_u)==0){
      #       dfk_selected_legend <- rbind(dfk_selected_legend,tmp_dfk)
      #     } else {
      #       dfk_selected <- rbind(dfk_selected,tmp_dfk)
      #     }
      #
      #
      #
      #     tmp_dfk=data.frame(x = test$Dim_1[tmp_index[tmp_tmp_index_counterpart]],
      #                        y=test$Dim_2[tmp_index[tmp_tmp_index_counterpart]],
      #                        Text_ID=test$Text_ID[tmp_index[tmp_tmp_index_counterpart]],
      #                        Width_border=test$Width_border[tmp_index[tmp_tmp_index_counterpart]],
      #                        Color=test$Color[tmp_tmp_index_counterpart[tmp_tmp_index_counterpart]],
      #                        Symbol=test$Symbol[tmp_index[tmp_tmp_index_counterpart]])
      #
      #     dfk_selected_counterpart=rbind(dfk_selected_counterpart, tmp_dfk)
      #
      #   } else {
      #     tmp_dfk <- data.frame(x = test$Dim_1[tmp_index],
      #                       y=test$Dim_2[tmp_index],
      #                       Text_ID=test$Text_ID[tmp_index],
      #                       Color=test$Color[tmp_index],
      #                       Symbol=test$Symbol[tmp_index])
      #
      #     dfk_unselected=rbind(dfk_unselected, tmp_dfk)
      #
      #   }
      # }
      #
      # if(nrow(dfk_selected)>0){
      #   fig <- fig %>%
      #     add_trace(data=dfk_selected,
      #               x = ~x,
      #               y = ~y,
      #               opacity = 0.9,
      #               text =  ~ Text_ID,
      #               key = ~ Text_ID,
      #               color = ~Color,
      #               symbol= ~Symbol,
      #               # size = ~Width_border,
      #               type = 'scatter',
      #               mode = 'markers',
      #               fill = ~'',
      #               hovertemplate = paste('<b>%{text}</b>'),
      #               marker = list(
      #                 line=list(color="black",
      #                           width=dfk_selected$Width_border[1])),
      #               name= dfk_selected$Color, showlegend=FALSE,
      #               legendgroup=dfk_selected$Color
      #     )
      #
      #
      # }
      #
      # if(nrow(dfk_selected_legend)>0) {
      #   fig <- fig %>%
      #     add_trace(data=dfk_selected_legend,
      #               x = ~x,
      #               y = ~y,
      #               opacity = 0.9,
      #               text =  ~ Text_ID,
      #               key = ~ Text_ID,
      #               color = ~Color,
      #               symbol= ~Symbol,
      #               # size = ~Width_border,
      #               type = 'scatter',
      #               mode = 'markers',
      #               fill = ~'',
      #               hovertemplate = paste('<b>%{text}</b>'),
      #               marker = list(
      #                 line=list(color="black",
      #                           width=dfk_selected_legend$Width_border[1])),
      #               name= dfk_selected_legend$Color, showlegend=TRUE,
      #               legendgroup=dfk_selected_legend$Color
      #     )
      # }
      #
      # if(nrow(dfk_selected_legend)>0 || nrow(dfk_selected)>0) {
      #   fig <- fig %>%
      #     add_trace(data=dfk_selected_counterpart,
      #               x = ~x,
      #               y = ~y,
      #               opacity = 0.9,
      #               text =  ~ Text_ID,
      #               key = ~ Text_ID,
      #               color = ~Color,
      #               colors = Sel_colors,
      #               # size = ~Width_border,
      #               symbol= ~Symbol,
      #               type = 'scatter',
      #               mode = 'markers',
      #               fill = ~'',
      #               hovertemplate = paste('<b>%{text}</b>'),
      #               marker =  list(
      #                 line=list(color="#ff4d6d",
      #                           width=dfk_selected_counterpart$Width_border[1])),
      #               name= dfk_selected_counterpart$Color, showlegend=FALSE,
      #               legendgroup=dfk_selected_counterpart$Color
      #     )
      # }
      #
      # if(nrow(dfk_unselected)>0) {
      #   fig <- fig %>%
      #     add_trace(data=dfk_unselected,
      #               x = ~x,
      #               y = ~y,
      #               opacity = 0.9,
      #               text =  ~ Text_ID,
      #               key = ~ Text_ID,
      #               color = ~Color,
      #               colors = Sel_colors,
      #               # size = ~Width_border,
      #               type = 'scatter',
      #               mode = 'markers',
      #               fill = ~'',
      #               hovertemplate = paste('<b>%{text}</b>'),
      #               marker = list(sizemode = 'diameter'),
      #               name= dfk_unselected$Color,
      #               legendgroup= dfk_unselected$Color
      #     )
      # }

      fig <- fig %>%
        layout(title = 'Selection plot', plot_bgcolor = "#e5ecf6",legend = list(orientation = 'v',y=0),showlegend=T, xaxis = list(title =  if(input$Selection_plot_type == "PCA") {paste('Dim ', dim1,' (', 100*Selection_values$Variance_explained[dim1], "%)", sep="")} else if(input$Selection_plot_type == "UMAP") {"Dim  1"}  ),
               yaxis = list(title = if(input$Selection_plot_type == "PCA") {paste('Dim ', dim2,' (', 100*Selection_values$Variance_explained[dim2], "%)", sep="")} else if(input$Selection_plot_type == "UMAP") {"Dim 2"} ,dragmode = "lasso")   )%>%
        config(displaylogo = FALSE, modeBarButtonsToRemove = c("zoomIn2d", "zoomOut2d", "autoScale2d","pan2d", "hoverCompareCartesian"))

      ID_selected_values$subclones=event_data("plotly_selected")

      print("Figure")
      fig%>% toWebGL()
    }

  })


  output$Exploration_plot <- renderPlotly({
    if(!is.null(Exploration_values$rows) && length(Exploration_values$rows)>0 && !is.null(Exploration_values$columns) && length(Exploration_values$columns)>0) {
      print("Exploration plot")
      # Ex_colors=unlist(branded_colors[c(1:length(unique(Exploration_values$Scores$Color)))])
      # names(Ex_colors)=sort(unique(Exploration_values$Scores$Color))
      # print(Ex_colors)
      # Ex_colors_BORDER=c(Ex_colors, "#050609", "#ffb703")
      # names(Ex_colors_BORDER)=c(names(Ex_colors),"User selected", "Counterpart")
      # print(Ex_colors_BORDER)
      # Ex_colors_BORDER_width=c(rep(1, length(Ex_colors)),4,4)
      # names(Ex_colors_BORDER_width)=c(names(Ex_colors),"User selected", "Counterpart" )
      # Ex_colors_BORDER_symbol=c("circle", "star-diamond")
      # names(Ex_colors_BORDER_symbol)=c("Repertoire", "Reconstructed_germline")
      #
      # print(dim(Exploration_values$Scores))
      # print(colnames(Exploration_values$Scores))
      # print("Render expl values")
      # print(Exploration_values$Scores$Selected[1:20])
      # print(Exploration_values$Scores$Color[1:20])
      # print(Ex_colors[match(Exploration_values$Scores$Color[1:20], names(Ex_colors))])
      # print(Ex_colors_BORDER[match(Exploration_values$Scores$Selected[1:20], names(Ex_colors_BORDER))])
      # print("Colors the same?")
      # print(length(Ex_colors_BORDER[match(Exploration_values$Scores$Selected, names(Ex_colors_BORDER))]) == length(Ex_colors[match(Exploration_values$Scores$Color, names(Ex_colors))]) )
      # print(all(Ex_colors_BORDER[match(Exploration_values$Scores$Selected, names(Ex_colors_BORDER))] == Ex_colors[match(Exploration_values$Scores$Color, names(Ex_colors))]))
      # print(all(Ex_colors_BORDER[base::match(Exploration_values$Scores$Selected, names(Ex_colors_BORDER))] == Ex_colors[match(Exploration_values$Scores$Color, names(Ex_colors))]))
      # print(which(Ex_colors_BORDER[base::match(Exploration_values$Scores$Selected, names(Ex_colors_BORDER))] != Ex_colors[match(Exploration_values$Scores$Color, names(Ex_colors))]))
      # print(Ex_colors_BORDER[base::match(Exploration_values$Scores$Selected, names(Ex_colors_BORDER))][which(Ex_colors_BORDER[base::match(Exploration_values$Scores$Selected, names(Ex_colors_BORDER))] != Ex_colors[match(Exploration_values$Scores$Color, names(Ex_colors))])])
      # Exploration_values$Scores$Colors_all=Ex_colors[match(Exploration_values$Scores$Color, names(Ex_colors))]
      # Exploration_values$Scores$Colors_selection=Ex_colors_BORDER[base::match(Exploration_values$Scores$Selected, names(Ex_colors_BORDER))]
      # Exploration_values$Scores$Colors_selection_width=Ex_colors_BORDER_width[base::match(Exploration_values$Scores$Selected, names(Ex_colors_BORDER_width))]
      # print(head(Exploration_values$Scores))
      #
      # fig_ex=plot_ly(
      #   data = Exploration_values$Scores ,
      #   x =  ~ Dim_1,
      #   y = ~ Dim_2,
      #   color = ~ Color,
      #   colors = Ex_colors_BORDER,
      #   opacity = 0.7,
      #   text =  ~ Text_ID,
      #   key = ~ Text_ID,
      #   type = 'scatter',
      #   mode = 'markers',
      #   hovertemplate = paste('<b>%{text}</b>'),
      #
      #   marker = list(
      #     size = 7,
      #     symbol = Ex_colors_BORDER_symbol[base::match(Exploration_values$Scores$Seq_type,
      #                                                  names(Ex_colors_BORDER_symbol))],
      #     line= list(color=(Ex_colors_BORDER[match(Exploration_values$Scores$Selected, names(Ex_colors_BORDER))]),
      #                width=( Ex_colors_BORDER_width[match(Exploration_values$Scores$Selected, names(Ex_colors_BORDER_width))]))
      #
      #   )
      # )  %>%
      #   layout(title = 'Exploration plot', plot_bgcolor = "#e5ecf6",legend = list(orientation = 'h',y=-0.3),showlegend=T, xaxis = list(title = paste('Dim 1 (', Exploration_values$Variance_explained[1], "%)", sep="") ),
      #          yaxis = list(title = paste('Dim 2 (', Exploration_values$Variance_explained[2], "%)", sep=""),dragmode = "lasso")   )%>%
      #   config(displaylogo = FALSE, modeBarButtonsToRemove = c("zoomIn2d", "zoomOut2d", "autoScale2d","pan2d", "hoverCompareCartesian"))
      #
      #
      # ID_selected_values$subclones=event_data("plotly_selected")
      # fig_ex

      if(input$plot_color_expl == "Best_V"){
        print("VVV")
        Ex_colors=Big_mem_color_values$V
      } else if (input$plot_color_expl == "Best_J") {
        Ex_colors=Big_mem_color_values$J
      } else if(input$plot_color_expl == "Best_D") {

        Ex_colors=Big_mem_color_values$D

      } else if(input$plot_color_expl == "V_and_D_and_J") {
        print("V_and_D_and_J")
        Ex_colors=Big_mem_color_values$VDJ

      } else if(input$plot_color_expl == "V_and_J"){
        print("VV_and_JVV")
        Ex_colors=Big_mem_color_values$VJ
      } else {
        print("ADAD")
        print(input$plot_color_expl)
        # Sel_colors=unlist(branded_colors[c(1:length(unique(Selection_values$Scores$Color)))])
        Ex_colors=Ab_palette(list_values=unique(Exploration_values$Scores$Color),
                              vect_genes_comb=NA,
                              type_values="cualitative",
                              colorblind=F)
        names(Ex_colors)=sort(unique(Exploration_values$Scores$Color))
      }

      print(input$plot_color_expl)
      # if(input$plot_color_expl =="V_gene"){
      #   print(input$plot_color_expl)
      #   Ex_colors=Big_mem_color_values$V
      # } else if(input$plot_color_expl =="J_gene") {
      #   Ex_colors=Big_mem_color_values$J
      # } else if(input$plot_color_expl =="D_gene") {
      #   Ex_colors=Big_mem_color_values$D
      # }else if(input$plot_color_expl =="V_and_J") {
      #   Ex_colors=Big_mem_color_values$VJ
      # }else if(input$plot_color_expl =="V_and_D_and_J") {
      #   Ex_colors=Big_mem_color_values$VDJ
      # }



      Ex_colors_BORDER=c(Ex_colors, "#050609", "#ff4d6d", "#FFAD05")
      names(Ex_colors_BORDER)=c(names(Ex_colors),"User selected", "Counterpart", "Clones")
      Ex_colors_BORDER_width=c(rep(1, length(Ex_colors)),4,4,4)
      names(Ex_colors_BORDER_width)=c(names(Ex_colors),"User selected", "Counterpart", "Clones")
      Ex_colors_BORDER_symbol=c("circle", "star-diamond")
      names(Ex_colors_BORDER_symbol)=c("Repertoire", "Reconstructed_germline")
      test=Exploration_values$Scores

      print(input$Exploration_plot_type)
      if(input$Exploration_plot_type == "PCA") {

        test=Exploration_values$Scores
        print(colnames(test))
      } else if(input$Exploration_plot_type == "UMAP") {
        test=Exploration_values$UMAP
        print(dim(test))
        print(colnames(test))
        print(head(test))
      }

      print("Here")
      test$Colors=Ex_colors[match(test$Color, names(Ex_colors))]
      print(test$Colors)
      print("Here2")
      print(Ex_colors_BORDER)
      print(names(Ex_colors_BORDER))
      print(test$Selected)
      print("Here2.5")
      test$Color_border=Ex_colors_BORDER[match(test$Selected, names(Ex_colors_BORDER))]
      print("Here3")
      test$Width_border=Ex_colors_BORDER_width[match(test$Selected, names(Ex_colors_BORDER_width))]
      test$Symbol=Ex_colors_BORDER_symbol[match(test$Seq_type, names(Ex_colors_BORDER_symbol))]


      if(length(input$Exploration_plot_type_dim)<2){
        dim1=1
        dim2=2
      }else{
        dim1=as.numeric(input$Exploration_plot_type_dim[1])
        dim2=as.numeric(input$Exploration_plot_type_dim[2])
      }
      tmpDim_1=test[,which(colnames(test)==paste0("Dim_", dim1))]
      tmpDim_2=test[,which(colnames(test)==paste0("Dim_", dim2))]
      test$Dim_1=tmpDim_1
      test$Dim_2=tmpDim_2

      fig <- plot_ly(data=test, type = 'scatter', mode = 'markers',
                     colors= Ex_colors)   %>%
        config(
          toImageButtonOptions = list(
            format = "png",
            filename = "PCA_exploration_plot",
            width = 1400,
            height = 1000
          )
        )

      ##NEW
      test=test[(order(test$Color)),]

      tmp_test_ns=test[intersect(intersect(which(test$Selected != "Clones"),
                                           which(test$Selected != "User selected")),
                                 which(test$Selected != "Counterpart")),]

      fig <- fig %>%
        add_trace(data=tmp_test_ns,
                  x = ~Dim_1,
                  y = ~Dim_2,
                  opacity = 0.7,
                  text =  ~ Text_ID,
                  key = ~ Text_ID,
                  color = ~Color,
                  colors = Ex_colors,
                  # size = ~Width_border,
                  type = 'scatter',
                  mode = 'markers',
                  fill = ~'',
                  hovertemplate = paste('<b>%{text}</b>'),
                  marker = list(sizemode = 'diameter'),
                  name= tmp_test_ns$Color,
                  legendgroup= tmp_test_ns$Color
        )


      tmp_test=test[which(test$Selected == "Counterpart"),]

      if(nrow(tmp_test)>0) {
        fig <- fig %>%
          add_trace(data=tmp_test,
                    x = ~Dim_1,
                    y = ~Dim_2,
                    opacity = 0.7,
                    text =  ~ Text_ID,
                    key = ~ Text_ID,
                    color = ~Color,
                    colors = Ex_colors,
                    # size = ~Width_border,
                    symbol= ~Symbol,
                    type = 'scatter',
                    mode = 'markers',
                    fill = ~'',
                    hovertemplate = paste('<b>%{text}</b>'),
                    marker =  list(
                      line=list(color="#ff4d6d",
                                width=tmp_test$Width_border[1])),
                    name= tmp_test$Color, showlegend=FALSE,
                    legendgroup=tmp_test$Color
          )
      }

      tmp_test=test[which(test$Selected == "Clones"),]

      if(nrow(tmp_test)>0) {
        not_in_leg=which(unique(tmp_test$Color) %!in% unique(tmp_test_ns$Color))
        if(length(not_in_leg)>0) {
          tmp_tmp_test=tmp_test[which(tmp_test$Color %in% not_in_leg),]

          fig <- fig %>%
            add_trace(data=tmp_tmp_test,
                      x = ~x,
                      y = ~y,
                      opacity = 0.9,
                      text =  ~ Text_ID,
                      key = ~ Text_ID,
                      color = ~Color,
                      symbol= ~Symbol,
                      # size = ~Width_border,
                      type = 'scatter',
                      mode = 'markers',
                      fill = ~'',
                      hovertemplate = paste('<b>%{text}</b>'),
                      marker = list(
                        line=list(color="black",
                                  width=tmp_tmp_test$Width_border[1])),
                      name= tmp_tmp_test$Color, showlegend=TRUE,
                      legendgroup=tmp_tmp_test$Color
            )
        }
        tmp_tmp_test=tmp_test[which(tmp_test$Color %!in% not_in_leg),]
        if(nrow(tmp_tmp_test)>0) {
          fig <- fig %>%
            add_trace(data=tmp_tmp_test,
                      x = ~x,
                      y = ~y,
                      opacity = 0.9,
                      text =  ~ Text_ID,
                      key = ~ Text_ID,
                      color = ~Color,
                      symbol= ~Symbol,
                      # size = ~Width_border,
                      type = 'scatter',
                      mode = 'markers',
                      fill = ~'',
                      hovertemplate = paste('<b>%{text}</b>'),
                      marker = list(
                        line=list(color="black",
                                  width=tmp_tmp_test$Width_border[1])),
                      name= tmp_tmp_test$Color, showlegend=FALSE,
                      legendgroup=tmp_tmp_test$Color
            )
        }

      }

      tmp_test=test[which(test$Selected == "User selected"),]

      if(nrow(tmp_test)>0) {
        not_in_leg=which(unique(tmp_test$Color) %!in% unique(tmp_test_ns$Color))
        if(length(not_in_leg)>0) {
          tmp_tmp_test=tmp_test[which(tmp_test$Color %in% not_in_leg),]

          fig <- fig %>%
            add_trace(data=tmp_tmp_test,
                      x = ~x,
                      y = ~y,
                      opacity = 0.9,
                      text =  ~ Text_ID,
                      key = ~ Text_ID,
                      color = ~Color,
                      symbol= ~Symbol,
                      # size = ~Width_border,
                      type = 'scatter',
                      mode = 'markers',
                      fill = ~'',
                      hovertemplate = paste('<b>%{text}</b>'),
                      marker = list(
                        line=list(color="black",
                                  width=tmp_tmp_test$Width_border[1])),
                      name= tmp_tmp_test$Color, showlegend=TRUE,
                      legendgroup=tmp_tmp_test$Color
            )
        }
        tmp_tmp_test=tmp_test[which(tmp_test$Color %!in% not_in_leg),]
        if(nrow(tmp_tmp_test)>0) {
          fig <- fig %>%
            add_trace(data=tmp_tmp_test,
                      x = ~x,
                      y = ~y,
                      opacity = 0.9,
                      text =  ~ Text_ID,
                      key = ~ Text_ID,
                      color = ~Color,
                      symbol= ~Symbol,
                      # size = ~Width_border,
                      type = 'scatter',
                      mode = 'markers',
                      fill = ~'',
                      hovertemplate = paste('<b>%{text}</b>'),
                      marker = list(
                        line=list(color="black",
                                  width=tmp_tmp_test$Width_border[1])),
                      name= tmp_tmp_test$Color, showlegend=FALSE,
                      legendgroup=tmp_tmp_test$Color
            )
        }

      }
      # for (Color in unique(test$Color)){
      #   tmp_index=which(test$Color == Color)
      #   SELECTED=F
      #
      #   if(any("User selected" %in% c(unique(test$Selected[tmp_index])))) {
      #     tmp_tmp_index=which(test$Selected[tmp_index] == "User selected")
      #     tmp_tmp_index_counterpart=which(test$Selected[tmp_index] == "Counterpart")
      #     tmp_tmp_index_no=intersect(which(test$Selected[tmp_index] != "User selected"), which(test$Selected[tmp_index] != "Counterpart"))
      #
      #     SELECTED=T
      #   }
      #
      #   # if(any("Clones" %in% c(unique(test$Selected[tmp_index])))) {
      #   #   tmp_tmp_index=which(test$Selected[tmp_index] == "Clones")
      #   #   tmp_tmp_index_counterpart=which(test$Selected[tmp_index] == "Counterpart")
      #   #   tmp_tmp_index_no=intersect(which(test$Selected[tmp_index] != "User selected"), which(test$Selected[tmp_index] != "Counterpart"))
      #   #
      #   #   SELECTED=T
      #   # }
      #   if(SELECTED) {
      #     dfk <- data.frame(x = test$Dim_1[tmp_index[tmp_tmp_index]], y=test$Dim_2[tmp_index[tmp_tmp_index]], Text_ID=test$Text_ID[tmp_index[tmp_tmp_index]],
      #                       Width_border=test$Width_border[tmp_index[tmp_tmp_index]], Color=test$Color[tmp_index[tmp_tmp_index]],
      #                       Symbol=test$Symbol[tmp_index[tmp_tmp_index]])
      #     fig <- fig %>%
      #       add_trace(data=dfk,
      #                 x = ~x,
      #                 y = ~y,
      #                 opacity = 0.9,
      #                 text =  ~ Text_ID,
      #                 key = ~ Text_ID,
      #                 color = ~Color,
      #                 colors = Ex_colors,
      #                 symbol= ~Symbol,
      #                 # size = ~Width_border,
      #                 type = 'scatter',
      #                 mode = 'markers',
      #                 fill = ~'',
      #                 hovertemplate = paste('<b>%{text}</b>'),
      #                 marker = list(
      #                   line=list(color="black",
      #                             width=dfk$Width_border[1])),
      #                 name= Color, showlegend=FALSE, legendgroup=Color
      #       )
      #
      #     dfk <- data.frame(x = test$Dim_1[tmp_index[tmp_tmp_index_no]], y=test$Dim_2[tmp_index[tmp_tmp_index_no]], Text_ID=test$Text_ID[tmp_index[tmp_tmp_index_no]],
      #                       Width_border=test$Width_border[tmp_index[tmp_tmp_index_no]], Color=test$Color[tmp_index[tmp_tmp_index_no]],
      #                       Symbol=test$Symbol[tmp_index[tmp_tmp_index_no]])
      #     fig <- fig %>%
      #       add_trace(data=dfk,
      #                 x = ~x,
      #                 y = ~y,
      #                 opacity = 0.9,
      #                 text =  ~ Text_ID,
      #                 key = ~ Text_ID,
      #                 color = ~Color,
      #                 colors = Ex_colors,
      #                 # size = ~Width_border,
      #                 symbol= ~Symbol,
      #                 type = 'scatter',
      #                 mode = 'markers',
      #                 fill = ~'',
      #                 hovertemplate = paste('<b>%{text}</b>'),
      #                 marker = list(sizemode = 'diameter'),
      #                 name= Color, legendgroup=Color
      #       )
      #
      #     dfk <- data.frame(x = test$Dim_1[tmp_index[tmp_tmp_index_counterpart]], y=test$Dim_2[tmp_index[tmp_tmp_index_counterpart]], Text_ID=test$Text_ID[tmp_index[tmp_tmp_index_counterpart]],
      #                       Width_border=test$Width_border[tmp_index[tmp_tmp_index_counterpart]], Color=test$Color[tmp_tmp_index_counterpart[tmp_tmp_index_counterpart]],
      #                       Symbol=test$Symbol[tmp_index[tmp_tmp_index_counterpart]])
      #     fig <- fig %>%
      #       add_trace(data=dfk,
      #                 x = ~x,
      #                 y = ~y,
      #                 opacity = 0.9,
      #                 text =  ~ Text_ID,
      #                 key = ~ Text_ID,
      #                 color = ~Color,
      #                 colors = Ex_colors,
      #                 # size = ~Width_border,
      #                 symbol= ~Symbol,
      #                 type = 'scatter',
      #                 mode = 'markers',
      #                 fill = ~'',
      #                 hovertemplate = paste('<b>%{text}</b>'),
      #                 marker =  list(
      #                   line=list(color="#ff4d6d",
      #                             width=dfk$Width_border[1])),
      #                 name= Color, showlegend=FALSE, legendgroup=Color
      #       )
      #   } else {
      #     dfk <- data.frame(x = test$Dim_1[tmp_index], y=test$Dim_2[tmp_index], Text_ID=test$Text_ID[tmp_index],
      #                       Width_border=test$Width_border[tmp_index], Color=test$Color[tmp_index],
      #                       Symbol=test$Symbol[tmp_index])
      #     fig <- fig %>%
      #       add_trace(data=dfk,
      #                 x = ~x,
      #                 y = ~y,
      #                 opacity = 0.9,
      #                 text =  ~ Text_ID,
      #                 key = ~ Text_ID,
      #                 color = ~Color,
      #                 colors = Ex_colors,
      #                 # size = ~Width_border,
      #                 symbol= ~Symbol,
      #                 type = 'scatter',
      #                 mode = 'markers',
      #                 fill = ~'',
      #                 hovertemplate = paste('<b>%{text}</b>'),
      #                 marker = list(sizemode = 'diameter'),
      #                 name= Color, legendgroup=Color
      #       )
      #
      #
      #   }
      # }
      fig <- fig %>%
        layout(title = 'Exploration plot', plot_bgcolor = "#e5ecf6",legend = list(orientation = 'v',y=-0),showlegend=T, xaxis = list(title = if(input$Exploration_plot_type == "PCA") {paste('Dim ', dim1,' (', 100*Exploration_values$Variance_explained[dim1], "%)", sep="")} else if(input$Exploration_plot_type == "UMAP") {"Dim  1"} ),
               yaxis = list(title = if(input$Exploration_plot_type == "PCA") {paste('Dim ', dim2,' (', 100*Exploration_values$Variance_explained[dim2], "%)", sep="")} else if(input$Exploration_plot_type == "UMAP") {"Dim  2"},dragmode = "lasso")   )%>%
        config(displaylogo = FALSE, modeBarButtonsToRemove = c("zoomIn2d", "zoomOut2d", "autoScale2d","pan2d", "hoverCompareCartesian"))

      ID_selected_values$subclones=event_data("plotly_selected")

      fig%>% toWebGL()

    }

  })
  output$Selection_plot_error <- renderUI({
    if(!is.null(Selection_values$rows) && length(Selection_values$rows)>0 && !is.null(Selection_values$columns) && length(Selection_values$columns)>0) {

    } else{
      HTML('<p style="color:red">ERROR! <br>
              Your selection lead to no rows and/or columns with this analysis. <br>
              Change some parameters!
             </p>')
    }

  })
  output$Exploration_plot_error <- renderUI({
    if(!is.null(Exploration_values$rows) && length(Exploration_values$rows)>0 && !is.null(Exploration_values$columns) && length(Exploration_values$columns)>0) {

    } else{
      HTML('<p style="color:red">ERROR! <br>
              Select any rows and/or columns<br>
              Or perhaps your selection lead to an error with this analysis
             </p>')
    }

  })

  output$Plot_color <- renderUI({
    if(is.null(Selection_values$rows) && length(Selection_values$rows)>0 && !is.null(Selection_values$columns) && length(Selection_values$columns)>0) {
      selectInput("plot_color","Color by", choices = c("Sample",	"Patient",	"Group",	"Subgroup",	"Chain"), selected="Group")
    } else{
      sample_info<-sample_info_react$table
      choices_color=c(colnames(sample_info))
      names(choices_color)=c(colnames(sample_info))
      if ("Chain" %in% colnames(Big_mem_values$Short_DF) && "Chain" %!in% choices_color) {
        tmp_names_color=names(choices_color)
        choices_color=c(choices_color, "Chain")
        names(choices_color)=c(tmp_names_color, "Ig Chain")
      }
      tmp_names_choices_color=names(choices_color)
      choices_color=c(choices_color, "Best_V", "Best_J", "V_and_J")
      names(choices_color)=c(tmp_names_choices_color, "Main V gene isoform", "Main J gene isoform", "Main V and J gene isoforms")
      choices_color=choices_color[which(choices_color != "Additional_info")]
      choices_color=choices_color[which(choices_color != "Folder")]
      choices_color=choices_color[which(choices_color != "Filename")]
      # choices_color=choices_color[which(choices_color != "V_and_J")]



      if ("D_region" %in% colnames(Big_mem_values$Short_DF)) {
        tmp_names_color=names(choices_color)
        choices_color=c(choices_color, "Best_D")
        names(choices_color)=c(tmp_names_color, "Main D gene isoform")
      }
      if ("V_and_D_and_J" %in% colnames(Big_mem_values$Short_DF)) {
        tmp_names_color=names(choices_color)
        choices_color=c(choices_color, "V_and_D_and_J")
        names(choices_color)=c(tmp_names_color, "Main V, D and J gene isoforms")
      }
      if ("C_region" %in% colnames(Big_mem_values$Short_DF)) {
        tmp_names_color=names(choices_color)
        choices_color=c(choices_color, "C_region")
        names(choices_color)=c(tmp_names_color, "Constant region")
      }
      if ("Clone_ID" %in% colnames(Big_mem_values$Short_DF)) {
        tmp_names_color=names(choices_color)
        choices_color=c(choices_color, "Clone_ID")
        names(choices_color)=c(tmp_names_color, "Clone ID")
      }

      if ("Dominance" %in% colnames(Big_mem_values$Short_DF)) {
        tmp_names_color=names(choices_color)
        choices_color=c(choices_color, "Dominance")
        names(choices_color)=c(tmp_names_color, "Dominance")
      }
      # if ("Cell_ID" %in% colnames(Big_mem_values$Short_DF)) {
      #   tmp_names_color=names(choices_color)
      #   choices_color=c(choices_color, "Cell_ID")
      #   names(choices_color)=c(tmp_names_color, "Cell ID")
      # }
      # print(choices_color)
      # choices_color=choices_color[match(choices_color,c("Sample",	"Patient",	"Group",	"Subgroup",	"Chain", "Best_V", "Best_J", "V_and_J", "D_region","C_region","Clone_ID","Cell_ID"))]
      # print(choices_color)

      selectInput("plot_color","Color by", choices = choices_color, selected="Group")
    }

  })


  output$Plot_color_expl <- renderUI({
    if(is.null(Exploration_values$rows) && length(Exploration_values$rows)>0 && !is.null(Exploration_values$columns) && length(Exploration_values$columns)>0) {
      selectInput("plot_color_expl","Color by", choices = c("Sample",	"Patient",	"Group",	"Subgroup",	"Chain"), selected="Group")
    } else{
      sample_info<-sample_info_react$table

      choices_color=c(colnames(sample_info))
      names(choices_color)=c(colnames(sample_info))
      if ("Chain" %in% colnames(Big_mem_values$Short_DF) && "Chain" %!in% choices_color) {
        tmp_names_color=names(choices_color)
        choices_color=c(choices_color, "Chain")
        names(choices_color)=c(tmp_names_color, "Ig Chain")
      }
      tmp_names_choices_color=names(choices_color)
      choices_color=c(choices_color, "Best_V", "Best_J", "V_and_J")
      names(choices_color)=c( tmp_names_choices_color, "Main V gene isoform", "Main J gene isoform", "Main V and J gene isoforms")
      choices_color=choices_color[which(choices_color != "Additional_info")]
      choices_color=choices_color[which(choices_color != "Folder")]
      choices_color=choices_color[which(choices_color != "Filename")]
      # choices_color=choices_color[which(choices_color != "V_and_J")]


      if ("D_region" %in% colnames(Big_mem_values$Short_DF)) {
        tmp_names_color=names(choices_color)
        choices_color=c(choices_color, "Best_D")
        names(choices_color)=c(tmp_names_color, "Main D gene isoform")
      }
      if ("V_and_D_and_J" %in% colnames(Big_mem_values$Short_DF)) {
        tmp_names_color=names(choices_color)
        choices_color=c(choices_color, "V_and_D_and_J")
        names(choices_color)=c(tmp_names_color, "Main V, D and J gene isoforms")
      }
      if ("C_region" %in% colnames(Big_mem_values$Short_DF)) {
        tmp_names_color=names(choices_color)
        choices_color=c(choices_color, "C_region")
        names(choices_color)=c(tmp_names_color, "Constant region")
      }
      if ("Clone_ID" %in% colnames(Big_mem_values$Short_DF)) {
        tmp_names_color=names(choices_color)
        choices_color=c(choices_color, "Clone_ID")
        names(choices_color)=c(tmp_names_color, "Clone ID")
      }
      if ("Dominance" %in% colnames(Big_mem_values$Short_DF)) {
        tmp_names_color=names(choices_color)
        choices_color=c(choices_color, "Dominance")
        names(choices_color)=c(tmp_names_color, "Dominance")
      }
      # if ("Cell_ID" %in% colnames(Big_mem_values$Short_DF)) {
      #   tmp_names_color=names(choices_color)
      #   choices_color=c(choices_color, "Cell_ID")
      #   names(choices_color)=c(tmp_names_color, "Cell ID")
      # }



      # choices_color=choices_color[match(choices_color,c("Sample",	"Patient",	"Group",	"Subgroup",	"Chain", "Best_V", "Best_J", "V_and_J", "D_region","C_region","Clone_ID","Cell_ID"))]
      # tst[match(c("D","AD","A","C")[which(c("D","AD","A","C")  %in%tst)],tst)]
      selectInput("plot_color_expl","Color by", choices = choices_color, selected="Group")
    }

  })

  toListenColorsSelection <- reactive({
    list(input$plot_color,input$plot_color_expl,ID_selected_values$subclones,ID_selected_values$clones,
         Exploration_values$rows, Exploration_values$columns,
         Selection_values$rows, Selection_values$columns, Big_mem_values$Run)
  })


  observeEvent(toListenColorsSelection(),{
    if(!is.null(input$plot_color_expl) && !is.null(input$plot_color)){
      print("Evento")
      print(input$plot_color_expl)
      tmp_sp=(unlist(Big_mem_values$Short_DF[Exploration_values$rows,get(input$plot_color_expl)]))

      if(input$plot_color_expl %in% c("Best_V", "Best_J")) {
        # print("Test_1")
        # print(tmp_sp[1:10])
        # tmp_sp=sapply(tmp_sp, function(z) strsplit(z,split="*", fixed=T)[[1]][1])
        # print("Test_2")
        print(tmp_sp[1:10])
      } else if (input$plot_color_expl == "V_and_J") {
        # tmp_sp=sapply(tmp_sp, function(z) paste(strsplit(z,split="*", fixed=T)[[1]][1],
        #                                         strsplit(strsplit(z,split="_", fixed=T)[[1]][2], split="*", fixed=T)[[1]][1],
        #                                         sep=" & "))
      }

      print(tmp_sp[1:10])
      tmp_sp[which(is.na(tmp_sp))]="NotAvailable"
      tmp_sp[which((tmp_sp)=="  ")]="NotAvailable"
      tmp_sp[which((tmp_sp)==" ")]="NotAvailable"
      tmp_sp[which((tmp_sp)=="")]="NotAvailable"
      tmp_sp[which((tmp_sp)=="NotAvailable")]="NotAvailable"

      tmp_sp=unname(tmp_sp)
      Exploration_values$Scores$Color =tmp_sp

      if(as.numeric(unname(strsplit(as.character(print(get_ram(), unit_system = "metric")), split=" ")[[1]][1])) >   round(length(Exploration_values$rows)*length(Exploration_values$columns)*8/2^{20}/1024, 2)) {
        Exploration_values$UMAP$Color=tmp_sp
      }


      old_sp=NULL
      all_sp=NULL
      if (!is.null(ID_selected_values$subclones)) {
        print("t3esting")
        print(tmp_sp)
        old_sp=tmp_sp

        tmp_sp[Exploration_values$Scores$Text_ID %in% ID_selected_values$subclones$key] = "User selected"



        counterpart_start=paste(sapply(ID_selected_values$subclones$key, function(z) strsplit(z, split="_&_")[[1]][1]), "_&_", sep="")
        tmp_txt_id=paste(sapply(Exploration_values$Scores$Text_ID, function(z) strsplit(z, split="_&_")[[1]][1]), "_&_", sep="")
        tmp_sp[intersect(which(tmp_txt_id %in% counterpart_start), which(tmp_sp != "User selected"))] = "Counterpart"
        print(tmp_sp)
      }
      print("t4esting")
      print(tmp_sp)
      tmp_sp=unname(tmp_sp)
      print(tmp_sp)
      Exploration_values$Scores$Selected=tmp_sp

      if(as.numeric(unname(strsplit(as.character(print(get_ram(), unit_system = "metric")), split=" ")[[1]][1])) >   round(length(Exploration_values$rows)*length(Exploration_values$columns)*8/2^{20}/1024, 2)) {
        Exploration_values$UMAP$Selected=tmp_sp
      }

      print("A ver")
      print(      all(Exploration_values$Scores$Color == Exploration_values$Scores$Selected))
      print("IDs")
      print(ID_selected_values$subclones$key)
      print("Text_ID")
      print(Exploration_values$Scores$Text_ID[which(Exploration_values$Scores$Text_ID %in% ID_selected_values$subclones$key)])
      print(old_sp[which(Exploration_values$Scores$Text_ID %in% ID_selected_values$subclones$key)])
      print(tmp_sp[which(Exploration_values$Scores$Text_ID %in% ID_selected_values$subclones$key)])

      tmp_sp=(unlist(Big_mem_values$Short_DF[Selection_values$rows,get(input$plot_color)]))
      print("Testete")
      print(length(tmp_sp))
      # if(input$plot_color %in% c("Best_V", "Best_J")) {
      #   tmp_sp=sapply(tmp_sp, function(z) strsplit(z,split="*", fixed=T)[[1]][1])
      # } else if (input$plot_color == "V_and_J") {
      #   tmp_sp=sapply(tmp_sp, function(z) paste(strsplit(z,split="*", fixed=T)[[1]][1],
      #                                           strsplit(strsplit(z,split="_", fixed=T)[[1]][2], split="*", fixed=T)[[1]][1],
      #                                           sep=" & "))
      # }
      print("Mmmmmmmmmmmmmmmmmmmmmmmmmmeow")
      tmp_sp[which(is.na(tmp_sp))]="Not specified"
      tmp_sp[which((tmp_sp)=="  ")]="Not specified"
      tmp_sp[which((tmp_sp)==" ")]="Not specified"
      tmp_sp[which((tmp_sp)=="")]="Not specified"
      tmp_sp[which((tmp_sp)=="NotAvailable")]="Not specified"

      print(length(tmp_sp))
      print(nrow(Selection_values$Scores))
      print(colnames(Selection_values$Scores))
      print(length(Selection_values$Scores$Color))
      Selection_values$Scores$Color=tmp_sp

      if(as.numeric(unname(strsplit(as.character(print(get_ram(), unit_system = "metric")), split=" ")[[1]][1])) >   round(length(Selection_values$rows)*length(Selection_values$columns)*8/2^{20}/1024, 2)) {
        Selection_values$UMAP$Color=tmp_sp
      }

      if (!is.null(ID_selected_values$subclones)) {
        tmp_sp[Selection_values$Scores$Text_ID %in% ID_selected_values$subclones$key] = "User selected"

        all_sp=rep("Non-selected", nrow(Big_mem_values$Short_DF))
        all_sp[Big_mem_values$Short_DF$Text_ID %in% ID_selected_values$subclones$key] = "User selected"

        counterpart_start=paste(sapply(ID_selected_values$subclones$key, function(z) strsplit(z, split="_&_")[[1]][1]), "_&_", sep="")
        tmp_txt_id=paste(sapply(Selection_values$Scores$Text_ID, function(z) strsplit(z, split="_&_")[[1]][1]), "_&_", sep="")
        tmp_sp[intersect(which(tmp_txt_id %in% counterpart_start), which(tmp_sp != "User selected"))] = "Counterpart"

      }
      print("Mmmmmmmmmmmmmmmmmmmmmmmmmmeow23")
      all_sp=unname(all_sp)
      Big_mem_values$Short_DF$Selected=all_sp
      print("Mmmmmmmmmmmmmmmmmmmmmmmmmmeow232")
      if (!is.null(ID_selected_values$clones)) {

        print("Wow clones working")
        print(ID_selected_values$clones)
        if(input$clonal_group=="Clone_ID") {
          clonepart_start=(Big_mem_values$Short_DF)[intersect(Selection_values$rows, which((Big_mem_values$Short_DF)[,
                                                                                                                     get("Clone_ID")] %in% ID_selected_values$clones$key)),get("ID")]

        } else {
          clonepart_start=(Big_mem_values$Short_DF)[intersect(Selection_values$rows, which((Big_mem_values$Short_DF)[,
                                                                                                                     get(paste("Clone",input$clonal_level_group,input$clonal_region_group, input$clonal_group,  input$identity_clonal_group,"simil", sep="_"))] %in% ID_selected_values$clones$key)),get("ID")]

        }
        tmp_txt_id=sapply(Selection_values$Scores$Text_ID, function(z) strsplit(z, split="_&_")[[1]][1])


        tmp_sp[which(tmp_txt_id %in% clonepart_start)] = "Clones"

        print(length(which(tmp_sp == "Clones")))
      }

      print(length(Selection_values$Scores$Selected))
      Selection_values$Scores$Selected=tmp_sp
      if(as.numeric(unname(strsplit(as.character(print(get_ram(), unit_system = "metric")), split=" ")[[1]][1])) >   round(length(Selection_values$rows)*length(Selection_values$columns)*8/2^{20}/1024, 2)) {
        Selection_values$UMAP$Selected=tmp_sp

      }

    }

  })

  observeEvent(ID_selected_values$subclones, {
    if (!is.null(ID_selected_values$subclones)) {
      # ID_selected_values$clones=NULL

    }
  })
  observeEvent(ID_selected_values$clones, {
    if (!is.null(ID_selected_values$clones)) {
      # ID_selected_values$subclones=NULL
    }
  })

  merged_data <- eventReactive(ID_selected_values$subclones, {

    if (!is.null(ID_selected_values$subclones) && !is.null(ID_selected_values$clones)) {
      keep_cols = colnames(Big_mem_values$Short_DF)
      print("a ver..1.")
      print(which(Selection_values$Scores$Selected == "Clones"))
      Big_mem_values$Short_DF[unique(c(match(ID_selected_values$subclones$key, paste(unlist(Big_mem_values$Short_DF[,get("ID")]),unlist(Big_mem_values$Short_DF[,get("Sequence_type")]), sep="_&_")),
                                       Selection_values$rows[which(Selection_values$Scores$Selected == "Clones")]
      )),..keep_cols]

    } else if (!is.null(ID_selected_values$subclones)) {
      keep_cols = colnames(Big_mem_values$Short_DF)
      print("a ver..2.")
      print(which(Selection_values$Scores$Selected == "Clones"))
      Big_mem_values$Short_DF[match(ID_selected_values$subclones$key, paste(unlist(Big_mem_values$Short_DF[,get("ID")]),unlist(Big_mem_values$Short_DF[,get("Sequence_type")]), sep="_&_")),..keep_cols]
    } else if (!is.null(ID_selected_values$clones)) {
      keep_cols = colnames(Big_mem_values$Short_DF)
      print("a ver..3.")
      print(which(Selection_values$Scores$Selected == "Clones"))
      Big_mem_values$Short_DF[ Selection_values$rows[which(Selection_values$Scores$Selected == "Clones")],..keep_cols]
    }


  }) #eventReactive


  output$Ab_table =DT::renderDataTable({
    rendered_table <- merged_data()
    DT::datatable(rendered_table,extensions = 'Buttons',options = list(
      dom = 'Bfrtip',
      buttons = c('copy', 'csv')))
  }, options = list(
    dom = 'Bfrtip',
    buttons = c('copy', 'csv'), scrollX = TRUE))


  observeEvent(input$Update_selection,{
    tmp_rows_cols=filter_merged(Big_mem_values$Big_DF,
                                Big_mem_values$Short_DF,
                                Big_mem_values$Header,
                                "Reconstructed germline" %in%  input$use_what,
                                "Repertoire" %in%  input$use_what,
                                "Productive" %in% input$use_productive_or_not,
                                "Non-productive" %in% input$use_productive_or_not,
                                input$my_regions,
                                input$my_var_elements,
                                input$my_vars,
                                input$my_vartypes,
                                input$use_sharedVDJ,
                                input$VJ_included,
                                input$groups_selected,
                                input$group_A,
                                input$group_B,
                                input$group_C,
                                input$use_univlog,
                                input$samples_selected,
                                input$exclude_variables,
                                input$pval_type,
                                input$pval_cutoff,
                                input$estimate_cutoff,
                                input$number_selected_vars,
                                input$VJ_deselected,
                                input$VDJ_normalized_per_size,
                                input$Rmut_filter[1],
                                input$Rmut_filter[2],
                                input$work_as_categories,
                                input$VDJ_maximize_clones,
                                input$VDJ_normalized_per_sample,
                                input$my_clone_def)
    Selection_values$rows=tmp_rows_cols$ROWS
    Selection_values$columns=tmp_rows_cols$COLUMNS
    Selection_values$Parameters = list(input$use_UMAP)
  })

  observeEvent(input$Update_exploration,{
    tmp_rows_cols=filter_merged(Big_mem_values$Big_DF,  Big_mem_values$Short_DF, Big_mem_values$Header,"Reconstructed germline" %in%  input$use_what, "Repertoire" %in%  input$use_what, "Productive" %in% input$use_productive_or_not,"Non-productive" %in% input$use_productive_or_not,
                                input$my_regions, input$my_var_elements, input$my_vars, input$my_vartypes, input$use_sharedVDJ, input$VJ_included, input$groups_selected, input$group_A, input$group_B,input$group_C, input$use_univlog,input$samples_selected, input$exclude_variables, input$pval_type, input$pval_cutoff, input$estimate_cutoff, input$number_selected_vars,
                                input$VJ_deselected, input$VDJ_normalized_per_size, input$Rmut_filter[1], input$Rmut_filter[2], input$work_as_categories, input$VDJ_maximize_clones, input$VDJ_normalized_per_sample, input$my_clone_def)
    Exploration_values$rows=tmp_rows_cols$ROWS
    Exploration_values$columns=tmp_rows_cols$COLUMNS
    Exploration_values$Parameters = list(input$use_UMAP)
  })






  output$Group_selection=renderUI({
    choices = c("Group", "Patient")
    if(!is.null(Exploration_values$rows) && length(Exploration_values$rows)>0 && !is.null(Exploration_values$columns) && length(Exploration_values$columns)>0) {
      if("Dominance" %in% colnames(Big_mem_values$Short_DF)) {
        choices=c(choices, "Dominance")
      }
      if("Selected" %in% colnames(Big_mem_values$Short_DF)) {
        choices=c(choices, "Selected")
      }
    }
    selectInput("groups_selected","Work only with", choices = choices, selected="Group")
  })
  # observeEvent(input$Move_to_analysis_real, {
  #
  #   info=merge_bigmems(folder_values$Featured)
  #   Big_mem_values$Header=info[[1]]
  #   Big_mem_values$Short_DF=info[[2]]
  #   Big_mem_values$Big_DF=info[[3]]
  #   rm(info)
  #   Exploration_values$rows=c(1:nrow(Big_mem_values$Big_DF))
  #   Exploration_values$columns=c(1:ncol(Big_mem_values$Big_DF))
  #   Selection_values$rows=c(1:nrow(Big_mem_values$Big_DF))
  #   Selection_values$columns=c(1:ncol(Big_mem_values$Big_DF))
  #
  # })


  # in server.R create reactiveVal
  Sample_selection_current_selection <- reactiveVal(NULL)

  # now store your current selection in the reactive value
  observeEvent(ignoreInit = TRUE,list(input$samples_selected, input$Rmut_filter), {

    filter_muts=intersect(which(Big_mem_values$Big_DF[,which(Big_mem_values$Header == "AA_Whole_Replacement_muts")] >= input$Rmut_filter[1]),
                          which(Big_mem_values$Big_DF[,which(Big_mem_values$Header == "AA_Whole_Replacement_muts")] <= input$Rmut_filter[2]))
    Big_mem_values$VJs=sort(unique(Big_mem_values$Short_DF$V_and_J[intersect(intersect(which(Big_mem_values$Short_DF$Patient_Sample %in% input$samples_selected),
                                                                                       which(Big_mem_values$Short_DF$Sequence_type %in% input$use_what)),
                                                                             filter_muts)]))
    counts_VJs=table(Big_mem_values$Short_DF$V_and_J[intersect(intersect(which(Big_mem_values$Short_DF$Patient_Sample %in% input$samples_selected), which(Big_mem_values$Short_DF$Sequence_type %in% input$use_what)),
                                                               filter_muts)])

    # counts_VJs=table(Big_mem_values$Short_DF$V_and_J[(intersect(which(Big_mem_values$Short_DF$Patient_Sample %in% input$samples_selected), which(Big_mem_values$Short_DF$Sequence_type %in% input$use_what))
    #                                                            )])
    names(Big_mem_values$VJs)= sapply(Big_mem_values$VJs, function(z) paste(z,  counts_VJs[which(names(counts_VJs)==z)], sep=" - Counts:") )
    Sample_selection_current_selection(input$samples_selected)
  })

  output$Sample_selection=renderUI({
    choices = c("RANDOMLETTERS")
    if(!is.null(Exploration_values$rows) && length(Exploration_values$rows)>0 && !is.null(Exploration_values$columns) && length(Exploration_values$columns)>0) {
      choices=Big_mem_values$Patient_Sample
    }
    if(!is.null(Sample_selection_current_selection()) && Sample_selection_current_selection()[1] != "") {
      tmp_selection=Sample_selection_current_selection()
    } else {
      tmp_selection=choices
    }
    if(length(tmp_selection)==0){
      tmp_selection=choices[1]
    }
    selectizeInput(
      'samples_selected', 'Select the sample(s)', choices = choices, selected = tmp_selection, multiple = TRUE,
      options = list(plugins= list('remove_button')))
  })


  # in server.R create reactiveVal
  VJ_selection_current_selection <- reactiveVal(NULL)

  # now store your current selection in the reactive value
  observeEvent(input$VJ_deselected, {
    VJ_selection_current_selection(input$VJ_deselected)
  })



  output$VJ_selection=renderUI({
    choices = c("RANDOMLETTERS")
    if(!is.null(Exploration_values$rows) && length(Exploration_values$rows)>0 && !is.null(Exploration_values$columns) && length(Exploration_values$columns)>0) {
      choices=Big_mem_values$VJs
    }
    if(!is.null(VJ_selection_current_selection()) && VJ_selection_current_selection()[1] != "") {
      tmp_selection=VJ_selection_current_selection()
      tmp_selection=tmp_selection[which(tmp_selection %in% Big_mem_values$VJs)]
    } else {
      tmp_selection=NULL
    }

    # selectizeInput(
    #   'VJ_selected', 'Include combination(s)', choices = choices, selected = tmp_selection, multiple = TRUE,
    #   options = list(plugins= list('remove_button')))
    if(length(tmp_selection)==length(choices)) {
      tmp_selection=choices[2:length(choices)]
    }
    pickerInput(
      inputId = "VJ_deselected",
      label = "Exclude the following VJ combinations:",
      choices = choices,
      selected= tmp_selection,
      options = pickerOptions(
        actionsBox = TRUE,
        liveSearch = TRUE,
        maxOptions=(length(choices))-1,
        size = 10,
        selectedTextFormat = "count > 3",
        selectAllText="Select All but 1"
      ),
      multiple = TRUE
    )

  })

  observeEvent(ignoreInit = TRUE,
               Big_mem_values$Short_DF,
               {
                 if(!is.null(Big_mem_values$Short_DF) && any(startsWith(colnames(Big_mem_values$Short_DF),"Clone_"))) {
                   updateSelectizeInput(session, "my_clone_def",
                                        choices = colnames(Big_mem_values$Short_DF)[which(startsWith(colnames(Big_mem_values$Short_DF),"Clone_"))]
                   )
                 }


               })



  output$Group_comparison=renderUI({

    labels=sort(unique(Big_mem_values$Short_DF[[input$groups_selected]][which(Big_mem_values$Short_DF$Patient_Sample %in% (input$samples_selected))]))
    bucket_list(
      header = "Include at least one element on each group",
      orientation="horizontal",
      group_name = "group_classification",
      add_rank_list(
        text = "Group A",
        labels = labels,
        input_id="group_A"
      ),
      add_rank_list(
        text = "Group B",
        labels = NULL,
        input_id="group_B"
      ),

      add_rank_list(
        text = "Ignore",
        labels = NULL,
        input_id="group_C"
      )
    )
  })


  observeEvent(ignoreInit = TRUE, list(
    input$Sample_selection,
    input$use_what,
    input$use_productive_or_not
  ) , {

    rows=c( )
    if("Repertoire" %in% input$use_what) {
      rows=c(rows,which(Big_mem_values$Short_DF$Sequence_type=="Repertoire"))
    }
    if("Reconstructed germline" %in% input$use_what) {
      rows=c(rows,which(Big_mem_values$Short_DF$Sequence_type!="Reconstructed germline"))
    }
    rows=intersect(rows, which(Big_mem_values$Short_DF$Patient_Sample %in% (input$samples_selected)))
    print(rows)
    updateSliderInput(inputId="Rmut_filter", min=min(Big_mem_values$Big_DF[rows,which(Big_mem_values$Header == "AA_Whole_Replacement_muts")]), max= max(Big_mem_values$Big_DF[rows,which(Big_mem_values$Header == "AA_Whole_Replacement_muts")]))

  })


  shared_VJ_list <- reactiveVal(NULL)

  observeEvent(ignoreInit = TRUE, list(
    input$samples_selected,
    input$use_what,
    input$groups_selected,
    input$VJ_deselected,
    input$use_sharedVDJ,
    input$VDJ_normalized_per_size,
    input$my_clone_def,
    input$Rmut_filter,
    input$VDJ_normalized_per_sample,
    input$VDJ_maximize_clones
  ) , {
    rows=c()

    if(input$use_sharedVDJ){
      if("Repertoire" %in% input$use_what) {
        rows=c(rows,which(Big_mem_values$Short_DF$Sequence_type=="Repertoire"))
      }
      if("Reconstructed germline" %in% input$use_what) {
        rows=c(rows,which(Big_mem_values$Short_DF$Sequence_type=="Reconstructed germline"))
      }
      print(length(rows))
      rows=sort(rows)
      intersect_mut=intersect(which(Big_mem_values$Big_DF[,which(Big_mem_values$Header == "AA_Whole_Replacement_muts")] >= input$Rmut_filter[1]),
                              which(Big_mem_values$Big_DF[,which(Big_mem_values$Header == "AA_Whole_Replacement_muts")] <= input$Rmut_filter[2]))

      intersect_groupA_samples=intersect(which(Big_mem_values$Short_DF[[input$groups_selected]] %in% input$group_A),
                                         which(Big_mem_values$Short_DF$Patient_Sample %in% (input$samples_selected)))
      intersect_all_A=intersect(intersect(intersect_groupA_samples, rows),
                                intersect_mut)
      int_group_A=Big_mem_values$Short_DF$V_and_J[intersect_all_A]
      table_int_group_A=table(int_group_A)
      intersect_groupB_samples=intersect(which(Big_mem_values$Short_DF[[input$groups_selected]] %in% input$group_B),
                                         which(Big_mem_values$Short_DF$Patient_Sample %in% (input$samples_selected)))
      intersect_all_B=intersect(intersect(intersect_groupB_samples, rows),
                                intersect_mut)
      int_group_B=Big_mem_values$Short_DF$V_and_J[intersect_all_B]
      table_int_group_B=table(int_group_B)

      labels=sort(unique(intersect(int_group_A,
                                   int_group_B)))
      # labels=intersect(labels, input$VJ_deselected)
      labels=labels[which(labels %!in% input$VJ_deselected)]
      samples_A=(unique(Big_mem_values$Short_DF$Patient_Sample[intersect_groupA_samples]))
      samples_B=(unique(Big_mem_values$Short_DF$Patient_Sample[intersect_groupB_samples]))
      if(input$VDJ_normalized_per_size && input$VDJ_normalized_per_sample) {

        min_sizes_A=sapply(labels,
                           function(z)  min(sapply(samples_A,
                                                   function(y) length(which(Big_mem_values$Short_DF$V_and_J[intersect(intersect_all_A, which(Big_mem_values$Short_DF$Patient_Sample %in% y )  )] ==z  )) ))
        )
        min_sizes_B=sapply(labels,
                           function(z)  min(sapply(samples_B,
                                                   function(y) length(which(Big_mem_values$Short_DF$V_and_J[intersect(intersect_all_B, which(Big_mem_values$Short_DF$Patient_Sample %in% y )  )] ==z  )) ))
        )
        min_total_sizes= sapply(1:length(labels),
                                function(z) min(min_sizes_A[z]*length(samples_A),min_sizes_B[z]*length(samples_B)) )
        min_total_sizes_per_sample_groupA=floor(min_total_sizes/length(samples_A))
        min_total_sizes_per_sample_groupB=floor(min_total_sizes/length(samples_B))

        labels=labels[which(min_total_sizes !=0)]
        min_total_sizes_per_sample_groupA=min_total_sizes_per_sample_groupA[which(min_total_sizes !=0)]
        min_total_sizes_per_sample_groupB=min_total_sizes_per_sample_groupB[which(min_total_sizes !=0)]
        min_total_sizes=min_total_sizes[which(min_total_sizes !=0)]
        min_sizes_B=min_sizes_B[which(min_total_sizes !=0)]
        min_sizes_A=min_sizes_A[which(min_total_sizes !=0)]
        if( !is.null(input$my_clone_def) && input$VDJ_maximize_clones){

          names(labels)= sapply(c(1:length(labels)), function(z) paste(labels[z],  paste(min_total_sizes[z], paste(" (",
                                                                                                                   paste(sapply(c(1:length(samples_A)), function(y)  min(min_total_sizes_per_sample_groupA[z], length(unique(Big_mem_values$Short_DF[[input$my_clone_def]][intersect(intersect(intersect_all_A, which(Big_mem_values$Short_DF$Patient_Sample %in% samples_A[y] )), which(Big_mem_values$Short_DF$V_and_J ==labels[z]) )]))) )  , collapse=", "),
                                                                                                                   " clones on each sample)",
                                                                                                                   sep=""),
                                                                                         " & group B: ",
                                                                                         min_total_sizes[z],
                                                                                         paste(" (",
                                                                                               paste(sapply(c(1:length(samples_B)), function(y)  min(min_total_sizes_per_sample_groupB[z], length(unique(Big_mem_values$Short_DF[[input$my_clone_def]][intersect(intersect(intersect_all_B, which(Big_mem_values$Short_DF$Patient_Sample %in% samples_B[y] )), which(Big_mem_values$Short_DF$V_and_J ==labels[z]) )]))) ) , collapse=", "),
                                                                                               " clones on each sample)",
                                                                                               sep=""),
                                                                                         sep=""),
                                                                       sep=" - Sequence counts: group A: ") )

        } else {
          names(labels)= paste(labels, " - Counts: group A - ", min_total_sizes," (",min_total_sizes_per_sample_groupA ," per sample) and ",min_total_sizes-min_total_sizes_per_sample_groupA*length(samples_A) ,
                               " random & group B - ", min_total_sizes," (",min_total_sizes_per_sample_groupB ," per sample) and ",min_total_sizes-min_total_sizes_per_sample_groupB*length(samples_B), " random", sep="")
        }

      } else if(input$VDJ_normalized_per_size) {

        if( !is.null(input$my_clone_def) && input$VDJ_maximize_clones){

          names(labels)= sapply(labels, function(z) paste(z,  paste(min(table_int_group_A[which(names(table_int_group_A)==z)], table_int_group_B[which(names(table_int_group_B)==z)]), paste(" (",
                                                                                                                                                                                             paste(sapply(samples_A, function(y)  min(min(table_int_group_A[which(names(table_int_group_A)==z)], table_int_group_B[which(names(table_int_group_B)==z)]), length(unique(Big_mem_values$Short_DF[[input$my_clone_def]][intersect(intersect(intersect_all_A, which(Big_mem_values$Short_DF$Patient_Sample %in% y )), which(Big_mem_values$Short_DF$V_and_J ==z) )]))) )  , collapse=", "),
                                                                                                                                                                                             " clones on each sample)",
                                                                                                                                                                                             sep=""),
                                                                    " & group B: ",
                                                                    min(table_int_group_A[which(names(table_int_group_A)==z)], table_int_group_B[which(names(table_int_group_B)==z)]),
                                                                    paste(" (",
                                                                          paste(sapply(samples_B, function(y)  min(min(table_int_group_A[which(names(table_int_group_A)==z)], table_int_group_B[which(names(table_int_group_B)==z)]), length(unique(Big_mem_values$Short_DF[[input$my_clone_def]][intersect(intersect(intersect_all_B, which(Big_mem_values$Short_DF$Patient_Sample %in% y )), which(Big_mem_values$Short_DF$V_and_J ==z) )]))) ) , collapse=", "),
                                                                          " clones on each sample)",
                                                                          sep=""),
                                                                    sep=""),
                                                          sep=" - Sequence counts: group A: ") )

        } else {
          names(labels)= sapply(labels, function(z) paste(z,  paste(min(table_int_group_A[which(names(table_int_group_A)==z)], table_int_group_B[which(names(table_int_group_B)==z)]), min(table_int_group_A[which(names(table_int_group_A)==z)], table_int_group_B[which(names(table_int_group_B)==z)]), sep=" & group B - "), sep=" - Counts: group A -") )

        }



      } else {

        if( !is.null(input$my_clone_def) && input$VDJ_maximize_clones) {
          names(labels)= sapply(labels, function(z) paste(z,  paste(table_int_group_A[which(names(table_int_group_A)==z)], paste(" (",
                                                                                                                                 paste(sapply(samples_A, function(y)  length(unique(Big_mem_values$Short_DF[[input$my_clone_def]][intersect(intersect(intersect_all_A, which(Big_mem_values$Short_DF$Patient_Sample %in% y )), which(Big_mem_values$Short_DF$V_and_J ==z) )])) )  , collapse=", "),
                                                                                                                                 " clones on each sample)",
                                                                                                                                 sep=""),
                                                                    " & group B - ",
                                                                    table_int_group_B[which(names(table_int_group_B)==z)],
                                                                    paste(" (",
                                                                          paste(sapply(samples_B, function(y)  length(unique(Big_mem_values$Short_DF[[input$my_clone_def]][intersect(intersect(intersect_all_B, which(Big_mem_values$Short_DF$Patient_Sample %in% y )), which(Big_mem_values$Short_DF$V_and_J ==z) )])) ) , collapse=", "),
                                                                          " clones on each sample)",
                                                                          sep=""),
                                                                    sep=""),
                                                          sep=" - Sequence counts: group A -") )

        } else{
          names(labels)= sapply(labels, function(z) paste(z,  paste(table_int_group_A[which(names(table_int_group_A)==z)], table_int_group_B[which(names(table_int_group_B)==z)], sep=" & group B - "), sep=" - Counts: group A -") )

        }
      }
      shared_VJ_list(labels)
    } else {
      shared_VJ_list(NULL)
    }

  })

  output$VDJ_subsetting=renderUI({
    print("VDJ_subsetting")



    labels=shared_VJ_list()
    # if(length(labels)==0) {
    #   labels="No shared VJ"
    # }
    # print(labels[1])
    # print(intersect(which(Big_mem_values$Short_DF[[input$groups_selected]] %in% input$group_A),
    #                 which(Big_mem_values$Short_DF$Patient_Sample %in% (input$samples_selected))))
    # print(Big_mem_values$Short_DF$V_and_J[1:3])
    # print(rows[1:10])
    # print(Big_mem_values$Short_DF$V_and_J[intersect(intersect(which(Big_mem_values$Short_DF[[input$groups_selected]] %in% input$group_A),
    #                                                           which(Big_mem_values$Short_DF$Patient_Sample %in% (input$samples_selected))), rows)])
    # print(Big_mem_values$Short_DF$V_and_J[intersect(intersect(which(Big_mem_values$Short_DF[[input$groups_selected]] %in% input$group_B),
    #                                                           which(Big_mem_values$Short_DF$Patient_Sample %in% (input$samples_selected))), rows)])
    # print(labels)
    # bucket_list(
    #   header = "Select shared V-J combinations between groups",
    #   orientation="horizontal",
    #   group_name = "V_J_classification",
    #   add_rank_list(
    #     text = "To include",
    #     labels = NULL,
    #     input_id="VJ_included"
    #   ),
    #   add_rank_list(
    #     text = "To exclude",
    #     labels = labels,
    #     input_id="VJ_excluded"
    #   )
    # )

    pickerInput(
      inputId = "VJ_included",
      label = "To include only the following VJ combinations:",
      choices = labels,
      options = pickerOptions(
        actionsBox = TRUE,
        liveSearch = TRUE,
        size = 10,
        selectedTextFormat = "count > 3"
      ),
      multiple = TRUE
    )

  })




  individual_variables_current_list <- reactiveVal(NULL)

  # now store your current selection in the reactive value
  observeEvent(ignoreInit = TRUE, list(
    input$samples_selected,
    input$my_vars,
    input$use_what,
    input$use_productive_or_not,
    input$my_regions,
    input$my_var_elements,
    input$my_vartypes
  ) , {

    choices=list()

    if(!is.null(input$my_vars)){
      columns=c(1:length(Big_mem_values$Header))
      columns=columns[which(columns %in% which((sapply(Big_mem_values$Header, function(x) strsplit(x, split="_")[[1]][2])) %in% input$my_regions) )]
      columns=columns[which(columns %in% which(grepl(paste0(paste("^", input$my_var_elements, sep=""), collapse ="|"),Big_mem_values$Header)))]
      if("Germline diff" %!in% input$my_vartypes) {
        columns=columns[which(columns %!in% which(endsWith(Big_mem_values$Header, "_Diff")))]
      }
      if("Baseline" %!in% input$my_vartypes) {
        columns=columns[which(columns %in% which(endsWith(Big_mem_values$Header, "_Diff")))]
      }

      merged_header_splitted_at_3=sapply( Big_mem_values$Header, function(x) strsplit(x, split="_")[[1]][3])
      merged_header_splitted_at_4=sapply( Big_mem_values$Header, function(x) strsplit(x, split="_")[[1]][4])
      merged_header_splitted_at_3_and_4=sapply( Big_mem_values$Header, function(x) paste(strsplit(x, split="_")[[1]][3], strsplit(x, split="_")[[1]][4], sep="_"))


      if ("Length" %in% input$my_vars) {
        tmp_index_var=c(columns[which(columns %in% which(endsWith(Big_mem_values$Header, "length")))],
                        columns[which(columns %in% which(endsWith(Big_mem_values$Header, "length_Diff")))])
        choices[["Length"]] = as.list(sort(Big_mem_values$Header[tmp_index_var]))
      }

      if ("Composition" %in% input$my_vars) {
        nucleotides=c("A","G","T","C")
        peptides=Peptides::aaList()

        tmp_index_var=c(columns[which(columns %in% which(merged_header_splitted_at_3_and_4 %in%
                                                           unique(c(paste(c(nucleotides, peptides), "count", sep="_"), paste(c(nucleotides, peptides), "norm", sep="_") )) ))],
                        columns[which(columns %in% which(merged_header_splitted_at_4 %in% c("counts")   ))])
        choices[["Composition"]] = as.list(sort(Big_mem_values$Header[tmp_index_var]))
      }

      if ("Hot/Cold motifs" %in% input$my_vars) {
        tmp_index_var=c(columns[which(columns %in% which(merged_header_splitted_at_3 %in% c("hot","cold","potential")   ))])
        choices[["Hot/Cold motifs"]] = as.list(sort(Big_mem_values$Header[tmp_index_var]))
      }

      if ("Substitutions" %in% input$my_vars) {
        tmp_index_var=c(columns[which(columns %in% which(merged_header_splitted_at_3_and_4 %in%
                                                           c("Sub", "Sub_prc", "SIDT_sum", "SID_sum")   ))],
                        columns[which(columns %in% which(merged_header_splitted_at_3 %in%
                                                           c("Sub")   ))])
        choices[["Substitutions"]] = as.list(sort(Big_mem_values$Header[tmp_index_var]))
      }

      if ("Insertions" %in% input$my_vars) {
        tmp_index_var=c(columns[which(columns %in% which(merged_header_splitted_at_3_and_4 %in%
                                                           c("Ins", "Ins_prc", "SIDT_sum", "SID_sum")   ))],
                        columns[which(columns %in% which(merged_header_splitted_at_3 %in%
                                                           c("Ins")   ))])
        choices[["Insertions"]] = as.list(sort(Big_mem_values$Header[tmp_index_var]))
      }
      if ("Deletions" %in% input$my_vars) {
        tmp_index_var=c(columns[which(columns %in% which(merged_header_splitted_at_3_and_4 %in%
                                                           c("Del", "Del_prc", "SIDT_sum", "SID_sum")   ))],
                        columns[which(columns %in% which(merged_header_splitted_at_3 %in%
                                                           c("Del")   ))])
        choices[["Deletions"]] = as.list(sort(Big_mem_values$Header[tmp_index_var]))
      }
      if ("Translocations" %in% input$my_vars) {

        tmp_index_var=c(columns[which(columns %in% which(merged_header_splitted_at_3_and_4 %in%
                                                           c("Trasl", "Trasl_prc", "SIDT_sum")   ))],
                        columns[which(columns %in% which(merged_header_splitted_at_3 %in%
                                                           c("Trasl")   ))])
        choices[["Translocations"]] = as.list(sort(Big_mem_values$Header[tmp_index_var]))
      }

      if("Leveshtein distance" %in% input$my_vars) {
        tmp_index_var=c(columns[which(columns %in% which(merged_header_splitted_at_3 %in% c("lv")   ))])
        choices[["Leveshtein distance"]] = as.list(sort(Big_mem_values$Header[tmp_index_var]))
      }

      if ("Transitions and transversions" %in% input$my_vars) {
        tmp_index_var=c(columns[which(columns %in% which(merged_header_splitted_at_3 %in% c("Transitions","Transversions")   ))],
                        columns[which(columns %in% which(merged_header_splitted_at_3_and_4 %in%
                                                           c("Ratio_Transitions-Transversions")   ))])
        choices[["Transitions and transversions"]] = as.list(sort(Big_mem_values$Header[tmp_index_var]))
      }
      if ("Replacement and silent mutations" %in% input$my_vars) {
        tmp_index_var=c(columns[which(columns %in% which(merged_header_splitted_at_3 %in% c("Replacement","Silent")   ))],
                        columns[which(columns %in% which(merged_header_splitted_at_3_and_4 %in%
                                                           c("Ratio_Silent-Replacement")   ))])
        choices[["Replacement and silent mutations"]] = as.list(sort(Big_mem_values$Header[tmp_index_var]))
      }

      if ("Mutations from X to Y" %in% input$my_vars) {
        # columns=columns[which(columns %!in% which(sapply(merged_header, function(x) paste(strsplit(x, split="_")[[1]][3], strsplit(x, split="_")[[1]][4], sep="_")) %in%
        #                                             c("Ratio_Silent")   ))]
        tmp_index_var=c(columns[which(columns %in% which(merged_header_splitted_at_4 %in% c("to")   ))])
        choices[["Mutations from X to Y"]] = as.list(sort(Big_mem_values$Header[tmp_index_var]))
      }

      if ("NGly sites" %in% input$my_vars) {
        tmp_index_var=c(columns[which(columns %in% which(merged_header_splitted_at_3 %in% c("NGly")   ))])
        choices[["NGly sites"]] = as.list(sort(Big_mem_values$Header[tmp_index_var]))
      }
      if ("Peptide features" %in% input$my_vars) {
        tmp_index_var=c(columns[which(columns %in% which(merged_header_splitted_at_3 %in% c("Peptides", "alkzm",
                                                                                            "Small","Tiny", "Aliphatic",
                                                                                            "Charged","Polar", "Basic","NonPolar",
                                                                                            "Aromatic", "Acidic")   ))])
        choices[["Peptide features"]] = as.list(sort(Big_mem_values$Header[tmp_index_var]))
      }


    }

    individual_variables_current_list(choices)
  })




  output$individual_variables=renderUI({
    choices = c("")
    if(!is.null(Exploration_values$rows) && length(Exploration_values$rows)>0 && !is.null(Exploration_values$columns) && length(Exploration_values$columns)>0) {
      choices=individual_variables_current_list()

    }
    # selectizeInput("exclude_variables", "Excluding individual variables (start writing the category):",
    #                choices=choices, multiple = TRUE,
    #                options = list(
    #                  # onInitialize = I(onInitialize),
    #                  onDropdownOpen = I(onDropdownOpen)
    #                )
    # )

    virtualSelectInput(
      inputId = "exclude_variables",
      label = "Remove individual variables:",
      choices = choices,
      showValueAsTags = TRUE,
      search = TRUE,
      multiple = TRUE
    )
  })


  output$RAM_S_Memory_Box <- renderInfoBox({
    infoBox(
      title="Selection dataset",value=tags$p(
        if (as.numeric(unname(strsplit(as.character(print(get_ram(), unit_system = "metric")), split=" ")[[1]][1])) >   round(length(Selection_values$rows)*length(Selection_values$columns)*8/2^{20}/1024, 2)) { "Fits in RAM, more ML methods are available. If both datasets fit in RAM, UMAP and PCA biplot will be available" } else { "Does not fit in RAM"}, style = "font-size: 50%;"),
      icon = icon("info-sign", lib = "glyphicon"),
      color = if (as.numeric(unname(strsplit(as.character(print(get_ram(), unit_system = "metric")), split=" ")[[1]][1])) >   round(length(Selection_values$rows)*length(Selection_values$columns)*8/2^{20}/1024, 2)) { "aqua" } else { "red" }
    )
  })
  output$RAM_E_Memory_Box <- renderInfoBox({
    infoBox(
      title="Exploration dataset",value=tags$p(
        if (as.numeric(unname(strsplit(as.character(print(get_ram(), unit_system = "metric")), split=" ")[[1]][1])) >   round(length(Exploration_values$rows)*length(Exploration_values$columns)*8/2^{20}/1024, 2)) { "Fits in RAM. If both datasets fit in RAM, UMAP and PCA biplot will be available." } else { "Does not fit in RAM"}, style = "font-size: 50%;"),
      icon = icon("info-sign", lib = "glyphicon"),
      color = if (as.numeric(unname(strsplit(as.character(print(get_ram(), unit_system = "metric")), split=" ")[[1]][1])) >   round(length(Exploration_values$rows)*length(Exploration_values$columns)*8/2^{20}/1024, 2)) { "aqua" } else { "red" }
    )
  })



  output$sunburst_selection_Whole <- renderUI({
    if(!is.null(Selection_values$rows) && length(Selection_values$rows)>0 && !is.null(Selection_values$columns) && length(Selection_values$columns)>0) {
      tmp_sunburst=unname(unlist(Big_mem_values$Header[Selection_values$columns]))
      # regions=c("Whole","FWR1","CDR1","FWR2","CDR2","FWR3","CDR3","FWR4")
      region="Whole"
      sunburst_df=show_selected_features(tmp_sunburst, region=region)
      # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel= "Variables analyzed per each region in the selection set") )

      if(dim(sunburst_df)[1]==0) {
        HTML(paste(c("There is no", region, "sequence variables analyzed in the selection set"), collapse=" "))
      } else {
        sund2b(sunburst_df,showLabels =F, rootLabel= paste(c(region, "sequence variables analyzed in the selection set"), collapse=" "))
      }

    }


  })
  output$sunburst_selection_FWR1 <- renderUI({
    if(!is.null(Selection_values$rows) && length(Selection_values$rows)>0 && !is.null(Selection_values$columns) && length(Selection_values$columns)>0) {
      tmp_sunburst=unname(unlist(Big_mem_values$Header[Selection_values$columns]))
      # regions=c("Whole","FWR1","CDR1","FWR2","CDR2","FWR3","CDR3","FWR4")
      region="FWR1"
      sunburst_df=show_selected_features(tmp_sunburst, region=region)
      # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel= "Variables analyzed per each region in the selection set") )
      if(dim(sunburst_df)[1]==0) {
        HTML(paste(c("There is no", region, "sequence variables analyzed in the selection set"), collapse=" "))
      } else {
        sund2b(sunburst_df,showLabels =F, rootLabel= paste(c(region, "sequence variables analyzed in the selection set"), collapse=" "))
      }
    }
  })
  output$sunburst_selection_CDR1 <- renderUI({
    if(!is.null(Selection_values$rows) && length(Selection_values$rows)>0 && !is.null(Selection_values$columns) && length(Selection_values$columns)>0) {
      tmp_sunburst=unname(unlist(Big_mem_values$Header[Selection_values$columns]))
      # regions=c("Whole","FWR1","CDR1","FWR2","CDR2","FWR3","CDR3","FWR4")
      region="CDR1"
      sunburst_df=show_selected_features(tmp_sunburst, region=region)
      # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel= "Variables analyzed per each region in the selection set") )
      if(dim(sunburst_df)[1]==0) {
        HTML(paste(c("There is no", region, "sequence variables analyzed in the selection set"), collapse=" "))
      } else {
        sund2b(sunburst_df,showLabels =F, rootLabel= paste(c(region, "sequence variables analyzed in the selection set"), collapse=" "))
      }
    }
  })
  output$sunburst_selection_FWR2 <- renderUI({
    if(!is.null(Selection_values$rows) && length(Selection_values$rows)>0 && !is.null(Selection_values$columns) && length(Selection_values$columns)>0) {
      tmp_sunburst=unname(unlist(Big_mem_values$Header[Selection_values$columns]))
      # regions=c("Whole","FWR1","CDR1","FWR2","CDR2","FWR3","CDR3","FWR4")
      region="FWR2"
      sunburst_df=show_selected_features(tmp_sunburst, region=region)
      # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel= "Variables analyzed per each region in the selection set") )
      if(dim(sunburst_df)[1]==0) {
        HTML(paste(c("There is no", region, "sequence variables analyzed in the selection set"), collapse=" "))
      } else {
        sund2b(sunburst_df,showLabels =F, rootLabel= paste(c(region, "sequence variables analyzed in the selection set"), collapse=" "))
      }
    }
  })
  output$sunburst_selection_CDR2 <- renderUI({
    if(!is.null(Selection_values$rows) && length(Selection_values$rows)>0 && !is.null(Selection_values$columns) && length(Selection_values$columns)>0) {
      tmp_sunburst=unname(unlist(Big_mem_values$Header[Selection_values$columns]))
      # regions=c("Whole","FWR1","CDR1","FWR2","CDR2","FWR3","CDR3","FWR4")
      region="CDR2"
      sunburst_df=show_selected_features(tmp_sunburst, region=region)
      # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel= "Variables analyzed per each region in the selection set") )
      if(dim(sunburst_df)[1]==0) {
        HTML(paste(c("There is no", region, "sequence variables analyzed in the selection set"), collapse=" "))
      } else {
        sund2b(sunburst_df,showLabels =F, rootLabel= paste(c(region, "sequence variables analyzed in the selection set"), collapse=" "))
      }
    }
  })
  output$sunburst_selection_FWR3 <- renderUI({
    if(!is.null(Selection_values$rows) && length(Selection_values$rows)>0 && !is.null(Selection_values$columns) && length(Selection_values$columns)>0) {
      tmp_sunburst=unname(unlist(Big_mem_values$Header[Selection_values$columns]))
      # regions=c("Whole","FWR1","CDR1","FWR2","CDR2","FWR3","CDR3","FWR4")
      region="FWR3"
      sunburst_df=show_selected_features(tmp_sunburst, region=region)
      # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel= "Variables analyzed per each region in the selection set") )
      if(dim(sunburst_df)[1]==0) {
        HTML(paste(c("There is no", region, "sequence variables analyzed in the selection set"), collapse=" "))
      } else {
        sund2b(sunburst_df,showLabels =F, rootLabel= paste(c(region, "sequence variables analyzed in the selection set"), collapse=" "))
      }
    }
  })
  output$sunburst_selection_CDR3 <- renderUI({
    if(!is.null(Selection_values$rows) && length(Selection_values$rows)>0 && !is.null(Selection_values$columns) && length(Selection_values$columns)>0) {
      tmp_sunburst=unname(unlist(Big_mem_values$Header[Selection_values$columns]))
      # regions=c("Whole","FWR1","CDR1","FWR2","CDR2","FWR3","CDR3","FWR4")
      region="CDR3"
      sunburst_df=show_selected_features(tmp_sunburst, region=region)
      # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel= "Variables analyzed per each region in the selection set") )
      if(dim(sunburst_df)[1]==0) {
        HTML(paste(c("There is no", region, "sequence variables analyzed in the selection set"), collapse=" "))
      } else {
        sund2b(sunburst_df,showLabels =F, rootLabel= paste(c(region, "sequence variables analyzed in the selection set"), collapse=" "))
      }
    }
  })
  output$sunburst_selection_FWR4 <- renderUI({
    if(!is.null(Selection_values$rows) && length(Selection_values$rows)>0 && !is.null(Selection_values$columns) && length(Selection_values$columns)>0) {
      tmp_sunburst=unname(unlist(Big_mem_values$Header[Selection_values$columns]))
      # regions=c("Whole","FWR1","CDR1","FWR2","CDR2","FWR3","CDR3","FWR4")
      region="FWR4"
      sunburst_df=show_selected_features(tmp_sunburst, region=region)
      # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel= "Variables analyzed per each region in the selection set") )
      if(dim(sunburst_df)[1]==0) {
        HTML(paste(c("There is no", region, "sequence variables analyzed in the selection set"), collapse=" "))
      } else {
        sund2b(sunburst_df,showLabels =F, rootLabel= paste(c(region, "sequence variables analyzed in the selection set"), collapse=" "))
      }
    }
  })


  output$sunburst_exploration_Whole <- renderUI({
    if(!is.null(Exploration_values$rows) && length(Exploration_values$rows)>0 && !is.null(Exploration_values$columns) && length(Exploration_values$columns)>0) {
      tmp_sunburst=unname(unlist(Big_mem_values$Header[Exploration_values$columns]))
      # regions=c("Whole","FWR1","CDR1","FWR2","CDR2","FWR3","CDR3","FWR4")
      region="Whole"
      sunburst_df=show_selected_features(tmp_sunburst, region=region)
      # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel= "Variables analyzed per each region in the exploration set") )

      if(dim(sunburst_df)[1]==0) {
        HTML(paste(c("There is no", region, "sequence variables analyzed in the exploration set"), collapse=" "))
      } else {
        sund2b(sunburst_df,showLabels =F, rootLabel= paste(c(region, "sequence variables analyzed in the exploration set"), collapse=" "))
      }

    }


  })
  output$sunburst_exploration_FWR1 <- renderUI({
    if(!is.null(Exploration_values$rows) && length(Exploration_values$rows)>0 && !is.null(Exploration_values$columns) && length(Exploration_values$columns)>0) {
      tmp_sunburst=unname(unlist(Big_mem_values$Header[Exploration_values$columns]))
      # regions=c("Whole","FWR1","CDR1","FWR2","CDR2","FWR3","CDR3","FWR4")
      region="FWR1"
      sunburst_df=show_selected_features(tmp_sunburst, region=region)
      # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel= "Variables analyzed per each region in the exploration set") )
      if(dim(sunburst_df)[1]==0) {
        HTML(paste(c("There is no", region, "sequence variables analyzed in the exploration set"), collapse=" "))
      } else {
        sund2b(sunburst_df,showLabels =F, rootLabel= paste(c(region, "sequence variables analyzed in the exploration set"), collapse=" "))
      }
    }
  })
  output$sunburst_exploration_CDR1 <- renderUI({
    if(!is.null(Exploration_values$rows) && length(Exploration_values$rows)>0 && !is.null(Exploration_values$columns) && length(Exploration_values$columns)>0) {
      tmp_sunburst=unname(unlist(Big_mem_values$Header[Exploration_values$columns]))
      # regions=c("Whole","FWR1","CDR1","FWR2","CDR2","FWR3","CDR3","FWR4")
      region="CDR1"
      sunburst_df=show_selected_features(tmp_sunburst, region=region)
      # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel= "Variables analyzed per each region in the exploration set") )
      if(dim(sunburst_df)[1]==0) {
        HTML(paste(c("There is no", region, "sequence variables analyzed in the exploration set"), collapse=" "))
      } else {
        sund2b(sunburst_df,showLabels =F, rootLabel= paste(c(region, "sequence variables analyzed in the exploration set"), collapse=" "))
      }
    }
  })
  output$sunburst_exploration_FWR2 <- renderUI({
    if(!is.null(Exploration_values$rows) && length(Exploration_values$rows)>0 && !is.null(Exploration_values$columns) && length(Exploration_values$columns)>0) {
      tmp_sunburst=unname(unlist(Big_mem_values$Header[Exploration_values$columns]))
      # regions=c("Whole","FWR1","CDR1","FWR2","CDR2","FWR3","CDR3","FWR4")
      region="FWR2"
      sunburst_df=show_selected_features(tmp_sunburst, region=region)
      # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel= "Variables analyzed per each region in the exploration set") )
      if(dim(sunburst_df)[1]==0) {
        HTML(paste(c("There is no", region, "sequence variables analyzed in the exploration set"), collapse=" "))
      } else {
        sund2b(sunburst_df,showLabels =F, rootLabel= paste(c(region, "sequence variables analyzed in the exploration set"), collapse=" "))
      }
    }
  })
  output$sunburst_exploration_CDR2 <- renderUI({
    if(!is.null(Exploration_values$rows) && length(Exploration_values$rows)>0 && !is.null(Exploration_values$columns) && length(Exploration_values$columns)>0) {
      tmp_sunburst=unname(unlist(Big_mem_values$Header[Exploration_values$columns]))
      # regions=c("Whole","FWR1","CDR1","FWR2","CDR2","FWR3","CDR3","FWR4")
      region="CDR2"
      sunburst_df=show_selected_features(tmp_sunburst, region=region)
      # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel= "Variables analyzed per each region in the exploration set") )
      if(dim(sunburst_df)[1]==0) {
        HTML(paste(c("There is no", region, "sequence variables analyzed in the exploration set"), collapse=" "))
      } else {
        sund2b(sunburst_df,showLabels =F, rootLabel= paste(c(region, "sequence variables analyzed in the exploration set"), collapse=" "))
      }
    }
  })
  output$sunburst_exploration_FWR3 <- renderUI({
    if(!is.null(Exploration_values$rows) && length(Exploration_values$rows)>0 && !is.null(Exploration_values$columns) && length(Exploration_values$columns)>0) {
      tmp_sunburst=unname(unlist(Big_mem_values$Header[Exploration_values$columns]))
      # regions=c("Whole","FWR1","CDR1","FWR2","CDR2","FWR3","CDR3","FWR4")
      region="FWR3"
      sunburst_df=show_selected_features(tmp_sunburst, region=region)
      # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel= "Variables analyzed per each region in the exploration set") )
      if(dim(sunburst_df)[1]==0) {
        HTML(paste(c("There is no", region, "sequence variables analyzed in the exploration set"), collapse=" "))
      } else {
        sund2b(sunburst_df,showLabels =F, rootLabel= paste(c(region, "sequence variables analyzed in the exploration set"), collapse=" "))
      }
    }
  })
  output$sunburst_exploration_CDR3 <- renderUI({
    if(!is.null(Exploration_values$rows) && length(Exploration_values$rows)>0 && !is.null(Exploration_values$columns) && length(Exploration_values$columns)>0) {
      tmp_sunburst=unname(unlist(Big_mem_values$Header[Exploration_values$columns]))
      # regions=c("Whole","FWR1","CDR1","FWR2","CDR2","FWR3","CDR3","FWR4")
      region="CDR3"
      sunburst_df=show_selected_features(tmp_sunburst, region=region)
      # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel= "Variables analyzed per each region in the exploration set") )
      if(dim(sunburst_df)[1]==0) {
        HTML(paste(c("There is no", region, "sequence variables analyzed in the exploration set"), collapse=" "))
      } else {
        sund2b(sunburst_df,showLabels =F, rootLabel= paste(c(region, "sequence variables analyzed in the exploration set"), collapse=" "))
      }
    }
  })
  output$sunburst_exploration_FWR4 <- renderUI({
    if(!is.null(Exploration_values$rows) && length(Exploration_values$rows)>0 && !is.null(Exploration_values$columns) && length(Exploration_values$columns)>0) {
      tmp_sunburst=unname(unlist(Big_mem_values$Header[Exploration_values$columns]))
      # regions=c("Whole","FWR1","CDR1","FWR2","CDR2","FWR3","CDR3","FWR4")
      region="FWR4"
      sunburst_df=show_selected_features(tmp_sunburst, region=region)
      # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel= "Variables analyzed per each region in the exploration set") )
      if(dim(sunburst_df)[1]==0) {
        HTML(paste(c("There is no", region, "sequence variables analyzed in the exploration set"), collapse=" "))
      } else {
        sund2b(sunburst_df,showLabels =F, rootLabel= paste(c(region, "sequence variables analyzed in the exploration set"), collapse=" "))
      }
    }
  })





  output$correlated_variables <- renderUI({
    if(!is.null(Selection_values$rows) && length(Selection_values$rows)>0 && !is.null(Selection_values$columns) && length(Selection_values$columns)>0) {



      print("MMMMMmmmmmmmmm correlating")
      correlated_BGM <- bigstatsr::big_cor(Big_mem_values$Big_DF, ind.col=Selection_values$columns, ind.row = Selection_values$rows)

      print("MMMMMmmmmmmmmm")

      index_correlated=c()
      pairs_correlated=list()
      react_table=data.frame(Region_of_the_variable=as.character(),	Variable=as.character(),	Region_of_the_correlated_variable=as.character(),
                             Correlated_variable=as.character(),	Pearson_Correlation=as.numeric(),	Raw_P_value=as.numeric())
      for (column in c(1:ncol(correlated_BGM))) {
        if (any(correlated_BGM[c(1:ncol(correlated_BGM))[which(c(1:ncol(correlated_BGM)) != column)],column] >=0.8)) {

          tmp_real_corr=c()
          tmp_pval_corr=c()
          for (column_corr in c(1:ncol(correlated_BGM))[which(correlated_BGM[c(1:ncol(correlated_BGM)),column] >=0.8)]) {


            if(column!= column_corr) {
              non_zero_rows=unique(c(which(Big_mem_values$Big_DF[Selection_values$rows,Selection_values$columns[column] ] != 0), which(Big_mem_values$Big_DF[Selection_values$rows,Selection_values$columns[column_corr]  ]!=0)))


              tmp_corr=cor(Big_mem_values$Big_DF[Selection_values$rows[non_zero_rows],Selection_values$columns[column]  ], Big_mem_values$Big_DF[Selection_values$rows[non_zero_rows],Selection_values$columns[column_corr] ])


              if(!is.na(tmp_corr) && tmp_corr >=0.8) {

                if (gsub("_count","",Big_mem_values$Header[Selection_values$columns[column]]) != gsub("_norm","",Big_mem_values$Header[Selection_values$columns[column_corr]]) ) {
                  if (gsub("_norm","",Big_mem_values$Header[Selection_values$columns[column]]) != gsub("_count","",Big_mem_values$Header[Selection_values$columns[column_corr]]) ) {
                    pval_corr=cor.test(Big_mem_values$Big_DF[Selection_values$rows,Selection_values$columns[column]  ], Big_mem_values$Big_DF[Selection_values$rows,Selection_values$columns[column_corr] ])$p.value

                    if(pval_corr <= 0.05) {
                      tmp_real_corr=c(tmp_real_corr, column_corr)
                      tmp_pval_corr=c(tmp_pval_corr, pval_corr)
                      react_table=rbind(react_table,
                                        data.frame(Region_of_the_variable=strsplit(Big_mem_values$Header[Selection_values$columns[column]], split="_")[[1]][2],	Variable=Big_mem_values$Header[Selection_values$columns[column]],
                                                   Region_of_the_correlated_variable=strsplit(Big_mem_values$Header[Selection_values$columns[column_corr]], split="_")[[1]][2],
                                                   Correlated_variable=Big_mem_values$Header[Selection_values$columns[column_corr]],
                                                   Pearson_Correlation=correlated_BGM[column_corr,column],	Raw_P_value=pval_corr))

                    }

                  }
                }

              }

            }

          }

          if(length(tmp_real_corr)>0) {
            pairs_correlated[[length(pairs_correlated)+1]]=tmp_real_corr
            names(pairs_correlated)=c(names(pairs_correlated)[1:(length(pairs_correlated)-1)], column)
            index_correlated=unique(c(index_correlated, column))

          }

        }
      }

      # for (keyvar in names(pairs_correlated)) {
      #   for (keyvar2 in pairs_correlated[[keyvar]]) {
      #     if (as.numeric(keyvar) %!in% pairs_correlated[[as.character(keyvar2)]]) {
      #       print("MMMMMMMMMMMMWIRD")
      #       print(keyvar)
      #       print(pairs_correlated[[(keyvar)]])
      #       print(keyvar2)
      #       print(pairs_correlated[[as.character(keyvar2)]])
      #     }
      #   }
      # }
      print(paste("Correlated elements: ", length(index_correlated), sep=""))
      print(paste("Correlated elements bis: ", length(unique(unlist(pairs_correlated))), sep=""))
      if(length(pairs_correlated)==0) {
        HTML("There is no correlation equal or greater than 0.8 between the sequence variables analyzed in the selection set")
      } else {

        ##Former heatmap
        # reference_reduced_matrix=correlated_BGM[index_correlated,index_correlated]
        #
        # reduced_matrix=as.data.frame(matrix(0, nrow = nrow(reference_reduced_matrix), ncol = ncol(reference_reduced_matrix)))
        #
        #
        # for(column in names(pairs_correlated)) {
        #   reduced_matrix[which(index_correlated == column), which(index_correlated %in% pairs_correlated[[as.character(column)]])] =reference_reduced_matrix[which(index_correlated == column), which(index_correlated %in% pairs_correlated[[as.character(column)]])]
        #   reduced_matrix[which(index_correlated %in% pairs_correlated[[as.character(column)]]), which(index_correlated == column)] = reference_reduced_matrix[which(index_correlated %in% pairs_correlated[[as.character(column)]]), which(index_correlated == column)]
        # }
        #
        # regions=Big_mem_values$Header[Selection_values$columns][index_correlated]
        #
        # for (region in c("Whole","FWR1","CDR1","FWR2","CDR2","FWR3","CDR3","FWR4")) {
        #   regions[grepl(region, regions)]=region
        # }
        # heatmaply(reduced_matrix,
        #           labRow=Big_mem_values$Header[Selection_values$columns[index_correlated]],
        #           labCol=Big_mem_values$Header[Selection_values$columns[index_correlated]],
        #           xlab = "Features",
        #           ylab = "Features",
        #           main = "Pearson correlation >=0.8 between non-zero values of different variable features",
        #           row_side_colors=data.frame(Regions=regions),
        #           col_side_colors=data.frame(Regions=regions),
        #           dendrogram="none",
        #           # limits = c(0.9, 1)
        # ) %>% layout(height=1400,width=1400)

        reactable(
          react_table,
          groupBy = c("Region_of_the_variable", "Variable"),
          columns = list(
            # Variable = colDef(aggregate = "unique"),
            Region_of_the_correlated_variable = colDef(aggregate = "unique")
          ),
          bordered = TRUE,
          elementId = "Correlation-table",
          searchable = TRUE
        )

      }
    }
  })

  output$numerical_covariables_plot <- renderUI({
    if(!is.null(Selection_values$rows) && length(Selection_values$rows)>0 && !is.null(Selection_values$columns) && length(Selection_values$columns)>0 && input$numerical_covariable_file !="") {

      plot_heatmap=F

      # gene_table=fread("/home/rgarcia/Descargas/Dani_genes/expression_matrix.txt")
      # gene_table2=readRDS("/home/rgarcia/Descargas/Dani_genes/expression_matrix.rds", refhook = NULL)
      # index_null=which(sapply(1:nrow(gene_table), function(z)  sum(gene_table[z, 2:ncol(gene_table)])) ==0)
      # gene_table=t(gene_table)
      # colnames(gene_table)=gene_table[1,]
      # gene_table=gene_table[2:nrow(gene_table),]
      # gene_table$Cell_ID=rownames(gene_table)
      #
      # write.table(gene_table, file="/home/rgarcia/Descargas/Dani_genes/expression_matrix_cor.txt", sep = "\t",quote = F,row.names = F,col.names = T)
      if ("Sequence_ID" %in% colnames(input$numerical_covariable_file)) {


        plot_heatmap=T
      } else if ("Cell_ID" %in% colnames(input$numerical_covariable_file )) {




        plot_heatmap=T
      }


      if (plot_heatmap) {

      }


    }
  })

  # 4. Features ####

  observeEvent(Selection_values$columns,{
    tmp_choices=Big_mem_values$Header[Selection_values$columns]
    # updateSelectInput(session, inputId = "plot_feature", choices=tmp_choices, selected=tmp_choices[1])
    updateSelectizeInput(session, inputId = "plot_feature", choices=tmp_choices, selected=tmp_choices[1], server=T,
                         options = list(maxOptions = length(tmp_choices)))
  })

  output$Group_selection_for_feature <- renderUI({
    if(is.null(Selection_values$rows) && length(Selection_values$rows)>0 && !is.null(Selection_values$columns) && length(Selection_values$columns)>0) {
      selectInput("plot_color_feature","Color by", choices = c("Sample",	"Patient",	"Group",	"Subgroup",	"Chain"), selected="Group")
    } else{
      sample_info<-sample_info_react$table
      choices_color=c(colnames(sample_info))
      names(choices_color)=c(colnames(sample_info))
      if ("Chain" %in% colnames(Big_mem_values$Short_DF) && "Chain" %!in% choices_color) {
        tmp_names_color=names(choices_color)
        choices_color=c(choices_color, "Chain")
        names(choices_color)=c(tmp_names_color, "Ig Chain")
      }
      tmp_names_choices_color=names(choices_color)
      choices_color=c(choices_color, "Best_V", "Best_J", "V_and_J")
      names(choices_color)=c(tmp_names_choices_color, "Main V gene isoform", "Main J gene isoform", "Main V and J gene isoforms")
      choices_color=choices_color[which(choices_color != "Additional_info")]
      choices_color=choices_color[which(choices_color != "Folder")]
      choices_color=choices_color[which(choices_color != "Filename")]
      # choices_color=choices_color[which(choices_color != "V_and_J")]



      if ("D_region" %in% colnames(Big_mem_values$Short_DF)) {
        tmp_names_color=names(choices_color)
        choices_color=c(choices_color, "Best_D")
        names(choices_color)=c(tmp_names_color, "Main D gene isoform")
      }
      if ("V_and_D_and_J" %in% colnames(Big_mem_values$Short_DF)) {
        tmp_names_color=names(choices_color)
        choices_color=c(choices_color, "V_and_D_and_J")
        names(choices_color)=c(tmp_names_color, "Main V, D and J gene isoforms")
      }
      if ("C_region" %in% colnames(Big_mem_values$Short_DF)) {
        tmp_names_color=names(choices_color)
        choices_color=c(choices_color, "C_region")
        names(choices_color)=c(tmp_names_color, "Constant region")
      }
      if ("Clone_ID" %in% colnames(Big_mem_values$Short_DF)) {
        tmp_names_color=names(choices_color)
        choices_color=c(choices_color, "Clone_ID")
        names(choices_color)=c(tmp_names_color, "Clone ID")
      }

      if ("Dominance" %in% colnames(Big_mem_values$Short_DF)) {
        tmp_names_color=names(choices_color)
        choices_color=c(choices_color, "Dominance")
        names(choices_color)=c(tmp_names_color, "Dominance")
      }


      selectInput("plot_color_feature","Color by", choices = choices_color, selected="Group")
    }

  })

  output$Violin_feature_plot <- renderPlotly({
    if(!is.null(Selection_values$rows) && length(Selection_values$rows)>0 && !is.null(Selection_values$columns) && length(Selection_values$columns)>0) {


      fig_violin=draw_feature_violinplot(values=Big_mem_values$Big_DF[,which(Big_mem_values$Header == input$plot_feature)],
                                         name_values=input$plot_feature,
                                         sequence_info_df=Big_mem_values$Short_DF,
                                         group_info=input$groups_selected,
                                         additional_group_info=input$plot_color_feature,
                                         show_reconstructed=input$show_reconstructed,
                                         selected_rows=Selection_values$rows,
                                         selected_subclones=NULL,
                                         hide_dots=input$hide_points)
      ID_selected_values$subclones=event_data("plotly_selected")
      print(ID_selected_values$clones$key)
      fig_violin
    }

  })





  # 5. Clones ####


  # in server.R create reactiveVal
  clonal_group_current_selection <- reactiveVal(NULL)

  # now store your current selection in the reactive value
  observeEvent(input$clonal_group, {
    clonal_group_current_selection(input$clonal_group)
  })

  output$clonal_group_output=renderUI({
    choices = c("")
    tmp_selection=""
    if(!is.null(Exploration_values$rows) && length(Exploration_values$rows)>0 && !is.null(Exploration_values$columns) && length(Exploration_values$columns)>0) {
      choices=c(colnames(Big_mem_values$Short_DF)[grepl("Clone_", colnames(Big_mem_values$Short_DF))])
      if(!is.null(clonal_group_current_selection()) && clonal_group_current_selection()[1] != "") {
        tmp_selection=clonal_group_current_selection()
      } else {
        tmp_selection=choices[1]
      }
    }

    selectInput("clonal_group","Use this clonal definition:", choices = choices, selected=tmp_selection)
  })

  new_dominance_calculation <- reactive({
    list(input$clonal_group,input$dominance_threshold)
  })

  observeEvent(new_dominance_calculation(), {
    if(!is.null(Selection_values$rows) && length(Selection_values$rows)>0 && !is.null(Selection_values$columns) && length(Selection_values$columns)>0) {
      if (input$clonal_group != "" ) {
        Big_mem_values$Short_DF[, "Dominance"] <- ""

        for(pat_sample in unique(Big_mem_values$Short_DF$Patient_Sample)) {
          index_tmp=which(Big_mem_values$Short_DF[, "Patient_Sample"]==pat_sample)
          dom_table=table(Big_mem_values$Short_DF[[input$clonal_group]][index_tmp])/length(Big_mem_values$Short_DF[[input$clonal_group]][index_tmp])
          dom_table=sapply(dom_table, function(z) if(z >= input$dominance_threshold/100){"Dominant"}else{"Non-dominant"})
          print(names(dom_table[1:5]))
          Big_mem_values$Short_DF[index_tmp, "Dominance"] <- dom_table[match(Big_mem_values$Short_DF[[input$clonal_group]][index_tmp],names(dom_table))]
        }

      }
    }
  })


  observeEvent(input$calculate_new_clone, {
    if(!is.null(Selection_values$rows) && length(Selection_values$rows)>0 && !is.null(Selection_values$columns) && length(Selection_values$columns)>0) {
      if (!(paste("Clone", input$clonal_level_group, input$clonal_region_group, input$new_clonal_group,input$identity_clonal_group, "simil", sep="_")  %in% colnames(Big_mem_values$Short_DF))  ) {
        print("2/2")
        Big_mem_values$Short_DF=calculate_clone(seq_df=Big_mem_values$Short_DF,  clonotype=input$new_clonal_group, AA_or_NT=input$clonal_level_group, region= input$clonal_region_group , percentage=input$identity_clonal_group, calculate_shared_clones=input$calculate_shared_clones)
        print("DONE")
        print(colnames(Big_mem_values$Short_DF)[(length(Big_mem_values$Short_DF)-10):length(Big_mem_values$Short_DF)])
      }
    }
  })



  output$Violin_plot <- renderPlotly({
    if(!is.null(Selection_values$rows) && length(Selection_values$rows)>0 && !is.null(Selection_values$columns) && length(Selection_values$columns)>0 && !is.null(input$clonal_group) && input$clonal_group != "") {


      # TOFIX ####
      fig_violin=draw_violinplots(seq_df=Big_mem_values$Short_DF, group="Patient_Sample", selected_rows=Selection_values$rows, clonotype=input$clonal_group, AA_or_NT=input$clonal_level_group,
                                  region= input$clonal_region_group , percentage=input$identity_clonal_group, freq_filter=input$filter_clonal_group, Selected_clones=ID_selected_values$clones, dominance_threshold = input$dominance_threshold)
      ID_selected_values$clones=event_data("plotly_selected")
      print(ID_selected_values$clones$key)
      fig_violin
    }

  })

  output$upset_plot <- upsetjs::renderUpsetjs({
    if(!is.null(Selection_values$rows) && length(Selection_values$rows)>0 && !is.null(Selection_values$columns) && length(Selection_values$columns)>0 && !is.null(input$clonal_group) && input$clonal_group != "") {

      fig_upset=draw_upsetplot(seq_df=Big_mem_values$Short_DF, group="Patient_Sample", selected_rows=Selection_values$rows, clonotype=input$clonal_group, AA_or_NT=input$clonal_level_group,
                               region= input$clonal_region_group , percentage=input$identity_clonal_group, freq_filter=input$filter_clonal_group, Selected_clones=ID_selected_values$clones)
      fig_upset %>%
        upsetjs::interactiveChart('click', events_nonce = TRUE) %>%
        upsetjs::generateDistinctIntersections()
    }

  })

  observeEvent(ID_selected_values$intersection_samples, {
    upsetjs::upsetjsProxy('upset_plot', session) %>%
      upsetjs::setSelection(ID_selected_values$intersection_samples)
  })


  observeEvent(input$upset_plot_click, {
    if(isTRUE(input$upset_plot_click[['isSelected']])) {
      ID_selected_values$intersection_samples <- ''
    } else {
      ID_selected_values$intersection_samples <- input$upset_plot_click
      print(ID_selected_values$intersection_samples$setNames)
      print(ID_selected_values$intersection_samples$elems)
    }
  })

  output$Comparison_plot <- renderPlotly({
    if(!is.null(ID_selected_values$intersection_samples$setNames) && !is.null(ID_selected_values$intersection_samples$elems) && !is.null(Selection_values$rows) && length(Selection_values$rows)>0 && !is.null(Selection_values$columns) && length(Selection_values$columns)>0 && !is.null(input$clonal_group) && input$clonal_group != "") {

      if(length(unlist(ID_selected_values$intersection_samples$setNames))>1){
        fig_sharedclones=draw_sharedclonesplot(seq_df=Big_mem_values$Short_DF, sets=ID_selected_values$intersection_samples$setNames, group="Patient_Sample", selected_rows=Selection_values$rows, clonotype=input$clonal_group, AA_or_NT=input$clonal_level_group,
                                               region= input$clonal_region_group , percentage=input$identity_clonal_group, freq_filter=input$filter_clonal_group, Selected_clones=ID_selected_values$clones, dominance_threshold=input$dominance_threshold)
        ID_selected_values$clones=event_data("plotly_selected")
        fig_sharedclones%>% toWebGL()
      }

    }

  })

  # 6.ML ####

  numbers <- reactive({
    validate(
      need(is.numeric(input$seedML), "Please input a number")
    )
  })
  output$onlynumbers <- renderPrint({ numbers() })

  hide("ML_clone_def_to_use")
  observeEvent(input$ML_only_one_subclone, {
    toggle("ML_clone_def_to_use")
  })

  output$ML_clone_def_to_use=renderUI({
    if(!is.null(Selection_values$rows) && length(Selection_values$rows)>0 && !is.null(Selection_values$columns) && length(Selection_values$columns)>0) {
      selectInput("MLclone","Clone definitions yet to generate, go to 5.", choices = NULL, selected=NULL)
    } else{
      selectInput("MLclone","Clone definition", choices = colnames(Big_mem_values$Short_DF)[which(startsWith(colnames(Big_mem_values$Short_DF), "Clone"))], selected=NULL)
    }

  })
  output$ML_group_selection=renderUI({
    if(!is.null(Selection_values$rows) && length(Selection_values$rows)>0 && !is.null(Selection_values$columns) && length(Selection_values$columns)>0) {
      selectInput("groups_selected","Separate", choices = c("Group"), selected="Group")
    } else{
      selectInput("groups_selected","Separate", choices = c("Group", "Patients"), selected="Group")
    }

  })


  output$ML_method_selection=renderUI({
    if(!is.null(Selection_values$rows) && length(Selection_values$rows)>0 && !is.null(Selection_values$columns) && length(Selection_values$columns)>0) {
      selectInput("groups_selected_ML","Work only with", choices = c("Linear regression"), selected="Linear regression")
    } else{
      if((as.numeric(unname(strsplit(as.character(print(get_ram(), unit_system = "metric")), split=" ")[[1]][1])) >   round(length(Selection_values$rows)*length(Selection_values$columns)*8/2^{20}/1024, 2))) {
        selectInput("groups_selected_ML","Work only with", choices = c("Linear regression"), selected="Linear regression")
      } else {
        selectInput("groups_selected_ML","Work only with", choices = c("Random Forest", "Linear regression", "SVM"), selected="Random Forest")
      }

    }

  })

  output$ML_assign_sample_to_training=renderUI({
    if(!is.null(Selection_values$rows) && length(Selection_values$rows)>0 && !is.null(Selection_values$columns) && length(Selection_values$columns)>0) {

    } else{
      sample_info<-sample_info_react$table
      selectizeInput(
        'patients_to_exclude', 'Select the subjects to include only in the ML validation (if used, this overrides the % size of samples destined to training and validation sets as the other subjects will be used as training)', choices = sample_info$Patient, selected=NULL, multiple = TRUE,
        options = list(plugins= list('remove_button')))
    }

  })

  ML_values<-reactiveValues(models_generated=list(),current_training_rows=NULL,current_validation_rows=NULL, group_category=NULL)

  observeEvent(input$update_ML_group_parameters, {
    set.seed(input$seedML)
    groupselection=input$groups_selected_ML

    tmp_rows=Selection_values$rows
    training_rows=c()
    validation_rows=c()

    if(input$ML_only_one_subclone && !is.null(input$MLclone)){
      clonedef_to_use=input$MLclone
      tmp_rows=tmp_rows[which(Big_mem_values$Short_DF[Selection_values$rows,..clonedef_to_use])]
    }

    if(input$ML_use_percentage_selection_dataset != 100) {
      tmp_decimated_rows=c()

      for(subject in unique(Big_mem_values$Short_DF[tmp_rows,get("Patient")])) {

        tmp_intersection=tmp_rows[which(Big_mem_values$Short_DF[tmp_rows,get("Patient")]==subject)]

        tmp_decimated_rows=c(tmp_decimated_rows, sample(tmp_intersection, round(input$ML_use_percentage_selection_dataset * length(tmp_intersection)/100)) )
      }


      tmp_rows=tmp_decimated_rows
    }


    if(!is.null(input$patients_to_exclude)){

    } else {
      for(subject in input$patients_to_exclude) {
        validation_rows=c(validation_rows, tmp_rows[which(Big_mem_values$Short_DF[tmp_rows,get("Patient")]==subject)])
      }

    }

    if (input$ML_equal_VJs && input$ML_equal_set_sizes) {

      VJs=c()
      sizes=list()
      for (group in unique(Big_mem_values$Short_DF[tmp_rows,..groupselection])) {
        tmp_group_rows=tmp_rows[which(Big_mem_values$Short_DF[tmp_rows,..groupselection]==group)]
        tmp_VJ=unique(paste(Big_mem_values$Short_DF[tmp_group_rows,get("Best_V")], Big_mem_values$Short_DF[tmp_group_rows,get("Best_J")], sep="_&_"))
        if(length(VJs)==0){
          VJs=tmp_VJ
        } else {
          VJs=intersect(VJs, tmp_VJ)
        }
        sizes[[length(sizes)+1]]=table(paste(Big_mem_values$Short_DF[tmp_group_rows,get("Best_V")], Big_mem_values$Short_DF[tmp_group_rows,get("Best_J")], sep="_&_"))
      }

      VJ_sizes=rep(-100, length(VJs))
      names(VJ_sizes)=VJs
      for (i in c(1:length(sizes))) {
        tmp_sub_sel=sizes[[i]][which(names(sizes[[i]] %in% names(VJ_sizes)))]
        for (j in c(1:length(VJ_sizes))){
          if(VJ_sizes[j] == -100) {
            VJ_sizes[j]=tmp_sub_sel[names(VJ_sizes)][j]
          } else {
            VJ_sizes[j]=min(tmp_sub_sel[names(VJ_sizes)][j], VJ_sizes[j])
          }

        }

      }



      tmp_new_tmp_rows=c()
      for (group in unique(Big_mem_values$Short_DF[tmp_rows,..groupselection])) {
        tmp_group_rows=tmp_rows[which(Big_mem_values$Short_DF[tmp_rows,..groupselection]==group)]
        tmp_VJ=paste(Big_mem_values$Short_DF[tmp_group_rows,get("Best_V")], Big_mem_values$Short_DF[tmp_group_rows,get("Best_J")], sep="_&_")

        for (j in c(1:length(VJ_sizes))) {
          tmp_new_tmp_rows=c(tmp_new_tmp_rows, sample(tmp_group_rows[which(tmp_VJ %in% names(VJ_sizes)[j])],VJ_sizes[j]))

        }

      }


      for (group in unique(Big_mem_values$Short_DF[tmp_rows,..groupselection])) {
        tmp_group_rows=tmp_rows[which(Big_mem_values$Short_DF[tmp_rows,..groupselection]==group)]
        tmp_VJ=paste(Big_mem_values$Short_DF[tmp_group_rows,get("Best_V")], Big_mem_values$Short_DF[tmp_group_rows,get("Best_J")], sep="_&_")
        tmp_new_tmp_rows=c(tmp_new_tmp_rows, tmp_group_rows[which(tmp_VJ %in% VJs)])
      }
      tmp_rows=tmp_new_tmp_rows


    } else if(input$ML_equal_VJs){

      VJs=c()
      for (group in unique(Big_mem_values$Short_DF[tmp_rows,..groupselection])) {
        tmp_group_rows=tmp_rows[which(Big_mem_values$Short_DF[tmp_rows,..groupselection]==group)]
        tmp_VJ=unique(paste(Big_mem_values$Short_DF[tmp_group_rows,get("Best_V")], Big_mem_values$Short_DF[tmp_group_rows,get("Best_J")], sep="_&_"))
        if(length(VJs)==0){
          VJs=tmp_VJ
        } else {
          VJs=intersect(VJs, tmp_VJ)
        }
      }
      tmp_new_tmp_rows=c()
      for (group in unique(Big_mem_values$Short_DF[tmp_rows,..groupselection])) {
        tmp_group_rows=tmp_rows[which(Big_mem_values$Short_DF[tmp_rows,..groupselection]==group)]
        tmp_VJ=paste(Big_mem_values$Short_DF[tmp_group_rows,get("Best_V")], Big_mem_values$Short_DF[tmp_group_rows,get("Best_J")], sep="_&_")
        tmp_new_tmp_rows=c(tmp_new_tmp_rows, tmp_group_rows[which(tmp_VJ %in% VJs)])
      }
      tmp_rows=tmp_new_tmp_rows

    } else if(input$ML_equal_set_sizes){
      sizes=c()
      for (group in unique(Big_mem_values$Short_DF[tmp_rows,..groupselection])) {
        tmp_group_rows=tmp_rows[which(Big_mem_values$Short_DF[tmp_rows,..groupselection]==group)]
        sizes=c(sizes, length(tmp_group_rows))
      }
      sizes=min(sizes)

      tmp_new_tmp_rows=c()
      for (group in unique(Big_mem_values$Short_DF[tmp_rows,..groupselection])) {
        tmp_group_rows=tmp_rows[which(Big_mem_values$Short_DF[tmp_rows,..groupselection]==group)]
        stmp_new_tmp_rows=c(tmp_new_tmp_rows, sample(tmp_group_rows, sizes))
      }
      tmp_rows=tmp_new_tmp_rows
    }

    Big_mem_values$Short_DF[tmp_rows,..groupselection]


    ## ML_equal_VJs  ML_equal_set_sizes
    validation_rows=intersect(validation_rows, tmp_rows)
    size_training=round(input$ML_perc_training * (length(tmp_rows))/100)
    size_validation_compensate=length(tmp_rows)-size_training-length(validation_rows)
    if(size_validation_compensate<0) {
      size_validation_compensate=length(validation_rows)
    }
    percentage_to_compensate=(size_training-size_validation_compensate)/size_training
    training_rows=sample(tmp_rows[tmp_rows %!in% validation_rows], round(percentage_to_compensate* (length(tmp_rows[tmp_rows %!in% validation_rows]))/100))
    validation_rows=tmp_rows[tmp_rows  %!in% training_rows]
    ML_values$current_training_rows=training_rows
    ML_values$current_validation_rows=validation_rows
    ML_values$group_category=groupselection
  })



  output$ML_sunburst_dataset_ML <- renderUI({
    if(!is.null(Selection_values$rows) && length(Selection_values$rows)>0 && !is.null(Selection_values$columns) && length(Selection_values$columns)>0 && length(ML_values$current_validation_rows) >0 && length(ML_values$current_validation_rows) >0) {
      rows=c(ML_values$current_training_rows, ML_values$current_validation_rows)

      tmp_ML=c(rep("Training", length(ML_values$current_training_rows)), rep("Validation", length(ML_values$current_validation_rows)))
      tmp_groups=Big_mem_values$Short_DF[rows, ..ML_values$group_category]
      tmp_subjects=Big_mem_values$Short_DF[rows,get("Patient")]

      tmp_sunburst=paste(tmp_ML, tmp_subjects, tmp_groups, sep="THISISARANDOMSTRING")
      tmp_sunburst=gsub("-","_",tmp_sunburst)
      tmp_sunburst=gsub("THISISARANDOMSTRING","-",tmp_sunburst)
      tmp_sunburst=table(tmp_sunburst)
      df_sunburst=data.frame(V1=names(tmp_sunburst),V2=unname(tmp_sunburst))


      sund2b(df_sunburst,showLabels =F, rootLabel= "Training and validation set composition")

    } else {
      # HTML("ERROR: the set(s) are empty!")
    }


  })

  output$ML_sunburst_dataset_grouping <- renderUI({
    if(!is.null(Selection_values$rows) && length(Selection_values$rows)>0 && !is.null(Selection_values$columns) && length(Selection_values$columns)>0 && length(ML_values$current_validation_rows) >0 && length(ML_values$current_validation_rows) >0) {

      rows=c(ML_values$current_training_rows, ML_values$current_validation_rows)

      tmp_groups=Big_mem_values$Short_DF[rows, ..ML_values$group_category]
      tmp_subjects=Big_mem_values$Short_DF[rows,get("Patient")]
      tmp_VJ=paste(Big_mem_values$Short_DF[rows,get("Best_V")], Big_mem_values$Short_DF[rows,get("Best_J")], sep="_&_")

      tmp_sunburst=paste(tmp_groups, tmp_subjects, tmp_VJ, sep="THISISARANDOMSTRING")
      tmp_sunburst=gsub("-","_",tmp_sunburst)
      tmp_sunburst=gsub("THISISARANDOMSTRING","-",tmp_sunburst)
      tmp_sunburst=table(tmp_sunburst)
      df_sunburst=data.frame(V1=names(tmp_sunburst),V2=unname(tmp_sunburst))


      sund2b(df_sunburst,showLabels =F, rootLabel= "Total group composition")

    } else {
      # HTML("ERROR: the set(s) are empty!")
    }


  })



  # 7.Help ####

  HELP_values<-reactiveValues(stage=0)
  output$HELP_output <- renderUI({
    if(HELP_values$stage==0) {
      HTML("In this tab (0.Project information) we will define what type of project, samples and workspace will be used. You need to indicate:<br>
           <ul>
              <li>If your data is <b>TCR or BCR</b> repertoire-based.</li>
              <li>Tell us about your data. <b>Include the required descriptive table</b>.</li>
              <li>Specify <b>the folder where all the related-files will be created</b>.</li>
          </ul>
           ")
    } else if(HELP_values$stage==1) {
      HTML("In this tab (1.AIRR-Seq conversion) we will reconstruct the repertoire germlines based on your input. You need to indicate:<br>
           <ul>
              <li>The <b>folder where your data is located</b>. All files must be in the same folder. Filenames should match the filenames you indicated in the previous step. If they are found in the folder, it will say 'Yes' in Found_in_folder, otherwise 'No' and you need to either a) change the folder, b) move the file to the folder, c) change the input table or d) change the filename.</li>
              <li>Additional information to better<b>reconstruct the germlines</b>.</li>
          </ul>
           ")
    }else if(HELP_values$stage==2) {

    }else if(HELP_values$stage==2.5) {

    }else if(HELP_values$stage==3) {

    } else {
      "Strange."
    }
  })


  # 8.Background functions ####
  #Disable menuitem when the app loads

  hideTab(inputId = "menutabset", target = "1.AIRR-Seq conversion")
  hideTab(inputId = "menutabset", target = "2.Sequence feature determination")
  hideTab(inputId = "menutabset", target = "1&2. Select your bigmem files")
  hideTab(inputId = "menutabset", target = "3.Dataset exploration and variable selection")
  hideTab(inputId = "menutabset", target = "4.Feature exploration")
  hideTab(inputId = "menutabset", target = "5.Clonal exploration")
  hideTab(inputId = "menutabset", target = "5.Cell exploration")
  hideTab(inputId = "menutabset", target = "6.ML analysis")


  observeEvent(input$Move_to_1_airr, {
    HELP_values$stage=1
    hideTab(inputId = "menutabset", target = "0.Project information")
    showTab(inputId = "menutabset", target = "1.AIRR-Seq conversion", select=T)
  })

  observeEvent(input$Move_to_2, {
    HELP_values$stage=2
    hideTab(inputId = "menutabset", target = "1.AIRR-Seq conversion")
    showTab(inputId = "menutabset", target = "2.Sequence feature determination", select=T)
    can_show_button_to_step$step_2=FALSE
    can_show_button_to_step$step_2_AIRR=FALSE
  })

  observeEvent(input$Move_to_analysis, {
    HELP_values$stage=2.5
    hideTab(inputId = "menutabset", target = "0.Project information")
    hideTab(inputId = "menutabset", target = "2.Sequence feature determination")
    showTab(inputId = "menutabset", target = "1&2. Select your bigmem files", select=T)
  })

  observeEvent(input$Move_to_analysis_real, {
    HELP_values$stage=3
    hideTab(inputId = "menutabset", target = "0.Project information")
    hideTab(inputId = "menutabset", target = "2.Sequence feature determination")
    hideTab(inputId = "menutabset", target = "1&2. Select your bigmem files")
    showTab(inputId = "menutabset", target = "3.Dataset exploration and variable selection", select=T)
    showTab(inputId = "menutabset", target = "4.Feature exploration")
    showTab(inputId = "menutabset", target = "5.Clonal exploration")
    if (Selection_values$Cell) {
      showTab(inputId = "menutabset", target = "5.Cell exploration")
    }

    showTab(inputId = "menutabset", target = "6.ML analysis")
    if (!all(is.numeric(unlist(input$bigmem_folder))) ) {
      folder_values$Featured=parseDirPath(volumes(), input$bigmem_folder)
    }
  })

  observeEvent(input$Reset_all, {
    session$reload()
    rm(list = names(input))
    HELP_values$stage=0
    showTab(inputId = "menutabset", target = "0.Project information", select=T)
    hideTab(inputId = "menutabset", target = "1.AIRR-Seq conversion")
    hideTab(inputId = "menutabset", target = "2.Sequence feature determination")
    hideTab(inputId = "menutabset", target = "3.Dataset exploration and variable selection")
    hideTab(inputId = "menutabset", target = "4.Feature exploration")
    hideTab(inputId = "menutabset", target = "5.Clonal exploration")
    hideTab(inputId = "menutabset", target = "5.Cell exploration")
    hideTab(inputId = "menutabset", target = "6.ML analysis")
    can_show_button_to_step$step_2=FALSE
    can_show_button_to_step$step_2_AIRR=FALSE
    can_show_button_to_step$step_3=FALSE
    folder_values$Featured=""
    folder_values$AIRR_parsed=""
    Big_mem_values$Header=NULL
    Big_mem_values$Short_DF=NULL
    Big_mem_values$Big_DF=NULL
    Exploration_values$rows=NULL
    Exploration_values$columns=NULL
    Exploration_values$Scores=NULL
    Exploration_values$Variance_explained=NULL
    Exploration_values$UMAP=NULL
    Exploration_values$Parameters=NULL
    Selection_values$rows=NULL
    Selection_values$columns=NULL
    Selection_values$Scores=NULL
    Selection_values$Variance_explained=NULL
    Selection_values$UMAP=NULL
    Selection_values$Parameters=NULL
    Selection_values$Cell=F
    ID_selected_values$subclones=NULL
    ID_selected_values$clones=NULL


  })

  print(session$userData)


  session$allowReconnect(TRUE)
}
