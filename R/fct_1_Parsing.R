#' 1_Parsing
#'
#' @description A fct function
#'
#' @return The return value, if any, from executing the function.
#' @import stringr
#'
#' @keywords internal
#' @export
parse_AIRRSeq_file <- function(
    file, group, patient, subgroup, sample, input_path, C_region_included, FWR1partial,
    FWR4partial, D_gene, repertoire, output_path, is_example=F
) {


    ### Am I working with the insertions?

    AIRR_input <- read.table(
        file = paste(
            input_path, paste(file, ".tsv", sep = ""),
            sep = ""
        ), nrows=if(is_example){20} else {-1},
        sep = "\t", header = TRUE, fill = T
    )

    AIRR_input[is.na(AIRR_input)] = ""

    if(is_example){
      # AIRR_input=AIRR_input[c(1:min(nrow(AIRR_input),10)),]

    }
    # Removing Ns
    AIRR_input <- AIRR_input[which(!grepl("N", AIRR_input$sequence)),
        ]
    AIRR_input <- AIRR_input[which(!grepl("n", AIRR_input$sequence)),
        ]
    AIRR_input <- AIRR_input[which(AIRR_input$sequence_alignment != ""),
        ]

    AIRR_input <- AIRR_input[which(AIRR_input$germline_alignment != ""),
    ]
    AIRR_input <- AIRR_input[which(AIRR_input$cdr2 != ""),
    ]
    AIRR_input <- AIRR_input[which(AIRR_input$cdr1 != ""),
    ]
    AIRR_input <- AIRR_input[which(AIRR_input$fwr2 != ""),
    ]
    AIRR_input <- AIRR_input[which(AIRR_input$cdr3 != ""),
    ]
    AIRR_input <- AIRR_input[which(AIRR_input$fwr3 != ""),
    ]
    AIRR_input <- AIRR_input[which(AIRR_input$fwr4 != ""),
    ]
    ##some programmes use - for deletions, others dont, so we standarize this way
    AIRR_input$sequence_alignment=sapply(AIRR_input$sequence_alignment, function(x) gsub("-",".", x, fixed=T))
    AIRR_input$germline_alignment=sapply(AIRR_input$germline_alignment, function(x) gsub("-",".", x, fixed=T))
    AIRR_input$fwr1=sapply(AIRR_input$fwr1, function(x) gsub("-",".", x, fixed=T))
    AIRR_input$fwr2=sapply(AIRR_input$fwr2, function(x) gsub("-",".", x, fixed=T))
    AIRR_input$fwr3=sapply(AIRR_input$fwr3, function(x) gsub("-",".", x, fixed=T))
    AIRR_input$fwr4=sapply(AIRR_input$fwr4, function(x) gsub("-",".", x, fixed=T))
    AIRR_input$cdr1=sapply(AIRR_input$cdr1, function(x) gsub("-",".", x, fixed=T))
    AIRR_input$cdr2=sapply(AIRR_input$cdr2, function(x) gsub("-",".", x, fixed=T))
    AIRR_input$cdr3=sapply(AIRR_input$cdr3, function(x) gsub("-",".", x, fixed=T))

    regions <- c("NT_FWR1", "NT_CDR1", "NT_FWR2", "NT_CDR2", "NT_FWR3", "NT_CDR3", "NT_FWR4")

    ########## R: insertions are not included here lalalala
    IMGT_parsed_index <- data.frame(
        ID = AIRR_input$sequence_id, Patient = patient, Sample = sample, Group = group,
        Sequence_type = "Repertoire", Subgroup = subgroup, V_region = AIRR_input$v_call,
        J_region = AIRR_input$j_call, NT_FWR1 = AIRR_input$fwr1, NT_CDR1 = AIRR_input$cdr1,
        NT_FWR2 = AIRR_input$fwr2, NT_CDR2 = AIRR_input$cdr2, NT_FWR3 = AIRR_input$fwr3,
        NT_CDR3 = AIRR_input$cdr3, NT_FWR4 = AIRR_input$fwr4, stringsAsFactors = F
    )


    if ("c_call" %in% colnames(AIRR_input)) {
        IMGT_parsed_index$C_region <- AIRR_input$c_call
    }
    if ("d_call" %in% colnames(AIRR_input)) {
        IMGT_parsed_index$D_region <- AIRR_input$d_call
    }

    if ("clone_id" %in% colnames(AIRR_input) &&
        any(!is.na(unique(AIRR_input$clone_id))) &&
        any(
            unique(AIRR_input$clone_id) !=
                " "
        ) &&
        any(
            unique(AIRR_input$clone_id) !=
                ""
        )) {

        IMGT_parsed_index$Clone_ID <- AIRR_input$clone_id
    } else if ("raw_clonotype_id" %in% colnames(AIRR_input)) {
        IMGT_parsed_index$Clone_ID <- AIRR_input$raw_clonotype_id
    } else if ("clone_id" %in% colnames(AIRR_input)) {
        IMGT_parsed_index$Clone_ID <- AIRR_input$clone_id
    }
    if ("cell_id" %in% colnames(AIRR_input)) {
        IMGT_parsed_index$Cell_ID <- AIRR_input$cell_id
    }

    if ("chain" %in% colnames(AIRR_input)) {
        IMGT_parsed_index$Chain <- AIRR_input$chain
    }
    if ("locus" %in% colnames(AIRR_input)) {
        IMGT_parsed_index$Chain <- AIRR_input$locus
    }



    IMGT_parsed_index$ORF_begins=sapply(c(1:nrow(AIRR_input)), function(z)
      if(any(str_locate_all(AIRR_input$sequence_alignment[z],".")[[1]][1,1]==1) || any(str_locate_all(AIRR_input$sequence_alignment[z],"-")[[1]][1,1]==1)){
        if(FWR1partial){
          ((str_locate_all(tolower(AIRR_input$sequence_alignment[z]),"a|c|g|t")[[1]][1,2] - nchar(gsub(".","",AIRR_input$fwr1[z], fixed=T)))%%3)
        } else {
          (str_locate_all(tolower(AIRR_input$sequence_alignment[z]),"a|c|g|t")[[1]][1,2]%%3)
        }

      }else {
        1
      } )



    ### even if the format says aligned, not all programmes (e.g. IMGT,
    ### 02/2024) have these fields aligned. Let's add a checkup
    if (any(grepl(".", IMGT_parsed_index$NT_FWR1, fixed = T)) ||
        any(grepl(".", IMGT_parsed_index$NT_FWR2, fixed = T)) ||
        any(grepl(".", IMGT_parsed_index$NT_FWR3, fixed = T)) ||
        any(grepl(".", IMGT_parsed_index$NT_FWR4, fixed = T)) ||
        any(grepl(".", IMGT_parsed_index$NT_CDR1, fixed = T)) ||
        any(grepl(".", IMGT_parsed_index$NT_CDR2, fixed = T)) ||
        any(grepl(".", IMGT_parsed_index$NT_CDR3, fixed = T))) {
      aligned_regions <- T
    } else {
      aligned_regions <- F
    }

    ########## R: insertions are not included here lalalala
    if (aligned_regions){
      or_CDR1 <- 1 + nchar(AIRR_input$fwr1)
      or_FWR2 <- or_CDR1 + nchar(AIRR_input$cdr1)
      or_CDR2 <- or_FWR2 + nchar(AIRR_input$fwr2)
      or_FWR3 <- or_CDR2 + nchar(AIRR_input$cdr2)
      or_CDR3 <- or_FWR3 + nchar(AIRR_input$fwr3)
      or_FWR4 <- or_CDR3 + nchar(AIRR_input$cdr3)
    } else {
      # or_CDR1 <- 1 + as.numeric(sapply(AIRR_input$fwr1_end, function(z) if(z ==""){0}else{z}))
      # or_FWR2 <- 1 + as.numeric(sapply(AIRR_input$cdr1_end, function(z) if(z ==""){0}else{z}))
      # or_CDR2 <- 1 + as.numeric(sapply(AIRR_input$fwr2_end, function(z) if(z ==""){0}else{z}))
      # or_FWR3 <- 1 + as.numeric(sapply(AIRR_input$cdr2_end, function(z) if(z ==""){0}else{z}))
      # or_CDR3 <- 1 + as.numeric(sapply(AIRR_input$fwr3_end, function(z) if(z ==""){0}else{z}))
      # or_FWR4 <- 1 + as.numeric(sapply(AIRR_input$cdr3_end, function(z) if(z ==""){0}else{z}))


      or_CDR1 <- c()
      or_FWR2 <- c()
      or_CDR2 <- c()
      or_FWR3 <- c()
      or_CDR3 <- c()
      or_FWR4 <- c()
      end_FWR4 <- c()
      for (numrow in c(1:nrow(AIRR_input))) {
        numbers=c(1:nchar(AIRR_input$sequence_alignment[numrow]))
        characters=strsplit(AIRR_input$sequence_alignment[numrow], split="")[[1]]
        gap_positions=str_locate_all(gsub(".","!",AIRR_input$sequence_alignment[numrow], fixed=T),"!")[[1]][,1]
        characters=characters[which(numbers %!in% gap_positions)]
        names(characters)=numbers[which(numbers %!in% gap_positions)]

        tmp_pos=c()
        init_char=1
        for(region in c("fwr1", "cdr1","fwr2","cdr2","fwr3","cdr3", "fwr4")) {
          char_pos=str_locate( paste(characters[init_char:length(characters)], collapse=""),AIRR_input[numrow,region])[1,2] + 1
          if(length(tmp_pos)>0 && char_pos==2) {

            char_pos=1
          } else if (region == "fwr4") {
            char_pos=char_pos-1
          } else if (region == "fwr1" && AIRR_input[numrow,region]=="") {
            char_pos=1
          }
          tmp_pos=c(tmp_pos,as.numeric(names(characters[init_char:length(characters)])[char_pos]))

          init_char=c(init_char:length(characters))[char_pos]


        }
        or_CDR1 <- c(or_CDR1, tmp_pos[1])
        or_FWR2 <- c(or_FWR2, tmp_pos[2])
        or_CDR2 <- c(or_CDR2, tmp_pos[3])
        or_FWR3 <- c(or_FWR3,tmp_pos[4] )
        or_CDR3 <- c(or_CDR3, tmp_pos[5])
        or_FWR4 <- c(or_FWR4, tmp_pos[6])
        end_FWR4 <- c(end_FWR4, tmp_pos[7] )
      }
    }



    AIRR_input$germline_alignment_nospace <- gsub(".", "", AIRR_input$germline_alignment, fixed = T)

    AIRR_input$germline_alignment_nospace[which(startsWith(AIRR_input$sequence_alignment,"."))]=
      sapply(which(startsWith(AIRR_input$sequence_alignment,".")),
             function(z)
               gsub(".", "",
               paste(strsplit(AIRR_input$germline_alignment[z],
                              split="")[[1]][str_locate_all(tolower(AIRR_input$sequence_alignment[z]),
                                                            "a|c|g|t")[[1]][1,2]:
                                               nchar(AIRR_input$germline_alignment[z])],
               collapse=""),
              fixed = T)
             )


    if (aligned_regions) {
        NT_CDR3_germ <- strsplit(
            substring(AIRR_input$germline_alignment, or_CDR3, or_CDR3 - 1 + nchar(AIRR_input$cdr3)),
            split = ""
        )
    } else {
        NT_CDR3_germ <- strsplit(
            substring(
                AIRR_input$germline_alignment, or_CDR3, or_FWR4 - 1
            ),
            split = ""
        )
    }
    # NT_CDR3_germ=substring(AIRR_input$sequence_alignment,or_CDR3,or_CDR3-1+nchar(AIRR_input$cdr3))


    for (index_cdr3 in c(1:length(NT_CDR3_germ))) {
        tmp <- NT_CDR3_germ[[index_cdr3]]
        index_n <- c(
            which(tmp == "N"),
            which(tmp == "n")
        )
        if (D_gene == F) {
          if(length(index_n)!=0) {
            index_n <- c(min(index_n):max(index_n))
          }
        }
        if (aligned_regions) {
            tmp[index_n] <- strsplit(
                substring(
                  AIRR_input$sequence_alignment[index_cdr3], or_CDR3[index_cdr3],
                  or_CDR3[index_cdr3] - 1 + nchar(AIRR_input$cdr3[index_cdr3])
              ),
                split = ""
            )[[1]][index_n]
        } else {
            tmp[index_n] <- strsplit(
                substring(
                  AIRR_input$sequence_alignment[index_cdr3],
                  or_CDR3[index_cdr3], or_FWR4[index_cdr3] - 1
              ),
                split = ""
            )[[1]][index_n]
        }
        NT_CDR3_germ[[index_cdr3]] <- paste(tmp, collapse = "")
    }



    if (aligned_regions) {
        reconstructed_IMGT_parsed_index <- data.frame(
            ID = AIRR_input$sequence_id, Patient = patient, Sample = sample, Group = group,
            Sequence_type = "Reconstructed_germline", Subgroup = subgroup, V_region = AIRR_input$v_call,
            J_region = AIRR_input$j_call, NT_FWR1 = gsub(
                "-", "", gsub(
                  ".", "", substring(AIRR_input$germline_alignment, 1, nchar(AIRR_input$fwr1)),
                  fixed = T
              ),
                fixed = T
            ),
            NT_CDR1 = gsub(
                "-", "", gsub(
                  ".", "", substring(AIRR_input$germline_alignment, or_CDR1, or_CDR1 - 1 + nchar(AIRR_input$cdr1)),
                  fixed = T
              ),
                fixed = T
            ),
            NT_FWR2 = gsub(
                "-", "", gsub(
                  ".", "", substring(AIRR_input$germline_alignment, or_FWR2, or_FWR2 - 1 + nchar(AIRR_input$fwr2)),
                  fixed = T
              ),
                fixed = T
            ),
            NT_CDR2 = gsub(
                "-", "", gsub(
                  ".", "", substring(AIRR_input$germline_alignment, or_CDR2, or_CDR2 - 1 + nchar(AIRR_input$cdr2)),
                  fixed = T
              ),
                fixed = T
            ),
            NT_FWR3 = gsub(
                "-", "", gsub(
                  ".", "", substring(AIRR_input$germline_alignment, or_FWR3, or_FWR3 - 1 + nchar(AIRR_input$fwr3)),
                  fixed = T
              ),
                fixed = T
            ),
            NT_CDR3 = unlist(NT_CDR3_germ),
            NT_FWR4 = gsub(
                "-", "", gsub(
                  ".", "", substring(AIRR_input$germline_alignment, or_FWR4, or_FWR4 - 1 + nchar(AIRR_input$fwr4)),
                  fixed = T
              ),
                fixed = T
            ),
            stringsAsFactors = F
        )
    } else {
        reconstructed_IMGT_parsed_index <- data.frame(
            ID = AIRR_input$sequence_id, Patient = patient, Sample = sample, Group = group,
            Sequence_type = "Reconstructed_germline", Subgroup = subgroup, V_region = AIRR_input$v_call,
            J_region = AIRR_input$j_call, NT_FWR1 = gsub(
                "-", "", gsub(
                  ".", "", substring(AIRR_input$germline_alignment, 1, or_CDR1-1),
                  fixed = T
              ),
                fixed = T
            ),
            NT_CDR1 = gsub(
                "-", "", gsub(
                  ".", "", substring(
                    AIRR_input$germline_alignment, or_CDR1, or_FWR2 - 1
                ),
                  fixed = T
              ),
                fixed = T
            ),
            NT_FWR2 = gsub(
                "-", "", gsub(
                  ".", "", substring(
                    AIRR_input$germline_alignment, or_FWR2, or_CDR2 - 1
                ),
                  fixed = T
              ),
                fixed = T
            ),
            NT_CDR2 = gsub(
                "-", "", gsub(
                  ".", "", substring(
                    AIRR_input$germline_alignment, or_CDR2, or_FWR3 - 1
                ),
                  fixed = T
              ),
                fixed = T
            ),
            NT_FWR3 = gsub(
                "-", "", gsub(
                  ".", "", substring(
                    AIRR_input$germline_alignment, or_FWR3, or_CDR3 - 1
                ),
                  fixed = T
              ),
                fixed = T
            ),
            NT_CDR3 = unlist(NT_CDR3_germ),
            NT_FWR4 = gsub(
                "-", "", gsub(
                  ".", "", substring(
                    AIRR_input$germline_alignment, or_FWR4, end_FWR4
                ),
                  fixed = T
              ),
                fixed = T
            ),
            stringsAsFactors = F
        )
    }


    if ("c_call" %in% colnames(AIRR_input)) {
        reconstructed_IMGT_parsed_index$C_region <- AIRR_input$c_call
    }
    if ("d_call" %in% colnames(AIRR_input)) {
        reconstructed_IMGT_parsed_index$D_region <- AIRR_input$d_call
    }

    if ("raw_clonotype_id" %in% colnames(AIRR_input)) {
        reconstructed_IMGT_parsed_index$Clone_ID <- AIRR_input$raw_clonotype_id
    }
    if ("clone_id" %in% colnames(AIRR_input)) {
        reconstructed_IMGT_parsed_index$Clone_ID <- IMGT_parsed_index$Clone_ID
    }
    if ("cell_id" %in% colnames(AIRR_input)) {
        reconstructed_IMGT_parsed_index$Cell_ID <- AIRR_input$cell_id
    }
    if ("chain" %in% colnames(AIRR_input)) {
        reconstructed_IMGT_parsed_index$Chain <- AIRR_input$chain
    }

    if ("locus" %in% colnames(AIRR_input)) {
        reconstructed_IMGT_parsed_index$Chain <- AIRR_input$locus
    }


    if (C_region_included) {
        index_insertions <- which(
            !(sapply(
                c(1:nrow(AIRR_input)),
                function(z) grepl(
                  gsub("-", "", gsub(".", "", AIRR_input$sequence_alignment[z], fixed = T)),
                  AIRR_input$sequence[z]
              )
            ))
        )

    } else {
        index_insertions <- which(
            (sapply(
                c(1:nrow(AIRR_input)),
                function(z) (gsub("-", "", gsub(".", "", AIRR_input$sequence_alignment[z], fixed = T)) !=
                  AIRR_input$sequence[z])
            ))
        )
    }
    index_deletions_germ <- which(
      (sapply(
        c(1:nrow(AIRR_input)),
        function(z) any(str_locate_all(gsub(".","!",AIRR_input$germline_alignment[z], fixed=T),"!")[[1]][,1] %!in%
                          str_locate_all(gsub(".","!",AIRR_input$sequence_alignment[z], fixed=T),"!")[[1]][,1])
      ))
    )

    index_deletions_rep <- which(
      (sapply(
        c(1:nrow(AIRR_input)),
        function(z) any(str_locate_all(gsub(".","!",AIRR_input$sequence_alignment[z], fixed=T),"!")[[1]][,1] %!in%
                          str_locate_all(gsub(".","!",AIRR_input$germline_alignment[z], fixed=T),"!")[[1]][,1])
      ))
    )

    insertion_sizes <- sapply(
        index_insertions, function(z) drop(
            attr(
                adist(
                  gsub("-", "", gsub(".", "", AIRR_input$sequence_alignment[z], fixed = T)),
                  AIRR_input$sequence[z], counts = TRUE
              ),
                "counts"
            )
        )[which(
            names(
                drop(
                  attr(
                    adist(
                      gsub("-", "", gsub(".", "", AIRR_input$sequence_alignment[z], fixed = T)),
                      AIRR_input$sequence[z], counts = TRUE
                  ),
                    "counts"
                )
              )
            ) ==
                "ins"
        )]
    )


    for (region in regions) {
      IMGT_parsed_index[, region] <- toupper(
        gsub(
          "-", "", gsub(".", "", IMGT_parsed_index[, region], fixed = T),
          fixed = T
        )
      )
      reconstructed_IMGT_parsed_index[, region] <- toupper(
        gsub(
          "-", "", gsub(".", "", reconstructed_IMGT_parsed_index[, region], fixed = T),
          fixed = T
        )
      )
    }


    for (deletion_seq in index_deletions_germ) {
      all_germ_located=str_locate_all(gsub(".","!",AIRR_input$germline_alignment[deletion_seq], fixed=T),"!")[[1]][,1]
      all_rep_located= str_locate_all(gsub(".","!",AIRR_input$sequence_alignment[deletion_seq], fixed=T),"!")[[1]][,1]
      deletion_pos=which(all_germ_located %!in%all_rep_located)
      germ_located=all_germ_located[deletion_pos]


      for (pos in germ_located) {

        curr_pos=0
        previous_dels=length(which(all_rep_located<pos)) ############ THINK RODRIGO THINK

        is_this_on=T
        reg_number=1
        while(is_this_on) {
          region=regions[reg_number]
          old_curr_pos=curr_pos
          curr_pos=curr_pos+nchar(IMGT_parsed_index[deletion_seq, region])
          if(pos<=(curr_pos+previous_dels) && pos >= (old_curr_pos+previous_dels)){
            # print("Y")
            # break
            tmp_del=strsplit(IMGT_parsed_index[deletion_seq, region], split="")[[1]]
            tmp_del[pos-old_curr_pos-previous_dels]=tolower(tmp_del[pos-old_curr_pos-previous_dels])
            IMGT_parsed_index[deletion_seq, region]=paste(tmp_del, collapse="")
            is_this_on=F
          }
          reg_number=reg_number+1

        }
        is_this_on=T

      }
    }

    for (deletion_seq in index_deletions_rep) {
      all_germ_located=str_locate_all(gsub(".","!",AIRR_input$germline_alignment[deletion_seq], fixed=T),"!")[[1]][,1]
      all_rep_located= str_locate_all(gsub(".","!",AIRR_input$sequence_alignment[deletion_seq], fixed=T),"!")[[1]][,1]
      deletion_pos=which(all_rep_located %!in% all_germ_located)
      germ_located=all_rep_located[deletion_pos]

      for (pos in germ_located) {

        curr_pos=0
        previous_dels=length(which(all_germ_located<pos))

        is_this_on=T
        reg_number=1
        while(is_this_on) {

          region=regions[reg_number]
          old_curr_pos=curr_pos
          curr_pos=curr_pos+nchar(reconstructed_IMGT_parsed_index[deletion_seq, region])
          if(pos<=(curr_pos+previous_dels) && pos >= (old_curr_pos+previous_dels)){


            tmp_del=strsplit(reconstructed_IMGT_parsed_index[deletion_seq, region], split="")[[1]]
            tmp_del[pos-old_curr_pos-previous_dels]=tolower(tmp_del[pos-old_curr_pos-previous_dels])
            reconstructed_IMGT_parsed_index[deletion_seq, region]=paste(tmp_del, collapse="")
            is_this_on=F
          }
          reg_number=reg_number+1


        }
        is_this_on=T

      }
    }


    ### PARTIAL FWRs
    if (FWR1partial) {
      IMGT_parsed_index$NT_FWR1 <- ""
    }

    if (FWR4partial) {
      IMGT_parsed_index$NT_FWR4 <- ""
    }

    reconstructed_IMGT_parsed_index$ORF_begins=IMGT_parsed_index$ORF_begins

    ### PARTIAL FWRs
    if (FWR1partial) {
      reconstructed_IMGT_parsed_index$NT_FWR1 <- sapply(
        c(1:nrow(reconstructed_IMGT_parsed_index)),
        function(z) if (nchar(
          gsub(
            "-", "", gsub(".", "", IMGT_parsed_index$NT_FWR1[z], fixed = T),
            fixed = T
          )
        ) ==
        0) {
          ""
        } else {
          str_sub(
            gsub(
              "-", "", gsub(".", "", reconstructed_IMGT_parsed_index$NT_FWR1[z], fixed = T),
              fixed = T
            ),
            -nchar(
              gsub(
                "-", "", gsub(".", "", IMGT_parsed_index$NT_FWR1[z], fixed = T),
                fixed = T
              )
            ),
            -1
          )
        }
      )
      reconstructed_IMGT_parsed_index$NT_FWR1 <- ""

      tmp_lengths_pre=sapply(reconstructed_IMGT_parsed_index$NT_CDR1, function(x) nchar(gsub(
        "-", "", gsub(".", "", x, fixed = T),
        fixed = T
      )))

      tmp_ref=reconstructed_IMGT_parsed_index$NT_CDR1
      reconstructed_IMGT_parsed_index$NT_CDR1 <- sapply(
        c(1:nrow(reconstructed_IMGT_parsed_index)),
        function(z) if (nchar(
          gsub(
            "-", "", gsub(".", "", IMGT_parsed_index$NT_CDR1[z], fixed = T),
            fixed = T
          )
        ) ==
        0) {
          ""
        } else {
          str_sub(
            gsub(
              "-", "", gsub(".", "", reconstructed_IMGT_parsed_index$NT_CDR1[z], fixed = T),
              fixed = T
            ),
            -nchar(
              gsub(
                "-", "", gsub(".", "", IMGT_parsed_index$NT_CDR1[z], fixed = T),
                fixed = T
              )
            )-length(str_match_all(sub("^.*?([A-Z])", "\\1", reconstructed_IMGT_parsed_index$NT_CDR1[z]),
                                   "[a-z]")[[1]])+
              length(str_match_all(IMGT_parsed_index$NT_CDR1[z],
                                   "[a-z]")[[1]]),
            -1
          )
        }
      )

      tmp_lengths_pro=sapply(reconstructed_IMGT_parsed_index$NT_CDR1, function(x) nchar(x))

      tmp_ref_orf=reconstructed_IMGT_parsed_index$ORF_begins
      # reconstructed_IMGT_parsed_index$ORF_begins=(reconstructed_IMGT_parsed_index$ORF_begins+(tmp_lengths_pre-tmp_lengths_pro)%%3)%%3
    }

    if (FWR4partial) {
      reconstructed_IMGT_parsed_index$NT_FWR4 <- sapply(
        c(1:nrow(reconstructed_IMGT_parsed_index)),
        function(z) if (nchar(
          gsub(
            "-", "", gsub(".", "", IMGT_parsed_index$NT_FWR4[z], fixed = T),
            fixed = T
          )
        ) ==
        0) {
          ""
        } else {
          str_sub(
            gsub(
              "-", "", gsub(".", "", reconstructed_IMGT_parsed_index$NT_FWR4[z], fixed = T),
              fixed = T
            ),
            1, nchar(
              gsub(
                "-", "", gsub(".", "", IMGT_parsed_index$NT_FWR4[z], fixed = T),
                fixed = T
              )
            )
          )
        }
      )
      reconstructed_IMGT_parsed_index$NT_FWR4 <- ""

      reconstructed_IMGT_parsed_index$NT_CDR3 <- sapply(
        c(1:nrow(reconstructed_IMGT_parsed_index)),
        function(z) if (nchar(
          gsub(
            "-", "", gsub(".", "", IMGT_parsed_index$NT_CDR3[z], fixed = T),
            fixed = T
          )
        ) ==
        0) {
          ""
        } else {
          str_sub(
            gsub(
              "-", "", gsub(".", "", reconstructed_IMGT_parsed_index$NT_CDR3[z], fixed = T),
              fixed = T
            ),
            1, nchar(
              gsub(
                "-", "", gsub(".", "", IMGT_parsed_index$NT_CDR3[z], fixed = T),
                fixed = T
              )
            )
          )
        }
      )
    }

    for (index_ins in index_insertions) {


        tmp_alignment_object <- pwalign::pairwiseAlignment(
            pattern = toupper(AIRR_input[index_ins, ]$sequence),
            subject = toupper(
              gsub( "-", "",
                    gsub(".", "",
                         AIRR_input[index_ins, ]$sequence_alignment,
                         fixed = T),
                    fixed = T)
            )
        )
        # tmp_alignment_seq=Biostrings::pairwiseAlignment(pattern=AIRR_input[index_ins,]$sequence,
        # subject=sequence) ##no es fiable

        tmp_align <- list(pattern = as.character(), subject = as.character())
        tmp_align$pattern <- as.character(pwalign::alignedPattern(tmp_alignment_object))
        tmp_align$subject <- as.character(pwalign::alignedSubject(tmp_alignment_object))

        gap_1 <- stringr::str_locate_all(
            as.character(tmp_align$pattern),
            "-"
        )[[1]][,
            1]
        gap_2 <- stringr::str_locate_all(
            as.character(tmp_align$subject),
            "-"
        )[[1]][,
            1]
        splitted_gap_1 <- split(
            gap_1, cumsum(
                c(
                  1, diff(gap_1) !=
                    1
              )
            )
        )
        if(length(splitted_gap_1)>0){

        }
        splitted_gap_2 <- split(
            gap_2, cumsum(
                c(
                  1, diff(gap_2) !=
                    1
              )
            )
        )

        num_gap <- 0
        diff_with_unfiltered_presequence=0
        for (gap in splitted_gap_2) {
            init <- T
            finit <- T
            num_gap <- num_gap + 1


            if (C_region_included && num_gap == length(splitted_gap_2) &&
                (gap[length(gap)] == nchar(as.character(tmp_align$subject)) ||
                 startsWith(paste(strsplit(tmp_align$pattern, split="")[[1]][splitted_gap_2[[num_gap]]],
                                  collapse=""),
                            paste(strsplit(tmp_align$subject, split="")[[1]][(splitted_gap_2[[num_gap]][1]+1):nchar(tmp_align$subject)],
                                  collapse=""))
                 )

                )
              {

            } else if (C_region_included && num_gap == 1 &&
                       (gap[1]==1 ||
                        endsWith(paste(strsplit(tmp_align$pattern, split="")[[1]][splitted_gap_2[[1]]], collapse=""),
                                 paste(strsplit(tmp_align$subject, split="")[[1]][1:(splitted_gap_2[[1]][1]-1)], collapse=""))
                        ))
              {
                diff_with_unfiltered_presequence=nchar(paste(strsplit(tmp_align$pattern, split="")[[1]][splitted_gap_2[[1]]], collapse=""))
              } else if (FWR1partial && num_gap == 1 && any(gap == 1)) {

            } else {
                for (region in regions) {
                  tmp_length <- sum(
                    nchar(
                      gsub(
                        "-", "", gsub(
                          ".", "", IMGT_parsed_index[index_ins, regions[c(1:(which(regions == region)))]],
                          fixed = T
                      ),
                        fixed = T
                    )
                  )
                ) +
                    1
                  sequence <- gsub(
                    "-", "", gsub(".", "", IMGT_parsed_index[index_ins, region], fixed = T),
                    fixed = T
                )


                  # if (!init && finit){ pos_init=0 }



                  if (gap[1] <= (tmp_length+diff_with_unfiltered_presequence) && init) {
                    if (region == "NT_FWR1") {
                      pos_init <- gap[1] - 1
                    } else {
                      pos_init <- gap[1] - sum(
                        nchar(
                          gsub(
                            "-", "", gsub(
                              ".", "", IMGT_parsed_index[index_ins, regions[c(
                                1:(which(regions == region) -
                                  1)
                            )]],
                              fixed = T
                          ),
                            fixed = T
                        )
                      )
                    ) -
                        1-diff_with_unfiltered_presequence
                    }

                    what_to_insert <- strsplit(tmp_align$pattern, split = "")[[1]][gap[1]:gap[length(gap)]]
                    what_to_insert <- tolower(what_to_insert)
                    init <- F

                    test <- strsplit(sequence, split = "")[[1]]
                    test <- append(test, what_to_insert, pos_init)

                    IMGT_parsed_index[index_ins, region] <- paste(test, collapse = "")
                  }

                  # if(gap[length(gap)]<=tmp_length && finit) { finit=F }

                }
            }



        }

    }


    IMGT_parsed_index <- rbind(IMGT_parsed_index, reconstructed_IMGT_parsed_index)

    IMGT_parsed_index <- IMGT_parsed_index[with(IMGT_parsed_index,
                                                order(ID, Sequence_type)),
                                           ]



    IMGT_parsed_index$Repertoire <- repertoire
    if(!is_example){
      suppressWarnings(
        write.table(
          IMGT_parsed_index, file = paste(output_path, "IMGT_parsed_index.txt", sep = ""),
          append = T, row.names = F,
          col.names = !file.exists(paste(output_path, "IMGT_parsed_index.txt", sep = "")),
          sep = "\t", quote = F
        )
      )
    } else {
      return(IMGT_parsed_index[c(1:8),])
    }


}
