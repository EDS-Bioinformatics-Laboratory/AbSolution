#' 1_Parsing
#'
#' @description A fct function
#'
#' @return The return value, if any, from executing the function.
#' @import Biostrings
#' @import stringr
#'
#' @noRd
parse_AIRRSeq_file <- function(
    file, group, patient, subgroup, sample, input_path, C_region_included, FWR1partial,
    FWR4partial, D_gene, repertoire, output_path
) {


    ### Am I working with the insertions?
    print(file)
    AIRR_input <- read.table(
        file = paste(
            input_path, paste(file, ".tsv", sep = ""),
            sep = ""
        ),
        sep = "\t", header = TRUE, fill = T
    )

    # Removing Ns
    AIRR_input <- AIRR_input[which(!grepl("N", AIRR_input$sequence)),
        ]
    AIRR_input <- AIRR_input[which(!grepl("n", AIRR_input$sequence)),
        ]
    AIRR_input <- AIRR_input[which(AIRR_input$sequence_alignment != ""),
        ]
    regions <- c("NT_FWR1", "NT_CDR1", "NT_FWR2", "NT_CDR2", "NT_FWR3", "NT_CDR3", "NT_FWR4")

    print("Testing parse 1")
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
    print("Testing parse 2")
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
        print("Testing parse 2a")
        IMGT_parsed_index$Clone_ID <- AIRR_input$clone_id
    } else if ("raw_clonotype_id" %in% colnames(AIRR_input)) {
        print("Testing parse 2b")
        IMGT_parsed_index$Clone_ID <- AIRR_input$raw_clonotype_id
    } else if ("clone_id" %in% colnames(AIRR_input)) {
        print("Testing parse 2c")
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

    ### PARTIAL FWRs
    if (FWR1partial) {
        IMGT_parsed_index$NT_FWR1 <- ""
    }

    if (FWR4partial) {
        IMGT_parsed_index$NT_FWR4 <- ""
    }




    ########## R: insertions are not included here lalalala

    or_CDR1 <- 1 + nchar(AIRR_input$fwr1)
    or_FWR2 <- or_CDR1 + nchar(AIRR_input$cdr1)
    or_CDR2 <- or_FWR2 + nchar(AIRR_input$fwr2)
    or_FWR3 <- or_CDR2 + nchar(AIRR_input$cdr2)
    or_CDR3 <- or_FWR3 + nchar(AIRR_input$fwr3)
    or_FWR4 <- or_CDR3 + nchar(AIRR_input$cdr3)

    ### even if the format says aligned, not all programmes (e.g. IMGT,
    ### 02/2024) have these fields aligned. Let's add a checkup
    if (any(grepl(".", AIRR_input$fwr1, fixed = T)) ||
        any(grepl(".", AIRR_input$fwr2, fixed = T)) ||
        any(grepl(".", AIRR_input$fwr3, fixed = T)) ||
        any(grepl(".", AIRR_input$fwr4, fixed = T)) ||
        any(grepl(".", AIRR_input$cdr1, fixed = T)) ||
        any(grepl(".", AIRR_input$cdr2, fixed = T)) ||
        any(grepl(".", AIRR_input$cdr3, fixed = T))) {
        aligned_regions <- T
    } else {
        aligned_regions <- F
    }
    AIRR_input$germline_alignment_nospace <- gsub(".", "", AIRR_input$germline_alignment, fixed = T)

    if (aligned_regions) {
        NT_CDR3_germ <- strsplit(
            substring(AIRR_input$germline_alignment, or_CDR3, or_CDR3 - 1 + nchar(AIRR_input$cdr3)),
            split = ""
        )
    } else {
        NT_CDR3_germ <- strsplit(
            substring(
                AIRR_input$germline_alignment_nospace, or_CDR3, or_CDR3 - 1 + nchar(AIRR_input$cdr3)
            ),
            split = ""
        )
    }
    # NT_CDR3_germ=substring(AIRR_input$sequence_alignment,or_CDR3,or_CDR3-1+nchar(AIRR_input$cdr3))




    print("Testing parse 3")


    for (index_cdr3 in c(1:length(NT_CDR3_germ))) {
        tmp <- NT_CDR3_germ[[index_cdr3]]
        index_n <- c(
            which(tmp == "N"),
            which(tmp == "n")
        )
        if (D_gene == F) {
            index_n <- c(min(index_n):max(index_n))
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
                  gsub(".", "", AIRR_input$sequence_alignment[index_cdr3], fixed = T),
                  or_CDR3[index_cdr3], or_CDR3[index_cdr3] - 1 + nchar(AIRR_input$cdr3[index_cdr3])
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
                  ".", "", substring(AIRR_input$germline_alignment_nospace, 1, nchar(AIRR_input$fwr1)),
                  fixed = T
              ),
                fixed = T
            ),
            NT_CDR1 = gsub(
                "-", "", gsub(
                  ".", "", substring(
                    AIRR_input$germline_alignment_nospace, or_CDR1, or_CDR1 - 1 +
                      nchar(AIRR_input$cdr1)
                ),
                  fixed = T
              ),
                fixed = T
            ),
            NT_FWR2 = gsub(
                "-", "", gsub(
                  ".", "", substring(
                    AIRR_input$germline_alignment_nospace, or_FWR2, or_FWR2 - 1 +
                      nchar(AIRR_input$fwr2)
                ),
                  fixed = T
              ),
                fixed = T
            ),
            NT_CDR2 = gsub(
                "-", "", gsub(
                  ".", "", substring(
                    AIRR_input$germline_alignment_nospace, or_CDR2, or_CDR2 - 1 +
                      nchar(AIRR_input$cdr2)
                ),
                  fixed = T
              ),
                fixed = T
            ),
            NT_FWR3 = gsub(
                "-", "", gsub(
                  ".", "", substring(
                    AIRR_input$germline_alignment_nospace, or_FWR3, or_FWR3 - 1 +
                      nchar(AIRR_input$fwr3)
                ),
                  fixed = T
              ),
                fixed = T
            ),
            NT_CDR3 = unlist(NT_CDR3_germ),
            NT_FWR4 = gsub(
                "-", "", gsub(
                  ".", "", substring(
                    AIRR_input$germline_alignment_nospace, or_FWR4, or_FWR4 - 1 +
                      nchar(AIRR_input$fwr4)
                ),
                  fixed = T
              ),
                fixed = T
            ),
            stringsAsFactors = F
        )
    }


    print("Testing parse 4")
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

    print("Testing parse  5")
    print(length(index_insertions))
    for (index_ins in index_insertions) {


        tmp_alignment_object <- Biostrings::pairwiseAlignment(
            pattern = AIRR_input[index_ins, ]$sequence, subject = gsub(
                "-", "", gsub(".", "", AIRR_input[index_ins, ]$sequence_alignment, fixed = T),
                fixed = T
            )
        )
        # tmp_alignment_seq=Biostrings::pairwiseAlignment(pattern=AIRR_input[index_ins,]$sequence,
        # subject=sequence) ##no es fiable

        tmp_align <- list(pattern = as.character(), subject = as.character())
        tmp_align$pattern <- as.character(Biostrings::alignedPattern(tmp_alignment_object))
        tmp_align$subject <- as.character(Biostrings::alignedSubject(tmp_alignment_object))

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
        splitted_gap_2 <- split(
            gap_2, cumsum(
                c(
                  1, diff(gap_2) !=
                    1
              )
            )
        )

        num_gap <- 0
        for (gap in splitted_gap_2) {
            init <- T
            finit <- T
            num_gap <- num_gap + 1

            # if(FWR1partial && num_gap==1) { print('-----------')
            # print(index_ins) print(gap) print('-----------') }
            if (C_region_included && num_gap == length(splitted_gap_2) &&
                gap[length(gap)] ==
                  nchar(as.character(tmp_align$subject))) {

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



                  if (gap[1] <= tmp_length && init) {
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
                        1
                    }

                    what_to_insert <- strsplit(tmp_align$pattern, split = "")[[1]][gap[1]:gap[length(gap)]]
                    init <- F

                    test <- strsplit(sequence, split = "")[[1]]
                    test <- append(test, what_to_insert, pos_init)
                    print(region)
                    print(index_ins)
                    print(sequence)
                    print(paste(test, collapse = ""))
                    IMGT_parsed_index[index_ins, region] <- paste(test, collapse = "")
                  }

                  # if(gap[length(gap)]<=tmp_length && finit) { finit=F }

                }
            }



        }

    }


    IMGT_parsed_index <- rbind(IMGT_parsed_index, reconstructed_IMGT_parsed_index)

    IMGT_parsed_index <- IMGT_parsed_index[with(IMGT_parsed_index, order(ID, Sequence_type)),
        ]

    for (region in regions) {
        IMGT_parsed_index[, region] <- toupper(
            gsub(
                "-", "", gsub(".", "", IMGT_parsed_index[, region], fixed = T),
                fixed = T
            )
        )
    }

    IMGT_parsed_index$Repertoire <- repertoire
    write.table(
        IMGT_parsed_index, file = paste(output_path, "IMGT_parsed_index.txt", sep = ""),
        append = T, row.names = F, col.names = !file.exists(paste(output_path, "IMGT_parsed_index.txt", sep = "")),
        sep = "\t", quote = F
    )

}
