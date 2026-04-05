#' 2_Feature_determination_1
#'
#' @description A fct function
#'
#' @return The return value, if any, from executing the function.
#' @import dplyr
#' @importFrom Biostrings subseq GENETIC_CODE
#' @keywords internal
#' @export
Feature_1 <- function(path_base, grouping_by) {
    input_path <- paste(path_base, "/1.Files_parsed", sep = "")
    output_path <- paste(path_base, "/2.Feature_determination", sep = "")

    ######## 2. Preparing IMGT index and adding V and J

    IMGT_parsed_index <- data.table::fread(
        file.path(input_path, "IMGT_parsed_index.txt"),
        header = T, na.strings=NULL, sep = "\t"
    )
    IMGT_parsed_index[is.na(IMGT_parsed_index)] = "" ## in case NA strings doesnt work

    prefix_gene <- strsplit(IMGT_parsed_index$V_region[1], split = " ")[[1]][1]
    if (prefix_gene == IMGT_parsed_index$V_region[1]) {
        prefix_gene <- ""
    } else {
        prefix_gene <- paste0(prefix_gene, " ")
    }
    if ("D_region" %in% colnames(IMGT_parsed_index)) {
        IMGT_parsed_index <- IMGT_parsed_index %>%
            dplyr::mutate(Best_V = gsub(prefix_gene, "", .data$V_region)) %>%
            dplyr::mutate(Best_J = gsub(prefix_gene, "", .data$J_region)) %>%
            dplyr::mutate(Best_D = gsub(prefix_gene, "", .data$D_region)) %>%
            # dplyr::mutate(Best_V=gsub(' [[:graph:][:space:]]{1, }$','',Best_V)) %>%
            # What is this for ARHJ mutate(Best_J=gsub('
            # [[:graph:][:space:]]{1, }$','',Best_J)) %>% mutate(Best_D=gsub('
            # [[:graph:][:space:]]{1, }$','',Best_D)) %>%
        dplyr::mutate(Best_V = gsub(",.*", "", .data$Best_V)) %>%
            dplyr::mutate(Best_J = gsub(",.*", "", .data$Best_J)) %>%
            dplyr::mutate(Best_D = gsub(",.*", "", .data$Best_D)) %>%
            dplyr::mutate(V_and_J = paste0(.data$Best_V, "_", .data$Best_J)) %>%
            dplyr::mutate(V_and_D_and_J = paste0(.data$Best_V, "_", .data$Best_D, "_", .data$Best_J))

        IMGT_parsed_index$D_region[which(is.na(IMGT_parsed_index$D_region))] <- "NotAvailable"
        IMGT_parsed_index$D_region[which(IMGT_parsed_index$D_region == "")] <- "NotAvailable"
    } else {
        IMGT_parsed_index <- IMGT_parsed_index %>%
            dplyr::mutate(Best_V = gsub(prefix_gene, "", .data$V_region)) %>%
            dplyr::mutate(Best_J = gsub(prefix_gene, "", .data$J_region)) %>%
            dplyr::mutate(Best_V = gsub(",.*", "", .data$Best_V)) %>%
            dplyr::mutate(Best_J = gsub(",.*", "", .data$Best_J)) %>%
            # mutate(Best_V=gsub(' [[:graph:][:space:]]{1, }$','',Best_V)) %>%
            # mutate(Best_J=gsub(' [[:graph:][:space:]]{1, }$','',Best_J)) %>%
        dplyr::mutate(V_and_J = paste0(.data$Best_V, "_", .data$Best_J))
    }


    ##print("Constructing Gene Usage Table ...")  ###absolute and relative info ########## maybe its better for the Shiny part, also the Best_ thing? not really because we already decided the best V


    # colnames(IMGT_parsed_index)[which(colnames(IMGT_parsed_index)
    # %in%regions)]=NT_regions
    write.table(
        as.data.frame(IMGT_parsed_index),
        file.path(output_path, "IMGT_parsed_index_extended.txt"),
        sep = "\t", row.names = FALSE, quote = FALSE
    )

    return(
        length(
            split(
              data.table::as.data.table(IMGT_parsed_index),
                by = grouping_by
            )
        )
    )
}


#' 2_Feature_determination_2
#'
#' @description A fct function
#'
#' @return The return value, if any, from executing the function.
#' @import parallel
#' @import foreach
#' @import bigstatsr
#' @import bigparallelr
#' @import stringr
#' @importFrom data.table :=
#' @keywords internal
#' @export
Feature__dataset <- function(
    path_base, DF_to_parse, name_DF_to_parse, FWR1partial, FWR4partial, standard = T
) {
  DF_to_parse[is.na(DF_to_parse)] = "" ## in case NA strings doesnt work

  # To avoid issues
  Sequence_type <- ID <- AA_Whole <- ORFs_match <- ORF <- NULL
  subsetDF <- d_rep <- d_recgerm <- NULL

    if (standard) {
        input_path <- paste(path_base, "/1.Files_parsed", sep = "")
        output_path <- paste(path_base, "/2.Feature_determination", sep = "")
    } else {
        input_path <- path_base
        output_path <- path_base
    }

    regions <- c("FWR1", "CDR1", "FWR2", "CDR2", "FWR3", "CDR3", "FWR4")
    NT_regions <- paste0("NT_", regions)
    AA_regions <- paste0("AA_", regions)

    ######## 3. Subdivide into smaller dfs according group and chain





    ######## 4. Iterate over each smaller df to generate DF

    ### Generating NT_Whole and removing sequences with missing regions (< 2
    ### nts) and strange characters other than std nts

    "%!in%" <- function(x, y) !(x %in% y)



    # DF_to_parse$NT_Whole= apply( DF_to_parse[ , ..NT_regions ] , 1 , paste ,
    # collapse = '' )

    if (FWR1partial) {
        init <- 2
    } else {
        init <- 1
    }

    if (FWR4partial) {
        fini <- length(NT_regions) -
            1
    } else {
        fini <- length(NT_regions)
    }

    tmp_NT_regions <- NT_regions[init:fini]

    DF_to_parse$NT_Whole <- do.call(paste, c(DF_to_parse[, tmp_NT_regions, with = FALSE], sep = ""))
    if (any(duplicated(DF_to_parse[Sequence_type != "Reconstructed_germline"]$ID))) {

        DF_to_parse=DF_to_parse[-(which(duplicated(DF_to_parse[Sequence_type != "Reconstructed_germline"]$ID))),]
    }

    ### Check absent regions (or less than 2 nts so no AA) and strange
    ### characters (there should not be strange chars)

    seqs_with_strange_chars <- which(
        unname(
            sapply(
                toupper(DF_to_parse$NT_Whole), function(x) if (any(
                  unique(strsplit(x, split = "")[[1]]) %!in%
                    c("A", "G", "C", "T")
              )) {
                  T
                } else {
                  F
                }
            )
        )
    )



    seqs_with_missing_regions <- c()
    nt_cutoff_length <- 2
    ignore_FWR1 <- c()
    ignore_FWR4 <- c()
    for (NT_region in NT_regions) {
        if (FWR1partial && NT_region == NT_regions[1]) {
            ignore_FWR1 <- c(
                ignore_FWR1, which(
                  nchar(unlist(DF_to_parse[, NT_region, with = FALSE])) <=
                    nt_cutoff_length
              )
            )
        } else if (FWR4partial && NT_region == NT_regions[7]) {
            ignore_FWR4 <- c(
                ignore_FWR4, which(
                  nchar(unlist(DF_to_parse[, NT_region, with = FALSE])) <=
                    nt_cutoff_length
              )
            )
        } else {
            seqs_with_missing_regions <- c(
                seqs_with_missing_regions, which(
                  nchar(unlist(DF_to_parse[, NT_region, with = FALSE])) <=
                    nt_cutoff_length
              )
            )

        }
    }

    to_remove <- unique(c(seqs_with_strange_chars, seqs_with_missing_regions))
    if (length(to_remove) >
        0) {
        IDs_to_remove <- DF_to_parse$ID[to_remove]
        message(
            paste0(
                "----Removed ", length(DF_to_parse$ID[to_remove]) *
                  2, " because of ", length(unique(seqs_with_missing_regions)),
                " sequences with short/missing regions and ", length(unique(seqs_with_strange_chars)),
                " with wrong chars"
            )
        )

        if (length(to_remove) >
            0) {
          data.table::fwrite(
                DF_to_parse[to_remove, ], file = paste0(output_path, "/", name_DF_to_parse, "_DISCARDED_bc_length-chars.txt"),
                sep = "\t"
            )
        }
        DF_to_parse <- DF_to_parse[!(ID %in% IDs_to_remove)]
    } else {
      message("----No removed sequences")
    }






    ## Let's try with an iterator
    cl <- parallel::makePSOCKcluster(
        max(
            1, detectCores(logical = FALSE) -
                4
        ),
        outfile = paste(output_path, "/", name_DF_to_parse, "_ORF.outfile.txt", sep = "")
    )
    doParallel::registerDoParallel(cl)

    i_DF_to_parse <- iterors::iteror(DF_to_parse, by = "row")

    parallel::clusterExport(
        cl = cl, unclass(
            lsf.str(
                envir = asNamespace("AbSolution"),
                all = T
            )
        ),
        envir = as.environment(asNamespace("AbSolution"))
    )
    AA_output <- foreach(
        subsetDF = i_DF_to_parse, .export = c(
            # "NT_regions",
            "nt_changes", "subseq_func", "str_length", "translate_fun"
        ),
        .packages = c("data.table", "Biostrings"),
        .verbose = F
    ) %dopar%
        {
            as.data.frame(
                t(
                  as.data.frame(
                    (unlist(
                      subseq_func(
                        as.character(subsetDF$NT_Whole),
                        IRanges::width(
                          sapply(
                            as.character(subsetDF[, NT_regions, with = FALSE]),
                            function(z) if (z ==
                              "NA" || is.na(z)) {
                              ""
                            } else {
                              z
                            }
                        )
                      ),
                        FWR1partial, FWR4partial, subsetDF$ORF_begins
                    )
                  ))
                )
              )
            )
        }
    parallel::stopCluster(cl)
    gc()
    # AA_output=lapply(AA_output, function(z) as.data.frame(z))
    AA_output <- data.table::rbindlist(AA_output)
    colnames(AA_output) <- c(
        paste0("AA_", c(regions, "Whole")),
        "ORF"
    )

    # write.table(data.frame(A=123), file=file.path(output_path,
    # paste0(name_DF_to_parse, '.sasdadasdexample')), sep='\t', quote =
    # F,row.names =F,col.names =T )

    # start.time <- Sys.time() for (j in c(1:nrow(DF_to_parse))) { if(j == 1) {
    # test=unlist(subseq_func(as.character(DF_to_parse[j, 'NT_Whole']),
    # width(sapply(as.character(DF_to_parse[j,..NT_regions]), function(z)
    # if(z=='NA'){''}else{z})),FWR1partial,FWR4partial )) } else {
    # test=rbind(test, unlist(subseq_func(as.character(DF_to_parse[j,
    # 'NT_Whole']), width(sapply(as.character(DF_to_parse[j,..NT_regions]),
    # function(z) if(z=='NA'){''}else{z})),FWR1partial,FWR4partial ))) } }
    # end.time <- Sys.time() time.taken_Linear <- end.time - start.time
    # print(time.taken_Linear)

    DF_to_parse <- cbind(DF_to_parse, AA_output)
    #print("Removing no ORF")
    message(
        paste0(
            "Removed ", length(unique(DF_to_parse[AA_Whole == "NO_GOOD_ORF"]$ID)) *
                2, " sequences because ", length(which(DF_to_parse$AA_Whole == "NO_GOOD_ORF")),
            " (", nrow(
                DF_to_parse[AA_Whole == "NO_GOOD_ORF"][Sequence_type != "Reconstructed_germline"]
            ),
            " from repertoire)", " had no good ORF"
        )
    )
    if (length(which(DF_to_parse$AA_Whole == "NO_GOOD_ORF")) >
        0) {
      data.table::fwrite(
            DF_to_parse[AA_Whole == "NO_GOOD_ORF"], file = paste0(output_path, "/", name_DF_to_parse, "_DISCARDED_bc_ORF.txt"),
            sep = "\t"
        )
    }
    DF_to_parse <- DF_to_parse[which(!(DF_to_parse$ID %in% DF_to_parse[AA_Whole == "NO_GOOD_ORF"]$ID))]

    # Check if same ORF ##in current case ofc not same orf
    data.table::setDT(DF_to_parse)[,
        ORFs_match := c("Not Match", "Match")[(data.table::uniqueN(ORF) ==
            1) + 1], ID]

    DF_to_parse$ORF=as.numeric(DF_to_parse$ORF)
    ######## 5. Iterate over each FBM

    ### Calculate num of ab features that go into the FBM


    #print("Calculating and storing Mutations and  Levenshtein distance...")  ######### its 0 for reconstructed



    ID_tmp <- DF_to_parse$ID[1]

    tmp <- DF_to_parse[ID == ID_tmp]

    regions <- c("FWR1", "CDR1", "FWR2", "CDR2", "FWR3", "CDR3", "FWR4", "Whole")
    NT_regions <- paste0("NT_", regions)
    AA_regions <- paste0("AA_", regions)
    NGly_motif <- "N[^P]{1}[S|T]"
    names(NGly_motif) <- "NGly"
    NGly_motif_pos <- c(1)
    names(NGly_motif_pos) <- names(NGly_motif)
    hot_cold_mtfs <- edited_hot_cold_mtfs()
    mut_points_mtfs <- mut_points_hot_cold_mtfs(hot_cold_mtfs)



    test <- full_analysis(
        tmp, AA_regions, NT_regions, hot_cold_mtfs, mut_points_mtfs, NGly_motif,
        NGly_motif_pos, FWR1partial, FWR4partial
    )


    tmp_example <- (as.matrix(rbind(test[["Matrix"]][["Germ"]], test[["Matrix"]][["Rep"]])))
    colnames(tmp_example) <- names(test$Matrix$Germ)
    write.table(
        tmp_example, file = file.path(output_path, paste0(name_DF_to_parse, ".example")),
        sep = "\t", quote = F, row.names = F, col.names = T
    )

    # tmp_FBM=FBM(nrow=length(test$Matrix$Rep),ncol=nrow(DF_to_parse),
    # is_read_only = F)
    tmp_FBM <- FBM(
        ncol = length(test$Matrix$Rep),
        nrow = nrow(DF_to_parse),
        is_read_only = F, backingfile = file.path(output_path, name_DF_to_parse)
    )
    # block.size=2*(round(block_size(ncol(tmp_FBM), ncores = nb_cores())/2))

    #print("Feature calculation...")

    # end=0 step_size=20000 tmp_info_df=list() ###reduce DF_to_parse size while
    # ((end + 1) < ncol(tmp_FBM)) { start=end + 1 print(paste(start,
    # paste0(hour(Sys.time()), ':', minute(Sys.time())))) end = end + step_size
    # subset_DF_to_parse=DF_to_parse[c(start:min(end, ncol(tmp_FBM))),]
    # block.size=2*(round(block_size(ncol(subset_DF_to_parse), ncores =
    # nb_cores())/2)) block.size=2*(round(block_size(ncol(DF_to_parse), ncores
    # = nb_cores())/2)) tmp_tmp_info_df= big_applyAb(tmp_FBM, function(X, ind,
    # full_analysis, DF_to_parse,ind_vector, name_DF_to_parse,output_path,
    # AA_regions,
    # NT_regions,hot_cold_mtfs,mut_points_mtfs,NGly_motif,NGly_motif_pos,
    # FWR1partial, FWR4partial) { tmp_df=data.frame(matrix(NA, nrow = nrow(X),
    # ncol = length(ind))) library(data.table) library(bigassertr)
    # source('Scripts/0_Global_shiny.R') position=0 for (ID_tmp in
    # unique(DF_to_parse$ID[ind])) { position=position+1
    # results=full_analysis(DF_to_parse[ID==ID_tmp], AA_regions,
    # NT_regions,hot_cold_mtfs,mut_points_mtfs,NGly_motif,NGly_motif_pos,
    # FWR1partial, FWR4partial) if(position == 1) {
    # results_df=data.frame(matrix(NA, nrow = length(ind), ncol =
    # length(results$DF$Germ))) } tmp_df[, position] <- results$Matrix$Germ
    # results_df[position, ] <- unname(results$DF$Germ) position=position+1
    # tmp_df[, position] <- results$Matrix$Rep results_df[position, ] <-
    # unname(results$DF$Rep) } X[, ind_vector[ind]] <- tmp_df if (1 %in%
    # ind_vector[ind]) { tmp_example=t(as.matrix(tmp_df[,1:2]))
    # colnames(tmp_example)=names(results$Matrix$Germ) write.table(tmp_example,
    # file=file.path(output_path, paste0(name_DF_to_parse, '.example')),
    # sep='\t', quote = F,row.names =F,col.names =T ) } # if (ind %% 1000 == 0)
    # { # fwrite(paste(paste0(hour(Sys.time()), ':', minute(Sys.time())),ind,
    # sep=' '), file = paste0(output_path, '/', name_DF_to_parse,
    # '_feature_outfile.txt'), sep='\t') # }
    # colnames(results_df)=colnames(data.frame(as.list(results$DF$Germ)))
    # results_df }, full_analysis = full_analysis,
    # DF_to_parse=subset_DF_to_parse,
    # name_DF_to_parse=name_DF_to_parse,output_path=output_path,
    # AA_regions=AA_regions,
    # NT_regions=NT_regions,hot_cold_mtfs=hot_cold_mtfs,mut_points_mtfs=mut_points_mtfs,NGly_motif=NGly_motif,NGly_motif_pos=NGly_motif_pos,
    # FWR1partial=FWR1partial, FWR4partial=FWR4partial, a.combine = NULL,
    # ncores = max(1,nb_cores()-4), block.size=block.size,
    # ind_vector=c(start:min(end, ncol(tmp_FBM))), ind=
    # c(1:length(c(start:min(end, ncol(tmp_FBM))))) )
    # #min(ceiling((ncol(tmp_FBM)/block.size)), nb_cores()) print('done')
    # tmp_info_df[[length(tmp_info_df)+1]]=rbindlist(tmp_tmp_info_df)
    # print('rbindlist done') } info_df=rbindlist(tmp_info_df)
    # rm(tmp_tmp_info_df) rm(tmp_info_df) print('Feature calculation... Done')


    cl <- parallel::makePSOCKcluster(
        max(1, detectCores() - 5),
        outfile = paste(output_path, "/", name_DF_to_parse, "_feature_outfile.txt", sep = "")
    )
    doParallel::registerDoParallel(cl)
    parallel::clusterExport(
        cl = cl, unclass(
            lsf.str(
                envir = asNamespace("AbSolution"),
                all = T
            )
        ),
        envir = as.environment(asNamespace("AbSolution"))
    )
    ### perhaps use flock package better to update on-the-go the FBM but still
    ### I would need the DF back
    end <- 0
    step_size <- 2000

    tmp_info_df <- list()
    ### reduce DF_to_parse size
    while ((end + 1) < nrow(tmp_FBM)) {
        start <- end + 1
        #print(
        #     paste(
        #         start, paste0(
        #           hour(Sys.time()),
        #           ":", minute(Sys.time())
        #       )
        #     )
        # )
        end <- end + step_size

        subset_DF_to_parse <- DF_to_parse[c(start:min(end, nrow(tmp_FBM))),
            ]

        i_rep <- iterors::iteror(
            subset_DF_to_parse[seq(
                from = 1, to = nrow(subset_DF_to_parse),
                by = 2
            ),
                ], by = "row"
        )
        i_recgerm <- iterors::iteror(
            subset_DF_to_parse[seq(
                from = 2, to = nrow(subset_DF_to_parse),
                by = 2
            ),
                ], by = "row"
        )

        results <- foreach(
            d_rep = i_rep, d_recgerm = i_recgerm, .export = c("nt_changes", "lv_changes", "length_of_region",
                "count_elements", "locate_motif", "group_freqs", "alakazam_properties",
                "Peptides_properties", "compare_sequences", "instaIndex_improved",
                "full_analysis",
                # "NT_regions", "AA_regions", "hot_cold_mtfs", "mut_points_mtfs",
                # "NGly_motif", "NGly_motif_pos", "FWR1partial", "FWR4partial",
                "stringdist"
            ),
            .packages = c("alakazam", "Peptides", "Biostrings", "BiocGenerics", "data.table", "stringr"),
            .verbose = F
        ) %dopar%
            {
                full_analysis(
                  data.table::as.data.table(rbind(d_rep, d_recgerm)),
                  AA_regions, NT_regions, hot_cold_mtfs, mut_points_mtfs, NGly_motif,
                  NGly_motif_pos, FWR1partial, FWR4partial
              )
            }



        tmp_tmp_info_df <- do.call(
            rbind, lapply(results, function(z) rbind(z[["DF"]][["Germ"]], z[["DF"]][["Rep"]]))
        )

        tmp_info_df[[length(tmp_info_df) +
            1]] <- as.data.frame(tmp_tmp_info_df)

        results_df <- do.call(
            rbind, lapply(
                results, function(z) rbind(
                  unname(z[["Matrix"]][["Germ"]]),
                  unname(z[["Matrix"]][["Rep"]])
              )
            )
        )


        # tmp_FBM[,c(start:min(end, ncol(tmp_FBM)))] = results_df
        tmp_FBM[c(start:min(end, nrow(tmp_FBM))),
            ] <- results_df
        rm(results_df)
    }

    parallel::stopCluster(cl)

    info_df <- data.table::rbindlist(tmp_info_df)

    info_df <- cbind(DF_to_parse, info_df)
    write.table(
        info_df, file = file.path(output_path, paste0(name_DF_to_parse, ".info")),
        sep = "\t", quote = F, row.names = F, col.names = T
    )
    rm(info_df)
    gc()

    # def_FBM <-big_transpose(tmp_FBM, backingfile = file.path(output_path,
    # name_DF_to_parse)) def_FBM$save()
    tmp_FBM$save()
}


#' 2_Feature_determination
#'
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#'
#' @noRd
full_analysis <- function(
    tmp, AA_regions, NT_regions, hot_cold_mtfs, mut_points_mtfs, NGly_motif, NGly_motif_pos,
    FWR1partial, FWR4partial
) {
    germline_index <- which(tmp$Sequence_type == "Reconstructed_germline")
    repertoire_index <- which(tmp$Sequence_type != "Reconstructed_germline")

    if (length(c(germline_index, repertoire_index)) !=
        2) {
        #print(paste0("[ERROR] More sequences than foreseen with ID ", unique(tmp$ID)))
    }

    list_row_values <- list(Rep = c(), Germ = c())
    list_row_values_df <- list(Rep = c(), Germ = c())
    list_lengths <- list(Rep_NT = c(), Rep_AA = c(), Germ_NT = c(), Germ_AA = c())
    values_not_to_diff <- c()

    for (index in c(repertoire_index, germline_index)) {
        row_values <- c()
        row_values_df <- c()
        length_regions_index <- c()
        for (NT_region in NT_regions) {
            if ((NT_region == "NT_FWR1" && FWR1partial) || (NT_region == "NT_FWR4" &&
                FWR4partial)) {
                orig_NT_region <- NT_region
                NT_region <- "NT_FWR2"
                seq_rec.germline <- (as.character(tmp[germline_index, NT_region, with = FALSE]))
                seq_repertoire <- (as.character(tmp[repertoire_index, NT_region, with = FALSE]))
                sequence <- (as.character(tmp[index, NT_region, with = FALSE]))
                if (index == repertoire_index) {
                  values_not_to_diff <- c(
                    values_not_to_diff, nt_changes("", "", orig_NT_region)
                )
                }
                seq_rec.germline <- toupper(seq_rec.germline)
                seq_repertoire <- toupper(seq_repertoire)
                sequence <- toupper(sequence)
                tmp_partial <- length_of_region(sequence, NT_region)
                tmp_partial_value <- rep(0, length(tmp_partial))
                names(tmp_partial_value) <- gsub(NT_region, orig_NT_region, names(tmp_partial))
                length_regions_index <- c(length_regions_index, tmp_partial_value)

                tmp_partial <- count_elements(sequence, NT_region)
                tmp_partial_value <- rep(0, length(tmp_partial))
                names(tmp_partial_value) <- gsub(NT_region, orig_NT_region, names(tmp_partial))
                count_el <- tmp_partial_value

            } else {
              seq_rec.germline <- (as.character(tmp[germline_index, NT_region, with = FALSE]))
              seq_repertoire <- (as.character(tmp[repertoire_index, NT_region, with = FALSE]))
              sequence <- (as.character(tmp[index, NT_region, with = FALSE]))
                if (index == repertoire_index) {
                  values_not_to_diff <- c(values_not_to_diff, nt_changes(sequence, seq_rec.germline, NT_region))
                }
                seq_rec.germline <- toupper(seq_rec.germline)
                seq_repertoire <- toupper(seq_repertoire)
                sequence <- toupper(sequence)
                length_regions_index <- c(length_regions_index, length_of_region(sequence, NT_region))
                count_el <- count_elements(sequence, NT_region)
            }


            row_values <- c(row_values, length_regions_index[length(length_regions_index)])
            row_values <- c(row_values, count_el)
        }
        if (index == repertoire_index) {
            list_lengths$Rep_NT <- length_regions_index
        } else {
            list_lengths$Germ_NT <- length_regions_index
        }


        sequence <- toupper(as.character(tmp[index, "NT_Whole"]))
        for (motif_i in c(1:length(hot_cold_mtfs))) {
            tmp_results <- locate_motif(
                sequence, length_regions_index, hot_cold_mtfs[motif_i], names(hot_cold_mtfs[motif_i]),
                motif_ini = mut_points_mtfs[which(
                  names(mut_points_mtfs) ==
                    names(hot_cold_mtfs[motif_i])
              )],
                NT_regions
            )
            row_values_df <- c(row_values_df, tmp_results$index_locis)
            row_values <- c(row_values, tmp_results$motif_count, tmp_results$perc_locis)
        }


        length_regions_index <- c()
        for (AA_region in AA_regions) {

            if ((AA_region == "AA_FWR1" && FWR1partial) || (AA_region == "AA_FWR4" &&
                FWR4partial)) {
                orig_AA_region <- AA_region
                AA_region <- "AA_FWR2"
                seq_rec.germline <- as.character(tmp[germline_index, AA_region, with = FALSE ])
                seq_repertoire <- as.character(tmp[repertoire_index, AA_region, with = FALSE ])
                sequence <- as.character(tmp[index, AA_region, with = FALSE ])
                # if (index == repertoire_index) {
                #   values_not_to_diff <- c(values_not_to_diff, aa_changes(seq_repertoire, seq_rec.germline, AA_region))
                # }

                tmp_partial <- length_of_region(sequence, AA_region)
                tmp_partial_value <- rep(0, length(tmp_partial))
                names(tmp_partial_value) <- gsub(AA_region, orig_AA_region, names(tmp_partial))
                length_regions_index <- c(length_regions_index, tmp_partial_value)
                row_values <- c(row_values, length_regions_index[length(length_regions_index)])

                tmp_partial <- c(
                  count_elements(sequence, AA_region),
                  group_freqs(sequence, AA_region),
                  alakazam_properties(sequence, AA_region),
                  Peptides_properties(sequence, AA_region)
              )
                tmp_partial_value <- rep(0, length(tmp_partial))
                names(tmp_partial_value) <- gsub(AA_region, orig_AA_region, names(tmp_partial))

                row_values <- c(row_values, tmp_partial_value)
            } else {
                seq_rec.germline <- as.character(tmp[germline_index, AA_region, with = FALSE ])
                seq_repertoire <- as.character(tmp[repertoire_index, AA_region, with = FALSE ])
                sequence <- as.character(tmp[index, AA_region, with = FALSE ])
                # if (index == repertoire_index) {
                #   values_not_to_diff <- c(values_not_to_diff, aa_changes(seq_repertoire, seq_rec.germline, AA_region))
                # }
                length_regions_index <- c(length_regions_index, length_of_region(sequence, AA_region))
                row_values <- c(row_values, length_regions_index[length(length_regions_index)])
                row_values <- c(row_values, count_elements(sequence, AA_region))
                row_values <- c(row_values, group_freqs(sequence, AA_region))
                row_values <- c(row_values, alakazam_properties(sequence, AA_region))
                row_values <- c(row_values, Peptides_properties(sequence, AA_region))
            }

        }
        if (index == repertoire_index) {
            list_lengths$Rep_AA <- length_regions_index
        } else {
            list_lengths$Germ_AA <- length_regions_index
        }

        sequence <- as.character(tmp[index, "AA_Whole"])
        tmp_results <- locate_motif(
            sequence, length_regions_index, NGly_motif, names(NGly_motif),
            motif_ini = NGly_motif_pos, AA_regions
        )
        row_values_df <- c(row_values_df, tmp_results$index_locis)
        row_values <- c(row_values, tmp_results$motif_count, tmp_results$perc_locis)

        if ((AA_region == "AA_FWR1" && FWR1partial) || (AA_region == "NT_FWR4" &&
            FWR4partial)) {
            row_values[which(is.nan(row_values))] <- 0
        }

        if (index == repertoire_index) {
            list_row_values$Rep <- row_values
            list_row_values_df$Rep <- row_values_df
        } else {
            list_row_values$Germ <- row_values
            list_row_values_df$Germ <- row_values_df
        }
    }



    seq_aa_rec.germline <- as.character(tmp[germline_index, "AA_Whole"])
    seq_aa_repertoire <- as.character(tmp[repertoire_index, "AA_Whole"])
    values_not_to_diff <- c(
        values_not_to_diff, compare_sequences(
            rep_seq = seq_aa_repertoire, germ_seq = seq_aa_rec.germline, mode = "AA",
            aa_rep_lengths = list_lengths$Rep_AA, aa_germ_lengths = list_lengths$Germ_AA,
            AA_regions = AA_regions,
            aant_rep_seq=as.character(tmp[repertoire_index, "NT_Whole"]),
            aant_germ_seq=as.character(tmp[germline_index, "NT_Whole"])
        )
    )

    seq_nt_rec.germline <- (as.character(tmp[germline_index, "NT_Whole"]))
    seq_nt_repertoire <- (as.character(tmp[repertoire_index, "NT_Whole"]))
    values_not_to_diff <- c(
        values_not_to_diff, compare_sequences(
            rep_seq = seq_nt_repertoire, germ_seq = seq_nt_rec.germline, mode = "NT",
            nt_rep_lengths = list_lengths$Rep_NT, nt_germ_lengths = list_lengths$Germ_NT,
            NT_regions = NT_regions
        )
    )

    seq_nt_rec.germline <- (as.character(tmp[germline_index, "NT_Whole"]))
    seq_nt_repertoire <- (as.character(tmp[repertoire_index, "NT_Whole"]))

    codon_analysis <- compare_sequences(
        rep_seq = seq_nt_repertoire, germ_seq = seq_nt_rec.germline, mode = "Codon",
        nt_rep_lengths = list_lengths$Rep_NT, nt_germ_lengths = list_lengths$Germ_NT,
        NT_regions = NT_regions, aa_rep_lengths = list_lengths$Rep_AA, aa_germ_lengths = list_lengths$Germ_AA,
        AA_regions = AA_regions, ORF_rep = tmp$ORF[repertoire_index], ORF_germ = tmp$ORF[germline_index]
    )

    values_not_to_diff <- c(values_not_to_diff, codon_analysis$counts)
    list_row_values$Rep <- c(list_row_values$Rep, codon_analysis$rep_codon_count)
    list_row_values$Germ <- c(list_row_values$Germ, codon_analysis$germ_codon_count)

    diff_tmp <- list_row_values$Rep - list_row_values$Germ
    names(diff_tmp) <- paste0(
        names(diff_tmp),
        "_Diff"
    )
    diff_tmp_ger <- rep(0, length(diff_tmp))
    names(diff_tmp_ger) <- names(diff_tmp)
    list_row_values$Rep <- c(list_row_values$Rep, diff_tmp)
    list_row_values$Germ <- c(list_row_values$Germ, diff_tmp_ger)



    list_row_values$Rep <- c(list_row_values$Rep, values_not_to_diff)
    empty_values_not_to_diff <- rep(0, length(values_not_to_diff))
    names(empty_values_not_to_diff) <- names(values_not_to_diff)
    list_row_values$Germ <- c(list_row_values$Germ, empty_values_not_to_diff)

    return(list(DF = list_row_values_df, Matrix = list_row_values))

}
