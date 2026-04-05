#' Merge FBMs and meta files
#'
#' @description A function to join by row individual FBMs and .info files
#'     into a unique FBM and a unique .info file, merged_bm.
#' @param filepath Path to the individual FBMs.
#' @return A list with the merged FBM, the merged info file and the FBM header.
#' @importFrom utils  glob2rx
#' @import bigstatsr
#' @import tools
#' @import bigparallelr
#' @keywords internal
#' @export
merge_FBMs <- function(filepath) {
    basenames <- basename(list.files(path = filepath, pattern = glob2rx("*.rds")))
    prefixes <- tools::file_path_sans_ext(basenames)

    list_descriptions <- list()
    list_FBMs <- list()
    list_headers <- list()
    for (prefix in unique(prefixes)[which(
        unique(prefixes) !=
            "Merged_bm"
    )]) {
        list_descriptions[[length(list_descriptions) +
            1]] <- data.table::fread(file.path(filepath, paste0(prefix, ".info")))
        list_headers[[length(list_headers) +
            1]] <- data.table::fread(file.path(filepath, paste0(prefix, ".example")))
        if (file.exists(file.path(filepath, paste0(prefix, ".rds")))) {
            list_FBMs[[length(list_FBMs) +
                1]] <- bigstatsr::big_attach(file.path(filepath, paste0(prefix, ".rds")))
        }

    }


    merged_df <- data.table::rbindlist(list_descriptions, fill=TRUE)
    merged_df[is.na(merged_df)]<- "Undetermined"

    if("C_region" %in% colnames(merged_df)){
      merged_df$C_region[which(merged_df$C_region=="")]="Undetermined"
    }

    if("Chain" %in% colnames(merged_df)){
      merged_df$Chain[which(merged_df$Chain=="")]="Undetermined"
    }
    # merged_df[merged_df == ""] <- "Undetermined"

    merged_header <- colnames(list_headers[[1]])

    for (i in c(1:length(list_headers))) {
        header_tmp <- colnames(list_headers[[i]])
        if (all(
            length(merged_header) ==
                length(header_tmp)
        ) &&
            all(merged_header == header_tmp)) {

        } else {
            #print("[ERROR] different number of columns and/or colnames")
        }

    }



    m <- ncol(list_FBMs[[1]])  # assuming that all have the same number of columns
    n <- sum(sapply(list_FBMs, nrow))

    if (file.exists(file.path(filepath, paste0("Merged_bm", ".bin")))) {
        file.remove(file.path(filepath, paste0("Merged_bm", ".bin")))
        if (file.exists(file.path(filepath, paste0("Merged_bm", ".desc")))) {
            file.remove(file.path(filepath, paste0("Merged_bm", ".desc")))
        }

    }

    merged_del <- T
    if (file.exists(file.path(filepath, paste0("Merged_bm", ".bk")))) {

        if (file.exists(file.path(filepath, paste0("Merged_bm", ".info")))) {
            tmp_merged_df <- data.table::fread(file.path(filepath, paste0("Merged_bm", ".info")))
            if (nrow(merged_df) ==
                nrow(tmp_merged_df) &
                all(tmp_merged_df[, 1] == merged_df[, 1]) &
                all(tmp_merged_df[, 2] == merged_df[, 2]) &
                all(tmp_merged_df[, 3] == merged_df[, 3])) {
                merged_del <- F
                merged_bm <- bigstatsr::big_attach(file.path(filepath, paste0("Merged_bm", ".rds")))
                merged_df <- tmp_merged_df
                if (nrow(merged_bm) !=
                  nrow(tmp_merged_df)) {
                  merged_del <- T
                }
            }
        }
    }


    if (merged_del == T) {
        if (file.exists(file.path(filepath, paste0("Merged_bm", ".bk")))) {
            file.remove(file.path(filepath, paste0("Merged_bm", ".bk")))
            if (file.exists(file.path(filepath, paste0("Merged_bm", ".rds")))) {
                file.remove(file.path(filepath, paste0("Merged_bm", ".rds")))
            }
        }

        merged_bm <- bigstatsr::FBM(
            n, m, backingfile = file.path(filepath, "Merged_bm"),
            is_read_only = F
        )
        data.table::fwrite(
            merged_df, file.path(filepath, paste0("Merged_bm", ".info")),
            sep = "\t"
        )


        ### Contribution by F. Privé in
        ### https://github.com/privefl/bigstatsr/issues/176
        ## NOTE THIS NEEDS TO SCALE UP BETTER, otherwise it occupies too much RAM
        # bigstatsr::big_apply(
        #     merged_bm, function(X, ind, list_fbm) {
        #         X[, ind] <- do.call(
        #           "rbind", lapply(
        #             list_fbm, function(fbm) fbm[,
        #               ind, drop = FALSE]
        #         )
        #       )
        #         gc()
        #         NULL
        #     }, list_fbm = list_FBMs, a.combine = "c", ncores = nb_cores()
        # )


        # bigstatsr::big_apply(
        #   merged_bm, function(X, ind, list_fbm) {
        #
        #     offset_FBM=0
        #     for (i in c(1:length(list_fbm))) {
        #       for(j in c(1:nrow(list_fbm[[i]]))) {
        #         X[j+offset_FBM, ind]=list_fbm[[i]][j,ind, drop = FALSE]
        #       }
        #       offset_FBM=offset_FBM+nrow(list_fbm[[i]])
        #     }
        #     NULL
        #   }, list_fbm = list_FBMs, a.combine = "c", ncores = nb_cores()
        # )


        #slow but RAM-friendly
        # offset_FBM=0
        # for (i in c(1:length(list_FBMs))) {
        #   print(paste(i,length(list_FBMs), sep="/"))
        #   for(j in c(1:nrow(list_FBMs[[i]]))) {
        #     merged_bm[j+offset_FBM, c(1:m)]=list_FBMs[[i]][j,c(1:m), drop = FALSE]
        #   }
        #   offset_FBM=offset_FBM+nrow(list_FBMs[[i]])
        # }



        bigstatsr::big_apply(
            merged_bm, function(X, ind, list_fbm) {
                X[, ind] <- do.call(
                  "rbind", lapply(
                    list_fbm, function(fbm) fbm[,
                      ind, drop = FALSE]
                )
              )
                gc()
                NULL
            }, list_fbm = list_FBMs, a.combine = "c", ncores = 1
        )

        merged_bm$save()

    }

    return(list(Header = merged_header, Short_DF = merged_df, Big_DF = merged_bm))

}

#' 3_Data_processing
#'
#' @description A fct function
#'
#' @return The return value, if any, from executing the function.
#' @import bigstatsr
#' @import bigparallelr
#' @keywords internal
#' @export
filter_merged <- function(
    FBM, merged_df, merged_header, use_rgermline, use_repertoire, use_productive,
    use_nonproductive, my_regions, my_var_elements, my_vars, my_vartypes, use_sharedVDJ,
    V_J_to_use, groups, group_A, group_B, group_C, univlog = F, samples_to_keep,
    variables_to_remove, pval_type, pval_cutoff, estimate_cutoff, number_selected_vars,
    VJ_deselected, VDJ_normalized_per_size, R_mut_threshold_min, R_mut_threshold_max,
    to_compare_groups, VDJ_maximize_clones, VDJ_normalized_per_sample, my_clone_def,
    seed, chains, igsubtypes
) {
    "%!in%" <- function(x, y) !(x %in%  y)
    columns <- c(1:ncol(FBM))
    rows <- c(1:nrow(FBM))
    index_repertoire <- which(merged_df$Sequence_type == "Repertoire")


    # print(index_repertoire)

    if (!is.null(samples_to_keep) &&
        samples_to_keep[1] != "RANDOMLETTERS") {
        #print("Filtering samples")

        rows <- rows[which(rows %in% which(merged_df$Patient_Sample %in% samples_to_keep))]

    }

    if (!use_rgermline) {
        rows <- rows[which(rows %in% index_repertoire)]
    }
    if (!use_repertoire) {
        rows <- rows[which(rows %!in% index_repertoire)]
    }

    if(!is.null(chains) && chains[1] != "Unknown") {
      rows <- rows[which(rows %in% which(merged_df$Chain %in% chains))]
    }

    if(!is.null(igsubtypes) && igsubtypes[1] != "Unknown") {
      rows <- rows[which(rows %in% which(merged_df$C_region %in% c(igsubtypes,"")))]
    }

    rows <- rows[which(
        rows %in% intersect(
            which(
                FBM[, which(merged_header == "AA_Whole_Replacement_muts_counts")] >=
                  R_mut_threshold_min
            ),
            which(
                FBM[, which(merged_header == "AA_Whole_Replacement_muts_counts")] <=
                  R_mut_threshold_max
            )
        )
    )]

    if (length(VJ_deselected) !=
        0) {

        rows <- rows[which(rows %in% which(merged_df$V_and_J %!in% VJ_deselected))]
    }

    set.seed(seed)

    if (use_sharedVDJ && to_compare_groups && !is.null(group_B) &&
        !is.null(group_A)) {
        # rows=rows[which(rows%!in%index_repertoire)] if(is.null(groups)){
        # groups='Groups' } list_VJs=list() for(group in
        # unique(merged_df[,get(groups)])) {
        # list_VJs[[length(list_VJs)+1]]=unique(merged_df[which(merged_df[,get(groups)]
        # == group),get('V_and_J')]) } VJs_allowed=Reduce(intersect, list_VJs)
        # print('VJ reduction') print(length(VJs_allowed))
        # rows=rows[which(rows%in%which(merged_df$V_and_J %in% VJs_allowed))]
        # print(length(rows))

        #print("Filtering by VJ")

        if (VDJ_normalized_per_size) {
            tmp_index <- c()

            if (VDJ_normalized_per_sample) {

                id_samples_A <- unique(
                  merged_df$Patient_Sample[intersect(
                    which(merged_df$Patient_Sample %in% samples_to_keep),
                    which(merged_df[[groups]] %in% group_A)
                )]
              )
                n_samples_A <- length(id_samples_A)
                id_samples_B <- unique(
                  merged_df$Patient_Sample[intersect(
                    which(merged_df$Patient_Sample %in% samples_to_keep),
                    which(merged_df[[groups]] %in% group_B)
                )]
              )
                n_samples_B <- length(id_samples_B)
                id_samples_C <- unique(
                  merged_df$Patient_Sample[intersect(
                    which(merged_df$Patient_Sample %in% samples_to_keep),
                    which(merged_df[[groups]] %in% group_C)
                )]
              )
                n_samples_C <- length(id_samples_C)
            } else {
                n_samples_A <- 1
                n_samples_B <- 1
                n_samples_C <- 1
            }
            for (V_J_comb in V_J_to_use) {
                tmp_group_A <- rows[which(
                  rows %in% intersect(
                    which(merged_df$V_and_J %in% V_J_comb),
                    which(merged_df[[groups]] %in% group_A)
                )
              )]
                tmp_group_B <- rows[which(
                  rows %in% intersect(
                    which(merged_df$V_and_J %in% V_J_comb),
                    which(merged_df[[groups]] %in% group_B)
                )
              )]
                tmp_group_C <- rows[which(
                  rows %in% intersect(
                    which(merged_df$V_and_J %in% V_J_comb),
                    which(merged_df[[groups]] %in% group_C)
                )
              )]

                if (VDJ_normalized_per_sample) {
                  min_size_A <- min(
                    sapply(
                      unique(
                        merged_df$Patient_Sample[intersect(
                          which(merged_df$Patient_Sample %in% samples_to_keep),
                          which(merged_df[[groups]] %in% group_A)
                      )]
                    ),
                      function(z) length(
                        intersect(
                          which(merged_df$Patient_Sample == z),
                          tmp_group_A
                      )
                    )
                  )
                )
                  min_size_B <- min(
                    sapply(
                      unique(
                        merged_df$Patient_Sample[intersect(
                          which(merged_df$Patient_Sample %in% samples_to_keep),
                          which(merged_df[[groups]] %in% group_B)
                      )]
                    ),
                      function(z) length(
                        intersect(
                          which(merged_df$Patient_Sample == z),
                          tmp_group_B
                      )
                    )
                  )
                )


                  min_size_AB <- min(c(min_size_A * n_samples_A, min_size_B * n_samples_B))

                  if (!is.null(n_samples_C) &&
                    n_samples_C > 0) {
                    min_size_C <- min(
                      sapply(
                        unique(
                          merged_df$Patient_Sample[intersect(
                            which(merged_df$Patient_Sample %in% samples_to_keep),
                            which(merged_df[[groups]] %in% group_C)
                        )]
                      ),
                        function(z) length(
                          intersect(
                            which(merged_df$Patient_Sample == z),
                            tmp_group_C
                        )
                      )
                    )
                  )
                    min_size_ABC <- min(
                      c(
                        min_size_A * n_samples_A, min_size_B * n_samples_B, min_size_C *
                          n_samples_C
                    )
                  )
                  }

                } else {
                  min_size_A <- length(tmp_group_A)
                  min_size_B <- length(tmp_group_B)
                  min_size_C <- length(tmp_group_C)

                  min_size_AB <- min(c(min_size_A, min_size_B))

                  if (!is.null(n_samples_C)) {
                    min_size_ABC <- min(c(min_size_A, min_size_B, min_size_C))
                  }

                }





                if (VDJ_normalized_per_sample) {
                  min_size_A <- floor(min_size_AB/n_samples_A)
                  min_size_B <- floor(min_size_AB/n_samples_B)
                  if (!is.null(group_C) && length(group_C) > 0 && n_samples_C > 0) {
                    min_size_C <- floor(min_size_ABC / n_samples_C)
                  } else {
                    min_size_C <- 0
                  }
                }

                if (min_size_AB == 0) {

                } else {

                  max_clones_A <- NULL

                  if (VDJ_maximize_clones) {


                    if (VDJ_normalized_per_sample) {
                      max_clones_A <- c()
                      tmp_tmp_group_A <- c()
                      for (n in c(1:n_samples_A)) {

                        sub_tmp_group_A <- tmp_group_A[which(
                          tmp_group_A %in% intersect(
                            intersect(
                              which(merged_df$V_and_J %in% V_J_comb),
                              which(merged_df[[groups]] %in% group_A)
                          ),
                            which(merged_df$Patient_Sample %in% id_samples_A[n])
                        )
                      )]

                        sub_all_clones_A <- sub_tmp_group_A[match(
                          unique(merged_df[sub_tmp_group_A, get(my_clone_def)]),
                          merged_df[sub_tmp_group_A, get(my_clone_def)]
                      )]
                        sub_all_clones_A <- sub_all_clones_A[1:min(min_size_A, length(sub_all_clones_A))]
                        if(length(sub_tmp_group_A[which(
                          merged_df[sub_tmp_group_A, get(my_clone_def)] %!in%
                          merged_df[max_clones_A, get(my_clone_def)]
                        )] > 1)) {
                          tmp_tmp_group_A <- c(
                            tmp_tmp_group_A, c(
                              sub_all_clones_A, sample(
                                sub_tmp_group_A[which(
                                  merged_df[sub_tmp_group_A, get(my_clone_def)] %!in%
                                    merged_df[max_clones_A, get(my_clone_def)]
                                )],
                                min_size_A - length(sub_all_clones_A)
                              )
                            )
                          )
                        } else {
                          tmp_tmp_group_A <- c(
                            tmp_tmp_group_A, c(
                              sub_all_clones_A))
                        }

                      }
                      if ((min_size_AB - min_size_A * n_samples_A) > 0) {
                        if(length(tmp_group_A[which(tmp_group_A %!in% tmp_tmp_group_A)])>1){
                          tmp_tmp_group_A <- c(
                            tmp_tmp_group_A, sample(
                              tmp_group_A[which(tmp_group_A %!in% tmp_tmp_group_A)],
                              min_size_AB - min_size_A * n_samples_A
                            )
                          )
                        } else {
                          tmp_tmp_group_A <- c(tmp_tmp_group_A, tmp_group_A[which(tmp_group_A %!in% tmp_tmp_group_A)])
                        }

                      }

                      tmp_group_A <- tmp_tmp_group_A
                    } else {
                      all_clones_A <- tmp_group_A[match(
                        unique(merged_df[tmp_group_A, get(my_clone_def)]),
                        merged_df[tmp_group_A, get(my_clone_def)]
                        )]
                      if (length(all_clones_A) > 1){
                        max_clones_A <- sample(all_clones_A, min(min_size_AB, length(all_clones_A)))
                      } else {
                        max_clones_A <- all_clones_A
                      }

                      if(length(tmp_group_A[which(tmp_group_A %!in% max_clones_A)]) > 1) {
                        tmp_group_A <- c(
                          max_clones_A, sample(
                            tmp_group_A[which(tmp_group_A %!in% max_clones_A)],
                            min_size_AB - length(max_clones_A)
                          ))
                      } else {
                        if((min_size_AB - length(max_clones_A))>0) {
                          tmp_group_A <- c(
                            max_clones_A, tmp_group_A[which(tmp_group_A %!in% max_clones_A)])
                        } else {
                          tmp_group_A=max_clones_A
                        }

                      }



                    }

                  } else {

                    if (VDJ_normalized_per_sample) {

                      tmp_tmp_group_A <- c()
                      for (n in c(1:n_samples_A)) {

                        sub_tmp_group_A <- tmp_group_A[which(
                          tmp_group_A %in% intersect(
                            intersect(
                              which(merged_df$V_and_J %in% V_J_comb),
                              which(merged_df[[groups]] %in% group_A)
                          ),
                            which(merged_df$Patient_Sample %in% id_samples_A[n])
                        )
                      )]
                        tmp_tmp_group_A <- c(tmp_tmp_group_A,  sub_tmp_group_A[sample.int(length(sub_tmp_group_A), min_size_A)])
                      }
                      if ((min_size_AB - min_size_A * n_samples_A) > 0) {

                        if(length(tmp_group_A[which(tmp_group_A %!in% tmp_tmp_group_A)])>1){
                          tmp_tmp_group_A <- c(
                            tmp_tmp_group_A, sample(
                              tmp_group_A[which(tmp_group_A %!in% tmp_tmp_group_A)],
                              min_size_AB - min_size_A * n_samples_A
                            ))
                        } else {
                          tmp_tmp_group_A <- c(
                            tmp_tmp_group_A, tmp_group_A[which(tmp_group_A %!in% tmp_tmp_group_A)]
                            )
                        }


                      }


                      tmp_group_A <- tmp_tmp_group_A
                    } else {
                      if (length(tmp_group_A)>1){
                        tmp_group_A <- sample(tmp_group_A, min_size_AB)
                      }
                    }


                  }

                  max_clones_B <- NULL

                  if (VDJ_maximize_clones) {


                    if (VDJ_normalized_per_sample) {
                      max_clones_B <- c()
                      tmp_tmp_group_B <- c()
                      for (n in c(1:n_samples_B)) {

                        sub_tmp_group_B <- tmp_group_B[which(
                          tmp_group_B %in% intersect(
                            intersect(
                              which(merged_df$V_and_J %in% V_J_comb),
                              which(merged_df[[groups]] %in% group_B)
                          ),
                            which(merged_df$Patient_Sample %in% id_samples_B[n])
                        )
                      )]

                        sub_all_clones_B <- sub_tmp_group_B[match(
                          unique(merged_df[sub_tmp_group_B, get(my_clone_def)]),
                          merged_df[sub_tmp_group_B, get(my_clone_def)]
                      )]
                        sub_all_clones_B <- sub_all_clones_B[1:min(min_size_B, length(sub_all_clones_B))]

                        if(length(sub_tmp_group_B[which(
                          merged_df[sub_tmp_group_B, get(my_clone_def)] %!in%
                          merged_df[max_clones_B, get(my_clone_def)]
                        )])>1){
                          tmp_tmp_group_B <- c(
                            tmp_tmp_group_B, c(
                              sub_all_clones_B, sample(
                                sub_tmp_group_B[which(
                                  merged_df[sub_tmp_group_B, get(my_clone_def)] %!in%
                                    merged_df[max_clones_B, get(my_clone_def)]
                                )],
                                min_size_B - length(sub_all_clones_B)
                              )
                            )
                          )
                        } else {
                          if ((min_size_B - length(sub_all_clones_B))>0){
                            tmp_tmp_group_B <- c(
                              tmp_tmp_group_B, sub_all_clones_B,
                              sub_tmp_group_B[which(
                                merged_df[sub_tmp_group_B, get(my_clone_def)] %!in%
                                  merged_df[max_clones_B, get(my_clone_def)]
                              )]
                            )
                          }
                        }

                      }
                      if ((min_size_AB - min_size_B * n_samples_B) > 0) {
                        if(length(tmp_group_B[which(tmp_group_B %!in% tmp_tmp_group_B)])>1) {
                          tmp_tmp_group_B <- c(
                            tmp_tmp_group_B, sample(
                              tmp_group_B[which(tmp_group_B %!in% tmp_tmp_group_B)],
                              min_size_AB - min_size_B * n_samples_B
                            ))
                        } else {
                          tmp_tmp_group_B <- c(
                            tmp_tmp_group_B, tmp_group_B[which(tmp_group_B %!in% tmp_tmp_group_B)]
                            )
                        }
                      }

                      tmp_group_B <- tmp_tmp_group_B
                    } else {
                      all_clones_B <- tmp_group_B[match(
                        unique(merged_df[tmp_group_B, get(my_clone_def)]),
                        merged_df[tmp_group_B, get(my_clone_def)]
                      )]
                      if (length(all_clones_B) > 1){
                        max_clones_B <- sample(all_clones_B, min(min_size_AB, length(all_clones_B)))
                      } else {
                        max_clones_B <- all_clones_B
                      }

                      if(length(tmp_group_B[which(tmp_group_B %!in% max_clones_B)]) > 1) {
                        tmp_group_B <- c(
                          max_clones_B, sample(
                            tmp_group_B[which(tmp_group_B %!in% max_clones_B)],
                            min_size_AB - length(max_clones_B)
                          ))
                      } else {
                        if((min_size_AB - length(max_clones_B))>0) {
                          tmp_group_B <- c(
                            max_clones_B, tmp_group_B[which(tmp_group_B %!in% max_clones_B)])
                        } else {
                          tmp_group_B=max_clones_B
                        }

                      }
                    }

                  } else {

                    if (VDJ_normalized_per_sample) {

                      tmp_tmp_group_B <- c()
                      for (n in c(1:n_samples_B)) {

                        sub_tmp_group_B <- tmp_group_B[which(
                          tmp_group_B %in% intersect(
                            intersect(
                              which(merged_df$V_and_J %in% V_J_comb),
                              which(merged_df[[groups]] %in% group_B)
                          ),
                            which(merged_df$Patient_Sample %in% id_samples_B[n])
                        )
                      )]

                        tmp_tmp_group_B <- c(tmp_tmp_group_B, sub_tmp_group_B[sample.int(length(sub_tmp_group_B), min_size_B)])
                      }
                      if ((min_size_AB - min_size_B * n_samples_B) > 0) {
                        if(length(tmp_group_B[which(tmp_group_B %!in% tmp_tmp_group_B)])>1){
                          tmp_tmp_group_B <- c(
                            tmp_tmp_group_B, sample(
                              tmp_group_B[which(tmp_group_B %!in% tmp_tmp_group_B)],
                              min_size_AB - min_size_B * n_samples_B
                            ))
                        } else {
                          tmp_tmp_group_B <- c(tmp_tmp_group_B,
                                               tmp_group_B[which(tmp_group_B %!in% tmp_tmp_group_B)])
                        }

                      }


                      tmp_group_B <- tmp_tmp_group_B


                    } else {
                      if (length(tmp_group_B)>1){
                        tmp_group_B <- sample(tmp_group_B, min_size_AB)
                      }

                    }


                  }


                  # if(length(tmp_group_C)>0){ # tmp_index=c(tmp_index,
                  # tmp_group_C[1:min(c(length(tmp_group_A),length(tmp_group_B),
                  # length(tmp_group_C)))]) max_clones_C=NULL
                  # if(VDJ_maximize_clones) {
                  # max_clones_C=tmp_group_C[match(unique(merged_df[tmp_group_C,which(colnames(merged_df)==my_clone_def)]),
                  # merged_df[tmp_group_C,which(colnames(merged_df)==my_clone_def)])]
                  # max_clones_C=sample( tmp_group_C, min_size_ABC) }
                  # tmp_group_C= c(max_clones_C,
                  # sample(tmp_group_C[which(tmp_group_C %!in% max_clones_C)],
                  # min_size_ABC)) }

                  if (length(tmp_group_C) >
                    0) {
                    max_clones_C <- NULL

                    if (VDJ_maximize_clones) {


                      if (VDJ_normalized_per_sample) {
                        max_clones_C <- c()
                        tmp_tmp_group_C <- c()
                        for (n in c(1:n_samples_C)) {

                          sub_tmp_group_C <- tmp_group_C[which(
                            tmp_group_C %in% intersect(
                              intersect(
                                which(merged_df$V_and_J %in% V_J_comb),
                                which(merged_df[[groups]] %in% group_C)
                            ),
                              which(merged_df$Patient_Sample %in% id_samples_C[n])
                          )
                        )]

                          sub_all_clones_C <- sub_tmp_group_C[match(
                            unique(merged_df[sub_tmp_group_C, get(my_clone_def)]),
                            merged_df[sub_tmp_group_C, get(my_clone_def)]
                        )]
                          sub_all_clones_C <- sub_all_clones_C[1:min(min_size_C, length(sub_all_clones_C))]

                          if(length(sub_tmp_group_C[which(
                            merged_df[sub_tmp_group_C, get(my_clone_def)] %!in%
                            merged_df[max_clones_C, get(my_clone_def)]
                          )])>1) {
                            tmp_tmp_group_C <- c(
                              tmp_tmp_group_C, c(
                                sub_all_clones_C, sample(
                                  sub_tmp_group_C[which(
                                    merged_df[sub_tmp_group_C, get(my_clone_def)] %!in%
                                      merged_df[max_clones_C, get(my_clone_def)]
                                  )],
                                  min_size_C - length(sub_all_clones_C)
                                )
                              )
                            )
                          } else {
                            if((min_size_C - length(sub_all_clones_C))>0) {
                              tmp_tmp_group_C <- c(
                                tmp_tmp_group_C, c(
                                  sub_all_clones_C, sub_tmp_group_C[which(
                                    merged_df[sub_tmp_group_C, get(my_clone_def)] %!in%
                                      merged_df[max_clones_C, get(my_clone_def)]
                                  )]
                                )
                              )
                            }
                          }

                        }
                        if ((min_size_ABC - min_size_C * n_samples_C) > 0) {
                          tmp_tmp_group_C <- c(
                            tmp_tmp_group_C, sample(
                              tmp_group_C[which(tmp_group_C %!in% tmp_tmp_group_C)],
                              min_size_ABC - min_size_C * n_samples_C
                            )
                          )
                        } else {
                          tmp_tmp_group_C <- c(
                            tmp_tmp_group_C, c(
                              sub_all_clones_C))
                        }

                        tmp_group_C <- tmp_tmp_group_C
                      } else {
                        all_clones_C <- tmp_group_C[match(
                          unique(merged_df[tmp_group_C, get(my_clone_def)]),
                          merged_df[tmp_group_C, get(my_clone_def)]
                        )]
                        if (length(all_clones_C) > 1){
                          max_clones_C <- sample(all_clones_C, min(min_size_ABC, length(all_clones_C)))
                        } else {
                          max_clones_C <- all_clones_C
                        }

                        if(length(tmp_group_C[which(tmp_group_C %!in% max_clones_C)]) > 1) {
                          tmp_group_C <- c(
                            max_clones_C, sample(
                              tmp_group_C[which(tmp_group_C %!in% max_clones_C)],
                              min_size_ABC - length(max_clones_C)
                            ))
                        } else {
                          if((min_size_ABC - length(max_clones_C))>0) {
                            tmp_group_C <- c(
                              max_clones_C, tmp_group_C[which(tmp_group_C %!in% max_clones_C)])
                          } else {
                            tmp_group_C=max_clones_C
                          }

                        }
                      }

                    } else {

                      if (VDJ_normalized_per_sample) {

                        tmp_tmp_group_C <- c()
                        for (n in c(1:n_samples_C)) {

                          sub_tmp_group_C <- tmp_group_C[which(
                            tmp_group_C %in% intersect(
                              intersect(
                                which(merged_df$V_and_J %in% V_J_comb),
                                which(merged_df[[groups]] %in% group_C)
                            ),
                              which(merged_df$Patient_Sample %in% id_samples_C[n])
                          )
                        )]

                          tmp_tmp_group_C <- c(tmp_tmp_group_C,
                                               sub_tmp_group_C[sample.int(length(sub_tmp_group_C), min_size_C)])
                        }
                        if ((min_size_AB - min_size_C * n_samples_C) > 0) {
                          tmp_tmp_group_C <- c(
                            tmp_tmp_group_C,
                            tmp_group_C[which(tmp_group_C %!in% tmp_tmp_group_C)][sample.int(length(tmp_group_C[which(tmp_group_C %!in% tmp_tmp_group_C)]),
                                                                                             min_size_ABC - min_size_C * n_samples_C)]

                            )
                        }


                        tmp_group_C <- tmp_tmp_group_C


                      } else {
                        if (length(tmp_group_C)>1){
                          tmp_group_C <- sample(c(tmp_group_C), min_size_ABC)
                        }

                      }


                    }

                  } else {
                    tmp_group_C <- NULL
                  }

                  # No random tmp_index=c(tmp_index,
                  # c(tmp_group_A[1:min(c(length(tmp_group_A),length(tmp_group_B)))],
                  # tmp_group_B[1:min(c(length(tmp_group_A),length(tmp_group_B)))])
                  # )

                  tmp_index <- c(tmp_index, c(tmp_group_A, tmp_group_B, tmp_group_C))

                }



            }


            rows <- rows[which(rows %in% tmp_index)]
        } else {
            rows <- rows[which(rows %in% which(merged_df$V_and_J %in% V_J_to_use))]
        }


    }







    if (!is.null(my_regions)) {
        columns <- columns[which(
            columns %in% which(
                (sapply(merged_header, function(x) strsplit(x, split = "_")[[1]][2])) %in%
                  my_regions
            )
        )]
    } else {
        columns <- c()
    }

    if (!is.null(my_var_elements)) {
        columns <- columns[which(
            columns %in% which(
                grepl(
                  paste0(
                    paste("^", my_var_elements, sep = ""),
                    collapse = "|"
                ),
                  merged_header
              )
            )
        )]
    } else {
        columns <- c()
    }

    if (!is.null(my_vars)) {

        merged_header_splitted_at_3 <- sapply(merged_header, function(x) strsplit(x, split = "_")[[1]][3])
        merged_header_splitted_at_4 <- sapply(merged_header, function(x) strsplit(x, split = "_")[[1]][4])
        merged_header_splitted_at_3_and_4 <- sapply(
            merged_header, function(x) paste(
                strsplit(x, split = "_")[[1]][3],
                strsplit(x, split = "_")[[1]][4],
                sep = "_"
            )
        )

        if ("Length" %!in% my_vars) {
            columns <- columns[which(columns %!in% which(endsWith(merged_header, "length")))]
            columns <- columns[which(columns %!in% which(endsWith(merged_header, "length_Diff")))]
        }

        if ("Composition" %!in% my_vars) {
            nucleotides <- c("A", "G", "T", "C")
            peptides <- Peptides::aaList()
            columns <- columns[which(
                columns %!in% which(
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
            )]


            columns <- columns[which(columns %!in% which(merged_header_splitted_at_4 %in% c("counts")))]
        }

        if ("Hot/Cold motifs" %!in% my_vars) {
            columns <- columns[which(
                columns %!in% which(merged_header_splitted_at_3 %in% c("hot", "cold", "potential"))
            )]
        }

        if ("Substitutions" %!in% my_vars) { ###Outdated
            columns <- columns[which(
                columns %!in% which(
                  merged_header_splitted_at_3_and_4 %in% c("Sub", "Sub_prc", "SIDT_sum", "SID_sum")
              )
            )]
            columns <- columns[which(columns %!in% which(merged_header_splitted_at_3 %in% c("Sub")))]
            columns <- columns[which(columns %!in% which(endsWith(merged_header, "_S")))]  ##Old, remove later
            columns <- columns[which(columns %!in% which(endsWith(merged_header, "_S_prc")))]  ##Old, remove later
        }

        if ("Insertions" %!in% my_vars) {
            columns <- columns[which(
                columns %!in% which(
                  merged_header_splitted_at_3 %in% c("NumIns", "LenghtIns")
              )
            )]
            # columns <- columns[which(columns %!in% which(merged_header_splitted_at_3 %in% c("Ins")))]
            # columns <- columns[which(columns %!in% which(endsWith(merged_header, "_I")))]  ##Old, remove later
            # columns <- columns[which(columns %!in% which(endsWith(merged_header, "_I_prc")))]  ##Old, remove later
        }
        if ("Deletions" %!in% my_vars) {
            columns <- columns[which(
                columns %!in% which(
                  merged_header_splitted_at_3 %in% c("LenghtDels", "NumDels")
              )
            )]
            # columns <- columns[which(columns %!in% which(merged_header_splitted_at_3 %in% c("Del")))]
            # columns <- columns[which(columns %!in% which(endsWith(merged_header, "_D")))]  ##Old, remove later
            # columns <- columns[which(columns %!in% which(endsWith(merged_header, "_D_prc")))]  ##Old, remove later
        }
        if ("Translocations" %!in% my_vars) { ###Outdated
            columns <- columns[which(
                columns %!in% which(merged_header_splitted_at_3_and_4 %in% c("Trasl", "Trasl_prc", "SIDT_sum"))
            )]
            columns <- columns[which(columns %!in% which(merged_header_splitted_at_3 %in% c("Trasl")))]
            columns <- columns[which(columns %!in% which(endsWith(merged_header, "_T")))]  ##Old, remove later
            columns <- columns[which(columns %!in% which(endsWith(merged_header, "_T_prc")))]  ##Old, remove later
        }

        if ("Leveshtein distance" %!in% my_vars) {  ###Outdated
            columns <- columns[which(columns %!in% which(merged_header_splitted_at_3 %in% c("lv")))]
        }

        if ("Transitions and transversions" %!in% my_vars) {
            columns <- columns[which(
                columns %!in% which(merged_header_splitted_at_3 %in% c("Transitions", "Transversions"))
            )]
            columns <- columns[which(
                columns %!in% which(merged_header_splitted_at_3_and_4 %in% c("Ratio_Transitions-Transversions"))
            )]
        }
        if ("Replacement and silent mutations" %!in% my_vars) {
            columns <- columns[which(
                columns %!in% which(merged_header_splitted_at_3 %in% c("Replacement", "Silent", "Total"))
            )]
            columns <- columns[which(
                columns %!in% which(merged_header_splitted_at_3_and_4 %in% c("Ratio_Replacement-Silent"))
            )]
        }

        if ("Mutations from X to Y" %!in% my_vars) {
            columns <- columns[which(columns %!in% which(merged_header_splitted_at_4 %in% c("to")))]
            # columns=columns[which(columns %!in% which(sapply(merged_header,
            # function(x) paste(strsplit(x, split='_')[[1]][3], strsplit(x,
            # split='_')[[1]][4], sep='_')) %in% c('Ratio_Silent') ))]
        }

        if ("NGly sites" %!in% my_vars) {
            columns <- columns[which(columns %!in% which(merged_header_splitted_at_3 %in% c("NGly")))]
        }
        if ("Peptide features" %!in% my_vars) {

          tmp_columns <- columns[intersect(which(
            columns %in% which(
              merged_header_splitted_at_3 %in% c(
                "Peptides", "alkzm", "Small", "Tiny", "Aliphatic", "Charged",
                "Polar", "Basic", "NonPolar", "Aromatic", "Acidic"
              )
            )
          ),which(columns %in% which(merged_header_splitted_at_4 %in% c("to"))))]

            columns <- columns[which(
                columns %!in% which(
                  merged_header_splitted_at_3 %in% c(
                    "Peptides", "alkzm", "Small", "Tiny", "Aliphatic", "Charged",
                    "Polar", "Basic", "NonPolar", "Aromatic", "Acidic"
                )
              )
            )]

            columns=c(columns, tmp_columns)
            # columns=columns[which(columns %!in% which(sapply(merged_header,
            # function(x) strsplit(x, split='_')[[1]][3]) %in% c('alkzm') ))]
            # columns=columns[which(columns %!in% which(sapply(merged_header,
            # function(x) strsplit(x, split='_')[[1]][3]) %in% c('alkzm') ))]
        }


    } else {
        columns <- c()
    }

    if (!is.null(my_vartypes)) {
        if ("Germline diff" %!in% my_vartypes) {
            columns <- columns[which(columns %!in% which(endsWith(merged_header, "_Diff")))]
        }
        if ("Baseline" %!in% my_vartypes) {
            columns <- columns[which(columns %in% which(endsWith(merged_header, "_Diff")))]
        }
    } else {
        columns <- c()
    }

    if (!is.null(variables_to_remove)) {
        columns <- columns[which(columns %!in% which(merged_header %in% variables_to_remove))]
    }



    big_summary <- bigstatsr::big_colstats(FBM, ind.row = rows, ind.col = columns, ncores = bigstatsr::nb_cores())
    vars_to_remove <- columns[unique(
        c(
            which(big_summary$sum == 0),
            which(big_summary$var == 0),
            which(is.na(big_summary$sum)),
            which(is.infinite(big_summary$sum)),
            which(is.nan(big_summary$sum)),
            which(is.na(big_summary$var)),
            which(is.infinite(big_summary$var)),
            which(is.nan(sqrt(big_summary$var)))
        )
    )]

    # if(any(is.nan(sqrt(big_summary$var)))) { print(FBM[rows,
    # columns[which(is.nan(sqrt(big_summary$var)))]]) }
    columns <- columns[which(columns %!in% vars_to_remove)]

    if (univlog && !is.null(columns) &&
        !is.null(group_B) &&
        !is.null(group_A) &&
        to_compare_groups) {

        merged_df$Binary_classif <- 4
        merged_df$Binary_classif[which(merged_df[[groups]] %in% group_A)] <- 0
        merged_df$Binary_classif[which(merged_df[[groups]] %in% group_B)] <- 1
        tmp_rows <- rows[(which(rows %in% which(merged_df[[groups]] %in% c(group_B, group_A))))]

        # print(length(tmp_rows))
        # print(table(merged_df$Binary_classif[tmp_rows]))
        # if(any(merged_df$Binary_classif[tmp_rows] %!in% c(0,1))){
        #   print("Unusual")
        #   print(which(merged_df$Binary_classif[tmp_rows] %!in% c(0,1)))
        #   print(merged_df$Binary_classif[tmp_rows][which(merged_df$Binary_classif[tmp_rows] %!in% c(0,1))])
        # }
        testuniv <- big_univLogReg(
            FBM, merged_df$Binary_classif[tmp_rows], ind.train = tmp_rows, ind.col = columns,
            ncores = nb_cores()
        )

        testuniv$p.value <- stats::predict(testuniv, log10 = FALSE)


        if (pval_type == "Corrected by Benjamini & Hochberg method") {
            testuniv$p.value <- round(
                stats::p.adjust(testuniv$p.value, "BH"),
                6
            )
        } else if (pval_type == "Corrected by Bonferroni method") {
            testuniv$p.value <- round(
                stats::p.adjust(testuniv$p.value, "bonferroni"),
                6
            )
        } else if (pval_type == "Corrected by Benjamini & Yekutieli method") {
          testuniv$p.value <- round(
            stats::p.adjust(testuniv$p.value, "BY"),
            6
          )
        }


        intersection <- intersect(
          which(
            abs(testuniv$p.value) <=
              pval_cutoff
          ),
          which(
            abs(testuniv$estim) >=
              estimate_cutoff
          )
        )
        length(intersection)

        if (number_selected_vars == "All") {
            columns <- columns[intersection]
        } else if (number_selected_vars == "Top 10 (according to estimate)") {

            columns <- columns[intersection[sort(
                abs(testuniv$estim)[intersection],
                index.return = TRUE, decreasing = TRUE
            )$ix[1:min(10, length(intersection))]]]

        } else if (number_selected_vars == "Top 50 (according to estimate)") {

            columns <- columns[intersection[sort(
                abs(testuniv$estim)[intersection],
                index.return = TRUE, decreasing = TRUE
            )$ix[1:min(50, length(intersection))]]]
        }



    }

    columns=columns[which(!is.na(columns))]




    return(list(ROWS = rows, COLUMNS = columns))
}

#' Big PCA
#'
#' @description A function to perform PCA over a FBM object.
#' @param FBM A FBM object.
#' @param rows Index of rows of the FBM that will be used to calculate the PCA.
#' @param columns Index of columns of the FBM that will be used to calculate the PCA..
#' @return A list with the PCA scores and its explained variances
#' @import bigstatsr
#' @keywords internal
#' @export
big_PCA <- function(FBM, rows, columns) {

    if (length(columns) >
        0 && length(rows) >
        0) {

        PCA_error_perhaps <- tryCatch(
            {

              testing <- bigstatsr::big_randomSVD(
                FBM, fun.scaling = big_scale(), ind.row = rows, ind.col = columns, k = 5,
                verbose = F
              )

              scores <- stats::predict(testing)

                var_exp <- testing$d^2/big_norm(
                  FBM, ind.row = rows, ind.col = columns, center = testing$center,
                  scale = testing$scale
              )
                variance <- signif(var_exp, 2)

                cumulative_variance <- round(
                  cumsum(var_exp),
                  3
              )
                # Contribuciones (%): para estandarizado es simplemente 100 * v^2
                contrib <- 100 * testing$v^2
                # Comprobación: cada columna debe sumar 100
                colSums(contrib)  # ~ 100

                # Nombres útiles (si tienes nombres de variables/PCs)
                # rownames(contrib) <- colnames(mat)  # si partiste de 'mat'
                colnames(contrib) <- paste0("PCA", seq_len(ncol(contrib)))

                list(Scores = scores, Variance_explained = variance,
                     Variance=contrib)
            }, error = function(e) {
                # print("Caught an error!")
                list(Scores = NULL, Variance_explained = NULL,
                     Variance=NULL)
            }, warning = function(w) {
                # print("Caught an warning!")
                list(Scores = NULL, Variance_explained = NULL,
                     Variance=NULL)
            }
        )
        return(
            list(
                Scores = PCA_error_perhaps$Scores,
                Variance_explained = PCA_error_perhaps$Variance_explained,
                Variance=PCA_error_perhaps$Variance
            )
        )

    } else {
        # print("No PCA")
        return(list(Scores = NULL, Variance_explained = NULL,
                    Variance = NULL))
    }


}
