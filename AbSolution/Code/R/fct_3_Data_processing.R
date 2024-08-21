#' Merge FBMs and meta files
#'
#' @description A function to join by row individual FBMs and .info files
#'     into a unique FBM and a unique .info file, merged_bm.
#' @param filepath Path to the individual FBMs.
#' @return A list with the merged FBM, the merged info file and the FBM header.
#' @importFrom utils  glob2rx
#' @import bigstatsr
#' @import data.table
#' @importFrom tools file_path_sans_ext
#' @import bigparallelr
#' @noRd
merge_FBMs =function(filepath) {
  basenames=basename(list.files(path=filepath, pattern=glob2rx("*.rds")))
  prefixes=tools::file_path_sans_ext(basenames)

  list_descriptions=list()
  list_FBMs=list()
  list_headers=list()
  for (prefix in unique(prefixes)[which(unique(prefixes) != "Merged_bm")] ){
    list_descriptions[[length(list_descriptions)+1]] =data.table::fread(file.path(filepath, paste0(prefix, ".info")))
    list_headers[[length(list_headers)+1]] =data.table::fread(file.path(filepath, paste0(prefix, ".example")))
    if(file.exists(file.path(filepath, paste0(prefix, ".rds")))) {
      list_FBMs[[length(list_FBMs)+1]] =bigstatsr::big_attach(file.path(filepath,paste0(prefix, ".rds")))
    }

  }


  merged_df=rbindlist( list_descriptions )
  merged_header=colnames(list_headers[[1]])

  for (i in c(1:length(list_headers)) ) {
    header_tmp=colnames(list_headers[[i]])
    if(all(length(merged_header)==length(header_tmp)) && all(merged_header==header_tmp)) {

    } else {
      print("[ERROR] different number of columns and/or colnames")
    }

  }



  m <- ncol(list_FBMs[[1]]) # assuming that all have the same number of columns
  n <- sum(sapply(list_FBMs, nrow))

  if (file.exists(file.path(filepath, paste0("Merged_bm", ".bin")))) {
    file.remove(file.path(filepath, paste0("Merged_bm", ".bin")))
    if (file.exists(file.path(filepath, paste0("Merged_bm", ".desc")))) {
      file.remove(file.path(filepath, paste0("Merged_bm", ".desc")))
    }

  }

  merged_del=T
  if (file.exists(file.path(filepath, paste0("Merged_bm", ".bk")))) {

    if(file.exists(file.path(filepath, paste0("Merged_bm", ".info")))) {
      tmp_merged_df=data.table::fread(file.path(filepath, paste0("Merged_bm", ".info")))
      if (nrow(merged_df)==nrow(tmp_merged_df) & all(tmp_merged_df[,1]==merged_df[,1])& all(tmp_merged_df[,2]==merged_df[,2]) & all(tmp_merged_df[,3]==merged_df[,3])) {
        merged_del=F
        merged_bm=bigstatsr::big_attach(file.path(filepath,paste0("Merged_bm", ".rds")))
        if(nrow(merged_bm) != nrow(tmp_merged_df)) {
          merged_del=T
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

    merged_bm=bigstatsr::FBM(n,m,backingfile = file.path(filepath,"Merged_bm"), is_read_only = F)
    fwrite(merged_df, file.path(filepath, paste0("Merged_bm", ".info")), sep="\t")


    ### Contribution by F. Privé in https://github.com/privefl/bigstatsr/issues/176
    bigstatsr::big_apply(merged_bm, function(X, ind, list_fbm) {
      X[, ind] <- do.call("rbind", lapply(list_fbm, function(fbm) fbm[, ind, drop = FALSE]))
      NULL
    }, list_fbm = list_FBMs, a.combine  = "c", ncores = nb_cores())


    print("flush")
    merged_bm$save()
    print("flushed")
  }

  return(list(Header=merged_header, Short_DF=merged_df, Big_DF=merged_bm))

}

#' 3_Data_processing
#'
#' @description A fct function
#'
#' @return The return value, if any, from executing the function.
#' @import Peptides
#' @import bigstatsr
#' @import bigparallelr
#' @noRd
filter_merged=function(copy_FBM, merged_df, merged_header, use_rgermline,
                       use_repertoire, use_productive,use_nonproductive,
                       my_regions, my_var_elements, my_vars, my_vartypes,
                       use_sharedVDJ, V_J_to_use, groups, group_A, group_B,
                       group_C, univlog=F, samples_to_keep, variables_to_remove,
                       pval_type, pval_cutoff, estimate_cutoff,
                       number_selected_vars, VJ_deselected,
                       VDJ_normalized_per_size, R_mut_threshold_min,
                       R_mut_threshold_max, to_compare_groups,
                       VDJ_maximize_clones, VDJ_normalized_per_sample,
                       my_clone_def) {
  '%!in%' <- function(x,y)!('%in%'(x,y))
  columns=c(1:ncol(copy_FBM))
  rows=c(1:nrow(copy_FBM))
  index_repertoire=which(merged_df$Sequence_type=="Repertoire")

  print("Original sizes: ")
  print(length(rows))
  print(length(columns))
  # print(index_repertoire)

  print(samples_to_keep)
  if(!is.null(samples_to_keep) && samples_to_keep[1] != "RANDOMLETTERS") {
    print("Filtering samples")
    print(length(rows))
    rows=rows[which(rows%in%which(merged_df$Patient_Sample %in% samples_to_keep))]
    print(length(rows))
  }

  if(!use_rgermline){
    rows=rows[which(rows%in%index_repertoire)]
  }
  if(!use_repertoire){
    rows=rows[which(rows%!in%index_repertoire)]
  }

  rows=rows[which(rows %in% intersect(which(copy_FBM[,which(merged_header == "AA_Whole_Replacement_muts")]>=R_mut_threshold_min),
                                      which(copy_FBM[,which(merged_header == "AA_Whole_Replacement_muts")] <=R_mut_threshold_max)))]


  set.seed(1234)
  if(use_sharedVDJ && to_compare_groups && !is.null(group_B) && !is.null(group_A)){
    # rows=rows[which(rows%!in%index_repertoire)]
    # if(is.null(groups)){
    #   groups="Groups"
    # }
    #
    # list_VJs=list()
    # for(group in unique(merged_df[,get(groups)])) {
    #   list_VJs[[length(list_VJs)+1]]=unique(merged_df[which(merged_df[,get(groups)] == group),get("V_and_J")])
    # }
    # VJs_allowed=Reduce(intersect, list_VJs)
    # print("VJ reduction")
    # print(length(VJs_allowed))
    # rows=rows[which(rows%in%which(merged_df$V_and_J %in% VJs_allowed))]
    # print(length(rows))

    print("Filtering by VJ")
    print(length(rows))
    if(VDJ_normalized_per_size){
      tmp_index=c()

      if (VDJ_normalized_per_sample) {

        id_samples_A=unique(merged_df$Patient_Sample[intersect(which(merged_df$Patient_Sample %in% samples_to_keep), which(merged_df[[groups]] %in% group_A))])
        n_samples_A=length(id_samples_A)
        id_samples_B=unique(merged_df$Patient_Sample[intersect(which(merged_df$Patient_Sample %in% samples_to_keep), which(merged_df[[groups]] %in% group_B))])
        n_samples_B=length(id_samples_B)
        id_samples_C=unique(merged_df$Patient_Sample[intersect(which(merged_df$Patient_Sample %in% samples_to_keep), which(merged_df[[groups]] %in% group_C))])
        n_samples_C=length(id_samples_C)
      } else {
        n_samples_A=1
        n_samples_B=1
        n_samples_C=1
      }
      for (V_J_comb in V_J_to_use) {
        tmp_group_A=rows[which(rows %in% intersect(which(merged_df$V_and_J %in% V_J_comb),which(merged_df[[groups]] %in% group_A) ))]
        tmp_group_B=rows[which(rows %in% intersect(which(merged_df$V_and_J %in% V_J_comb), which(merged_df[[groups]] %in% group_B)))]
        tmp_group_C=rows[which(rows %in% intersect(which(merged_df$V_and_J %in% V_J_comb), which(merged_df[[groups]] %in% group_C)))]

        if (VDJ_normalized_per_sample) {
          min_size_A=min(sapply(unique(merged_df$Patient_Sample[intersect(which(merged_df$Patient_Sample %in% samples_to_keep), which(merged_df[[groups]] %in% group_A))]),
                                function(z) length(intersect(which(merged_df$Patient_Sample == z), tmp_group_A)) ))
          min_size_B=min(sapply(unique(merged_df$Patient_Sample[intersect(which(merged_df$Patient_Sample %in% samples_to_keep), which(merged_df[[groups]] %in% group_B))]),
                                function(z) length(intersect(which(merged_df$Patient_Sample == z), tmp_group_B)) ))


          min_size_AB=min(c(min_size_A * n_samples_A, min_size_B * n_samples_B))

          if(!is.null(n_samples_C) && n_samples_C >0){
            min_size_C=min(sapply(unique(merged_df$Patient_Sample[intersect(which(merged_df$Patient_Sample %in% samples_to_keep), which(merged_df[[groups]] %in% group_C))]),
                                  function(z) length(intersect(which(merged_df$Patient_Sample == z), tmp_group_B)) ))
            min_size_ABC=min(c(min_size_A * n_samples_A, min_size_B * n_samples_B,  min_size_C * n_samples_C ))
          }

        } else {
          min_size_A=length(tmp_group_A)
          min_size_B=length(tmp_group_B)
          min_size_C=length(tmp_group_C)

          min_size_AB=min(c(min_size_A, min_size_B ))

          if(!is.null(n_samples_C)) {
            min_size_ABC=min(c(min_size_A, min_size_B,  min_size_C ))
          }

        }





        if (VDJ_normalized_per_sample) {
          min_size_A=floor(min_size_AB/n_samples_A)
          min_size_B=floor(min_size_AB/n_samples_B)
          min_size_C=floor(min_size_AB/n_samples_C)
        }

        if (min_size_AB==0) {

        } else {

          max_clones_A=NULL

          if(VDJ_maximize_clones) {


            if(VDJ_normalized_per_sample) {
              max_clones_A=c()
              tmp_tmp_group_A=c()
              for (n in c(1:n_samples_A)) {

                sub_tmp_group_A=tmp_group_A[which(tmp_group_A %in% intersect(intersect(which(merged_df$V_and_J %in% V_J_comb),which(merged_df[[groups]] %in% group_A) ), which(merged_df$Patient_Sample %in% id_samples_A[n])))]

                sub_all_clones_A=sub_tmp_group_A[match(unique(merged_df[sub_tmp_group_A,get(my_clone_def)]),
                                                       merged_df[sub_tmp_group_A,get(my_clone_def)])]
                sub_all_clones_A=sub_all_clones_A[min(min_size_A, length(sub_all_clones_A))]
                tmp_tmp_group_A=c(tmp_tmp_group_A,
                                  c(sub_all_clones_A, sample(sub_tmp_group_A[which(merged_df[sub_tmp_group_A,get(my_clone_def)] %!in% merged_df[max_clones_A,get(my_clone_def)])], min_size_A-length(sub_all_clones_A))))
              }
              if( (min_size_AB - min_size_A*n_samples_A) > 0) {
                tmp_tmp_group_A=c(tmp_tmp_group_A,
                                  sample(tmp_group_A[which(tmp_group_A %!in% tmp_tmp_group_A)], min_size_AB - min_size_A*n_samples_A))
              }

              tmp_group_A=tmp_tmp_group_A
            } else {
              all_clones_A=tmp_group_A[match(unique(merged_df[tmp_group_A,get(my_clone_def)]),
                                             merged_df[tmp_group_A,get(my_clone_def)])]
              max_clones_A=sample( all_clones_A, min(min_size_AB, length(all_clones_A)))
              tmp_group_A=c(max_clones_A, sample(tmp_group_A[which(tmp_group_A %!in% max_clones_A)], min_size_AB-length(max_clones_A)))
            }

          } else {

            if(VDJ_normalized_per_sample) {

              tmp_tmp_group_A=c()
              for (n in c(1:n_samples_A)) {

                sub_tmp_group_A=tmp_group_A[which(tmp_group_A %in% intersect(intersect(which(merged_df$V_and_J %in% V_J_comb),which(merged_df[[groups]] %in% group_A) ), which(merged_df$Patient_Sample %in% id_samples_A[n])))]

                tmp_tmp_group_A=c(tmp_tmp_group_A,  sample(sub_tmp_group_A, min_size_A))
              }
              if( (min_size_AB - min_size_A*n_samples_A) > 0) {
                tmp_tmp_group_A=c(tmp_tmp_group_A,
                                  sample(tmp_group_A[which(tmp_group_A %!in% tmp_tmp_group_A)], min_size_AB - min_size_A*n_samples_A))
              }


              tmp_group_A=tmp_tmp_group_A
            } else {
              tmp_group_A=sample(tmp_group_A, min_size_AB)
            }


          }

          max_clones_B=NULL

          if(VDJ_maximize_clones) {


            if(VDJ_normalized_per_sample) {
              max_clones_B=c()
              tmp_tmp_group_B=c()
              for (n in c(1:n_samples_B)) {

                sub_tmp_group_B=tmp_group_B[which(tmp_group_B %in% intersect(intersect(which(merged_df$V_and_J %in% V_J_comb),which(merged_df[[groups]] %in% group_B) ), which(merged_df$Patient_Sample %in% id_samples_B[n])))]

                sub_all_clones_B=sub_tmp_group_B[match(unique(merged_df[sub_tmp_group_B,get(my_clone_def)]),
                                                       merged_df[sub_tmp_group_B,get(my_clone_def)])]
                sub_all_clones_B=sub_all_clones_B[min(min_size_B, length(sub_all_clones_B))]
                tmp_tmp_group_B=c(tmp_tmp_group_B,
                                  c(sub_all_clones_B, sample(sub_tmp_group_B[which(merged_df[sub_tmp_group_B,get(my_clone_def)] %!in% merged_df[max_clones_B,get(my_clone_def)])], min_size_B-length(sub_all_clones_B))))
              }
              if( (min_size_AB - min_size_B*n_samples_B) > 0) {
                tmp_tmp_group_B=c(tmp_tmp_group_B,
                                  sample(tmp_group_B[which(tmp_group_B %!in% tmp_tmp_group_B)], min_size_AB - min_size_B*n_samples_B))
              }

              tmp_group_B=tmp_tmp_group_B
            } else {
              all_clones_B=tmp_group_B[match(unique(merged_df[tmp_group_B,get(my_clone_def)]),
                                             merged_df[tmp_group_B,get(my_clone_def)])]
              max_clones_B=sample( all_clones_B, min(min_size_AB, length(all_clones_B)))
              tmp_group_B=c(max_clones_B, sample(tmp_group_B[which(tmp_group_B %!in% max_clones_B)], min_size_AB-length(max_clones_B)))
            }

          }  else {

            if(VDJ_normalized_per_sample) {

              tmp_tmp_group_B=c()
              for (n in c(1:n_samples_B)) {

                sub_tmp_group_B=tmp_group_B[which(tmp_group_B %in% intersect(intersect(which(merged_df$V_and_J %in% V_J_comb),which(merged_df[[groups]] %in% group_B) ), which(merged_df$Patient_Sample %in% id_samples_B[n])))]

                tmp_tmp_group_B=c(tmp_tmp_group_B,  sample(sub_tmp_group_B, min_size_B))
              }
              if( (min_size_AB - min_size_B*n_samples_B) > 0) {
                tmp_tmp_group_B=c(tmp_tmp_group_B,
                                  sample(tmp_group_B[which(tmp_group_B %!in% tmp_tmp_group_B)], min_size_AB - min_size_B*n_samples_B))
              }


              tmp_group_B=tmp_tmp_group_B


            } else {
              tmp_group_B=sample(tmp_group_B, min_size_AB)
            }


          }


          # if(length(tmp_group_C)>0){
          #   # tmp_index=c(tmp_index, tmp_group_C[1:min(c(length(tmp_group_A),length(tmp_group_B), length(tmp_group_C)))])
          #   max_clones_C=NULL
          #   if(VDJ_maximize_clones) {
          #     max_clones_C=tmp_group_C[match(unique(merged_df[tmp_group_C,which(colnames(merged_df)==my_clone_def)]),
          #                                    merged_df[tmp_group_C,which(colnames(merged_df)==my_clone_def)])]
          #     max_clones_C=sample( tmp_group_C, min_size_ABC)
          #   }
          #   tmp_group_C= c(max_clones_C, sample(tmp_group_C[which(tmp_group_C %!in% max_clones_C)], min_size_ABC))
          #
          # }

          if(length(tmp_group_C)>0){
            max_clones_C=NULL

            if(VDJ_maximize_clones) {


              if(VDJ_normalized_per_sample) {
                max_clones_C=c()
                tmp_tmp_group_C=c()
                for (n in c(1:n_samples_C)) {

                  sub_tmp_group_C=tmp_group_B[which(tmp_group_C %in% intersect(intersect(which(merged_df$V_and_J %in% V_J_comb),which(merged_df[[groups]] %in% group_C) ), which(merged_df$Patient_Sample %in% id_samples_C[n])))]

                  sub_all_clones_C=sub_tmp_group_C[match(unique(merged_df[sub_tmp_group_C,get(my_clone_def)]),
                                                         merged_df[sub_tmp_group_C,get(my_clone_def)])]
                  sub_all_clones_C=sub_all_clones_C[min(min_size_C, length(sub_all_clones_C))]
                  tmp_tmp_group_C=c(tmp_tmp_group_C,
                                    c(sub_all_clones_C, sample(sub_tmp_group_C[which(merged_df[sub_tmp_group_C,get(my_clone_def)] %!in% merged_df[max_clones_C,get(my_clone_def)])], min_size_C-length(sub_all_clones_C))))
                }
                if( (min_size_ABC - min_size_C*n_samples_C) > 0) {
                  tmp_tmp_group_C=c(tmp_tmp_group_C,
                                    sample(tmp_group_C[which(tmp_group_C %!in% tmp_tmp_group_C)], min_size_ABC - min_size_C*n_samples_C))
                }

                tmp_group_C=tmp_tmp_group_C
              } else {
                all_clones_C=tmp_group_C[match(unique(merged_df[tmp_group_C,get(my_clone_def)]),
                                               merged_df[tmp_group_C,get(my_clone_def)])]
                max_clones_C=sample( all_clones_C, min(min_size_AB, length(all_clones_C)))
                tmp_group_C=c(max_clones_C, sample(tmp_group_C[which(tmp_group_C %!in% max_clones_C)], min_size_ABC-length(max_clones_C)))
              }

            }  else {

              if(VDJ_normalized_per_sample) {

                tmp_tmp_group_C=c()
                for (n in c(1:n_samples_C)) {

                  sub_tmp_group_C=tmp_group_C[which(tmp_group_C %in% intersect(intersect(which(merged_df$V_and_J %in% V_J_comb),which(merged_df[[groups]] %in% group_C) ), which(merged_df$Patient_Sample %in% id_samples_C[n])))]

                  tmp_tmp_group_C=c(tmp_tmp_group_C,  sample(sub_tmp_group_C, min_size_C))
                }
                if( (min_size_AB - min_size_C*n_samples_C) > 0) {
                  tmp_tmp_group_C=c(tmp_tmp_group_C,
                                    sample(tmp_group_C[which(tmp_group_C %!in% tmp_tmp_group_C)], min_size_ABC - min_size_C*n_samples_C))
                }


                tmp_group_C=tmp_tmp_group_C


              } else {
                tmp_group_C=sample(tmp_group_C, min_size_ABC)
              }


            }

          } else {
            tmp_group_C=NULL
          }

          # No random
          # tmp_index=c(tmp_index, c(tmp_group_A[1:min(c(length(tmp_group_A),length(tmp_group_B)))], tmp_group_B[1:min(c(length(tmp_group_A),length(tmp_group_B)))]) )

          tmp_index=c(tmp_index, c(tmp_group_A,tmp_group_B,tmp_group_C))

        }



      }



      print(tmp_index)
      rows=rows[which(rows%in%tmp_index)]
    } else {
      rows=rows[which(rows%in%which(merged_df$V_and_J %in% V_J_to_use))]
    }

    print(length(rows))
  }

  if(length(VJ_deselected)!= 0){

    rows=rows[which(rows%in%which(merged_df$V_and_J %!in% VJ_deselected))]
  }





  if(!is.null(my_regions)){
    columns=columns[which(columns %in% which((sapply(merged_header, function(x) strsplit(x, split="_")[[1]][2])) %in% my_regions) )]
  } else {
    columns=c()
  }

  if(!is.null(my_var_elements)){
    columns=columns[which(columns %in% which(grepl(paste0(paste("^", my_var_elements, sep=""), collapse ="|"),merged_header)))]
  } else {
    columns=c()
  }

  if(!is.null(my_vars)){

    merged_header_splitted_at_3=sapply(merged_header, function(x) strsplit(x, split="_")[[1]][3])
    merged_header_splitted_at_4=sapply(merged_header, function(x) strsplit(x, split="_")[[1]][4])
    merged_header_splitted_at_3_and_4=sapply(merged_header, function(x) paste(strsplit(x, split="_")[[1]][3], strsplit(x, split="_")[[1]][4], sep="_"))

    if ("Length" %!in% my_vars) {
      columns=columns[which(columns %!in% which(endsWith(merged_header, "length")))]
      columns=columns[which(columns %!in% which(endsWith(merged_header, "length_Diff")))]
    }

    if ("Composition" %!in% my_vars) {
      nucleotides=c("A","G","T","C")
      peptides=Peptides::aaList()
      columns=columns[which(columns %!in% which(merged_header_splitted_at_3_and_4 %in%
                                                  unique(c(paste(c(nucleotides, peptides), "count", sep="_"), paste(c(nucleotides, peptides), "norm", sep="_") )) ))]


      columns=columns[which(columns %!in% which(merged_header_splitted_at_4 %in% c("counts")   ))]
    }

    if ("Hot/Cold motifs" %!in% my_vars) {
      columns=columns[which(columns %!in% which(merged_header_splitted_at_3 %in% c("hot","cold","potential")   ))]
    }

    if ("Substitutions" %!in% my_vars) {
      columns=columns[which(columns %!in% which(merged_header_splitted_at_3_and_4 %in%
                                                  c("Sub", "Sub_prc", "SIDT_sum", "SID_sum")   ))]
      columns=columns[which(columns %!in% which(merged_header_splitted_at_3 %in%
                                                  c("Sub")   ))]
      columns=columns[which(columns %!in% which(endsWith(merged_header, "_S") ))] ##Old, remove later
      columns=columns[which(columns %!in% which(endsWith(merged_header, "_S_prc") ))] ##Old, remove later
    }

    if ("Insertions" %!in% my_vars) {
      columns=columns[which(columns %!in% which(merged_header_splitted_at_3_and_4 %in%
                                                  c("Ins", "Ins_prc", "SIDT_sum", "SID_sum")   ))]
      columns=columns[which(columns %!in% which(merged_header_splitted_at_3 %in%
                                                  c("Ins")   ))]
      columns=columns[which(columns %!in% which(endsWith(merged_header, "_I") ))] ##Old, remove later
      columns=columns[which(columns %!in% which(endsWith(merged_header, "_I_prc") ))] ##Old, remove later
    }
    if ("Deletions" %!in% my_vars) {
      columns=columns[which(columns %!in% which(merged_header_splitted_at_3_and_4 %in%
                                                  c("Del", "Del_prc", "SIDT_sum", "SID_sum")   ))]
      columns=columns[which(columns %!in% which(merged_header_splitted_at_3 %in%
                                                  c("Del")   ))]
      columns=columns[which(columns %!in% which(endsWith(merged_header, "_D") ))] ##Old, remove later
      columns=columns[which(columns %!in% which(endsWith(merged_header, "_D_prc") ))] ##Old, remove later
    }
    if ("Translocations" %!in% my_vars) {
      columns=columns[which(columns %!in% which(merged_header_splitted_at_3_and_4 %in%
                                                  c("Trasl", "Trasl_prc", "SIDT_sum")   ))]
      columns=columns[which(columns %!in% which(merged_header_splitted_at_3 %in%
                                                  c("Trasl")   ))]
      columns=columns[which(columns %!in% which(endsWith(merged_header, "_T") ))] ##Old, remove later
      columns=columns[which(columns %!in% which(endsWith(merged_header, "_T_prc") ))] ##Old, remove later
    }

    if("Leveshtein distance" %!in% my_vars) {
      columns=columns[which(columns %!in% which(merged_header_splitted_at_3 %in% c("lv")   ))]
    }

    if ("Transitions and transversions" %!in% my_vars) {
      columns=columns[which(columns %!in% which(merged_header_splitted_at_3 %in% c("Transitions","Transversions")   ))]
      columns=columns[which(columns %!in% which(merged_header_splitted_at_3_and_4 %in%
                                                  c("Ratio_Transitions-Transversions")   ))]
    }
    if ("Replacement and silent mutations" %!in% my_vars) {
      columns=columns[which(columns %!in% which(merged_header_splitted_at_3 %in% c("Replacement","Silent")   ))]
      columns=columns[which(columns %!in% which(merged_header_splitted_at_3_and_4 %in%
                                                  c("Ratio_Silent-Replacement")   ))]
    }

    if ("Mutations from X to Y" %!in% my_vars) {
      columns=columns[which(columns %!in% which(merged_header_splitted_at_4 %in% c("to")   ))]
      # columns=columns[which(columns %!in% which(sapply(merged_header, function(x) paste(strsplit(x, split="_")[[1]][3], strsplit(x, split="_")[[1]][4], sep="_")) %in%
      #                                             c("Ratio_Silent")   ))]
    }

    if ("NGly sites" %!in% my_vars) {
      columns=columns[which(columns %!in% which(merged_header_splitted_at_3 %in% c("NGly")   ))]
    }
    if ("Peptide features" %!in% my_vars) {
      columns=columns[which(columns %!in% which(merged_header_splitted_at_3 %in% c("Peptides", "alkzm",
                                                                                   "Small","Tiny", "Aliphatic",
                                                                                   "Charged","Polar", "Basic","NonPolar",
                                                                                   "Aromatic", "Acidic")   ))]

      # columns=columns[which(columns %!in% which(sapply(merged_header, function(x) strsplit(x, split="_")[[1]][3]) %in% c("alkzm")   ))]
      # columns=columns[which(columns %!in% which(sapply(merged_header, function(x) strsplit(x, split="_")[[1]][3]) %in% c("alkzm")   ))]
    }


  } else {
    columns=c()
  }

  if(!is.null(my_vartypes)){
    if("Germline diff" %!in% my_vartypes) {
      columns=columns[which(columns %!in% which(endsWith(merged_header, "_Diff")))]
    }
    if("Baseline" %!in% my_vartypes) {
      columns=columns[which(columns %in% which(endsWith(merged_header, "_Diff")))]
    }
  } else {
    columns=c()
  }

  if(!is.null(variables_to_remove)) {
    columns=columns[which(columns %!in% which(merged_header %in% variables_to_remove))]
  }

  print("Big_summary for")
  print(length(merged_header[columns]))

  big_summary=bigstatsr::big_colstats(copy_FBM, ind.row = rows, ind.col=columns, ncores=nb_cores())
  vars_to_remove=columns[unique(c(which(big_summary$sum ==0), which(big_summary$var ==0),
                                  which(is.na(big_summary$sum)),which(is.infinite(big_summary$sum)),which(is.nan(big_summary$sum)),
                                  which(is.na(big_summary$var)),which(is.infinite(big_summary$var)),which(is.nan(sqrt(big_summary$var)))
  ))]

  # if(any(is.nan(sqrt(big_summary$var)))) {
  #   print(copy_FBM[rows, columns[which(is.nan(sqrt(big_summary$var)))]])
  # }
  columns=columns[which(columns%!in%vars_to_remove)]

  if(univlog && !is.null(columns) && !is.null(group_B) && !is.null(group_A) && to_compare_groups){

    merged_df$Binary_classif=0
    merged_df$Binary_classif[which(merged_df[[groups]] %in% group_B)]=1
    tmp_rows=rows[(which(rows %in% which(merged_df[[groups]] %in% c(group_B, group_A))))]
    print(table(merged_df$Binary_classif[tmp_rows]))
    testuniv <- big_univLogReg(copy_FBM, merged_df$Binary_classif[tmp_rows], ind.train = tmp_rows,
                               ind.col =columns,
                               # covar.train = covar[covar_training, ],
                               ncores=nb_cores())
    testuniv$p.value=predict(testuniv, log10=FALSE)
    print(min(testuniv$p.value))

    if (pval_type =="Corrected by B-H") {
      testuniv$p.value=round(p.adjust(testuniv$p.value, "BH"), 6)
    } else if(pval_type == "Corrected by Bonferroni") {
      testuniv$p.value=round(p.adjust(testuniv$p.value, "bonferroni"), 6)
    }

    print("P.values")

    print(min(testuniv$p.value))
    print(class(pval_cutoff))
    print(class(estimate_cutoff))
    if (number_selected_vars == "All") {
      columns=columns[intersect(which(abs(testuniv$p.value) <= pval_cutoff), which(abs(testuniv$estim) >estimate_cutoff))]
    } else if (number_selected_vars == "Top 10 (according to estimate)" ) {
      intersection=intersect(which(abs(testuniv$p.value) <= pval_cutoff), which(abs(testuniv$estim) >estimate_cutoff))
      columns=columns[intersection[sort(abs(testuniv$estim)[intersection], index.return=TRUE, decreasing=TRUE)$ix[1:min(10, length(intersection))]]]

    } else if(number_selected_vars == "Top 50 (according to estimate)") {
      intersection=intersect(which(abs(testuniv$p.value) < pval_cutoff), which(abs(testuniv$estim) >estimate_cutoff))
      columns=columns[intersection[sort(abs(testuniv$estim)[intersection], index.return=TRUE, decreasing=TRUE)$ix[1:min(50, length(intersection))]]]
    }

  }

  print("Filter passed by: ")
  print(length(rows))
  print(length(columns))


  return(list(ROWS=rows, COLUMNS=columns))
}

#' 3_Data_processing
#'
#' @description A fct function
#'
#' @return The return value, if any, from executing the function.
#' @import bigstatsr
#' @import stats
#' @noRd
big_PCA_plot=function(copy_FBM, rows, columns){

  if(length(columns)>0 && length(rows)>0) {

    testing= bigstatsr::big_randomSVD(copy_FBM,
                                      fun.scaling = big_scale(),
                                      ind.row = rows,
                                      ind.col = columns,
                                      k = 5,
                                      verbose=T)

    scores <- predict(testing)

    PCA_error_perhaps=tryCatch({

      var_exp <- testing$d^2 / big_norm(copy_FBM,ind.row=rows,
                                                   ind.col=columns,
                                                   center = testing$center,
                                                   scale = testing$scale)
      variance=signif(var_exp, 2)

      cumulative_variance=round(cumsum(var_exp), 3)

      list(Scores=scores, Variance_explained=variance)
    },error = function(e){
      print('Caught an error!')
      list(Scores=NULL, Variance_explained=NULL)
    },
    warning = function(w){
      print('Caught an warning!')
      list(Scores=NULL, Variance_explained=NULL)
    })
    return(list(Scores=PCA_error_perhaps$Scores,
                Variance_explained=PCA_error_perhaps$Variance_explained))

  } else {
    return (list(Scores=NULL, Variance_explained=NULL))
  }


}
