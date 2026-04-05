#' translate_fun
#'
#' @description A function to translate NT sequences into AA sequences and select
#' the optimal translation (if it exists). Prioritizes mclosest ORF to 1.
#'
#' @return A list with either the optimal-without-stop-codon translation and its
#' ORF, or a list indicating there is no a correct translation.
#' @noRd
translate_fun <- function(x, ORF_begins) {

    if (nchar(x)<=2) {
        return(list("", ""))
    } else {
      seq_to_translate= Biostrings::subseq(
         Biostrings::DNAString(as.character(x)),
        start = if(ORF_begins==1){1}else if(ORF_begins==2){3} else if(ORF_begins==0) {2})

      AA_seq=Biostrings::translate(seq_to_translate,
                                   genetic.code= Biostrings::GENETIC_CODE,
                                   no.init.codon=T)

      is_there_stop=sum(length(Biostrings::matchPattern("*", AA_seq)))
      AA_seq=as.character(AA_seq)
      if(AA_seq == ""){
            return(list("", ""))
        } else if (is_there_stop == 0) {
            return(list(AA_seq, if(ORF_begins==1){1}else if(ORF_begins==2){3} else if(ORF_begins==0) {2}))
        } else {
            return(list("NO_GOOD_ORF", 0))  #If the ORF has stop codons. This value is used for filtering afterwards.
        }
        # ORFs <- lapply(
        #     1:3, function(pos) Biostrings::subseq(
        #         Biostrings::DNAString(as.character(x)),
        #         start = pos
        #     )
        # )
        # suppressWarnings(
        #     translated_ORFs <- lapply(
        #         ORFs, Biostrings::translate, genetic.code = Biostrings::GENETIC_CODE,
        #         no.init.codon = TRUE, if.fuzzy.codon = "error"
        #     )
        # )
        #
        # # Search for the best ORF
        # Stop_per_ORF <- vector()
        # for (i in 1:3) {
        #     Stop_per_ORF[i] <- sum(length(Biostrings::matchPattern("*", unlist(translated_ORFs[[i]]))))
        # }
        #
        # # Find STOP codons per ORF and choose the ORF with the less
        # BestORF <- which.min(Stop_per_ORF)
        # aa_seq <- as.character(translated_ORFs[[BestORF]])
        #
        # if(aa_seq == ""){
        #     return(list("", ""))
        # } else if (Stop_per_ORF[BestORF] == 0) {
        #     return(list(aa_seq, BestORF))
        # } else {
        #     return(list("NO_GOOD_ORF", 0))  #If every ORF has stop codons. This value is used for filtering afterwards.
        # }
    }

}

#' Peptides_aaCheck
#'
#' @description Original code from package Peptides, internal function aaCheck.
#' Original authors: Daniel Osorio ORCID iD [aut, cre], Paola Rondon-Villarreal ORCID iD [aut, ths], Rodrigo Torres ORCID iD [aut, ths], J. Sebastian Paez ORCID iD [ctb], Luis Pedro Coelho ORCID iD [ctb], Richèl J.C. Bilderbeek ORCID iD [ctb], Florian C. Sigloch ORCID iD [ctb]
#' Original source: https://github.com/dosorio/Peptides/blob/master/R/aaCheck.R
#' Original citation: Osorio, D., Rondon-Villarreal, P. & Torres, R. Peptides: A package for data mining of antimicrobial peptides. The R Journal. 7(1), 4–14 (2015).
#' Original license: GPL-2 (see https://cran.r-project.org/package=Peptides)
#' @return A list of aminoacids in uppercase
#' @noRd
Peptides_aaCheck <- function (seq)
{
  if (!any(lengths(seq) > 1)) {
    seq <- toupper(seq)
    seq <- gsub(pattern = "[[:space:]]+", replacement = "",
                x = seq)
    seq <- strsplit(x = seq, split = "")
  }
  else {
    seq <- lapply(seq, toupper)
  }
  check <- unlist(lapply(seq, function(sequence) {
    !all(sequence %in% c("A", "C", "D", "E", "F", "G", "H",
                         "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T",
                         "V", "W", "Y", "-"))
  }))
  if (sum(check) > 0) {
    sapply(which(check == TRUE), function(sequence) {
      warning(paste0("Sequence ", sequence, " has unrecognized amino acid types. Output value might be wrong calculated"),
              call. = FALSE)
    })
  }
  return(seq)
}

#' 2_Feature_determination
#'
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#' @import stringr
#' @noRd
subseq_func <- function(NT_seq, size_NT_regions, FWR1partial = F, FWR4partial = F, ORF_begins) {
    ########### TOCAR
    AA_info <- translate_fun(toupper(NT_seq), ORF_begins)
    AA_seq <- AA_info[[1]]
    seq_list <- list()
    ORF <- AA_info[[2]]

    if ((AA_seq == "NO_GOOD_ORF") == TRUE) {
        seq_list[1:(length(size_NT_regions) +
            2)] <- "NO_GOOD_ORF"
    } else {
        # initialize
        if (FWR1partial) {
            size_NT_regions <- size_NT_regions[2:length(size_NT_regions)]
        }
        if (FWR4partial) {
            size_NT_regions <- size_NT_regions[1:c(
                length(size_NT_regions) -
                  1
            )]
        }
        size_NT_regions[1] <- size_NT_regions[1] - ORF + 1

        size_AA_regions=c()
        off_set=0
        for(i in c(1:(length(size_NT_regions)))){
          if((size_NT_regions[i]-off_set)%%3 == 0) {
            size_AA_regions=c(size_AA_regions,( size_NT_regions[i]-off_set)/3)
            off_set=0
          } else {
            if(i != (length(size_NT_regions))) {
              onset=3-(size_NT_regions[i]-off_set)%%3

              if(nchar(size_NT_regions[i+1])>= onset) {
                if(onset==2){
                  size_AA_regions=c(size_AA_regions, floor((size_NT_regions[i]-off_set)/3))
                  off_set=-1
                } else {
                  size_AA_regions=c(size_AA_regions, ceiling((size_NT_regions[i]-off_set)/3))
                  off_set=onset
                }

              } else {
                size_AA_regions=c(size_AA_regions, floor(size_NT_regions[i]-off_set)/3)
                off_set=0
              }
            } else {
              size_AA_regions=c(size_AA_regions, floor(size_NT_regions[i]-off_set)/3)
              off_set=0
            }


          }
        }
        # size_AA_regions <- round(size_NT_regions/3)
        # size_AA_regions[length(size_AA_regions) -
        #     1] <- floor(
        #     (size_NT_regions[length(size_AA_regions) -
        #         1]/3)
        # )
        # size_AA_regions[length(size_AA_regions)] <- floor((size_NT_regions[length(size_AA_regions)]/3))
        end <- 0


        for (i in 1:(length(size_AA_regions) -
            1)) {
            start <- end + 1
            end <- start - 1 + size_AA_regions[i]
            seq_list[i] <- Biostrings::subseq(AA_seq, start, end)
        }

        seq_list[[length(size_AA_regions)]] <- Biostrings::subseq(
            AA_seq, sum(
                size_AA_regions[1:length(size_AA_regions) -
                  1]
            ) +
                1, stringr::str_length(AA_seq)
        )

        if (FWR1partial) {
            seq_list <- append(seq_list, "", after = 0)
        }

        if (FWR4partial) {
            seq_list <- append(seq_list, "", after = length(seq_list))
        }


        if(paste(unlist(seq_list), collapse="") != AA_seq){
          strange=T
        } else {
          strange=F
        }
        seq_list[length(seq_list) +
            1] <- AA_seq
        seq_list[length(seq_list) +
            1] <- ORF

        if(strange){
          seq_list[length(seq_list) +
                     1] <- ("MMMMMMMMMMMmmmmmmmmmmm strange")
        }




    }
    return(seq_list)
}

#' 2_Feature_determination
#'
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#'
#' @noRd
length_of_region <- function(sequence, region) {
    region_length <- nchar(sequence)
    names(region_length) <- paste0(region, "_length")
    return(region_length)
}

#' 2_Feature_determination
#'
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#' @import stringdist
#' @noRd
lv_changes <- function(seq_repertoire, seq_rec.germline, region) {
    lv_dist <- stringdist::stringdist(seq_rec.germline, seq_repertoire, method = "lv")
    if (Biostrings::width(as.character(seq_repertoire)) ==
        0) {
        lv_dist_norm <- 0
    } else {
        lv_dist_norm <- (lv_dist/Biostrings::width(as.character(seq_repertoire))) *
            100
    }

    changes <- c(lv_dist, lv_dist_norm)
    names(changes) <- c("lv_dist", "lv_dist_norm")
    names(changes) <- paste(
        as.character(region),
        names(changes),
        sep = "_"
    )
    return(changes)
}

#' 2_Feature_determination
#'
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#' @import stringdist
#' @import utils
#' @noRd
nt_changes <- function(seq_repertoire, seq_rec.germline, region) {
    # lv_info <- lv_changes(seq_repertoire, seq_rec.germline, region)
    # lv_dist <- lv_info[1]
    # osa_dist <- stringdist::stringdist(seq_rec.germline, seq_repertoire, method = "osa")
    # no.trans <- lv_dist - osa_dist
    # tab <- drop(
    #     attr(
    #         adist(seq_rec.germline, seq_repertoire, count = TRUE),
    #         "counts"
    #     )
    # )
    # no.ins <- tab[["ins"]] - no.trans
    # no.del <- tab[["del"]] - no.trans
    # no.subs <- tab[["sub"]]
    # changes <- c(no.subs, no.ins, no.del, no.trans)
    # if (sum(changes) !=
    #     0) {
    #     changes_plus <- c(
    #         changes, sum(changes),
    #         (changes[1]/sum(changes)) *
    #             100, (changes[2]/sum(changes)) *
    #             100, (changes[3]/sum(changes)) *
    #             100, (changes[4]/sum(changes)) *
    #             100
    #     )
    # } else {
    #     changes_plus <- c(
    #         changes, sum(changes),
    #         0, 0, 0, 0
    #     )
    # }
    # names(changes_plus) <- c(
    #     "Sub", "Ins", "Del", "Trasl", "SIDT_sum", "Sub_prc", "Ins_prc", "Del_prc",
    #     "Trasl_prc"
    # )
    # names(changes_plus) <- gsub(
    #     "^", paste0(
    #         as.character(region),
    #         "_"
    #     ),
    #     names(changes_plus)
    # )

  if(nchar(seq_repertoire) > 0 ) {
    insertions=as.data.frame(str_locate_all(seq_repertoire,"[a-z]+")[[1]])
    n_ins=nrow(insertions)
    if(nrow(insertions) >0) {
      insertions$length=insertions$end-insertions$start+1
      length_ins=sum(insertions$length)
    } else {
      length_ins=0
    }

    deletions=as.data.frame(str_locate_all(seq_rec.germline,"[a-z]+")[[1]])
    n_dels=nrow(deletions)
    if(nrow(deletions) >0) {
      deletions$length=deletions$end-deletions$start+1
      length_del=sum(deletions$length)
    } else {
      length_del=0
    }

  } else {
    n_ins=0
    n_dels=0
    length_ins=0
    length_del=0
  }

    changes=c(n_ins, n_dels,length_ins,length_del)
    names(changes)=c("NumIns","NumDels","LenghtIns","LenghtDels")

    names(changes) <- gsub(
        "^", paste0(
            as.character(region),
            "_"
        ),
        names(changes)
    )
    return(changes)
}

#' 2_Feature_determination
#'
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#' @noRd
show_selected_features <- function(tmp_sunburst = tmp_sunburst, region = region) {

    # for (region in regions) {
    index_regions <- which(grepl(region, tmp_sunburst, fixed = T))
    tmp_sunburst <- tmp_sunburst[index_regions]
    tmp_sunburst <- gsub(
        paste(region, "_", sep = ""),
        "", tmp_sunburst, fixed = T
    )

    # tmp_sunburst[index_regions]=paste(region, tmp_sunburst[index_regions],
    # sep='-') }


    tmp_sunburst <- gsub("NT_", "NT-", tmp_sunburst, fixed = T)
    tmp_sunburst <- gsub("AA_", "AA-", tmp_sunburst, fixed = T)
    tmp_sunburst <- gsub("_Diff", "-Diff_with_germline", tmp_sunburst, fixed = T)
    tmp_index <- which(!grepl("-Diff_with_germline", tmp_sunburst))
    tmp_sunburst[tmp_index] <- paste(tmp_sunburst[tmp_index], "-Sequence_values", sep = "")




    nucleotides <- c("A", "G", "T", "C")
    peptides <- Peptides::aaList()

    tmp_index <- which(
        grepl(
            paste(
                c(
                  paste(
                    "^NT-", paste(
                      c(nucleotides),
                      "count", sep = "_"
                  ),
                    sep = ""
                ),
                  paste(
                    "^NT-", paste(
                      c(nucleotides),
                      "norm", sep = "_"
                  ),
                    sep = ""
                ),
                  paste(
                    "^AA-", paste(
                      c(peptides),
                      "count", sep = "_"
                  ),
                    sep = ""
                ),
                  paste(
                    "^AA-", paste(
                      c(peptides),
                      "norm", sep = "_"
                  ),
                    sep = ""
                ),
                  "counts"
              ),
                collapse = "|"
            ),
            tmp_sunburst
        )
    )
    tmp_sunburst[tmp_index] <- gsub("NT-", "NT-Composition-", tmp_sunburst[tmp_index])
    tmp_sunburst[tmp_index] <- gsub("^AA-", "AA-Composition-Aminoacids-", tmp_sunburst[tmp_index], perl = T)
    tmp_index <- which(grepl("counts", tmp_sunburst))
    tmp_sunburst[tmp_index] <- gsub("^AA-", "AA-Composition-Codons-", tmp_sunburst[tmp_index], perl = T)
    # print(tmp_sunburst[tmp_index])

    tmp_index <- which(
        grepl(
            paste(
                c("hot", "cold", "potential"),
                collapse = "|"
            ),
            tmp_sunburst
        )
    )
    tmp_sunburst[tmp_index] <- gsub("cold_", "cold-", tmp_sunburst[tmp_index], fixed = T)
    tmp_sunburst[tmp_index] <- gsub("hot_", "hot-", tmp_sunburst[tmp_index], fixed = T)
    tmp_sunburst[tmp_index] <- gsub("_R_", "-R_", tmp_sunburst[tmp_index], fixed = T)
    tmp_sunburst[tmp_index] <- gsub("_F_", "-F_", tmp_sunburst[tmp_index], fixed = T)
    tmp_sunburst[tmp_index] <- gsub("_norm", "-norm", tmp_sunburst[tmp_index], fixed = T)
    tmp_sunburst[tmp_index] <- gsub("_count", "-count", tmp_sunburst[tmp_index], fixed = T)
    tmp_sunburst[tmp_index] <- gsub("NT-", "NT-Hot/Cold motifs-", tmp_sunburst[tmp_index])
    tmp_sunburst[tmp_index] <- gsub("^AA-", "AA-Hot/Cold motifs-", tmp_sunburst[tmp_index], perl = T)

    tmp_index <- which(
        grepl(
            paste(
                c("Sub", "Sub_prc", "SIDT", "SID_sum"),
                collapse = "|"
            ),
            tmp_sunburst
        )
    )
    tmp_sunburst[tmp_index] <- gsub("NT-", "NT-Substitutions-", tmp_sunburst[tmp_index])
    tmp_sunburst[tmp_index] <- gsub("AA-", "AA-Substitutions-", tmp_sunburst[tmp_index])


    tmp_index <- which(
        grepl(
            paste(
                c("Ins", "Ins_prc", "SIDT", "SID_sum"),
                collapse = "|"
            ),
            tmp_sunburst
        )
    )
    tmp_sunburst[tmp_index] <- gsub("NT-", "NT-Insertions-", tmp_sunburst[tmp_index])
    tmp_sunburst[tmp_index] <- gsub("AA-", "AA-Insertions-", tmp_sunburst[tmp_index])

    tmp_index <- which(
        grepl(
            paste(
                c("Del", "Del_prc", "SIDT", "SID_sum"),
                collapse = "|"
            ),
            tmp_sunburst
        )
    )
    tmp_sunburst[tmp_index] <- gsub("NT-", "NT-Deletions-", tmp_sunburst[tmp_index])
    tmp_sunburst[tmp_index] <- gsub("AA-", "AA-Deletions-", tmp_sunburst[tmp_index])


    tmp_index <- which(
        grepl(
            paste(
                c("Trasl", "Trasl_prc", "SIDT"),
                collapse = "|"
            ),
            tmp_sunburst
        )
    )
    tmp_sunburst[tmp_index] <- gsub("NT-", "NT-Translocations-", tmp_sunburst[tmp_index])
    tmp_sunburst[tmp_index] <- gsub("AA-", "AA-Translocations-", tmp_sunburst[tmp_index])


    tmp_index <- which(grepl("lv_dist", tmp_sunburst))
    tmp_sunburst[tmp_index] <- gsub(
        "Transitions-Transversions", "TransitionsToTransversions", tmp_sunburst[tmp_index],
        fixed = T
    )
    tmp_sunburst[tmp_index] <- gsub("_norm", "-norm", tmp_sunburst[tmp_index], fixed = T)
    tmp_sunburst[tmp_index] <- gsub("_count", "-count", tmp_sunburst[tmp_index], fixed = T)
    tmp_sunburst[tmp_index] <- gsub("NT-", "NT-Leveshtein distance-", tmp_sunburst[tmp_index])
    tmp_sunburst[tmp_index] <- gsub("AA-", "AA-Leveshtein distance-", tmp_sunburst[tmp_index])


    tmp_index <- which(grepl("Transitions|Transversions", tmp_sunburst))
    tmp_sunburst[tmp_index] <- gsub(
        "Transitions-Transversions", "TransitionsToTransversions", tmp_sunburst[tmp_index],
        fixed = T
    )
    tmp_sunburst[tmp_index] <- gsub("_norm", "-norm", tmp_sunburst[tmp_index], fixed = T)
    tmp_sunburst[tmp_index] <- gsub("_count", "-count", tmp_sunburst[tmp_index], fixed = T)
    tmp_sunburst[tmp_index] <- gsub("NT-", "NT-Transitions and transversions-", tmp_sunburst[tmp_index])
    tmp_sunburst[tmp_index] <- gsub("AA-", "AA-Transitions and transversions-", tmp_sunburst[tmp_index])

    tmp_index <- which(grepl("Replacement|Silent", tmp_sunburst))
    tmp_sunburst[tmp_index] <- gsub(
        "Silent-Replacement", "SilentToReplacement", tmp_sunburst[tmp_index], fixed = T
    )
    tmp_sunburst[tmp_index] <- gsub("_norm", "-norm", tmp_sunburst[tmp_index], fixed = T)
    tmp_sunburst[tmp_index] <- gsub("NT-", "NT-Replacement and silent mutations-", tmp_sunburst[tmp_index])
    tmp_sunburst[tmp_index] <- gsub("AA-", "AA-Replacement and silent mutations-", tmp_sunburst[tmp_index])

    tmp_index <- which(grepl("_to_", tmp_sunburst))
    tmp_sunburst[tmp_index] <- gsub("_count", "-count", tmp_sunburst[tmp_index], fixed = T)
    tmp_sunburst[tmp_index] <- gsub("_norm", "-norm", tmp_sunburst[tmp_index], fixed = T)
    tmp_sunburst[tmp_index] <- gsub("NT-", "NT-Mutations, from A to B-", tmp_sunburst[tmp_index])
    tmp_sunburst[tmp_index] <- gsub("AA-", "AA-Mutations, from A to B-", tmp_sunburst[tmp_index])


    tmp_index <- which(
        grepl(
            paste(
                c(
                  "Peptides", "alkzm", "Small", "Tiny", "Aliphatic", "Charged", "Polar",
                  "Basic", "NonPolar", "Aromatic", "Acidic"
              ),
                collapse = "|"
            ),
            tmp_sunburst
        )
    )
    tmp_sunburst[tmp_index] <- gsub("Peptides_", "Peptides-", tmp_sunburst[tmp_index], fixed = T)
    tmp_sunburst[tmp_index] <- gsub("alkzm_", "alkzm-", tmp_sunburst[tmp_index], fixed = T)
    tmp_sunburst[tmp_index] <- gsub("_count", "-count", tmp_sunburst[tmp_index], fixed = T)
    tmp_sunburst[tmp_index] <- gsub("_norm", "-norm", tmp_sunburst[tmp_index], fixed = T)
    tmp_sunburst[tmp_index] <- gsub("NT-", "NT-Peptide features-", tmp_sunburst[tmp_index])
    tmp_sunburst[tmp_index] <- gsub("AA-", "AA-Peptide features-", tmp_sunburst[tmp_index])


    sunburst_df <- data.frame(V1 = tmp_sunburst)
    sunburst_df$V2 <- 1

    return(sunburst_df)
}

#' 2_Feature_determination
#'
#' @description A utils function, not used currently
#'
#' @return The return value, if any, from executing the utility.
#' @importFrom utils adist
#' @noRd
#'
aa_changes <- function(seq_repertoire, seq_rec.germline, region) {
    tab <- drop(
        attr(
            adist(seq_rec.germline, seq_repertoire, counts = TRUE),
            "counts"
        )
    )
    no.ins <- tab[["ins"]]
    no.del <- tab[["del"]]
    no.subs <- tab[["sub"]]
    changes <- c(no.subs, no.ins, no.del)
    if (sum(changes) !=
        0) {
        changes_plus <- c(
            changes, sum(changes),
            (changes[1]/sum(changes)) *
                100, (changes[2]/sum(changes)) *
                100, (changes[3]/sum(changes)) *
                100
        )
    } else {
        changes_plus <- c(
            changes, sum(changes),
            0, 0, 0
        )
    }
    names(changes_plus) <- c("Sub", "Ins", "Del", "SID_sum", "Sub_prc", "Ins_prc", "Del_prc")
    names(changes_plus) <- gsub(
        "^", paste0(
            as.character(region),
            "_"
        ),
        names(changes_plus)
    )

    lv_info <- lv_changes(seq_repertoire, seq_rec.germline, region)
    return(c(changes_plus, lv_info))
}

#' 2_Feature_determination
#'
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#' @import stringr
#' @noRd
count_elements <- function(sequence, region) {
    if (startsWith(region, "NT")) {
        elements <- c("A", "G", "C", "T")
    } else if (startsWith(region, "AA")) {
        aa_names <- c(
            "Phe", "Met", "Leu", "Ile", "Val", "Pro", "Tyr", "Trp", "Cys", "Ala",
            "Gly", "Ser", "Thr", "His", "Glu", "Gln", "Asp", "Asn", "Lys", "Arg"
        )
        elements <- seqinr::a(aa_names)
    } else {
        # print("[ERROR] Wrong region included")
    }
    length_seq <- stringr::str_length(sequence)
    results <- c()

    for (i in elements) {
        counts <- str_count(sequence, i)
        previous_names <- names(results)
        if (length_seq == 0) {
            results <- c(results, counts, 0)
        } else {
            results <- c(results, counts, 100 * counts/length_seq)
        }

        names(results) <- c(
            previous_names, paste0(region, "_", i, "_count"),
            paste0(region, "_", i, "_norm")
        )
    }
    return(results)
}

#' 2_Feature_determination
#'
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#' @import Peptides
#' @noRd
instaIndex_improved <- function(seq) {
    #### fake instability index of 0 to single aminoacids changed to be able to
    #### work with dipeptides
    guruprasad <- c(
        WW = 1, WC = 1, WM = 24.68, WH = 24.68, WY = 1, WF = 1, WQ = 1, WN = 13.34,
        WI = 1, WR = 1, WD = 1, WP = 1, WT = -14.03, WK = 1, WE = 1, WV = -7.49,
        WS = 1, WG = -9.37, WA = -14.03, WL = 13.34, CW = 24.68, CC = 1, CM = 33.6,
        CH = 33.6, CY = 1, CF = 1, CQ = -6.54, CN = 1, CI = 1, CR = 1, CD = 20.26,
        CP = 20.26, CT = 33.6, CK = 1, CE = 1, CV = -6.54, CS = 1, CG = 1, CA = 1,
        CL = 20.26, MW = 1, MC = 1, MM = -1.88, MH = 58.28, MY = 24.68, MF = 1, MQ = -6.54,
        MN = 1, MI = 1, MR = -6.54, MD = 1, MP = 44.94, MT = -1.88, MK = 1, ME = 1,
        MV = 1, MS = 44.94, MG = 1, MA = 13.34, ML = 1, HW = -1.88, HC = 1, HM = 1,
        HH = 1, HY = 44.94, HF = -9.37, HQ = 1, HN = 24.68, HI = 44.94, HR = 1, HD = 1,
        HP = -1.88, HT = -6.54, HK = 24.68, HE = 1, HV = 1, HS = 1, HG = -9.37, HA = 1,
        HL = 1, YW = -9.37, YC = 1, YM = 44.94, YH = 13.34, YY = 13.34, YF = 1, YQ = 1,
        YN = 1, YI = 1, YR = -15.91, YD = 24.68, YP = 13.34, YT = -7.49, YK = 1,
        YE = -6.54, YV = 1, YS = 1, YG = -7.49, YA = 24.68, YL = 1, FW = 1, FC = 1,
        FM = 1, FH = 1, FY = 33.6, FF = 1, FQ = 1, FN = 1, FI = 1, FR = 1, FD = 13.34,
        FP = 20.26, FT = 1, FK = -14.03, FE = 1, FV = 1, FS = 1, FG = 1, FA = 1,
        FL = 1, QW = 1, QC = -6.54, QM = 1, QH = 1, QY = -6.54, QF = -6.54, QQ = 20.26,
        QN = 1, QI = 1, QR = 1, QD = 20.26, QP = 20.26, QT = 1, QK = 1, QE = 20.26,
        QV = -6.54, QS = 44.94, QG = 1, QA = 1, QL = 1, NW = -9.37, NC = -1.88, NM = 1,
        NH = 1, NY = 1, NF = -14.03, NQ = -6.54, NN = 1, NI = 44.94, NR = 1, ND = 1,
        NP = -1.88, NT = -7.49, NK = 24.68, NE = 1, NV = 1, NS = 1, NG = -14.03,
        NL = 1, IW = 1, IC = 1, IM = 1, IH = 13.34, IY = 1, IF = 1, IQ = 1, IN = 1,
        II = 1, IR = 1, ID = 1, IP = -1.88, IT = 1, IK = -7.49, IE = 44.94, IV = -7.49,
        IS = 1, IG = 1, IA = 1, IL = 20.26, RW = 58.28, RC = 1, RM = 1, RH = 20.26,
        RY = -6.54, RF = 1, RQ = 20.26, RN = 13.34, RI = 1, RR = 58.28, RD = 1, RP = 20.26,
        RT = 1, RK = 1, RE = 1, RV = 1, RS = 44.94, RG = -7.49, RA = 1, RL = 1, DW = 1,
        DC = 1, DM = 1, DH = 1, DY = 1, DF = -6.54, DQ = 1, DN = 1, DI = 1, DR = -6.54,
        DD = 1, DP = 1, DT = -14.03, DK = -7.49, DE = 1, DV = 1, DS = 20.26, DG = 1,
        DA = 1, DL = 1, PW = -1.88, PC = -6.54, PM = -6.54, PH = 1, PY = 1, PF = 20.26,
        PQ = 20.26, PN = 1, PI = 1, PR = -6.54, PD = -6.54, PP = 20.26, PT = 1, PK = 1,
        PE = 18.38, PV = 20.26, PS = 20.26, PG = 1, PA = 20.26, PL = 1, TW = -14.03,
        TC = 1, TM = 1, TH = 1, TY = 1, TF = 13.34, TQ = -6.54, TN = -14.03, TI = 1,
        TR = 1, TD = 1, TP = 1, TT = 1, TK = 1, TE = 20.26, TV = 1, TS = 1, TG = -7.49,
        TA = 1, TL = 1, KW = 1, KC = 1, KM = 33.6, KH = 1, KY = 1, KF = 1, KQ = 24.68,
        KN = 1, KI = -7.49, KR = 33.6, KD = 1, KP = -6.54, KT = 1, KK = 1, KE = 1,
        KV = -7.49, KS = 1, KG = -7.49, KA = 1, KL = -7.49, EW = -14.03, EC = 44.94,
        EM = 1, EH = -6.54, EY = 1, EF = 1, EQ = 20.26, EN = 1, EI = 20.26, ER = 1,
        ED = 20.26, EP = 20.26, ET = 1, EK = 1, EE = 33.6, EV = 1, ES = 20.26, EG = 1,
        EA = 1, EL = 1, VW = 1, VC = 1, VM = 1, VH = 1, VY = -6.54, VF = 1, VQ = 1,
        VN = 1, VI = 1, VR = 1, VD = -14.03, VP = 20.26, VT = -7.49, VK = -1.88,
        VE = 1, VV = 1, VS = 1, VG = -7.49, VA = 1, VL = 1, SW = 1, SC = 33.6, SM = 1,
        SH = 1, SY = 1, SF = 1, SQ = 20.26, SN = 1, SI = 1, SR = 20.26, SD = 1, SP = 44.94,
        ST = 1, SK = 1, SE = 20.26, SV = 1, SS = 20.26, SG = 1, SA = 1, SL = 1, GW = 13.34,
        GC = 1, GM = 1, GH = 1, GY = -7.49, GF = 1, GQ = 1, GN = -7.49, GI = -7.49,
        GR = 1, GD = 1, GP = 1, GT = -7.49, GK = -7.49, GE = -6.54, GV = 1, GS = 1,
        GG = 13.34, GA = -7.49, GL = 1, AW = 1, AC = 44.94, AM = 1, AH = -7.49, AY = 1,
        AF = 1, AQ = 1, AN = 1, AI = 1, AR = 1, AD = -7.49, AP = 20.26, AT = 1, AK = 1,
        AE = 1, AV = 1, AS = 1, AG = 1, AA = 1, AL = 1, LW = 24.68, LC = 1, LM = 1,
        LH = 1, LY = 1, LF = 1, LQ = 33.6, LN = 1, LI = 1, LR = 20.26, LD = 1, LP = 20.26,
        LT = 1, LK = -7.49, LE = 1, LV = 1, LS = 1, LG = 1, LA = 1, LL = 1, `NA` = 1
    )
    aa <- Peptides_aaCheck(seq)
    dp <- lapply(
        aa, function(aa) {
            if (length(aa) ==
                1) {
                0
            } else if (length(aa) ==
                2) {
                unname(
                  apply(
                    t(
                      as.data.frame(
                        stats::embed(aa, 2)[,
                          c(2:1)]
                    )
                  ),
                    1, paste0, collapse = ""
                )
              )
            } else {
                apply(
                  stats::embed(aa, 2)[,
                    2:1], 1, paste0, collapse = ""
              )
            }
        }
    )
    gp <- lapply(
        dp, function(dp) {
            (10/(length(dp) +
                1)) * sum(guruprasad[dp], na.rm = TRUE)
        }
    )
    return(unlist(gp))
}

#' 2_Feature_determination
#'
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#' @importFrom stringr str_split
#' @noRd
group_freqs <- function(peptide, region) {
    # define groups
    group_list <- list()
    group_list[[1]] <- c("A", "C", "G", "S", "T")
    group_list[[2]] <- c("A", "B", "C", "D", "G", "N", "P", "S", "T", "V")
    group_list[[3]] <- c("A", "I", "L", "V")
    group_list[[4]] <- c("F", "H", "W", "Y")
    group_list[[5]] <- c("A", "C", "F", "G", "I", "L", "M", "P", "V", "W", "Y")
    group_list[[6]] <- c("D", "E", "H", "K", "N", "Q", "R", "S", "T", "Z")
    group_list[[7]] <- c("B", "D", "E", "H", "K", "R", "Z")
    group_list[[8]] <- c("H", "K", "R")
    group_list[[9]] <- c("B", "D", "E", "Z")
    aa_groups <- c(
        "Tiny", "Small", "Aliphatic", "Aromatic", "NonPolar", "Polar", "Charged",
        "Basic", "Acidic"
    )
    names(group_list) <- aa_groups
    # Prepare peptide
    split_peptide <- unlist(str_split(peptide, ""))
    number <- c()
    perc <- c()
    for (group in aa_groups) {
        number <- c(
            number, sum(
                table(split_peptide)[group_list[[group]]],
                na.rm = TRUE
            )
        )
        if (length(split_peptide) ==
            0) {
            perc <- c(perc, 0)
        } else {
            perc <- c(
                perc, (sum(
                  table(split_peptide)[group_list[[group]]],
                  na.rm = TRUE
              )/length(split_peptide)) *
                  100
            )
        }

    }
    names(number) <- paste0(region, "_", aa_groups, "_count")
    names(perc) <- paste0(region, "_", aa_groups, "_norm")
    stats <- c(number, round(perc, 3))
    return(stats)
}

#' 2_Feature_determination
#'
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#' @noRd
alakazam_properties <- function(peptide, region) {
    bulkiness <- alakazam::bulk(peptide)
    hydrophobicity <- alakazam::gravy(peptide)
    aliphatic_index <- alakazam::aliphatic(peptide)
    avg_polarity <- alakazam::polar(peptide)
    charge <- alakazam::charge(peptide)
    normalized_charge <- alakazam::charge(peptide, normalize = TRUE)
    alkzm_prprts <- c(
        bulkiness, hydrophobicity, aliphatic_index, avg_polarity, charge, normalized_charge
    )
    names(alkzm_prprts) <- c(
        "bulkiness", "hydrophobicity", "aliphatic_index", "avg_polarity", "charge",
        "normalized_charge"
    )
    names(alkzm_prprts) <- gsub(
        "^", paste0(region, "_", "alkzm_"),
        names(alkzm_prprts)
    )
    return(alkzm_prprts)
}

#' 2_Feature_determination
#'
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#' @noRd
Peptides_properties <- function(peptide, region) {
    aliphatic_index <- Peptides::aIndex(peptide)
    boman_index <- Peptides::boman(peptide)  #A protein has high binding potential if the index value is higher than 2.48
    charge <- Peptides::charge(peptide, pH = 7.4, pKscale = "EMBOSS")  #ph blood 7.4
    hmoment <- Peptides::hmoment(peptide)
    hydrophobicity <- Peptides::hydrophobicity(peptide)
    instability_index <- instaIndex_improved(peptide)
    mol_weight <- Peptides::mw(peptide)
    isoelectric_point <- Peptides::pI(peptide)

    # these need names
    names(aliphatic_index) <- "aliphatic_index"
    names(boman_index) <- "boman_index"
    names(charge) <- "charge"
    names(hmoment) <- "hmoment"
    names(hydrophobicity) <- "hydrophobicity"
    names(instability_index) <- "instability_index"
    names(mol_weight) <- "mol_weight"
    names(isoelectric_point) <- "isoelectric_point"

    # these are named already
    cruciani_properties <- unlist(Peptides::crucianiProperties(peptide))
    fasgai_vecs <- unlist(Peptides::fasgaiVectors(peptide))
    blosum_indicces <- unlist(Peptides::blosumIndices(peptide))
    kidera_factors <- unlist(Peptides::kideraFactors(peptide))
    mswhim_scores <- unlist(Peptides::mswhimScores(peptide))
    st_scales <- unlist(Peptides::stScales(peptide))
    t_scales <- unlist(Peptides::tScales(peptide))
    vhse_scales <- unlist(Peptides::vhseScales(peptide))
    z_scales <- unlist(Peptides::zScales(peptide))

    pptds_prprts <- c(
        aliphatic_index, boman_index, charge, hmoment, hydrophobicity, instability_index,
        mol_weight, isoelectric_point, cruciani_properties, fasgai_vecs, blosum_indicces,
        kidera_factors, mswhim_scores, st_scales, t_scales, vhse_scales, z_scales
    )
    names(pptds_prprts) <- gsub(
        "^", paste0(region, "_", "Peptides_"),
        names(pptds_prprts)
    )
    return(pptds_prprts)
}

#' 2_Feature_determination
#'
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#'
#' @noRd
edited_hot_cold_mtfs <- function() {
    hot_1_F <- "WRC"
    hot_1_R <- "GYW"
    hot_2_F <- "WA"
    hot_2_R <- "TW"
    cold_F <- "SYC"
    cold_R <- "GRS"
    potential_1_F <- "CRCY"
    potential_1_R <- "RGYG"
    potential_2_F <- "ATCT"
    potential_2_R <- "AGAT"
    hot_cold_motif_vec <- c(
        hot_1_F, hot_1_R, hot_2_F, hot_2_R, cold_F, cold_R, potential_1_F, potential_1_R,
        potential_2_F, potential_2_R
    )
    edit_mtf <- function(mtf) {
        mtf <- gsub("S", "[GC]", mtf)
        mtf <- gsub("W", "[AT]", mtf)
        mtf <- gsub("R", "[GA]", mtf)
        mtf <- gsub("Y", "[CT]", mtf)
        return(mtf)
    }
    edited_hot_cold_mtfs <- sapply(hot_cold_motif_vec, edit_mtf)
    names(edited_hot_cold_mtfs) <- paste0(
        c(
            rep("hot_", 4),
            rep("cold_", 2),
            rep("potential_hot_", 4)
        ),
        hot_cold_motif_vec
    )
    names(edited_hot_cold_mtfs) <- paste0(
        names(edited_hot_cold_mtfs),
        c(
            rep("_1", 2),
            rep("_2", 2),
            rep("", 2),
            rep("_1", 2),
            rep("_2", 2)
        )
    )
    names(edited_hot_cold_mtfs) <- paste0(
        names(edited_hot_cold_mtfs),
        c("_F", "_R")
    )
    return(edited_hot_cold_mtfs)
}

#' 2_Feature_determination
#'
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#'
#' @noRd
mut_points_hot_cold_mtfs <- function(edited_hot_cold_mtfs) {
    hot_1_F <- 3
    hot_1_R <- 1
    hot_2_F <- 2
    hot_2_R <- 1
    cold_F <- 3
    cold_R <- 1
    potential_1_F <- 3
    potential_1_R <- 2
    potential_2_F <- 3
    potential_2_R <- 2
    mut_points <- c(
        hot_1_F, hot_1_R, hot_2_F, hot_2_R, cold_F, cold_R, potential_1_F, potential_1_R,
        potential_2_F, potential_2_R
    )
    names(mut_points) <- names(edited_hot_cold_mtfs)
    return(mut_points)
}

#' 2_Feature_determination
#'
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#' @importFrom stringr str_locate_all
#' @noRd
locate_motif <- function(full_sequence, lengths, motif, motif_ID, motif_ini, regions) {
    loci <- unlist(data.frame(str_locate_all(full_sequence, motif))["start"])
    loci_mut_point <- loci + (motif_ini - 1)
    index_locis <- rep(
        "", length(lengths) -
            1
    )
    motif_count <- rep(
        0, length(lengths) -
            1
    )
    if (length(loci) >
        0) {
        for (i in c(
            1:(length(lengths) -
                1)
        )) {
            if (i != 1) {
                min <- sum(lengths[1:(i - 1)])
            } else {
                min <- 0
            }
            index_locis[i] <- paste(
                loci[which(loci_mut_point > min & loci_mut_point <= min + lengths[i])] -
                  min, collapse = "/"
            )
            motif_count[i] <- length(which(loci_mut_point > min & loci_mut_point <= min + lengths[i]))
        }
    }

    index_locis <- c(index_locis, paste(loci, collapse = "/"))
    names(index_locis) <- paste0(regions, "_", motif_ID, "_Pos")

    motif_count <- c(motif_count, sum(motif_count))
    names(motif_count) <- paste0(regions, "_", motif_ID, "_count")
    perc_locis <- 100 * motif_count/lengths
    names(perc_locis) <- paste0(regions, "_", motif_ID, "_norm")

    perc_locis[which(is.nan(perc_locis))] <- 0

    results <- list(index_locis = index_locis, motif_count = motif_count, perc_locis = perc_locis)
    return(results)
}

#' 2_Feature_determination
#'
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#' @import stringr
#' @import utils
#' @noRd
compare_sequences <- function(
    rep_seq, germ_seq, mode, nt_rep_lengths = NULL, aa_rep_lengths = NULL, nt_germ_lengths = NULL,
    aa_germ_lengths = NULL, NT_regions = NULL, AA_regions = NULL, ORF_rep = NULL,
    ORF_germ = NULL, aant_rep_seq=NULL, aant_germ_seq=NULL
) {

  ##lets find the indel coordinates and their respective ones
  if(mode=="AA") {
    germ_indel_seq=aant_germ_seq
    rep_indel_seq=aant_rep_seq
  } else {
    if(mode=="Codon") {
      germ_seq=substr(germ_seq, ORF_germ, nchar(germ_seq))
      rep_seq=substr(rep_seq, ORF_rep, nchar(rep_seq))
    }
    germ_indel_seq=germ_seq
    rep_indel_seq=rep_seq
  }

  indel=F
  indel_g=F
  indel_r=F
  indel_table=data.frame()

  indels_germ=as.data.frame(str_locate_all(germ_indel_seq, "[a-z]+")[[1]])
  if(nrow(indels_germ)>0){
    indels_germ=cbind(indels_germ, rep("Germ", nrow(indels_germ)))
    colnames(indels_germ)[3]="Seq"
    indel_g=T
  }

  indels_rep=as.data.frame(str_locate_all(rep_indel_seq, "[a-z]+")[[1]])
  if(nrow(indels_rep)>0){
    indels_rep=cbind(indels_rep, rep("Rep", nrow(indels_rep)))
    colnames(indels_rep)[3]="Seq"
    indel_r=T
  }
  if(indel_g && indel_r) {
    indel_table=rbind(indels_rep, indels_germ)
    indel=T
  } else if(indel_g){
    indel_table=indels_germ
    indel=T
  } else if (indel_r) {
    indel_table=indels_rep
    indel=T
  }

  if(indel){
    indel_table=indel_table[order(indel_table[,"start"]),]
    indel_table$length=indel_table$end-indel_table$start+1

    if(indel_g) {
      for (row in  which(indel_table[,"Seq"]=="Germ")){
        index_lower_start=intersect(which(indel_table$start < indel_table[row, "start"]),
                                    which(indel_table[,"Seq"]!="Germ"))
        indel_table[row, c("start", "end")]=indel_table[row, c("start", "end")]-sum(indel_table$length[index_lower_start])
      }

    }

    indel_positions=data.frame(Seq=as.character(), Position=as.numeric(),AA_Codon_position=as.numeric())

    for(row in c(1:nrow(indel_table))) {
      start=indel_table$start[row]
      for(nt in c(1:indel_table$length[row])) {
        indel_positions=rbind(indel_positions, data.frame(Seq=indel_table$Seq[row], Position=start,AA_Codon_position=ceiling(start/3)))
        start=start+1
      }
    }
    indel_positions$AA_Codon_full=sapply(c(1:nrow(indel_positions)), function(z) length(intersect(which(indel_positions$Seq ==indel_positions$Seq[z]),
                                                                                                  which(indel_positions$AA_Codon_position ==indel_positions$AA_Codon_position[z]))) )
  }

  if(mode=="Codon") {

    if(indel) {
      mut_germ_seq=strsplit(germ_seq, split="")[[1]]
      mut_rep_seq=strsplit(rep_seq, split="")[[1]]
      mut_germ_seq=mut_germ_seq[which(c(1:length(mut_germ_seq)) %!in% indel_positions[which(indel_positions[,"Seq"]=="Germ"),"Position"])  ]

      for(i in indel_positions[which(indel_positions[,"Seq"]!="Germ"),"Position"])  {
        mut_germ_seq <- append(mut_germ_seq, mut_rep_seq[i], i-1)
      }

      mut_germ_seq=toupper(paste(mut_germ_seq, collapse=""))
      mut_rep_seq=toupper(paste(mut_rep_seq, collapse=""))
    }
  }

  rep_seq <- toupper(rep_seq)
  germ_seq <- toupper(germ_seq)


  if (mode == "AA") {
    regions <- AA_regions
    aa_names <- c(
      "Phe", "Met", "Leu", "Ile", "Val", "Pro", "Tyr", "Trp", "Cys", "Ala",
      "Gly", "Ser", "Thr", "His", "Glu", "Gln", "Asp", "Asn", "Lys", "Arg"
    )
    elements <- seqinr::a(aa_names)

    group_list <- list()
    group_list[[1]] <- c("A", "C", "G", "S", "T")
    group_list[[2]] <- c("A", "B", "C", "D", "G", "N", "P", "S", "T", "V")
    group_list[[3]] <- c("A", "I", "L", "V")
    group_list[[4]] <- c("F", "H", "W", "Y")
    group_list[[5]] <- c("A", "C", "F", "G", "I", "L", "M", "P", "V", "W", "Y")
    group_list[[6]] <- c("D", "E", "H", "K", "N", "Q", "R", "S", "T", "Z")
    group_list[[7]] <- c("B", "D", "E", "H", "K", "R", "Z")
    group_list[[8]] <- c("H", "K", "R")
    group_list[[9]] <- c("B", "D", "E", "Z")
    aa_groups <- c(
      "Tiny", "Small", "Aliphatic", "Aromatic", "NonPolar", "Polar", "Charged",
      "Basic", "Acidic"
    )
    names(group_list) <- aa_groups

    rep_seq=strsplit(rep_seq, split="")[[1]]
    germ_seq=strsplit(germ_seq, split="")[[1]]

    rep_lengths <- aa_rep_lengths
    germ_lengths <- aa_germ_lengths

  } else if (mode == "NT") {
    regions <- NT_regions
    elements <- c("A", "C", "G", "T")
    transition_muts <- c("A_to_G", "G_to_A", "T_to_C", "C_to_T")
    rep_lengths <- nt_rep_lengths
    germ_lengths <- nt_germ_lengths

    rep_seq=strsplit(rep_seq, split="")[[1]]
    germ_seq=strsplit(germ_seq, split="")[[1]]

  } else if (mode == "Codon") {
    elements <- sort(names(Biostrings::GENETIC_CODE))
    regions <- AA_regions

    rep_seq=as.character(Biostrings::codons(Biostrings::DNAString(rep_seq)))
    germ_seq=as.character(Biostrings::codons(Biostrings::DNAString(germ_seq)))

    rep_codon_count <- c()
    germ_codon_count <- c()


    rep_lengths <- aa_rep_lengths
    germ_lengths <- aa_germ_lengths
  } else {
    # print("There is an error with the mode")
  }


  ###coordinates of each region
  coord_rep_regions=c(1:rep_lengths[8])
  coord_germ_regions=c(1:germ_lengths[8])

  tmp_rep_lengths=rep_lengths[1:7]
  tmp_germ_lengths=germ_lengths[1:7]
  names(tmp_rep_lengths)=gsub("_length","", names(tmp_rep_lengths))
  names(tmp_germ_lengths)=gsub("_length","", names(tmp_germ_lengths))
  for(i in c(2:7)) {
    tmp_rep_lengths[i]=tmp_rep_lengths[i]+tmp_rep_lengths[i-1]
    tmp_germ_lengths[i]=tmp_germ_lengths[i]+tmp_germ_lengths[i-1]
  }
  names(coord_rep_regions)=sapply(coord_rep_regions, function(z) names(tmp_rep_lengths)[ which(tmp_rep_lengths>=z)[1]])
  names(coord_germ_regions)=sapply(coord_germ_regions, function(z) names(tmp_germ_lengths)[ which(tmp_germ_lengths>=z)[1]])

  coord_regions=list(coord_rep_regions, coord_germ_regions)
  counts_per_sequence=list(Rep=c(), Germ=c())
  ###count codons
  if (mode == "Codon") {
    for (region in regions) {


      if(region == "AA_Whole") {
        location=c(regions)
      } else {
        location=region
      }

      for(i in c(1:2)) {
        coords_region=which(names(coord_regions[[i]]) %in% location)
        codons_present=table(rep_seq[coords_region])

        tmp_counts=rep(0, length(elements[which(!(elements %in% names(codons_present)))]))
        names(tmp_counts)= elements[which(!(elements %in% names(codons_present)))]
        tmp_counts=c(tmp_counts, codons_present)
        tmp_counts=tmp_counts[sort(names(tmp_counts))]
        names(tmp_counts)=paste(region, names(tmp_counts), sep="_")

        if(length(coords_region)!=0) {
          tmp_counts_norm=tmp_counts/length(coords_region)
        } else {
          tmp_counts_norm=tmp_counts
        }
        names(tmp_counts)=paste(names(tmp_counts), "counts", sep="_")
        names(tmp_counts_norm)=paste(names(tmp_counts), "norm", sep="_")


        counts_per_sequence[[i]]=c(counts_per_sequence[[i]], tmp_counts, tmp_counts_norm)
      }

    }

  }

  ##confirmed with only tmp_counts        sum(counts_per_sequence[[i]][which(startsWith(names(counts_per_sequence[[i]]), "AA_Whole"))])


  ## NT to NT
  if (mode == "NT" || mode =="AA" || mode =="Codon") {

    counts=c()
    for (region in regions) {


      if(region == "AA_Whole" || region == "NT_Whole") {
        location=c(regions)
      } else {
        location=region
      }
      coords_region=which(names(coord_regions[[1]]) %in% location)

      coords_rep=which(names(coord_regions[[1]]) %in% location)
      coords_germ=which(names(coord_regions[[2]]) %in% location)


      if(indel) {
        if (mode == "NT") {
          coords_rep=setdiff(coords_rep,indel_positions$Position[which(indel_positions$Seq=="Rep")] )
          coords_germ=setdiff(coords_germ,indel_positions$Position[which(indel_positions$Seq=="Germ")] )
        } else if (mode == "AA"){
          coords_rep=setdiff(coords_rep,
                             unique(indel_positions$AA_Codon_position[intersect(which(indel_positions$Seq=="Rep"),
                                                                                which(indel_positions$AA_Codon_full>1))]))
          coords_germ=setdiff(coords_germ,
                              unique(indel_positions$AA_Codon_position[intersect(which(indel_positions$Seq=="Germ"),
                                                                                 which(indel_positions$AA_Codon_full>1))]))
        } else if (mode=="Codon") {
          coords_rep=setdiff(coords_rep,
                             unique(indel_positions$AA_Codon_position[intersect(which(indel_positions$Seq=="Rep"),
                                                                                which(indel_positions$AA_Codon_full>1))]))
          coords_germ=setdiff(coords_germ,
                              unique(indel_positions$AA_Codon_position[intersect(which(indel_positions$Seq=="Germ"),
                                                                                 which(indel_positions$AA_Codon_full>1))]))
        }
      }




      # if(length(coords_rep) != length(coords_germ)) {
      #   print("Unusual")
      # }

      differences=which(rep_seq[coords_rep] != germ_seq[coords_germ])

      differences=differences[1:(length(differences)-length(which(c(is.na(rep_seq[coords_rep][differences]), is.na(germ_seq[coords_germ][differences])))))]
      tmp_counts=compare_diffs(elements_rep=rep_seq[coords_rep][differences],
                               elements_germ=germ_seq[coords_germ][differences],
                               elements=elements,
                               region=region)



      ### lets add transversions and transitions
      if (mode == "NT") {
        if (sum(tmp_counts) != 0) {
          transition_counts <- sum(
            tmp_counts[which(
              grepl(paste(transition_muts, collapse="|") , names(tmp_counts))
            )]
          )
          transversion_counts <- sum(tmp_counts) - transition_counts
          if (transversion_counts == 0) {
            Ratio <- 10
          } else {
            Ratio <- transition_counts/transversion_counts
          }

        } else {
          transition_counts <- 0
          transversion_counts <- 0
          Ratio <- 0
        }
        tr_tmp <- c(transition_counts, transversion_counts, Ratio)
        names(tr_tmp) <- paste(
          region, c("Transitions", "Transversions", "Ratio_Transitions-Transversions"),
          sep = "_"
        )

        if(length(coords_rep)!=0 && length(coords_germ)!=0) {
          tmp_counts_norm=tmp_counts/length(coords_rep)
        } else {
          tmp_counts_norm=tmp_counts
        }

        names(tmp_counts) =paste(names(tmp_counts), "count", sep="_")
        names(tmp_counts_norm) =paste(names(tmp_counts_norm), "norm", sep="_")


        counts <- c(counts, tmp_counts, tmp_counts_norm, tr_tmp)

      } else if(mode =="AA") {
        tmp_counts_2=tmp_counts
        combinations <- expand.grid(names(group_list), names(group_list))
        combinations <- combinations[which(combinations$Var1 != combinations$Var2),
        ]

        combinations$Counts=0
        for (i in which(tmp_counts_2 !=0)) {
          from=strsplit(names(tmp_counts_2)[i], split="_")[[1]][5]
          to=strsplit(names(tmp_counts_2)[i], split="_")[[1]][3]

          from_type=names(group_list)[grep(from, group_list)]
          to_type=names(group_list)[grep(to, group_list)]

          indexes_combinations=intersect(which(combinations$Var1 %in% from_type),
                                         which(combinations$Var2 %in% to_type))
          combinations[indexes_combinations, "Counts"]=combinations[indexes_combinations, "Counts"]+tmp_counts_2[i]



        }
        tmp_type_counts=combinations$Counts
        names(tmp_type_counts)=paste(region, paste(combinations$Var1, combinations$Var2, sep="_to_"), sep="_")

        tmp_type_counts=tmp_type_counts[sort(names(tmp_type_counts))]

        if(length(coords_region)!=0) {
          tmp_counts_norm=tmp_counts/length(coords_region)
          tmp_type_counts_norm=tmp_type_counts/length(coords_rep)
        } else {
          tmp_counts_norm=tmp_counts
          tmp_type_counts_norm=tmp_type_counts
        }
        names(tmp_counts)=paste(names(tmp_counts), "count", sep="_")
        names(tmp_counts_norm)=paste(names(tmp_counts_norm), "norm", sep="_")
        names(tmp_type_counts)=paste(names(tmp_type_counts), "count", sep="_")
        names(tmp_type_counts_norm)=paste(names(tmp_type_counts_norm), "norm", sep="_")
        # counts <- c(counts, tmp_counts, tmp_counts_norm, tmp_type_counts, tmp_type_counts_norm)

        counts <- c(counts, tmp_counts, tmp_type_counts)
      }  else if(mode =="Codon") {

        coords_rep=which(names(coord_regions[[1]]) %in% location)
        coords_germ=which(names(coord_regions[[2]]) %in% location)

        if(indel){
          # coords_rep=setdiff(coords_rep,
          #                    unique(indel_positions$AA_Codon_position[intersect(which(indel_positions$Seq=="Rep"),
          #                                                                       which(indel_positions$AA_Codon_full>1))]))
          # coords_germ=setdiff(coords_germ,
          #                     unique(indel_positions$AA_Codon_position[intersect(which(indel_positions$Seq=="Germ"),
          #                                                                        which(indel_positions$AA_Codon_full>1))]))
          #
          # index_1_rep=c(indel_positions$AA_Codon_position[intersect(which(indel_positions$Seq=="Rep"),
          #                                                           which(indel_positions$AA_Codon_full==1))],
          #               coords_rep[which(coords_germ==unique(indel_positions$AA_Codon_position[intersect(which(indel_positions$Seq=="Germ"),
          #                                                                                     which(indel_positions$AA_Codon_full==1))]))]
          # )
          # index_1_germ=c(indel_positions$AA_Codon_position[intersect(which(indel_positions$Seq=="Germ"),
          #                                                           which(indel_positions$AA_Codon_full==1))],
          #                coords_germ[which(coords_rep==unique(indel_positions$AA_Codon_position[intersect(which(indel_positions$Seq=="Rep"),
          #                                                                                     which(indel_positions$AA_Codon_full==1))]))]
          # )
          #
          # coords_rep=setdiff(coords_rep,
          #                    index_1_rep)
          # coords_germ=setdiff(coords_germ,
          #                     index_1_germ)

          rep_seq=as.character(Biostrings::codons(Biostrings::DNAString(mut_rep_seq)))
          germ_seq=as.character(Biostrings::codons(Biostrings::DNAString(mut_germ_seq)))
          coords_germ=coords_rep
        }



        # if(length(coords_rep) != length(coords_germ)) {
        #   print("Unusual")
        # }

        differences=which(rep_seq[coords_rep] != germ_seq[coords_germ])

        if(length(differences)==0) {
          silent=0
          replacement=0
          Ratio=0
          total=0
        } else {
          differences=differences[1:(length(differences)-length(which(c(is.na(rep_seq[coords_rep][differences]), is.na(germ_seq[coords_germ][differences])))))]

          type_mut=sapply(differences, function(z) if(Biostrings::GENETIC_CODE[match(rep_seq[coords_rep][z],
                                                                         names(Biostrings::GENETIC_CODE))] ==
                                                      Biostrings::GENETIC_CODE[match(germ_seq[coords_germ][z],
                                                                         names(Biostrings::GENETIC_CODE))]) {
            "Silent"
          } else {
            "Replacement"
          })
          num_muts=  sapply(differences,
                            function(z)  adist(rep_seq[coords_rep][z],germ_seq[coords_germ][z])[1,1] )

          names(num_muts)=paste(rep_seq[coords_rep][differences],germ_seq[coords_germ][differences], sep="_to_")

          if(length(which(type_mut=="Silent"))>0) {
            silent=sum(sapply(which(type_mut=="Silent"),
                              function(z)
                                sum(num_muts[z]) ))
          } else {
            silent=0
          }

          if(length(which(type_mut!="Silent"))>0) {
            replacement=sum(sapply(which(type_mut!="Silent"),
                                   function(z)
                                     sum(num_muts[z]) ))
          } else{
            replacement=0
          }


          total=replacement+silent
          if (silent == 0) {
            Ratio <- 10
          } else {
            Ratio <- replacement/silent
          }
        }


        count_tmp <- c(silent, replacement, Ratio, total)
        names(count_tmp) <- paste(
          region, c("Silent_muts", "Replacement_muts",
                    "Ratio_Replacement-Silent", "Total_muts"),
          sep = "_"
        )

        if(length(coords_rep)!=0) {
          count_tmp_norm=count_tmp/(length(coords_rep)*3)
        } else {
          count_tmp_norm=count_tmp
        }

        names(count_tmp)=paste(names(count_tmp), "counts", sep="_")
        names(count_tmp_norm)=paste(names(count_tmp_norm), "norm", sep="_")

        names(tmp_counts)=paste(names(tmp_counts), "count", sep="_")

        if(region=="AA_Whole") {
          counts <- c(counts, tmp_counts, count_tmp, count_tmp_norm)
        } else {
          counts <- c(counts, count_tmp, count_tmp_norm)
        }

      }

    }

  }



  ##to_check
  if (mode == "Codon") {
    counts_full <- list(
      counts = counts, rep_codon_count = counts_per_sequence$Rep, germ_codon_count = counts_per_sequence$Germ
    )
  } else {
    counts_full <- counts
  }
  return(counts_full)
}

#' 2_Feature_determination
#'
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#' @noRd
compare_diffs<- function(elements_rep, elements_germ, elements, region){

  m_vec <- paste0(elements_rep, elements_germ)
  m_vec <- m_vec[which(!grepl("NA",m_vec))]
  combinations <- expand.grid(elements, elements)
  combinations <- combinations[which(combinations$Var1 != combinations$Var2),
  ]

  counts <- rep(0, nrow(combinations))
  names(counts) <- paste(combinations[, 1], combinations[, 2], sep = "_to_")
  combinations <- paste0(combinations[, 1], combinations[, 2])

  table_counts <- table(m_vec)
  indexes <- match(combinations, names(table_counts))
  counts[which(!is.na(indexes))] <- table_counts[indexes[which(!is.na(indexes))]]
  names(counts) <- paste(
    region, names(counts),
    sep = "_"
  )

  counts=counts[sort(names(counts))]
  return(counts)

}
