#' 2_Feature_determination
#'
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#' @import Biostrings
#' @noRd
translate_fun <- function(x) {

    if (x == "") {
        return(list("", ""))
    } else {
        ORFs <- lapply(
            1:3, function(pos) Biostrings::subseq(
                Biostrings::DNAString(as.character(x)),
                start = pos
            )
        )
        suppressWarnings(
            translated_ORFs <- lapply(
                ORFs, Biostrings::translate, genetic.code = Biostrings::GENETIC_CODE,
                no.init.codon = TRUE, if.fuzzy.codon = "error"
            )
        )

        # Search for the best ORF
        Stop_per_ORF <- vector()
        for (i in 1:3) {
            Stop_per_ORF[i] <- sum(length(Biostrings::matchPattern("*", unlist(translated_ORFs[[i]]))))
        }

        # Find STOP codons per ORF and choose the ORF with the less
        BestORF <- which.min(Stop_per_ORF)
        aa_seq <- as.character(translated_ORFs[[BestORF]])

        if (Stop_per_ORF[BestORF] == 0) {
            return(list(aa_seq, BestORF))
        } else {
            return(list("NO_GOOD_ORF", 0))  #If every ORF has stop codons. This value is used for filtering afterwards.
        }
    }

}

#' 2_Feature_determination
#'
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#' @import stringr
#' @import Biostrings
#' @noRd
subseq_func <- function(NT_seq, size_NT_regions, FWR1partial = F, FWR4partial = F) {
    ########### TOCAR
    AA_info <- translate_fun(NT_seq)
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
        size_AA_regions <- round(size_NT_regions/3)
        size_AA_regions[length(size_AA_regions) -
            1] <- floor(
            (size_NT_regions[length(size_AA_regions) -
                1]/3)
        )
        size_AA_regions[length(size_AA_regions)] <- floor((size_NT_regions[length(size_AA_regions)]/3))
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

        seq_list[length(seq_list) +
            1] <- AA_seq
        seq_list[length(seq_list) +
            1] <- ORF

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
    if (width(as.character(seq_repertoire)) ==
        0) {
        lv_dist_norm <- 0
    } else {
        lv_dist_norm <- (lv_dist/width(as.character(seq_repertoire))) *
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
    lv_info <- lv_changes(seq_repertoire, seq_rec.germline, region)
    lv_dist <- lv_info[1]
    osa_dist <- stringdist::stringdist(seq_rec.germline, seq_repertoire, method = "osa")
    no.trans <- lv_dist - osa_dist
    tab <- drop(
        attr(
            adist(seq_rec.germline, seq_repertoire, count = TRUE),
            "counts"
        )
    )
    no.ins <- tab[["ins"]] - no.trans
    no.del <- tab[["del"]] - no.trans
    no.subs <- tab[["sub"]]
    changes <- c(no.subs, no.ins, no.del, no.trans)
    if (sum(changes) !=
        0) {
        changes_plus <- c(
            changes, sum(changes),
            (changes[1]/sum(changes)) *
                100, (changes[2]/sum(changes)) *
                100, (changes[3]/sum(changes)) *
                100, (changes[4]/sum(changes)) *
                100
        )
    } else {
        changes_plus <- c(
            changes, sum(changes),
            0, 0, 0, 0
        )
    }
    names(changes_plus) <- c(
        "Sub", "Ins", "Del", "Trasl", "SIDT_sum", "Sub_prc", "Ins_prc", "Del_prc",
        "Trasl_prc"
    )
    names(changes_plus) <- gsub(
        "^", paste0(
            as.character(region),
            "_"
        ),
        names(changes_plus)
    )
    return(c(changes_plus, lv_info))
}

#' 2_Feature_determination
#'
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#' @import Peptides
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
#' @description A utils function
#'
#' @return The return value, if any, from executing the utility.
#' @importFrom utils adist
#' @noRd
aa_changes <- function(seq_repertoire, seq_rec.germline, region) {
    tab <- drop(
        attr(
            adist(seq_rec.germline, seq_repertoire, count = TRUE),
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
#' @import seqinr
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
        print("[ERROR] Wrong region included")
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
#' @import stats
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
    aa <- Peptides:::aaCheck(seq)
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
                        embed(aa, 2)[,
                          c(2:1)]
                    )
                  ),
                    1, paste0, collapse = ""
                )
              )
            } else {
                apply(
                  embed(aa, 2)[,
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
#' @import alakazam
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
#' @import Peptides
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
#' @import seqinr
#' @import Biostrings
#' @import stringr
#' @import utils
#' @noRd
compare_sequences <- function(
    rep_seq, germ_seq, mode, nt_rep_lengths = NULL, aa_rep_lengths = NULL, nt_germ_lengths = NULL,
    aa_germ_lengths = NULL, NT_regions = NULL, AA_regions = NULL, ORF_rep = NULL,
    ORF_germ = NULL
) {

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

        rep_lengths <- aa_rep_lengths
        germ_lengths <- aa_germ_lengths

    } else if (mode == "NT") {
        regions <- NT_regions
        elements <- c("A", "C", "G", "T")
        transition_muts <- c("A_to_G", "G_to_A", "T_to_C", "C_to_T")
        rep_lengths <- nt_rep_lengths
        germ_lengths <- nt_germ_lengths
    } else if (mode == "Codon") {
        elements <- sort(names(GENETIC_CODE))
        regions <- AA_regions
        ORF_rep <- as.numeric(ORF_rep)
        ORF_germ <- as.numeric(ORF_germ)

        codons <- elements
        codon_dict <- as.character(c(0:9, letters, LETTERS, "+", "-"))
        names(codon_dict) <- codons
        rep_seq <- substr(
            rep_seq, ORF_rep, ORF_rep - 1 + 3 * aa_rep_lengths[which(
                startsWith(
                  names(aa_rep_lengths),
                  "AA_Whole"
              )
            )]
        )
        rep_seq_codon <- as.character(codons(DNAString(rep_seq)))
        rep_seq <- paste(
            codon_dict[match(rep_seq_codon, names(codon_dict))],
            collapse = ""
        )
        germ_seq <- substr(
            germ_seq, ORF_germ, ORF_germ - 1 + 3 * aa_germ_lengths[which(
                startsWith(
                  names(aa_germ_lengths),
                  "AA_Whole"
              )
            )]
        )
        germ_seq_codon <- as.character(codons(DNAString(germ_seq)))
        germ_seq <- paste(
            codon_dict[match(germ_seq_codon, names(codon_dict))],
            collapse = ""
        )

        rep_codon_count <- c()
        germ_codon_count <- c()



        rep_lengths <- aa_rep_lengths
        germ_lengths <- aa_germ_lengths
    } else {
        print("There is an error with the mode")
    }

    counts_full <- c()
    for (region in regions) {
        germ_start <- 1
        rep_start <- 1
        region_i <- which(
            startsWith(
                names(rep_lengths),
                region
            )
        )
        if (endsWith(region, "Whole")) {
            rep_end <- nchar(rep_seq)
            germ_end <- nchar(germ_seq)
        } else {
            if (region_i != 1) {
                rep_start <- sum(rep_lengths[1:(region_i - 1)]) +
                  1
                germ_start <- sum(germ_lengths[1:(region_i - 1)]) +
                  1
            }
            rep_end <- rep_lengths[region_i] + rep_start - 1
            germ_end <- germ_lengths[region_i] + germ_start - 1

        }
        rep_seq_tmp <- substr(rep_seq, rep_start, rep_end)
        if (nchar(rep_seq_tmp) !=
            rep_lengths[region_i]) {
            print(
                paste(
                  "ERROR: calculating rep_seq length", "Mode:", mode, "Rep:", rep_seq,
                  "Germ:", germ_seq, "Region: ", region, sep = " "
              )
            )
        }
        germ_seq_tmp <- substr(germ_seq, germ_start, germ_end)
        if (nchar(germ_seq_tmp) !=
            germ_lengths[region_i]) {
            print(
                paste(
                  "ERROR: calculating germ_seq length:", "Mode:", mode, "Rep:", rep_seq,
                  "Germ:", germ_seq, "Region: ", region, sep = " "
              )
            )

        }


        if (mode == "Codon") {

            tmp_rep_codon_count <- rep(0, length(elements))
            names(tmp_rep_codon_count) <- elements
            tmp_codon_counts <- table(rep_seq_codon[rep_start:rep_end])
            tmp_rep_codon_count[match(
                names(tmp_codon_counts),
                names(tmp_rep_codon_count)
            )] <- tmp_codon_counts
            names(tmp_rep_codon_count) <- paste0(
                paste(
                  region, names(tmp_rep_codon_count),
                  sep = "_"
              ),
                "_counts"
            )
            rep_codon_count <- c(rep_codon_count, tmp_rep_codon_count)

            tmp_germ_codon_count <- rep(0, length(elements))
            names(tmp_germ_codon_count) <- elements
            tmp_codon_counts <- table(germ_seq_codon[germ_start:germ_end])
            tmp_germ_codon_count[match(
                names(tmp_codon_counts),
                names(tmp_germ_codon_count)
            )] <- tmp_codon_counts
            names(tmp_germ_codon_count) <- paste0(
                paste(
                  region, names(tmp_germ_codon_count),
                  sep = "_"
              ),
                "_counts"
            )
            germ_codon_count <- c(germ_codon_count, tmp_germ_codon_count)
        }

        comp_unr_to_rec1 <- drop(
            attr(
                adist(rep_seq_tmp, germ_seq_tmp, counts = T),
                "trafos"
            )
        )
        comp_rec_to_unr1 <- drop(
            attr(
                adist(germ_seq_tmp, rep_seq_tmp, counts = T),
                "trafos"
            )
        )

        comp_unr_to_rec <- gsub("[^MSD]", "", comp_unr_to_rec1)
        comp_rec_to_unr <- gsub("[^MSD]", "", comp_rec_to_unr1)
        #
        if (!(stringr::str_length(comp_unr_to_rec) ==
            stringr::str_length(rep_seq_tmp) &&
            stringr::str_length(comp_rec_to_unr) ==
                stringr::str_length(germ_seq_tmp))) {
            print(
                paste(
                  "[ERROR] Comparing mutations went wrong:", "Mode:", mode, "Rep:",
                  rep_seq, "Germ:", germ_seq, "Values:", print(stringr::str_length(comp_unr_to_rec)),
                  stringr::str_length(rep_seq_tmp),
                  stringr::str_length(comp_rec_to_unr),
                  stringr::str_length(germ_seq_tmp),
                  comp_unr_to_rec, "vs", rep_seq_tmp, "AND", comp_rec_to_unr, "vs",
                  germ_seq_tmp, "AND", comp_unr_to_rec1, "vs", comp_rec_to_unr1,
                  sep = " "
              )
            )
        }

        subs_in_unrec <- str_locate_all(comp_unr_to_rec, "S")[[1]][,
            "start"]
        subs_in_rec <- str_locate_all(comp_rec_to_unr, "S")[[1]][,
            "start"]
        #
        if (mode == "Codon") {
            m_from <- rep_seq_codon[subs_in_rec]
            m_to <- germ_seq_codon[subs_in_unrec]
        } else {
            m_from <- str_split(germ_seq, "")[[1]][subs_in_rec]
            m_to <- str_split(rep_seq, "")[[1]][subs_in_unrec]
        }

        m_vec <- paste0(m_from, m_to)

        combinations <- expand.grid(elements, elements)
        combinations <- combinations[which(combinations$Var1 != combinations$Var2),
            ]

        counts <- rep(0, nrow(combinations))
        names(counts) <- paste(combinations[, 1], combinations[, 2], sep = "_to_")
        combinations <- paste0(combinations[, 1], combinations[, 2])
        table_counts <- table(m_vec)
        ############### reconvertimos variables a codones

        if (mode == "Codon") {
            names(table_counts) <- sapply(
                names(table_counts),
                function(x) paste(
                  names(codon_dict)[match(
                    strsplit(x, split = "")[[1]],
                    (codon_dict)
                )],
                  collapse = ""
              )
            )

        }
        indexes <- match(combinations, names(table_counts))
        counts[which(!is.na(indexes))] <- table_counts[indexes[which(!is.na(indexes))]]
        names(counts) <- paste(
            region, names(counts),
            sep = "_"
        )


        if (mode == "NT") {
            total_counts <- sum(counts)
            if (total_counts != 0) {
                transition_counts <- sum(
                  counts[which(
                    endsWith(
                      names(counts),
                      transition_muts
                  )
                )]
              )
                transversion_counts <- total_counts - transition_counts
                if (transversion_counts == 0) {
                  Ratio <- 1000
                } else {
                  Ratio <- transition_counts/transversion_counts
                }

            } else {
                transition_counts <- 0
                transversion_counts <- 0
                Ratio <- 0
            }


            count_tmp <- c(transition_counts, transversion_counts, Ratio)
            names(count_tmp) <- paste(
                region, c("Transitions", "Transversions", "Ratio_Transitions-Transversions"),
                sep = "_"
            )
            counts <- c(counts, count_tmp)
        } else if (mode == "AA") {
            ### Changes between aa groups: This is a bit more complex and
            ### maaaaaaaaaaaaaaaaybe not super informative, so lets leave it
            ### like this for now
        } else if (mode == "Codon") {
            ### R/S mutations go in the next if because for now we dont want to
            ### include the codon mutations in counts_full, its too much
        }

        counts_norm <- counts/(1 + rep_end - rep_start)
        counts_norm[which(is.nan(counts_norm))] <- 0
        names(counts_norm) <- paste0(
            names(counts_norm),
            "_norm"
        )
        names(counts) <- paste0(
            names(counts),
            "_count"
        )
        counts <- c(counts, counts_norm)
        if (mode != "Codon") {
            counts_full <- c(counts_full, counts)
        } else {
            total_muts <- sum(table_counts)
            if (total_muts != 0) {
                silent_muts <- sum(
                  sapply(
                    names(table_counts),
                    function(x) if (GENETIC_CODE[match(
                      substr(x, 1, 3),
                      names(GENETIC_CODE)
                  )] ==
                      GENETIC_CODE[match(
                        substr(x, 4, 6),
                        names(GENETIC_CODE)
                    )]) {
                      1
                    } else {
                      0
                    }
                )
              )
                replacement_muts <- total_muts - silent_muts
                if (replacement_muts == 0) {
                  Ratio <- 1000
                } else {
                  Ratio <- silent_muts/replacement_muts
                }

            } else {
                silent_muts <- 0
                replacement_muts <- 0
                Ratio <- 0
            }

            count_tmp <- c(silent_muts, replacement_muts, Ratio)
            names(count_tmp) <- paste(
                region, c("Silent_muts", "Replacement_muts", "Ratio_Silent-Replacement"),
                sep = "_"
            )
            ### We could normalize by length of the region... repertoire region
            ### I guess
            count_tmp_norm <- count_tmp/(1 + rep_end - rep_start)
            count_tmp_norm[which(is.nan(count_tmp_norm))] <- 0
            names(count_tmp_norm) <- paste0(
                names(count_tmp_norm),
                "_norm"
            )
            counts_full <- c(counts_full, count_tmp, count_tmp_norm)
        }

    }

    if (mode == "Codon") {
        counts_full <- list(
            counts = counts_full, rep_codon_count = rep_codon_count, germ_codon_count = germ_codon_count
        )
    } else {
        counts_full <- counts_full
    }
    return(counts_full)
}

