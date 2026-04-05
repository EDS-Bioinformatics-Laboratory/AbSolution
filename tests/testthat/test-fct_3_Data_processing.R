test_that("Filtering comparing groups works", {

  FBM=bigstatsr::as_FBM(matrix(c(rep(c(0:5),2),4,rep(14,4),15,rep(15,8)),
                              13)
                        )

  Short_DF=data.frame(
    Sequence_type=rep(c("Repertoire"), 13),
    Patient_Sample=c(rep(c("Pat_a","Pat_b"), each=5), "Pat_b", rep("Pat_c",2)),
    Chain=rep("IGH",13),
    V_and_J=c(rep("VJ",3),rep("VJ2",2),rep("VJ",5),"VJ2",rep("VJ",2)),
    Groupname=c(rep("A",4),"C",rep("B",6),"C","A"),
    Clonedef=c(rep("Clone_1a", 3), "Clone_2a", "Clone_3a",
               rep("Clone_1b", 3), "Clone_2b", "Clone_3b","Clone_4b",
               "Clone_1c","Clone_2c")
  )
  Short_DF=data.table::as.data.table(Short_DF)

  Header=c("AA_Whole_Replacement_muts_counts",
           "NT_CDR3_length")


  expect_equal(c(1:13),
               AbSolution:::filter_merged(
                 FBM=FBM,
                 merged_df=Short_DF,
                 merged_header=Header,
                 use_rgermline=FALSE,
                 use_repertoire=TRUE,
                 use_productive=TRUE,
                 use_nonproductive= FALSE,
                 my_regions= c("CDR3"),
                 my_var_elements=c("NT"),
                 my_vars=c(
                   "Length"
                 ),
                 my_vartypes="Baseline",
                 groups="Groupname",
                 group_A="A",
                 group_B="B",
                 group_C=c("C"),
                 univlog=T,
                 samples_to_keep=c("Pat_a",
                                   "Pat_b",
                                   "Pat_c"),
                 variables_to_remove=NULL,
                 pval_type="Bonferroni",
                 pval_cutoff=0.05,
                 estimate_cutoff=1,
                 number_selected_vars="All",
                 VJ_deselected= NULL,
                 R_mut_threshold_min=0,
                 R_mut_threshold_max=5,
                 V_J_to_use=c("VJ","VJ2"),
                 to_compare_groups=T,
                 use_sharedVDJ=F,
                 VDJ_normalized_per_size=F,
                 VDJ_maximize_clones=F,
                 VDJ_normalized_per_sample=F,
                 my_clone_def="Clonedef",
                 seed=1234,
                 chains="IGH",
                 igsubtypes=NULL
               )$ROWS)
})



test_that("Filtering used shared-VJ works", {

  FBM=bigstatsr::as_FBM(matrix(c(rep(c(0:5),2),4,rep(14,4),15,rep(15,8)),
                               13)
  )

  Short_DF=data.frame(
    Sequence_type=rep(c("Repertoire"), 13),
    Patient_Sample=c(rep(c("Pat_a","Pat_b"), each=5), "Pat_b", rep("Pat_c",2)),
    Chain=rep("IGH",13),
    V_and_J=c(rep("VJ",3),rep("VJ2",2),rep("VJ",5),"VJ2",rep("VJ",2)),
    Groupname=c(rep("A",4),"C",rep("B",6),"C","A"),
    Clonedef=c(rep("Clone_1a", 3), "Clone_2a", "Clone_3a",
               rep("Clone_1b", 3), "Clone_2b", "Clone_3b","Clone_4b",
               "Clone_1c","Clone_2c")
  )
  Short_DF=data.table::as.data.table(Short_DF)

  Header=c("AA_Whole_Replacement_muts_counts",
           "NT_CDR3_length")


  expect_equal(c(1:13),
               AbSolution:::filter_merged(
                 FBM=FBM,
                 merged_df=Short_DF,
                 merged_header=Header,
                 use_rgermline=FALSE,
                 use_repertoire=TRUE,
                 use_productive=TRUE,
                 use_nonproductive= FALSE,
                 my_regions= c("CDR3"),
                 my_var_elements=c("NT"),
                 my_vars=c(
                   "Length"
                 ),
                 my_vartypes="Baseline",
                 groups="Groupname",
                 group_A="A",
                 group_B="B",
                 group_C=c("C"),
                 univlog=T,
                 samples_to_keep=c("Pat_a",
                                   "Pat_b",
                                   "Pat_c"),
                 variables_to_remove=NULL,
                 pval_type="Bonferroni",
                 pval_cutoff=0.05,
                 estimate_cutoff=1,
                 number_selected_vars="All",
                 VJ_deselected= NULL,
                 R_mut_threshold_min=0,
                 R_mut_threshold_max=5,
                 V_J_to_use=c("VJ","VJ2"),
                 to_compare_groups=T,
                 use_sharedVDJ=T,
                 VDJ_normalized_per_size=F,
                 VDJ_maximize_clones=F,
                 VDJ_normalized_per_sample=F,
                 my_clone_def="Clonedef",
                 seed=1234,
                 chains="IGH",
                 igsubtypes=NULL
               )$ROWS)
})


test_that("Filtering used shared-VJ normalized per group works", {

  FBM=bigstatsr::as_FBM(matrix(c(rep(c(0:5),2),4,rep(14,4),15,rep(15,8)),
                               13)
  )

  Short_DF=data.frame(
    Sequence_type=rep(c("Repertoire"), 13),
    Patient_Sample=c(rep(c("Pat_a","Pat_b"), each=5), "Pat_b", rep("Pat_c",2)),
    Chain=rep("IGH",13),
    V_and_J=c(rep("VJ",3),rep("VJ2",2),rep("VJ",5),"VJ2",rep("VJ",2)),
    Groupname=c(rep("A",4),"C",rep("B",6),"C","A"),
    Clonedef=c(rep("Clone_1a", 3), "Clone_2a", "Clone_3a",
               rep("Clone_1b", 3), "Clone_2b", "Clone_3b","Clone_4b",
               "Clone_1c","Clone_2c")
  )
  Short_DF=data.table::as.data.table(Short_DF)

  Header=c("AA_Whole_Replacement_muts_counts",
           "NT_CDR3_length")


  expect_equal(c(c(1:6),c(8:13)),
               AbSolution:::filter_merged(
                 FBM=FBM,
                 merged_df=Short_DF,
                 merged_header=Header,
                 use_rgermline=FALSE,
                 use_repertoire=TRUE,
                 use_productive=TRUE,
                 use_nonproductive= FALSE,
                 my_regions= c("CDR3"),
                 my_var_elements=c("NT"),
                 my_vars=c(
                   "Length"
                 ),
                 my_vartypes="Baseline",
                 groups="Groupname",
                 group_A="A",
                 group_B="B",
                 group_C=c("C"),
                 univlog=T,
                 samples_to_keep=c("Pat_a",
                                   "Pat_b",
                                   "Pat_c"),
                 variables_to_remove=NULL,
                 pval_type="Bonferroni",
                 pval_cutoff=0.05,
                 estimate_cutoff=1,
                 number_selected_vars="All",
                 VJ_deselected= NULL,
                 R_mut_threshold_min=0,
                 R_mut_threshold_max=5,
                 V_J_to_use=c("VJ","VJ2"),
                 to_compare_groups=T,
                 use_sharedVDJ=T,
                 VDJ_normalized_per_size=T,
                 VDJ_maximize_clones=F,
                 VDJ_normalized_per_sample=F,
                 my_clone_def="Clonedef",
                 seed=1234,
                 chains="IGH",
                 igsubtypes=NULL
               )$ROWS)
})


test_that("Filtering used shared-VJ normalized per group and sample works", {

  FBM=bigstatsr::as_FBM(matrix(c(rep(c(0:5),2),4,rep(14,4),15,rep(15,8)),
                               13)
  )

  Short_DF=data.frame(
    Sequence_type=rep(c("Repertoire"), 13),
    Patient_Sample=c(rep(c("Pat_a","Pat_b"), each=5), "Pat_b", rep("Pat_c",2)),
    Chain=rep("IGH",13),
    V_and_J=c(rep("VJ",3),rep("VJ2",2),rep("VJ",5),"VJ2",rep("VJ",2)),
    Groupname=c(rep("A",4),"C",rep("B",6),"C","A"),
    Clonedef=c(rep("Clone_1a", 3), "Clone_2a", "Clone_3a",
               rep("Clone_1b", 3), "Clone_2b", "Clone_3b","Clone_4b",
               "Clone_1c","Clone_2c")
  )
  Short_DF=data.table::as.data.table(Short_DF)

  Header=c("AA_Whole_Replacement_muts_counts",
           "NT_CDR3_length")


  expect_equal(c(c(2,9),c(10,13)),
               AbSolution:::filter_merged(
                 FBM=FBM,
                 merged_df=Short_DF,
                 merged_header=Header,
                 use_rgermline=FALSE,
                 use_repertoire=TRUE,
                 use_productive=TRUE,
                 use_nonproductive= FALSE,
                 my_regions= c("CDR3"),
                 my_var_elements=c("NT"),
                 my_vars=c(
                   "Length"
                 ),
                 my_vartypes="Baseline",
                 groups="Groupname",
                 group_A="A",
                 group_B="B",
                 group_C=c("C"),
                 univlog=T,
                 samples_to_keep=c("Pat_a",
                                   "Pat_b",
                                   "Pat_c"),
                 variables_to_remove=NULL,
                 pval_type="Bonferroni",
                 pval_cutoff=0.05,
                 estimate_cutoff=1,
                 number_selected_vars="All",
                 VJ_deselected= NULL,
                 R_mut_threshold_min=0,
                 R_mut_threshold_max=5,
                 V_J_to_use=c("VJ","VJ2"),
                 to_compare_groups=T,
                 use_sharedVDJ=T,
                 VDJ_normalized_per_size=T,
                 VDJ_maximize_clones=F,
                 VDJ_normalized_per_sample=T,
                 my_clone_def="Clonedef",
                 seed=1234,
                 chains="IGH",
                 igsubtypes=NULL
               )$ROWS)
})


test_that("Filtering used shared-VJ normalized per group and clone works", {

  FBM=bigstatsr::as_FBM(matrix(c(rep(c(0:5),2),4,rep(14,4),15,rep(15,8)),
                               13)
  )

  Short_DF=data.frame(
    Sequence_type=rep(c("Repertoire"), 13),
    Patient_Sample=c(rep(c("Pat_a","Pat_b"), each=5), "Pat_b", rep("Pat_c",2)),
    Chain=rep("IGH",13),
    V_and_J=c(rep("VJ",3),rep("VJ2",2),rep("VJ",5),"VJ2",rep("VJ",2)),
    Groupname=c(rep("A",4),"C",rep("B",6),"C","A"),
    Clonedef=c(rep("Clone_1a", 3), "Clone_2a", "Clone_3a",
               rep("Clone_1b", 3), "Clone_2b", "Clone_3b","Clone_4b",
               "Clone_1c","Clone_2c")
  )
  Short_DF=data.table::as.data.table(Short_DF)

  Header=c("AA_Whole_Replacement_muts_counts",
           "NT_CDR3_length")


  expect_equal(c(c(1:7),c(9:13)),
               AbSolution:::filter_merged(
                 FBM=FBM,
                 merged_df=Short_DF,
                 merged_header=Header,
                 use_rgermline=FALSE,
                 use_repertoire=TRUE,
                 use_productive=TRUE,
                 use_nonproductive= FALSE,
                 my_regions= c("CDR3"),
                 my_var_elements=c("NT"),
                 my_vars=c(
                   "Length"
                 ),
                 my_vartypes="Baseline",
                 groups="Groupname",
                 group_A="A",
                 group_B="B",
                 group_C=c("C"),
                 univlog=T,
                 samples_to_keep=c("Pat_a",
                                   "Pat_b",
                                   "Pat_c"),
                 variables_to_remove=NULL,
                 pval_type="Bonferroni",
                 pval_cutoff=0.05,
                 estimate_cutoff=1,
                 number_selected_vars="All",
                 VJ_deselected= NULL,
                 R_mut_threshold_min=0,
                 R_mut_threshold_max=5,
                 V_J_to_use=c("VJ","VJ2"),
                 to_compare_groups=T,
                 use_sharedVDJ=T,
                 VDJ_normalized_per_size=T,
                 VDJ_maximize_clones=T,
                 VDJ_normalized_per_sample=F,
                 my_clone_def="Clonedef",
                 seed=1234,
                 chains="IGH",
                 igsubtypes=NULL
               )$ROWS)
})


test_that("Filtering used shared-VJ normalized per group, sample and clone works", {

  FBM=bigstatsr::as_FBM(matrix(c(rep(c(0:5),2),4,rep(14,4),15,rep(15,8)),
                               13)
  )

  Short_DF=data.frame(
    Sequence_type=rep(c("Repertoire"), 13),
    Patient_Sample=c(rep(c("Pat_a","Pat_b"), each=5), "Pat_b", rep("Pat_c",2)),
    Chain=rep("IGH",13),
    V_and_J=c(rep("VJ",3),rep("VJ2",2),rep("VJ",5),"VJ2",rep("VJ",2)),
    Groupname=c(rep("A",4),"C",rep("B",6),"C","A"),
    Clonedef=c(rep("Clone_1a", 3), "Clone_2a", "Clone_3a",
               rep("Clone_1b", 3), "Clone_2b", "Clone_3b","Clone_4b",
               "Clone_1c","Clone_2c")
  )
  Short_DF=data.table::as.data.table(Short_DF)

  Header=c("AA_Whole_Replacement_muts_counts",
           "NT_CDR3_length")


  expect_equal(c(c(1,6),c(9,12), 13),
               AbSolution:::filter_merged(
                 FBM=FBM,
                 merged_df=Short_DF,
                 merged_header=Header,
                 use_rgermline=FALSE,
                 use_repertoire=TRUE,
                 use_productive=TRUE,
                 use_nonproductive= FALSE,
                 my_regions= c("CDR3"),
                 my_var_elements=c("NT"),
                 my_vars=c(
                   "Length"
                 ),
                 my_vartypes="Baseline",
                 groups="Groupname",
                 group_A="A",
                 group_B="B",
                 group_C=c("C"),
                 univlog=T,
                 samples_to_keep=c("Pat_a",
                                   "Pat_b",
                                   "Pat_c"),
                 variables_to_remove=NULL,
                 pval_type="Bonferroni",
                 pval_cutoff=0.05,
                 estimate_cutoff=1,
                 number_selected_vars="All",
                 VJ_deselected= NULL,
                 R_mut_threshold_min=0,
                 R_mut_threshold_max=5,
                 V_J_to_use=c("VJ","VJ2"),
                 to_compare_groups=T,
                 use_sharedVDJ=T,
                 VDJ_normalized_per_size=T,
                 VDJ_maximize_clones=T,
                 VDJ_normalized_per_sample=T,
                 my_clone_def="Clonedef",
                 seed=1234,
                 chains="IGH",
                 igsubtypes=NULL
               )$ROWS)
})

