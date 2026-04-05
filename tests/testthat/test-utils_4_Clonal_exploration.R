test_that("Clonal grouping works for VCDR3J grouping with non-shared clones", {
  expect_equal(c("VJ1_______A__Pat1","VJ1_______A__Pat2","VJ2_______A__Pat1"),
               AbSolution:::calculate_clone(
                 seq_df=data.table::as.data.table(data.frame(V_and_J=c("VJ1","VJ1","VJ2"),
                                                             AA_CDR3=c("A","A","A"),
                                                             Patient_Sample=c("Pat1","Pat2","Pat1"))),
                 clonotype="VCDR3J",
                 AA_or_NT = "AA",
                 region="CDR3",
                 percentage = 100,
                 calculate_shared_clones = F)$'Clone_AA_CDR3_VCDR3J_100_non-shared_simil')
})


test_that("Clonal grouping works for VCDR3J grouping with shared clones", {
  expect_equal(c("VJ1_______A","VJ1_______A","VJ2_______A"),
               AbSolution:::calculate_clone(
                 seq_df=data.table::as.data.table(data.frame(V_and_J=c("VJ1","VJ1","VJ2"),
                                                             AA_CDR3=c("A","A","A"),
                                                             Patient_Sample=c("Pat1","Pat2","Pat1"))),
                 clonotype="VCDR3J",
                 AA_or_NT = "AA",
                 region="CDR3",
                 percentage = 100,
                 calculate_shared_clones = T)$'Clone_AA_CDR3_VCDR3J_100_shared_simil')
})


test_that("Clonal grouping works for Reconstructed_germline grouping with shared clones", {
  expect_equal(rep(c("VJ1_______B","VJ1_______B","VJ2_______B"), each=2),
               AbSolution:::calculate_clone(
                 seq_df=data.table::as.data.table(data.frame(ID=rep(c("1","2","3"), each = 2),
                                                             Sequence_type=rep(c("Repertoire","Reconstructed_germline"),3),
                                                             V_and_J=rep(c("VJ1","VJ1","VJ2"), each = 2),
                                                             AA_Whole=rep(c("A","B"), 3),
                                                             Patient_Sample=rep(c("Pat1","Pat2","Pat1"), each = 2))),
                 clonotype="Reconstructed_germline",
                 AA_or_NT = "AA",
                 region="Whole",
                 percentage = 100,
                 calculate_shared_clones = T)$'Clone_AA_Whole_Reconstructed_germline_100_shared_simil')
})


test_that("Clonal grouping works for Reconstructed_germline grouping with shared clones", {
  expect_equal(rep(c("VJ1_______B__Pat1","VJ1_______B__Pat2","VJ2_______B__Pat1"), each=2),
               AbSolution:::calculate_clone(
                 seq_df=data.table::as.data.table(data.frame(ID=rep(c("1","2","3"), each = 2),
                                                             Sequence_type=rep(c("Repertoire","Reconstructed_germline"),3),
                                                             V_and_J=rep(c("VJ1","VJ1","VJ2"), each = 2),
                                                             AA_Whole=rep(c("A","B"), 3),
                                                             Patient_Sample=rep(c("Pat1","Pat2","Pat1"), each = 2))),
                 clonotype="Reconstructed_germline",
                 AA_or_NT = "AA",
                 region="Whole",
                 percentage = 100,
                 calculate_shared_clones = F)$'Clone_AA_Whole_Reconstructed_germline_100_non-shared_simil')
})

