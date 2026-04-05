test_that("Ab_palette works for VDJ combinations", {

  library(alakazam)
  data("ExampleDb")
  vect_V_g=ExampleDb$v_call[1:10]
  vect_D_g=ExampleDb$d_call[1:10]
  vect_J_g=ExampleDb$j_call[1:10]

  vect_genes_comb=paste(vect_V_g,
                        vect_D_g,
                        vect_J_g,
                        sep="_")
  ##Example: Homsap IGHV3-11*05 F_Homsap IGHD3-10*01 F_Homsap IGHJ4*03 F

  list_values=list()
  list_values[["V"]]=list()
  list_values[["D"]]=list()
  list_values[["J"]]=list()
  for (com in c(1:length(vect_V_g))) {

    list_values[["V"]][[strsplit(strsplit(strsplit(vect_V_g[com],
                                                   split="*",
                                                   fixed=T)[[1]][1],
                                          split="/")[[1]][1],
                                 split="-")[[1]][1]]]= sort(unique(c(list_values[["V"]][[strsplit(strsplit(strsplit(vect_V_g[com], split="*", fixed=T)[[1]][1], split="/")[[1]][1], split="-")[[1]][1]]], vect_V_g[com])))
    list_values[["D"]][[strsplit(strsplit(strsplit(vect_D_g[com], split="*", fixed=T)[[1]][1], split="/")[[1]][1], split="-")[[1]][1]]]= sort(unique(c(list_values[["D"]][[strsplit(strsplit(strsplit(vect_D_g[com], split="*", fixed=T)[[1]][1], split="/")[[1]][1], split="-")[[1]][1]]], vect_D_g[com])))
    list_values[["J"]][[strsplit(strsplit(strsplit(vect_J_g[com], split="*", fixed=T)[[1]][1], split="/")[[1]][1], split="-")[[1]][1]]]= sort(unique(c(list_values[["J"]][[strsplit(strsplit(strsplit(vect_J_g[com], split="*", fixed=T)[[1]][1], split="/")[[1]][1], split="-")[[1]][1]]], vect_J_g[com])))
  }

  expect_equal(unname(Ab_palette(list_values=list_values,
                                 vect_genes_comb=vect_genes_comb,
                                 type_values="VDJ", colorblind=F)),
               c("#7F837A", "#75A9BF", "#87ACFE", "#87D4FE", "#FFFA79",
                 "#7E8EFE", "#7EFCFE", "#7ED4FE", "#7EACFE", "#7EA9FE"))

  expect_equal(names(Ab_palette(list_values=list_values,
                                vect_genes_comb=vect_genes_comb,
                                type_values="VDJ",
                                colorblind=F)),
               vect_genes_comb)
})




test_that("Ab_palette works for VJ combinations", {

  library(alakazam)
  data("ExampleDb")
  vect_V_g=ExampleDb$v_call[1:10]
  vect_J_g=ExampleDb$j_call[1:10]

  vect_genes_comb=paste(vect_V_g,
                        vect_J_g,
                        sep="_")
  ##Example: Homsap IGHV3-11*05 F_Homsap IGHJ4*03 F

  list_values=list()
  list_values[["V"]]=list()
  list_values[["J"]]=list()
  for (com in c(1:length(vect_V_g))) {
    list_values[["V"]][[strsplit(strsplit(strsplit(vect_V_g[com], split="*", fixed=T)[[1]][1], split="/")[[1]][1], split="-")[[1]][1]]]= sort(unique(c(list_values[["V"]][[strsplit(strsplit(strsplit(vect_V_g[com], split="*", fixed=T)[[1]][1], split="/")[[1]][1], split="-")[[1]][1]]], vect_V_g[com])))
    list_values[["J"]][[strsplit(strsplit(strsplit(vect_J_g[com], split="*", fixed=T)[[1]][1], split="/")[[1]][1], split="-")[[1]][1]]]= sort(unique(c(list_values[["J"]][[strsplit(strsplit(strsplit(vect_J_g[com], split="*", fixed=T)[[1]][1], split="/")[[1]][1], split="-")[[1]][1]]], vect_J_g[com])))
  }

  expect_equal(unname(Ab_palette(list_values=list_values,
                                 vect_genes_comb=vect_genes_comb,
                                 type_values="VJ",
                                 colorblind=F)),
               c( "#7F7A8E", "#757ABF", "#877AF8", "#877AF8", "#FF7A83",
                  "#7E7AF8", "#7E7AF8", "#7E7AF8", "#7E7AF8", "#7E7AF8"))

  expect_equal(names(Ab_palette(list_values=list_values,
                                vect_genes_comb=vect_genes_comb,
                                type_values="VJ",
                                colorblind=F)),
               vect_genes_comb)
})

test_that("Ab_palette works for V (ergo also D or J) gene sets", {

  library(alakazam)
  data("ExampleDb")
  vect_V_g=ExampleDb$v_call[1:10]


  list_values=list()
  list_values[["V"]]=list()
  for (com in c(1:length(vect_V_g))) {
    list_values[["V"]][[strsplit(strsplit(strsplit(vect_V_g[com], split="*", fixed=T)[[1]][1], split="/")[[1]][1], split="-")[[1]][1]]]= sort(unique(c(list_values[["V"]][[strsplit(strsplit(strsplit(vect_V_g[com], split="*", fixed=T)[[1]][1], split="/")[[1]][1], split="-")[[1]][1]]], vect_V_g[com])))
  }

  expect_equal(unname(Ab_palette(list_values=list_values[["V"]],
                                 vect_genes_comb=NA,
                                 type_values="V",
                                 colorblind=F)),
               c("#7FAA7F","#7F9F7F", "#7F8A7F", "#7FA97F","#4DFFFF" ))

  expect_equal(names(Ab_palette(list_values=list_values[["V"]],
                                vect_genes_comb=NA,
                                type_values="V",
                                colorblind=F)),
               unname(unlist(list_values[["V"]])))
})

test_that("Ab_palette works for quantitative", {
  list_values=c(1,10,4,5,0)
  expect_equal(unname(Ab_palette(list_values=list_values,
                                 vect_genes_comb=NA,
                                 type_values="cuantitative",
                                 colorblind=F)),
               c("#D1495B", "#D55958", "#DA6A55", "#DF7C52", "#E38C4F",
                 "#E89D4C", "#EDAE49", "#C0A95B", "#93A56E", "#66A182"))
})

test_that("Ab_palette works for qualitative", {
  list_values=c("Group_1","Group_2")
  expect_equal(unname(Ab_palette(list_values=list_values,
                                 vect_genes_comb=NA,
                                 type_values="cualitative",
                                 colorblind=F)),
               c("#00798c", "#9e51bd"))
})

test_that("Ab_palette colorblind mode works for V(D)J combinations/quantitative/qualitative", {
  list_values=c("Group_1","Group_2")
  expect_equal(unname(Ab_palette(list_values=list_values,
                                 vect_genes_comb=NA,
                                 type_values="cualitative",
                                 colorblind=T)),
               c("#440154FF", "#FDE725FF"))

  list_values=c(1,10,4,5,0)
  expect_equal(unname(Ab_palette(list_values=list_values,
                                 vect_genes_comb=NA,
                                 type_values="cuantitative",
                                 colorblind=T)),
               c("#440154FF", "#482878FF", "#3E4A89FF", "#31688EFF",
                 "#26828EFF", "#1F9E89FF", "#35B779FF", "#6DCD59FF",
                 "#B4DE2CFF", "#FDE725FF"))


  library(alakazam)
  data("ExampleDb")
  vect_V_g=ExampleDb$v_call[1:10]
  vect_D_g=ExampleDb$d_call[1:10]
  vect_J_g=ExampleDb$j_call[1:10]

  vect_genes_comb=paste(vect_V_g,
                        vect_D_g,
                        vect_J_g,
                        sep="_")
  ##Example: Homsap IGHV3-11*05 F_Homsap IGHD3-10*01 F_Homsap IGHJ4*03 F

  list_values=list()
  list_values[["V"]]=list()
  list_values[["D"]]=list()
  list_values[["J"]]=list()
  for (com in c(1:length(vect_V_g))) {

    list_values[["V"]][[strsplit(strsplit(strsplit(vect_V_g[com],
                                                   split="*",
                                                   fixed=T)[[1]][1],
                                          split="/")[[1]][1],
                                 split="-")[[1]][1]]]= sort(unique(c(list_values[["V"]][[strsplit(strsplit(strsplit(vect_V_g[com], split="*", fixed=T)[[1]][1], split="/")[[1]][1], split="-")[[1]][1]]], vect_V_g[com])))
    list_values[["D"]][[strsplit(strsplit(strsplit(vect_D_g[com], split="*", fixed=T)[[1]][1], split="/")[[1]][1], split="-")[[1]][1]]]= sort(unique(c(list_values[["D"]][[strsplit(strsplit(strsplit(vect_D_g[com], split="*", fixed=T)[[1]][1], split="/")[[1]][1], split="-")[[1]][1]]], vect_D_g[com])))
    list_values[["J"]][[strsplit(strsplit(strsplit(vect_J_g[com], split="*", fixed=T)[[1]][1], split="/")[[1]][1], split="-")[[1]][1]]]= sort(unique(c(list_values[["J"]][[strsplit(strsplit(strsplit(vect_J_g[com], split="*", fixed=T)[[1]][1], split="/")[[1]][1], split="-")[[1]][1]]], vect_J_g[com])))
  }

  expect_equal(unname(Ab_palette(list_values=list_values,
                                 vect_genes_comb=vect_genes_comb,
                                 type_values="VDJ", colorblind=T)),
               c("#440154FF", "#482878FF", "#3E4A89FF", "#31688EFF",
                 "#26828EFF",  "#1F9E89FF", "#35B779FF",  "#6DCD59FF",
                 "#B4DE2CFF", "#FDE725FF"))


  vect_V_g=ExampleDb$v_call[1:10]


  list_values=list()
  list_values[["V"]]=list()
  for (com in c(1:length(vect_V_g))) {
    list_values[["V"]][[strsplit(strsplit(strsplit(vect_V_g[com], split="*", fixed=T)[[1]][1], split="/")[[1]][1], split="-")[[1]][1]]]= sort(unique(c(list_values[["V"]][[strsplit(strsplit(strsplit(vect_V_g[com], split="*", fixed=T)[[1]][1], split="/")[[1]][1], split="-")[[1]][1]]], vect_V_g[com])))
  }

  expect_equal(unname(Ab_palette(list_values=list_values[["V"]],
                                 vect_genes_comb=NA,
                                 type_values="V",
                                 colorblind=T)),
               c("#440154FF", "#3B528BFF", "#21908CFF",
                 "#5DC863FF", "#FDE725FF"))
})


