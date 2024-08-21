#' Ab_palette: color ranges adapted to the BCR/TCR gene structure
#'
#' @description Returns a vector of colors according to the input.
#'    For V/D/J/VJ/VDJ inputs, the colors are produced in a way that genes and
#'    combination of genes from the same family are assigned a similar color.
#' @param list_values For V(D)J combinations, a list of lists of the families of
#'     each segment and the genes within them. For V/D/J segments, a list of the
#'      families and the genes within them. For cuantitative or cualitative
#'      values, a list of values.
#' @param vect_genes_comb The vector of the V(D)J combinations present in the
#'    dataset if type_values is "VJ" or "VDJ". Otherwise, NA.
#' @param type_values One of "V","D","J","VJ","VDJ", "cuantitative" or
#'    "cualitative".  "V","D","J","VJ","VDJ",
#' @param colorblind If TRUE, the output is a colorblind-friendly vector of
#'    colors from the viridis package. The similarity of the V-D-J colors is
#'    lost.
#'
#' @return palette_colors: a (named) vector of colors
#' @import stats
#' @import grDevices
#' @import viridis
#' @noRd
Ab_palette=function(list_values, vect_genes_comb=NA, type_values=c("V","D","J","VJ","VDJ", "cuantitative","cualitative"), colorblind=F){

  palette_colors=c()
  set.seed(123)

  if(colorblind){
    if (type_values=="cuantitative") {
      palette_colors=viridis(max(list_values)-min(list_values))
    } else if (type_values=="cualitative" || type_values %in% c("V","D","J")) {
      palette_colors=viridis(length(unlist(unlist(list_values))))
      names(palette_colors)=sort(unlist(unlist(list_values)))
    } else {
      palette_colors=viridis(length(unlist(unlist(vect_genes_comb))))
      names(palette_colors)=sort(unlist(unlist(vect_genes_comb)))
    }
  } else {
    if(type_values =="VJ" || type_values == "VDJ") {
      n_genes=length(list_values)
      RGB_list=list()

      for (k in c(1:n_genes)) {
        n_family = length(list_values[[k]])
        divisions=seq(0, 255, length.out=n_family)
        # RGB_list[[names(list_values)[k]]]=list()

        for (i in c(1:n_family)){
          # RGB_list[[names(list_values)[k]]][[as.character(i)]]=list()
          for (j in c(1:length(list_values[[k]][[i]]))) {
            RGB_combn_exists=T
            RGB_val=divisions[i]
            tmp_RGB_val=RGB_val
            while(RGB_combn_exists){

              tmp_inc=rnorm(1, mean=5.5, sd=10)
              if(tmp_RGB_val %in% unlist(RGB_list)) {
                tmp_RGB_val=RGB_val
                if(sample(1:2)[1] == 1 ) {
                  tmp_RGB_val=max(min(tmp_RGB_val+tmp_inc, 255), 0)

                } else {
                  tmp_RGB_val=min(max(tmp_RGB_val-tmp_inc, 0), 255)
                }

              } else {
                RGB_combn_exists=F
              }

            }

            RGB_list[[as.character(sort(list_values[[k]][[i]])[[j]])]]=tmp_RGB_val
          }

        }

      }

      if(type_values =="VJ") {
        palette_colors=sapply(vect_genes_comb, function(z) grDevices::rgb(RGB_list[[strsplit(z, split="_")[[1]][1]]],
                                                               122,
                                                               RGB_list[[strsplit(z, split="_")[[1]][2]]],
                                                               max=255) )
        # seecol(usecol(pal = palette_colors, n = "all"))

      } else {
        palette_colors=sapply(vect_genes_comb,
                              function(z) if(grepl("__",z)){
                                grDevices::rgb(RGB_list[[strsplit(z, split="_")[[1]][1]]],
                                    RGB_list[[ which(names(RGB_list) =="")]],
                                    RGB_list[[strsplit(z, split="_")[[1]][3]]],
                                    max=255)
                              } else {
                                grDevices::rgb(RGB_list[[strsplit(z, split="_")[[1]][1]]],
                                    RGB_list[[strsplit(z, split="_")[[1]][2]]],
                                    RGB_list[[strsplit(z, split="_")[[1]][3]]],
                                    max=255)
                              } )
        # seecol(usecol(pal = palette_colors, n = "all"))

      }


      # sapply(vect_genes_comb, function(z) RGB_list[[strsplit(z, split="__")[[1]][1]]] )


    }  else if (type_values == "cualitative"){
      some_colors <- list(
        "blue"   = "#00798c",
        "lavender" ="#9e51bd",
        "rosedust" = "#9D5568",
        "red"    = "#d1495b",
        "orange" = "#DF7C52"  ,
        "yellow" = "#edae49",
        "green"  = "#66a182",
        "greenmoss" ="#AAA866",
        "greeeeen"   = "#2e853b",
        "grey"   = "#333840",
        "pinkputi" = "#F72585",
        "azulon" = "#3A0CA3",
        "azulin" = "#4CC9F0",
        "raro" = "#D81159",
        "naranja" ="#f86624",
        "gold" ="#ff7d00",
        "green" = "#007200",
        "darkgreen" = "#004b23",
        "d81159"="#d81159",
        "verypink"="#ff0a54",
        "goldie"="#fbb13c",
        "softmalva"="#9381ff",
        "softyellow"="#fee440",
        "lightgreen"="#abff4f",
        "keyred"="#e5383b",
        "brown"="#78290f",
        "brownred"="#ad2831",
        "verybrown"="#941b0c",
        "darkbrown"="#220901",
        "opaqueblack"="#242423",
        "lightlessblue"="#003d5b",
        "verylightbrown"="#d69f7e",
        "veryldarkblue"="#04052e",
        "clodublue"="#50d8d7",
        "purpura"="#4f0147",
        "aqua"="#037971",
        "odoperdido"="#76520e",
        "cloudgreen"="#7bf1a8",
        "darkmagenta"="#290025",
        "powerpink"="#ff70a6",
        "lightgrey"="#d5bdaf"
      )

      for (element in c(1:length(list_values))) {
        if(element<=length(some_colors)) {

          new_palette=some_colors[[element]]
        } else {
          is_in_palette=T
          while(is_in_palette) {
            RGB=round(runif(3, min=-0.49, max=255.49))
            new_palette=grDevices::rgb(RGB[1],RGB[2],RGB[3], maxColorValue = 255)
            if(new_palette %!in% palette_colors && sum(sapply(RGB, function(z) 255-z))>15 && sum (RGB) > 20) {
              is_in_palette=F
            }
          }


        }
        palette_colors=c(palette_colors, new_palette)
      }
      names(palette_colors)=sort(list_values)

    } else if (type_values=="cuantitative"){
      some_colors <- list(
        "red"    = "#d1495b",
        "orange" = "#DF7C52"  ,
        "yellow" = "#edae49",
        "green"  = "#66a182")
      palette_colors=colorRampPalette(some_colors)(max(list_values)-min(list_values))
    } else {
      n_family = length(list_values)
      divisions=seq(0, 255, length.out=n_family)
      RGB_list=list()
      for (i in c(1:n_family)){
        RGB_combn_exists=T

        while(RGB_combn_exists){
          RGB_combn=round(runif(3, min=0.51, max=(n_family+0.49)))
          RGB_combn[1:3]=RGB_combn[1]
          RGB_combn[round(runif(1, min=0.51, max=3.49))]=NA
          RGB_combn=paste(RGB_combn, collapse = "_")
          if(RGB_combn %!in% RGB_list) {
            RGB_combn_exists=F
          }
        }

        RGB_list[[length(RGB_list)+1]]=RGB_combn
      }
      suppressWarnings({
        for (i in c(1:n_family)) {
          tmp_palette_colors=c()

          # fam_divisions=seq(0, 255, length.out=length(list_values[[i]]))
          fam_divisions=round(runif(1, min=-0.49, max=255.49))
          for (j in c(1:length(list_values[[i]]))) {
            tmp_RGB=as.numeric(strsplit(RGB_list[[i]], split="_")[[1]])
            tmp_RGB=divisions[tmp_RGB]
            NA_index=which(is.na(tmp_RGB))
            tmp_RGB[NA_index]=fam_divisions


            tmp_RGB_exists=T

            while(tmp_RGB_exists) {

              if(sum(sapply(tmp_RGB, function(z) 255-z))<15 || sum(tmp_RGB)<45 ||  grDevices::rgb(tmp_RGB[1],tmp_RGB[2],tmp_RGB[3], maxColorValue = 255) %in% palette_colors || grDevices::rgb(tmp_RGB[1],tmp_RGB[2],tmp_RGB[3], maxColorValue = 255) %in% tmp_palette_colors ){

                tmp_inc=rnorm(1, mean=10, sd=15)
                if(sample(1:2)[1] == 1 ) {
                  tmp_RGB[NA_index]=max(min(fam_divisions+tmp_inc, 255), 0)

                } else {
                  tmp_RGB[NA_index]=min(max(fam_divisions-tmp_inc, 0), 255)
                }
              } else {
                tmp_RGB_exists=F
                tmp_palette_colors=c(tmp_palette_colors, grDevices::rgb(tmp_RGB[1],tmp_RGB[2],tmp_RGB[3], maxColorValue = 255))

              }
            }
          }

          names(tmp_palette_colors)=list_values[[i]]

          palette_colors=c(palette_colors, tmp_palette_colors)
        }
      })



    }
  }

  return(palette_colors)
}



#' Total variance of the PCA
#'
#' @description (squared) Frobenius norm (which is the total variance of the
#'    matrix) https://github.com/privefl/bigstatsr/issues/83
#' @author Florian Privé
#' @param X FBM object used in the bigstatsr::big_randomSVD calculation
#' @param ind.row FBM rows used in the bigstatsr::big_randomSVD calculation
#' @param ind.col FBM cols used in the bigstatsr::big_randomSVD calculation
#' @param center centering vector used in bigstatsr::big_randomSVD
#' @param scale scaling vector used in bigstatsr::big_randomSVD
#' @return (squared) Frobenius norm
#' @import bigstatsr
#' @import bigassertr
#' @noRd
big_norm <- function(X, ind.row = rows_along(X), ind.col = cols_along(X),
                     center = rep(0, length(ind.col)),
                     scale = rep(1, length(ind.col))) {
  bigassertr::assert_lengths(center, scale, ind.col)
  stats <- bigstatsr::big_colstats(X, ind.row, ind.col)
  n <- length(ind.row)
  return(sum(((n - 1) * stats$var + n * (stats$sum / n - center)^2) / scale^2))
}


#' TEMPORAL
#'
#' @description colorss
#' @return color
#' @noRd

branded_colors <- list(
  "blue"   = "#00798c",
  "lavender" ="#9e51bd",
  "rosedust" = "#9D5568",
  "red"    = "#d1495b",
  "orange" = "#DF7C52"  ,
  "yellow" = "#edae49",
  "green"  = "#66a182",
  "greenmoss" ="#AAA866",
  "navy"   = "#2e853b",
  "grey"   = "#333840",
  "pinkputi" = "#F72585",
  "azulon" = "#3A0CA3",
  "azulin" = "#4CC9F0",
  "raro" = "#D81159",
  "naranja" ="#f86624",
  "gold" ="#ff7d00",
  "green" = "#007200",
  "darkgreen" = "#004b23",
  "d81159"="#d81159",
  "verypink"="#ff0a54",
  "goldie"="#fbb13c",
  "softmalva"="#9381ff",
  "softyellow"="#fee440",
  "lightgreen"="#abff4f",
  "keyred"="#e5383b",
  "brown"="#78290f",
  "brownred"="#ad2831",
  "verybrown"="#941b0c",
  "darkbrown"="#220901",
  "opaqueblack"="#242423",
  "lightlessblue"="#003d5b",
  "verylightbrown"="#d69f7e",
  "veryldarkblue"="#04052e",
  "clodublue"="#50d8d7",
  "purpura"="#4f0147",
  "aqua"="#037971",
  "odoperdido"="#76520e",
  "cloudgreen"="#7bf1a8",
  "darkmagenta"="#290025",
  "powerpink"="#ff70a6",
  "lightgrey"="#d5bdaf"
  # "black" ="#050609",
  # "pinkie"="#ff4d6d"
)
