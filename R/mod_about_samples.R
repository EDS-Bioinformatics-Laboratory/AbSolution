#' about_samples UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @import shiny
#' @importFrom shinycssloaders withSpinner
mod_about_samples_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
      column(
        5, htmlOutput("sunburst_gene") %>% withSpinner(type=6, color="#cf395c")
      ),
      column(7, DT::dataTableOutput(ns("table_gene"))%>%
               withSpinner(type = 6, color = "#cf395c"))
    )
  )
}

#' about_samples Server Functions
#' @import  sunburstR
#' @noRd
mod_about_samples_server <- function(id, data, group_by, show, gene_combination){
  moduleServer( id, function(input, output, session){
    ns <- session$ns


    output$sunburst_gene <- renderUI({
      if(!is.null(Selection_values$rows) && length(Selection_values$rows)>0 && !is.null(Selection_values$columns) && length(Selection_values$columns)>0) {
        tmp_sunburst=unname(unlist(Big_mem_values$Header[Selection_values$columns]))
        # regions=c("Whole","FWR1","CDR1","FWR2","CDR2","FWR3","CDR3","FWR4")
        region="FWR1"
        sunburst_df=show_selected_features(tmp_sunburst, region=region)
        # add_shiny(sund2b(sunburst_df,showLabels =F, rootLabel= "Variables analyzed per each region in the selection set") )
        if(dim(sunburst_df)[1]==0) {
          HTML(paste(c("There is no", region, "sequence variables analyzed in the selection set"), collapse=" "))
        } else {
          sund2b(sunburst_df,showLabels =F, rootLabel= paste(c(region, "sequence variables analyzed in the selection set"), collapse=" "))
        }
      }
    })

  })
}

## To be copied in the UI
# mod_about_samples_ui("about_samples_1")

## To be copied in the server
# mod_about_samples_server("about_samples_1")
