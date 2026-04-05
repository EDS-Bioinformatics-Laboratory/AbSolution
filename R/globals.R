# Suppress R CMD check notes for data.table NSE and reactive variables
utils::globalVariables(c(
  # data.table
  ".", ".N", ".SD", "N",
  # columnas creadas on-the-fly
  "Clone", "Clone_ID", "color",
  # reactiveValues usados en módulos
  "Selection_values", "Big_mem_values"
))
