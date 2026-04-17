# 4_Clonal_exploration

\*\*Internal function.\*\* Not intended for direct use. Exported only
for \`shinymeta\` report rendering via \`::\` access. Use \[run_app()\]
instead.

## Usage

``` r
draw_violinplots(
  seq_df,
  group = "Patient_Sample",
  selected_rows,
  clonotype,
  AA_or_NT = "AA",
  region = "CDR3",
  percentage = 100,
  freq_filter = 0,
  Selected_clones = NULL,
  dominance_threshold,
  seed = 1234,
  really_hide_dots = FALSE,
  width = 1400,
  height = 1000,
  img_type = "png",
  scale = 4,
  source = "clone_violinplot"
)
```

## Value

The return value, if any, from executing the utility.

## Examples

``` r
# \donttest{
  # Internal function exported for shinymeta :: access during report rendering.
  # Requires a live Shiny reactive context and real AIRR-seq data.
  # Use run_app() as the user-facing entry point.
# }
```
