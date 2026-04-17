# 4_Clonal_exploration

\*\*Internal function.\*\* Not intended for direct use. Exported only
for \`shinymeta\` report rendering via \`::\` access. Use \[run_app()\]
instead.

## Usage

``` r
calculate_clone(
  seq_df,
  clonotype,
  AA_or_NT = "NT",
  region = "CDR3",
  percentage = 100,
  calculate_shared_clones
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
