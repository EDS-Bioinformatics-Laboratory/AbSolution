# 1_Parsing

\*\*Internal function.\*\* Not intended for direct use. Exported only
for \`shinymeta\` report rendering via \`::\` access. Use \[run_app()\]
instead.

## Usage

``` r
parse_AIRRSeq_file(
  file,
  group,
  patient,
  subgroup,
  sample,
  input_path,
  C_region_included,
  FWR1partial,
  FWR4partial,
  D_gene,
  repertoire,
  output_path,
  is_example = FALSE
)
```

## Value

The return value, if any, from executing the function.

## Examples

``` r
# \donttest{
  # Internal function exported for shinymeta :: access during report rendering.
  # Requires a live Shiny reactive context and real AIRR-seq data.
  # Use run_app() as the user-facing entry point.
# }
```
