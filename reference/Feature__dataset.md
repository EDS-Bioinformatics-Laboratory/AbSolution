# 2_Feature_determination_2

\*\*Internal function.\*\* Not intended for direct use. Exported only
for \`shinymeta\` report rendering via \`::\` access. Use \[run_app()\]
instead.

## Usage

``` r
Feature__dataset(
  path_base,
  DF_to_parse,
  name_DF_to_parse,
  FWR1partial,
  FWR4partial,
  standard = TRUE
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
