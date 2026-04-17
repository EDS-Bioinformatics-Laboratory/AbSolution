# Merge FBMs and meta files

\*\*Internal function.\*\* Not intended for direct use. Exported only
for \`shinymeta\` report rendering via \`::\` access. Use \[run_app()\]
instead. A function to join by row individual FBMs and .info files into
a unique FBM and a unique .info file, merged_bm.

## Usage

``` r
merge_FBMs(filepath)
```

## Arguments

- filepath:

  Path to the individual FBMs.

## Value

A list with the merged FBM, the merged info file and the FBM header.

## Examples

``` r
# \donttest{
  # Internal function exported for shinymeta :: access during report rendering.
  # Requires a live Shiny reactive context and real AIRR-seq data.
  # Use run_app() as the user-facing entry point.
# }
```
