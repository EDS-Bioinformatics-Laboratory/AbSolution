# Big PCA

\*\*Internal function.\*\* Not intended for direct use. Exported only
for \`shinymeta\` report rendering via \`::\` access. Use \[run_app()\]
instead. A function to perform PCA over a FBM object.

## Usage

``` r
big_PCA(FBM, rows, columns)
```

## Arguments

- FBM:

  A FBM object.

- rows:

  Index of rows of the FBM that will be used to calculate the PCA.

- columns:

  Index of columns of the FBM that will be used to calculate the PCA..

## Value

A list with the PCA scores and its explained variances

## Examples

``` r
# \donttest{
  # Internal function exported for shinymeta :: access during report rendering.
  # Requires a live Shiny reactive context and real AIRR-seq data.
  # Use run_app() as the user-facing entry point.
# }
```
