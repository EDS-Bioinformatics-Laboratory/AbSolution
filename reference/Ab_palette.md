# Ab_palette: color ranges adapted to the BCR/TCR gene structure

Returns a vector of colors according to the input. For V/D/J/VJ/VDJ
inputs, the colors are produced in a way that genes and combination of
genes from the same family are assigned a similar color.

## Usage

``` r
Ab_palette(
  list_values,
  vect_genes_comb = NA,
  type_values = c("V", "D", "J", "VJ", "VDJ", "cuantitative", "cualitative"),
  colorblind = F,
  seed = 1234
)
```

## Arguments

- list_values:

  For V(D)J combinations, a list of lists of the families of each
  segment and the genes within them. For V/D/J segments, a list of the
  families and the genes within them. For cuantitative or cualitative
  values, a list of values.

- vect_genes_comb:

  The vector of the V(D)J combinations present in the dataset if
  type_values is 'VJ' or 'VDJ'. Otherwise, NA.

- type_values:

  One of 'V','D','J','VJ','VDJ', 'cuantitative' or 'cualitative'.
  'V','D','J','VJ','VDJ',

- colorblind:

  If TRUE, the output is a colorblind-friendly vector of colors from the
  viridis package. The similarity of the V-D-J colors is lost.

## Value

palette_colors: a (named) vector of colors
