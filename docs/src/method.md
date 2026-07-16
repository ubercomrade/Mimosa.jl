# Method

Mimosa compares heterogeneous motif models by the score profiles they produce
on the same DNA sequences. It does not convert their internal parameters to a
common matrix representation.

## Profile comparison pipeline

1. Scan the comparison sequences with each motif to obtain score profiles.
2. Convert raw scores to empirical `-log10(FPR)` values. A separate background
   sequence set can be used for calibration; otherwise the comparison sequences
   provide the empirical background.
3. Select one best anchor per sequence when `min_logfpr == 0`, or every anchor
   above the requested threshold.
4. Compare site-centered profile windows over the configured shifts and the
   four strand orientations `++`, `+-`, `-+`, and `--`.
5. Return the best score together with its offset, orientation, and number of
   contributing sites.

The offset is the query displacement relative to the target. A positive offset
means that the query is shifted to the right. Ties are resolved by score,
contributing site count, smaller absolute shift, and orientation order `++`,
`+-`, `-+`, `--`.

## Similarity metrics

For two aligned, non-negative profile windows `v1` and `v2`, the pooled
continuous overlap and Dice metrics are

```math
CO(v_1,v_2) =
\frac{\sum_i \min(v_{1i},v_{2i})}
     {\min(\sum_i v_{1i},\sum_i v_{2i})}
```

```math
Dice(v_1,v_2) =
\frac{2\sum_i \min(v_{1i},v_{2i})}
     {\sum_i v_{1i}+\sum_i v_{2i}}.
```

`co_rowwise` and `dice_rowwise` calculate the corresponding value for each
selected window and average finite values, so each anchor contributes equally.
`cosine` is also calculated per window and averaged:

```math
Cosine(v_1,v_2) =
\frac{v_1 \cdot v_2}{\lVert v_1 \rVert_2\lVert v_2 \rVert_2}.
```

All metrics are oriented so that a larger value means a stronger match. See
[Data Layout](data_layout.md) for coordinate and orientation conventions.

## Statistical significance

`build_null` compares eligible cross-group model pairs and pools their scores
into an empirical null sample. Mimosa fits a generalized extreme value (GEV)
distribution in `Float64`; the upper-tail survival probability of an observed
score is its p-value. `annotate_results` additionally computes
Benjamini-Hochberg adjusted p-values and E-values.

Null bundles are reusable only when their profile metric, comparison settings,
sequence and background fingerprints, and format contract are compatible with
the observed comparison. See [Storage Format](storage.md) for the null bundle
contract. Serial and threaded execution preserve result order and exact
comparison fields.

## References

1. Gupta S. et al. (2007). Quantifying similarity between motifs. *Genome
   Biology*, 8, R24. <https://doi.org/10.1186/gb-2007-8-2-r24>
2. Lambert S. A. et al. (2016). Motif comparison based on similarity of binding
   affinity profiles. *Bioinformatics*, 32(22), 3504-3506.
   <https://doi.org/10.1093/bioinformatics/btw489>
3. Siebert M. and Soding J. (2016). Bayesian Markov models consistently
   outperform PWMs at predicting motifs in nucleotide sequences. *Nucleic
   Acids Research*, 44(13), 6055-6069.
   <https://doi.org/10.1093/nar/gkw521>
