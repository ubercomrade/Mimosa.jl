# Julia API Reference

## Model types

```@docs
AbstractProfileSource
AbstractMotifModel
PWM
BaMM
SiteGA
Dimont
Slim
ScoreProfile
is_scannable
pfm_to_pwm
pwm_from_pfm
extend_pwm_with_n
```

## Model geometry and extension contract

```@docs
modelname
motif_length
left_context
right_context
window_size
npositions
site_start_offset
scan_kernel!
validate_model
```

## Model I/O

```@docs
readmodel
writemodel
read_scores
read_bamm
read_sitega
write_sitega
read_dimont
read_slim
readsequences
read_fasta
```

## Scanning

```@docs
scan
scan!
scan_forward!
scan_reverse!
best_hits!
scan_both!
scan_result_lengths
scorebounds
scorematrix
scoretype
StrandPair
```

## Sequence representation

```@docs
EncodedSequenceBatch
nsequences
seqlength
sequence
empty_sequence_batch
N_CODE
encode_base
encode_sequence
reverse_complement
reverse_complement!
to_padded
from_padded
make_random_sequences
RaggedArray
nrows
rowlength
row
build_ragged
empty_ragged
```

## Strand policies

```@docs
StrandPolicy
ForwardOnly
ReverseOnly
BestStrand
BothStrands
```

## Profile comparison

```@docs
compare
ComparisonResult
metric_name
AbstractProfileMetric
OverlapCoefficient
OverlapCoefficientRowwise
DiceSimilarity
DiceSimilarityRowwise
CosineSimilarityProfile
ProfileConfig
PreparedProfile
prepare_profile
profile_bundle
parse_profile_metric
LogTailTable
EmpiricalLogTail
fit
lookup_score
transform_scores
flatten_bundle
normalize_bundle
AnchorCSR
build_anchor_csr
collect_best_anchors
collect_threshold_anchors
collect_anchors
score_shift
profile_compare
```

## Site extraction

```@docs
SiteSelector
BestPerSequence
ThresholdHits
TopFractionHits
SiteHit
SiteCollection
selectsites
reconstruct_pfm
extract_site_matrix
site_strings
build_pcm
pcm_to_pfm
sort_hits!
select_top_fraction
empty_site_collection
```

## Statistics

```@docs
GEVFit
GEVFitFailure
GEVFitResult
fit_gev
survival
cdf
pvalue
scipy_params
BenjaminiHochberg
adjusted_pvalues
evalue
NullPair
NullDistribution
ProfileComparisonContract
NullBuildConfig
NullBuildResult
build_null
annotate_results
AnnotatedResult
eligible_targets
parse_group_relations
GroupRelations
savenull
loadnull
MODEL_FORMAT_VERSION
NULL_FORMAT_VERSION
ANNOTATED_RESULT_SCHEMA_VERSION
```

## Execution policies

```@docs
ExecutionPolicy
SerialExecution
ThreadedExecution
```

## Cache

```@docs
Cache
cache_key
cache_has
cache_get
cache_get_meta
cache_set
clearcache
content_fingerprint
model_fingerprint
model_collection_fingerprint
sequence_fingerprint
```

## Serialization

```@docs
to_json
to_dict
main
```

## Errors

```@docs
MimosaError
ModelFormatError
ModelDimensionError
InvariantError
ModelInterfaceError
```
