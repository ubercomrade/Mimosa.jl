# MotifHORDE Downstream Contract

Mimosa.jl owns motif representation, I/O, scanning, profile comparison, site
workflows, statistics, portable storage, cache primitives, and the CLI adapter.
Downstream orchestration owns discovery-tool execution, parameter grids, model
selection, reruns, and pipeline directory layout.

## Stable Workflows

```julia
model = readmodel(path)
writemodel(path, model)
batch, names = readsequences(fasta_path)

scores = scan(model, batch; strands=BestStrand(), execution=SerialExecution())
result = compare(query, target, batch; metric=:co)

prepared = prepare_profile(score_profile)
results = compare(prepared, targets; execution=ThreadedExecution(4))

sites = selectsites(model, batch, selector; execution=SerialExecution())
pfm = reconstruct_pfm(model, batch, selector; execution=SerialExecution())

null_result = build_null(
    models; sequences=batch, metric=:co, n_samples=2000, shuffle=true, seed=127
)
savenull(path, null_result.distribution)
dist = loadnull(path)
annotated = annotate_results(results, dist)
```

Supported comparison metrics are `co`, `co_rowwise`, `dice`, `dice_rowwise`,
and `cosine`. Model comparison always requires sequences. Null strategy is
always `"profile"`.

Downstream code must not import private names, depend on source-file layout,
mutate model structs, bypass validation, assume a global cache, or infer retired
matrix-comparison APIs from historical documents.

Downstream custom models subtype `AbstractMotifModel` and implement
`modelname`, `motif_length`, and `scan_kernel!` (plus
`left_context`/`right_context` for context models). They must not rely
on fields named `name`, `representation`, `weights`, `order`, or
`span`. See [Extending Mimosa](extending.md) for the full contract.

The separate package under `test/downstream/` verifies this contract
using public exports only, including a custom model defined in a separate
downstream module.
