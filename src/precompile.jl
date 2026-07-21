# Precompilation workload using PrecompileTools.
#
# This workload exercises representative code paths (model construction, PWM
# scanning, reverse complement, motif comparison, GEV fit, JSON serialization)
# during package precompilation — NOT at `using Mimosa` time.
#
# PrecompileTools ensures the workload code runs only during precompilation
# and has zero runtime cost when the module is loaded.

using PrecompileTools: @compile_workload

function _precompile_workload()
    # ── PWM construction and scanning ────────────────────────────────────
    # Small synthetic PWM: 5 rows (A,C,G,T,N) × 8 columns
    weights = Float32[
        0.5 -0.3 0.8 -0.2 0.1 0.6 -0.4 0.3
        -0.2 0.7 -0.5 0.3 0.8 -0.1 0.2 -0.6
        0.1 -0.4 0.2 0.6 -0.3 0.5 0.7 -0.2
        0.3 0.1 -0.6 0.4 0.2 -0.5 -0.1 0.8
        -0.2 -0.3 -0.5 -0.2 -0.3 -0.5 -0.4 -0.6
    ]
    bg = (Float32(0.25), Float32(0.25), Float32(0.25), Float32(0.25))
    pwm = PWM("precompile_motif", weights, bg)

    # Score bounds
    scorebounds(pwm)

    # Reverse complement
    rc = reverse_complement(weights)

    # Encode and scan a short sequence
    seq_str = "ACGTACGTACGTACGTACGT"
    seq_enc = encode_sequence(seq_str)
    scan(pwm, seq_enc; strands=ForwardOnly())
    scan(pwm, seq_enc; strands=ReverseOnly())
    scan(pwm, seq_enc; strands=BestStrand())
    scan(pwm, seq_enc; strands=BothStrands())

    # In-place scan
    dest = Vector{Float32}(undef, npositions(length(seq_enc), 8))
    scan!(dest, pwm, seq_enc; strands=ForwardOnly())

    # Batch scan
    batch = make_random_sequences(5, 100; seed=42)
    scan(pwm, batch; strands=BestStrand(), execution=Execution())

    # ── Motif comparison ────────────────────────────────────────────────
    compare(pwm, pwm, batch; metric=OverlapCoefficient(), min_logfpr=0.0f0)

    # ── Site extraction ─────────────────────────────────────────────────
    sites = selectsites(pwm, batch, BestPerSequence(); strands=BestStrand())

    # ── PFM reconstruction ───────────────────────────────────────────────
    if nsequences(batch) > 0
        reconstruct_pfm(pwm, batch, BestPerSequence(); pseudocount=Float32(1e-4))
    end

    # ── GEV fit ──────────────────────────────────────────────────────────
    samples = Float32[
        0.1,
        0.2,
        0.3,
        0.4,
        0.5,
        0.6,
        0.7,
        0.8,
        0.9,
        1.0,
        0.15,
        0.25,
        0.35,
        0.45,
        0.55,
        0.65,
        0.75,
        0.85,
        0.95,
        1.05,
    ]
    fit_gev(samples)

    # ── BH FDR ───────────────────────────────────────────────────────────
    pvals = Float32[0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.5, 0.9]
    adjusted_pvalues(pvals; method=BenjaminiHochberg())

    # ── JSON serialization ──────────────────────────────────────────────
    result = compare(pwm, pwm, batch; metric=:co, min_logfpr=0.0f0)
    to_json(result)
    to_dict(result)

    # ── Higher-order model scan (BaMM geometry) ────────────────────────
    # Small BaMM: order=0, 5×4 representation
    bamm_weights = Float32[
        0.1 -0.2 0.3 -0.1
        0.2 0.1 -0.3 0.2
        -0.1 0.3 0.1 -0.2
        0.3 -0.1 0.2 0.1
        -0.1 -0.1 -0.1 -0.1
    ]
    bamm = BaMM("precompile_bamm", bamm_weights, 0, 4)
    scorebounds(bamm)
    scan(bamm, seq_enc; strands=ForwardOnly())
    scan(bamm, batch; strands=BestStrand())

    # ── Cache fingerprint ───────────────────────────────────────────────
    model_fingerprint(pwm)
    return sequence_fingerprint(batch)
end

@compile_workload begin
    _precompile_workload()
end
