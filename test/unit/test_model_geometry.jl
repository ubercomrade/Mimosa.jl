# Characterization tests for the model geometry contract (ADR 0003).
#
# These tests pin the public geometry formulas for all five built-in model
# families before the extensibility API migration. They must keep passing
# unchanged through Stages 1-5: the migration is only allowed to *provide*
# these formulas, not to alter existing results.

using Test
using Mimosa

const EXAMPLES = joinpath(dirname(dirname(@__DIR__)), "examples")

# ── Geometry formula pinning ──────────────────────────────────────────────
#
# For every built-in model the following identities must hold:
#   window_size(model) == left_context(model) + motif_length(model) + right_context(model)
#   npositions(model, L) == max(L - window_size(model) + 1, 0)
#   site_start_offset(model) == left_context(model)
# Built-in models may have their own overrides of window_size / site_start_offset,
# but the identities above are the contract characterization.

function _check_geometry_identities(model::AbstractMotifModel, name::AbstractString)
    L = 200
    ml = motif_length(model)
    lc = left_context(model)
    rc = right_context(model)
    @test ml == motif_length(model)
    @test ml > 0
    @test lc >= 0
    @test rc >= 0
    @test window_size(model) == lc + ml + rc
    @test site_start_offset(model) == lc
    @test npositions(model, 0) == 0
    @test npositions(model, window_size(model) - 1) == 0
    @test npositions(model, window_size(model)) == 1
    @test npositions(model, L) == L - window_size(model) + 1
    # left/right defaults: model without overrides must still satisfy identities
    # (defaults are 0 for AbstractMotifModel)
    return nothing
end

@testset "ADR 0003: built-in geometry identity" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    bamm = readmodel(joinpath(EXAMPLES, "foxa2.ihbcp"))

    # Construct one model of each remaining family from minimal matrices.
    sitega = SiteGA("sitega", reshape(Float32.(1:100), 25, 4), 4)
    dimont = Dimont("dimont", reshape(Float32.(1:75), 25, 3), 1, 3)
    slim = Slim("slim", reshape(Float32.(1:75), 25, 3), 1, 3)

    _check_geometry_identities(pwm, "PWM")
    _check_geometry_identities(bamm, "BaMM")
    _check_geometry_identities(sitega, "SiteGA")
    _check_geometry_identities(dimont, "Dimont")
    _check_geometry_identities(slim, "Slim")
end

@testset "ADR 0003: built-in geometry per-family values" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    @test motif_length(pwm) == length(pwm)
    @test left_context(pwm) == 0
    @test right_context(pwm) == 0
    @test window_size(pwm) == motif_length(pwm)
    @test site_start_offset(pwm) == 0

    bamm = readmodel(joinpath(EXAMPLES, "foxa2.ihbcp"))
    @test motif_length(bamm) == bamm.motif_length
    @test left_context(bamm) == bamm.order
    @test right_context(bamm) == 0
    @test window_size(bamm) == bamm.motif_length + bamm.order
    @test site_start_offset(bamm) == bamm.order

    sitega = SiteGA("sitega", reshape(Float32.(1:100), 25, 4), 4)
    @test motif_length(sitega) == sitega.motif_length
    @test left_context(sitega) == 0
    @test right_context(sitega) == 0
    @test window_size(sitega) == sitega.motif_length
    @test site_start_offset(sitega) == 0

    dimont = Dimont("dimont", reshape(Float32.(1:75), 25, 3), 1, 3)
    @test motif_length(dimont) == dimont.motif_length
    @test left_context(dimont) == dimont.span
    @test right_context(dimont) == 0
    @test window_size(dimont) == dimont.motif_length + dimont.span
    @test site_start_offset(dimont) == dimont.span

    slim = Slim("slim", reshape(Float32.(1:75), 25, 3), 1, 3)
    @test motif_length(slim) == slim.motif_length
    @test left_context(slim) == slim.span
    @test right_context(slim) == 0
    @test window_size(slim) == slim.motif_length + slim.span
    @test site_start_offset(slim) == slim.span
end

# ── Scan-result invariants must remain stable ──────────────────────────────

@testset "ADR 0003: built-in scan result lengths unchanged" begin
    pwm = readmodel(joinpath(EXAMPLES, "pif4.meme"))
    bamm = readmodel(joinpath(EXAMPLES, "foxa2.ihbcp"))
    seq = UInt8[0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4]

    for model in (pwm, bamm)
        n_pos = npositions(model, length(seq))
        @test n_pos == max(length(seq) - window_size(model) + 1, 0)
        # Single-sequence scan returns a flat Vector{Float32} of length n_pos.
        forward = scan(model, seq; strands=ForwardOnly())
        reverse = scan(model, seq; strands=ReverseOnly())
        best = scan(model, seq; strands=BestStrand())
        both = scan(model, seq; strands=BothStrands())
        @test length(forward) == n_pos
        @test length(reverse) == n_pos
        @test length(best) == n_pos
        @test length(both.forward) == n_pos
        @test length(both.reverse) == n_pos
        @test best == max.(forward, reverse)
    end
end
