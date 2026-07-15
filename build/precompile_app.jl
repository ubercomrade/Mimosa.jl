using Mimosa

const ROOT = normpath(joinpath(@__DIR__, ".."))
const EXAMPLES = joinpath(ROOT, "examples")

Mimosa.main(["--version"])
Mimosa.main(["inspect-model", joinpath(EXAMPLES, "pif4.meme"), "--type", "pwm"])
Mimosa.main([
    "profile",
    joinpath(EXAMPLES, "scores_1.fasta"),
    joinpath(EXAMPLES, "scores_2.fasta"),
    "--model1-type",
    "scores",
    "--model2-type",
    "scores",
    "--metric",
    "cosine",
])
