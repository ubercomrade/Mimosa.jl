#!/usr/bin/env julia
# Mimosa CLI entry point for use as a standalone script or package app.
#
# Usage:
#   julia --project=. app/mimosa.jl profile scores1.fasta scores2.fasta --model1-type scores --model2-type scores
#   julia --project=. app/mimosa.jl inspect-model model.ihbcp --type bamm
#   julia --project=. app/mimosa.jl --help

using Mimosa

exit(Mimosa.main(ARGS))
