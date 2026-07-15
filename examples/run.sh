#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd -- "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_DIR"

julia --project=. app/mimosa.jl \
  profile examples/myog.ihbcp examples/pif4.meme \
  --model1-type bamm --model2-type pwm \
  --fasta examples/foreground.fa --metric dice

julia --project=. app/mimosa.jl \
  profile examples/sitega.mat examples/pif4.meme \
  --model1-type sitega --model2-type pwm \
  --fasta examples/foreground.fa --metric co

julia --project=. app/mimosa.jl \
  profile examples/scores_1.fasta examples/scores_2.fasta \
  --model1-type scores --model2-type scores --metric cosine
