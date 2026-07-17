$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest
Set-Location (Join-Path $PSScriptRoot "..")

julia -m Mimosa `
  profile examples/myog.ihbcp examples/pif4.meme `
  --model1-type bamm --model2-type pwm `
  --fasta examples/foreground.fa --metric dice

julia -m Mimosa `
  profile examples/sitega.mat examples/pif4.meme `
  --model1-type sitega --model2-type pwm `
  --fasta examples/foreground.fa --metric co

julia -m Mimosa `
  profile examples/scores_1.fasta examples/scores_2.fasta `
  --model1-type scores --model2-type scores --metric cosine
