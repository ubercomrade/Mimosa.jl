$ErrorActionPreference = "Stop"
Set-StrictMode -Version Latest
Set-Location $PSScriptRoot

mimosa profile ./myog.ihbcp ./pif4.meme `
  --model1-type bamm `
  --model2-type pwm `
  --metric dice

mimosa profile ./gata2.ihbcp ./gata4.ihbcp `
  --model1-type bamm `
  --model2-type bamm `
  --metric co

mimosa profile ./sitega_stat6.mat ./pif4.meme `
  --model1-type sitega `
  --model2-type pwm `
  --metric co

mimosa motif ./sitega_gata2.mat ./pif4.meme `
  --model1-type sitega `
  --model2-type pwm `
  --metric ed

mimosa motif ./pif4.meme ./pif4.meme `
  --model1-type pwm `
  --model2-type pwm `
  --metric ed

mimosa motif ./sitega_stat6.mat ./pif4.meme `
  --model1-type sitega `
  --model2-type pwm `
  --metric pcc `
  --pfm-mode

mimosa profile ./sitega.mat ./pif4.meme `
  --model1-type sitega `
  --model2-type pwm `
  --metric co

mimosa motif ./sitega_stat6.mat ./sitega_gata2.mat `
  --model1-type sitega `
  --model2-type sitega `
  --metric pcc `
  --pfm-mode

mimosa motif ./sitega_stat6.mat ./sitega_gata2.mat `
  --model1-type sitega `
  --model2-type sitega `
  --metric ed `
  --pfm-mode

mimosa motif ./sitega_stat6.mat ./sitega_stat6.mat `
  --model1-type sitega `
  --model2-type sitega `
  --metric ed `
  --pfm-mode

mimosa motif ./gata2.meme ./sitega_gata2.mat `
  --model1-type pwm `
  --model2-type sitega `
  --metric ed `
  --pfm-mode

mimosa motif ./gata2.ihbcp ./gata4.ihbcp `
  --model1-type bamm `
  --model2-type bamm `
  --metric ed `
  -v

mimosa profile ./scores_1.fasta ./scores_2.fasta `
  --model1-type scores `
  --model2-type scores
