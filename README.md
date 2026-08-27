[![Build and Push Docker Image](https://github.com/bixBeta/sartools/actions/workflows/docker-build.yml/badge.svg)](https://github.com/bixBeta/sartools/actions/workflows/docker-build.yml)
[![Docker Image](https://img.shields.io/badge/ghcr.io-bixbeta%2Fsartools-blue?logo=docker)](https://github.com/bixBeta/sartools/pkgs/container/sartools)

`nextflow pull bixBeta/sartools -r g2 `




`nextflow run bixBeta/sartools -r g2 --help `

```
S A R  T O O L S     W O R K F L O W  -  @bixBeta
=======================================================================================================================================================================
Usage:
    nextflow run https://github.com/bixbeta/sartools -r main < args ... >

Args:
    * --id             : TREx Project ID 
    * --ref            : Base Level (Denominator for log2FC calcs, must be in the group column of targetFile)
    * --target         : targetFile.txt (tab delim file with label, files and group mandatory columns)
    * --genome         : Reference genome (GRCh38, GRCm38 etc.)
    * --quarto         : < default: k3 > (render params used in makeReport.sh)
    * --annots         : Annotations source (NCBI or < default: ENSEMBL >)
    * --trimmer        : Read pre-processing tool named in the report (< default: fastp >)
```
