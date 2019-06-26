# The genomic landscape of divergence across the speciation continuum in an island-colonising bird

This GitHub repository contains custom scripts used in the manuscript: 'Sendell-Price et al. The genomic landscape of divergence across the speciation continuum in an island-colonising bird' submitted to Evolution. Custom scripts include:

### R script used to compare distributional skew (and kurtosis) of FST values:
1) CompareSkew.R  

### R scripts used to detect genomic island of divergence and genomic valleys of simularity; extract statistics from islands and valleys (mean FST, mean dxy, mean pi etc); get genes contained within each island/valley. 
1) Detect_Islands_Valleys_BiomartQuery.R - main script for island/valley detection.
2) Island_Detecting_Functions.R - script defining further functions required by Detect_Islands_Valleys_BiomartQuery.R.
note: island/valley detection scripts were originally developed by Benjamin Van Doren (benjamin.vandoren@zoo.ox.ac.uk) for the manuscript: Van Doren et al. (2017) Correlated patterns of genetic diversity and differentiation across an avian family.

### R script used to detect outlier SNPs (using the package PCAdapt)
1) Detect_outlier_SNPs_using_PCAdapt.R
note: this script requires use of pairwise VCF files (one file per population comparison).

### R script used to reorder VCF file into Zebra Finch chromosomes.
1) Reorder_VCF.R
note: in addition to Zosterops lateralis mapped VCF files, this script requires the following files: ZLat_scaffold_order_from_ZFinch.csv (output from satsuma synteny) and Zost_Scaffold_Lengths.txt 

### R scripts used to conduct individual-based models of population divergence:
1) simulations.R - main script for conducting simulations.
2) MainFunction.R - fuctions required 
note: simulation scripts were developed by Claudio Quilodrán (claudio.quilodran@zoo.ox.ac.uk) and Eric Anderson (eric.anderson@noaa.gov). For queries regarding simulating two-population divergence please contact either Claudio or Eric. 

### Data access:
Resequencing data from this study have been submitted to the National Center for Biotechnology Information (NCBI; https://www.ncbi.nlm.nih.gov) under accession number PRJNA489169. 

A Zebra Finch orientated VCF file, which formed the basis of our analyses, is available via Dropbox: "https://www.dropbox.com/sh/ivqj19l7qyrv5fm/AAB5_mEfV8njcf2AsyO_0aVVa?dl=0". Assignment of samples to populations are detailed in the file 'Samples_in_pops.txt'. The five comparisons used in this study are:

SI vs. CI = Early stage with gene flow, SI vs. FP = Early stage no gene flow, HI vs. ML = Mid stage with gene flow, LF vs. GT = Late stage with gene flow, LH vs. VN = Late stage no gene flow.
