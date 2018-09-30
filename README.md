# The genomic landscape of divergence across the speciation continuum in an island-colonising bird

This GitHub repository contains custom scripts used in the manuscript 'The genomic landscape of divergence across the speciation continuum in an island-colonising bird' submitted to Molecular Ecology. Custom scripts include:

### R script used to compare distributional skew (and kurtosis) of FST values:
1) CompareSkew.R  

### R scripts used to detect genomic island of divergence and genomic valleys of simularity:
1) Detect_Islands_Valleys.R - main script for island/valley detection.
2) Island_Detecting_Functions.R - fuctions required 

note: island/valley detection scripts were developed by Benjamin Van Doren (benjamin.vandoren@zoo.ox.ac.uk).

### R scripts used to conduct individual-based models of population divergence:
1) simulations.R - main script for conducting simulations.
2) MainFunction.R - fuctions required 

note: simulation scripts were developed by Claudio Quilodrán (claudio.quilodran@zoo.ox.ac.uk) and Eric Anderson (eric.anderson@noaa.gov). For queries regarding simulating two-population divergence please contact either Claudio or Eric. 

### Data access:
Resequencing data from this study have been submitted to the National Center for Biotechnology Information (NCBI; https://www.ncbi.nlm.nih.gov) under accession number PRJNA489169. 

A Zebra Finch orientated VCF file, which formed the basis for our analyses, is available via Dropbox: "ADD LINK". Assignment of samples to populations are detailed in the file 'Samples_in_pops.txt'. The five comparisons used in this study are:

SI vs. CI = Early stage with gene flow
SI vs. FP = Early stage no gene flow
HI vs. ML = Mid stage with gene flow
LF vs. GT = Late stage with gene flow
LH vs. VN = Late stage no gene flow
