# Phylogic 1d clustering

This directory contains the scripts for the Terra method mleventh/phylogic\_1d\_cluster,
which implements Dirichlet clustering as described in Leschiner et al. 2019 for an individual
sample to identify clusters of somatic mutations. The script formats oncotated MAF files
for use with Phylogic, runs the clustering in one dimension without a seg file, and then
outputs two files: one that includes all clones, and one that selects the largest clone
with at least four mutations. If running phylogic\_preprocess\_genome.py, the expected 
input is a VCF, but whole genome mutation calls from a maf can be run with 
phylogic\_preprocess.py


It is recommended that as much filtering of somatic mutations be performed before clustering
and before deciding on what largest clonal cluster to analyze.