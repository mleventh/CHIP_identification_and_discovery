# Filtering out rare germline mutations with SVM regression (SNPKebab)

This directory contains the notebooks used to prepare and train SVM regression using
20,000 germline SNPs from a subset of the ExAC consortium. 

The proximal_afs notebook details how mutation annotation files and call stats files were
annotated with information about the allele fractions of the nearest SNPs and the distribution
of SNP variant allele fractions at given sites. In addition, the notebook also annotates these
files for the distance to the nearest SNP and bait boundary. In addition to using this notebook,
MappabilityFilter.R uses a GenomicRanges implementation to identify sites that overlap with
common germline CNVs, segmental duplications, and the UCSC mappability tracks. In addition,
this script can be used to filter out sites that fail Hardy-Weinberg Equilibrium and sites
excluded by the 1000 Genomes strict mask.

## Requirements

Hardy-Weinberg failures: https://github.com/lh3/varcmp/blob/master/scripts/LCRhs37d5.bed.gz 
and https://github.com/lh3/varcmp/blob/master/scripts/1000g.hwe-bad.bed

Segmental duplications: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz

Strict mask of 1000 Genomes: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis\_results/supporting/accessible\_genome\_mask
s/20120824\_strict_mask.bed

Common germline CNVs: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated\_sv_map/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz

Mappability track: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign75mer.bigWig

Bait boundaries: gs://getzlab-workflows-reference\_files-oa/hg19/agilent/whole\_exome\_agilent\_1.1\_refseq\_plus\_3\_boosters.Homo\_sapiens\_assembly19.baits.interval_list

