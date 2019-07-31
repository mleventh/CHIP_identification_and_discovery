# Filtering out rare germline mutations with SVM regression (SNPKebab)

This directory contains the notebooks used to prepare and train SVM regression using
20,000 germline SNPs from a subset of the ExAC consortium to predict the expected allele 
fraction of a variant at a given site in the genome. This tool is motivated by the observation
that rare germline mutations, due to poor mappability, bait bias or common germline copy
number events can be observed at low allele fractions and cluster with those of somatic
variants. These features are intrinsic to the genomic region in which these variants are
found. SNPKebab predicts what the expected allele fraction of germline mutations are
at a given site to see if the observed allele fraction of a putative somatic mutation
is consistent with the observed sample's purity or with the site-specific germline SNPs

## Data preparation
The proximal_afs notebook details how mutation annotation files and call stats files were
annotated with information about the allele fractions of the nearest SNPs and the distribution
of SNP variant allele fractions at given sites. In addition, the notebook also annotates these
files for the distance to the nearest SNP and bait boundary. In addition to using this notebook,
MappabilityFilter.R uses a GenomicRanges implementation to identify sites that overlap with
common germline CNVs, segmental duplications, and the UCSC mappability tracks. In addition,
this script can be used to filter out sites that fail Hardy-Weinberg Equilibrium and sites
excluded by the 1000 Genomes strict mask.

Note that since we are most concerned with the allele fraction skews of heterozygous mutations,
our training set only includes variants with allele fractions of at most 0.6, which is higher
than any individual sample's germline skew (quantified by the "germline\_cluster" output in 
the Phylogic-1D clustering output maf). The mutations are selected from the call stats by 
including sites with ExAC frequency greater than 1% and are rejected for the failure reasons
"germline_risk", "alt\_allele\_in\_normal" and "normal\_lod". To prevent class imbalance,
we make sure to include an equal number of mutations per allele fraction decile.

### Requirements

Hardy-Weinberg failures: https://github.com/lh3/varcmp/blob/master/scripts/LCRhs37d5.bed.gz 
and https://github.com/lh3/varcmp/blob/master/scripts/1000g.hwe-bad.bed

Segmental duplications: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz

Strict mask of 1000 Genomes: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase1/analysis\_results/supporting/accessible\_genome\_mask
s/20120824\_strict_mask.bed

Common germline CNVs: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated\_sv_map/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz

Mappability track: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign75mer.bigWig

Bait boundaries: gs://getzlab-workflows-reference\_files-oa/hg19/agilent/whole\_exome\_agilent\_1.1\_refseq\_plus\_3\_boosters.Homo\_sapiens\_assembly19.baits.interval_list

## Making predictions with SNPKebab

The SNPKebab notebook contains scripts for k-fold validation of the SVM regression model and 
applications to the test data using both a pre-determined training set ("small\_case\_shuffled\_medians\_segdup\_CNV\_proxy\_map\_dist\_bait\_gc.txt") and a pre-trained model
("svr\_germline.pickle"). For those interested, the notebook also includes attempts to train
a fully-connected Neural Net, whose performance may improve with the implementation of Kullback-
Leibler divergence

## Using SNPKebab predictions to filter mutations out

The Phylo_merger notebook details steps needed to filter based on the SNPKebab predictions. 
If the purity estimates come from phylogic clustering, one needs to merge with the pre-phylogic 
clustering maf in order to get the "t\_alt\_count" and "t\_ref\_count" annotations. Others can be added as well:
the "SwissProt\_acc\_Id","SwissProt\_entry\_Id" are highlighted because they are necessary for
Mutsig2CV stick figure plots. 

After the annotation merge, a likelihood ratio is calculated between two beta probability density 
functions shaped by the reference and alternate allele counts, one whose probability is the 
purity and another whose probability is the SNPKebab-predicted germline probability for that
variant. If the likelihood ratio is less than 0.9, the mutaiton is filtered out

Once the sites are filtered out, clone sizes are identified, filtering out clones with
fewer than 5 mutations and keeping the largest clone per sample. After this selection,
we test if, given the median depth of the mutations in the sample, the sample's germline
skew and the sample's purity the sample's true positive rate is above 0.9 and the 
sample's false discovery rate is under 0.1, keeping samples that fit this criteria

### Testing

Files test_annotated.txt and test.txt are provided for users to run proximal_afs->SNPKebab->Phylo_merger on their own
