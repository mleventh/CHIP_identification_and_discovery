# Blat realignment using simulated reads

The docker is for the Terra method mleventh/blat-reads-only, which uses the methodology employed
in the blat filtering step of the CGA pipeline. However, instead of using the reads from the 
sequencing data, simulated reads of 75 base pairs are generated from the hg19 reference sequence. 
The variant allele is tested in a sliding window of 5 base pairs with BLAT to identify to what other
sites in the genome the mutation maps. If mutations have higher scores with respect to mapping to 
other sites, then the variant is filtered out.

## Requirements
A 2-bit hg19 reference genome is necessary to run the tool, which can be found at 
gs://getzlab-workflows-reference_files-oa/hg19/hg19.2bit