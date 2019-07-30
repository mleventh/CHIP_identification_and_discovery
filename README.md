# CHIP identification and discovery
Scripts to identify CHIP mutations in mutation annotation files and for identifying somatic mutations in normal blood samples for the identification of recurrently mutated genes in blood

# CHIP identification

Presented are CHIP whitelists and scripts that parse them. The lists come from Jaiswal et al. 2014, Jaiswal et al. 2017, Lindsley et al. 2017, and a list that integrates the genes from the three publications. All of the scripts run on oncotated MAFs except for rules_filter.py, which runs on varscan outputs from the mleventh/ebert-pipeline method in Terra.

# CHIP discovery

In this directory are scripts used for filtering and post-processing of somatic mutation calls in order to identify recurrent somatic mutations. For bleed_through_filter, blat_realignment_reads and phylogic_1d_clustering, these directories contain the dockerfiles for the Terra methods mleventh/bleed_through_context, mleventh/blat-reads-only, and mleventh/phylogic_1d_cluster respectively. 
