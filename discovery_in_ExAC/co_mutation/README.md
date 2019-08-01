# Statistical validation of significantly mutated genes

The notebook stat\_validation takes two lists of significantly mutated genes and two
mafs. To generate the mafs, first divide the total set of samples in half. The mutations
in the discovery set are then split into two sets, depending on which samples those mutations
are found in. Mutsig2CV is run independently on each of these two sets. The stat\_validation
notebook then considers the significantly mutated genes in one set, determines how often 
mutations are found in those genes in the second set of mutations, and then calculates
how many genes have mutation rates within binomial expectation set by the first set in
95% and 99% confidence intervals. We perform the same analysis for the set of significantly
mutated genes in the second set, and using the mutations in the first set

## Using the Cancer Gene Census (CGC) to suggest biological plausibility

cgc\_permutation performs a permutation test to identify how many times a set of randomly selected
genes equal to the number of significantly mutated genes has more genes in the CGC than
the current list of significantly mutated genes. The notebook also performs Fisher's exact
tests to determine if the set of significantly mutated genes is enriched for genes found in
the Cancer Gene Census v89
