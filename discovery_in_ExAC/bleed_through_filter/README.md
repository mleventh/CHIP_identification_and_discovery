# Bleed-Through filter

The filter identifies if mutations whose alternate alleles are in the trinucleotide reference
context occur at a frequency greater than empirical expectation of a set of 1000 samples
whose predominant mutational signature is C>T CpG aging. This includes the scripts found 
in the Terra method mleventh/bleed\_through_context, including code for a lego plotter.

In addition, experimental code for filtering out bleed through based on the clustering of
mutations whose alternate alleles are found in the reference trinucleotide context in a 
call stats file found in the Terra method mleventh/bleed-through-position is included, but
has not been thoroughly benchmarked
