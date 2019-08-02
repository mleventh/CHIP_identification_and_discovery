# Indelocator filters

Indelocator, due to its lack of local realignment, miscounting of split reads and ignorance of homopolymers, 
needs additional filtering in order to get high confidence mutation calls. The jupyter notebook includes steps 
to eliminate indels found in long homopolymer runs, long indels, and indels found in sites of extraordinary high depth 
