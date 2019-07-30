library(data.table)
library(optparse)
options(stringsAsFactors=FALSE)

### Identify genes in the mutsig output that have CPM>=2 in hematopoietic stem cells (HSCs) 
### of myeloid lineage according to the Buenrostro et al. 2016 (GSE74246-GPL11154)

option.list <- list(
  make_option(c("--genelist"), action="store", default="sig_genes.txt",type="character",help="mutsig sig genes list"),
  make_option(c("-s","--samplenames"),action="store",default="",type="character",help="sample name"),
  make_option(c("--expresslist"),action="store",default="RNA_integrated_targets.txt",type="character",help="genes with expression in HSCs")
)

opt  <- parse_args(OptionParser(option_list=option.list, usage = "Rscript %prog [options]"), print_help_and_exit=FALSE)

genelist<-opt$genelist
expresslist<-opt$expresslist
samplenames<-opt$samplenames

key.combine.target.noDup<-fread(expresslist,sep='\t',
                                header=T,stringsAsFactors = F)
sig.gene.review<-fread(genelist,sep='\t',
                       header=T,stringsAsFactors=F)

sig.gene.review<-sig.gene.review[sig.gene.review$gene%in%key.combine.target.noDup$gene,]

sig.gene.review$q<-p.adjust(sig.gene.review$p,method='fdr')

sig.gene.review<-sig.gene.review[sig.gene.review$q<0.1,]

write.table(sig.gene.review,paste0(samplenames,"sig_gene_RNA.txt"),sep='\t',quote=F,
            row.names=F)