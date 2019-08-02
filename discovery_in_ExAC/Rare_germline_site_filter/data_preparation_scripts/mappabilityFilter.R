library(GenomicRanges)
library(data.table)

library(optparse)

options(stringsAsFactors=FALSE)

option.list <- list(
  make_option(c("-m","--maf"), action="store", default="maf.maf",type="character",help="specify maf here"),
  make_option(c("-s","--samplenames"),action="store",default="CHIP_",type="character",help="sample name"),
  make_option(c("-w","--annotation"),action="store",default="age_annotation.txt",type="character",help="Annotations associated with genomic intervals"),
)

opt  <- parse_args(OptionParser(option_list=option.list, usage = "Rscript %prog [options]"), print_help_and_exit=FALSE)

maf<-opt$maf
annotation<-opt$annotation
samplenames<-opt$samplenames

#Note: the segdup file from UCSC does not come with a header, so set header=FALSE for that file
downsize2<-fread(maf,sep='\t',header=TRUE,stringsAsFactors=F)
wig.iter<-fread(annotation,header=T,stringsAsFactors=F)

#set headers for segdup file and get rid of chr contigs
colnames(wig.iter)[2:4]<-c("Chromosome","Start_position","End_position")
wig.iter$Chromosome<-gsub("chr","",wig.iter$Chromosome)

##CNV annotations##

wig.iter<-wig.iter[,c(1:9,2514,2515)]
colnames(wig.iter)[1]<-"Chromosome"
wig.iter$Start_position<-as.numeric(wig.iter$Start_position)
wig.iter<-wig.iter[is.na(wig.iter$Start_position)==FALSE,]
wig.iter<-wig.iter[is.na(wig.iter$End_position)==FALSE,]

####

downsize.range<-with(downsize2, GRanges(Chromosome, IRanges(start=Start_position,end=End_position)))
wig.range<-with(wig.iter, GRanges(Chromosome, IRanges(start=Start_position,end=End_position)))

match<-findOverlaps(downsize.range,wig.range)
ranges(downsize.range)[queryHits(match)]<-ranges(wig.range)[subjectHits(match)]

start<-gsub("[0-9,A-Z]+:(.*?)-.*", "\\1", downsize.range)
start<-gsub("[0-9,A-Z]+:(.*?)","\\1",start)
end<-gsub("[0-9,A-Z]+:.*?-(.*)", "\\1", downsize.range)
end<-gsub("[0-9,A-Z]+:(.*?)","\\1",end)

downsize.score<-downsize2

downsize.score$start<-start 
downsize.score$end<-end

#ensure vectors are of the same type
wig.iter$Start_position<-as.character(wig.iter$Start_position)

scores<-merge(downsize.score,wig.iter,by.x="start",
              by.y="Start_position",all.x=TRUE)

scores.blacklisted<-scores[is.na(scores$POS),c(2:21)]

write.table(scores.blacklisted[,c(2:25,29)],"all_clone/blacklist/exac_cluster_allClone_hwe_segdup_CNV_map.maf",sep='\t',quote=F,
            row.names=F)

score.order<-scores$score[order(scores$score)]
plot(seq(1,length(scores$score),by=1),score.order)
good.align<-scores[scores$score>0.15,]
good.align.freq<-as.data.frame(table(good.align$Hugo_Symbol))
good.align.freq<-good.align.freq[order(-good.align.freq$Freq),]
good.align.freq[1:40,]

varids.non3 <- as.data.frame(table(good.align$Genome_Change))

varids.non4 <- varids.non3[varids.non3[,2]>=6,]
varids.non5 <- as.vector(varids.non4[,1])

non6 <- good.align[!good.align$Genome_Change%in%varids.non5,]

#sig.genes<-read.delim("mutsig/exac_purity_sig_genes.txt.sig_genes.txt",sep='\t',header=T,stringsAsFactors = F)
#non.sig.genes<-sig.genes[which(sig.genes$q>=0.1),]
full.mut<-non6
prob <- function (x,y){
  pbinom(x,x+y,prob=0.48)
}

#cumulative binomial probability
x<- full.mut$t_alt_count
y<- full.mut$t_ref_count

full.mut$prob <- prob(x,y)

#get q values from fdr
full.mut$qval <- p.adjust(full.mut$prob, method="fdr")

valid.genes<-read.delim("Sidd_valid_genes.txt",sep='\t',header=T)
normsns <- full.mut[(full.mut$qval<0.05&full.mut$Hugo_Symbol%in%valid.genes$x),]
somatic.freq<-as.data.frame(table(normsns$Hugo_Symbol))
somatic.freq<-somatic.freq[order(-somatic.freq$Freq),]

somatic.freq[1:40,]
write.table(somatic.freq[1:40,],"whitelist_all_genes_fixed_depth_purity_filter_3_22_beta=2.txt",row.names=F,quote=F,sep='\t')

somatic.freq$Var1 <- factor(somatic.freq$Var1, levels = somatic.freq$Var1[order(-somatic.freq$Freq)])

ggplot(somatic.freq[1:40,],aes(x=Var1,y=Freq))+ geom_col()+
  ggtitle("Genes with frequent mutations in ExAC")+
  theme(plot.title=element_text(lineheight=1.0,face="bold",size=15))+
  theme(axis.text.x=element_text(size=10,angle=90))+
  theme(axis.text.y=element_text(size=20))+
  theme(axis.title=element_text(size=20))+
  theme(legend.text=element_text(size=12))+
  theme(legend.title=element_text(size=16))+
  xlab("gene")+
  ylab("frequency")
ggsave("whitelist_all_genes_old_filter_fixed_depth_purity_filter_3_22_beta=2.png")

###HWE###

scores$hwe_fail<-0
scores$hwe_fail[is.na(scores$End_position.y)==FALSE]<-1
scores.blacklisted<-scores[scores$hwe_fail==0,2:25]
colnames(scores.blacklisted)<-gsub("\\.x","",colnames(scores.blacklisted))
write.table(scores.blacklisted,
            paste0(samplenames,"_blacklisted.txt",sep="")
            sep='\t',quote=F,row.names=F)
  
###segdup###

scores<-scores[!duplicated(scores[,2:25])]
scores$has_segdup<-0
scores$has_segdup[is.na(scores$End_position.y)==FALSE]<-1
scores.blacklisted<-scores[scores$has_segdup==0,c(2:25)]
colnames(scores.blacklisted)<-gsub("\\.x","",colnames(scores.blacklisted))
write.table(scores.blacklisted,
            paste0(samplenames,"_blacklisted.txt",sep="")
            sep='\t',quote=F,row.names=F)

### dealing with CN overlaps ###

scores$gain_overlap<-0
scores$loss_overlap<-0

scores$gain_overlap[grep("CN3",scores$ALT)]<-1
scores$loss_overlap[grep("CN0",scores$ALT)]<-1

scores.blacklisted<-scores[scores$gain_overlap==0&scores$loss_overlap==0,c(2:25)]
colnames(scores.blacklisted)<-gsub("\\.x","",colnames(scores.blacklisted))
write.table(scores.blacklisted,
            paste0(samplenames,"_blacklisted.txt",sep="")
            sep='\t',quote=F,row.names=F)


### 1kg mask ###

scores<-scores[!duplicated(scores[,2:25])]
scores$mask<-0
scores$mask[is.na(scores$End_position.y)==TRUE]<-1
scores.blacklisted<-scores[scores$mask==0,c(2:25)]
colnames(scores.blacklisted)<-gsub("\\.x","",colnames(scores.blacklisted))
write.table(scores.blacklisted,
            paste0(samplenames,"_blacklisted.txt",sep="")
            sep='\t',quote=F,row.names=F)

