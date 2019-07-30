source("https://bioconductor.org/biocLite.R")
#biocLite("GenomicRanges")
library(GenomicRanges)
library(data.table)

downsize<-fread("genovese_blacklists/phylogic_clustered_exacMuts_1kg_segdup_CNV_CHIP.txt",sep='\t',header=TRUE,stringsAsFactors=F)
downsize2<-downsize[is.na(downsize$Start_position)==FALSE,]
wig.iter<-fread("~/Documents/ExAC_filtering/germline_CNVs/kg_common_germline_CNVs.txt",header=T,stringsAsFactors=F)
colnames(wig.iter)[2:4]<-c("Chromosome","Start_position","End_position")
wig.iter$Chromosome<-gsub("chr","",wig.iter$Chromosome)

##CNV##

wig.iter<-wig.iter[,c(1:9,2514,2515)]
colnames(wig.iter)[1]<-"Chromosome"
wig.iter$Start_position<-as.numeric(wig.iter$Start_position)
wig.iter<-wig.iter[is.na(wig.iter$Start_position)==FALSE,]
wig.iter<-wig.iter[is.na(wig.iter$End_position)==FALSE,]

####

##for refseq and ucsc references##
wig.iter$Chromosome_refseq<-gsub("X","23",wig.iter$Chromosome_refseq)
wig.iter$Chromosome_refseq<-gsub("Y","24",wig.iter$Chromosome_refseq)
wig.iter<-wig.iter[is.na(wig.iter$Chromosome_refseq)==FALSE,]

#agilent<-fread("gaf_20111020_broad_wex_1.1_hg19.bed.interval_list",sep='\t',header=TRUE,stringsAsFactors=F)

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

#setkey(downsize2,Start_position)
#setkey(wig.iter,Start_position)
#scores<-downsize2[wig.iter]
scores<-merge(downsize.score,wig.iter,by.x="start",
              by.y="Start_position",all.x=TRUE)

scores.blacklisted<-scores[is.na(scores$POS),c(2:21)]
colnames(scores.blacklisted)<-gsub("\\.x","",colnames(scores.blacklisted))

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
  
###segdup###

scores<-scores[!duplicated(scores[,2:25])]
scores$has_segdup<-0
scores$has_segdup[is.na(scores$End_position.y)==FALSE]<-1
scores.blacklisted<-scores[scores$has_segdup==0,c(2:25)]


### dealing with CN overlaps ###

scores$gain_overlap<-0
scores$loss_overlap<-0

scores$gain_overlap[grep("CN3",scores$ALT)]<-1
scores$loss_overlap[grep("CN0",scores$ALT)]<-1

scores.blacklisted<-scores[scores$gain_overlap==0&scores$loss_overlap==0,c(2:25)]
write.table(scores.blacklisted,
            "all_clone/blacklist/exac_cluster_allClone_hwe_segdup_CNV.maf",
            sep='\t',quote=F,row.names=F)


### 1kg mask ###

scores<-scores[!duplicated(scores[,2:25])]
scores$mask<-0
scores$mask[is.na(scores$End_position.y)==TRUE]<-1
scores.blacklisted<-scores[scores$mask==0,c(2:25)]

