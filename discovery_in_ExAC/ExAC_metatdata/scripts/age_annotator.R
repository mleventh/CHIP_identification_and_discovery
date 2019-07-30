library(data.table)
library(ggrepel)
library(weights)
library(optparse)

options(stringsAsFactors=FALSE)

option.list <- list(
  make_option(c("-m","--maf"), action="store", default="maf.maf",type="character",help="specify maf here"),
  make_option(c("-s","--samplenames"),action="store",default="CHIP_",type="character",help="sample name"),
  make_option(c("-w","--annotation"),action="store",default="age_annotation.txt",type="character",help="Age of samples in ExAC"),
  make_option(c("--known"),action="store",default="canon_sig_genes.txt",type="character",help="Canonical CHIP drivers"),
  make_option(c("--uncharacterized"),action="store",default="unchar_sig_genes.txt",type="character",help="Previously uncharacterized Sig Genes")
)

opt  <- parse_args(OptionParser(option_list=option.list, usage = "Rscript %prog [options]"), print_help_and_exit=FALSE)

maf<-opt$maf
annotation<-opt$annotation
samplenames<-opt$samplenames
known<-opt$known
uncharacterized<-opt$uncharacterized

#annotate the somatic maf with information from the metadata
downsample.maf<-fread(maf,sep='\t',header=T,stringsAsFactors=F)

binomial<-function(a, b, p) {binom.test(a, b, p, alternative=
                                          c("two.sided"), conf.level = 0.95)$conf.int}

all.age<-fread(annotation,sep='\t',header=T,
               stringsAsFactors=F)


#annotate with age and remove duplicates
exac.age<-merge(downsample.maf,all.age,by="Tumor_Sample_Barcode")
exac.age<-exac.age[!duplicated(exac.age[,1:17]),]

#### plot number of samples per age ####

unique.samples<-exac.age[!duplicated(exac.age$Tumor_Sample_Barcode),]

bins=c("20-29","30-39","40-49","50-59","60-69","70-79","80+")
age.df<-data.frame(age=bins,counts=0,total=0,stringsAsFactors=F)

#age.df$counts[age.df$age=="<20"]<-length(which(unique.samples$age<20))
age.df$counts[age.df$age=="20-29"]<-length(which(unique.samples$age>=20&
                                            unique.samples$age<30))
age.df$counts[age.df$age=="30-39"]<-length(which(unique.samples$age>=30&
                                                   unique.samples$age<40))
age.df$counts[age.df$age=="40-49"]<-length(which(unique.samples$age>=40&
                                                   unique.samples$age<50))
age.df$counts[age.df$age=="50-59"]<-length(which(unique.samples$age>=50&
                                                   unique.samples$age<60))
age.df$counts[age.df$age=="60-69"]<-length(which(unique.samples$age>=60&
                                                   unique.samples$age<70))
age.df$counts[age.df$age=="70-79"]<-length(which(unique.samples$age>=70&
                                                   unique.samples$age<80))
age.df$counts[age.df$age=="80+"]<-length(which(unique.samples$age>=80))
#age.df$counts[age.df$age=="90+"]<-length(which(unique.samples$age>=90))


#age.df$total[age.df$age=="<20"]<-length(which(all.age$age<20))
age.df$total[age.df$age=="20-29"]<-length(which(all.age$age>=20&
                                                  all.age$age<30))
age.df$total[age.df$age=="30-39"]<-length(which(all.age$age>=30&
                                                  all.age$age<40))
age.df$total[age.df$age=="40-49"]<-length(which(all.age$age>=40&
                                                  all.age$age<50))
age.df$total[age.df$age=="50-59"]<-length(which(all.age$age>=50&
                                                  all.age$age<60))
age.df$total[age.df$age=="60-69"]<-length(which(all.age$age>=60&
                                                  all.age$age<70))
age.df$total[age.df$age=="70-79"]<-length(which(all.age$age>=70&
                                                  all.age$age<80))
age.df$total[age.df$age=="80+"]<-length(which(all.age$age>=80))
#age.df$total[age.df$age=="90+"]<-length(which(all.age$age>=90))

age.df$ratio<-age.df$counts/age.df$total

age.df$sig_gene_count<-0

#subdivide based on whether to plot trends for known or uncharacterized sig genes 
#with respect to previous CHIP literature

canon<-read.delim(known,sep='\t',header=T,stringsAsFactors=F)
novel<-read.delim(uncharacterized,sep='\t',header=T,stringsAsFactors = F)

whitelisted.muts<-exac.age[which(exac.age$Hugo_Symbol%in%novel$gene),]
sig.samples<-whitelisted.muts[!duplicated(whitelisted.muts$id),]

age.df$sig_gene_count[age.df$age=="20-29"]<-length(which(sig.samples$age>=20&
                                                           sig.samples$age<30))
age.df$sig_gene_count[age.df$age=="30-39"]<-length(which(sig.samples$age>=30&
                                                           sig.samples$age<40))
age.df$sig_gene_count[age.df$age=="40-49"]<-length(which(sig.samples$age>=40&
                                                           sig.samples$age<50))
age.df$sig_gene_count[age.df$age=="50-59"]<-length(which(sig.samples$age>=50&
                                                           sig.samples$age<60))
age.df$sig_gene_count[age.df$age=="60-69"]<-length(which(sig.samples$age>=60&
                                                           sig.samples$age<70))
age.df$sig_gene_count[age.df$age=="70-79"]<-length(which(sig.samples$age>=70&
                                                           sig.samples$age<80))
age.df$sig_gene_count[age.df$age=="80+"]<-length(which(sig.samples$age>=80))

age.df$sig_gene_ratio<-age.df$sig_gene_count/age.df$total

count.CI<-mapply(binomial,age.df$count,age.df$total,
                age.df$ratio)

age.df$low_CI<-(count.CI[1,])
age.df$high_CI<-(count.CI[2,])

ggplot(data = age.df[age.df$age!="20-29",],aes(age,ratio)) +
  geom_point()+
  geom_errorbar(aes(ymin=low_CI, ymax=high_CI))+
  #ggtitle("Frequency of mutations in genes by age")+
  theme_classic()+
  theme(plot.title=element_text(lineheight=1.0,face="bold",size=15))+
  theme(axis.text.x=element_text(size=20,angle=90))+
  theme(axis.text.y=element_text(size=20))+
  theme(axis.title=element_text(size=16))+
  theme(legend.text=element_text(size=14))+
  theme(legend.title=element_text(size=16))+
  #facet_grid(rows = vars(Hugo_Symbol),scale="free_y")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x="age",y="percent samples with uncharacterized significantly mutated gene",legend="gene")
ggsave(paste0(samplename,".png"))
ggsave(paste0(samplename,".pdf"))

