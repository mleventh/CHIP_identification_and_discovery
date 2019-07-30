library(data.table)
library(optparse)
options(stringsAsFactors=FALSE)

option.list <- list(
  make_option(c("-m","--maf"), action="store", default="maf.maf",type="character",help="specify maf here"),
  make_option(c("-s","--samplenames"),action="store",default="CHIP_",type="character",help="sample name"),
  make_option(c("-w","--whitelist"),action="store",default="chip_gene_list_Sidd_Coleman_combined.txt",type="character",help="CHIP mutation whitelist")
)

opt  <- parse_args(OptionParser(option_list=option.list, usage = "Rscript %prog [options]"), print_help_and_exit=FALSE)

maf<-opt$maf
whitelist<-opt$whitelist
samplenames<-opt$samplenames

#maf to whitelist
non4<-fread(maf,sep='\t',header=T,stringsAsFactors=F)

#whitelist
chip.whitelist<-read.delim(whitelist,sep='\t',header=T,stringsAsFactors=F)

#truncate protein change names to eliminate the 'p.' at the start 
miss.chip.protein.changes<-gsub("p\\.","",non4$Protein_Change)
whitelisted.chip.missense<-non4[which(miss.chip.protein.changes%in%chip.whitelist$Protein_Change&
                                        non4$Hugo_Symbol%in%chip.whitelist$Hugo_Symbol&non4$Variant_Classification=="Missense_Mutation"),]

#nonsense mutations but save ASXL1 for its own specific subset
chip.nonsense<-chip.whitelist[which(substring(chip.whitelist$Protein_Change,12,19)=="nonsense"&
                                      chip.whitelist$Hugo_Symbol!="ASXL1"),]
whitelisted.chip.nonsense<-non4[which(non4$Hugo_Symbol%in%chip.nonsense$Hugo_Symbol&non4$Variant_Classification=="Nonsense_Mutation"),]

chip.splice<-chip.whitelist[which(substring(chip.whitelist$Protein_Change,21,31)=="splice-site"),]
whitelisted.chip.splice<-non4[which(non4$Hugo_Symbol%in%chip.splice$Hugo_Symbol&non4$Variant_Classification=="Splice_Site"),]

#save ASXL1 as a special case
chip.indel<-chip.whitelist[which(substring(chip.whitelist$Protein_Change,1,10)=="Frameshift"&
                                   chip.whitelist$Hugo_Symbol!="ASXL1"),]
whitelisted.chip.indel<-non4[which(non4$Hugo_Symbol%in%chip.indel$Hugo_Symbol&(substring(non4$Variant_Classification,1,11)=="Frame_Shift")),]

chip.CBL<-non4[non4$Hugo_Symbol=="CBL"&nchar(miss.chip.protein.changes)==5&
                 (as.numeric(substring(miss.chip.protein.changes,2,4))>=381&
                    as.numeric(substring(miss.chip.protein.changes,2,4))<=421)&
                 non4$Variant_Classification=="Missense_Mutation",]

#special cases
chip.ASXL1<-non4[non4$Hugo_Symbol=="ASXL1"&(((as.numeric(substring(non4$Genome_Change,9,16))>=31021087&
                                                as.numeric(substring(non4$Genome_Change,9,16))<=31021720)|
                                               (as.numeric(substring(non4$Genome_Change,9,16))>=31021087&
                                                  as.numeric(substring(non4$Genome_Change,18,25))<=31021720))|
                                              ((as.numeric(substring(non4$Genome_Change,9,16))>=31022235&
                                                  as.numeric(substring(non4$Genome_Change,9,16))<=31025141)|
                                                 (as.numeric(substring(non4$Genome_Change,18,25))<=31025141)
                                               &as.numeric(substring(non4$Genome_Change,18,25))<=31025141)
)&
  (non4$Variant_Classification%in%c("Nonsense_Mutation","Frame_Shift_Ins","Frame_Shift_Del")),]

chip.ASXL1.locus<-non4[non4$Hugo_Symbol=="ASXL1"&nchar(miss.chip.protein.changes)>=5&
                         (as.numeric(substring(miss.chip.protein.changes,2,5))>=400&
                            as.numeric(substring(miss.chip.protein.changes,2,5))<=1540)&
                         (non4$Variant_Classification%in%c("Nonsense_Mutation","Frame_Shift_Ins","Frame_Shift_Del")),]

chip.TET2<-non4[non4$Hugo_Symbol=="TET2"&nchar(miss.chip.protein.changes)==6&
                  ((as.numeric(substring(miss.chip.protein.changes,2,5))>=1104&
                      as.numeric(substring(miss.chip.protein.changes,2,5))<=1481)|
                     (as.numeric(substring(miss.chip.protein.changes,2,5))>=1843&
                        as.numeric(substring(miss.chip.protein.changes,2,5))<=2002))&
                  (non4$Variant_Classification=="Missense_Mutation"|
                     non4$Variant_Classification=="In_Frame_Del"|
                     non4$Variant_Classification=="In_Frame_Ins"),]

chip.CBLB<-non4[non4$Hugo_Symbol=="CBLB"&nchar(miss.chip.protein.changes)==5&
                  (as.numeric(substring(miss.chip.protein.changes,2,4))>=360&
                     as.numeric(substring(miss.chip.protein.changes,2,4))<=430)&
                  non4$Variant_Classification=="Missense_Mutation",]

#aggregation
whitelisted.overall<-rbind(whitelisted.chip.missense,whitelisted.chip.nonsense,whitelisted.chip.splice,
                           whitelisted.chip.indel,chip.ASXL1,chip.CBL,chip.CBLB,chip.TET2)

#elimination of dups
whitelisted.overall<-whitelisted.overall[!duplicated(whitelisted.overall[,1:17]),]

write.table(whitelisted.overall,paste0(samplenames,"_whitelisted.txt"),sep='\t',quote=F,row.names=F)




