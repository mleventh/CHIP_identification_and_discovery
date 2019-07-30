library(data.table)
library(optparse)
options(stringsAsFactors=FALSE)

option.list <- list(
  make_option(c("-m","--maf"), action="store", default="maf.maf",type="character",help="specify maf here"),
  make_option(c("-s","--samplenames"),action="store",default="CHIP_",type="character",help="sample name"),
  make_option(c("-w","--whitelist"),action="store",default="integrative_chip_gene_list.txt",type="character",help="CHIP mutation whitelist")
)

opt  <- parse_args(OptionParser(option_list=option.list, usage = "Rscript %prog [options]"), print_help_and_exit=FALSE)

maf<-opt$maf
whitelist<-opt$whitelist
samplenames<-opt$samplenames

#read in maf file
non4<-fread(maf,sep='\t',header=T,stringsAsFactors=F)

#read in whitelist
chip.whitelist<-read.delim(whitelist,sep='\t',header=T,stringsAsFactors=F)

#identify missense mutation matches
miss.chip.protein.changes<-gsub("p\\.","",non4$Protein_Change)
whitelisted.chip.missense<-non4[which(miss.chip.protein.changes%in%chip.whitelist$Protein_Change&
                                        non4$Hugo_Symbol%in%chip.whitelist$Hugo_Symbol&non4$Variant_Classification=="Missense_Mutation"),]

#identify nonsense mutation changes except ASXL1, which is handled later
chip.nonsense<-chip.whitelist[which(substring(chip.whitelist$Protein_Change,12,19)=="nonsense"&
                                      chip.whitelist$Hugo_Symbol!="ASXL1"),]
whitelisted.chip.nonsense<-non4[which(non4$Hugo_Symbol%in%chip.nonsense$Hugo_Symbol&non4$Variant_Classification=="Nonsense_Mutation"),]

#identify splice site mtuations
chip.splice<-chip.whitelist[which(substring(chip.whitelist$Protein_Change,21,31)=="splice-site"),]
whitelisted.chip.splice<-non4[which(non4$Hugo_Symbol%in%chip.splice$Hugo_Symbol&non4$Variant_Classification=="Splice_Site"),]

#identify indels except ASXL1
chip.indel<-chip.whitelist[which(substring(chip.whitelist$Protein_Change,1,10)=="Frameshift"&
                                   chip.whitelist$Hugo_Symbol!="ASXL1"),]
whitelisted.chip.indel<-non4[which(non4$Hugo_Symbol%in%chip.indel$Hugo_Symbol&(substring(non4$Variant_Classification,1,11)=="Frame_Shift")),]

#call specific amino acid ranges and loci with low ExAC frequencies
chip.CBL<-non4[non4$Hugo_Symbol=="CBL"&nchar(miss.chip.protein.changes)==5&
                 (as.numeric(substring(miss.chip.protein.changes,2,4))>=381&
                    as.numeric(substring(miss.chip.protein.changes,2,4))<=421)&
                 non4$Variant_Classification=="Missense_Mutation",]

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

chip.PPM1D<-non4[non4$Hugo_Symbol=="PPM1D"&((as.numeric(substring(non4$Genome_Change,9,16))>=58733960|
                                               as.numeric(substring(non4$Genome_Change,9,16))<=58734202)|
                                              (as.numeric(substring(non4$Genome_Change,9,16))>=58740356|
                                                 as.numeric(substring(non4$Genome_Change,9,16))<=58740913))&
                                              (non4$Variant_Classification=="Nonsense_Mutation"|
                                              substring(non4$Variant_Classification,1,11)=="Frame_Shift")]

chip.CBLB<-non4[non4$Hugo_Symbol=="CBLB"&nchar(miss.chip.protein.changes)==5&
                 (as.numeric(substring(miss.chip.protein.changes,2,4))>=360&
                    as.numeric(substring(miss.chip.protein.changes,2,4))<=430)&
                  non4$Variant_Classification=="Missense_Mutation",]

chip.CEBPA<-non4[non4$Hugo_Symbol=="CEBPA"&non4$Protein_Change%in%c("p.R327W","p.G122E",
                                                                    "p.G223S","p.G242S","p.G102S"),]

chip.CSNK1A1<-non4[non4$Hugo_Symbol=="CSNK1A1"&non4$Protein_Change%in%c("p.T321A",
                                                                        "p.E82K","p.R282H"),]

chip.DDX41<-non4[non4$Hugo_Symbol=="DDX41"&non4$Protein_Change%in%c("p.P258L","p.G67R","p.R219C","p.Y33H",
                                                                    "p.P38L","p.H202Q","p.H274R"),]

chip.DNMT3A<-non4[non4$Hugo_Symbol=="DNMT3A"&nchar(miss.chip.protein.changes)==5&
                    ((as.numeric(substring(miss.chip.protein.changes,2,4))>=290&
                       as.numeric(substring(miss.chip.protein.changes,2,4))<=374)|
                    (as.numeric(substring(miss.chip.protein.changes,2,4))>=626&
                       as.numeric(substring(miss.chip.protein.changes,2,4))<=910))&
                    non4$Variant_Classification=="Missense_Mutation",]

chip.ETV6<-non4[non4$Hugo_Symbol=="ETV6"&((nchar(miss.chip.protein.changes)==5&
                  ((as.numeric(substring(miss.chip.protein.changes,2,4))>=338&
                      as.numeric(substring(miss.chip.protein.changes,2,4))<=424)))|
                     (nchar(miss.chip.protein.changes)==4&
                        (as.numeric(substring(miss.chip.protein.changes,2,4))>=58&
                        as.numeric(substring(miss.chip.protein.changes,2,4))<=99))|
                    (nchar(miss.chip.protein.changes)==5&
                       (as.numeric(substring(miss.chip.protein.changes,2,4))>=100&
                          as.numeric(substring(miss.chip.protein.changes,2,4))<=123)))&
                  non4$Variant_Classification=="Missense_Mutation",]

chip.ETV6.locus<-non4[non4$Hugo_Symbol=="ETV6"&non4$Protein_Change%in%c("p.R210H","p.R211H","p.R49L","p.S30L","p.P25L")]

chip.EZH2<-non4[non4$Hugo_Symbol=="EZH2"&nchar(non4$Protein_Change)==7&
                  (as.numeric(substring(miss.chip.protein.changes,2,4))>=500&
                     as.numeric(substring(miss.chip.protein.changes,2,4))<=738)&
                  non4$Variant_Classification=="Missense_Mutation",]

chip.GATA2<-non4[non4$Hugo_Symbol=="GATA2"&
                   substring(non4$Genome_Change,nchar(non4$Genome_Change)-6,
                             nchar(non4$Genome_Change)-3)>=2118&
                   substring(non4$Genome_Change,nchar(non4$Genome_Change)-6,
                             nchar(non4$Genome_Change)-3)<=2197,]

chip.NOTCH1.review<-non4[non4$Hugo_Symbol=="NOTCH1"&nchar(miss.chip.protein.changes)==5&
                    ((as.numeric(substring(miss.chip.protein.changes,2,4))>=617&
                        as.numeric(substring(miss.chip.protein.changes,2,4))<=738)|
                       (as.numeric(substring(miss.chip.protein.changes,2,4))>=500&
                          as.numeric(substring(miss.chip.protein.changes,2,4))<=616))&
                    non4$Genome_Change!="g.chr9:139409976C>T"&
                    non4$Variant_Classification=="Missense_Mutation",]

chip.NOTCH1<-non4[non4$Hugo_Symbol=="NOTCH1"&nchar(miss.chip.protein.changes)==6&
                    ((as.numeric(substring(miss.chip.protein.changes,2,4))>=2065&
                        as.numeric(substring(miss.chip.protein.changes,2,4))<=2555))
                    &non4$Variant_Classification%in%c("Nonsense_Mutation",
                                                      "Frame_Shift_Ins","Frame_Shift_Del"),]

chip.NOTCH2<-non4[non4$Hugo_Symbol=="NOTCH2"&nchar(miss.chip.protein.changes)==6&
                    ((as.numeric(substring(miss.chip.protein.changes,2,4))>=2010&
                        as.numeric(substring(miss.chip.protein.changes,2,4))<=2041))
                  &non4$Variant_Classification%in%c("Nonsense_Mutation",
                                                    "Frame_Shift_Ins","Frame_Shift_Del"),]

chip.NRAS<-non4[non4$Hugo_Symbol=="NRAS"&non4$Protein_Change%in%c("p.Q61L","p.P185S","p.A134T",
                                                                  "p.A91V"),]

chip.PIGA<-non4[non4$Hugo_Symbol=="PIGA"&(non4$Protein_Change%in%c("p.H11R","p.E476G","p.G474R","p.K293E",
                                                                   "p.H419Q","p.L56P")),]

chip.PTEN<-non4[non4$Hugo_Symbol=="PTEN"&non4$Protein_Change%in%c("p.N288S","p.I224M","p.N117S","p.G143D",
                                                                  "p.M270I"),]

chip.PTPN11<-non4[non4$Hugo_Symbol=="PTPN11"&(non4$Protein_Change%in%c("p.E97Q","p.P561L","p.S548P",
                                                                       "p.T566M","p.P424S")),]

chip.RPL11<-non4[non4$Hugo_Symbol=="RPL11"&non4$Protein_Change%in%c("p.A79V","p.R18H"),]

chip.RUNX1<-non4[non4$Hugo_Symbol=="RUNX1"&non4$Protein_Change%in%c("p.S410L","p.Q286K",
                                                                   "p.M240I","p.P395L",
                                                                   "p.G367D","p.G340D",
                                                                   "p.R223H","p.P395R","p.S12L",
                                                                   "p.P76R","p.R205L","p.H78Y","p.L29S","p.E53D"),]

chip.RPS7<-non4[non4$Hugo_Symbol=="RPS7"&non4$Protein_Change=="p.P104R",]

chip.SH2B3<-non4[non4$Hugo_Symbol=="SH2B3"&nchar(miss.chip.protein.changes)==5&
                    ((as.numeric(substring(miss.chip.protein.changes,2,4))>=364&
                        as.numeric(substring(miss.chip.protein.changes,2,4))<=441)|
                       (as.numeric(substring(miss.chip.protein.changes,2,4))>=195&
                          as.numeric(substring(miss.chip.protein.changes,2,4))<=307))&
                    non4$Variant_Classification=="Missense_Mutation",]

chip.SRSF2<-non4[non4$Hugo_Symbol=="SRSF2"&substring(non4$Protein_Change,1,5)=="p.P95",]

chip.EZH2<-non4[non4$Hugo_Symbol=="STAT5B"&nchar(miss.chip.protein.changes)==5&
                  (as.numeric(substring(miss.chip.protein.changes,2,4))>=593&
                     as.numeric(substring(miss.chip.protein.changes,2,4))<=670)&
                  non4$Variant_Classification=="Missense_Mutation",]

chip.TP53<-non4[non4$Hugo_Symbol=="TP53"&non4$Protein_Change%in%c("p.R110H","p.N131S",
                                                                  "p.G244A","p.M246L",
                                                                  "p.P8T","p.N131T",
                                                                  "p.L194H","p.G244A",
                                                                  "p.I50T","p.E287D","p.V203L",
                                                                  "p.Y107H","p.G112S",
                                                                  "p.V216L","p.D352V","p.R158G"),]

chip.VPS45<-non4[non4$Hugo_Symbol=="VPS45"&non4$Protein_Change%in%c("p.R433Q","p.K67R","p.V262A","p.I247V",
                                                                    "p.R386H","p.I62M","p.R176L","p.P253S","p.M63I"),]

chip.JAK3<-non4[non4$Hugo_Symbol=="JAK3"&non4$Protein_Change%in%c("p.R582Q","p.M558V","p.G1101R","p.A919V",
                                                                  "p.V619E","p.Q827E","p.A573V","p.F1123L",
                                                                  "p.E1072K","p.G642R","p.D47N","p.A853V"),]

chip.KRAS<-non4[non4$Hugo_Symbol=="KRAS"&non4$Protein_Change%in%c("p.G12R","p.G12C",
                                                                  "p.G13C",
                                                                  "p.G13D","p.G13R",
                                                                  "p.Q61R","p.Q61L",
                                                                  "p.Q22K","p.Q61H","p.Y64D",
                                                                  "p.A59E","p.L19F"),]

chip.GATA1<-non4[non4$Hugo_Symbol=="GATA1"&non4$Protein_Change%in%c("p.L386P","p.R317Q","p.V32I","p.L194P","p.L180V"),]

chip.PHF6<-non4[non4$Hugo_Symbol=="PHF6"&nchar(miss.chip.protein.changes)==5&
                  ((as.numeric(substring(miss.chip.protein.changes,2,3))>=42&
                      as.numeric(substring(miss.chip.protein.changes,2,4))<=132)|
                     (as.numeric(substring(miss.chip.protein.changes,2,4))>=239&
                        as.numeric(substring(miss.chip.protein.changes,2,4))<=330))&
                  non4$Variant_Classification=="Missense_Mutation"]

chip.SF3B1<-non4[non4$Hugo_Symbol=="SF3B1"&nchar(miss.chip.protein.changes)==5&
                   (as.numeric(substring(miss.chip.protein.changes,2,4))>=550&
                       as.numeric(substring(miss.chip.protein.changes,2,4))<=800)&
                      non4$Variant_Classification=="Missense_Mutation",]

chip.STAT5B<-non4[non4$Hugo_Symbol=="STAT5B"&non4$Protein_Change=="p.T628S",]

chip.TERT<-non4[non4$Hugo_Symbol=="TERT"&non4$Protein_Change%in%c("p.E441del","p.E439Q",
                                                                  "p.R962S","p.P173S",
                                                                  "p.R962C","p.G804V","p.R230Q",
                                                                  "p.A1117S","p.R646C","p.R743Q",
                                                                  "p.Q53H","p.T1111M","p.R756C","p.V1070M",
                                                                  "p.K1050N","p.L394Q","p.A745V","p.A789V"),]

chip.U2AF2<-non4[non4$Hugo_Symbol=="U2AF2"&non4$Protein_Change=="p.INQDK191del",]

#aggregate CHIP mutations
whitelisted.overall<-rbind(whitelisted.chip.missense,whitelisted.chip.nonsense,whitelisted.chip.splice,
                           whitelisted.chip.indel,chip.ASXL1,chip.CBL,chip.CBLB,chip.DNMT3A,
                           chip.ETV6,chip.EZH2,chip.GATA1,chip.GATA2,chip.JAK3,chip.KRAS,
                           chip.NOTCH1,chip.NRAS,chip.PIGA,chip.PTEN,chip.PTEN,chip.PTPN11,
                           chip.RPL11,chip.SH2B3,chip.SRSF2,chip.TET2,chip.TP53,chip.VPS45,
                           chip.ASXL1.locus,chip.DDX41,chip.RUNX1,chip.SF3B1,chip.STAT5B,chip.TERT,
                           chip.CEBPA,chip.PPM1D,chip.CSNK1A1,chip.ETV6.locus,chip.RPS7,chip.NOTCH2,
                           chip.U2AF2,fill=TRUE)

whitelisted.overall<-whitelisted.overall[!duplicated(whitelisted.overall[,1:17]),]

write.table(whitelisted.overall,paste0(samplenames,"_whitelisted.txt"),sep='\t',quote=F,row.names=F)
