## Nicolas Stransky
## The Broad Institute of MIT and Harvard / Cancer program.
## stransky@broadinstitute.org

## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.

## You should have received a copy of the GNU Lesser General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

suppressPackageStartupMessages(require(optparse))

option_list <- list(                                                                                                                         
		make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
				help="Print more output [default %default]"),
		make_option(c("--full.verbose"), action="store_true", default=FALSE,
				help="Print even more output [default %default]"),
		make_option(c("-o", "--outputdir"), action="store", type="character", default=".",
				help="Output directory", metavar="dir.path"),
		make_option(c("--firehose.mode"), action="store_true", default=FALSE,
				help="Run coMut as a FH task [default FALSE]"),
		make_option(c("-f", "--firehose.output"), action="store", type="character",
				help="Firehose output directory [default %default]", default="/local/cga-fh/cga/", metavar="dir.path"),
		make_option(c("-w", "--firehose.workspace"), action="store", type="character", default=NULL,
				help="Firehose workspace"),
		make_option(c("-r", "--firehose.mutsigrun"), action="store", type="character",
				help="MutSig version. e.g. mutsig1.0 mutsig1.5 etc. [default %default]", default="mutsig1.5"),
		make_option(c("--firehose.mutsig.mutcategs"), action="store", type="character", default=NULL,
				help="Full path to the mutcategs.txt file from a MutSig run. Only required in firehose.mode. [default FALSE]"),
		make_option(c("-p", "--mutsig.full.path"), action="store", type="character",
				help="Full path of a MutSig output directory (overrides automated guessing of Firehose directory)", default=NULL, metavar="dir.path"),
		make_option(c("-a", "--analysis.set"), action="store", type="character",
				help="Name of the analysis sample set"),
		make_option(c("-q", "--qthresh"), action="store", type="numeric", default=.1,
				help="q-value threshold of genes to display [default %default]", metavar="threshold"),
		make_option(c("-b", "--blacklist.file"), action="store", type="character",
				help="List of genes to ignore", default=NULL, metavar="file.path"),
		make_option(c("-l", "--whitelist.file"), action="store", type="character",
				help="Force the list of genes to appear on the plot", default=NULL, metavar="file.path"),
		make_option(c("-s", "--significance.file"), action="store", type="character",
				help="Provide the list of significant genes and q-values. Tab-delimited genes and q-values, no header. (overrides MutSig sig_genes)", 
				default=NULL, metavar="file.path"),
		make_option(c("-c", "--coverage.file"), action="store", type="character",
				help="Provide coverage for all samples. Tab-delimited sample and number of bases covered, no header. (overrides MutSig patients.counts_and_rates)", 
				default=NULL, metavar="file.path"),
		make_option(c("-m", "--maf"), action="store", type="character",
				help="MAF file containing all the required fields. (overrides MutSig maf)", 
				default=NULL, metavar="file.path", dest="maf.file"),
		make_option(c("--reference.genome"), action="store", type="character",
				help="Path to an alternate reference genome [default %default]", 
				default="/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta", metavar="file.path"),
		make_option(c("--no.allelic.fraction.boxplot"), action="store_false", default=TRUE,
				help="Suppress allelic fraction boxplots [default FALSE]", dest="allelic.fraction.boxplot"),
		make_option(c("--no.mutation.spectrum.plot"), action="store_false", default=TRUE,
				help="Suppress mutation spectrum plots [default FALSE]", dest="mutation.spectrum.plot"),
		make_option(c("--sort.by.mutation.status"), action="store_true", default=FALSE,
				help="Sort by mutation status of each gene instead of mutation frequency of each sample [default %default]"),
		make_option(c("--sort.genes.by.prevalence"), action="store_true", default=FALSE,
				help="Sort genes by gene mutation prevalence instead of q-value [default %default]"),
		make_option(c("--png"), action="store_true", default=FALSE,
				help="Output a PNG image (rasterized format) instead of a PDF (vectorial, editable) [default %default]")
)

opt <- parse_args(OptionParser(option_list=option_list, usage = "Rscript %prog [options]"), print_help_and_exit=FALSE)

if (opt[["help"]]) {
	print_help(OptionParser(option_list=option_list, usage = "Rscript %prog [options]"))
	cat("There are three modes to run coMut.R:
					- MutSig mode: coMut needs to know where MutSig results are located in the Firehose directory structure. This is if you want to run coMut on the 
					command line, to visualize MutSig results. To do so, the FH workspace (-w) and the FH analysis set (-a) are needed.
					- Firehose mode: Like the above, but don't guess anything, use all paths provided.
					- Standalone mode: coMut does not assume that some results have been precomputed in MutSig and recalculates them. 
					For that mode, you need to tell coMut where the maf file (-m), the coverage file (-c) and the list of q-values (-s) files are located, 
					as well as the name of your analysis set (-a, can be anything). If you don't provide a workspace or the complete path to the MutSig output, 
					it is assumed that you are running the standalone mode.
					
					Examples
					- Firehose mode
					Rscript coMut.R -v -o <output_dir> --firehose.mode -a <analysis_set> -s sig_genes.txt -c patients.counts_and_rates.txt -m final_analysis_set.maf --firehose.mutsig.mutcategs mutcategs.txt
					
					- MutSig mode
					Rscript coMut.R -v -o <output_dir> -a <analysis_set> -w <workspace>
					
					- Standalone mode
					Rscript coMut.R -v -o <output_dir> -a <analysis_set> -s sig_genes.txt -c patients.counts_and_rates.txt -m final_analysis_set.maf
					
					")
	quit(status = 1)
}

## Options parsing
if (TRUE) {
	mutsigrun <- "mutsig1.5"
	verbose <- FALSE
	full.verbose <- FALSE
	outputdir <- "coMut"
	firehose.output <- ""
	firehose.mode <- FALSE
	mutsig.full.path <- NULL
	workspace <- NULL
	qthresh <- 0.1
	allelic.fraction.boxplot <- FALSE
	mutation.spectrum.plot <- FALSE
	sort.by.mutation.status <- TRUE
	sort.by.disease.status <- FALSE
	sort.genes.by.prevalence <- FALSE
	reference.genome <- "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta"
	png.output <- FALSE
	analysis.set <- "coMut"
	significance.file <- "sig_genes.txt"
	maf.file <- "final_analysis_set.maf"
	coverage.file <- "patient_counts_and_rates.txt"
	blacklist.file <- ""
	whitelist.file <- ""
	firehose.mutsig.mutcategs <- ""
}

min.reads.to.count.allelic.fraction <- 20
max.mutation.rate <- 100

# verbose <- opt$verbose
# full.verbose <- opt$full.verbose
# outputdir <- opt$outputdir
# firehose.output <- opt$firehose.output
# workspace <- opt$firehose.workspace
# mutsigrun <- opt$firehose.mutsigrun
# firehose.mode <- opt$firehose.mode
# firehose.mutsig.mutcategs <- opt$firehose.mutsig.mutcategs
# analysis.set <- opt$analysis.set
# qthresh <- opt$qthresh
# maf.file <- ifelse(is.null(opt$maf.file), "", opt$maf.file)
# blacklist.file <- ifelse(is.null(opt$blacklist.file), "", opt$blacklist.file)
# whitelist.file <- ifelse(is.null(opt$whitelist.file), "", opt$whitelist.file)
# significance.file <- ifelse(is.null(opt$significance.file), "", opt$significance.file)
# coverage.file <- ifelse(is.null(opt$coverage.file), "", opt$coverage.file)
# mutsig.full.path <- opt$mutsig.full.path
# reference.genome <- opt$reference.genome
# allelic.fraction.boxplot <- opt$allelic.fraction.boxplot
# mutation.spectrum.plot <- opt$mutation.spectrum.plot
# sort.by.mutation.status <- opt$sort.by.mutation.status
# sort.genes.by.prevalence <- opt$sort.genes.by.prevalence
# png.output <- opt$png

if(full.verbose) verbose <- TRUE

if (verbose) cat("coMut v2.2 -- co-ocurrence of Mutations\n")

suppressPackageStartupMessages(require(RColorBrewer))
suppressPackageStartupMessages(require(gplots))

## Check arguments
if (is.null(analysis.set)) stop("analysis.set is missing")
if (firehose.mode & is.null(firehose.mutsig.mutcategs)) stop("firehose.mutsig.mutcategs is required with --firehose.mode")

if (!is.null(mutsig.full.path) | !is.null(workspace)) {
	if (verbose) cat("MutSig mode.\n")
	if (!is.null(workspace)) {  
		mutsig.output <- paste(firehose.output,"/", workspace, "/Individual_Set/", analysis.set, "/jobs/", mutsigrun, "/", sep="")
	} else {
		mutsig.output <- mutsig.full.path
	}
	if (!file.exists(mutsig.output)) stop(paste(mutsig.output, "is not readable"))
} else {
	mutsig.output <- ""
	if (!file.exists(significance.file) | !file.exists(coverage.file) | !file.exists(maf.file)) {stop("Check arguments. significance.file, coverage.file and maf.file are required.")}
}

if (file.exists(blacklist.file)) {
	if (verbose) cat("Reading blacklist")
	blacklist <- read.table(blacklist.file)[,1]
	if (verbose) cat(".\n")
} else {
	if (blacklist.file!="") warning("blacklist not found!")
	blacklist <- c()
}

if (file.exists(whitelist.file)) {
	if (verbose) cat("Reading whitelist")
	forceSelectedGenes <- read.table(whitelist.file)[,1]
	if (verbose) cat(".\n")
} else {
	if (whitelist.file!="") warning("whitelist not found!")
}

if (verbose) cat("Reading significance values")
if (file.exists(significance.file)) {
	if(any(c("gene", "name", "q", "q-value", "qvalue")%in%read.delim(significance.file, nrow=1, header=FALSE, colClasses="character")[1,])) {
		significant.gene.list <- read.delim(significance.file, header=TRUE, quote="")
	} else {
		significant.gene.list <- read.delim(significance.file, header=FALSE, quote="")
		colnames(significant.gene.list) <- c("gene", "q")
	}
} else {
	if (significance.file!="") stop("significance list not found!")
	significant.gene.list <- read.delim(paste(mutsig.output, "/", analysis.set,".sig_genes.txt", sep=""), header=TRUE, quote="")
}
if (verbose) cat(".\n")

if (verbose) cat("Reading coverage values")
if (file.exists(coverage.file)) {
	if(all(c("name", "N_tot")%in%read.delim(coverage.file, nrow=1, header=FALSE, colClasses="character")[1,])) {
		patients.counts_and_rates <- read.delim(coverage.file, header=TRUE, quote="")
	} else {
		patients.counts_and_rates <- read.delim(coverage.file, header=FALSE, quote="")
		colnames(patients.counts_and_rates) <- c("name", "N_tot")
	}
} else {
	if (coverage.file!="") stop("coverage values not found!")
	patients.counts_and_rates <- read.delim(paste(mutsig.output, "/", analysis.set,".patients.counts_and_rates.txt", sep=""), header=TRUE, quote="")
}
rownames(patients.counts_and_rates) <- patients.counts_and_rates$name
if (verbose) cat(".", nlevels(patients.counts_and_rates$name), "samples.\n")

if (verbose) cat("Reading maf file")
if (file.exists(maf.file)) {
	final_analysis_set <- read.delim(maf.file, header=TRUE, quote="", comment.char = "#")
} else {
	if (maf.file!="") stop("MAF file not found!")
	final_analysis_set <- read.delim(paste(mutsig.output, "/", analysis.set,".final_analysis_set.maf", sep=""), header=TRUE, quote="")
}


if (verbose) cat(".", nrow(final_analysis_set), "mutations in", nlevels(final_analysis_set$Tumor_Sample_Barcode), "samples.")
tmp <- lapply(list(c("_Mutation", ""), c("Silent", "Synonymous"), c("Splice_Site.*", "Splice_Site"), 
				c("Nonstop|De_novo_Start.*|Start_.*|Translation_.*|Read\\-through.*", "Other_non_syn."), c("In_frame.*", "In_frame_Indel"), c("Frame_Shift.*", "Frame_Shift"),
				c("3'\\-?UTR|5'\\-?UTR|3'\\-?Flank|5'\\-?Flank|IGR|Intron", "Silent")), 
		function(x) final_analysis_set$Variant_Classification <<- sub(x[1], x[2], as.character(final_analysis_set$Variant_Classification), ignore.case=TRUE))

## HERE modification for the new MAF from MutSig_2CV 
## :final_analysis_set$Variant_Classification now has 'RNA' 
final_analysis_set$Variant_Classification_Num <- as.numeric(factor(final_analysis_set$Variant_Classification, 
				levels=c("Synonymous", "In_frame_Indel", "Other_non_syn.", "Missense", "Splice_Site", "Frame_Shift", "Nonsense", "miRNA", "Silent", "RNA")))
final_analysis_set$Chromosome <- factor(sub("^M$", "MT", final_analysis_set$Chromosome))
if (verbose) cat("\n")

## Annotate MAF with the mutation categories
if (mutation.spectrum.plot && is.null(final_analysis_set$categ) == FALSE) {
	if (mutsig.output!="" & !firehose.mode) {
		mutcategs <- read.delim(paste(mutsig.output, "/", analysis.set, ".mutcategs.txt", sep=""))
	} else if (firehose.mode) {
		mutcategs <- read.delim(firehose.mutsig.mutcategs)
	} else {
		if (verbose) cat("Fetching genome reference")
		suppressPackageStartupMessages(require(Rsamtools))
		final_analysis_set$seq_context <- as.character(scanFa(reference.genome, GRanges(seqnames=final_analysis_set[, "Chromosome"], ranges=IRanges(start=final_analysis_set[, "Start_position"]-1, end=final_analysis_set[, "End_position"]+1), strand=final_analysis_set[, "Strand"])))
		if (verbose) cat(".\nDefining mutation categories")
		final_analysis_set$categ <- as.numeric(NA)
		final_analysis_set$categ[final_analysis_set$Variant_Type=="SNP" & grepl("CG$", final_analysis_set$seq_context) & final_analysis_set$Tumor_Seq_Allele1%in%c("T")] <- 4 # CpG transition
		final_analysis_set$categ[final_analysis_set$Variant_Type=="SNP" & is.na(final_analysis_set$categ) & grepl("^CG", final_analysis_set$seq_context) & final_analysis_set$Tumor_Seq_Allele1%in%c("A")] <- 4 # CpG transition
		final_analysis_set$categ[final_analysis_set$Variant_Type=="SNP" & is.na(final_analysis_set$categ) & grepl("CG$", final_analysis_set$seq_context) & final_analysis_set$Tumor_Seq_Allele1%in%c("G", "A")] <- 3 # CpG transversion
		final_analysis_set$categ[final_analysis_set$Variant_Type=="SNP" & is.na(final_analysis_set$categ) & grepl("^CG", final_analysis_set$seq_context) & final_analysis_set$Tumor_Seq_Allele1%in%c("C", "T")] <- 3 # CpG transversion
		final_analysis_set$categ[final_analysis_set$Variant_Type=="SNP" & is.na(final_analysis_set$categ) & grepl(".C.", final_analysis_set$seq_context) & final_analysis_set$Tumor_Seq_Allele1%in%c("T")] <- 2 # non CpG transition
		final_analysis_set$categ[final_analysis_set$Variant_Type=="SNP" & is.na(final_analysis_set$categ) & grepl(".G.", final_analysis_set$seq_context) & final_analysis_set$Tumor_Seq_Allele1%in%c("A")] <- 2 # non CpG transition
		final_analysis_set$categ[final_analysis_set$Variant_Type=="SNP" & is.na(final_analysis_set$categ) & grepl(".T.", final_analysis_set$seq_context) & final_analysis_set$Tumor_Seq_Allele1%in%c("C")] <- 2 # non CpG transition
		final_analysis_set$categ[final_analysis_set$Variant_Type=="SNP" & is.na(final_analysis_set$categ) & grepl(".A.", final_analysis_set$seq_context) & final_analysis_set$Tumor_Seq_Allele1%in%c("G")] <- 2 # non CpG transition
		final_analysis_set$categ[final_analysis_set$Variant_Type=="SNP" & is.na(final_analysis_set$categ) & grepl(".[CT].", final_analysis_set$seq_context) & final_analysis_set$Tumor_Seq_Allele1%in%c("G", "A")] <- 1 # non-CpG transversion
		final_analysis_set$categ[final_analysis_set$Variant_Type=="SNP" & is.na(final_analysis_set$categ) & grepl(".[GA].", final_analysis_set$seq_context) & final_analysis_set$Tumor_Seq_Allele1%in%c("C", "T")] <- 1 # non-CpG transversion
		final_analysis_set$categ[is.na(final_analysis_set$categ)] <- 5 # Indels, DNP and others.
		mutcategs <- data.frame(name=c("non-CpG tv (C/T<->A/G)", "non-CpG ts (A<->G, C<->T)", "CpG tv (C<->G)", "CpG ts (C->T, G->A)", "Indels, DNP, others"))
		if (verbose) cat(".\n")
	}
	
	if (max(final_analysis_set$categ, na.rm=TRUE) > nrow(mutcategs)) {
		if (verbose) warning("MutSig 1.5+ workaround: fixing mutcategs")
		fix.mutcategs <- mutcategs[NULL,]
		fix.mutcategs[1,"name"] <- NA
		fix.mutcategs$name <- factor("double_null")
		mutcategs <- rbind(mutcategs, fix.mutcategs)
	}
	
	if (max(final_analysis_set$categ, na.rm=TRUE) < nrow(mutcategs)) {
		final_analysis_set$categ <- factor(final_analysis_set$categ, levels=1:nrow(mutcategs))
	}
}

if (all(!c("patient_name", "patient")%in%colnames(final_analysis_set))) {
	final_analysis_set$patient_name <- factor(gsub("-Tumor", "", final_analysis_set$Tumor_Sample_Barcode))
} else if ("patient"%in%colnames(final_analysis_set)) {
	final_analysis_set$patient_name <- final_analysis_set$patient
}

final_analysis_set$patient_name <- factor(final_analysis_set$patient_name, levels=patients.counts_and_rates$name)

## Sanity checks
if (length(which(is.na(final_analysis_set$Variant_Classification_Num)))>0) {
	if (full.verbose) {
		cat(length(which(is.na(final_analysis_set$Variant_Classification_Num))),"mutations were not recognized (\"Variant_Classification\" does not fall into an expected category):\n")
		final_analysis_set[which(is.na(final_analysis_set$Variant_Classification_Num)), 
				intersect(c("Hugo_Symbol", "Tumor_Sample_Barcode", "Genome_Change", "Protein_Change", "Variant_Classification"), colnames(final_analysis_set))]
	} else warning(length(which(is.na(final_analysis_set$Variant_Classification_Num))), " SNP mutations were not recognized. Use --full.verbose to see the list.")
}

## Build gene list
if (!is.numeric(significant.gene.list$q)) {
	significant.gene.list$p <- as.numeric(sub("<", "", significant.gene.list$p))
	significant.gene.list$q <- as.numeric(sub("<", "", significant.gene.list$q))
}
if (mutsig.output!="" | firehose.mode) {
	all_sig_genes <- significant.gene.list
} else {
	all_sig_genes <- significant.gene.list
	final_analysis_set_nonsil <- subset(final_analysis_set, Variant_Classification_Num%in%c(2:7))
	final_analysis_set_sil <- subset(final_analysis_set, Variant_Classification_Num%in%c(1))
	all_sig_genes$n <- as.integer(tapply(final_analysis_set_nonsil$Tumor_Sample_Barcode, final_analysis_set_nonsil$Hugo_Symbol, length)[as.character(all_sig_genes$gene)])
	all_sig_genes$n[is.na(all_sig_genes$n)] <- 0
	all_sig_genes$npat <- as.integer(tapply(final_analysis_set_nonsil$Tumor_Sample_Barcode, final_analysis_set_nonsil$Hugo_Symbol, function(x) length(unique(x)))[as.character(all_sig_genes$gene)])
	all_sig_genes$npat[is.na(all_sig_genes$npat)] <- 0
	all_sig_genes$nsil <- as.integer(tapply(final_analysis_set_sil$Tumor_Sample_Barcode, final_analysis_set_sil$Hugo_Symbol, length)[as.character(all_sig_genes$gene)])
	all_sig_genes$nsil[is.na(all_sig_genes$nsil)] <- 0
	if(ncol(patients.counts_and_rates)==2) {
		patients.counts_and_rates$rate_non <- as.integer(tapply(final_analysis_set_nonsil$Variant_Classification, final_analysis_set_nonsil$patient_name, length)[as.character(patients.counts_and_rates$name)]) / patients.counts_and_rates$N_tot
		patients.counts_and_rates$rate_non[is.na(patients.counts_and_rates$rate_non)] <- 0
		patients.counts_and_rates$rate_sil <- as.integer(tapply(final_analysis_set_sil$Variant_Classification, final_analysis_set_sil$patient_name, length)[as.character(patients.counts_and_rates$name)]) / patients.counts_and_rates$N_tot
		patients.counts_and_rates$rate_sil[is.na(patients.counts_and_rates$rate_sil)] <- 0
		## For WGS:
#    patients.counts_and_rates$rate_non <- as.integer(tapply(final_analysis_set$Variant_Classification, final_analysis_set$patient_name, length)[as.character(patients.counts_and_rates$name)]) / patients.counts_and_rates$N_tot
#    patients.counts_and_rates$rate_sil <- 0
	}
}

if (!exists("forceSelectedGenes")) {
	selectedGenes <- subset(all_sig_genes, q <= qthresh)$gene
} else {
	selectedGenes <- intersect(all_sig_genes$gene,forceSelectedGenes)
}

# To prevent failure, select top 15 if there are no selected genes
if (length(selectedGenes) < 1) {
	selectedGenes <- all_sig_genes$gene[1:15]
}

selectedGenes <- setdiff(selectedGenes, blacklist)
sig_genes <- subset(all_sig_genes, gene%in%selectedGenes)
if ('n' %in% names(sig_genes)) {
	nonsilentColumnName = "n"
	print(nonsilentColumnName)
} else if ('nnon' %in% names(sig_genes)) {
	nonsilentColumnName = "nnon"
	print(nonsilentColumnName)
} else {
	stop("Nonsilent mutations column missing from significant genes file. Require 'n' or 'nnon' column.")
}


if(sort.genes.by.prevalence) {
	if ('n' %in% names(sig_genes)) {
		sig_genes <- sig_genes[order(sig_genes$n, -log10(sig_genes$q), sig_genes$gene, decreasing=TRUE),]
	} else {
		sig_genes <- sig_genes[order(sig_genes$nnon, -log10(sig_genes$q), sig_genes$gene, decreasing=TRUE),]
	}
}

if (verbose) cat(length(selectedGenes),"genes to display...\n q-values:\n")
if (verbose) print(summary(sig_genes$q))

## Build mutation matrix
final_analysis_subset <- subset(final_analysis_set, Hugo_Symbol%in%sig_genes$gene)

if ("miRNA"%in%(final_analysis_subset$Variant_Classification)) {
	warning("coMut does not handle miRNA mutations in significant genes yet. Please report this!")
}
final_analysis_subset <- subset(final_analysis_subset, Variant_Classification_Num<=7)
final_analysis_subset$Hugo_Symbol <- factor(final_analysis_subset$Hugo_Symbol, levels=sig_genes$gene)
final_analysis_subset <- aggregate(Variant_Classification_Num ~ Hugo_Symbol + patient_name, max, na.rm=TRUE, data=final_analysis_subset)

genematrix <- xtabs(Variant_Classification_Num ~ Hugo_Symbol + patient_name, data=final_analysis_subset, drop.unused.levels = FALSE)

## Build mutation categories matrix
if (mutation.spectrum.plot && is.null(final_analysis_set$categ) == FALSE) {
    # modification on mutation category 0 for MutSig2CV
    w.rm_indel_categ0 <- which(final_analysis_set$categ == 0)
    if(length(w.rm_indel_categ0) > 0){final_analysis_set <- final_analysis_set[-w.rm_indel_categ0, ]}

	final_analysis_set$count <- 1
	total_mut <- aggregate(count ~ patient_name, data=final_analysis_set, sum)
	mutational_signatures <- aggregate(count ~ patient_name + categ, data=final_analysis_set, sum)
	mutational_signatures <- merge(mutational_signatures, total_mut, by="patient_name")
	mutational_signatures$rate <- mutational_signatures$count.x/mutational_signatures$count.y
	mutational_signatures <- xtabs(rate ~ categ + patient_name, mutational_signatures)
	
	if (nrow(mutational_signatures) == 0) {
		stop("No mutational signatures found")
	}
	new_rownames <- mutcategs$name[as.numeric(rownames(mutational_signatures))] 
	if (length(new_rownames) != nrow(mutational_signatures)) {
		stop("Mutational categories do not match the categories provided in the MAF")
	}
	rownames(mutational_signatures) <- new_rownames
	mutational_signatures <- mutational_signatures[, order(mutational_signatures[rownames(mutational_signatures)[1],], decreasing=T), drop=FALSE]
	present_mut_categs <- rownames(mutational_signatures)
}

if (sort.by.mutation.status) {
	bin.genematrix <- rbind((genematrix[,as.character(patients.counts_and_rates$name), drop=FALSE]>1)+0, rate=rowSums(patients.counts_and_rates[,c("rate_non", "rate_sil"), drop=FALSE]))
	clustlabels <- as.character(patients.counts_and_rates$name)[do.call("order", c(lapply(1:nrow(bin.genematrix), function(x) bin.genematrix[x,]), decreasing=TRUE))]
}else if(sort.by.disease.status){
  bin.genematrix <- rbind((genematrix[,as.character(patients.counts_and_rates$name), drop=FALSE]>1)+0, rate=rowSums(patients.counts_and_rates[,c("rate_non", "rate_sil"), drop=FALSE]))
  status<-rep(0,ncol(bin.genematrix))
  
  status[as.numeric(substring(colnames(bin.genematrix),3,4))<30]<-1
  bin.genematrix<-rbind(status,bin.genematrix)
  
  clustlabels <- as.character(patients.counts_and_rates$name)[do.call("order", c(lapply(1:nrow(bin.genematrix), function(x) bin.genematrix[x,]), decreasing=TRUE))]
}else {
	clustlabels <- names(sort(rowSums(patients.counts_and_rates[,c("rate_non", "rate_sil"), drop=FALSE]), decreasing=TRUE, na.last = TRUE))
}

## Calculate allelic fractions
if (allelic.fraction.boxplot) {
	if (length(unlist(lapply(c("t_alt_count$", "t_ref_count$"), grep, colnames(final_analysis_set))))==2) {
		final_analysis_set$t_alt_count <- suppressWarnings(as.numeric(as.character(final_analysis_set[,grep("t_alt_count$", colnames(final_analysis_set))])))
		final_analysis_set$t_ref_count <- suppressWarnings(as.numeric(as.character(final_analysis_set[,grep("t_ref_count$", colnames(final_analysis_set))])))
		allelic.fraction.subset <- subset(final_analysis_set, (t_ref_count+t_alt_count)>=min.reads.to.count.allelic.fraction)
		allelic.fraction.subset$patient_name <- factor(allelic.fraction.subset$patient_name, levels=clustlabels)
	} else {
		allelic.fraction.boxplot <- FALSE
		warning("t_alt_count and t_ref_count not found in the maf file, or ambiguous, disabling the allelic fraction bloxplots!")
	}
}

mutperc <- rev(sig_genes$npat)/nrow(patients.counts_and_rates)*100

# remove silent
#genematrix[which(genematrix==1, arr.ind=TRUE)] <- 0

vsize <- nrow(sig_genes)/10+3

if (allelic.fraction.boxplot) {
	extra.vsize <- .8
} else {
	extra.vsize <- 0
}

mutsigrunsuffix <- paste(ifelse(firehose.mode | mutsig.output=="", "", "_"), ifelse(firehose.mode | mutsig.output=="", "", mutsigrun), sep="")

plot.coMut <- function() {
# X11("", 8.5, vsize+1+extra.vsize)
	
	if (verbose) cat("  initial image plot\n")
	## Mutation matrix
	layout(matrix(c(6, 4, 7, 2, 1, 3, 0, 9, 0, 8, 5, 0), ncol=3, byrow=T), widths=c(1,3,1), heights=c(1, vsize-2, extra.vsize, 2))
	par(mar=c(0,0,0,0), las=1)
	image(1:ncol(genematrix), 1:nrow(genematrix), t(genematrix[rev(as.character(sig_genes$gene)), clustlabels, drop=FALSE]), col=c("grey94", brewer.pal(7, "Set1")[c(3,6,7,2,4,5,1)]), zlim=c(0,7), xlab="", ylab="", axes=FALSE)
	segments(.5 + 1:(ncol(genematrix)-1), .5, .5 + 1:(ncol(genematrix)-1), .5 + nrow(genematrix), col="grey96", lwd=ifelse(ncol(genematrix)>200, .2, .5))
	segments(.5, .5 + 1:(nrow(genematrix)-1), .5 + ncol(genematrix), .5 + 1:(nrow(genematrix)-1), col="grey96", lwd=ifelse(ncol(genematrix)>200, .2, .5))
	
	## q-values
	par(mar=c(0,6,0,.5))
	if (verbose) cat("  q values\n")
	plot(0, xlim=c(max(-log10(sig_genes$q+0.0001)), min(.4, min(-log10(sig_genes$q+0.0001)))), type="n", yaxt="n", frame.plot=FALSE, xlab="", ylab="")
	par(usr=c(par("usr")[1:2], 0, nrow(sig_genes)), lwd=.8)
	mypos <- barplot(rev(-log10(sig_genes$q+0.001)), horiz=TRUE, axes=FALSE, add=TRUE, names.arg=rep("", nrow(sig_genes)), border="grey94", space=0, cex.names=.9,
			col=c("grey70", "grey55")[factor((rev(sig_genes$q)<=.1)+1, levels=c(1,2))])
	axis(2, at=mypos, labels=rev(sig_genes$gene), lwd=0, cex.axis=1.1, line=-.5)
	mtext("-log10(q-value)", side=1, cex=.7, line=2.2, adj=1.1)
	abline(v=-log10(.01), col="red", lwd=1)
	abline(v=-log10(.1), col="purple", lty=2, lwd=1)
	
	## Mutation counts
	par(mar=c(0,1.75, 0, 2))
	
	if (max(rowSums(sig_genes[,c(nonsilentColumnName, "nsil")]))>60) {
	  maxCountSigGenes <- as.integer(10^(quantile(log10(rowSums(sig_genes[,c(nonsilentColumnName, "nsil")])), 3/4) + 1.5*IQR(log10(rowSums(sig_genes[,c(nonsilentColumnName, "nsil")])))))
	} else {
	  maxCountSigGenes <- max(rowSums(sig_genes[,c(nonsilentColumnName, "nsil")]))
	}
	if (verbose) cat("  counts\n")
	plot(0, xlim=c(0, maxCountSigGenes), type="n", axes=FALSE, frame.plot=FALSE, xlab="", ylab="")
	par(usr=c(par("usr")[1:2], 0, nrow(sig_genes)), lwd=.8)
	mypos <- barplot(t(sig_genes[nrow(sig_genes):1,c(nonsilentColumnName, "nsil")]), col=c("dodgerblue4", "#4DAF4A"), 
	                 horiz=TRUE, names.arg=rep("", nrow(sig_genes)), add=TRUE, border="grey94", space=0, axes=FALSE)
	axis(2, at=mypos, labels=paste(round(mutperc,1), "%", sep=""), line=-1.25, cex.axis=1.0, tick=FALSE)
	if (maxCountSigGenes < max(rowSums(sig_genes[,c(nonsilentColumnName, "nsil")])))
	  text(rep(maxCountSigGenes*.9, 2), which(rev(rowSums(sig_genes[,c(nonsilentColumnName, "nsil")]))>maxCountSigGenes)-.5, 
	       cex=.8, col="white", labels=rev(rowSums(sig_genes[,c(nonsilentColumnName, "nsil")]))[which(rev(rowSums(sig_genes[,c(nonsilentColumnName, "nsil")]))>maxCountSigGenes)])
	par(new=TRUE, lwd=1)
	plot(0, xlim=c(0, maxCountSigGenes), type="n", yaxt="n", frame.plot=FALSE, xlab="", ylab="")
	mtext("# mutations", side=1, cex=.7, line=2)
	
	## Mutation rates
	par(mar=c(.5,0,.5,0), las=1)
#plot(1, ylim=range(rowSums(patients.counts_and_rates[clustlabels,c("rate_non", "rate_sil")], na.rm=TRUE)*1e6), type="n", axes=F, xlab="", ylab="")
	if (verbose) cat("  mutations rates\n")
	plot(1, ylim=c(0,min(max.mutation.rate, max(rowSums(patients.counts_and_rates[clustlabels,c("rate_non", "rate_sil")], na.rm=TRUE)*1e6))), type="n", axes=F, xlab="", ylab="")
	par(usr=c(0, nrow(patients.counts_and_rates), par("usr")[3:4]), lwd=ifelse(ncol(genematrix)>200, .2, .5))
	barplot(t(patients.counts_and_rates[clustlabels,c("rate_non", "rate_sil")])*1e6, col=c("dodgerblue4", "#4DAF4A"), axes=FALSE, add=TRUE, names.arg=rep("", nrow(patients.counts_and_rates)), border="grey94", space=0)
	par(lwd=1)
	axis(2, cex.axis=1, line=.3)
	if (any(rowSums(patients.counts_and_rates[clustlabels,c("rate_non", "rate_sil")]*1e6)>max.mutation.rate))
		text(which(rowSums(patients.counts_and_rates[clustlabels,c("rate_non", "rate_sil")]*1e6)>max.mutation.rate)-.5, max.mutation.rate*.9, 
				cex=.8, col="white", labels=round(rowSums(patients.counts_and_rates[clustlabels,c("rate_non", "rate_sil")]*1e6))[which(rowSums(patients.counts_and_rates[clustlabels,c("rate_non", "rate_sil")]*1e6)>max.mutation.rate)],
				srt=90)
	mtext("# mutations/Mb", side=2, cex=.7, line=2.7, las=0)
	
	## Mutation categories
	if (verbose) cat("  mutations categories\n")
	plot(0, ylim=c(0,100), type="n", axes=FALSE, xlab="", ylab="")
	if (mutation.spectrum.plot && is.null(final_analysis_set$categ) == FALSE) {
		par(mar=c(7,0,.5,0), las=1)
		par(usr=c(0, ncol(mutational_signatures), 0,100))
		box()
		mutation_spectrum_plot_brewer <- brewer.pal(length(present_mut_categs) +1, "Spectral")[seq_along(present_mut_categs)]
		mypos <- barplot(mutational_signatures[,clustlabels]*100, beside=FALSE, col=mutation_spectrum_plot_brewer, names.arg=rep("", ncol(mutational_signatures)), add=TRUE, border=NA, space=0, axes=FALSE)
		axis(4, cex.axis=1, line=.3)
		mtext("%", side=4, cex=.8, line=2.5, las=1)
		par(las=2)
		axis(1, at=mypos, labels=clustlabels, tick=FALSE, cex.axis=min((ifelse(length(clustlabels)<=130, .2, 0) + .7/log10(length(clustlabels))), 1), line=-.7)
	}
	
	## Display legends
	if (verbose) cat("  legends\n")
	par(mar=c(.3,0,.5,0), las=1)
	plot(0, type="n", axes=FALSE, xlab="", ylab="")
	legend("center", "center", c("Syn.", "Non syn."), fill=c("#4DAF4A", "dodgerblue4"), cex=.9, bty = "n")
	
	plot(0, type="n", axes=FALSE, xlab="", ylab="")
	par(xpd=NA)
	legend("right", "top", c("Syn.", "Missense", "Splice site", "Nonsense", "Frame shift", "In frame indel", "Other non syn."),
			fill=brewer.pal(7, "Set1")[c(3,2,4,1,5,6,7)], cex=.9, ncol=2, bty = "n", inset=0)
	par(xpd=FALSE)
	
	par(mar=c(.5, 0, .5, 0.5), las=1)
	plot(0, type="n", axes=FALSE, xlab="", ylab="")
	if (mutation.spectrum.plot && is.null(final_analysis_set$categ) == FALSE) {
		mutation.spectrum.legend <- rev(rownames(mutational_signatures))
		legend("right", ifelse(allelic.fraction.boxplot,"top", "center"), 
				mutation.spectrum.legend,
				fill=rev(mutation_spectrum_plot_brewer), cex=.9, ncol=1, inset=0, bty = "n")
	}
	
	## Allelic fraction boxplot
	if (allelic.fraction.boxplot) {
		if (verbose) cat("  allelic fraction boxplot\n")
		plot(0, ylim=c(0,100), type="n", axes=FALSE, xlab="", ylab="")
		par(mar=c(0.1, 0, .5, 0), las=1)
		par(usr=c(0, length(clustlabels), 0,100))
		box()
		b <- boxplot(allelic.fraction.subset$t_alt_count/(allelic.fraction.subset$t_alt_count+allelic.fraction.subset$t_ref_count)*100 ~ allelic.fraction.subset$patient_name, 
				pch=20, cex=.3, range=0, lty=1, col=brewer.pal(3, "Set1")[2], ylim=c(0,1), cex.axis=.6, axes=FALSE, add=TRUE, at=(1:length(clustlabels))-.5, plot=FALSE)
		segments(x0=(1:length(clustlabels))-.5, y0=b$stats[1,], y1=b$stats[5,], col=brewer.pal(4, "Purples")[3], lwd=ifelse(ncol(genematrix)>200, .8, 1.5))
		segments(x0=(1:length(clustlabels))-.5, y0=b$stats[2,], y1=b$stats[4,], col=brewer.pal(4, "Purples")[4], lwd=ifelse(ncol(genematrix)>200, 1.5, 3))
		abline(h=median(allelic.fraction.subset$t_alt_count/(allelic.fraction.subset$t_alt_count+allelic.fraction.subset$t_ref_count)*100), col="red", lwd=1, lty=2)
		axis(4, cex.axis=1, line=.3)
		mtext("Allelic\nfraction", side=4, cex=.7, line=3.6, las=0)
	}
	
}


if (verbose) cat("Plotting...\n")
if (!png.output | firehose.mode) {
	if (as.numeric(sub("(.*)\\..*","\\1",R.Version()$minor))>=14) {
		if (full.verbose) cat("  Using Cairo\n")
		cairo_pdf(paste(outputdir, "/", analysis.set, mutsigrunsuffix, "_coMut.pdf", sep=""), 8.5, vsize+1+extra.vsize) ## Only in R 2.14 and up
		plot.coMut()
		dev.off()
	} else {
		pdf(paste(outputdir, "/", analysis.set, mutsigrunsuffix, "_coMut.pdf", sep=""), 8.5, vsize+1+extra.vsize)
		plot.coMut()
		dev.off()
	}
} 

if (png.output | firehose.mode) {
	png(paste(outputdir, "/", analysis.set, mutsigrunsuffix, "_coMut.png", sep=""), 8.5*120, (vsize+1+extra.vsize)*120, type="cairo", pointsize=20)
	plot.coMut()
	dev.off()
}

