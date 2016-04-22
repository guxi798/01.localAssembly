###########################################################################################
###  Clean and release the memory before start R										###
###########################################################################################

rm(list=ls())
gc()

###########################################################################################
###  Load necessary R libraries or functions for later usage.							###
###########################################################################################


###########################################################################################
###  Read in arguments. 																###
###########################################################################################

cat("\n\n########## Reading Arguments ##########\n", sep="")
args = commandArgs(TRUE)
for (arg in args) cat("  ", arg, "\n", sep="")
exprFolder = args[1] 		# "13.abundance/run.3"
outFile = args[2]		# "13.abundance/run.3/gene.fpkm.txt"
geneFile = args[3]		# "01.data/05.splitGenes/02.Transcript/run.3/contig2gene.txt"
tissueFile = args[4]		# "01.data/00.PriorData/tissue_record.txt"

###########################################################################################
###  Read in annotation file. 															###
###########################################################################################
#workdir = "/escratch4/guxi798/guxi798_Sep_17/14.localAssembly/12.Parasite/01.OrAe"
#setwd(workdir)

genes = read.table(geneFile, header=FALSE, sep="", as.is=TRUE)
subs = list.files(exprFolder, pattern=".+", full.names=FALSE)
tissues = read.table(tissueFile, header=TRUE, as.is=TRUE)

i = 0
data = NULL
for (sub in subs){
	file = paste(exprFolder, sub, "RSEM.isoforms.results", sep="/")
	temp = read.table(file, header=TRUE, sep="\t", as.is=TRUE)
	if(i == 0){
		data = temp[,1]
		pos = match(data, genes[,1])
		data = cbind(data, genes[pos,2])
		data = cbind(data, temp[,7])
	}else{
		data = cbind(data, temp[,7])
	}
	i = i+1
}

pos = match(subs, tissues[,1])
colnames(data) = c("gene_id", "gene_name", tissues[pos, 2])

###########################################################################################
###  Write to output file. 																###
###########################################################################################

write.table(data, file=outFile, quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)

###########################################################################################
###  Clean and release the memory before start R										###
###########################################################################################

rm(list=ls())
gc()
