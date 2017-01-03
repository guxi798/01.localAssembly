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
#geneFile = args[3]		# "01.data/05.splitGenes/02.Transcript/run.3/contig2gene.txt"
tissueFile = args[3]		# "01.data/00.PriorData/tissue_record.txt"
method = args[4]		# method of estimate expression

###########################################################################################
###  Read in annotation file. 															###
###########################################################################################
#workdir = "/escratch4/guxi798/guxi798_Sep_17/14.localAssembly/12.Parasite/01.OrAe"
#setwd(workdir)
system(paste("rm -f", outFile))
#genes = read.table(geneFile, header=FALSE, sep="", as.is=TRUE)
subs = list.files(exprFolder, pattern=".+", full.names=FALSE)
tissues = read.table(tissueFile, header=TRUE, as.is=TRUE)

i = 0
data = NULL
count = NULL
if(method == "RSEM"){
	for (sub in subs){
		file = paste(exprFolder, sub, "RSEM.isoforms.results", sep="/")
		temp = read.table(file, header=TRUE, sep="\t", as.is=TRUE)
		if(i == 0){
			data = temp[,c(1,3,7)]
			count = temp[,6]
			#pos = match(data, genes[,1])
			#data = cbind(data, genes[pos,2])
			#data = cbind(data, temp[,7])
		}else{
			pos = match(data[,1], temp[,1])
			data = cbind(data, temp[,7])
			count = cbind(count, temp[,6])
		}
		i = i+1
	}
	data = cbind(data, count)
	pos = match(subs, tissues[,1])
	name1 = paste(tissues[pos,2], "FPKM", sep="_")
	name2 = paste(tissues[pos,2], "TPM", sep="_")
	colnames(data) = c("transcript_id", "transcript_len", name1, name2)
}else if(method == "eXpress"){
	for (sub in subs){
		file = paste(exprFolder, sub, "results.xprs", sep="/")
		temp = read.table(file, header=TRUE, sep="\t", as.is=TRUE)
		if(i == 0){
			data = temp[,c(2,3,11)]
			count = temp[,5]
			#pos = match(data, genes[,1])
			#data = cbind(data, genes[pos,2])
			#data = cbind(data, temp[,7])
		}else{
			pos = match(data[,1], temp[,2])
			data = cbind(data, temp[pos,11])
			count = cbind(count, temp[pos, 5])
		}
		i = i+1
	}
	data = cbind(data, count)
	pos = match(subs, tissues[,1])
	name1 = paste(tissues[pos,2], "FPKM", sep="_")
	name2 = paste(tissues[pos,2], "Count", sep="_")
	colnames(data) = c("transcript_id", "transcript_len", name1, name2)
}

###########################################################################################
###  Write to output file. 																###
###########################################################################################

write.table(data, file=outFile, quote=FALSE, sep="\t", col.names=TRUE, row.names=FALSE)

###########################################################################################
###  Clean and release the memory before start R										###
###########################################################################################

rm(list=ls())
gc()
