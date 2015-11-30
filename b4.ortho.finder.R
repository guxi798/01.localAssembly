###########################################################################################
###  Clean and release the memory before start R										###
###########################################################################################

rm(list=ls())
gc()

###########################################################################################
###  Read in arguments. 																###
###########################################################################################

cat("\n\n########## Reading Arguments ##########\n", sep="")
args = commandArgs(TRUE)
for (arg in args) cat("  ", arg, "\n", sep="")
groupFile = args[1] # "01.data/04.GeneOfInterest/GeneID.txt"
outFile = args[2] # "01.data/04.GeneOfInterest/GeneID.ref.txt"
refFile = args[3] # "01.data/00.PriorData/ortholog_groups_edit.tab"
size = args[4] # 1000
sp = args[5] # "AT"
sp2 = args[5] # "Potri"

system(paste(":> ", outFile, sep=""))

###########################################################################################
###  Combine paralog groups into meta-group of preset size.								###
###########################################################################################

con = file(refFile, "r", blocking = FALSE)
data = readLines(con)
close(con)

data.gene = lapply(data, strsplit, "\\s+")
group = NULL
refgroup = list()
j=0
for(i in 1:length(data.gene)){
	j = j + 1
	temp = unlist(data.gene[i])
	temp1 = temp[grep(sp, temp)]
	temp2 = temp[grep(sp2, temp)]
	temp1
	match(x, temp1)
	if(length(temp1) == 0){
		j = j - 1
		next
	}
	group = rbind(group, cbind(temp1, j))
	refgroup = c(refgroup, list(temp2))
}

con = file(groupFile, "r", blocking = FALSE)
ortho = readLines(con)
close(con)

i = 0
while(length(ortho) > 0){
	i = i + 1
	g = NULL
	gene = 0
	while(gene != "" && length(ortho) > 0){
		gene = ortho[1]
		ortho = ortho[-1]
		g = c(g, gene)
	}
	loc = match(g, group[,1])
	gID = unique(group[loc, 2])
	gID = as.numeric(gID[!is.na(gID)])
	orthogenes = unlist(refgroup[gID])
	match(x, g)
	write(orthogenes, file=outFile, sep="\n", append=TRUE)
	write("", file=outFile, sep="\n", append=TRUE)
}

###########################################################################################
###  Clean and release the memory before start R										###
###########################################################################################

rm(list=ls())
gc()
