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
groupFile = args[1] # "01.data/03.MCL/02.mcl/mcl.out.txt"
outFile = args[2] # "01.data/04.GeneOfInterest/GeneID.txt"
size = args[3] # 1000
sp = args[4] # "AT"

###########################################################################################
###  Combine paralog groups into meta-group of preset size.								###
###########################################################################################

con = file(groupFile, "r", blocking = FALSE)
data = readLines(con)
close(con)

data.gene = lapply(data, strsplit, "\\s+")
group = list()
i = 0
while(length(data.gene) > 0){
	i = i + 1
	temp = unlist(data.gene[1])
	temp = temp[grep(sp, temp)]
	data.gene = data.gene[-1]
	#n = sample(1:length(data.gene))
	n = 1:length(data.gene)
	m = NULL
	for(j in n){
		x = unlist(data.gene[j])
		x = x[grep(sp, x)]
		if(length(temp) + length(x) <= size){
			temp = c(temp, x)
			m = c(m, j)
		}
	}
	if(!is.null(m)){
		data.gene = data.gene[-m]
	}
	group[[i]] = temp
}

###########################################################################################
###  Write to output file. 																###
###########################################################################################

temp.dir = unlist(strsplit(outFile, "/"))
outDir = paste(temp.dir[-length(temp.dir)], collapse="/")
system(paste("mkdir -p ", outDir, sep=""))
system(paste(":> ", outFile))

for(i in 1:length(group)){
	write(group[[i]], file=outFile, sep="\n", append=TRUE)
	write("", file=outFile, append=TRUE)
}

###########################################################################################
###  Clean and release the memory before start R										###
###########################################################################################

rm(list=ls())
gc()
