rm(list=ls())
gc()

if(!require("devEMF")){
	install.packages("devEMF", dependencies = TRUE)
	library(devEMF)
}

args = commandArgs(TRUE)
for (arg in args) cat("  ", arg, "\n", sep="")

blsFile = tolower(args[1])				## 08.evaluate/03.bowtie.nucl/run.1/1000/1000.compare.tab
tgtfolder = args[2]
figureHead = tolower(args[3])

data = read.table(file = blsFile, header=TRUE, as.is=TRUE)

figureFile = paste(tgtfolder, "/", figureHead, ".identity.emf", sep="")
emf(file = figureFile, width = 6, height = 6)
#png(file = figureFile, width = 480, height = 480, family="Helvetica")
plot(x=data$OldRun_identity, y=data$NewRun_identity, 
		main="Blast Identity (%)",
		xlab="Old Run",	ylab="New Run")
dev.off()

figureFile = paste(tgtfolder, "/", figureHead, ".length.emf", sep="")
emf(file = figureFile, width = 6, height = 6, family="Helvetica")
#png(file = figureFile, width = 480, height = 480, family="Helvetica")
plot(x=data$OldRun_len, y=data$NewRun_len, 
		main="Blast Length",
		xlab="Old Run",	ylab="New Run")
dev.off()

figureFile = paste(tgtfolder, "/", figureHead, ".evalue.emf", sep="")
emf(file = figureFile, width = 6, height = 6, family="Helvetica")
#png(file = figureFile, width = 480, height = 480, family="Helvetica")
plot(x=log(data$OldRun_evalue), y=log(data$NewRun_evalue), 
		main="Blast Evalue (log)",
		xlab="Old Run",	ylab="New Run")
dev.off()
