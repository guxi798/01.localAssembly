###########################################################################################
###  Clean and release the memory before start R										###
###########################################################################################

rm(list=ls())
gc()

###########################################################################################
###  Load necessary R libraries or functions for later usage.							###
###########################################################################################
if(!require("devEMF")){
	install.packages("devEMF", dependencies = TRUE)
	library(devEMF)
}

###########################################################################################
###  Read in arguments. 																###
###########################################################################################

cat("\n\n########## Reading Arguments ##########\n", sep="")
args = commandArgs(TRUE)
for (arg in args) cat("  ", arg, "\n", sep="")
groupFile = "01.data/03.MCL/02.Transcript/02.mcl/wu.mcl.out.txt"
exprFile = "01.data/00.PriorData/genes.fpkm_tracking.deseq"
tissueFile = "01.data/00.PriorData/tissue_record.txt"
RDataDir = "02.expr.pattern"
size = 500

system(paste("mkdir -p ", RDataDir))
###########################################################################################
###  Read in annotation file. 															###
###########################################################################################
## this function is for fpkm file from cufflinks
read.fpkm <- function(file, sample.file=NULL){
    data = read.table(file=file, header=TRUE, sep='\t', as.is=T)
    ncol = dim(data)[2] #get column no
    col = seq(from=10, to=ncol, by=4) #1:tracking_id, 4:gene_id, 10-ncol: fpkm value
    fpkm = data[,col]
    rownames(fpkm) = data$tracking_id

    if(!is.null(sample.file)){
        tissues = read.csv(sample.file, header=TRUE, sep="\t", as.is=TRUE)
        names(fpkm) = tissues[,2]
    }else{
            print("Warning: No tissue names specified.")
        }

    return(fpkm)
}

## for deseq, fpkm file can directly be read
## read in gene family data
con = file(groupFile, "r", blocking = FALSE)
data = readLines(con)
close(con)

gene = NULL
group = NULL
for(i in 1:length(data)){
	temp = unlist(strsplit(data[[i]], "\t"))
	group.temp = rep(i, times = length(temp))
	gene = c(gene, temp)
	group = c(group, group.temp)
}

### read in expression file 
#fpkm = read.fpkm(exprFile, tissueFile)
fpkm = read.table(exprFile, as.is=TRUE)

loc = match(gene,rownames(fpkm))
loc.na = which(is.na(loc))
# no expression for "Potri.004G166300" "Potri.004G166100" "Potri.T143800"
if(length(loc.na) > 0){
	gene.new = gene[-loc.na]
	group.new = group[-loc.na]
}else{
	gene.new = gene
	group.new = group
}
loc = match(gene.new, rownames(fpkm))
fpkm.new = fpkm[loc, ]
#fpkm.new[fpkm.new==0] = 1e-7
#fpkm.new = log10(fpkm.new)
aver.expr = apply(fpkm.new, 1, mean)

### plot expression vs. group
file = paste(RDataDir, paste("expr.group.local100",".emf",sep=""), sep="/")
emf(file, height=6, width=6)
par(mar = c(5,5,2,2))
plot(y=aver.expr, x=group.new, xlab="Group ID", ylab="FPKM", main="", 
		ylim=c(0, 100), pch=20, cex=0.5, cex.axis=1.5, cex.lab=1.5)
dev.off()

### normalized by gene number
mean.expr = NULL
retain = NULL
for(i in 1:length(data)){
	pos = which(group.new == i)
	mean.expr =  c(mean.expr, mean(aver.expr[pos]))
	retain = c(retain, length(pos))
}

file = paste(RDataDir, paste("aver.expr.group.local100",".emf",sep=""), sep="/")
emf(file, height=6, width=6)
par(mar = c(5,5,2,2))
plot(y=mean.expr, x=1:length(data), xlab="Group ID", ylab="FPKM", main="", 
		ylim=c(0, 100), cex=0.5, cex.axis=1.5, cex.lab=1.5)
dev.off()

file = paste(RDataDir, paste("aver.logexpr.group",".emf",sep=""), sep="/")
emf(file, height=6, width=6)
par(mar = c(5,5,2,2))
plot(y=log(mean.expr), x=1:length(data), xlab="Group ID", ylab="Log FPKM", main="", 
		cex=0.5, cex.axis=1.5, cex.lab=1.5)
dev.off()

### relationship between gene number and expression
num = floor(length(data)/size)
frac = NULL
for(i in c(1, 5, 10, 20)){
	frac.temp = NULL
	for(j in 1:num){
		low = sum(mean.expr[seq(((j-1)*size+1),j*size)] > i)
		frac.temp = c(frac.temp, low/size)
	}
	frac = cbind(frac, frac.temp)
}

file = paste(RDataDir, paste("high.expr.fraction",".emf",sep=""), sep="/")
emf(file, height=6, width=6)
par(mar = c(5,5,2,2))
plot(x=1:num, y=frac[,1], xlab="Bin ID", ylab="High Expression Fraction", main="", 
		ylim=c(0,1), type="n", cex.axis=1.5, cex.lab=1.5)
	points(x=1:num, y=frac[,1], type="b", pch = 1)
	points(x=1:num, y=frac[,2], type="b", pch = 2)
	points(x=1:num, y=frac[,3], type="b", pch = 3)
	points(x=1:num, y=frac[,4], type="b", pch = 4)
legend("topright", legend=c("cutoff = 1", "cutoff = 5", "cutoff = 10", "cutoff = 20"),
		pch = c(1,2,3,4), lty=c(1,1,1), bty="n")
dev.off()


###########################################################################################
###  Clean and release the memory before start R										###
###########################################################################################

rm(list=ls())
gc()
