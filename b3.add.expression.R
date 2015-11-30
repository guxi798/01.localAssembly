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
geneFile = "01.data/04.GeneOfInterest/GeneID.v1.txt"
exprFile = "01.data/00.PriorData/genes.fpkm_tracking.deseq"
outFile = "01.data/04.GeneOfInterest/GeneID.v2.txt"

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
### read in expression file 
fpkm = read.table(exprFile, as.is=TRUE)
sub = c(47,50,53,56,59,62,63,66,69,72,76,78)
aver.expr = apply(fpkm[, sub], 1, mean)

## read in gene info data
data = read.csv(geneFile, header=TRUE, as.is=TRUE, sep="\t")
loc = match(data[,1], rownames(fpkm))
expr = round(aver.expr[loc], digits=2)
data = cbind(data, expr)

## write data
write.table(data, file=outFile, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

###########################################################################################
###  Clean and release the memory before start R										###
###########################################################################################

rm(list=ls())
gc()
