library("devEMF")

globfile = "05.Full.Length/global.uniq.gene.txt"
locfile = "05.Full.Length/local.uniq.gene.txt"
bothfile = "05.Full.Length/both.full.gene.txt"

## read in data
data.glob = read.table(globfile, header=FALSE, as.is=TRUE)
data.loc = read.table(locfile, header=FALSE, as.is=TRUE)
data.both = read.table(bothfile, header=FALSE, as.is=TRUE)

## properties
## protein length
summary(data.glob[,2])
summary(data.loc[,2])
## transcript length
summary(data.glob[,3])
summary(data.loc[,3])
## gc content
summary(data.glob[,4])
summary(data.loc[,4])
## expression
summary(data.glob[,5])
summary(data.loc[,5])
x = data.glob[data.glob[,5]!=0,5]
y = data.loc[data.loc[,5]!=0,5]
z = data.both[data.both[,5]!=0,5]
m = c(x, z, y)
name = factor(c(rep("Global",length(x)), rep("Common", length(z)), rep("Local", length(y))),
			levels = c("Global", "Common", "Local"))

tgtfolder = "13.figures/02.4Kand2K"
system(paste("mkdir -p ", tgtfolder, sep=""))
figureFile = paste(tgtfolder, "/", "uniqueset.tranlength.emf", sep="")
	emf(file = figureFile, width = 4, height = 4, family="Helvetica")
	boxplot(data.glob[,3], data.loc[,3], ylab="Transcript Length", names=c("Global", "Local"))
dev.off()

figureFile = paste(tgtfolder, "/", "uniqueset.gc.emf", sep="")
	emf(file = figureFile, width = 4, height = 4, family="Helvetica")
	boxplot(data.glob[,4], data.loc[,4], ylab="Transcript Length", names=c("Global", "Local"))
dev.off()

figureFile = paste(tgtfolder, "/", "uniqueset.expression.emf", sep="")
	emf(file = figureFile, width = 4, height = 4, family="Helvetica")
	boxplot(log(x), log(y), ylab="Log-transformed FPKM", names=c("Global", "Local"))
dev.off()

figureFile = paste(tgtfolder, "/", "set.expression.emf", sep="")
	emf(file = figureFile, width = 4, height = 4, family="Helvetica")
	boxplot(log(m) ~ name, ylab="Log-transformed FPKM")
dev.off()

## both properties
globgene = unique(c(data.glob[,1], data.both[,1]))
locgene = unique(c(data.loc[,1], data.both[,1]))

file = "../03.Protein.Subset/01.data/04.GeneOfInterest/GeneID.v2.txt"
data = read.table(file, header=TRUE, as.is=TRUE)

set = list()
set = c(set, list(data$Gene_ID[data$expr > 0 & data$expr <= 1]))
set = c(set, list(data$Gene_ID[data$expr > 1 & data$expr <= 10]))
set = c(set, list(data$Gene_ID[data$expr > 10 & data$expr <= 50]))
set = c(set, list(data$Gene_ID[data$expr > 50 & data$expr <= 100]))
set = c(set, list(data$Gene_ID[data$expr > 100 & data$expr <= 200]))
set = c(set, list(data$Gene_ID[data$expr > 200 & data$expr <= 500]))
set = c(set, list(data$Gene_ID[data$expr > 500]))

all = NULL
glob.num = NULL
loc.num = NULL
n = 7
for(i in 1:n){
	all = c(all, length(set[[i]]))
	glob.num = c(glob.num, sum(!is.na(match(globgene, set[[i]]))))
	loc.num = c(loc.num, sum(!is.na(match(locgene, set[[i]]))))
}

tgtfolder = "13.figures/03.full.length.recover"
system(paste("mkdir -p ", tgtfolder, sep=""))

figureFile = paste(tgtfolder, "/", "fraction.against.fpkm.emf", sep="")
	emf(file = figureFile, width = 4, height = 4, family="Helvetica")
	plot(x=1:n, y=loc.num/all, type="b", col="red", xaxt="n", xlab="FPKM range", ylab="Recovery Rate")
	points(x=1:n, y=glob.num/all, type="b")
	text(1:n-0.05, par("usr")[3]-0.03, labels= c("0-1", "1-10", "10-50", "50-100", "100-200", "200-500", ">500"),
		srt = 45, pos=1, xpd=TRUE)
	legend("bottomright", legend=c("local", "global"), col=c("red", "black"), bty="n", lty=1, pch=1)
dev.off()

################################################
tgtfolder = "13.figures/01.identity"
system(paste("mkdir -p ", tgtfolder, sep=""))

## compare common set between global and local
## identity
file = "05.Full.Length/both.full.identity.txt"
data = read.table(file, header=FALSE, as.is=TRUE, row.names=1)
data1= data[data[,2]!=0 & data[,4]!=0,]		## 0 usually means blastx is consistent with blastn

data = data1[,c(2,4)]
colnames(data) = c("global", "local")
figureFile = paste(tgtfolder, "/", "commonset.identity.emf", sep="")
	emf(file = figureFile, width = 4, height = 4, family="Helvetica")
	boxplot(data, ylab="Identity")
dev.off()

data = data1[,c(1,3)]

wilcox.test(data1[,1], data1[,2], alternative="less") # p-value < 2.2e-16
t.test(data1[,1], data1[,2], alternative="less") # p-value < 2.2e-16

## ordering
data1 = data1[order(data1[,2]),]
plot(data1[,1], col="gray")
points(data1[,2])
