rm(list=ls())
gc()

if(!require("devEMF")){
	install.packages("devEMF", dependencies = TRUE)
	library(devEMF)
}
if(!require("leaps")){
	install.packages("leaps", dependencies = TRUE)
	library(leaps)
}
if(!require("RColorBrewer")){
	install.packages("RColorBrewer", dependencies = TRUE)
	library(RColorBrewer)
}


args = commandArgs(TRUE)
for (arg in args) cat("  ", arg, "\n", sep="")

################################################
fullSeqFolder = "01.data/05.splitGenes/03.Full.Length/"
allSeqFolder = "01.data/05.splitGenes/02.Transcript/"
blastFolder = "07.map.back/02.blastn/"
tgtfolder = "13.figures"

## local
files = list.files(fullSeqFolder)
runs = matrix(unlist(strsplit(files[grep("run.*", files)], "\\.")), ncol=2, byrow=TRUE)[,2]
runs = sort(as.numeric(runs))

full.count = NULL
all.count = NULL
all.len = NULL
for(run in runs){
	if(run == 0){
		next
	}
	print(run)
	file = paste(fullSeqFolder, "run.", run, "/full.length.contigs.fasta", sep="")
	full.temp = system(paste("grep '>' ", file, " | wc -l", sep=""), intern = TRUE)
	full.count = c(full.count, as.numeric(full.temp))
	
	file = paste(allSeqFolder, "run.", run, "/all.length.fasta", sep="")
	all.temp = system(paste("grep '>' ", file, " | wc -l", sep=""), intern = TRUE)
	all.count = c(all.count, as.numeric(all.temp))
	
	data = read.table(file, as.is = TRUE)
	len.temp = mean(as.numeric(gsub("(len=)", "", data[,2], perl=TRUE)))
	all.len = c(all.len, len.temp)
}

figureFile = paste(tgtfolder, "/", "numoffull.run.emf", sep="")
emf(file = figureFile, width = 6, height = 6, family="Helvetica")
	plot(runs, full.count, xlab="Iteration", ylab="Count", type="b", main="Fully Assembled Transcripts")
dev.off()

figureFile = paste(tgtfolder, "/", "alllen.run.emf", sep="")
emf(file = figureFile, width = 6, height = 6, family="Helvetica")
	plot(all.len, xlab="Iteration", ylab="Base Pairs", type="b", main="Length of Transcripts")
dev.off()

####################################
## local
run = 6
subs = list.files(paste(allSeqFolder, "run.", run, sep=""))
subs = subs[grep("[0-9]+", subs)]

len.local = NULL
frac.local = NULL
prob = 0.75
for(sub in subs){
	file = paste(allSeqFolder, "run.", run, "/", sub, "/header.txt", sep="")
	seqinfo = read.csv(file, sep=" ", header=FALSE, as.is=TRUE)
	seqinfo[,1] = gsub("(>)", "", seqinfo[,1], perl=TRUE)
	seqinfo[,2] = as.numeric(gsub("(len=)", "", seqinfo[,2], perl=TRUE))
	seqinfo[,3] = as.numeric(gsub("(GC=)", "", seqinfo[,3], perl=TRUE))
					
	file = paste(blastFolder, "run.", run-1, "/", sub, "/", sub, ".new.contigs.blast.out", sep="")
	blast = read.table(file, sep="\t", header=FALSE, as.is=TRUE)
	loc = match(seqinfo[,1], blast[,1])
	cut = quantile(seqinfo[,2], prob=prob)
	pos = seqinfo[,2] >= cut
	
	temp = blast[loc[pos],4]/seqinfo[pos,2]
	temp = temp[!is.na(temp)]
	len.local = c(len.local, seqinfo[pos,2])
	frac.local = c(frac.local, temp)
}

summary(len.local)
round(summary(frac.local), digits = 4)

breaks = 50
interval = 1/breaks

x.local = NULL
y.local = NULL
for(i in 1:breaks){
	lower = (i-1)*interval
	upper = i*interval
	c = sum(frac>lower & frac<=upper)/length(frac)
	x.local = c(x.local, (lower+upper)/2)
	y.local = c(y.local, c)
}

plot(y.local ~ x.local, type="l", col="red")

###############################

## global
file = "../08.Trinity/05.Full.Length/header.txt"
seqinfo = read.csv(file, sep=" ", header=FALSE, as.is=TRUE)
seqinfo[,1] = gsub("(>)", "", seqinfo[,1], perl=TRUE)
seqinfo[,2] = as.numeric(gsub("(len=)", "", seqinfo[,2], perl=TRUE))
#seqinfo[,3] = as.numeric(gsub("(GC=)", "", seqinfo[,3], perl=TRUE))

prob = 0.75
file = "../08.Trinity/04.blastn/blast.new.out"
blast = read.table(file, sep="\t", header=FALSE, as.is=TRUE)
loc = match(seqinfo[,1], blast[,1])
cut = quantile(seqinfo[,2], prob=prob)
pos = seqinfo[,2] >= cut
temp = blast[loc[pos],4]/seqinfo[pos,2]

len.global = seqinfo[pos,2]
frac.global = temp[!is.na(temp)]

summary(len.global)
round(summary(frac.global), digits = 4)

x.global = NULL
y.global = NULL
for(i in 1:breaks){
	lower = (i-1)*interval
	upper = i*interval
	c = sum(frac>lower & frac<=upper)/length(frac)
	x.global = c(x.global, (lower+upper)/2)
	y.global = c(y.global, c)
}

lines(y.global ~ x.local)

hist(len.local)
hist(len.global)
########################################

fpkm.cls2 = fpkm <= 0 & fpkm > -4
all.full = data[full.loc,11:20]*100
all.cls2 = all.full[fpkm.cls2, ]
plot(seq(0, 9), data[1,11:20], ylim=c(0, 100), type="n")
for(i in 1:sum(fpkm.cls2)){
	points(seq(0, 9), all.cls2[i,], type="b")
}

nclass = 10
col = colorRampPalette(c("gray", "black"))(n=nclass)
col = palette(rainbow(n=nclass))

bound1 = quantile(fpkm, probs=seq(0,1,1/nclass))
bound2 = quantile(seq(min(fpkm), max(fpkm)+0.01, by=0.01), probs=seq(0,1,1/nclass))
class = cut(fpkm, bound2, label=1:nclass, include.lowest=TRUE)
plot(seq(0, 9), data[1,11:20], ylim=c(40, 100), type="n", xlab="Run No", ylab="Recovered Percentage")
for(i in 1:nclass){
	all.cls = all.full[class == i,]
	all.avg = apply(all.cls, 2, mean)
	points(seq(0,9), all.avg, type="b", col=col[i])
}
legend("bottomright", legend=paste("Class", 1:nclass, sep="_"), bty="n", col=col, lty=1, pch=1, ncol=2)

########################################
#fpkm = as.factor(round(log(data$FPKM)[full.loc]))
#ProLen = as.factor(round(log(data$Pro_len[full.loc])))
#GeneLen = as.factor(round(log(data$Gene_len[full.loc])))
#Iden = as.factor(round(data$Identity[full.loc]/10)*10)

null = lm(data$Stable_run[full.loc] ~ 1)
full = lm(data$Stable_run[full.loc] ~ fpkm + ProLen + GeneLen + Iden
				+ fpkm*ProLen + fpkm*Iden + ProLen*Iden
				+ fpkm*ProLen*Iden)
step(null, scope = list(lower=null,upper=full), direction="both", k = log(sum(full.loc)))
model = lm(data$Stable_run[full.loc] ~ fpkm + ProLen + Iden + ProLen*Iden)
summary(model)
anova(model)

new = data.frame(fpkm=2,ProLen=3,Iden=90)
predict(model, new)

############################## compare fully assembled with partially assembled
############ expression
fpkm = log(data$FPKM)
figureFile = paste(tgtfolder, "/", "fullvspart.fpkm.boxplot.emf", sep="")
	emf(file = figureFile, width = 6, height = 6, family="Helvetica")
	boxplot(fpkm[full.loc], fpkm[part.loc], names=c("Fully Assembled", "Partially Assembled"), ylab="Log FPKM", cex.lab=1.3)
dev.off()
figureFile = paste(tgtfolder, "/", "fullvspart.fpkm.distribution.emf", sep="")
	emf(file = figureFile, width = 6, height = 6, family="Helvetica")
	fpkm.full = density(fpkm[full.loc], adjust=2)
	fpkm.part = density(fpkm[part.loc], adjust=2)
	plot(x=fpkm[full.loc], y=fpkm[full.loc], ylim=c(0,max(c(fpkm.full$y, fpkm.part$y))), type="n", xlab="Log FPKM", ylab="Density", cex.lab=1.3)
	lines(fpkm.full, col="black", lty="solid", lwd=3)
	lines(fpkm.part, col="darkgrey", lty="dashed", lwd=3)
	legend("topleft", legend=c("Fully Assembled", "Partially Assembled"), col=c("black", "darkgrey"), lty=c("solid", "dashed"), bty="n", lwd=3)
dev.off()
	## statistic test
	t.test(fpkm ~ full.loc, alternative="less")
	wilcox.test(fpkm ~ full.loc, alternative="less")
############ interaction
FpkmLen = log(data$Pro_len*data$FPKM)
figureFile = paste(tgtfolder, "/", "fullvspart.prolen.boxplot.emf", sep="")
	emf(file = figureFile, width = 6, height = 6, family="Helvetica")
	boxplot(FpkmLen[full.loc], FpkmLen[part.loc], names=c("Fully Assembled", "Partially Assembled"), ylab="Log FPKM*Length", cex.lab=1.3)
dev.off()
figureFile = paste(tgtfolder, "/", "fullvspart.prolen.distribution.emf", sep="")
	emf(file = figureFile, width = 6, height = 6, family="Helvetica")
	FpkmLen.full = density(FpkmLen[full.loc], adjust=2)
	FpkmLen.part = density(FpkmLen[part.loc], adjust=2)
	plot(x=FpkmLen[full.loc], y=FpkmLen[full.loc], ylim=c(0,max(c(FpkmLen.full$y, FpkmLen.part$y))), type="n", xlab="Log FPKM*Length", ylab="Density", cex.lab=1.3)
	lines(FpkmLen.full, col="black", lty="solid", lwd=3)
	lines(FpkmLen.part, col="darkgrey", lty="dashed", lwd=3)
	legend("topleft", legend=c("Fully Assembled", "Partially Assembled"), col=c("black", "darkgrey"), lty=c("solid", "dashed"), bty="n", lwd=3)
dev.off()
	## statistic test
	t.test(FpkmLen ~ full.loc, alternative="greater")
	wilcox.test(FpkmLen ~ full.loc, alternative="greater")
############ identity
Iden = data$Identity
figureFile = paste(tgtfolder, "/", "fullvspart.prolen.boxplot.emf", sep="")
	emf(file = figureFile, width = 6, height = 6, family="Helvetica")
	boxplot(Iden[full.loc], Iden[part.loc], names=c("Fully Assembled", "Partially Assembled"), ylab="Identity with P. trichocarpa", cex.lab=1.3)
dev.off()
figureFile = paste(tgtfolder, "/", "fullvspart.prolen.distribution.emf", sep="")
	emf(file = figureFile, width = 6, height = 6, family="Helvetica")
	Iden.full = density(Iden[full.loc], adjust=2)
	Iden.part = density(Iden[part.loc], adjust=2)
	plot(x=Iden[full.loc], y=Iden[full.loc], ylim=c(0,max(c(Iden.full$y, Iden.part$y))), type="n", xlab="Identity", ylab="Density", cex.lab=1.3)
	lines(Iden.full, col="black", lty="solid", lwd=3)
	lines(Iden.part, col="darkgrey", lty="dashed", lwd=3)
	legend("topleft", legend=c("Fully Assembled", "Partially Assembled"), col=c("black", "darkgrey"), lty=c("solid", "dashed"), bty="n", lwd=3)
dev.off()
	## statistic test
	t.test(Iden ~ full.loc, alternative="less")
	wilcox.test(Iden ~ full.loc, alternative="less")
	






