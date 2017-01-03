data = read.table("12.Local.Trinity.comparison/03.CDS/Trinity_Local_intersect_quality.v2.txt", as.is=TRUE, header=FALSE)
set.seed(12345)
s = sample(1:nrow(data))

local_track = NULL
local_mismatch = 0
local_bits = 0
local_num = 0
trinity_track = NULL
trinity_mismatch = 0
trinity_bits = 0
trinity_num = 0

for(i in s){
	local_num = local_num + data[i,2]
	local_mismatch = local_mismatch + data[i,3]
	local_bits = local_bits + data[i,4]
	local_track = rbind(local_track, c(local_mismatch/local_num, local_bits/local_num))
	
	trinity_num = trinity_num + data[i,5]
	trinity_mismatch = trinity_mismatch + data[i,6]
	trinity_bits = trinity_bits + data[i,7]
	trinity_track = rbind(trinity_track, c(trinity_mismatch/trinity_num, trinity_bits/trinity_num))
}

plot(x=local_track[,1], y=trinity_track[,1], xlim=c(0,10), ylim=c(0,10))
abline(a=0, b=1, col="red")

plot(local_track[,1], type="l")
lines(trinity_track[,1], col="red")
local_track[5940,]
trinity_track[5940,]

data = read.csv("12.Local.Trinity.comparison/03.CDS/Trinity_Local_intersect.txt", as.is=TRUE, header=FALSE)

loc = grep("^>", data[,1])
data = data[-loc, ]
data = data.frame(matrix(unlist(strsplit(data, "\\t|\\s+")), ncol=13, byrow=TRUE))

set.seed(12345)
s = sample(1:nrow(data))

data = read.table("12.Local.Trinity.comparison/03.CDS/Trinity_Local_intersect_quality.v2.txt", as.is=TRUE, header=FALSE)
size = 0.7
w = 4
h = 4

emf(file="12.Local.Trinity.comparison/03.CDS/01.mismatch.comparison2.emf", width=w, height=h)
	plot(data[,2]-data[,5], xlab="Gene", ylab="Local - Trinity", main="Mismatch Number", cex=size)
dev.off()
data1 = data
data1[data1==0] = 0.1
emf(file="12.Local.Trinity.comparison/03.CDS/01.mismatch.comparison.emf", width=w, height=h)
	plot(x=log(data1[,2]), y=log(data1[,5]), xlab="Local assembly", ylab="Trinity assembly", main="Mismatch Number", cex=size)
dev.off()

emf(file="12.Local.Trinity.comparison/03.CDS/02.bitscore.comparison.emf", width=w, height=h)
	plot(x=data[,3], y=data[,6], xlab="Local assembly",ylab="Trinity assembly", main="Bitscore", cex=size)
dev.off()
emf(file="12.Local.Trinity.comparison/03.CDS/02.bitscore.comparison1.emf", width=w, height=h)
	plot(x=data[,3], y=data[,6], xlim=c(0,7000), ylim=c(0,7000), xlab="Local assembly",ylab="Trinity assembly", main="Bitscore", cex=size)
dev.off()

emf(file="12.Local.Trinity.comparison/03.CDS/02.bitscore.comparison2.emf", width=w, height=h)
	plot(data[,3] - data[,6], xlab="Gene", ylab="Local - Trinity", main="Bitscore", cex=size)
dev.off()
emf(file="12.Local.Trinity.comparison/03.CDS/02.bitscore.comparison3.emf", width=w, height=h)
	plot(data[,3] - data[,6], ylim=c(-1000, 1000), xlab="Gene", ylab="Local - Trinity", main="Bitscore", cex=size)
dev.off()

emf(file="12.Local.Trinity.comparison/03.CDS/03.length.comparison.emf", width=w, height=h)
	plot(x=data[,4], y=data[,7], xlab="Local assembly",ylab="Trinity assembly", main="Bitscore", cex=size)
dev.off()
emf(file="12.Local.Trinity.comparison/03.CDS/03.length.comparison1.emf", width=w, height=h)
	plot(x=data[,4], y=data[,7], xlim=c(0,4000), ylim=c(0,4000), xlab="Local assembly",ylab="Trinity assembly", main="Bitscore", cex=size)
dev.off()

emf(file="12.Local.Trinity.comparison/03.CDS/03.length.comparison2.emf", width=w, height=h)
	plot(data[,4] - data[,7], xlab="Gene", ylab="Local - Trinity", main="Transcript Length", cex=size)
dev.off()

blast.l = read.table("01.data/06.TargetTranscriptome/transcriptome.v1.ptr.blastn.cds.out", header=FALSE)
blast.t = read.table("05.Trinity/Trinity.new.ptr.blastn.cds.out", header=FALSE)
summary(blast.l[,4])
summary(blast.t[,4])

summary(blast.l[,5])
summary(blast.t[,5])

summary(blast.l[,5]/blast.l[,4])
summary(blast.t[,5]/blast.t[,4])

################################################################
## compare with public available data
public_data = read.table("14.evaluation/04.TransRate/public_data_set.txt", header=TRUE)
sum(public_data$score <= 0.35)/nrow(public_data)
sum(public_data$score <= 0.2833)/nrow(public_data)

data.l = read.csv("14.evaluation/04.TransRate/05.Local/transcriptome.v1/contigs.csv")
data.t = read.csv("14.evaluation/04.TransRate/06.Trinity/Trinity.new/contigs.csv")

contig_score = c(data.l$score, data.t$score)
method = factor(c(rep("Local", nrow(data.l)), rep("Trinity", nrow(data.t))))
contig_score_dfs = data.frame(Score = contig_score, Method = method)

emf(file="14.evaluation/04.TransRate/03.Local/transcriptome.v1/contig_score.emf", width=3.5, height=1.5)
ggplot(contig_score_dfs, aes(x=Score)) + geom_density(aes(group=Method, color=Method)) + xlab("Contig Score")
dev.off()

ref_cov = c(data.l$reference_coverage, data.t$reference_coverage)
method = factor(c(rep("Local", nrow(data.l)), rep("Trinity", nrow(data.t))))
ref_cov_dfs = data.frame(Coverage = ref_cov, Method = method)

emf(file="14.evaluation/04.TransRate/03.Local/transcriptome.v1/reference_coverage_cds.emf", width=3.5, height=1.5)
ggplot(ref_cov_dfs, aes(x=Coverage)) + geom_density(aes(group=Method, color=Method)) + xlab("Reference Coverage")#+ coord_cartesian(xlim = c(0.5,1))
dev.off()

p_good = c(data.l$p_good, data.t$p_good)
method = factor(c(rep("Local", nrow(data.l)), rep("Trinity", nrow(data.t))))
p_good_dfs = data.frame(Coverage = p_good, Method = method)

emf(file="14.evaluation/04.TransRate/03.Local/transcriptome.v1/p_good.emf", width=3.5, height=1.5)
ggplot(p_good_dfs, aes(x=Coverage)) + geom_density(aes(group=Method, color=Method)) + xlab("P Good")
dev.off()

l = c(data.l$length, data.t$length)
method = factor(c(rep("Local", nrow(data.l)), rep("Trinity", nrow(data.t))))
length_dfs = data.frame(Coverage = l, Method = method)

emf(file="14.evaluation/04.TransRate/03.Local/transcriptome.v1/length.emf", width=3.5, height=1.5)
ggplot(length_dfs, aes(x=Coverage)) + geom_density(aes(group=Method, color=Method)) + xlab("Length") + xlim(0,1500)
dev.off()

segment = c(data.l$p_not_segmented, data.t$p_not_segmented)
method = factor(c(rep("Local", nrow(data.l)), rep("Trinity", nrow(data.t))))
segment_dfs = data.frame(Not_Segmented = fragment, Method = method)

emf(file="14.evaluation/04.TransRate/03.Local/transcriptome.v1/chimerism2.emf", width=3.5, height=1.5)
ggplot(segment_dfs, aes(x=Not_Segmented)) + geom_density(aes(group=Method, color=Method)) + xlab("P Not Segmented")+ coord_cartesian(xlim = c(0.5,1))
dev.off()

p_seq_true = c(data.l$p_seq_true, data.t$p_seq_true)
method = factor(c(rep("Local", nrow(data.l)), rep("Trinity", nrow(data.t))))
p_seq_true_dfs = data.frame(p_seq_true = p_seq_true, Method = method)

emf(file="14.evaluation/04.TransRate/03.Local/transcriptome.v1/gene_collapse2.emf", width=3.5, height=1.5)
ggplot(p_seq_true_dfs, aes(x=p_seq_true)) + geom_density(aes(group=Method, color=Method)) + xlab("P Seq True")+ coord_cartesian(xlim = c(0.8,1))
dev.off()

p_bases_covered = c(data.l$p_bases_covered, data.t$p_bases_covered)
method = factor(c(rep("Local", nrow(data.l)), rep("Trinity", nrow(data.t))))
p_bases_covered_dfs = data.frame(p_bases_covered = p_bases_covered, Method = method)

emf(file="14.evaluation/04.TransRate/03.Local/transcriptome.v1/p_bases_covered2.emf", width=3.5, height=1.5)
ggplot(p_bases_covered_dfs, aes(x=p_bases_covered)) + geom_density(aes(group=Method, color=Method)) + xlab("P Bases Covered") + coord_cartesian(xlim = c(0.95,1))
dev.off()

## compare features of unique and common genes
expr = read.table("01.data/12.deseq/genes.fpkm_tracking")
expr = apply(expr, 1, mean)
g1 = expr[expr >=0 & expr <1]
g2 = expr[expr >=1 & expr <10]
g3 = expr[expr >=10 & expr <50]
g4 = expr[expr >=50 & expr <100]
g5 = expr[expr >=100 & expr <500]
g6 = expr[expr >=500 & expr <1000]
g7 = expr[expr >=1000]


data.l = read.csv("12.Local.Trinity.comparison/03.CDS/Local_uniq.txt", as.is=TRUE, header=FALSE)
loc = grep("^>", data.l[,1])
data.l = data.l[loc, ]
data.l = data.frame(matrix(unlist(strsplit(data.l, "\\t|\\s+")), ncol=2, byrow=TRUE), stringsAsFactors=FALSE)
data.l[,1] = gsub(">", "", data.l[,1])
loc = match(data.l[,1], names(expr))
data.l.expr = expr[loc]
data.l.expr[data.l.expr == 0] = 0.01
data.l.expr = log(data.l.expr)

data.t = read.csv("12.Local.Trinity.comparison/03.CDS/Trinity_uniq.txt", as.is=TRUE, header=FALSE)
loc = grep("^>", data.t[,1])
data.t = data.t[loc, ]
data.t = data.frame(matrix(unlist(strsplit(data.t, "\\t|\\s+")), ncol=2, byrow=TRUE), stringsAsFactors=FALSE)
data.t[,1] = gsub(">", "", data.t[,1])
loc = match(data.t[,1], names(expr))
data.t.expr = expr[loc]
data.t.expr[data.t.expr == 0] = 0.01
data.t.expr = log(data.t.expr)

data.i = read.csv("12.Local.Trinity.comparison/03.CDS/Trinity_Local_intersect.txt", as.is=TRUE, header=FALSE)
loc = grep("^>", data.i[,1])
data.i = data.i[loc, ]
data.i = data.frame(matrix(unlist(strsplit(data.i, "\\t|\\s+")), ncol=2, byrow=TRUE), stringsAsFactors=FALSE)
data.i[,1] = gsub(">", "", data.i[,1])
loc = match(data.i[,1], names(expr))
data.i.expr = expr[loc]
data.i.expr[data.i.expr == 0] = 0.01
data.i.expr = log(data.i.expr)

emf(file="12.Local.Trinity.comparison/03.CDS/05.expression.trinity_intersect_local.emf", width=w, height=h)
	boxplot(data.t.expr, data.i.expr, data.l.expr, names=c("Trinity Unique", "Intersect", "Local Unique"))
dev.off()
emf(file="12.Local.Trinity.comparison/03.CDS/05.expression.trinity_local_uniq.emf", width=w, height=h)
	boxplot(data.t.expr, data.l.expr, names=c("Trinity Unique", "PLAS Unique"))
dev.off()
emf(file="12.Local.Trinity.comparison/03.CDS/05.expression.trinity_local.emf", width=w, height=h)
	boxplot(c(data.t.expr, data.i.expr), c(data.i.expr,data.l.expr), names=c("Trinity", "Local"))
dev.off()

len.t = log(as.numeric(data.t[,2]))
len.i = log(as.numeric(data.i[,2]))
len.l = log(as.numeric(data.l[,2]))
emf(file="12.Local.Trinity.comparison/03.CDS/06.length.trinity_intersect_local.emf", width=w, height=h)
	boxplot(len.t, len.i, len.l, names=c("Trinity Unique", "Intersect", "Local Unique"))
dev.off()
emf(file="12.Local.Trinity.comparison/03.CDS/06.length.trinity_local_uniq.emf", width=w, height=h)
	boxplot(len.t, len.l, names=c("Trinity Unique", "PLAS Unique"))
dev.off()
emf(file="12.Local.Trinity.comparison/03.CDS/06.length.trinity_local.emf", width=w, height=h)
	boxplot(c(len.t, len.i), c(len.i,len.l), names=c("Trinity", "Local"))
dev.off()

l = unique(c(data.l[,1], data.i[,1]))
t = unique(c(data.t[,1], data.i[,1]))
rec_rate.l = c(sum(l %in% names(g1)), sum(l %in% names(g2)), sum(l %in% names(g3)), sum(l %in% names(g4)),
			sum(l %in% names(g5)), sum(l %in% names(g6)), sum(l %in% names(g7)))
rec_rate.t = c(sum(t %in% names(g1)), sum(t %in% names(g2)), sum(t %in% names(g3)), sum(t %in% names(g4)),
			sum(t %in% names(g5)), sum(t %in% names(g6)), sum(t %in% names(g7)))
rec_rate.l = rec_rate.l / c(length(g1),length(g2),length(g3),length(g4),length(g5),length(g6),length(g7))
rec_rate.t = rec_rate.t / c(length(g1),length(g2),length(g3),length(g4),length(g5),length(g6),length(g7))

lng1 = names(g1)[!(names(g1) %in% l)]
lng2 = names(g2)[!(names(g2) %in% l)]
lng3 = names(g3)[!(names(g3) %in% l)]
lng4 = names(g4)[!(names(g4) %in% l)]
lng5 = names(g5)[!(names(g5) %in% l)]
lng6 = names(g6)[!(names(g6) %in% l)]
lng7 = names(g7)[!(names(g7) %in% l)]

emf(file="12.Local.Trinity.comparison/03.CDS/04.recovery_rate.emf", width=w, height=h)
plot(rec_rate.l, type="b", col="red", xlab="FPKM Range", ylab="Recovery Rate", xaxt="n", cex=size)
lines(rec_rate.t, type="b", col="blue")
legend("bottomright", legend=c("PLAS", "Trinity"), col=c("red", "blue"), lty=1, pch=1, bty="n")
text(1:8, par("usr")[3]-0.035, labels=c("0-1","1-10","10-50","50-100","100-500", "500-1000",">1000"), srt=45, pos = 1, xpd = TRUE)
dev.off()

copy = read.table("12.Local.Trinity.comparison/03.CDS/copy.number.txt", as.is=TRUE)
c1 = copy[copy[,2]==1,1]
c2 = copy[copy[,2]==2,1]
c3 = copy[copy[,2]>=3 & copy[,2] <5,1]
c4 = copy[copy[,2]>=5 & copy[,2] <10,1]
c5 = copy[copy[,2]>=10 & copy[,2] <50,1]
c6 = copy[copy[,2]>=50,1]

l = unique(c(data.l[,1], data.i[,1]))
t = unique(c(data.t[,1], data.i[,1]))
rec_rate.l = c(sum(l %in% c1), sum(l %in% c2), sum(l %in% c3), sum(l %in% c4),
			sum(l %in% c5), sum(l %in% c6))
rec_rate.t = c(sum(t %in% c1), sum(t %in% c2), sum(t %in% c3), sum(t %in% c4),
			sum(t %in% c5), sum(t %in% c6))
rec_rate.l = rec_rate.l / c(length(c1),length(c2),length(c3),length(c4),length(c5),length(c6))
rec_rate.t = rec_rate.t / c(length(c1),length(c2),length(c3),length(c4),length(c5),length(c6))

emf(file="12.Local.Trinity.comparison/03.CDS/04.recovery_rate_copynumber1.emf", width=w, height=h)
	plot(rec_rate.l, type="b", col="red", xlab="Copy Number", ylab="Recovery Rate", xaxt="n", cex=size)
	lines(rec_rate.t, type="b", col="blue")
	legend("topright", legend=c("PLAS", "Trinity"), col=c("red", "blue"), lty=1, pch=1, bty="n")
	text(1:8, par("usr")[3]-0.006, labels=c("1","2","3-5","5-10","10-50", ">50"), srt=45, pos = 1, xpd = TRUE)
dev.off()

l = unique(data.l[,1])
t = unique(data.t[,1])
i = unique(data.i[,1])
rec_rate.l = c(sum(l %in% c1), sum(l %in% c2), sum(l %in% c3), sum(l %in% c4), sum(l %in% c5), sum(l %in% c6))
rec_rate.t = c(sum(t %in% c1), sum(t %in% c2), sum(t %in% c3), sum(t %in% c4), sum(t %in% c5), sum(t %in% c6))
rec_rate.i = c(sum(i %in% c1), sum(i %in% c2), sum(i %in% c3), sum(i %in% c4), sum(i %in% c5), sum(i %in% c6))
rec_rate.l = rec_rate.l / c(length(c1),length(c2),length(c3),length(c4),length(c5),length(c6))
rec_rate.t = rec_rate.t / c(length(c1),length(c2),length(c3),length(c4),length(c5),length(c6))
rec_rate.i = rec_rate.i / c(length(c1),length(c2),length(c3),length(c4),length(c5),length(c6))

emf(file="12.Local.Trinity.comparison/03.CDS/04.recovery_rate_copynumber2.emf", width=w, height=h)
	plot(rec_rate.i, type="b", col="black", xlab="Copy Number", ylab="Recovery Rate", xaxt="n", cex=size)
	lines(rec_rate.t, type="b", col="blue")
	lines(rec_rate.l, type="b", col="red")
	legend("topright", legend=c("PLAS Unique", "Trinity Unique", "Common Set"), col=c("red", "blue", "black"), lty=1, pch=1, bty="n")
	text(1:8, par("usr")[3]-0.006, labels=c("1","2","3-5","5-10","10-50", ">50"), srt=45, pos = 1, xpd = TRUE)
dev.off()

source("FASTKs/bin/run_mclust.r")
library(mclust)
pdf("ks_distn.pdf")
#runmclust("ptr_ptr_ks/ptr_ptr.kaks.txt",4, "ptr_PLAS")
runmclust("app_app_ks/app_app.kaks.txt",4, "app_PLAS")
dev.off()
