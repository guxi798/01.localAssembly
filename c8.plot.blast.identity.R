library(devEMF)

both = read.table("08.GlobLocal.Comparison/both.full.identity.txt", header=TRUE)
glob = read.table("08.GlobLocal.Comparison/global.uniq.identity.txt", header=TRUE)
loc = read.table("08.GlobLocal.Comparison/local.uniq.identity.txt", header=TRUE)

emf(file = "08.GlobLocal.Comparison/Boxplot.identity.emf", width = 5, height = 5, family="Helvetica")
par(mar=c(5,4,2,2))
m = boxplot(glob$Glob_ident, both$Glob_ident, both$Loc_ident, loc$Loc_ident,
        xaxt="n")
text(1:4, par("usr")[3], labels=c("Glob_Uniq", "Glob_Comm", "Loc_Comm", "Loc_Uniq"),
     srt=45, adj=c(1,1,1,1),xpd=TRUE, cex=1)
dev.off()

emf(file = "08.GlobLocal.Comparison/Boxplot.identityE.emf", width = 5, height = 5, family="Helvetica")
par(mar=c(5,4,2,2))
m = boxplot(glob$Glob_identE, both$Glob_identE, both$Loc_identE, loc$Loc_identE,
            xaxt="n")
text(1:4, par("usr")[3], labels=c("Glob_Uniq", "Glob_Comm", "Loc_Comm", "Loc_Uniq"),
     srt=45, adj=c(1,1,1,1),xpd=TRUE, cex=1)
dev.off()

emf(file = "08.GlobLocal.Comparison/Boxplot.conserved.emf", width = 5, height = 5, family="Helvetica")
par(mar=c(5,4,2,2))
m = boxplot(glob$Glob_cons, both$Glob_cons, both$Loc_cons, loc$Loc_cons,
            xaxt="n")
text(1:4, par("usr")[3], labels=c("Glob_Uniq", "Glob_Comm", "Loc_Comm", "Loc_Uniq"),
     srt=45, adj=c(1,1,1,1),xpd=TRUE, cex=1)
dev.off()

emf(file = "08.GlobLocal.Comparison/Boxplot.conservedE.emf", width = 5, height = 5, family="Helvetica")
par(mar=c(5,4,2,2))
m = boxplot(glob$Glob_consE, both$Glob_consE, both$Loc_consE, loc$Loc_consE,
            xaxt="n")
text(1:4, par("usr")[3], labels=c("Glob_Uniq", "Glob_Comm", "Loc_Comm", "Loc_Uniq"),
     srt=45, adj=c(1,1,1,1),xpd=TRUE, cex=1)
dev.off()
