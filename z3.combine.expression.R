annot = read.csv("14.annotation/annotation.txt", header=FALSE, sep="\t")
rsem = read.table("13.abundance/genome/RSEM/gene.fpkm.rsem.txt", header=TRUE)
expr = read.table("13.abundance/genome/eXpress/gene.fpkm.express.txt", header=TRUE)

colnames(annot) = c("transcript_id", "ptr_id", "ptr_GO", "ptr2ath_abbr", "ptr2ath_desc", "ath_id", "ath_GO", "ath2osa_desc", "uniprot")
pos1 = match(annot[,1], rsem$transcript_id)
pos2 = match(annot[,1], expr$transcript_id)

data = cbind(annot, rsem[pos1,-c(1,7:10)], expr[pos2,-1])

out = "13.abundance/gene.fpkm.all.txt"
write.table(data, file=out, quote=FALSE)
