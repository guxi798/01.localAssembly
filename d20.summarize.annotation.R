ortho.data = readLines("01.data/02.Transcript/Results_Dec15/WorkingDirectory/Orthogroups.txt")

## read in annotation file
# ath
ath.annot = read.csv("../01.Ptrichocarpa/01.data/00.PriorData/02.TAIR10_functional_descriptions",sep="\t", as.is=TRUE)
ath.annot = ath.annot[,-2]
ath.annot[,1] = gsub("\\.[0-9]*", "", ath.annot[,1])
# app
app.annot = read.csv("../15.Appalachian/16.annotation/annotation.txt", sep="\t", as.is=TRUE, header=FALSE)
app.annot[,1] = paste("app", app.annot[,1], sep=".")
# crp
crp.annot = read.csv("../21.Crepemyrtle/16.annotation/annotation.txt", sep="\t", as.is=TRUE, header=FALSE)
crp.annot[,1] = paste("crp", crp.annot[,1], sep=".")
# oak
oak.annot = read.csv("../07.Oak/16.annotation/annotation.txt", sep="\t", as.is=TRUE, header=FALSE)
oak.annot[,1] = paste("oak", oak.annot[,1], sep=".")
# pch
pch.annot = read.csv("../05.Peach/16.annotation/annotation.txt", sep="\t", as.is=TRUE, header=FALSE)
pch.annot[,1] = paste("pch", pch.annot[,1], sep=".")
# plm
plm.annot = read.csv("../20.Plum/16.annotation/annotation.txt", sep="\t", as.is=TRUE, header=FALSE)
plm.annot[,1] = paste("plm", plm.annot[,1], sep=".")
# ptr
ptr.annot = read.csv("../01.Ptrichocarpa/16.annotation/annotation.txt", sep="\t", as.is=TRUE, header=FALSE)
ptr.annot[,1] = paste("ptr", ptr.annot[,1], sep=".")
# wil 
wil.annot = read.csv("../22.Willow/16.annotation/annotation.txt", sep="\t", as.is=TRUE, header=FALSE)
wil.annot[,1] = paste("wil", wil.annot[,1], sep=".")

## read in expression file
# ath
ath.expr = read.table("../16.Arabidopsis/13.abundance/genome/genes.fpkm_tracking", header=TRUE, as.is=FALSE)
ath.expr = ath.expr[,-c(3,4)]
# app
app.expr = read.table("../15.Appalachian/13.abundance/genome/RSEM/gene.fpkm.rsem.txt", header=TRUE, as.is=FALSE)
app.expr = app.expr[,-c(2,(ncol(app.expr)/2+2):ncol(app.expr))]
app.expr$transcript_id = paste("app", app.expr$transcript_id, sep=".")
# crp
crp.expr = read.table("../21.Crepemyrtle/13.abundance/genome/RSEM/gene.fpkm.rsem.txt", header=TRUE, as.is=FALSE)
crp.expr = crp.expr[,-c(2,(ncol(crp.expr)/2+2):ncol(crp.expr))]
crp.expr$transcript_id = paste("crp", crp.expr$transcript_id, sep=".")
# oak
oak.expr = read.table("../07.Oak/13.abundance/genome/RSEM/gene.fpkm.rsem.txt", header=TRUE, as.is=FALSE)
oak.expr = oak.expr[,-c(2,(ncol(oak.expr)/2+2):ncol(oak.expr))]
oak.expr$transcript_id = paste("oak", oak.expr$transcript_id, sep=".")
# pch
pch.expr = read.table("../05.Peach/13.abundance/genome/RSEM/gene.fpkm.rsem.txt", header=TRUE, as.is=FALSE)
pch.expr = pch.expr[,-c(2,(ncol(pch.expr)/2+2):ncol(pch.expr))]
pch.expr$transcript_id = paste("pch", pch.expr$transcript_id, sep=".")
# plm
plm.expr = read.table("../20.Plum/13.abundance/genome/RSEM/gene.fpkm.rsem.txt", header=TRUE, as.is=FALSE)
plm.expr = plm.expr[,-c(2,(ncol(plm.expr)/2+2):ncol(plm.expr))]
plm.expr$transcript_id = paste("plm", plm.expr$transcript_id, sep=".")
# ptr
ptr.expr = read.table("../01.Ptrichocarpa/13.abundance/genome/RSEM/gene.fpkm.rsem.txt", header=TRUE, as.is=FALSE)
ptr.expr = ptr.expr[,-c(2,(ncol(ptr.expr)/2+2):ncol(ptr.expr))]
ptr.expr$transcript_id = paste("ptr", ptr.expr$transcript_id, sep=".")
# wil
wil.expr = read.table("../22.Willow/13.abundance/genome/RSEM/gene.fpkm.rsem.txt", header=TRUE, as.is=FALSE)
wil.expr = wil.expr[,-c(2,(ncol(wil.expr)/2+2):ncol(wil.expr))]
wil.expr$transcript_id = paste("wil", wil.expr$transcript_id, sep=".")

x = strsplit(ortho.data, "\\s+|:\\s+")
n0 = 4
n1 = 1
n2 = 1
n3 = ncol(app.annot) - 2
n4 = ncol(crp.annot) - 2
n5 = ncol(oak.annot) - 2
n6 = ncol(pch.annot) - 3
n7 = ncol(plm.annot) - 3
n8 = ncol(ptr.annot) - 2
n9 = ncol(wil.annot) - 2

ortho = NULL

m = lapply(x, summarizeAnnotation)
#ortho = do.call(rbind, m)

summarizeAnnotation <- function(x){
  #for(i in 1:100){
  #if(i %% 100 == 0){print(i)}
  temp = x
  group = temp[1]
  sp0 = temp[grep("^AT", temp)]
  sp1 = temp[grep("^Glyma", temp)]
  sp2 = temp[grep("^Migut", temp)]
  sp3 = temp[grep("^app", temp)]
  sp4 = temp[grep("^crp", temp)]
  sp5 = temp[grep("^oak", temp)]
  sp6 = temp[grep("^pch", temp)]
  sp7 = temp[grep("^plm", temp)]
  sp8 = temp[grep("^ptr", temp)]
  sp9 = temp[grep("^wil", temp)]
  
  l0 = length(sp0)
  l1 = length(sp1)
  l2 = length(sp2)
  l3 = length(sp3)
  l4 = length(sp4)
  l5 = length(sp5)
  l6 = length(sp6)
  l7 = length(sp7)
  l8 = length(sp8)
  l9 = length(sp9)
  l = max(c(l0, l1, l2, l3, l4, l5, l6, l7, l8, l9))
  
  if(l0 > 0){
    loc = match(sp0, ath.annot[,1])
    sp0.annot = cbind(sp0, ath.annot[loc, -1], stringsAsFactors=FALSE)
  }
  if(l3 > 0){
    loc = match(sp3, app.annot[,1])
    sp3.annot = cbind(sp3, app.annot[loc, -c(1,3,8)], stringsAsFactors=FALSE)
  }
  if(l4 > 0){
    loc = match(sp4, crp.annot[,1])
    sp4.annot = cbind(sp4, crp.annot[loc, -c(1,3,8)], stringsAsFactors=FALSE)
  }
  if(l5 > 0){
    loc = match(sp5, oak.annot[,1])
    sp5.annot = cbind(sp5, oak.annot[loc, -c(1,3,8)], stringsAsFactors=FALSE)
  }
  if(l6 > 0){
    loc = match(sp6, pch.annot[,1])
    sp6.annot = cbind(sp6, pch.annot[loc, -c(1,3,8,13)], stringsAsFactors=FALSE)
  }
  if(l7 > 0){
    loc = match(sp7, plm.annot[,1])
    sp7.annot = cbind(sp7, plm.annot[loc, -c(1,3,8,13)], stringsAsFactors=FALSE)
  }
  if(l8 > 0){
    loc = match(sp8, ptr.annot[,1])
    sp8.annot = cbind(sp8, ptr.annot[loc, -c(1,3,8)], stringsAsFactors=FALSE)
  }
  if(l9 > 0){
    loc = match(sp9, wil.annot[,1])
    sp9.annot = cbind(sp9, wil.annot[loc, -c(1,3,8)], stringsAsFactors=FALSE)
  }
  
  g = NULL
  for(j in 1:l){
    if(l0 == 0){g0 = rep(".",n0)}else if(j <= l0){g0 = sp0.annot[j,]}else{g0 = rep(".",n0)}
    if(l1 == 0){g1 = rep(".",n1)}else if(j <= l1){g1 = sp1[j]}else{g1 = rep(".",n1)}
    if(l2 == 0){g2 = rep(".",n2)}else if(j <= l2){g2 = sp2[j]}else{g2 = rep(".",n2)}
    if(l3 == 0){g3 = rep(".",n3)}else if(j <= l3){g3 = sp3.annot[j,]}else{g3 = rep(".",n3)}
    if(l4 == 0){g4 = rep(".",n4)}else if(j <= l4){g4 = sp4.annot[j,]}else{g4 = rep(".",n4)}
    if(l5 == 0){g5 = rep(".",n5)}else if(j <= l5){g5 = sp5.annot[j,]}else{g5 = rep(".",n5)}
    if(l6 == 0){g6 = rep(".",n6)}else if(j <= l6){g6 = sp6.annot[j,]}else{g6 = rep(".",n6)}
    if(l7 == 0){g7 = rep(".",n7)}else if(j <= l7){g7 = sp7.annot[j,]}else{g7 = rep(".",n7)}
    if(l8 == 0){g8 = rep(".",n8)}else if(j <= l8){g8 = sp8.annot[j,]}else{g8 = rep(".",n8)}
    if(l9 == 0){g9 = rep(".",n9)}else if(j <= l9){g9 = sp9.annot[j,]}else{g9 = rep(".",n8)}
    g = rbind(g, c(g0, g1, g2, g3, g4, g5, g6, g7, g8, g9))
  }
  g = data.frame(rep(group, l), g, stringsAsFactors = FALSE)
  names(g) = c("ortho_group","Ath_gene","Short_description","Curator_summary","Computational_description","Gma_gene","Mgu_gene",
			"Dogwood_gene","Best_hit_ptr_name","Ptr2Ath","Ptr_ShortName","Ptr_Description","Best_hit_ath_name","Ath_Description","UniProt_Description",
			"CrepeMyrtle_gene","Best_hit_ptr_name","Ptr2Ath","Ptr_ShortName","Ptr_Description","Best_hit_ath_name","Ath_Description","UniProt_Description",
			"Oak_gene","Best_hit_ptr_name","Ptr2Ath","Ptr_ShortName","Ptr_Description","Best_hit_ath_name","Ath_Description","UniProt_Description",
			"Peach_gene","Best_hit_peach_name","Peach2Ath","Peach_ShortName","Peach_Description","Best_hit_ptr_name","Ptr2Ath","Ptr_ShortName","Ptr_Description","Best_hit_ath_name","Ath_Description","UniProt_Description",
			"Plum_gene","Best_hit_peach_name","Peach2Ath","Peach_ShortName","Peach_Description","Best_hit_ptr_name","Ptr2Ath","Ptr_ShortName","Ptr_Description","Best_hit_ath_name","Ath_Description","UniProt_Description",
			"Ptr_gene","Best_hit_ptr_name","Ptr2Ath","Ptr_ShortName","Ptr_Description","Best_hit_ath_name","Ath_Description","UniProt_Description",
			"Willow_gene","Best_hit_ptr_name","Ptr2Ath","Ptr_ShortName","Ptr_Description","Best_hit_ath_name","Ath_Description","UniProt_Description")
  return(g)
}


for(i in 1:round(length(m)/100)){
	print(i)
	if(i < round(length(m)/100)){
		ortho = do.call(rbind, m[((i-1)*100+1):(i*100)])
		ortho = apply(ortho, 2, as.character)
		if(i == 1){
			write.table(ortho, file="01.data/02.Transcript/annotation/ortholog_annotation.txt", row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t", append=FALSE)
		}else{
			write.table(ortho, file="01.data/02.Transcript/annotation/ortholog_annotation.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", append=TRUE)
		}
	}else{
		ortho = do.call(rbind, m[((i-1)*100+1):length(m)])
		ortho = apply(ortho, 2, as.character)
		write.table(ortho, file="01.data/02.Transcript/annotation/ortholog_annotation.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t", append=TRUE)
	}
}

save(m, file="01.data/02.Transcript/annotation/ortholog_annotation.RData")
#write.table(ortho, "01.data/02.Transcript/annotation/ortholog_annotation.txt", row.names=FALSE, quote=FALSE, sep="\t")

ortho.data = read.csv("01.data/02.Transcript/annotation/ortholog_annotation.txt", as.is=TRUE, header=TRUE, sep="\t")
ath.loc = match(ortho.data$Ath_gene, ath.expr$transcript_id)
pos = which(names(ortho.data) == "Ath_gene")
ortho.data = data.frame(ortho.data[1:pos], ath.expr[ath.loc, -1], ortho.data[(pos+1):ncol(ortho.data)], stringsAsFactors = FALSE)
names(ortho.data)[pos+1] = "ath_pollen"

app.loc = match(ortho.data$Dogwood_gene, app.expr$transcript_id)
pos = which(names(ortho.data) == "Dogwood_gene")
ortho.data = data.frame(ortho.data[1:pos], app.expr[app.loc, -1], ortho.data[(pos+1):ncol(ortho.data)], stringsAsFactors = FALSE)

crp.loc = match(ortho.data$CrepeMyrtle_gene, crp.expr$transcript_id)
pos = which(names(ortho.data) == "CrepeMyrtle_gene")
ortho.data = data.frame(ortho.data[1:pos], crp.expr[crp.loc, -1], ortho.data[(pos+1):ncol(ortho.data)], stringsAsFactors = FALSE)

oak.loc = match(ortho.data$Oak_gene, oak.expr$transcript_id)
pos = which(names(ortho.data) == "Oak_gene")
ortho.data = data.frame(ortho.data[1:pos], oak.expr[oak.loc, -1], ortho.data[(pos+1):ncol(ortho.data)], stringsAsFactors = FALSE)

pch.loc = match(ortho.data$Peach_gene, pch.expr$transcript_id)
pos = which(names(ortho.data) == "Peach_gene")
ortho.data = data.frame(ortho.data[1:pos], pch.expr[pch.loc, -1], ortho.data[(pos+1):ncol(ortho.data)], stringsAsFactors = FALSE)

plm.loc = match(ortho.data$Plum_gene, plm.expr$transcript_id)
pos = which(names(ortho.data) == "Plum_gene")
ortho.data = data.frame(ortho.data[1:pos], plm.expr[plm.loc, -1], ortho.data[(pos+1):ncol(ortho.data)], stringsAsFactors = FALSE)

ptr.loc = match(ortho.data$Ptr_gene, ptr.expr$transcript_id)
pos = which(names(ortho.data) == "Ptr_gene")
ortho.data = data.frame(ortho.data[1:pos], ptr.expr[ptr.loc, -1], ortho.data[(pos+1):ncol(ortho.data)], stringsAsFactors = FALSE)

wil.loc = match(ortho.data$Willow_gene, wil.expr$transcript_id)
pos = which(names(ortho.data) == "Willow_gene")
ortho.data = data.frame(ortho.data[1:pos], wil.expr[wil.loc, -1], ortho.data[(pos+1):ncol(ortho.data)], stringsAsFactors = FALSE)

ortho.data[is.na(ortho.data)] = "."

write.table(ortho.data, "01.data/02.Transcript/annotation/ortholog_annotation_expr.txt", row.names=FALSE, quote=FALSE, sep="\t")


ath_og = ortho.data$ortho_group[ortho.data$ath_pollen >= 5 & !is.na(ortho.data$ath_pollen)]
length(unique(ath_og))
app_fpkm = apply(ortho.data[,10:12], 1, mean)
app_og = ortho.data$ortho_group[ortho.data$Dogwood_gene != "." & app_fpkm >= 5 & !is.na(app_fpkm)]
length(unique(app_og))
ptr_fpkm = apply(ortho.data[,70:72], 1, mean)
ptr_og = ortho.data$ortho_group[ortho.data$Ptr_gene != "." & ptr_fpkm >= 5 & !is.na(ptr_fpkm)]
length(unique(ptr_og))
plm_og = ortho.data$ortho_group[ortho.data$Plum_gene != "."]
length(unique(plm_og))

length(intersect(unique(ath_og), unique(app_og)))
length(intersect(unique(ath_og), unique(ptr_og)))
length(intersect(unique(ptr_og), unique(app_og)))
length(intersect(intersect(unique(ath_og), unique(ptr_og)), unique(app_og)))

length(intersect(unique(ath_og), unique(plm_og)))
length(intersect(unique(ptr_og), unique(plm_og)))
length(intersect(unique(app_og), unique(plm_og)))
