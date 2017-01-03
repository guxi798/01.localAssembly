#!/usr/bin/perl -w
# run the script: time 00.script/c13.compare.Trinity.local.pl 05.Trinity/count1 08.full.length/Final.fasta 10.unmapped.reads.trinity/full.length.contigs.fasta

trinity_data = read.table("trinity_full_length_info.txt", sep=" ")
local_data = read.table("local_full_length_info.txt", sep=" ")
trinity_blast = read.table("trinity_blast", sep="\t")
local_blast = read.table("local_blast", sep="\t")

loc = trinity_data[,3] %in% local_data[,3]
trinity_uniq = trinity_data[!loc,]
loc = local_data[,3] %in% trinity_data[,3]
local_uniq = local_data[!loc,]

## find blast info for trinity unique gene
trinity_blast_subset = trinity_blast[trinity_blast[,2] %in% trinity_uniq[,3],]
local_blast_subset = local_blast[local_blast[,2] %in% trinity_uniq[,3],]
trinity_blast_subset = cbind("Trinity", trinity_blast_subset)
local_blast_subset = cbind("Local", local_blast_subset)
blast_subset = rbind(trinity_blast_subset, local_blast_subset)

write.table(blast_subset, "trinity_uniq_blast_info.txt", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

## find blast info for local unique gene
trinity_blast_subset = trinity_blast[trinity_blast[,2] %in% local_uniq[,3],]
local_blast_subset = local_blast[local_blast[,2] %in% local_uniq[,3],]
trinity_blast_subset = cbind("Trinity", trinity_blast_subset)
local_blast_subset = cbind("Local", local_blast_subset)
blast_subset = rbind(trinity_blast_subset, local_blast_subset)

write.table(blast_subset, "local_uniq_blast_info.txt", col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")
