#!/bin/bash
#PBS -q batch
#PBS -N local_StHe
#PBS -l nodes=1:ppn=1:AMD
#PBS -l walltime=4800:00:00

module load perl/5.20.2
cd $PBS_O_WORKDIR

platform="Sapelo"

####################### this starts with proteome #######################
# WU-BLAST
#/usr/local/wublast/latest/xdformat -p -o 01.data/03.MCL/01.blast/database 01.data/00.PriorData/proteome.fa
#/usr/local/wublast/latest/blastp 01.data/03.MCL/01.blast/database 01.data/00.PriorData/proteome.fa -o 01.data/03.MCL/01.blast/blast.all.out -e 1e-5 -mformat 2 -cpus 4 -wordmask seg

# NCBI-BLAST
#/usr/local/ncbiblast+/latest/bin/makeblastdb -in 01.data/00.PriorData/proteome.fa prot
#time /usr/local/ncbiblast+/latest/bin/blastp -num_threads 4 -db 01.data/00.PriorData/proteome.fa -query 01.data/00.PriorData/proteome.fa -out 01.data/03.MCL/01.blast/blast.all.out -evalue 1e-5 -outfmt 6

# mcl
#time perl 00.script/a1.mcl.prepare.graph.pl 01.data/03.MCL/01.blast/ncbi.blast.all.out 01.data/03.MCL/02.mcl/mcl.graph.txt
#time /usr/local/mcl/latest/bin/mcl 01.data/03.MCL/02.mcl/mcl.graph.txt --abc -o 01.data/03.MCL/02.mcl/mcl.out.txt -I 1.5

#time perl 00.script/a1.mcl.prepare.graph.pl 01.data/03.MCL/01.blast/wu.blast.all.out 01.data/03.MCL/02.mcl/mcl.graph.txt wu
#time /usr/local/mcl/latest/bin/mcl 01.data/03.MCL/02.mcl/wu.mcl.graph.txt --abc -o 01.data/03.MCL/02.mcl/mcl.out.txt -I 1.5

# construct meta-group
#module load R/3.1.2
#time /usr/local/apps/R/3.1.2/bin/Rscript 00.script/a3.geneSelection.R 01.data/03.MCL/02.mcl/mcl.out.txt 01.data/04.GeneOfInterest/GeneID.txt 1000 Potri

# split gene
#time perl 00.script/a4.splitGene.pl 01.data/00.PriorData/proteome.fa 01.data/04.GeneOfInterest/GeneID.txt 01.data/05.splitGenes/01.Protein/run.0 1000
#time perl 00.script/a4.splitGene.pl 01.data/00.PriorData/transcriptome.fa 01.data/04.GeneOfInterest/GeneID.txt 01.data/05.splitGenes/02.Transcript/run.0 1000

#time perl 00.script/a5.releventInfo.pl 01.data/04.GeneOfInterest/GeneID.txt 01.data/00.PriorData/proteome.fa 01.data/00.PriorData/transcriptome.fa 01.data/00.PriorData/gene.gff3 01.data/04.GeneOfInterest/GeneID.v1.txt 1000

###### preprocess data
#time perl 00.script/02.makeblastdb.folder.pl 01.data/05.splitGenes/01.Protein/run.0 prot DIAMOND $platform
#time perl 00.script/02.makeblastdb.folder.pl 01.data/05.splitGenes/02.Transcript/run.0 nucl FALSE $platform
#time perl 00.script/01.folder.fastaCombinePairedEnd.pl  01.data/01.Fastq $platform
#time perl 00.script/01.fastq2fasta.folder.pl 01.data/01.Fastq 01.data/02.Fasta $platform
#time perl 00.script/01.folder.IDConverter.pl 01.data/02.Fasta $platform

###################################################l
a=5
while [ $a -le 11 ]
do
    b=`expr $a + 1`
    echo "The run: $b" >> job.monitor.txt

    if [ $b -eq 0 ];then
#        time perl 00.script/03.diamond.folder.pl 01.data/02.Fasta 01.data/05.splitGenes/01.Protein/run.$b 03.blast/03.bowtie.nucl/run.$b nucl bowtie.log/bowtie.run.$b 1e-10 $platform 
#        time perl 00.script/040.folder.retrievebowtie.reads.pl 03.blast/03.bowtie.nucl/run.$b 04.retrieve.reads/03.bowtie.nucl/run.$b nucl bowtie.log/bowtie.run.$b 1000 genome $platform 20
        echo "Have been here"
    else
        time perl 00.script/021.makebowtiedb.folder.pl 01.data/05.splitGenes/02.Transcript/run.$b $platform 10
        time perl 00.script/03.bowtie.folder.pl 01.data/02.Fasta 01.data/05.splitGenes/02.Transcript/run.$b 03.blast/03.bowtie.nucl/run.$b nucl local bowtie.log/bowtie.run.$b $platform 10
        time perl 00.script/04.folder.retrievebowtie.reads.pl 03.blast/03.bowtie.nucl/run.$b 04.retrieve.reads/03.bowtie.nucl/run.$b nucl bowtie.log/bowtie.run.$b $platform 20
    fi

    time perl 00.script/06.assembly.trinity.folder.pl 04.retrieve.reads/03.bowtie.nucl/run.$b 06.assembly/03.bowtie.nucl/run.$b genome $platform 20
    time perl 00.script/06.truncate.header.folder.pl 06.assembly/03.bowtie.nucl/run.$b $platform 20
    time perl 00.script/07.blastx.back.pl 06.assembly/03.bowtie.nucl/run.$b 01.data/05.splitGenes/01.Protein/run.0 07.map.back/03.bowtie.nucl/run.$b $platform 10
    time perl 00.script/07.blastn.back.pl 06.assembly/03.bowtie.nucl/run.$b 01.data/05.splitGenes/02.Transcript/run.0 07.map.back/02.blastn/run.$b $platform 10    

    c=`expr $b + 1`
    time perl 00.script/100.transfer.saturate.seq.pl 06.assembly/03.bowtie.nucl/run.$b 07.map.back/03.bowtie.nucl/run.$b 01.data/05.splitGenes/02.Transcript/run.$c 01.data/04.GeneOfInterest/GeneID.v1.txt pct 0.02 10 
    a=`expr $a + 1`
done

#### summarize all runs
#time perl 00.script/11.summarize.run.pl 01.data/04.GeneOfInterest/GeneID.v2.txt 11.stable.run/all.summary.run.txt 07.map.back/03.bowtie.nucl 9
#time perl 00.script/b2.extract.gene.of.interest.pl 01.data/00.PriorData/PhQ.txt 07.map.back/03.bowtie.nucl 06.assembly/03.bowtie.nucl 01.data/04.GeneOfInterest/GeneID.v1.txt 11.stable.run/01.blastx/PhQ.blast.summary.txt 11.stable.run/01.blastx/PhQ.contig.seq.fasta

#### assemble unmapped reads
time perl 00.script/c9.get.full.length.seq.pl 08.full.length/count3 08.full.length/full.length.contigs.fasta 08.full.length/Final.fasta

#module load bowtie2/2.2.4
#time bowtie2-build -f -q 08.full.length/Final.fasta 08.full.length/Final
#time perl 00.script/c10.folder.bowtie.full.length.pl 01.data/02.Fasta 08.full.length/Final 09.bowtie.full.length Sapelo
#time perl 00.script/06.truncate.header.pl 10.unmapped.reads.trinity/Trinity.fasta 10.unmapped.reads.trinity/Trinity.new.fasta

#time /usr/local/apps/ncbiblast+/2.2.29/bin/makeblastdb -in 01.data/00.PriorData/proteome.fa -dbtype prot
#time /usr/local/apps/ncbiblast+/2.2.29/bin/makeblastdb -in 01.data/00.PriorData/transcriptome.fa -dbtype nucl
#time /usr/local/apps/ncbiblast+/2.2.29/bin/makeblastdb -in 01.data/00.PriorData/Ath.proteome.fa -dbtype prot
## run examine.sh in 10.unmapped.reads.trinity
#time /usr/local/apps/ncbiblast+/2.2.29/bin/blastx -db ../01.data/00.PriorData/proteome.fa -query Trinity.new.fasta -out blastx.out -evalue 1e-5  -outfmt 6 -num_threads 8 -max_target_seqs 1
#time /usr/local/apps/ncbiblast+/2.2.29/bin/blastx -db ../01.data/00.PriorData/proteome.fa -query Trinity.new.fasta -out blastx.xml.out -evalue 1e-5  -outfmt 5 -num_threads 8 -max_target_seqs 1

### get full length contigs of unmapped reads assembly
#time perl 00.script/101.transfer.saturate.seq.pl 10.unmapped.reads.trinity/Trinity.new.fasta 10.unmapped.reads.trinity/blastx.out 10.unmapped.reads.trinity 01.data/04.GeneOfInterest/GeneID.v1.txt pct 0.02 
cat 10.unmapped.reads.trinity/incomplete.contigs.fasta 10.unmapped.reads.trinity/unmapped.contigs.fasta > 01.data/06.TargetTranscriptome/transcriptome.part2.v1.fa
#time /usr/local/apps/ncbiblast+/2.2.29/bin/blastx -db 01.data/00.PriorData/ath.proteome.fa -query 01.data/06.TargetTranscriptome/transcriptome.part2.v1.fa -out 01.data/06.TargetTranscriptome/transcriptome.part2.v1.Ath.blastx.out -evalue 1e-5 -outfmt 6 -num_threads 8 -max_target_seqs 1
mkdir -p 11.blast.Ath
time perl 00.script/101.transfer.saturate.seq.pl 01.data/06.TargetTranscriptome/transcriptome.part2.v1.fa 01.data/06.TargetTranscriptome/transcriptome.part2.v1.blastx.out 11.blast.Ath 01.data/04.GeneOfInterest/Ath.GeneID.v1.txt pct 0.02

### combine full length contigs from local assembly and assembly of unmapped reads
time /usr/local/apps/ncbiblast+/2.2.29/bin/makeblastdb -in 08.full.length/Final.fasta -dbtype nucl
time /usr/local/apps/ncbiblast+/2.2.29/bin/blastn -db 08.full.length/Final.fasta -query 10.unmapped.reads.trinity/full.length.contigs.fasta -out 10.unmapped.reads.trinity/full.length.contigs.blastn.out -evalue 1e-30 -outfmt 6 -num_threads 4 -max_target_seqs 1
time perl 00.script/c11.combine.full.length.pl 10.unmapped.reads.trinity/full.length.contigs.blastn.xml.out 10.unmapped.reads.trinity/full.length.contigs.fasta 10.unmapped.reads.trinity/temp.fasta
cat 08.full.length/Final.fasta 10.unmapped.reads.trinity/temp.fasta > 01.data/06.TargetTranscriptome/transcriptome.part1.v1.fa
rm 10.unmapped.reads.trinity/temp.fasta

## compare Ath full contigs vs Mimulus full contigs
time /usr/local/apps/ncbiblast+/2.2.29/bin/makeblastdb -in 01.data/06.TargetTranscriptome/transcriptome.part1.v1.fa -dbtype nucl
time /usr/local/apps/ncbiblast+/2.2.29/bin/blastn -db 01.data/06.TargetTranscriptome/transcriptome.part1.v1.fa -query 11.blast.Ath/full.length.contigs.fasta -out 11.blast.Ath/full.length.contigs.blastn.out -evalue 1e-30 -outfmt 6 -num_threads 4 -max_target_seqs 1
time perl 00.script/c11.combine.full.length.pl 11.blast.Ath/full.length.contigs.blastn.xml.out 11.blast.Ath/full.length.contigs.fasta 11.blast.Ath/temp.fasta
cp 11.blast.Ath/temp.fasta 01.data/06.TargetTranscriptome/transcriptome.part2.v2.fa

## collect incomplete and unmapped contigs for length filtering, contamination cleaning and annotation
cat 11.blast.Ath/incomplete.contigs.fasta 11.blast.Ath/unmapped.contigs.fasta > 01.data/06.TargetTranscriptome/transcriptome.part3.v0.fa
time perl 00.script/c12.filter.by.length.pl 01.data/06.TargetTranscriptome/transcriptome.part3.v0.fa 01.data/06.TargetTranscriptome/transcriptome.part3.v1.fa 500

## remove redundant contigs
time perl 00.script/d11.remove.redundant.contigs.pl  01.data/06.TargetTranscriptome/transcriptome.part1.v1.fa  01.data/06.TargetTranscriptome/transcriptome.part1.v2.fa
time perl 00.script/d11.remove.redundant.contigs.pl  01.data/06.TargetTranscriptome/transcriptome.part3.v1.fa  01.data/06.TargetTranscriptome/transcriptome.part3.v2.fa
time perl 00.script/d11.remove.redundant.contigs.pl  01.data/06.TargetTranscriptome/transcriptome.part2.v2.fa  01.data/06.TargetTranscriptome/transcriptome.part2.v3.fa
mkdir -p 12.MCL.remove.redundancy
cd 12.MCL.remove.redundancy
ln -sf ../01.data/06.TargetTranscriptome/transcriptome.part3.v2.fa transcriptome.part3.v2.fa
cd ../
/usr/local/wublast/latest/xdformat -n -o 12.MCL.remove.redundancy/transcriptome.part3.v2 12.MCL.remove.redundancy/transcriptome.part3.v2.fa
/usr/local/wublast/latest/blastn 12.MCL.remove.redundancy/transcriptome.part3.v2 12.MCL.remove.redundancy/transcriptome.part3.v2.fa -o 12.MCL.remove.redundancy/transcriptome.part3.v2.blastn.out -e 1e-5 -mformat 2 -cpus 4 -wordmask seg
time perl 00.script/d12.mcl.prepare.graph.pl 12.MCL.remove.redundancy/transcriptome.part3.v2.blastn.out 12.MCL.remove.redundancy/mcl.graph.txt wu
time /usr/local/mcl/latest/bin/mcl  12.MCL.remove.redundancy/mcl.graph.txt --abc -o 12.MCL.remove.redundancy/mcl.out.txt -I 5.0
time perl 00.script/d13.remove.redundant.mcl.pl 12.MCL.remove.redundancy/mcl.out.txt 12.MCL.remove.redundancy/transcriptome.part3.v2.fa 01.data/06.TargetTranscriptome/transcriptome.part3.v3.fa

## remove contamination from non-plant sequences
mkdir -p 13.remove.contamination/01.plant
mkdir -p 13.remove.contamination/02.bacteria
mkdir -p 13.remove.contamination/03.fungi
mkdir -p 13.remove.contamination/04.unplant
cd 13.remove.contamination/01.plant
ln -sf ../../01.data/06.TargetTranscriptome/transcriptome.part3.v3.fa transcriptome.part3.v3.fa
/usr/local/apps/ncbiblast+/2.2.29/bin/makeblastdb -in transcriptome.part3.v3.fa -dbtype nucl
# download sequences from NCBI RefSeq from three folders
download.sh
# blast assembled contigs against plant proteins
time perl ../../00.script/d14.blast.batch.pl tblastn protein.fa transcriptome.part3.v3.fa tblastn.batch 6 5000 zcluster
time perl ../../00.script/d15.get.unplant.seq.pl tblastn.batch transcriptome.part3.v3.fa ../04.unplant/transcriptome.part3.vNA.fa
cd ../
cat 02.bacteria/RNA.fna 03.fungi/RNA.fna 04.unplant/RNA.fa
# blast non-plant RNA
cd 04.unplant/
/usr/local/apps/ncbiblast+/2.2.29/bin/makeblastdb -in transcriptome.part3.vNA.fa -dbtype nucl
time perl ../../00.script/d14.blast.batch.pl blastn RNA.fa transcriptome.part3.vNA.fa blastn.batch 5 5000 zcluster


## remove contamination from the host

### blast to Ath and Mus for GO annotation, can put into other script to use multiple processors
#time /usr/local/apps/ncbiblast+/2.2.29/bin/blastx -db 01.data/00.PriorData/proteome.fa -query 01.data/06.TargetTranscriptome/transcriptome.part1.v1.fa -out Mus.blastx.out -evalue 1e-30 -outfmt 6 -num_threads 8 -max_target_seqs 1
#time /usr/local/apps/ncbiblast+/2.2.29/bin/blastx -db 01.data/00.PriorData/ath.proteome.fa -query 01.data/06.TargetTranscriptome/transcriptome.part1.v1.fa -out Ath.blastx.out -evalue 1e-30 -outfmt 6 -num_threads 8 -max_target_seqs 1

