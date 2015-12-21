#!/bin/bash
#PBS -q batch
#PBS -N localAssembly_allsample
#PBS -l nodes=1:ppn=1:AMD
#PBS -l walltime=4800:00:00

#PBS -M guxi798@hotmail.com
#PBS -m ae

cd $PBS_O_WORKDIR

####################### this starts with proteome #######################
# WU-BLAST
#/usr/local/wublast/latest/xdformat -p -o 01.data/03.MCL/01.Protein/01.blast/database 01.data/00.PriorData/Ptrichocarpa_210_v3.0.protein_primaryTranscriptOnly.fa
#/usr/local/wublast/latest/blastp 01.data/03.MCL/01.Protein/01.blast/database 01.data/00.PriorData/Ptrichocarpa_210_v3.0.protein_primaryTranscriptOnly.fa -o 01.data/03.MCL/01.Protein/01.blast/blast.all.out -e 1e-5 -mformat 2 -cpus 4 -wordmask seg

# NCBI-BLAST
#/usr/local/ncbiblast+/latest/bin/makeblastdb -in 01.data/00.PriorData/Ptrichocarpa_210_v3.0.protein_primaryTranscriptOnly.fa -dbtype prot
#time /usr/local/ncbiblast+/latest/bin/blastp -num_threads 4 -db 01.data/00.PriorData/Ptrichocarpa_210_v3.0.protein_primaryTranscriptOnly.fa -query 01.data/00.PriorData/Ptrichocarpa_210_v3.0.protein_primaryTranscriptOnly.fa -out 01.data/03.MCL/01.Protein/01.blast/ncbi.blast.all.out -evalue 1e-5 -outfmt 6

# mcl
#time perl 00.script/a1.mcl.prepare.graph.pl 01.data/03.MCL/01.Protein/01.blast/ncbi.blast.all.out 01.data/03.MCL/01.Protein/02.mcl/mcl.graph.txt
#time /usr/local/mcl/latest/bin/mcl 01.data/03.MCL/01.Protein/02.mcl/mcl.graph.txt --abc -o 01.data/03.MCL/01.Protein/02.mcl/mcl.out.txt -I 1.5

#time perl 00.script/a1.mcl.prepare.graph.pl 01.data/03.MCL/01.Protein/01.blast/wu.blast.all.out 01.data/03.MCL/01.Protein/02.mcl/wu.mcl.graph.txt wu
#time /usr/local/mcl/latest/bin/mcl 01.data/03.MCL/01.Protein/02.mcl/wu.mcl.graph.txt --abc -o 01.data/03.MCL/01.Protein/02.mcl/wu.mcl.out.txt -I 1.5

# split gene
#time perl 00.script/a4.splitGene.pl 01.data/00.PriorData/proteome.fa 01.data/04.GeneOfInterest/GeneID.txt 01.data/05.splitGenes/01.Protein 1000

#time perl 00.script/a5.releventInfo.pl 01.data/04.GeneOfInterest/GeneID.txt 01.data/00.PriorData/proteome.fa 01.data/00.PriorData/transcriptome.fa 01.data/00.PriorData/Ptrichocarpa_210_v3.0.gene.gff3 01.data/04.GeneOfInterest/GeneID.v1.txt

###### preprocess data
time perl 00.script/02.makeblastdb.folder.pl 01.data/05.splitGenes/01.Protein/run.0 prot
time perl 00.script/01.fastq2fasta.folder.pl 01.data/01.Fastq 01.data/02.Fasta 10
time perl 00.script/01.folder.fastaCombinePairedEnd.pl  01.data/02.Fasta 60
time perl 00.script/01.folder.IDConverter.pl  01.data/02.Fasta

###################################################
platform="Sapelo"

a=-1
while [ $a -le 20 ]
do
    b=`expr $a + 1`
    echo "The run: $b" >> job.monitor.txt

    if [ $b -eq 0 ];then
        time perl 00.script/03.diamond.folder.pl 01.data/02.Fasta 01.data/05.splitGenes/01.Protein/run.$b 03.blast/03.bowtie.nucl/run.$b nucl bowtie.log/bowtie.run.$b $platform
        time perl 00.script/040.folder.retrievebowtie.reads.pl 03.blast/03.bowtie.nucl/run.$b 04.retrieve.reads/03.bowtie.nucl/run.$b nucl bowtie.log/bowtie.run.$b 1000 $platform 20
    else
        time perl 00.script/021.makebowtiedb.folder.pl 01.data/05.splitGenes/02.Transcript/run.$b 10
        time perl 00.script/03.bowtie.folder.pl 01.data/02.Fasta 01.data/05.splitGenes/02.Transcript/run.$b 03.blast/03.bowtie.nucl/run.$b nucl local bowtie.log/bowtie.run.$b $platform 10
        time perl 00.script/04.folder.retrievebowtie.reads.pl 03.blast/03.bowtie.nucl/run.$b 04.retrieve.reads/03.bowtie.nucl/run.$b nucl bowtie.log/bowtie.run.$b $platform 20
    fi

    time perl 00.script/06.assembly.trinity.folder.pl 04.retrieve.reads/03.bowtie.nucl/run.$b 06.assembly/03.bowtie.nucl/run.$b $platform 60
    time perl 00.script/06.truncate.header.folder.pl 06.assembly/03.bowtie.nucl/run.$b/01.trinity $platform 600
    time perl 00.script/07.blastx.back.pl 06.assembly/03.bowtie.nucl/run.$b/01.trinity 01.data/05.splitGenes/01.Protein/run.0 07.map.back/03.bowtie.nucl/run.$b $platform 600
    
    if [ $b -gt 0 ];then
        time perl 00.script/08.folder.summarize.blastx.pl 07.map.back/03.bowtie.nucl/run.$a 07.map.back/03.bowtie.nucl/run.$b 08.evaluate/03.bowtie.nucl/run.$b 10
        time perl 00.script/09.folder.plot.blastx.pl 08.evaluate/03.bowtie.nucl/run.$b 09.plot/03.bowtie.nucl/run.$b 10
    else
        time perl 00.script/09.folder.plot.touch.pl 07.map.back/03.bowtie.nucl/run.$b 09.plot/03.bowtie.nucl/run.$b 10
    fi

    c=`expr $b + 1`
    time perl 00.script/10.transfer.saturate.seq.pl 06.assembly/03.bowtie.nucl/run.$b/01.trinity 07.map.back/03.bowtie.nucl/run.$b 01.data/05.splitGenes/02.Transcript/run.$c 01.data/04.GeneOfInterest/GeneID.v1.txt 10
    a=`expr $a + 1`
done

#### summarize all runs
#time perl 00.script/11.summarize.run.pl 01.data/04.GeneOfInterest/GeneID.v1.txt 11.stable.run/all.summary.run.txt 07.map.back/03.bowtie.nucl 9
