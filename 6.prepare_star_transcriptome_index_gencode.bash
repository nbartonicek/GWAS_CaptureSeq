#!/bin/bash

#load all modules
module load gi/samtools/1.0
module load nenbar/star/2.4.0d
#number of cores
ncore=15

#project directory
homedir="/share/ClusterScratch/nenbar"

sample="gencode"
#genome directory
genomeDir="/share/ClusterShare/biodata/contrib/nenbar/genomes/star/hg19_ercc_$sample"
ln -s /share/ClusterShare/biodata/contrib/nenbar/genomes/star/hg19_ercc/hg19_ercc.fa $genomeDir/hg19_ercc.fa
gtfFile="/share/ClusterShare/biodata/contrib/nenbar/genomes/human/hg19/gencode/gencode_ercc.v19.annotation.gtf"
mkdir -p $genomeDir
genome="hg19"

starLine="STAR --runMode genomeGenerate --genomeDir $genomeDir 
	--sjdbGTFfile $gtfFile 
	--sjdbOverhang 99 
	--genomeFastaFiles $genomeDir/hg19_ercc.fa 
	--runThreadN $ncore"
qsub -N star_trin -b y -cwd -j y -R y -pe smp $ncore -V $starLine 



