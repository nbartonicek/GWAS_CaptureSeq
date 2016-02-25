#!/bin/bash
module load gi/gcc/4.8.2
module load nenbar/star/2.4.0d
module load gi/samtools/1.1
module load fabbus/python/2.7.3
module load gi/novosort/1.02.02

numcores=15

#directory hierarchy
#raw files directory
homedir="/share/ClusterShare/biodata/contrib/nenbar"
projectDir="$homedir/projects/melanoma"
resultsDir="$projectDir/project_results/"
tempDir="/share/ClusterShare/biodata/contrib/nenbar"

projectname=( "melanoma" )

genomeDir="/share/ClusterShare/biodata/contrib/nenbar/genomes/star/hg19_ercc_"$projectname
rsem_index="/share/ClusterShare/biodata/contrib/nenbar/genomes/rsem/cuffmerged_"$projectname"_gencode/melanoma_combined"
#extension of the files to be used
inExt="fq.gz"

#scripts directory
scriptsPath="$homedir/projects/melanoma/scripts/trancsript_assembly"
logDir=$projectDir"/scripts/transcript_assembly/logs"

#name to append to projectname and create a folder
inType="benchmark"


#out directory
inPath="$homedir/projects/melanoma/project_results/$projectname.test/raw_files"
#inPath="$tempDir/capseq/trimmed_files/$projectname"
outTool1="star"
outPath1="$tempDir/capseq/project_results/$projectname.$outTool1"

outTool2="rsem"
outPath2="$projectDir/project_results/$projectname.$outTool2"

#log and command files for bsub
logPath="logs"
commandPath="commands"
#make the directory structure   
mkdir -p $outPath1
mkdir -p $outPath2
mkdir -p $logPath
mkdir -p $commandPath

rm -f $commandFile


#get the name of the script for the logs
scriptName=`basename $0`
i=0
echo $inPath
files=`ls $inPath/*$inExt`

for file in ${files[@]};do
                #echo The file used is: $file
                filesTotal[i]=$file;
                let i++;
done 

j=0
echo ${#filesTotal[@]}
while [ $j -lt ${#filesTotal[@]} ]; do

        
	inFile1=${filesTotal[$j]}
	inFile2=${filesTotal[$(($j+1))]}
        #uniqueID=`basename $inFile1 | sed s/_R1.fastq.gz//`
        uniqueID=`basename $inFile1 | sed s/_1.fq//`
        name=$uniqueID
        outDir1=$outPath1/$uniqueID/
        outDir2=$outPath2/$uniqueID
	mkdir -p $outDir1
        mkdir -p $outDir2       
        echo $name
	#echo $command_line


	starJobName="star."$name
	rsemJobName="rsem."$name
        sortJobName="sort."$name
        filterJobName="filter."$name
        htseqJobName="htseq."$name

	outBam=$outDir1"Aligned.toTranscriptome.out.bam"
        outSortedBam=$outDir1"Aligned.toTranscriptome.sorted.out"
        outFilteredBam=$outDir1"Aligned.toTranscriptome.filtered.out.bam"
        outNameFilteredBam=$outDir1"Aligned.toTranscriptome.filtered.names.sorted"
        fixmate=$outDir1"Aligned.toTranscriptome.fixmate.bam"
        nongapped=$outDir1"Aligned.toTranscriptome.nongapped2.bam"

        #--readFilesCommand zcat \
	star_line="STAR --runMode alignReads \
                --readFilesCommand zcat \
		--genomeDir $genomeDir \
                --outFilterType BySJout \
                --outSAMattributes NH HI AS NM MD\
                --outFilterMultimapNmax 20 \
                --outFilterMismatchNmax 999 \
                --outFilterMismatchNoverReadLmax 0.04 \
                --alignIntronMin 20 \
                --alignIntronMax 1500000 \
                --alignMatesGapMax 1500000 \
                --alignSJoverhangMin 6 \
                --alignSJDBoverhangMin 1 \
                --readFilesIn $inFile1 $inFile2 \
                --outFileNamePrefix $outDir1 \
                --runThreadN $numcores \
                --quantMode TranscriptomeSAM \
                --outFilterMatchNmin 101" 
        #sort_line="samtools sort -n -m 16G -@ $numcores $outBam $outSortedBam"
        filter_line="~/local/lib/samtools-1.1/samtools view -m 16G -@ $numcores $outBam -f 3 -b > $outFilteredBam"
        #sortname_line="~/local/lib/samtools-1.1/samtools sort -n -m 16G -@ $numcores $outFilteredBam $outSortedBam"
        sortname_line="novosort -n -m 16G -c $numcores $outFilteredBam >$outSortedBam.bam"
        index_line="samtools index $outSortedBam.bam"
        #nongapped_line="samtools view -m 16G -@ $numcores $outBam -f 2 -b | samtools sort -n -m 16G -@ $numcores - $nongapped"


        rsem_line="~/local/lib/rsem-1.2.15/rsem-calculate-expression \
                --paired-end \
                --bam \
                --no-bam-output \
                --seed 12345 -p $numcores \
                --forward-prob 0 \
                $outSortedBam.bam \
                $rsem_index  \
                $outDir2"
        echo $index_line
        
        #echo $rsem_line
        #echo $htseq_line
        #echo "samtools idxstats $outSortedBam | head -n1 | cut -f3 >$outDir$name.tmp;"
	#qsub -N $starJobName -b y -wd $logDir -j y -R y -pe smp $numcores -V $star_line 
        #qsub -N $filterJobName -b y -hold_jid $starJobName -wd $logDir -j y -R y -pe smp $numcores -l h_vmem=17G,mem_requested=16G -V $filter_line 
        #qsub -N $sortJobName -b y -hold_jid $filterJobName -wd $logDir -j y -R y -pe smp $numcores -l h_vmem=17G,mem_requested=16G -V $sortname_line 
        qsub -N $rsemJobName -hold_jid $sortJobName -b y -wd $logDir -j y -R y -pe smp $numcores -V $rsem_line

	#qsub -N $indexJobName -hold_jid $sortJobName -b y -wd $logDir -j y -R y -pe smp 1 -V index_line
        #qsub -N $htseqJobName -hold_jid $starJobName -b y -wd $logDir -j y -R y -pe smp $numcores -V $htseq_line
        
        j=$(($j+2))


done;


