#$ -S /bin/bash

module load gi/boost/1.53.0
module load gi/samtools/0.1.19
module load nenbar/rsem/1.2.18

numcores=15

#directory hierarchy
#raw files directory
homedir="/share/ClusterShare/biodata/contrib/nenbar"
projectDir="$homedir/projects/melanoma"
resultsDir="$projectDir/project_results/"
tempDir="/share/ClusterScratch/nenbar"

genomeDir="/share/ClusterShare/biodata/contrib/nenbar/genomes/star/hg19_ercc_trinity"
rsem_index="/share/ClusterShare/biodata/contrib/nenbar/genomes/rsem/cuffmerged_alltissues_gencode"
#extension of the files to be used
inExt="gz"

#scripts directory
scriptsPath="$homedir/projects/melanoma/scripts/trancsript_assembly"
logDir=$projectDir"/scripts/transcript_assembly/logs"

#name to append to projectname and create a folder
inType="benchmark"

projectnames=( "tissues" )
		
rsem_line="rsem-calculate-expression -p 12 \
	--bam --paired-end 
	--no-bam-output \
	--estimate-rspd AlignedToTranscriptome.out.bam \
	/share/ClusterShare/biodata/contrib/nenbar/genomes/rsem/cuffmerged_alltissues/ $outDir/$name"
echo $rsem_line                

#qsub -N RSEM_count_$outPathBit -hold_jid RSEM_index_$outPathBit -wd $logDir -b y -j y -R y -pe smp $numcores -V $rsem_line

