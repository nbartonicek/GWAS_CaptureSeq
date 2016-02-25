#$ -S /bin/bash

module load gi/boost/1.53.0
module load gi/samtools/0.1.19
module load gi/cufflinks/2.2.1

numcores=10


homedir="/share/ClusterShare/biodata/contrib/nenbar"
projectDir="$homedir/projects/melanoma"
resultsDir="$projectDir/project_results/"
projectname="transcript_assembly"

scriptsPath="$projectDir/scripts/$projectname"
logDir=$scriptsPath"/logs"
mkdir -p $logDir

sample="all_tissues"
genomeFile="/share/ClusterShare/biodata/contrib/shinyRseq/hg19_ercc/genome.fa"
annotationFile="/share/ClusterShare/biodata/contrib/shinyRseq/hg19_ercc/gencode19_ercc.gtf"


#input/output
inTool="gsnap"
inDir="$projectDir/project_results/$projectname.$inTool/all_tissues/"
outTool="cuffcompare"


manifestFile="$inDir/assembly_list.txt"
ls -d $inDir/* > $manifestFile

outDir="$projectDir/project_results/$projectname.$outTool/all_tissues_nongenome/"
mkdir -p $outDir
cuffcompare_line="cuffcompare -G -r $annotationFile -s $genomeFile -p $outDir -i $manifestFile"
echo $cuffcompare_line
#qsub -N cuffcompare -wd $logDir -b y -j y -R y -pe smp $numcores -V $cuffcompare_line


