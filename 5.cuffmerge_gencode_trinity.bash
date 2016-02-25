#$ -S /bin/bash

module load gi/boost/1.53.0
module load gi/samtools/0.1.19
module load gi/cufflinks/2.2.1

numcores=5


homedir="/share/ClusterShare/biodata/contrib/nenbar"
projectDir="$homedir/projects/melanoma"
resultsDir="$projectDir/project_results/"
projectname="transcript_assembly"

scriptsPath="$projectDir/scripts/$projectname"
logDir=$scriptsPath"/logs"
mkdir -p $logDir

annotationFile="/share/ClusterShare/biodata/contrib/shinyRseq/hg19_ercc/genome.fa"
gencodeFile="/share/ClusterShare/biodata/contrib/nenbar/genomes/human/hg19/gencode/gencode_ercc.v19.annotation.gtf"

sample="melanoma"
#input/output
inTool="cuffmerge"
outTool="cuffmerge"
inDir="$projectDir/project_results/$projectname.$inTool/$sample/"
inFile="$projectDir/project_results/$projectname.$inTool/$sample/merged.gtf"
manifestFile="$inDir/assembly_list.txt"
outDir="$projectDir/project_results/$projectname.$outTool/$sample.gencode_trinity/"
#outDir="$projectDir/project_results/$projectname.$outTool/gencode/"

mkdir -p $outDir

#create manifest file
echo $gencodeFile > $manifestFile
echo $inFile >> $manifestFile


cuffmerge_line="cuffmerge -p $numcores -s $annotationFile -o $outDir $manifestFile"
qsub -N cuffmerge_gencode -wd $logDir -hold_jid cuffmerge_guided -b y -j y -R y -pe smp $numcores -V $cuffmerge_line


