#$ -S /bin/bash

#module load marsmi/trinityrnaseq/2014-04-13p1
#module load marsmi/gmap-gsnap/2014-06-10
module unload fabbus/perl/5.14.2
module load fabbus/perl/5.14.2
module load marsmi/bowtie/0.12.8
module load gi/samtools/0.1.19
module load marsmi/java/1.6.0_37
#module load marsmi/java/1.6.0_37

numcores=5 

homedir="/share/ClusterShare/biodata/contrib/nenbar"
projectDir="$homedir/projects/melanoma"
resultsDir="$projectDir/project_results/"
projectname="transcript_assembly"

scriptsPath="$projectDir/scripts/$projectname"
logDir=$scriptsPath"/logs"
mkdir -p $logDir

sample="capseq"
annotation="hg19_ercc"

#input/output
inTool="trinity"
inPath=$projectDir"/project_results/$projectname.$sample"
outTool="gsnap"
outDir="$projectDir/project_results/$projectname.$outTool/all_tissues/"
mkdir -p $outDir

files=($(ls -d $inPath/pbmc.fasta))



for inFile in ${files[@]};do

		outPathBit=`basename $inFile | sed 's/.fasta//'`
		echo $outPathBit
		outGff=$outDir/$outPathBit".gff"
		gsnap_params="-K 2000000 -w 4000000 -L 5000000 --canonical-mode=0 "
		gmap_gff_line="/share/ClusterShare/software/contrib/nenbar/gmap-gsnap/2014-02-28/gmap \
			-t $numcores -n 0 -d $annotation $inFile -f gff3_gene > $outGff"
		qsub -N gff\_$sample\_$outPathBit -wd $logDir -b y -j y -R y -pe smp $numcores -V $gmap_gff_line

done;
