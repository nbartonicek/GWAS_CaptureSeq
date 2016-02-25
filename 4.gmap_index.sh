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
annotation="Rn4_UCSC"

qsub -N gmapidx_rat -wd $logDir -b y -j y -R y -pe smp $numcores -V /share/ClusterShare/software/contrib/nenbar/gmap-gsnap/2014-02-28/gmap_build -D . -d $annotation /share/ClusterShare/biodata/contrib/nenbar/genomes/gsnap/Rn4_UCSC.fa

