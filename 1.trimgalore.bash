#!/bin/bash

module load hugfre/trimgalore/0.2.8
module load gi/fastqc/0.10.1

############## directory hierarchy ##############
#raw files directory
homedir="/share/ClusterShare/biodata/contrib/nenbar"
#inPath="$homedir/projects/MRTA/raw_files/"
inPath="/share/ClusterScratch/nenbar/MRTA/raw_files/"

#extension of the files to be used
inExt="gz"

#scripts directory
scriptsPath="$homedir/projects/MRTA/scripts/QC"

#name to append to projectname and create a folder
inType="trimgalore"
projectname="ATT_1"

#out directory
outPath="/$homedir/projects/MRTA/project_results/"$projectname.$inType
outPath="/share/ClusterScratch/nenbar/MRTA/project_results/"$projectname.$inType

#log and command files for bsub
logDir=$scriptsPath/"logs"
commandPath="commands"

#make the directory structure   
mkdir -p $outPath
mkdir -p $logDir
mkdir -p $commandPath
rm -f $commandFile

############## fetch file names ##############

i=0   
files=( $(ls $inPath/*) )
for file in ${files[@]};do
        echo The file used is: $file
        filesTotal[i]=$file;
        let i++;
done;

############## perform analysis in pairs ##############


j=0
echo -e "The total number of files is:"
echo ${#filesTotal[@]}
echo -e

while [ $j -lt ${#filesTotal[@]} ]; do

        inFile1=${files[$j]}
        inFile2=${files[$(($j+1))]}

        uniqueID=`basename $inFile1 | sed 's/-RNA.*_L00//' | sed s/_R1_001.fastq.gz// | sed s/ATT1-//`
        echo $uniqueID
        name=$uniqueID
        outDir=$outPath/$uniqueID/
        mkdir -p $outDir
        echo $name
        #echo $command_line

        command_line="trim_galore $inFile1 $inFile2 --fastqc --paired --retain_unpaired --length 16 -o $outDir"
        #echo $command_line
        qsub -b y -wd $logDir -j y -N trimgalore -R y -pe smp 1 -V $command_line
        j=$(($j+2))

done;

