#!/bin/bash

#This script uses fastqc to for QC control 

module load marcow/fastqc/prebuilt/0.10.1
#directory hierarchy
#raw files directory
homedir="/share/ClusterShare/biodata/contrib/nenbar"

#extension of the files to be used
inExt="fq.gz"

#scripts directory
scriptsPath="$homedir/projects/melanoma/scripts/QC"

#name to append to projectname and create a folder
projectnames=( "melanomaGWAS_1" )
analyses=( "trimgalore" )
imageTypes=( "duplication_levels" "sequence_length_distribution" "kmer_profiles" "per_base_quality" "per_base_sequence_content" )

for projectname in ${projectnames[@]}; do

        outType="figures"
        outPath="/$homedir/projects/melanoma/project_results/"$projectname.$outType
        mkdir -p $outPath
        for analysis in ${analyses[@]};do

                for imageType in ${imageTypes[@]};do

                        montage -label '%d' `ls -d ../../project_results/$projectname.$analysis/**/**/Images/$imageType*` -tile 2x16 -geometry +15+15 $outPath/$imageType.$analysis.png
                        echo "$imageType.$analysis.png"
                done;
        done;
done;

