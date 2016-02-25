#$ -S /bin/bash

#module load marsmi/trinityrnaseq/2014-04-13p1
#module load marsmi/gmap-gsnap/2014-06-10
module unload fabbus/perl/5.14.2
module load fabbus/perl/5.14.2
module load marsmi/bowtie/0.12.8
module load gi/samtools/0.1.19
module load marsmi/java/1.6.0_37
#module load marsmi/java/1.6.0_37

#genomeFile="/share/ClusterShare/biodata/contrib/genomeIndices_garvan/pwbc/genomes/Hsapiens/hg19/bowtie2/hg19.fa"
genomeFile="/share/ClusterShare/biodata/contrib/shinyRseq/hg19_ercc/genome.fa"
homedir="/share/ClusterShare/biodata/contrib/nenbar"
projectDir="$homedir/projects/melanoma"
resultsDir="$projectDir/project_results/"
projectname="transcript_assembly"

scriptsPath="$projectDir/scripts/$projectname"
logDir=$scriptsPath"/logs"
mkdir -p $logDir

numcores=20

inPathBit="liver"
inPath=$projectDir"/project_results/"$projectname".benchmark/raw_files/$inPathBit"

#position of output on clusterScratch
tempDir="/share/ClusterScratch/nenbar/shinyRseq/assembled/trinity"
mkdir tempDir

#Get Files
#files=`ls $inPath`
files=($(ls -d $inPath/*.fq))

#Set up conditions
versions=( "old" "new" )
genomes=( "genome" )
#genomes=( "nogenome" "genome" )
#grids=( "nogrid" "grid" )
grids=( "nogrid" )
for version in ${versions[@]};do

	for genome in ${genomes[@]};do

		#for grid in ${grid[@]};do

			#output
			outPathBit="capseq_"$version\_$genome
			tempOutPath="$tempDir/$inPathBit/$outPathBit/"
			echo $outPathBit 
			mkdir -p $tempOutPath
			outDir=$projectDir"/project_results/"$projectname".trinity/$inPathBit/$outPathBit"
			mkdir -p $outDir

			if [[ ($version =~ "old") ]]; then 
				trinity_version="/home/nenbar/local/lib/trinityrnaseq_r20140413p1/Trinity"
			elif [[ ($version =~ "new") ]]; then 
				trinity_version="/home/nenbar/local/lib/trinityrnaseq_r20140710beta/Trinity"
			fi

			if [[ ($genome =~ "nogenome") ]]; then 
				genome_line=""
			else
				genome_line="--genome $genomeFile \
					--genome_guided_max_intron 1000000 \
					--genome_guided_sort_buffer 16G \
					--genome_guided_CPU 2 \
					--GMAP_CPU $numcores"
			fi

			#if [[ ($grid =~ "nogrid") ]]; then 
			#	grid_line=""
			#else
			#	grid_line="--grid_conf /home/nenbar/conf/wolfpack.conf"
			#fi

			trinity_line="$trinity_version --seqType fq \
				--JM 8G \
				--left ${files[0]} --right ${files[1]} \
				--CPU $numcores \
				--normalize_reads \
				--normalize_max_read_cov 50 \
				--SS_lib_type FR \
				--output $tempOutPath \
				--full_cleanup
				$genome_line" 
				#\"
				#$grid_line"
			trinity_command="trinity_command_$outPathBit"
			echo "#$ -S /bin/bash" >$trinity_command
			echo $trinity_line >>$trinity_command
			qsub -N tri_$outPathBit -wd $logDir -j y -R y -pe smp $numcores -V ./$trinity_command
		#done;
	done;
done;

#--full_cleanup \
# \/home/nenbar/local/lib/trinityrnaseq_r20140413p1/
#	--grid_conf /home/nenbar/conf/wolfpack.conf
#trinity_line_genome="Trinity --seqType fq /
#	--JM 2G /
#	----genome genome.fasta
#	--left ${files[0]} --right ${files[1]} /
#	--CPU $numcores /
#	--normalize_max_read_cov 50 /
#	--SS_lib_type FR /
#	--output $outDir /
#	--tmp_dir_name /share/ClusterScratch/nenbar/capseq/trinity /
#	--PARALLE:wL_STATS"

#Trinity --genome genome.fasta --genome_guided_max_intron 10000 --genome_guided_sort_buffer 10G \
#        --genome_guided_CPU 4 --GMAP_CPU 10 \
#        --seqType fq --JM 2G --left reads_1.fq  --right reads_2.fq --CPU 10

#echo "#$ -S /bin/bash" >trinity_line
#echo $trinity_line >>trinity_line
#qsub -N trinity_bm -wd $logDir -j y -R y -pe smp $numcores -V ./trinity_line
