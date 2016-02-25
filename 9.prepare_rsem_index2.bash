#$ -S /bin/bash

module load gi/boost/1.53.0
module load gi/bowtie/1.0.0 
module load borgue/rsem/1.2.12

numcores=1 
sample="melanoma"
type="melanoma"
inDir="../../project_results/transcript_assembly.cuffmerge/$sample.gencode_trinity"
outDir="/share/ClusterShare/biodata/contrib/nenbar/genomes/rsem/cuffmerged_melanoma_gencode/"
mkdir -p $outDir
cp $inDir/*.fa $outDir/
cp $inDir/*.txt $outDir/


#echo $file
rsem_index_line="~/local/lib/rsem-1.2.15/rsem-prepare-reference --no-polyA \
	$outDir/$type.fa \
	--transcript-to-gene-map $outDir/$type.txt \
    $outDir$type"
#echo $rsem_index_line
qsub -N RSEM_index -cwd -b y -j y -R y -pe smp 1 -V $rsem_index_line

