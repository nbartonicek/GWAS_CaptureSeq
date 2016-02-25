#$ -S /bin/bash
module load gi/boost/1.53.0
module load gi/samtools/0.1.19
module load gi/cufflinks/2.2.1

sample="melanoma"
inDir="/share/ClusterShare/biodata/contrib/nenbar/genomes/rsem/hg19_ercc_gencode19/"
genomeFile="/share/ClusterShare/biodata/contrib/nenbar/genomes/star/hg19_ercc/hg19_ercc.fa"

name="gencode_ercc.v19.annotation"
#remove the regions without strands
strandedFile=$name"_stranded.gtf"
#rm $strandedFile
file=`ls -d $inDir/$name.gtf`
more $file | perl -ne '@cols=split("\t",$_);if($cols[6]!~/\./){print $_}' >$inDir$strandedFile

gffread -w $inDir/$name.fa -g $genomeFile $inDir/$strandedFile


files=`ls -d $inDir/$name.fa`
for file in ${files[@]};do 
	
	echo $name;
	cat $file | grep ">" | sed s/\>//| sort | uniq | perl -ne 'chomp;$gene = $_;$transcript=$_;$transcript=~s/\s.*//;$gene=~s/.*gene=//;$gene=~s/\s.*//;print "$gene\t$transcript\n"' >$inDir$name.txt 
	#cat $file | cut -f 9 | sed s/\\..*// | sed s/.*=// | sort | uniq | perl -ne 'chomp;$gene = $_;$transcript=$_;$gene=~s/_i.*//;print "$gene\t$transcript\n"' >$name.txt
done; 
