
library(ShortRead)
library(rtracklayer)
library(ggplot2)
library(reshape2)
library(scales)
library(RUVSeq)
library(edgeR)
library(limma)

homedir="/share/ClusterShare/biodata/contrib/nenbar"
#enable for Rstudio
#homedir="../../../../"

timeStamp <- format(Sys.time(), "%Y_%m_%d")


######## directory structure #######
projectDir=paste0(homedir,"/projects/melanoma")
resultsDir=paste0(projectDir,"/project_results")
imageDir=paste0(resultsDir,"/figures")
annotationDir=paste0(projectDir,"/annotation/ERCC")
projectname="tissues"
robjectsDir = paste(resultsDir,"/gwas_capseq.Robjects/",sep="")
scriptsPath=paste(projectDir,"/scripts/",projectname)
logDir=paste0(scriptsPath,"/logs")
inDir=paste0(projectDir,"/project_results/",projectname,".rsem")
 

#Get the ERCC information to remove before normalization
ercc_transcripts_file=paste0(homedir,"/genomes/rsem/cuffmerged_",projectname,"_gencode/ercc.gtf")
ercc_txs=import(ercc_transcripts_file)
ercc_txs_df=as.data.frame(values(ercc_txs))
ercc_txs_df$ercc<-as.character(seqnames(ercc_txs))

############### NORMALIZATION ##################

#Fetch RSEM output for genes
files<-list.files(inDir,full.names=T)
files<-files[!grepl("pbmc",files)]

resultsCounts<-list()
resultsCountsFile=paste0("../../project_results/transcript_assembly.Rdata/resultsCounts_",projectname,".Rdata")
#if(!file.exists(resultsCountsFile)){
	for(sampleFile in files){
	  sampleName=basename(sampleFile)
	  cat(sampleName)
	  cat("\t")
	  fileName<-paste0(sampleFile,"/",sampleName,".genes.results")
	  data=read.table(fileName,header=T,sep="\t") 
	  resultsCounts[[sampleName]]<-as.integer(data$expected_count)
	}
#} else {load(resultsCountsFile)}
save(resultsCounts,file=resultsCountsFile)
alltissue_counts=as.data.frame(resultsCounts)
alltissue_counts=alltissue_counts[!grepl("pbmc",names(alltissue_counts))]
alltissue_counts$gene_id=data$gene_id
alltissue_counts$transcript_id=data$transcript_id.s.

#remove ERCC
alltissue_counts=alltissue_counts[!alltissue_counts$gene_id %in% ercc_txs_df$gene_id,]
alltissue_counts<-alltissue_counts[1:22]
row.names(alltissue_counts)=alltissue_counts[,22]
alltissue_counts=alltissue_counts[,-22]
filter <- apply(alltissue_counts, 1, function(x) length(x[x>10])>=1)
filtered <- alltissue_counts[filter,]

#x <- as.factor(colnames(alltissue_counts))
#set <- newSeqExpressionSet(as.matrix(filtered), phenoData = data.frame(x, row.names=colnames(filtered)))

########## Normalize means by DESeq normalization factors, and dispersion with quantile normalization
library(DESeq)
countTable=filtered
condition=as.factor(rep(c("1", "2", "1",  "3", "4", "5"), c(3,1,2,5,5,5)))
cds = newCountDataSet( countTable, condition )
cds = estimateSizeFactors( cds )
cds = estimateDispersions( cds, method="blind", sharingMode="fit-only" )
normalized.count <- counts(cds, normalized = TRUE)
normalized.count<-normalizeQuantiles(normalized.count)

################ Make transcripts for testing ##################

#get width to calculate RPKMs
load("../../project_results/transcript_assembly.Rdata/RSEM_tissues.Rdata")
gene_width<-alltissue_fpkms[,c("gene_id","width")]
merged<-merge(normalized.count,gene_width,by.x=0,by.y="gene_id")
rpkms=rpkm(merged[2:22],gene.length=merged$width, lib.size=mean(colSums(normalized.count)))
row.names(rpkms)=merged$Row.names


######### Load chromosomal locations ###############
load(paste0("../../project_results/transcript_assembly.Rdata/merged_gtf_",projectname,".Rdata"))

######### Get normalized counts and coords ##############
normGR<-tl[values(tl)$gene_id %in% row.names(rpkms)]

#### annotate normGR
sampleNumber=dim(rpkms)[2]

data=as.data.frame(values(normGR))
merged<-merge(data,rpkms,by.x="gene_id",by.y=0,all.x=T)
merged$maxFPKM<-apply(merged[,11:(11+sampleNumber-1)],1,max)
annotGR<-normGR
values(annotGR)=merged

#Load the capture regions
load(paste0("../../project_results/transcript_assembly.Rdata/capture_regions_",projectname,".Rdata"))
load(paste0("../../project_results/transcript_assembly.Rdata/capture_regions_short_",projectname,".Rdata"))

##################### Find overlaps with regions #######################


final_regions=capture_regions_short[[1]]
mat<-findOverlaps(annotGR,final_regions)

#make the total, single and multi objects
captureTranscripts<-annotGR[unique(queryHits(mat))]
captureGeneIds=unique(values(captureTranscripts)$gene_id)
captureTranscriptIds=unique(values(captureTranscripts)$transcript_id)


#TRANSCRIPTS
captureTxExonCount<-table(captureTranscripts$transcript_id)
captureMultiexonTranscriptIds<-names(captureTxExonCount)[captureTxExonCount>1]
captureSingleTranscriptIds=captureTranscriptIds[!captureTranscriptIds %in% captureMultiexonTranscriptIds]

#GENES
captureMultiexonGeneIds<-unique(values(captureTranscripts)$gene_id[values(captureTranscripts)$transcript_id %in% captureMultiexonTranscriptIds])
captureSingleGeneIds=captureGeneIds[!captureGeneIds %in% captureMultiexonGeneIds]
save(captureMultiexonGeneIds,file="../../project_results/transcript_assembly.Rdata/captureMultiexonGeneIds_tissues.Rdata")
save(captureSingleGeneIds,file="../../project_results/transcript_assembly.Rdata/captureSingleGeneIds_tissues.Rdata")

#34797 Total genes
#Total of 11766 genes overlap the capture regions
#Total of 1832 Multiexonic, 9934 Singletons 
sampleNumber=dim(rpkms)[2]
#maxFPKM<-apply(rpkms[,1:(sampleNumber)],1,max)

save(rpkms,file="../../project_results/transcript_assembly.Rdata/all_rpkms_tissues.Rdata")

total_rpkms=rpkms[row.names(rpkms) %in% captureGeneIds,]
#save(total_rpkms,file="../../project_results/transcript_assembly.Rdata/total_rpkms_tissues.Rdata")
save(total_rpkms,file="../../project_results/transcript_assembly.Rdata/capture_rpkms_tissues.Rdata")

maxFPKM<-apply(total_rpkms[,1:(sampleNumber)],1,max)

tempRpkms=rpkms
rpkms<-rpkms[row.names(rpkms) %in% captureMultiexonGeneIds,]
#save(rpkms,file="../../project_results/transcript_assembly.Rdata/rpkms_tissues.Rdata")
save(rpkms,file="../../project_results/transcript_assembly.Rdata/capture_multi_rpkms_tissues.Rdata")

rpkms<-tempRpkms[row.names(tempRpkms) %in% captureSingleGeneIds,]
save(rpkms,file="../../project_results/transcript_assembly.Rdata/capture_single_rpkms_tissues.Rdata")


#OBJECTS
captureMulti_annotGR<-captureTranscripts[values(captureTranscripts)$gene_id %in% captureMultiexonGeneIds]
captureSingle_annotGR<-captureTranscripts[!values(captureTranscripts)$gene_id %in% captureMultiexonGeneIds]
save(captureSingle_annotGR,file=paste0(robjectsDir,"captureSingle_annotGR.Rdata"))
captureSingle_annotGR_RPKM1=captureSingle_annotGR[values(captureSingle_annotGR)$maxFPKM>=1]



capture_RPKM_none=captureMulti_annotGR[values(captureMulti_annotGR)$maxFPKM<1]
capture_RPKM1=captureMulti_annotGR[values(captureMulti_annotGR)$maxFPKM>=1&values(captureMulti_annotGR)$maxFPKM<10]
capture_RPKM10=captureMulti_annotGR[values(captureMulti_annotGR)$maxFPKM>=10&values(captureMulti_annotGR)$maxFPKM<100]
capture_RPKM100=captureMulti_annotGR[values(captureMulti_annotGR)$maxFPKM>=100]

capture_RPKM=capture_RPKM_none
save(capture_RPKM,file=paste0(robjectsDir,"capture_RPKM_none.Rdata"))
capture_RPKM=capture_RPKM1
save(capture_RPKM,file=paste0(robjectsDir,"capture_RPKM1.Rdata")) 
capture_RPKM=capture_RPKM10
save(capture_RPKM,file=paste0(robjectsDir,"capture_RPKM10.Rdata")) 
capture_RPKM=capture_RPKM100
save(capture_RPKM,file=paste0(robjectsDir,"capture_RPKM100.Rdata")) 

capture_confident=c(capture_RPKM1,capture_RPKM10,capture_RPKM100)
save(capture_confident,file=paste0(robjectsDir,"capture_confident.Rdata")) 
save(captureMulti_annotGR,file=paste0(robjectsDir,"all_capture_transcripts.Rdata"))


############ Full length transcripts

#Take the transcripts, split them by gene
multiexonicRPKM1plus=captureMulti_annotGR[values(captureMulti_annotGR)$maxFPKM>=1]
save(captureMulti_annotGR,file=paste0(robjectsDir,"captureMulti_annotGR.Rdata"))
save(multiexonicRPKM1plus,file=paste0(robjectsDir,"multiexonicRPKM1plus.Rdata"))
genes<-split(captureMulti_annotGR,values(captureMulti_annotGR)$gene_id)
gene_start<-sapply(genes,function(x){min(start(x))})
gene_end<-sapply(genes,function(x){max(end(x))})
gene_chromosome<-sapply(genes,function(x){as.character(seqnames(x[1]))})
gene_strand<-sapply(genes,function(x){as.character(strand(x[1]))})
gene_maxvalue<-sapply(genes,function(x){values(x[1])$maxFPKM})

gr<-GRanges(seqnames=gene_chromosome,IRanges(start=gene_start,end=gene_end),strand=gene_strand,gene_id=names(genes),seqlengths=seqlengths(tl),maxFPKM=gene_maxvalue)
save(gr,file=paste0(robjectsDir,"capture_transcripts_start2end.Rdata"))

#Save singletons that overlap multiexonic regions
mat<-findOverlaps(captureSingle_annotGR,gr)
captureSingle_annotGR<-captureSingle_annotGR[-unique(queryHits(mat))]
save(captureSingle_annotGR,file=paste0(robjectsDir,"captureSingle_all.Rdata"))


#eliminate the multiexonic regions from the singletons
mat<-findOverlaps(captureSingle_annotGR_RPKM1,gr)
captureSingle_annotGR_RPKM1<-captureSingle_annotGR_RPKM1[-unique(queryHits(mat))]
save(captureSingle_annotGR_RPKM1,file=paste0(robjectsDir,"captureSingle_RPKM1.Rdata"))
captureSingle_annotGR_confident=captureSingle_annotGR_RPKM1[values(captureSingle_annotGR_RPKM1)$maxFPKM>=10]
save(captureSingle_annotGR_confident,file=paste0(robjectsDir,"captureSingle_confident.Rdata"))


#take normalized counts into voom to make a linear model fit
v <- voom(normalized.count)
#1. fetch values
voomCounts<-v$E

#2. extract genes from RPKM1-100
RPKMnogenes<-unique(values(capture_RPKM_none)$gene_id)
RKPM1genes<-unique(values(capture_RPKM1)$gene_id)
RKPM10genes<-unique(values(capture_RPKM10)$gene_id)
RKPM100genes<-unique(values(capture_RPKM100)$gene_id)

capture_RPKM=as.data.frame(voomCounts[RPKMnogenes,])
save(capture_RPKM,file=paste0(robjectsDir,"voom_RPKM_none.Rdata"))
capture_RPKM=as.data.frame(voomCounts[RKPM1genes,])
save(capture_RPKM,file=paste0(robjectsDir,"voom_RPKM1.Rdata")) 
capture_RPKM=as.data.frame(voomCounts[RKPM10genes,])
save(capture_RPKM,file=paste0(robjectsDir,"voom_RPKM10.Rdata")) 
capture_RPKM=as.data.frame(voomCounts[RKPM100genes,])
save(capture_RPKM,file=paste0(robjectsDir,"voom_RPKM100.Rdata")) 

x <- as.factor(colnames(alltissue_counts))
set <- newSeqExpressionSet(as.matrix(normalized.count), phenoData = data.frame(x, row.names=colnames(filtered)))
imageFile=paste0("../../project_results/figure1/all_tissues_normalized_",projectname,"_",timeStamp,".pdf")
pdf(imageFile,width=13,height=8)
par(mfrow=c(1,2))
plotRLE(set, outline=FALSE, ylim=c(-4, 4),las=2)
pca<-princomp(normalized.count)
plot(pca$loading,pch=19,col="white")
text(pca$loading,colnames(normalized.count))
dev.off()












