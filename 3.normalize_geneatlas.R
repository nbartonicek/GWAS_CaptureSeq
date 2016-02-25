
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
inDirGA =paste0(projectDir,"/project_results/geneatlas.rsem")

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


filesGA<-list.files(inDirGA,full.names=T,pattern="genes")

for(sampleFile in filesGA){
  sampleName=gsub(".genes.results","",basename(sampleFile))
  cat(sampleName)
  cat("\t")
  dataGA=read.table(sampleFile,header=T,sep="\t") 
  resultsCounts[[sampleName]]<-as.integer(dataGA$expected_count)
}


############save(resultsCounts,file=resultsCountsFile)
alltissue_counts=as.data.frame(resultsCounts)
alltissue_counts$gene_id=data$gene_id
alltissue_counts$transcript_id=data$transcript_id.s.

#remove ERCC
alltissue_counts=alltissue_counts[!alltissue_counts$gene_id %in% ercc_txs_df$gene_id,]
alltissue_counts<-alltissue_counts[1:35]
row.names(alltissue_counts)=alltissue_counts[,35]
alltissue_counts=alltissue_counts[,-35]
filter <- apply(alltissue_counts, 1, function(x) length(x[x>10])>=1)
filtered <- alltissue_counts[filter,]

#x <- as.factor(colnames(alltissue_counts))
#set <- newSeqExpressionSet(as.matrix(filtered), phenoData = data.frame(x, row.names=colnames(filtered)))

########## Normalize means by DESeq normalization factors, and dispersion with quantile normalization
library(DESeq)
countTable=filtered
condition=as.factor(rep(c("1", "2", "1",  "3", "4", "5","6"), c(3,1,2,5,5,5,13)))
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
rpkms=rpkm(merged[2:35],gene.length=merged$width, lib.size=mean(colSums(normalized.count)))
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
############save(captureMultiexonGeneIds,file="../../project_results/transcript_assembly.Rdata/captureMultiexonGeneIds_tissues.Rdata")
############save(captureSingleGeneIds,file="../../project_results/transcript_assembly.Rdata/captureSingleGeneIds_tissues.Rdata")

#34797 Total genes
#Total of 11766 genes overlap the capture regions
#Total of 1832 Multiexonic, 9934 Singletons 
sampleNumber=dim(rpkms)[2]
#maxFPKM<-apply(rpkms[,1:(sampleNumber)],1,max)

############save(rpkms,file="../../project_results/transcript_assembly.Rdata/all_rpkms_tissues.Rdata")

total_rpkms=rpkms[row.names(rpkms) %in% captureGeneIds,]
#save(total_rpkms,file="../../project_results/transcript_assembly.Rdata/total_rpkms_tissues.Rdata")
######save(total_rpkms,file="../../project_results/transcript_assembly.Rdata/capture_rpkms_tissues.Rdata")

maxFPKM<-apply(total_rpkms[,1:(sampleNumber)],1,max)

tempRpkms=rpkms
rpkms<-rpkms[row.names(rpkms) %in% captureMultiexonGeneIds,]
#save(rpkms,file="../../project_results/transcript_assembly.Rdata/rpkms_tissues.Rdata")
############save(rpkms,file="../../project_results/transcript_assembly.Rdata/capture_multi_rpkms_tissues.Rdata")

rpkms<-tempRpkms[row.names(tempRpkms) %in% captureSingleGeneIds,]
############save(rpkms,file="../../project_results/transcript_assembly.Rdata/capture_single_rpkms_tissues.Rdata")


#OBJECTS
captureMulti_annotGR<-captureTranscripts[values(captureTranscripts)$gene_id %in% captureMultiexonGeneIds]
captureSingle_annotGR<-captureTranscripts[!values(captureTranscripts)$gene_id %in% captureMultiexonGeneIds]
############save(captureSingle_annotGR,file=paste0(robjectsDir,"captureSingle_annotGR.Rdata"))
captureSingle_annotGR_RPKM1=captureSingle_annotGR[values(captureSingle_annotGR)$maxFPKM>=1]



capture_RPKM_none=captureMulti_annotGR[values(captureMulti_annotGR)$maxFPKM<1]
capture_RPKM1=captureMulti_annotGR[values(captureMulti_annotGR)$maxFPKM>=1&values(captureMulti_annotGR)$maxFPKM<10]
capture_RPKM10=captureMulti_annotGR[values(captureMulti_annotGR)$maxFPKM>=10&values(captureMulti_annotGR)$maxFPKM<100]
capture_RPKM100=captureMulti_annotGR[values(captureMulti_annotGR)$maxFPKM>=100]

capture_RPKM=capture_RPKM_none
############save(capture_RPKM,file=paste0(robjectsDir,"capture_RPKM_none.Rdata"))
capture_RPKM=capture_RPKM1
############save(capture_RPKM,file=paste0(robjectsDir,"capture_RPKM1.Rdata")) 
capture_RPKM=capture_RPKM10
############save(capture_RPKM,file=paste0(robjectsDir,"capture_RPKM10.Rdata")) 
capture_RPKM=capture_RPKM100
############save(capture_RPKM,file=paste0(robjectsDir,"capture_RPKM100.Rdata")) 

capture_confident=c(capture_RPKM1,capture_RPKM10,capture_RPKM100)
############save(capture_confident,file=paste0(robjectsDir,"capture_confident.Rdata")) 
############save(captureMulti_annotGR,file=paste0(robjectsDir,"all_capture_transcripts.Rdata"))

geneatlasAnno<-read.table("E-MTAB-513.sdrf.txt",header=T,sep="\t")
conversion<-unique(geneatlasAnno[,c("Comment.ENA_RUN.","Characteristics.organism.part.")])

forPlot<-tempRpkms

oldcolnames<-colnames(forPlot)
conversionNew<-conversion[conversion$Comment.ENA_RUN. %in% oldcolnames,]
conversionNew<-conversionNew[order(conversionNew$Comment.ENA_RUN.),]
colnames(forPlot)[22:34]=as.character(conversionNew$Characteristics.organism.part.)
colnames(forPlot)[25]<-"skmusc"
colnames(forPlot)[22]<-"testes"
colnames(forPlot)[22:34]<-paste0(colnames(forPlot)[22:34],"_bodyatlas")

############# Separate multiexonic transcripts between annotations ##########


#1. There will be 3. categories
#a) captured, non-captured protein, non-captured lncRNA

#First load the regions
haploblocks=read.table("../../project_results/transcript_assembly.tabFiles/GWAS_nocontig_linked_traits.txt",header=T,sep="\t")
haploblockGR<-GRanges(seqnames=haploblocks$Chr,IRanges(start=haploblocks$LD.region.start,end=haploblocks$LD.region.end),strand="*",seqlengths=seqlengths(annotGR))
haploblocksPilot=read.table("../../project_results/transcript_assembly.tabFiles/GWAS_nocontig_linked_traits_pilot.txt",header=F,sep="\t")
haploblockPilotGR<-GRanges(seqnames=haploblocksPilot[,1],IRanges(start=haploblocksPilot[,2],end=haploblocksPilot[,3]),strand="*",seqlengths=seqlengths(annotGR))
#Make the regions unique
haploblockGR<-reduce(c(haploblockGR,haploblockPilotGR))



#Extract all capture
final_regions=haploblockGR
mat<-findOverlaps(annotGR,final_regions)
#make the total, single and multi objects
captureTranscripts<-annotGR[unique(queryHits(mat))]
captureGeneIdsCapture=unique(values(captureTranscripts)$gene_id)

#### Get all the captured
rpkmsDF<-as.data.frame(forPlot)
idsCapture<-row.names(rpkmsDF)[row.names(rpkmsDF) %in% captureGeneIdsCapture]

#### Get all the noncaptured
idsNC<-row.names(rpkmsDF)[!row.names(rpkmsDF) %in% captureGeneIdsCapture]


#Extract protein coding ids
final_regions=capture_regions_short[[3]]
mat<-findOverlaps(annotGR,final_regions)
#make the total, single and multi objects
captureTranscripts<-annotGR[unique(queryHits(mat))]
capturePCGeneIds=unique(values(captureTranscripts)$gene_id)

final_regions=capture_regions_short[[4]]
mat<-findOverlaps(annotGR,final_regions)
#make the total, single and multi objects
captureTranscripts<-annotGR[unique(queryHits(mat))]
captureNCGeneIds=unique(values(captureTranscripts)$gene_id)


##### Protein coding
rpkmsNC=rpkmsDF[!row.names(rpkmsDF) %in% captureGeneIdsCapture,]
idsNC_PC<-row.names(rpkmsNC)[row.names(rpkmsNC) %in% capturePCGeneIds]
idsNC_NC<-row.names(rpkmsNC)[row.names(rpkmsNC) %in% captureNCGeneIds]


df<-data.frame(id=row.names(rpkmsDF),capture=rpkmsDF$testes,bodyatlas=rpkmsDF$testes_bodyatlas)
#Turn the dataset for ggplot
#samples<-c("ovary","skmusc","prostate","lung","brain","breast","kidney","heart","liver")
#for(sample in samples){
#	cat(sample)
#	cat("\n")
#	temp<-data.frame(id=row.names(rpkmsDF),capture=rpkmsDF[,sample],bodyatlas=rpkmsDF[,paste0(sample,"_bodyatlas")])
#	df<-rbind(df,temp)
#}
df$capture<-log10(df$capture+1)
df$bodyatlas<-log10(df$bodyatlas+1)

dfGencode<-df[df$id %in% idsNC,]

figureName=paste0("../../project_results/figure1/bodyatlas_gencode.pdf")
pdf(figureName,width=8,height=8)
p<-ggplot(dfGencode,aes(bodyatlas,capture)) 
p<-p+geom_point(size=3,alpha=0.5,col='black')
p<-p+xlim(c(0,6)) 
p<-p+ylim(c(0,6)) 
p<-p+theme_bw() +
  theme(plot.background = element_blank()
        ,panel.border = element_blank()
        ,axis.text.x  = element_text(size=12)
        ,axis.text.y  = element_text(size=12)
  )
 p
dev.off()

dfCapture<-df[df$id %in% idsCapture,]
figureName=paste0("../../project_results/figure1/bodyatlas_capture.pdf")
pdf(figureName,width=8,height=8)
p<-ggplot(dfCapture,aes(bodyatlas,capture)) 
p<-p+geom_point(size=3,alpha=0.5,,col='darkred') 
p<-p+xlim(c(0,6)) 
p<-p+ylim(c(0,6)) 
p<-p+theme_bw() +
  theme(plot.background = element_blank()
        ,panel.border = element_blank()
        ,axis.text.x  = element_text(size=12)
        ,axis.text.y  = element_text(size=12)
  )
 p
dev.off()


#Now, for each 
pdf("test_all.pdf",width=8,height=8)
smoothScatter(dfGencode$capture~dfGencode$bodyatlas,colramp = colorRampPalette(c(blues9,"white")),xlim=c(0,6),ylim=c(0,6),col="darkblue",transformation=function(x){x*1000})
#points(log10(rpkmsDF$testes[row.names(rpkmsDF) %in% idsCapture]+1)~log10(rpkmsDF$testes_bodyatlas[row.names(rpkmsDF) %in% idsCapture]+1),xlim=c(0,6),ylim=c(0,6),col="red",pch=".")
dev.off()




#Now, for each 
pdf("test.pdf",width=8,height=8)
rpkmsDF<-as.data.frame(forPlot)
smoothScatter(log10(rpkmsDF$testes[row.names(rpkmsDF) %in% idsNC]+1)~log10(rpkmsDF$testes_bodyatlas[row.names(rpkmsDF) %in% idsNC]+1),xlim=c(0,6),ylim=c(0,6),col="darkblue")
points(log10(rpkmsDF$testes[row.names(rpkmsDF) %in% idsCapture]+1)~log10(rpkmsDF$testes_bodyatlas[row.names(rpkmsDF) %in% idsCapture]+1),xlim=c(0,6),ylim=c(0,6),col="red",pch=".")
dev.off()








