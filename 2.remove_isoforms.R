
library(ShortRead)
library(rtracklayer)
library(ggplot2)
library(reshape2)
library(scales)
#library(RUVSeq)
library(edgeR)
library(limma)
library(stringr)

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
tracksDir = paste(resultsDir,"/gwas_capseq.tracks/",sep="")


#Get the ERCC information to remove before normalization
ercc_transcripts_file=paste0(homedir,"/genomes/rsem/cuffmerged_",projectname,"_gencode/ercc.gtf")
ercc_txs=import(ercc_transcripts_file)
ercc_txs_df=as.data.frame(values(ercc_txs))
ercc_txs_df$ercc<-as.character(seqnames(ercc_txs))
#clean up the isomers that never get more than 1% of isomer count.
#Fetch RSEM output for genes
filesIso<-list.files(inDir,full.names=T)
filesIso<-filesIso[!grepl("pbmc",filesIso)]

resultsCountsIso<-list()
resultsCountsIsoPct<-list()
resultsCountsIsoFile=paste0("../../project_results/transcript_assembly.Rdata/resultsCountsIso_",projectname,".Rdata")
#if(!file.exists(resultsCountsFile)){
	for(sampleFile in filesIso){
	  sampleName=basename(sampleFile)
	  cat(sampleName)
	  cat("\t")
	  fileName<-paste0(sampleFile,"/",sampleName,".isoforms.results")
	  data=read.table(fileName,header=T,sep="\t") 
	  resultsCountsIso[[sampleName]]<-as.integer(data$expected_count)
  	  resultsCountsIsoPct[[sampleName]]<-data$IsoPct

	}
#} else {load(resultsCountsFile)}

#1. First normalize all the isoforms. 
#2. Then eliminate those that do not add more than 1% to isomer count
#3. Correct the gtf files for nice images
#4. Then add this to the rose chart and the heatmap


########## Eliminate those isoforms with less than 1% contribution

alltissue_countsIsoPct=as.data.frame(resultsCountsIsoPct)
alltissue_countsIsoPct=alltissue_countsIsoPct[!grepl("pbmc",names(alltissue_countsIsoPct))]

alltissue_countsIsoPct$gene_id=data$gene_id
alltissue_countsIsoPct$transcript_id=data$transcript_id
alltissue_countsIsoPct=alltissue_countsIsoPct[!alltissue_countsIsoPct$gene_id %in% ercc_txs_df$gene_id,]

#find genes with several equally strong isoforms?
maxIsoPct<-apply(alltissue_countsIsoPct[,c(-22,-23)],1,max)
alltissue_countsIsoPct=alltissue_countsIsoPct[maxIsoPct>=1,]


################ load the data with all isoforms

load(paste0(robjectsDir,"capture_confident.Rdata")) 

capture_confident=capture_confident[capture_confident$transcript_id %in% alltissue_countsIsoPct$transcript_id]

forCuffcompare<-capture_confident
tempCC=as.data.frame(values(forCuffcompare))
tempCC<-tempCC[,-c(11:32)]


genes<-data.frame(ids=sort(unique(capture_confident$gene_id)))
genes$gene_name=paste0("garvan_",str_pad(1:length(genes$ids), 4, pad = "0"))

temp=as.data.frame(values(capture_confident))
merged<-merge(genes,temp,by.x="ids",by.y="gene_id")
mergedCC<-merge(genes,tempCC,by.x="ids",by.y="gene_id")
colnames(mergedCC)[2]<-"gene_id"
values(capture_confident)=merged
values(forCuffcompare)=mergedCC
export.gff2(forCuffcompare,"forCuffcompare.gtf")
#rename transcripts
tempShort=unique(as.data.frame(values(capture_confident)[,c("gene_name","transcript_id")]))
genesL<-split(tempShort,tempShort$gene_name)
out<-lapply(genesL,function(x){x$transcript_id_new<-paste0(x$gene_name,"_",1:length(x[,2]));x})

outDF<-do.call("rbind",out)
colnames(outDF)<-c("gene_name_new","transcript_id","transcript_id_new")
merged<-merge(merged,outDF,by.x="transcript_id",by.y="transcript_id")

merged<-merged[,c(3,35,5,8,12:33)]
values(capture_confident)=merged

colnames(values(capture_confident))[1]="gene_id"
colnames(values(capture_confident))[2]="transcript_id"


######### Add cuffcompare data

data=read.table("../figure1/cuffcmp.forCuffcompare.gtf.tmap",header=T)
code<-unique(data[,c(3,4)])





write.table(as.data.frame(capture_confident),"GWAS_transcripts.txt",row.names=F,quote=F,sep="\t")
export.gff2(capture_confident,"GWAS_transcripts.gtf")

save(capture_confident,file=paste0(robjectsDir,"capture_confident_iso.Rdata")) 
save(captureMulti_annotGR,file=paste0(robjectsDir,"all_capture_transcripts_iso.Rdata"))


############ Full length transcripts

##Take the transcripts, split them by gene
#multiexonicRPKM1plus=captureMulti_annotGR[values(captureMulti_annotGR)$maxFPKM>=1]
#save(captureMulti_annotGR,file=paste0(robjectsDir,"captureMulti_annotGR.Rdata"))
#save(multiexonicRPKM1plus,file=paste0(robjectsDir,"multiexonicRPKM1plus.Rdata"))
genes<-split(capture_confident,values(capture_confident)$gene_id)
gene_start<-sapply(genes,function(x){min(start(x))})
gene_end<-sapply(genes,function(x){max(end(x))})
gene_chromosome<-sapply(genes,function(x){as.character(seqnames(x[1]))})
gene_strand<-sapply(genes,function(x){as.character(strand(x[1]))})
gene_maxvalue<-sapply(genes,function(x){values(x[1])$maxFPKM})

gr<-GRanges(seqnames=gene_chromosome,IRanges(start=gene_start,end=gene_end),strand=gene_strand,gene_id=names(genes),seqlengths=seqlengths(tl),maxFPKM=gene_maxvalue)
save(gr,file=paste0(robjectsDir,"capture_transcripts_start2end_iso.Rdata"))
#
##Save singletons that overlap multiexonic regions
#mat<-findOverlaps(captureSingle_annotGR,gr)
#captureSingle_annotGR<-captureSingle_annotGR[-unique(queryHits(mat))]
#save(captureSingle_annotGR,file=paste0(robjectsDir,"captureSingle_all.Rdata"))



#Clean up the gtfs
load(paste0(robjectsDir,"capture_RPKM1.Rdata"))
capture_RPKM=capture_RPKM[capture_RPKM$transcript_id %in% captureIsoMultiexonTranscriptIds]
export.gff2(capture_RPKM,paste0(tracksDir,"capture_RPKM1.gtf"))

load(paste0(robjectsDir,"capture_RPKM10.Rdata"))
capture_RPKM=capture_RPKM[capture_RPKM$transcript_id %in% captureIsoMultiexonTranscriptIds]
export.gff2(capture_RPKM,paste0(tracksDir,"capture_RPKM10.gtf"))

load(paste0(robjectsDir,"capture_RPKM100.Rdata"))
capture_RPKM=capture_RPKM[capture_RPKM$transcript_id %in% captureIsoMultiexonTranscriptIds]
export.gff2(capture_RPKM,paste0(tracksDir,"capture_RPKM100.gtf"))

load(paste0(robjectsDir,"capture_confident.Rdata"))
capture_RPKM=capture_RPKM[capture_RPKM$transcript_id %in% captureIsoMultiexonTranscriptIds]
export.gff2(capture_RPKM,paste0(tracksDir,"capture_confident.gtf"))

load(paste0(robjectsDir,"all_capture_transcripts.Rdata"))
capture_RPKM=capture_RPKM[capture_RPKM$transcript_id %in% captureIsoMultiexonTranscriptIds]
export.gff2(capture_RPKM,paste0(tracksDir,"all_capture_transcripts.gtf"))

















