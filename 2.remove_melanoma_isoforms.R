
library(ShortRead)
library(rtracklayer)
library(ggplot2)
library(reshape2)
library(scales)
#library(RUVSeq)
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
projectname="melanoma_combined"
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
filesIso<-list.files(inDir,full.names=T,pattern="isoforms")
filesIso<-filesIso[grepl("Melanoma",filesIso)]

resultsCountsIso<-list()
resultsCountsIsoPct<-list()
resultsCountsIsoFile=paste0("../../project_results/transcript_assembly.Rdata/resultsCountsIso_",projectname,".Rdata")
#if(!file.exists(resultsCountsFile)){
	for(sampleFile in filesIso){
	  sampleName=basename(sampleFile)
	  cat(sampleName)
	  cat("\t")
	  #fileName<-paste0(sampleFile,"/",sampleName,".isoforms.results")
	  fileName<-sampleFile
	  data=read.table(fileName,header=T,sep="\t") 
	  resultsCountsIso[[sampleName]]<-as.integer(data$expected_count)
  	  resultsCountsIsoPct[[sampleName]]<-data$IsoPct

	}
#} else {load(resultsCountsFile)}
alltissue_IsoPct=as.data.frame(resultsCountsIsoPct)
row.names(alltissue_IsoPct)=data$transcript_id
alltissue_IsoPct$maxIso<-unlist(apply(alltissue_IsoPct,1,max))
alltissue_IsoPct_short<-alltissue_IsoPct[alltissue_IsoPct$maxIso>1,]

#Simply get the isoforms that are never below 1%
#Eliminate them from the tl or annotGR 
load(paste0("../../project_results/transcript_assembly.Rdata/merged_gtf_",projectname,".Rdata"))
load(paste0("../../project_results/transcript_assembly.Rdata/annotGR_",projectname,".Rdata"))
isoClean<-annotGR[annotGR$transcript_id %in% row.names(alltissue_IsoPct_short)]

#eliminate single and small maxFPKM
annotGR_confident=annotGR[annotGR$maxFPKM>1]
idList<-table(annotGR_confident$transcript_id)
txList=idList[idList>1]

#multiAnnotGR
multiAnnotGR<-annotGR_confident[annotGR_confident$transcript_id %in% names(txList)]
multiAnnotGR$score=1

export.gff2(multiAnnotGR,"/share/ClusterShare/biodata/contrib/nenbar/projects/melanoma/project_results/cuffcompare/iso_melanoma/melanoma_isoclean.gtf")

