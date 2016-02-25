library(DESeq)
library(ShortRead)
library(rtracklayer)
library(ggplot2)
library(reshape2)
library(scales)
library(RUVSeq)
library(edgeR)
library(limma)
library(gplots)
library(RColorBrewer)


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
 

#Get the ERCC information to remove before normalization
ercc_transcripts_file=paste0(homedir,"/genomes/rsem/cuffmerged_",projectname,"_gencode/ercc.gtf")
ercc_txs=import(ercc_transcripts_file)
ercc_txs_df=as.data.frame(values(ercc_txs))
ercc_txs_df$ercc<-as.character(seqnames(ercc_txs))

############### NORMALIZATION ##################

#Fetch RSEM output for genes
files<-list.files(inDir,pattern="genes",full.names=T)
files<-files[!grepl("pbmc",files)]

resultsCounts<-list()
results<-list()

if(!file.exists(paste0(robjectsDir,"raw_counts_melanoma.Rdata"))){
  for(sampleFile in files){
    sampleName=basename(sampleFile)
    sampleName=gsub("\\..*","",sampleName)
    if(grepl("Melanoma",sampleName)){
       sampleName=gsub(".*Melanoma_","melanoma.",sampleName)
       sampleName=gsub("_.*","",sampleName)
    }
    cat(sampleName)
    cat("\t")
    data=read.table(sampleFile,header=T,sep="\t") 
    resultsCounts[[sampleName]]<-as.integer(data$expected_count)
    results[[sampleName]]<-data$FPKM
  }
  save(results,file=paste0(robjectsDir,"raw_counts_melanoma.Rdata"))
  save(resultsCounts,file=paste0(robjectsDir,"sum_raw_counts_melanoma.Rdata"))
  save(data,file=paste0(robjectsDir,"sample_counts_melanoma.Rdata"))

} else {
  load(paste0(robjectsDir,"raw_counts_melanoma.Rdata"))
  load(paste0(robjectsDir,"sum_raw_counts_melanoma.Rdata"))
  load(paste0(robjectsDir,"sample_counts_melanoma.Rdata"))
}

alltissue_counts=as.data.frame(resultsCounts)
alltissue_counts=alltissue_counts[!grepl("pbmc",names(alltissue_counts))]
alltissue_counts$gene_id=data$gene_id
alltissue_counts$transcript_id=data$transcript_id.s.

alltissue_fpkms=as.data.frame(results)
alltissue_fpkms$gene_id=data$gene_id
alltissue_fpkms$transcript_id=data$transcript_id.s.
alltissue_fpkms$width=data$effective_length
save(alltissue_fpkms,file=paste0("../../project_results/transcript_assembly.Rdata/RSEM_",projectname,".Rdata"))

#remove ERCC
alltissue_counts=alltissue_counts[!alltissue_counts$gene_id %in% ercc_txs_df$gene_id,]
alltissue_counts<-alltissue_counts[1:35]
row.names(alltissue_counts)=alltissue_counts[,35]
alltissue_counts=alltissue_counts[,-35]
#take more than 10 reads
filter <- apply(alltissue_counts, 1, function(x) length(x[x>10])>=1)
filtered <- alltissue_counts[filter,]
filteredNew <- filtered[,c(1:11,25:34,12:24)]
filtered=filteredNew
#x <- as.factor(colnames(alltissue_counts))
#set <- newSeqExpressionSet(as.matrix(filtered), phenoData = data.frame(x, row.names=colnames(filtered)))

########## Normalize means by DESeq normalization factors, and dispersion with quantile normalization
countTable=filtered
condition=as.factor(rep(c("1", "2", "1",  "3", "4", "5", "melanoma"), c(3,1,2,5,5,5,13)))
cds = newCountDataSet( countTable, condition )
cds = estimateSizeFactors( cds )
cds = estimateDispersions( cds, method="blind", sharingMode="fit-only" )
normalized.count <- counts(cds, normalized = TRUE)
normalized.count<-normalizeQuantiles(normalized.count)


#differential expression: only looking for clear examples where melanoma tissues are different
if(!file.exists(paste0(robjectsDir,"melanoma_res.Rdata"))){
  conditionDiff=as.factor(rep(c("tissues", "melanoma"), c(21,13)))
  cds = newCountDataSet( countTable, conditionDiff )
  cds = estimateSizeFactors( cds )
  cds = estimateDispersions( cds, method="blind")
  res = nbinomTest( cds, "tissues", "melanoma" )
  save(res,file=paste0(robjectsDir,"melanoma_res.Rdata"))
} else {load(paste0(robjectsDir,"melanoma_res.Rdata"))}
resSig = res[ res$padj < 0.1, ]
#head( resSig[ order(resSig$pval), ],n=50)

#1. reads overlapping melanoma GWAS
load(paste0("../../project_results/transcript_assembly.Rdata/capture_regions_short_tissues.Rdata"))
final_regions=capture_regions_short[[1]]

melanomaGR<-GRanges(seqnames="chr10",IRanges(start=107516352,end=107516352),seqlengths=seqlengths(final_regions))

mat=findOverlaps(melanomaGR,final_regions)
melanomaGWAS=final_regions[unique(subjectHits(mat))]
#Find trancsripts that overlap this region

 

################ Make transcripts for testing ##################

#get width to calculate RPKMs
gene_width<-alltissue_fpkms[,c("gene_id","width")]
merged<-merge(normalized.count,gene_width,by.x=0,by.y="gene_id")
rpkms=rpkm(merged[2:35],gene.length=merged$width, lib.size=mean(colSums(normalized.count)))
row.names(rpkms)=merged$Row.names
save(rpkms,file=paste0(robjectsDir,"rpkms_melanoma.Rdata"))
rpkms<-as.data.frame(rpkms)
rpkms$maxFPKM<-apply(rpkms[,22:34],1,max)


######### Load chromosomal locations ###############
load(paste0("../../project_results/transcript_assembly.Rdata/merged_gtf_",projectname,".Rdata"))

######### Get normalized counts and coords ##############
normGR<-tl[values(tl)$gene_id %in% row.names(rpkms)]


#### annotate normGR
sampleNumber=dim(rpkms)[2]

data=as.data.frame(values(normGR))
merged<-merge(data,rpkms,by.x="gene_id",by.y=0,all.x=T)
merged$maxFPKM<-apply(merged[,15:(15+sampleNumber-1)],1,max)
annotGR<-normGR
values(annotGR)=merged
save(annotGR,file=paste0("../../project_results/transcript_assembly.Rdata/annotGR_",projectname,".Rdata"))


##################################################################
##################################################################
##################################################################

#There are three analyses that have to be done
#First investigate whether the intergenic GWAS locus associated with melanoma contains anything
#There are no transcripts in that haploblock
#Second: investigate only melanoma-specific regions


##################### Find overlap with Melanoma GWAS SNP intergenic region #######################


final_regions=capture_regions_short[[1]]
melanomaGR<-GRanges(seqnames="chr10",IRanges(start=107516352,end=107516352),seqlengths=seqlengths(final_regions))
mat=findOverlaps(melanomaGR,final_regions)
melanomaGWAS=final_regions[unique(subjectHits(mat))]
mat<-findOverlaps(annotGR,melanomaGWAS)
melanomaGWAStranscripts=annotGR[unique(queryHits(mat))]
save(melanomaGWAStranscripts,file=paste0("../../project_results/transcript_assembly.Rdata/melanomaGWAStranscripts_",projectname,".Rdata"))

#plot the profile
data=as.data.frame(values(melanomaGWAStranscripts))[,c(1,15:(15+sampleNumber-1))]
data=unique(data)

dataM<-melt(data)
names(dataM)=c("gene_id","tissue","FPKM")
imageFile=paste0("../../project_results/figure5/melanomaGWAS_",timeStamp,".pdf")
pdf(imageFile,width=11,height=8)
p<-ggplot(dataM,aes(tissue,FPKM,group=gene_id, colour = gene_id))
p<-p+geom_line()
p<-p+ylab("FPKM")
p<-p+xlab("Tissues")
p<-p+theme(text = element_text(size=20))
p<-p+theme_bw() +
  theme(plot.background = element_blank()
   ,panel.border = element_blank()
   ,axis.text.x=element_text(angle=90)
  )
p
dev.off()



############### find melanoma-regions multiexonic genes ##############
load(paste0("../../project_results/transcript_assembly.Rdata/capture_regions_short_melanoma.Rdata"))


#melanoma_haploblocks_clean<-import("../captureSpace/S1E.melanoma_haploblocks_captured_regions.bed")


final_regions=c(capture_regions_short[[2]])
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
save(captureMultiexonGeneIds,file=paste0("../../project_results/transcript_assembly.Rdata/captureMultiexonGeneIds_",projectname,".Rdata"))

multiexonicMelanoma<-annotGR[annotGR$gene_id %in% captureMultiexonGeneIds]
multiexonicMelanoma<-multiexonicMelanoma[multiexonicMelanoma$gene_id %in% row.names(rpkms[rpkms$maxFPKM>=1,])]

############## The number of novel genes is:
length(unique(multiexonicMelanoma$gene_id))
#get the start2end of the transcripts
genes<-split(multiexonicMelanoma,values(multiexonicMelanoma)$gene_id)
gene_start<-sapply(genes,function(x){min(start(x))})
gene_end<-sapply(genes,function(x){max(end(x))})
gene_chromosome<-sapply(genes,function(x){as.character(seqnames(x[1]))})
gene_strand<-sapply(genes,function(x){as.character(strand(x[1]))})
gene_maxvalue<-sapply(genes,function(x){values(x[1])$maxFPKM})

gr<-GRanges(seqnames=gene_chromosome,IRanges(start=gene_start,end=gene_end),strand=gene_strand,gene_id=names(genes),seqlengths=seqlengths(tl),maxFPKM=gene_maxvalue)
save(gr,file=paste0(robjectsDir,"melanoma_DEmelanomacapture_multiexonic_start2end.Rdata"))



###################  plot the profile of 78 genes #####################
data=as.data.frame(values(multiexonicMelanoma))[,c(1,15:(15+sampleNumber-1))]
data=unique(data)

row.names(data)=data[,1]
data=data[,-1]
data=data[,-35]

data=as.matrix(data)
dataLog=log10(data+1)
imageFile=paste0("../../project_results/figure5/melanomaHeatmap_",timeStamp,".pdf")
pdf(imageFile,width=11,height=15)
palette<-brewer.pal(11,"RdBu")
heatmap.2(as.matrix(dataLog),col=rev(palette),
  trace="none",cexRow=0.5,cexCol=1.0,
  lmat=rbind(c(0, 3, 0), c(2,1,0),c(0,0,0),c(0,4,0)), 
  lhei=c(0.01, 0.08, 0.005,0.018 ),lwid=c(0.1, 0.8, 0.5), density.info = "none", Rowv=T, Colv=T)
dev.off()


############## The number of novel, DE genes is:
ids<-row.names(rpkms[rpkms$maxFPKM>=1,])
resSigShort<-resSig[resSig$id %in% ids,]
DEids<-resSigShort$id
melanomaDE<-multiexonicMelanoma[multiexonicMelanoma$gene_id %in% DEids]
multiexonicMelanomaDEgwas=melanomaDE
length(unique(melanomaDE$gene_id))
export.gff2(multiexonicMelanomaDEgwas,"melanoma_regions_multiexonic_DE.gtf")
save(multiexonicMelanomaDEgwas,file=paste0("../../project_results/transcript_assembly.Rdata/multiexonicMelanomaDEgwas",projectname,".Rdata"))
finalIDs<-unique(multiexonicMelanomaDEgwas$gene_id)

############### find all multiexonic genes that are differentially expressed ##############

final_regions=c(capture_regions_short[[1]])
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
#save(captureMultiexonGeneIds,file=paste0("../../project_results/transcript_assembly.Rdata/captureMultiexonGeneIds_",projectname,".Rdata"))


#get the number of genes that have FPKM>1 in melanoma
ids<-row.names(rpkms[rpkms$maxFPKM>=1,])

#Get the differentially expressed genes
resSigShort<-resSig[resSig$id %in% ids,]
DEids<-resSigShort$id
melanomaDEgwas<-annotGR[annotGR$gene_id %in% resSig$id]
multiexonicMelanomaDEgwas<-melanomaDEgwas[melanomaDEgwas$gene_id %in% captureMultiexonGeneIds]
finalIDs<-unique(multiexonicMelanomaDEgwas$gene_id)
save(multiexonicMelanomaDEgwas,file=paste0("../../project_results/transcript_assembly.Rdata/multiexonicMelanomaDEgwas",projectname,".Rdata"))

break()
#plot the profile
data=as.data.frame(values(multiexonicMelanomaDEgwas))[,c(1,15:(15+sampleNumber-1))]
data=unique(data)
#ids<-data[,c("gene_id","maxFPKM")]
#annot<-unique(values(annotGR)[,c(1,8)])
#annot<-annot[!is.na(annot$gene_name),]
#annot$gene_name=gsub("\\..*","",annot$gene_name)
#annot=unique(annot)
#merged<-merge(ids,annot,by.x="gene_id",by.y="gene_id")
#mergedL<-split(merged,merged$gene_id)
#gene_names<-sapply(mergedL,function(x){x$gene_name[1]})
#names(gene_names)<-names(mergedL)

row.names(data)=data[,1]
data=data[,-1]
data=data[,-35]
data=as.matrix(data)
dataLog=log10(data+1)
imageFile=paste0("../../project_results/figure5/melanomaHeatmapCaptureOnly_",timeStamp,".pdf")
pdf(imageFile,width=11,height=15)
palette<-brewer.pal(11,"RdBu")
heatmap.2(as.matrix(dataLog),col=rev(palette),
  trace="none",cexRow=0.5,cexCol=1.0,
  lmat=rbind(c(0, 3, 0), c(2,1,0),c(0,0,0),c(0,4,0)), 
  lhei=c(0.01, 0.08, 0.005,0.018 ),lwid=c(0.1, 0.8, 0.5), density.info = "none", Rowv=T, Colv=T)
dev.off()

#plot the protein_coding ones
PCids<-DEids[!grepl("XLOC","",DEids)]
melanomaDEgwas<-annotGR[annotGR$gene_id %in% PCids]
multiexonicMelanomaDEgwas<-melanomaDEgwas[melanomaDEgwas$gene_id %in% captureMultiexonGeneIds]
finalIDsPC<-unique(multiexonicMelanomaDEgwas$gene_id)
#save(multiexonicMelanomaDEgwas,file=paste0("../../project_results/transcript_assembly.Rdata/multiexonicMelanomaDEgwas",projectname,".Rdata"))

break()
#plot the profile
data=as.data.frame(values(multiexonicMelanomaDEgwas))[,c(1,15:(15+sampleNumber-1))]
data=unique(data)
#ids<-data[,c("gene_id","maxFPKM")]
#annot<-unique(values(annotGR)[,c(1,8)])
#annot<-annot[!is.na(annot$gene_name),]
#annot$gene_name=gsub("\\..*","",annot$gene_name)
#annot=unique(annot)
#merged<-merge(ids,annot,by.x="gene_id",by.y="gene_id")
#mergedL<-split(merged,merged$gene_id)
#gene_names<-sapply(mergedL,function(x){x$gene_name[1]})
#names(gene_names)<-names(mergedL)

row.names(data)=data[,1]
data=data[,-1]
data=data[,-35]
data=as.matrix(data)
dataLog=log10(data+1)
imageFile=paste0("../../project_results/figure5/melanomaHeatmapPC_",timeStamp,".pdf")
pdf(imageFile,width=11,height=15)
palette<-brewer.pal(11,"RdBu")
heatmap.2(as.matrix(dataLog),col=rev(palette),
  trace="none",cexRow=0.5,cexCol=1.0,
  lmat=rbind(c(0, 3, 0), c(2,1,0),c(0,0,0),c(0,4,0)), 
  lhei=c(0.01, 0.08, 0.005,0.018 ),lwid=c(0.1, 0.8, 0.5), density.info = "none", Rowv=T, Colv=T)
dev.off()












#Get the multiexonic regions first
annotDE<-annotGR[annotGR$gene_id %in% DEids]
MEtx<-unique(annotDE$transcript_id[table(annotDE$transcript_id)>1])
MEgenes<-annotDE$gene_id[annotDE$transcript_id %in% MEtx]
annotDE_ME<-annotDE[annotDE$gene_id %in% MEgenes]

export.gff2(annotDE_ME,"all_melanoma_novel.gtf")




#get the number of genes that have FPKM>1 in melanoma
ids<-row.names(rpkms[rpkms$maxFPKM>=1,])
resSigShort<-resSig[resSig$id %in% ids,]
DEids<-resSigShort$id

#find the pathways for noncapture
noncapture=DEids[!grepl("XLOC",DEids)]
write.table(noncapture,"noncapture.txt",row.names=F,quote=F,sep="\t")


#Load the capture regions
load(paste0("../../project_results/transcript_assembly.Rdata/capture_regions_",projectname,".Rdata"))
load(paste0("../../project_results/transcript_assembly.Rdata/capture_regions_short_",projectname,".Rdata"))
#load(paste0("../../project_results/transcript_assembly.Rdata/capture_regions_tissues.Rdata"))
#load(paste0("../../project_results/transcript_assembly.Rdata/capture_regions_short_tissues.Rdata"))


##################### Find overlap with Capture gwas ################

#Instead of using the n=50 I will take the multiexonic regions (caluclated in the section above)
#data=head( resSig[ order(resSig$pval), ],n=50)

data=resSig[resSig$id %in% captureMultiexonGeneIds, ]

data=data[data$foldChange>1,]
ids=data$id

gr=annotGR[annotGR$gene_id %in% ids]
save(gr,file=paste0(robjectsDir,"melanoma_diff.Rdata"))

melanomaTranscripts=annotGR[annotGR$gene_id %in% ids]
melanomaTranscripts2=annotGR[annotGR$gene_name %in% ids]
melanomaTranscripts=c(melanomaTranscripts,melanomaTranscripts2)

data=as.data.frame(values(melanomaTranscripts))[,c(1,15:(15+sampleNumber-1))]
data=unique(data)
data[,2:35]=log10(data[,2:35]+1)
dataM<-melt(data)
names(dataM)=c("gene_id","tissue","FPKM")
imageFile=paste0("../../project_results/figure5/melanomaCaptureDiff_",timeStamp,".pdf")
pdf(imageFile,width=11,height=8)
p<-ggplot(dataM,aes(tissue,FPKM,group=gene_id, colour = gene_id))
p<-p+geom_line()
p<-p+ylab("Log 10 FPKM")
p<-p+xlab("Tissues")
p<-p+theme(text = element_text(size=20))
p<-p+theme_bw() +
  theme(plot.background = element_blank()
   ,panel.border = element_blank()
   ,axis.text.x=element_text(angle=90)
  )
p
dev.off()



##################### Find overlap with targeted Melanoma capture space #######################

final_regions=capture_regions_short[[2]]
annotGRshort=annotGR[annotGR$maxFPKM>1]
mat<-findOverlaps(annotGRshort,final_regions)

melanomaTranscripts=annotGRshort[unique(queryHits(mat))]
#plot the profile
data=as.data.frame(values(melanomaTranscripts))[,c(1,15:(15+sampleNumber-1))]
data=unique(data)

#plot the profile
data=as.data.frame(values(melanomaTranscripts))[,c(1,15:(15+sampleNumber-1))]
data=unique(data)

row.names(data)=data[,1]
data=data[,-1]
data=as.matrix(data)
dataLog=log10(data+1)
imageFile=paste0("../../project_results/figure5/melanomaHeatmap_",timeStamp,".pdf")
pdf(imageFile,width=11,height=15)
palette<-brewer.pal(11,"RdBu")
heatmap.2(as.matrix(dataLog),col=rev(palette),
  trace="none",cexRow=0.5,cexCol=1.0,
  lmat=rbind(c(0, 3, 0), c(2,1,0),c(0,0,0),c(0,4,0)), 
  lhei=c(0.01, 0.08, 0.005,0.018 ),lwid=c(0.1, 0.8, 0.5), density.info = "none", Rowv=T, Colv=T)
dev.off()
















row.names(data)=data[,1]
data=data[,-1]
dataLog=as.matrix(data)
imageFile=paste0("../../project_results/figure5/melanomaCaptureMelanomaHeatmap_",timeStamp,".pdf")
pdf(imageFile,width=11,height=15)
palette<-brewer.pal(11,"RdBu")
cols=rev(palette)
cols=cols[c(1,2,7:11)]
breaks=c(0,0.5,1,5,10,50,100,100000)
breakslog=log10(breaks+1)

heatmap.2(as.matrix(dataLog),col=cols,
  trace="none",cexRow=0.5,cexCol=1.0,
  breaks = breakslog, 
  lmat=rbind(c(0, 3, 0), c(2,1,0),c(0,0,0),c(0,4,0)), 
  lhei=c(0.01, 0.08, 0.005,0.018 ),lwid=c(0.1, 0.8, 0.5), density.info = "none", Rowv=T, Colv=T)
dev.off()







#########################

final_regions=c(capture_regions_short[[1]],capture_regions_short[[2]])
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
save(captureMultiexonGeneIds,file=paste0("../../project_results/transcript_assembly.Rdata/captureMultiexonGeneIds_",projectname,".Rdata"))



##################### Find overlap with all capture regions that are significantly diff. from tissues ############

















#34797 Total genes
#Total of 11766 genes overlap the capture regions
#Total of 1832 Multiexonic, 9934 Singletons 
sampleNumber=dim(rpkms)[2]
#maxFPKM<-apply(rpkms[,1:(sampleNumber)],1,max)

save(rpkms,file=paste0("../../project_results/transcript_assembly.Rdata/all_rpkms_",projectname,".Rdata"))

total_rpkms=rpkms[row.names(rpkms) %in% captureGeneIds,]
#save(total_rpkms,file="../../project_results/transcript_assembly.Rdata/total_rpkms_tissues.Rdata")
save(total_rpkms,file=paste0("../../project_results/transcript_assembly.Rdata/capture_rpkms_melanoma.Rdata"))













