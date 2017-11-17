#This script takes the results of bam Robjects and makes FPKMs

library(GenomicRanges)
library(BSgenome)
library("BSgenome.Hsapiens.UCSC.hg19")
library(rtracklayer)
library(ggplot2)
library(RColorBrewer)
library(data.table)
library(gwascat)

chrs=seqlengths(Hsapiens)[!grepl("_",names(seqlengths(Hsapiens)))]
#chrs=seqlengths(Hsapiens)
#annotation_directory
homedir="../../../../"

projectDir =paste0(homedir,"projects/melanoma/project_results/")
robjectsDir = paste(projectDir,"gwas_capseq.Robjects/",sep="")
annotDir = paste0(homedir,"projects/melanoma/annotation/")
tracksDir = paste(annotDir,"captureSpace/",sep="")

#### annotate normGR
#data=as.data.frame(values(normGR))
#merged<-merge(data,rpkms,by.x="gene_id",by.y=0,all.x=T)
#merged$maxFPKM<-apply(merged[,11:21],1,max)
#annotGR<-normGR
#values(annotGR)=merged

#get the multiexons
#multiexons<-table(annotGR$transcript_id)
#multi_annotGR<-annotGR[multiexons>1]

############Find their overlaps with GWAS regions



#All NHGRI
#data(gwrngs)
#gwasGR<-as(gwrngs,"GRanges")
data(gwrngs19)
gwasGR<-as(gwrngs19,"GRanges")
gwasGR<-gwasGR[values(gwasGR)$p.Value<5e-08]
#6266 GWAS 

#Capture space
captureSpaces<-GRangesList()

for(file in list.files(tracksDir,pattern="mergeoverlapping",full.names=T)){
	
	data<-read.table(file,header=F,sep="\t")
	data<-data[!grepl("_",data[,1]),]
	temp<-GRanges(seqnames=as.character(data[,1]),IRanges(start=data[,2],end=data[,3]),seqlengths=seqlengths(gwasGR),type=data[,4])
	captureSpaces[[basename(file)]]<-temp

}

#Annotate NHGRI with capture space
names(captureSpaces)<-c("targetRegion")
targetRegion<-captureSpaces[["targetRegion"]]

types<-as.character(values(targetRegion)$type)
types[grepl("mela",types)]="melanoma"
types[grepl("exon",types)]="pilot_exons"
types[grepl("GWAS_loci",types)]="GWAS_loci"
types[grepl("breast",types)]="breast_cancer"
types[grepl("GWAS_intronic_control",types)]="pilot_introns"

values(targetRegion)$type=types

targetShort<-targetRegion[types %in% c("pilot_exons","pilot_introns","GWAS_loci","melanoma")]
temp<-targetShort
values(targetShort)<-NULL
capture_type<-split(targetShort,values(temp)$type)
#Transcripts are split into genes that cover coding, noncoding, gwas, exon_pilot, intron_pilot 
# and regular gencode annotation is added as a control

#Add gencode regions to the annotation
load("../../annotation/Robjects/gencodeAnnotation.Rdata")
protein_coding=gencodeAnnotation[values(gencodeAnnotation)$gene_type=="protein_coding"&values(gencodeAnnotation)$type=="exon"]
lncRNA=gencodeAnnotation[values(gencodeAnnotation)$gene_type=="lincRNA"&values(gencodeAnnotation)$type=="exon"]
protein_coding_clean<-protein_coding
lncRNA_clean<-lncRNA
values(protein_coding_clean)<-NULL
values(lncRNA_clean)<-NULL

#add the id column to the annotation
values(protein_coding_clean)$id=1:length(protein_coding_clean)
values(lncRNA_clean)$id=1:length(lncRNA_clean)

############### Separate multiexonic transcripts between annotations ##########

#First load the regions
haploblocks=read.table("../../project_results/transcript_assembly.tabFiles/GWAS_nocontig_linked_traits.txt",header=T,sep="\t")
haploblockGR<-GRanges(seqnames=haploblocks$Chr,IRanges(start=haploblocks$LD.region.start,end=haploblocks$LD.region.end),strand="*",seqlengths=seqlengths(capture_type[[1]]))
haploblocksPilot=read.table("../../project_results/transcript_assembly.tabFiles/GWAS_nocontig_linked_traits_pilot.txt",header=T,sep="\t")
haploblockPilotGR<-GRanges(seqnames=haploblocksPilot[,1],IRanges(start=haploblocksPilot[,2],end=haploblocksPilot[,3]),strand="*",seqlengths=seqlengths(capture_type[[1]]))

#Make the regions unique
haploblockGR<-reduce(haploblockGR)
haploblockPilotGR<-reduce(haploblockPilotGR)
combined<-reduce(c(haploblockGR,haploblockPilotGR))
#If we combine the GWAS 296 and pilot 329, there 560 unique Regions. They should be cleaned from introns and identifiable
#Another option is to keep them separate. There is overlap though. Guys should make a decision. I'll make both plots.


#Ensure that after cleaning of the regions the blocks remain identifiable
haploblockGR$id=1:length(haploblockGR)
haploblockPilotGR$id=1:length(haploblockPilotGR)
introns<-capture_type[[3]]
values(introns)$id=1:length(introns)
values(combined)$id=1:length(combined)
melanoma<-capture_type[["melanoma"]]
values(melanoma)$id=1:length(melanoma)

capture_regions<-GRangesList()
capture_regions_short<-GRangesList()
capture_regions[["lower_significance_GWAS"]]<-haploblockPilotGR
capture_regions[["higher_significance_GWAS"]]<-haploblockGR
capture_regions[["melanoma"]]<-melanoma
capture_regions[["pilot_introns"]]<-introns
capture_regions[["protein_coding"]]<-protein_coding_clean
capture_regions[["lncRNA"]]<-lncRNA_clean
strand(capture_regions[[5]])<-"*"
strand(capture_regions[[6]])<-"*"
#Then make sure that the captured regions do not overlap (sequential elimination of lncRNAs, protein, introns, pilot regions from GWAS capture regions)

capture_types<-capture_regions
capture_regions_short[["GWAS_loci"]]<-combined
capture_regions_short[["melanoma"]]<-melanoma
capture_regions_short[["pilot_introns"]]<-introns
capture_regions_short[["protein_coding"]]<-protein_coding_clean
capture_regions_short[["lncRNA"]]<-lncRNA_clean
strand(capture_regions_short[[4]])<-"*"
strand(capture_regions_short[[5]])<-"*"

cleanRegions<-function(regionList){
  types=rev(names(regionList))
  for(type in types){
    cat(type)
    cat("\t")
    typestoClean<-types[!types %in% type]
    for(typetoClean in typestoClean){
      cat("\t")
      cat(typetoClean)
      #setdiff cleans out the id, make it keep it
      diffRegion<-setdiff(regionList[[typetoClean]],regionList[[type]])
      mat<-findOverlaps(diffRegion,regionList[[typetoClean]])
      cat(paste0(" (",sum(width(regionList[[typetoClean]])),") "))
      values(diffRegion)$id="NA"
      values(diffRegion)$id[queryHits(mat)]<-values(regionList[[typetoClean]])$id[subjectHits(mat)]
      regionList[[typetoClean]]<-diffRegion
      if(sum(is.na(values(diffRegion)$id))>0){break()}
    }  
    cat("\n")
  }
  return(regionList)
}
capture_regions<-cleanRegions(capture_regions)
capture_regions_short<-cleanRegions(capture_regions_short)

sapply(capture_types,length)
sapply(capture_regions,length)


melanoma_haploblocks=read.table("../../project_results/transcript_assembly.tabFiles/hg19_melanoma_regions_annotated.bed",header=F,sep="\t")
melanoma_haploblocksGR<-GRanges(seqnames=melanoma_haploblocks[,1],IRanges(start=melanoma_haploblocks[,2],end=melanoma_haploblocks[,3]),strand="*",seqlengths=seqlengths(capture_type[[1]]))
write.table(melanoma_haploblocks[,c(1,2,3,5,4)],"S1D.melanoma_haploblocks.txt",quote=F,row.names=F,sep="\t")
export(melanoma_haploblocksGR,"S1D.melanoma_haploblocks.bed")

#Annotate melanoma haploblocks
mat<-findOverlaps(melanoma_haploblocksGR,capture_regions_short[["melanoma"]])
melanomaL<-split(capture_regions_short[["melanoma"]][subjectHits(mat)],queryHits(mat))
for(i in 1:8){
  melanomaL[[i]]$id=i
}
melanoma_clean<-unlist(melanomaL)

projectname="melanoma"
names(melanoma_clean)=1:length(melanoma_clean)
melanoma_clean_table<-as.data.frame(melanoma_clean)
save(capture_regions,file=paste0("../../project_results/transcript_assembly.Rdata/capture_regions_",projectname,".Rdata"))
save(capture_regions_short,file=paste0("../../project_results/transcript_assembly.Rdata/capture_regions_short_",projectname,".Rdata"))
write.table(melanoma_clean_table,"S1E.melanoma_haploblocks_captured_regions.txt",quote=F,row.names=F,sep="\t")
export(melanoma_clean,"S1E.melanoma_haploblocks_captured_regions.bed")
##################### Find overlaps with regions
