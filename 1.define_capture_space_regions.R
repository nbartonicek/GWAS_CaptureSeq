#This script defines the capture space
#A modification has to be made, since some of the regions were intronic or overlapped known exons


library(GenomicRanges)
library(BSgenome)
library("BSgenome.Hsapiens.UCSC.hg19")
library(rtracklayer)
library(ggplot2)
library(RColorBrewer)
library(data.table)
library(gwascat)

chrs=seqlengths(Hsapiens)[!grepl("_",names(seqlengths(Hsapiens)))]
homedir="../../../../"

projectDir =paste0(homedir,"projects/melanoma/project_results/")
robjectsDir = paste(projectDir,"gwas_capseq.Robjects/",sep="")
annotDir = paste0(homedir,"projects/melanoma/annotation/")
tracksDir = paste(annotDir,"captureSpace/",sep="")

projectname="tissues"

############Find their overlaps with GWAS regions



#All NHGRI
#data(gwrngs)
#gwasGR<-as(gwrngs,"GRanges")
data(gwrngs19)
gwasGR<-as(gwrngs19,"GRanges")
gwasGRSig<-gwasGR[values(gwasGR)$p.Value<=5e-08]
gwasGRSemiSig<-gwasGR[values(gwasGR)$p.Value<=1e-05]

captureSpaces<-GRangesList()
for(file in list.files(tracksDir,pattern="mergeoverlapping",full.names=T)){	
	data<-read.table(file,header=F,sep="\t")
	data<-data[!grepl("_",data[,1]),]
	temp<-GRanges(seqnames=as.character(data[,1]),IRanges(start=data[,2],end=data[,3]),type=data[,4])
	captureSpaces[[basename(file)]]<-temp
}

#Annotate NHGRI with capture space
names(captureSpaces)<-c("targetRegion")
targetRegion<-captureSpaces[["targetRegion"]]

types<-as.character(values(targetRegion)$type)
types[grepl("exon",types)]="pilot_exons"
types[grepl("GWAS_loci",types)]="GWAS_loci"
types[grepl("mela",types)]="melanoma"
types[grepl("breast",types)]="breast_cancer"
types[grepl("GWAS_intronic_control",types)]="pilot_introns"

values(targetRegion)$type=types

targetShort<-targetRegion[types %in% c("pilot_exons","pilot_introns","GWAS_loci")]
temp<-targetShort
values(targetShort)<-NULL
capture_type<-split(targetShort,values(temp)$type)
#Transcripts are split into genes that cover coding, noncoding, gwas, exon_pilot, intron_pilot 
# and regular gencode annotation is added as a control

#Get the gencode v12 genes. Get the regions. See how many of the GWAS capture ones are intronic for protein coding
#Make a table of how many of them are 
gencode12<-import("../../annotation/gencode_archive/gencode.v12.annotation.gtf.gz")
gencode12PC<-gencode12[gencode12$gene_type=="protein_coding"]
if(!file.exists(robjectsDir,"gencodeV12_proteincoding_start2end.Rdata")){
  genes<-split(gencode12PC,values(gencode12PC)$gene_id)
  gene_start<-sapply(genes,function(x){min(start(x))})
  gene_end<-sapply(genes,function(x){max(end(x))})
  gene_chromosome<-sapply(genes,function(x){as.character(seqnames(x[1]))})
  gene_strand<-sapply(genes,function(x){as.character(strand(x[1]))})

  gr<-GRanges(seqnames=gene_chromosome,IRanges(start=gene_start,end=gene_end),strand=gene_strand,gene_id=names(genes))
  save(gr,file=paste0(robjectsDir,"gencodeV12_proteincoding_start2end.Rdata"))
} else {load(paste0(robjectsDir,"gencodeV12_proteincoding_start2end.Rdata"))}
gencode12lncRNA<-gencode12[gencode12$gene_type=="lincRNA"]
if(!file.exists(robjectsDir,"gencodeV12_lncRNA_start2end.Rdata")){
  genes<-split(gencode12lncRNA,values(gencode12lncRNA)$gene_id)
  gene_start<-sapply(genes,function(x){min(start(x))})
  gene_end<-sapply(genes,function(x){max(end(x))})
  gene_chromosome<-sapply(genes,function(x){as.character(seqnames(x[1]))})
  gene_strand<-sapply(genes,function(x){as.character(strand(x[1]))})
  grLncRNA<-GRanges(seqnames=gene_chromosome,IRanges(start=gene_start,end=gene_end),strand=gene_strand,gene_id=names(genes))
  save(grLncRNA,file=paste0(robjectsDir,"gencodeV12_lncRNA_start2end.Rdata"))
} else {load(paste0(robjectsDir,"gencodeV12_lncRNA_start2end.Rdata"))}

#Now load refseq and define the overlaps
refseqExons<-import("../../annotation/captureSpace/Coding_refSeq_exons_12_6_12update_on_UCSC.bed")
refseqGenes<-import("../../annotation/captureSpace/Coding_refSeq_Genes_12_6_12update_on_UCSC.bed")


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
haploblocksPilot=read.table("../../project_results/transcript_assembly.tabFiles/GWAS_nocontig_linked_traits_pilot.txt",header=F,sep="\t")
haploblockPilotGR<-GRanges(seqnames=haploblocksPilot[,1],IRanges(start=haploblocksPilot[,2],end=haploblocksPilot[,3]),strand="*",seqlengths=seqlengths(capture_type[[1]]))

#Find how many of haploblocks and pilot haploblocks overlap a) known annotation b) significant SNPs c) insignificant SNPs
length(haploblockGR)
length(haploblockPilotGR)

#Make the regions unique
haploblockGR<-reduce(haploblockGR)
haploblockPilotGR<-reduce(haploblockPilotGR)
length(haploblockGR)
length(haploblockPilotGR)
sum(countOverlaps(haploblockGR,gwasGR)>0)
sum(countOverlaps(haploblockGR,gwasGRSemiSig)>0)

sum(countOverlaps(haploblockGR,gwasGRSig)>0)
sum(countOverlaps(haploblockPilotGR,gwasGR)>0)
sum(countOverlaps(haploblockPilotGR,gwasGRSig)>0)
sum(countOverlaps(haploblockPilotGR,gwasGRSemiSig)>0)

#Find how many overlap protein coding genes, lncRNAs, refseq: exons and total
genes<-GRangesList()
genes[["gencodePC"]]<-gr
genes[["gencodelncRNA"]]<-grLncRNA
refseqGenesShort<-refseqGenes[,"name"]
colnames(values(refseqGenesShort))<-"gene_id"
genes[["refseq"]]<-refseqGenesShort

exons<-GRangesList()
exons[["gencodePC"]]<-gencode12PC[,NULL]
exons[["gencodelncRNA"]]<-gencode12lncRNA[,NULL]
refseqexonsShort<-refseqExons
values(refseqexonsShort)<-NULL
exons[["refseq"]]<-refseqexonsShort

########## Count the overlaps with annotation
gene_overlap<-sum(countOverlaps(haploblockGR,genes)>0)
gene_overlap
exons_overlap<-sum(countOverlaps(haploblockGR,exons)>0)
exons_overlap

sapply(genes,function(x){sum(countOverlaps(haploblockGR,x)>0)})
sapply(exons,function(x){sum(countOverlaps(haploblockGR,x)>0)})



combined<-reduce(c(haploblockGR,haploblockPilotGR))
#If we combine the GWAS 296 and pilot 329, there 560 unique Regions. They should be cleaned from introns and identifiable
#Another option is to keep them separate. There is overlap though. Guys should make a decision. I'll make both plots.
sum(countOverlaps(combined,gwasGR)>0)
sum(countOverlaps(combined,gwasGRSemiSig)>0)
sum(countOverlaps(combined,gwasGRSig)>0)


#Ensure that after cleaning of the regions the blocks remain identifiable
haploblockGR$id=1:length(haploblockGR)
haploblockPilotGR$id=1:length(haploblockPilotGR)
introns<-capture_type[[3]]
values(introns)$id=1:length(introns)
values(combined)$id=1:length(combined)


capture_regions<-GRangesList()
capture_regions_short<-GRangesList()
capture_regions[["lower_significance_GWAS"]]<-haploblockPilotGR
capture_regions[["higher_significance_GWAS"]]<-haploblockGR
capture_regions[["pilot_introns"]]<-introns
capture_regions[["protein_coding"]]<-protein_coding_clean
capture_regions[["lncRNA"]]<-lncRNA_clean
strand(capture_regions[[4]])<-"*"
strand(capture_regions[[5]])<-"*"
#Then make sure that the captured regions do not overlap (sequential elimination of lncRNAs, protein, introns, pilot regions from GWAS capture regions)

capture_types<-capture_regions
capture_regions_short[["GWAS_loci"]]<-combined
capture_regions_short[["pilot_introns"]]<-introns
capture_regions_short[["protein_coding"]]<-protein_coding_clean
capture_regions_short[["lncRNA"]]<-lncRNA_clean
strand(capture_regions_short[[3]])<-"*"
strand(capture_regions_short[[4]])<-"*"

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


save(capture_regions,file=paste0("../../project_results/transcript_assembly.Rdata/capture_regions_",projectname,".Rdata"))
save(capture_regions_short,file=paste0("../../project_results/transcript_assembly.Rdata/capture_regions_short_",projectname,".Rdata"))

############# Create the Supplementary Table 1. with annotation of diseases and SNPs
############# Columns: chromosome, start, end, SNP(s), disease(s), 
annotated<-combined
mat<-findOverlaps(combined,gwasGRSemiSig)

SNPs<-sapply(regionAnnotationL,function(x){paste(unique(paste0(x$SNPs," (",x$p.Value,",",x$PUBMEDID,")")),collapse=";")})

regionAnnotationL<-split(gwasGRSemiSig[subjectHits(mat)],queryHits(mat))
annotated$SNP<-sapply(regionAnnotationL,function(x){paste(unique(x$SNPs),collapse=";")})
annotated$min.p.value<-sprintf("%g",sapply(regionAnnotationL,function(x){min(x$p.Value)}))
annotated$Disease.Trait<-sapply(regionAnnotationL,function(x){paste(unique(x$Disease.Trait),collapse=";")})
annotated$p.value<-SNPs

df=as.data.frame(annotated)
df=df[,-5]
colnames(df)<-c("Chromosome","Start","End","Width","SNPs","Minimal pvalue","Disease.Trait","SNP(P.value,Pubmed ID)")
write.table(df,file="Supplementary_table1.txt",quote=F,row.names=F,sep="\t")





