

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("methylKit")
BiocManager::install("genomation")
BiocManager::install("rtracklayer")

library("methylKit")
library("genomation")
library("rtracklayer")

setwd("/home/danielruzzante/BrookTroutEpigenomics2024/methyldata")


file.list = list( file.path( "23_CV_SFO_005_S24_CpG.txt"),
                  file.path("23_CV_SFO_010_S18_CpG.txt"),
                  file.path("23_CV_SFO_007_S25_CpG.txt"),
                  file.path("23_CV_SFO_012_S19_CpG.txt"),
                  file.path( "23_CV_SFO_014_S20_CpG.txt"),
                  file.path( "23_CV_SFO_015_S21_CpG.txt"),
                  file.path( "23_CV_SFO_017_S26_CpG.txt"),
                  file.path("23_CV_SFO_024_S22_CpG.txt"),
                  file.path("23_CV_SFO_025_S23_CpG.txt"),
                  file.path( "23_RCD_SFO_001_S1_CpG.txt"),
                  file.path("23_RCD_SFO_007_S2_CpG.txt"),
                  file.path("23_RCD_SFO_009_S3_CpG.txt"),
                  file.path("23_RCD_SFO_012_S4_CpG.txt"),
                  file.path( "23_RCD_SFO_015_S5_CpG.txt"),
                  file.path("23_RCD_SFO_017_S6_CpG.txt"),
                  file.path("23_RCD_SFO_021_S7_CpG.txt"),
                  file.path( "23_RCD_SFO_023_S8_CpG.txt"),
                  file.path( "23_RCD_SFO_024_S9_CpG.txt"),
                  file.path("23_WWD_SFO_003_S10_CpG.txt"),
                  file.path( "23_WWD_SFO_004_S11_CpG.txt"),
                  file.path( "23_WWD_SFO_008_S12_CpG.txt"),
                  file.path( "23_WWD_SFO_015_S13_CpG.txt"),
                  file.path( "23_WWD_SFO_016_S27_CpG.txt"),
                  file.path("23_WWD_SFO_018_S14_CpG.txt"),
                  file.path( "23_WWD_SFO_021_S15_CpG.txt"),
                  file.path("23_WWD_SFO_023_S16_CpG.txt"),
                  file.path("23_WWD_SFO_024_S17_CpG.txt"))
                  
                  myobj = methRead(file.list,
                                   sample.id=list("23_CV_SFO_005_S24_CpG.txt","23_CV_SFO_007_S25_CpG.txt", "23_CV_SFO_010_S18_CpG.txt",  "23_CV_SFO_012_S19_CpG.txt", 
                                                  "23_CV_SFO_014_S20_CpG.txt", "23_CV_SFO_015_S21_CpG.txt", "23_CV_SFO_017_S26_CpG.txt", "23_CV_SFO_024_S22_CpG.txt", 
                                                  "23_CV_SFO_025_S23_CpG.txt",  "23_RCD_SFO_001_S1_CpG.txt", "23_RCD_SFO_007_S2_CpG.txt", "23_RCD_SFO_009_S3_CpG.txt",
                                                  "23_RCD_SFO_012_S4_CpG.txt",  "23_RCD_SFO_015_S5_CpG.txt", "23_RCD_SFO_017_S6_CpG.txt",  "23_RCD_SFO_021_S7_CpG.txt", 
                                                  "23_RCD_SFO_023_S8_CpG.txt",  "23_RCD_SFO_024_S9_CpG.txt", "23_WWD_SFO_003_S10_CpG.txt", "23_WWD_SFO_004_S11_CpG.txt",
                                                  "23_WWD_SFO_008_S12_CpG.txt", "23_WWD_SFO_015_S13_CpG.txt", "23_WWD_SFO_016_S27_CpG.txt", "23_WWD_SFO_018_S14_CpG.txt",
                                                  "23_WWD_SFO_021_S15_CpG.txt", "23_WWD_SFO_023_S16_CpG.txt", "23_WWD_SFO_024_S17_CpG.txt"),
                                   assembly = "GCF_029448725.1_ASM2944872v1_genomic.fa",
                                   treatment = c(1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2),
                                   context= "CpG",
                                   mincov = 10)

#clean it up
 getMethylationStats(myobj[[4]], plot=TRUE, both.strands=FALSE)                 
 getCoverageStats(myobj[[4]], plot=TRUE, both.strands=FALSE)                 
 filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,hi.count=NULL,hi.perc=99.9)                
 normalized.myobj=normalizeCoverage(filtered.myobj) 

meth = unite(normalized.myobj, destrand=TRUE)  

pm = percMethylation(meth)
mds = matrixStats::rowSds(pm)
head(meth[mds>10,])
hist(mds, col="cornflowerblue", xlab="Std. dev. per CpG")


clusterSamples(meth, dist="correlation", method="ward.D2", plot=TRUE)

hc=clusterSamples(meth, dist="correlation", method="ward.D2", plot=FALSE)

pc=PCASamples(meth,  screeplot = TRUE, adj.lim=c(1,1))
pc=PCASamples(meth, obj.return = TRUE, adj.lim=c(1,1))

PCASamples(meth)




myDiff = calculateDiffMeth(meth)                  
                  
myDiff.hyper = getMethylDiff(myDiff,qvalue=0.01,difference=18,type="hyper")
bedgraph(myDiff.hyper, file.name = "hyper.CpG.bedGraph", col.name = "qvalue")

my.Diff.hypo = getMethylDiff(myDiff, qvalue=0.01, difference=18, type = "hypo")

# myDiff.hypo = getMethylDiff(myDiff,qvalue=0.05,difference=10,type="hypo")
bedgraph(myDiff.hypo, file.name = "hypo.CpG.bedGraph", col.name = "qvalue")                  
                  
tiles = tileMethylCounts(myobj,win.size=1000,step.size=1000,cov.bases = 10)
meth.tiles = unite(tiles, destrand=TRUE) 
                  
myDiff.tiles = calculateDiffMeth(meth.tiles)

myDiff.tiles.hyper = getMethylDiff(myDiff.tiles,qvalue=0.1,difference=10,type="hyper")
bedgraph(myDiff.tiles.hyper, file.name = "hyper.DMR.bedGraph", col.name = "qvalue")

myDiff.tiles.hypo = getMethylDiff(myDiff.tiles,qvalue=0.1,difference=10,type="hypo")
bedgraph(myDiff.tiles.hypo, file.name = "hypo.DMR.bedGraph", col.name = "qvalue")  
  
myDiff25p=getMethylDiff(myDiff,difference=15,qvalue=0.01, type = "all")

#PCR of CPG's                  
PCASamples(meth, screeplot = FALSE, adj.lim = c(1, 1), scale = TRUE, center = TRUE, comp = c(1,2), transpose = TRUE, sd.filter = TRUE, sd.threshold = 0.5, filterByQuantile = TRUE, obj.return = TRUE, chunk.size = 1e+06)                  
#pca of DMR's
PCASamples(meth.tiles, screeplot = FALSE, adj.lim = c(1, 1), scale = TRUE, center = TRUE, comp = c(1,2), transpose = TRUE, sd.filter = TRUE, sd.threshold = 0.5, filterByQuantile = TRUE, obj.return = TRUE, chunk.size = 1e+06)                

#i had to convert gtf to bed with $convert2bed -i gff < genomic.gff > ncbioutput.bed
#move bed genome to file "extdata"


gene.obj <- readTranscriptFeatures(system.file("extdata", "bedparseout.bed", package = "methylKit"), remove.unusual = FALSE)
annotateWithGeneParts(as(myDiff25p, "GRanges"), gene.obj)

cpg.obj <-  readFeatureFlank(system.file("extdata", "bedparseout.bed", package = "methylKit"), feature.flank.name = c("CpGi", "shores"))

diffCpGann = annotateWithFeatureFlank(as(myDiff25p, "GRanges"), cpg.obj$CpGi, cpg.obj$shores, feature.name = "CpGi", flank.name = "shoreS")


promoters = regionCounts(myobj, gene.obj$promoters)

head(promoters[[1]])
diffAnn = annotateWithGeneParts(as(myDiff25p, "GRanges"), gene.obj)
TSS = getAssociationWithTSS(diffAnn)
write.csv(TSS, "Sdev23_T13_12h_TSS.csv")

diffMethPerChr(myDiff,plot =FALSE,qvalue.cutoff=0.05, meth.cutoff=25)

plotTargetAnnotation(diffCpGann, col=c("green", "blue", "grey"), main = "Differential Methylation Annotation CpG")

#hyper
diffHyperAnn = annotateWithGeneParts(as(myDiff.hyper, "GRanges"), gene.obj)
TSS_Hyper = getAssociationWithTSS(diffHyperAnn)
plotTargetAnnotation(diffHyperAnn, precedence = TRUE, main = "Differential Methylation Annotation")

#hypo
diffHypoAnn = annotateWithGeneParts(as(my.Diff.hypo, "GRanges"), gene.obj)
TSS_Hyper = getAssociationWithTSS(diffHypoAnn)
plotTargetAnnotation(diffHypoAnn, precedence = TRUE, main = "Differential Methylation Annotation")


                                    