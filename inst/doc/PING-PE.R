### R code from vignette source 'PING-PE.Rnw'

###################################################
### code chunk number 1: PING-PE.Rnw:26-27
###################################################
options(continue=" ")


###################################################
### code chunk number 2: Loading-PING
###################################################
library(PING)


###################################################
### code chunk number 3: Read-data
###################################################
yeastBam<- system.file("extdata/yeastChrI.bam",package="PING")


###################################################
### code chunk number 4: bam2gr
###################################################
library(PICS)
gr<-bam2gr(bamFile=yeastBam, PE=TRUE)


###################################################
### code chunk number 5: subset-GR
###################################################
grI<-gr[seqnames(gr)=="chrI"]
seqlevels(grI)<-"chrI"


###################################################
### code chunk number 6: Genome-segmentation
###################################################
segPE<-segmentPING(grI, PE=TRUE)


###################################################
### code chunk number 7: Cluster-initialization
###################################################
library(parallel)


###################################################
### code chunk number 8: PING-analysis
###################################################
ping<-PING(segPE, nCores=2)


###################################################
### code chunk number 9: Post-process-PING-result
###################################################
PS=postPING(ping, segPE)


###################################################
### code chunk number 10: makeRangedDataOutput (eval = FALSE)
###################################################
## rdBed<-makeRangedDataOutput(PS, type="bed")
## library(rtracklayer)
## export(rdBed, "nucPrediction.bed")


###################################################
### code chunk number 11: plotSummary-PE
###################################################
plotSummary(PS, ping,  grI, chr="chrI", from=149000, to=153000)


