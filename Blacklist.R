library(GenomicRanges)
library(rtracklayer)
library(ChIPpeakAnno)

set.seed(123)

load('/media/usuario/Datos/eve/ifibyne/1_Genomic/2_OurProjects/1_GR-RAR_MicaSil/2_ATACseq/3_Analysis/Overlaps/maxgap10L/overlaps_all_ChIPpeakAnno.RData')
ls()
blacklist <- import.bed('/media/usuario/Datos/eve/ifibyne/1_Genomic/2_OurProjects/1_GR-RAR_MicaSil/2_ATACseq/3_Analysis/Blacklist/consensusBlacklist.bed')


findOverlapsOfPeaks(overlaps_CTRL_peaklist, blacklist)$mergedPeaks
findOverlapsOfPeaks(overlaps_DEX_peaklist, blacklist)$mergedPeaks
findOverlapsOfPeaks(overlaps_RA_peaklist, blacklist)$mergedPeaks
findOverlapsOfPeaks(overlaps_RAD_peaklist, blacklist)$mergedPeaks
