library(ChIPpeakAnno)

path = '/media/usuario/Datos/eve/ifibyne/1_Genomic/2_OurProjects/1_GR-RAR_MicaSil/2_ATACseq/'
load(paste0(path,'3_Analysis/Overlaps/maxgap10L/overlaps_all_ChIPpeakAnno.RData'))

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
annoData <- toGRanges(TxDb.Hsapiens.UCSC.hg19.knownGene, feature="gene")

library(org.Hs.eg.db)


set.seed(123)

AnnotationChIPpeakAnno <- function(file){
  
  annotatedPeaks <- annotatePeakInBatch(file, 
                                        AnnotationData = annoData, 
                                        PeakLocForDistance='middle',
                                        FeatureLocForDistance ='TSS',
                                        output="both", maxgap=1000000)
  
  annotatedPeaks_IDs <- addGeneIDs(annotatedPeaks, orgAnn="org.Hs.eg.db", IDs2Add=c("symbol"),
                                   feature_id_type = "entrez_id")
  annotatedPeaks_IDs$closestGene <- NULL
  print(annotatedPeaks_IDs)
  annotatedPeaks_IDs_df <- data.frame(annotatedPeaks_IDs)
  
  path <- '/media/usuario/Datos/eve/ifibyne/1_Genomic/2_OurProjects/1_GR-RAR_MicaSil/2_ATACseq/3_Analysis/Annotation/ChIPpeakAnno/'
  write.csv(annotatedPeaks_IDs_df, paste0(path,deparse(substitute(file)),'_AnnotationChIPpeakAnno.csv'))
  
  return(annotatedPeaks_IDs)
}

AnnotationChIPpeakAnnoCTRL <- AnnotationChIPpeakAnno(overlaps_CTRL_peaklist)
AnnotationChIPpeakAnnoDEX <- AnnotationChIPpeakAnno(overlaps_DEX_peaklist)
AnnotationChIPpeakAnnoRA <- AnnotationChIPpeakAnno(overlaps_RA_peaklist)
AnnotationChIPpeakAnnoRAD <- AnnotationChIPpeakAnno(overlaps_RAD_peaklist)

path <- '/media/usuario/Datos/eve/ifibyne/1_Genomic/2_OurProjects/1_GR-RAR_MicaSil/2_ATACseq/3_Analysis/Annotation/ChIPpeakAnno/'
save(AnnotationChIPpeakAnnoCTRL, 
     AnnotationChIPpeakAnnoDEX,
     AnnotationChIPpeakAnnoRA, 
     AnnotationChIPpeakAnnoRAD,
     file=paste0(path,"AnnotationChIPpeakAnno_ALL.RData"))
     
