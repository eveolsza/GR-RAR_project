library(ChIPpeakAnno)
library(plyranges)
library(VennDiagram)
library(rtracklayer)

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
annoData <- toGRanges(TxDb.Hsapiens.UCSC.hg19.knownGene, feature="gene")

library(org.Hs.eg.db)

source('functions.R')

set.seed(123)

path = '/media/usuario/Datos/eve/ifibyne/1_Genomic/2_OurProjects/1_GR-RAR_MicaSil/2_ATACseq/'
load(paste0(path,'1_Data_Transfromed/Diffpeaks_Matt_ChangedName/Diffpeaks_Matt_ChangedName.RData'))
load(paste0(path,'3_Analysis/Annotation/ChIPpeakAnno/AnnotationChIPpeakAnno_ALL.RData'))
ls()

A_C = read.csv('/media/usuario/Datos/eve/ifibyne/1_Genomic/2_OurProjects/1_GR-RAR_MicaSil/1_bulkRNAseq/tablas/A+C.csv', sep=',')
head(A_C)
genes <- A_C$gene_name


BuscarPicosEnPromotores=function(archivo){
  PicosEnPromotores <- GRanges()
  for (fila in seq_along(archivo)){
    x = archivo[fila]
    distStart = start(x) - x$start_position
    distEnd = end(x) - x$start_position
    if(x$feature_strand == "-") {
      distStart =  x$end_position - start(x)
      distEnd = x$end_position - end(x)
    }
    if ((-2000<=distStart & distStart<=500)|(-2000<=distEnd & distEnd<=500)){
      PicosEnPromotores <- append(PicosEnPromotores, x)
    }
  }
  return(PicosEnPromotores)
}



GenesConPicosDIFFdistalesyConPicosEnPromotor = function(GAINorLOSS){
  
  annotatedPeaks <- annotatePeakInBatch(GAINorLOSS, 
                                        AnnotationData = annoData, 
                                        PeakLocForDistance='middle',
                                        FeatureLocForDistance ='TSS',
                                        output="both", maxgap=1000000)
  
  annotatedPeaks_IDs <- addGeneIDs(annotatedPeaks, orgAnn="org.Hs.eg.db", IDs2Add=c("symbol"),
                                   feature_id_type = "entrez_id")
  
  # quearme con los picos asociados a los genes de interes
  annotatedPeaks_IDs_Filter <- plyranges::filter(annotatedPeaks_IDs, symbol %in% genes)
  annotatedPeaks_IDs_Filter$closestGene <- NULL
  annotatedPeaks_IDs_Filter_df <- data.frame(annotatedPeaks_IDs_Filter)
  
  # quedarme con los picos distales
  GAINorLOSS_AC_distal <- plyranges::filter(annotatedPeaks_IDs_Filter, ((distancetoFeature>=3000)|(distancetoFeature<=(-5000))))
  
  # Filtrar las anotaciones de todo el genoma por los genes que tienen picos diferenciales entre los tratamientos a analizar en los genes AC
  AnnotationChIPpeakAnnoCTRL_filterACdistal = plyranges::filter(AnnotationChIPpeakAnnoCTRL, symbol %in% GAINorLOSS_AC_distal$symbol)
  AnnotationChIPpeakAnnoRA_filterACdistal = plyranges::filter(AnnotationChIPpeakAnnoRA, symbol %in% GAINorLOSS_AC_distal$symbol)
  AnnotationChIPpeakAnnoRAD_filterACdistal = plyranges::filter(AnnotationChIPpeakAnnoRAD, symbol %in% GAINorLOSS_AC_distal$symbol)
  
  Promoters_CTRL_filterAC = BuscarPicosEnPromotores(AnnotationChIPpeakAnnoCTRL_filterACdistal)
  Promoters_RA_filterAC = BuscarPicosEnPromotores(AnnotationChIPpeakAnnoRA_filterACdistal)
  Promoters_RAD_filterAC = BuscarPicosEnPromotores(AnnotationChIPpeakAnnoRAD_filterACdistal)
  
  # picos distales que además tienen un pico en el promotor
  GenesCOnPicosEnPromotor <- unique(c(Promoters_CTRL_filterAC$symbol, Promoters_RA_filterAC$symbol, Promoters_RAD_filterAC$symbol))
  
  SI <- plyranges::filter(GAINorLOSS_AC_distal, symbol %in% GenesCOnPicosEnPromotor)
  
  return(SI)
  
}

RADandRAvsCTRL_gain <- c(RADvsCTRL_gain, RAvsCTRL_gain)
RADandRAvsCTRL_gain_unique <- unique(RADandRAvsCTRL_gain)
RADandRAvsCTRL_gain_reduced <- GenomicRanges::reduce(RADandRAvsCTRL_gain)
RADandRAvsCTRL_gain_picosPromotores = GenesConPicosDIFFdistalesyConPicosEnPromotor(RADandRAvsCTRL_gain_reduced)

RADandRAvsCTRL_loss <- c(RADvsCTRL_loss, RAvsCTRL_loss)
RADandRAvsCTRL_loss_unique <- unique(RADandRAvsCTRL_loss)
RADandRAvsCTRL_loss_reduced <- GenomicRanges::reduce(RADandRAvsCTRL_loss)
RADandRAvsCTRL_loss_picosPromotores = GenesConPicosDIFFdistalesyConPicosEnPromotor(RADandRAvsCTRL_loss_reduced)

RADandRAvsCTRL_gainYloss_picosPromotores = unique(c(RADandRAvsCTRL_gain_picosPromotores$symbol, RADandRAvsCTRL_loss_picosPromotores$symbol))
write.csv(RADandRAvsCTRL_gainYloss_picosPromotores, "RADandRAvsCTRL_gainYloss_picosPromotores.csv")


genes_RADandRAvsCTRL_gainYloss_AC_distal_withPeaksInPromoter_centroYcorta = read.csv("/media/usuario/Datos/eve/ifibyne/1_Genomic/2_OurProjects/1_GR-RAR_MicaSil/2_ATACseq/3_Analysis/Cruza_con_RNAseq/ChIPpeakAnno/A+C/genes_RADandRAvsCTRL_gainYloss_AC_distal_withPeaksInPromoter_centroYcorta.csv")

setdiff(genes_RADandRAvsCTRL_gainYloss_AC_distal_withPeaksInPromoter_centroYcorta$x, RADandRAvsCTRL_gainYloss_picosPromotores)

setdiff(RADandRAvsCTRL_gainYloss_picosPromotores, genes_RADandRAvsCTRL_gainYloss_AC_distal_withPeaksInPromoter_centroYcorta$x)

c("STK24", "SLC12A4", "PLD1", "B3GNT5", "SH3RF3-AS1", "TOB1", "PCNX1", "LOC100499489", "NETO2", "ETS2") %in% RADandRAvsCTRL_gainYloss_picosPromotores
