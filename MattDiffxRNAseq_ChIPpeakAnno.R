library(ChIPpeakAnno)
library(plyranges)
library(VennDiagram)
library(rtracklayer)

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
annoData <- toGRanges(TxDb.Hsapiens.UCSC.hg19.knownGene, feature="gene")

library(org.Hs.eg.db)



set.seed(123)

path = '/media/usuario/Datos/eve/ifibyne/1_Genomic/2_OurProjects/1_GR-RAR_MicaSil/2_ATACseq/'
#load(paste0(path,'3_Analysis/Overlaps/maxgap10L/overlaps_all_ChIPpeakAnno.RData'))
load(paste0(path,'1_Data_Transfromed/Diffpeaks_Matt_ChangedName/Diffpeaks_Matt_ChangedName.RData'))
ls()

A_C = read.csv('/media/usuario/Datos/eve/ifibyne/1_Genomic/2_OurProjects/1_GR-RAR_MicaSil/1_bulkRNAseq/tablas/A+C.csv', sep=',')
head(A_C)

CenterPeaks <- function(file, distancia){
  
  file = format_granges2(file)
  feature.recentered = file
  feature.centered = file
  
  start(feature.centered) <- 1
  end(feature.centered) <- 2
  
  start(feature.recentered) <- 1
  end(feature.recentered) <- 2
  
  for (i in seq(length(file))){
    if (file[i]$feature_strand == '+'){
      end(feature.centered)[i] = file[i]$start_position
      start(feature.centered)[i] = file[i]$start_position
    }
    if (file[i]$feature_strand == '-'){
      end(feature.centered)[i] = file[i]$end_position
      start(feature.centered)[i] = file[i]$end_position
    }
  }
  end(feature.recentered) <- end(feature.centered) + distancia
  start(feature.recentered) <- start(feature.centered) - distancia
  
  end = list(Center=feature.centered, Recentered=feature.recentered)
  
  return(end)
  
}



Annotate_Filter_VennDiagram_HeatMap <- function(GAINorLOSS, genes){
  
  annotatedPeaks <- annotatePeakInBatch(GAINorLOSS, 
                                        AnnotationData = annoData, 
                                        PeakLocForDistance='middle',
                                        FeatureLocForDistance ='TSS',
                                        output="both", maxgap=1000000)
  
  annotatedPeaks_IDs <- addGeneIDs(annotatedPeaks, orgAnn="org.Hs.eg.db", IDs2Add=c("symbol"),
                                   feature_id_type = "entrez_id")
  
  annotatedPeaks_IDs_Filter <- plyranges::filter(annotatedPeaks_IDs, symbol %in% genes)
  annotatedPeaks_IDs_Filter$closestGene <- NULL
  annotatedPeaks_IDs_Filter_df <- data.frame(annotatedPeaks_IDs_Filter)
  write.csv(annotatedPeaks_IDs_Filter_df, paste0(path,deparse(substitute(GAINorLOSS)),'_AnnotationChIPpeakAnno_FilterAC.csv'))
  
  ##################################################
  # DIAGRAMA DE VENN CON PICOS DISTALES Y PROXIMALES
  
  # calculo picos distales y los guardo
  distal <- plyranges::filter(annotatedPeaks_IDs_Filter, ((distancetoFeature>=3000)|(distancetoFeature<=(-5000))))
  # calculo picos proximales y los guardo  
  proximal <- plyranges::filter(annotatedPeaks_IDs_Filter, ((distancetoFeature<3000)&(distancetoFeature>(-5000))))
  
  # busco los genes que están en cada parte del diagrama de Venn  
  distal_only = setdiff(unique(distal$symbol), unique(proximal$symbol))
  intersect = intersect(unique(distal$symbol), unique(proximal$symbol))
  proximal_only = setdiff(unique(proximal$symbol), unique(distal$symbol))
  
  # los guardo en una tabla   
  max.len = max(length(distal_only), length(intersect), length(proximal_only))
  
  distal_only = c(distal_only, rep(NA, max.len - length(distal_only)))
  intersect = c(intersect, rep(NA, max.len - length(intersect)))
  proximal_only = c(proximal_only, rep(NA, max.len - length(proximal_only)))
  
  write.csv(data.frame(distal_only,
                       intersect,
                       proximal_only), 
            paste0(path,deparse(substitute(GAINorLOSS)),'_AC_VennDiagram.csv'))
  
  
  # grafico el diagrama de Venn propiamente dicho 
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  venn.diagram(list(Distal = unique(distal$symbol), Proximal = unique(proximal$symbol)), 
               fill = c("white", "white"), col = c('black', 'black'), lwd =1, main=deparse(substitute(file)),
               paste0(path,deparse(substitute(GAINorLOSS)),"_AC_VennDiagram.pdf"))
  
  
  #####################################################################
  ############################## HEATMAP ##############################
  
  #### Center Peaks
  
  features <- annotatedPeaks_IDs_Filter
  feature.recentered <- reCenterPeaks(features, width=2000)
  feature.center <- reCenterPeaks(features, width=1)
  
  # cargo los bigwigs centrados en las features a analizar
  path_bigwigs='/media/usuario/Datos/eve/ifibyne/1_Genomic/2_OurProjects/1_GR-RAR_MicaSil/2_ATACseq/1_Data_Transfromed/MeanReplicatesBigWigs/'
  
  CTRL_bw <- rtracklayer::import(paste0(path_bigwigs,'CTRL1_CTRL2_mean.bw'),
                                 format="BigWig",
                                 which=feature.recentered,
                                 as="RleList") 
  
  
  DEX_bw <- rtracklayer::import(paste0(path_bigwigs,'DEX1_DEX2_mean.bw'),
                                format="BigWig",
                                which=feature.recentered,
                                as="RleList") 
  
  RA_bw <- rtracklayer::import(paste0(path_bigwigs,'RA1_RA2_mean.bw'),
                               format="BigWig",
                               which=feature.recentered,
                               as="RleList") 
  
  RAD_bw <- rtracklayer::import(paste0(path_bigwigs,'RAD1_RAD2_mean.bw'),
                                format="BigWig",
                                which=feature.recentered,
                                as="RleList") 
  
  cvglists <- list(CTRL=CTRL_bw,
                   DEX=DEX_bw,
                   RA=RA_bw, 
                   RA_DEX=RAD_bw)
  
  sig <- featureAlignedSignal(cvglists, feature.center, 
                              upstream=2000, downstream=2000)
  
  # graficar el Heatmap
  pdf(file=paste0(path,paste0(deparse(substitute(GAINorLOSS)),'_Heatmap_centerPeaks.pdf')))
  heatmap <- featureAlignedHeatmap(sig, feature.center, 
                                   upstream=2000, downstream=2000,
                                   upper.extreme=15,
                                   color=colorRampPalette(c("blue", "white", "red"))(50))
  dev.off()                       
  
  
  #### Center TSS

  df = data.frame(annotatedPeaks_IDs_Filter)
  first <- df[match(unique(df$symbol), df$symbol),]
  first_GRanges <- makeGRangesFromDataFrame(first, keep.extra.columns=TRUE)
  
  centered = CenterPeaks(first_GRanges, 2000)
  feature.recentered = centered$Recentered
  feature.center = centered$Center

  # cargo los bigwigs centrados en las features a analizar
  path_bigwigs='/media/usuario/Datos/eve/ifibyne/1_Genomic/2_OurProjects/1_GR-RAR_MicaSil/2_ATACseq/1_Data_Transfromed/MeanReplicatesBigWigs/'
  
  CTRL_bw <- rtracklayer::import(paste0(path_bigwigs,'CTRL1_CTRL2_mean.bw'),
                                 format="BigWig",
                                 which=feature.recentered,
                                 as="RleList") 
  
  
  DEX_bw <- rtracklayer::import(paste0(path_bigwigs,'DEX1_DEX2_mean.bw'),
                                format="BigWig",
                                which=feature.recentered,
                                as="RleList") 
  
  RA_bw <- rtracklayer::import(paste0(path_bigwigs,'RA1_RA2_mean.bw'),
                               format="BigWig",
                               which=feature.recentered,
                               as="RleList") 
  
  RAD_bw <- rtracklayer::import(paste0(path_bigwigs,'RAD1_RAD2_mean.bw'),
                                format="BigWig",
                                which=feature.recentered,
                                as="RleList") 
  
  cvglists <- list(CTRL=CTRL_bw,
                   DEX=DEX_bw,
                   RA=RA_bw, 
                   RA_DEX=RAD_bw)
  
  sig <- featureAlignedSignal(cvglists, feature.center, 
                              upstream=2000, downstream=2000)
  
  # graficar el Heatmap
  pdf(file=paste0(path,paste0(deparse(substitute(GAINorLOSS)),'_Heatmap_centerTSS.pdf')))
  heatmap <- featureAlignedHeatmap(sig, feature.center, 
                                   upstream=2000, downstream=2000,
                                   upper.extreme=15,
                                   color=colorRampPalette(c("blue", "white", "red"))(50))
  dev.off()                       
  
  
}                    


RADandRAvsCTRL_gain <- c(RADvsCTRL_gain, RAvsCTRL_gain)
RADandRAvsCTRL_gain <- unique(RADandRAvsCTRL_gain)

Annotate_Filter_VennDiagram_HeatMap(RADandRAvsCTRL_gain, A_C$gene_name)

RADandRAvsCTRL_loss <- c(RADvsCTRL_loss, RAvsCTRL_loss)
RADandRAvsCTRL_loss <- unique(RADandRAvsCTRL_loss)

Annotate_Filter_VennDiagram_HeatMap(RADandRAvsCTRL_loss, A_C$gene_name)

Annotate_Filter_VennDiagram_HeatMap(RADvsRA_gain, A_C$gene_name)

Annotate_Filter_VennDiagram_HeatMap(RADvsRA_loss, A_C$gene_name)

