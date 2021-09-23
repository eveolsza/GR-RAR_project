path <- '/media/usuario/Datos/eve/ifibyne/1_Genomic/2_OurProjects/1_GR-RAR_MicaSil/2_ATACseq/3_Analysis/Annotation/ChIPpeakAnno/'
load(file=paste0(path,"AnnotationChIPpeakAnno_ALL.RData"))

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


Picos <- function(folder, GAINorLOSS) {
  
  name = substr(folder,1,nchar(folder)-1)
  
  path = '/media/usuario/Datos/eve/ifibyne/1_Genomic/2_OurProjects/1_GR-RAR_MicaSil/2_ATACseq/3_Analysis/Cruza_con_RNAseq/ChIPpeakAnno/'
  GAINorLOSS_AC <- read.csv(paste0(path,folder,GAINorLOSS,'.csv'))
  GAINorLOSS_AC <- makeGRangesFromDataFrame(GAINorLOSS_AC, keep.extra.columns=TRUE)
  
  # quedarme con los picos distales
  GAINorLOSS_AC_distal <- plyranges::filter(GAINorLOSS_AC, ((distancetoFeature>=3000)|(distancetoFeature<=(-5000))))

  # Filtrar las anotaciones de todo el genoma por los genes que tienen picos diferenciales entre los tratamientos a analizar en los genes AC
  AnnotationChIPpeakAnnoCTRL_filterACdistal = plyranges::filter(AnnotationChIPpeakAnnoCTRL, symbol %in% GAINorLOSS_AC_distal$symbol)
  AnnotationChIPpeakAnnoRA_filterACdistal = plyranges::filter(AnnotationChIPpeakAnnoRA, symbol %in% GAINorLOSS_AC_distal$symbol)
  AnnotationChIPpeakAnnoRAD_filterACdistal = plyranges::filter(AnnotationChIPpeakAnnoRAD, symbol %in% GAINorLOSS_AC_distal$symbol)
  
  # Quedarme con los picos en promotores
  Promoters_CTRL_filterAC = plyranges::filter(AnnotationChIPpeakAnnoCTRL_filterACdistal, (distancetoFeature > -2000 & distancetoFeature < 500))
  Promoters_RA_filterAC = plyranges::filter(AnnotationChIPpeakAnnoRA_filterACdistal, (distancetoFeature > -2000 & distancetoFeature < 500))
  Promoters_RAD_filterAC = plyranges::filter(AnnotationChIPpeakAnnoRAD_filterACdistal, (distancetoFeature > -2000 & distancetoFeature < 500))
  
  # picos distales que además tienen un pico en el promotor
  GenesCOnPicosEnPromotor <- unique(c(Promoters_CTRL_filterAC$symbol, Promoters_RA_filterAC$symbol, Promoters_RAD_filterAC$symbol))
  
  SI <- plyranges::filter(GAINorLOSS_AC_distal, symbol %in% GenesCOnPicosEnPromotor)
  write.csv(SI, paste0(path,folder,name,'_AC_distal_withPeaksInPromoter.csv'))
  
  # Creating a not in operator:
  `%notin%` <- Negate(`%in%`)
  # picos distales que NO tienen un pico en el promotor
  NO <- plyranges::filter(GAINorLOSS_AC_distal, symbol %notin% GenesCOnPicosEnPromotor)
  write.csv(NO, paste0(path,folder,name,'_AC_distal_withOUTPeaksInPromoter.csv'))

            
  #############################################
  ######### HEATMAP withPeaksInPromoter #######

  features <- SI
  df = data.frame(SI)
  first <- df[match(unique(df$symbol), df$symbol),]
  first_GRanges <- makeGRangesFromDataFrame(first, keep.extra.columns=TRUE)
  
  #### ZOOM IN

  centered = CenterPeaks(first_GRanges, 2000)
  feature.recentered = centered$Recentered
  feature.center = centered$Center
  
  
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
  
  pdf(file=paste0(path,folder,name,'_Heatmap_AC_distal_withPeaksInPromoter_ZOOMIN.pdf'))

  heatmap <- featureAlignedHeatmap(sig, feature.center, 
                                   upstream=2000, downstream=2000,
                                   upper.extreme=15, 
                                   color=colorRampPalette(c("black", "white", "red"))(50))
  dev.off()
  
  
  #### ZOOM OUT
  
  centered = CenterPeaks(first_GRanges, 1000000)
  feature.recentered = centered$Recentered
  feature.center = centered$Center
  
  
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
                              upstream=1000000, downstream=1000000)
  
  pdf(file=paste0(path,folder,name,'_Heatmap_AC_distal_withPeaksInPromoter_ZOOMOUT.pdf'))
  
  heatmap <- featureAlignedHeatmap(sig, feature.center, 
                                   upstream=1000000, downstream=1000000,
                                   upper.extreme=15, 
                                   color=colorRampPalette(c("black", "yellow",  "blue", "orange", "red"))(50))
  dev.off()
  
  
  #############################################
  ######### HEATMAP withOUTPeaksInPromoter ####
  
  features <- NO
  df = data.frame(NO)
  first <- df[match(unique(df$symbol), df$symbol),]
  first_GRanges <- makeGRangesFromDataFrame(first, keep.extra.columns=TRUE)
  
  #### ZOOM IN
  
  centered = CenterPeaks(first_GRanges, 2000)
  feature.recentered = centered$Recentered
  feature.center = centered$Center
  
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
  
  pdf(file=paste0(path,folder,name,'_Heatmap_AC_distal_withOUTPeaksInPromoter_ZOOMIN.pdf'))
  
  heatmap <- featureAlignedHeatmap(sig, feature.center, 
                                   upstream=2000, downstream=2000,
                                   upper.extreme=15, 
                                   color=colorRampPalette(c("black", "white", "red"))(50))
  dev.off()
  
  #### ZOOM OUT
  
  centered = CenterPeaks(first_GRanges, 1000000)
  feature.recentered = centered$Recentered
  feature.center = centered$Center
  
  
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
                              upstream=1000000, downstream=1000000)
  
  pdf(file=paste0(path,folder,name,'_Heatmap_AC_distal_withOUTPeaksInPromoter_ZOOMOUT.pdf'))
  
  heatmap <- featureAlignedHeatmap(sig, feature.center, 
                                   upstream=1000000, downstream=1000000,
                                   upper.extreme=15, 
                                   color=colorRampPalette(c("black", "yellow", "blue", "orange", "red"))(50))
  dev.off()
  
  
            
  
}
  
Picos('RADandRAvsCTRL_gain/','RADandRAvsCTRL_gain_AnnotationChIPpeakAnno_FilterAC')
Picos('RADandRAvsCTRL_loss/','RADandRAvsCTRL_loss_AnnotationChIPpeakAnno_FilterAC')
Picos('RADvsRA_gain/','RADvsRA_gain_AnnotationChIPpeakAnno_FilterAC')
Picos('RADvsRA_loss/','RADvsRA_loss_AnnotationChIPpeakAnno_FilterAC')
