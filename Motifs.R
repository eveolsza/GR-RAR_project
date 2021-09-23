library(memes)
library(magrittr)
library(ggplot2)
library(IRanges)
library(GenomicRanges)
library(ChIPpeakAnno)
library(tidyverse)
library(universalmotif)
library(Biostrings)

suppressPackageStartupMessages(library(GenomicRanges))

options(meme_bin = "/home/usuario/Programs/meme-5.4.1/src/")
check_meme_install()

HS.genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

meme_HS <- read_meme("/media/usuario/Datos/eve/ifibyne/1_Genomic/2_OurProjects/1_GR-RAR_MicaSil/2_ATACseq/Motivos_Mica/motifs/motif_databases.12.18/motif_databases/HUMAN/HOCOMOCOv11_full_HUMAN_mono_meme_format.meme") %>% 
  to_df()
head(meme_HS)

options(meme_HS = to_list(meme_HS, extrainfo = FALSE))

#############################################
#              PICOS DISTALES              #
#############################################

RADandRAvsCTRL_loss_AC_distal_withPeaksInPromoter = read.csv('/media/usuario/Datos/eve/ifibyne/1_Genomic/2_OurProjects/1_GR-RAR_MicaSil/2_ATACseq/3_Analysis/Cruza_con_RNAseq/ChIPpeakAnno/A+C/RADandRAvsCTRL_loss/RADandRAvsCTRL_loss_AC_distal_withPeaksInPromoter.csv')
RADandRAvsCTRL_loss_AC_distal_withPeaksInPromoter = makeGRangesFromDataFrame(RADandRAvsCTRL_loss_AC_distal_withPeaksInPromoter, keep.extra.columns=TRUE)
length(unique(RADandRAvsCTRL_loss_AC_distal_withPeaksInPromoter$symbol))
newStyle <- mapSeqlevels(seqlevels(RADandRAvsCTRL_loss_AC_distal_withPeaksInPromoter), "UCSC")
RADandRAvsCTRL_loss_AC_distal_withPeaksInPromoter <- renameSeqlevels(RADandRAvsCTRL_loss_AC_distal_withPeaksInPromoter, newStyle)

RADandRAvsCTRL_gain_AC_distal_withPeaksInPromoter = read.csv('/media/usuario/Datos/eve/ifibyne/1_Genomic/2_OurProjects/1_GR-RAR_MicaSil/2_ATACseq/3_Analysis/Cruza_con_RNAseq/ChIPpeakAnno/A+C/RADandRAvsCTRL_gain/RADandRAvsCTRL_gain_AC_distal_withPeaksInPromoter.csv')
RADandRAvsCTRL_gain_AC_distal_withPeaksInPromoter = makeGRangesFromDataFrame(RADandRAvsCTRL_gain_AC_distal_withPeaksInPromoter, keep.extra.columns=TRUE)
length(unique(RADandRAvsCTRL_gain_AC_distal_withPeaksInPromoter$symbol))
newStyle <- mapSeqlevels(seqlevels(RADandRAvsCTRL_gain_AC_distal_withPeaksInPromoter), "UCSC")
RADandRAvsCTRL_gain_AC_distal_withPeaksInPromoter <- renameSeqlevels(RADandRAvsCTRL_gain_AC_distal_withPeaksInPromoter, newStyle)

RADandRAvsCTRL_gainYloss_AC_distal_withPeaksInPromoter <- c(RADandRAvsCTRL_gain_AC_distal_withPeaksInPromoter, RADandRAvsCTRL_loss_AC_distal_withPeaksInPromoter)
RADandRAvsCTRL_gainYloss_AC_distal_withPeaksInPromoter_reduced <- GenomicRanges::reduce(RADandRAvsCTRL_gainYloss_AC_distal_withPeaksInPromoter)

RADandRAvsCTRL_gainYloss_AC_distal_withPeaksInPromoter_reduced_sequence <- RADandRAvsCTRL_gainYloss_AC_distal_withPeaksInPromoter_reduced %>%
  # Get a list of chip peaks 
  # look up the DNA sequence of each peak within each group
  get_sequence(HS.genome)
writeXStringSet(unique(RADandRAvsCTRL_gainYloss_AC_distal_withPeaksInPromoter_reduced_sequence), filepath = paste0("/media/usuario/Datos/eve/ifibyne/1_Genomic/2_OurProjects/1_GR-RAR_MicaSil/2_ATACseq/3_Analysis/Cruza_con_RNAseq/ChIPpeakAnno/A+C/RADandRAvsCTRL_gainYloss/RADandRAvsCTRL_gainYloss_AC_distal_withPeaksInPromoter_sequence.fasta"), format = "fasta")


#### CONTROL for STREME
# distal peaks that do not change between RADandRAvsCTRL
load('/media/usuario/Datos/eve/ifibyne/1_Genomic/2_OurProjects/1_GR-RAR_MicaSil/2_ATACseq/3_Analysis/Annotation/ChIPpeakAnno/AnnotationChIPpeakAnno_ALL.RData')
rm(AnnotationChIPpeakAnnoDEX)

AnnotationChIPpeakAnnoCTRL_distal <- plyranges::filter(AnnotationChIPpeakAnnoCTRL, ((distancetoFeature>=3000)|(distancetoFeature<=(-5000))))
AnnotationChIPpeakAnnoRA_distal <- plyranges::filter(AnnotationChIPpeakAnnoRA, ((distancetoFeature>=3000)|(distancetoFeature<=(-5000))))
AnnotationChIPpeakAnnoRAD_distal <- plyranges::filter(AnnotationChIPpeakAnnoRAD, ((distancetoFeature>=3000)|(distancetoFeature<=(-5000))))

AnnotationChIPpeakAnnoCTRL_distal_filterAC = plyranges::filter(AnnotationChIPpeakAnnoCTRL_distal, symbol %in% RADandRAvsCTRL_gainYloss_AC_distal_withPeaksInPromoter$symbol)
AnnotationChIPpeakAnnoRA_distal_filterAC = plyranges::filter(AnnotationChIPpeakAnnoRA_distal, symbol %in% RADandRAvsCTRL_gainYloss_AC_distal_withPeaksInPromoter$symbol)
AnnotationChIPpeakAnnoRAD_distal_filterAC = plyranges::filter(AnnotationChIPpeakAnnoRAD_distal, symbol %in% RADandRAvsCTRL_gainYloss_AC_distal_withPeaksInPromoter$symbol)

AnnotationChIPpeakAnnoCTRLRARAD_distal_filterAC <- c(AnnotationChIPpeakAnnoCTRL_distal_filterAC, AnnotationChIPpeakAnnoRA_distal_filterAC, AnnotationChIPpeakAnnoRAD_distal_filterAC)

diff = setdiff(AnnotationChIPpeakAnnoCTRLRARAD_distal_filterAC, RADandRAvsCTRL_gainYloss_AC_distal_withPeaksInPromoter)

#chequear que no haya solapamientos
findOverlaps(RADandRAvsCTRL_gainYloss_AC_distal_withPeaksInPromoter, diff)

diff_sequence <- diff %>%
  # Get a list of chip peaks 
  # look up the DNA sequence of each peak within each group
  get_sequence(HS.genome)

writeXStringSet(unique(diff_sequence), filepath = paste0("/media/usuario/Datos/eve/ifibyne/1_Genomic/2_OurProjects/1_GR-RAR_MicaSil/2_ATACseq/3_Analysis/Cruza_con_RNAseq/ChIPpeakAnno/A+C/RADandRAvsCTRL_gainYloss/Motifs/CONTROL.fasta"), format = "fasta")

# Run STREME -> no lo pude hacer por un error, lo corro en la web
# streme = runStreme(unique(by_binding), control = unique(diff))




#############################################
#           PICOS EN EL PROMOTOR            #
#############################################

# Quedarme con los picos en promotores
Promoters_CTRL = plyranges::filter(AnnotationChIPpeakAnnoCTRL, (distancetoFeature > -2000 & distancetoFeature < 500))
Promoters_RA = plyranges::filter(AnnotationChIPpeakAnnoRA, (distancetoFeature > -2000 & distancetoFeature < 500))
Promoters_RAD = plyranges::filter(AnnotationChIPpeakAnnoRAD, (distancetoFeature > -2000 & distancetoFeature < 500))

# Filtrar las anotaciones de todo el genoma por los genes que tienen picos diferenciales entre los tratamientos a analizar en los genes AC
Promoters_CTRL_filterAC = plyranges::filter(Promoters_CTRL, symbol %in% RADandRAvsCTRL_gainYloss_AC_distal_withPeaksInPromoter$symbol)
Promoters_RA_filterAC = plyranges::filter(Promoters_RA, symbol %in% RADandRAvsCTRL_gainYloss_AC_distal_withPeaksInPromoter$symbol)
Promoters_RAD_filterAC = plyranges::filter(Promoters_RAD, symbol %in% RADandRAvsCTRL_gainYloss_AC_distal_withPeaksInPromoter$symbol)

Promotores <- c(Promoters_CTRL_filterAC, Promoters_RA_filterAC, Promoters_RAD_filterAC)
Promotores_reduced <- GenomicRanges::reduce(Promotores)  

Promotores_reduced_sequence <- Promotores_reduced %>%
  # Get a list of chip peaks 
  # look up the DNA sequence of each peak within each group
  get_sequence(HS.genome)
writeXStringSet(Promotores_reduced_sequence, filepath = paste0("/media/usuario/Datos/eve/ifibyne/1_Genomic/2_OurProjects/1_GR-RAR_MicaSil/2_ATACseq/3_Analysis/Cruza_con_RNAseq/ChIPpeakAnno/A+C/RADandRAvsCTRL_gainYloss/Motifs/AME/Promotores/Promotores_reduced_sequence.fasta"), format = "fasta")


`%notin%` <- Negate(`%in%`)

# CONTROL

genesNOregulados <-read.csv('/media/usuario/Datos/eve/ifibyne/1_Genomic/2_OurProjects/1_GR-RAR_MicaSil/1_bulkRNAseq/genesNOregulados.csv', header = T)
length(genesNOregulados$name)
length(unique(genesNOregulados$name))

Promoters_CTRL_CONTROL = plyranges::filter(Promoters_CTRL, symbol %in% genesNOregulados$name)
Promoters_RA_CONTROL = plyranges::filter(Promoters_RA, symbol %in% genesNOregulados$name)
Promoters_RAD_CONTROL = plyranges::filter(Promoters_RAD, symbol %in% genesNOregulados$name)

Promotores_CONTROL <- c(Promoters_CTRL_CONTROL, Promoters_RA_CONTROL, Promoters_RAD_CONTROL)
length(Promotores_CONTROL$symbol)
length(unique(Promotores_CONTROL$symbol))
Promotores_CONTROL_reduced <- GenomicRanges::reduce(Promotores_CONTROL)  
length(Promotores_CONTROL_reduced)

setdiff(Promotores_CONTROL$symbol, genesNOregulados$name)
setdiff(genesNOregulados$name, Promotores_CONTROL$symbol)

Promotores_CONTROL_reduced_sequence <- Promotores_CONTROL_reduced %>%
  # Get a list of chip peaks 
  # look up the DNA sequence of each peak within each group
  get_sequence(HS.genome)
writeXStringSet(Promotores_CONTROL_reduced_sequence, filepath = paste0("/media/usuario/Datos/eve/ifibyne/1_Genomic/2_OurProjects/1_GR-RAR_MicaSil/2_ATACseq/3_Analysis/Cruza_con_RNAseq/ChIPpeakAnno/A+C/RADandRAvsCTRL_gainYloss/Motifs/AME/Promotores/Promotores_CONTROL_reduced_sequence.fasta"), format = "fasta")


