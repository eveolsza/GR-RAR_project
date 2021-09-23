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

RADandRAvsCTRL_loss_AC_distal_withPeaksInPromoter_sequence <- RADandRAvsCTRL_loss_AC_distal_withPeaksInPromoter %>%
  # Get a list of chip peaks 
  # look up the DNA sequence of each peak within each group
  get_sequence(HS.genome)
writeXStringSet(unique(RADandRAvsCTRL_loss_AC_distal_withPeaksInPromoter_sequence), filepath = paste0("/media/usuario/Datos/eve/ifibyne/1_Genomic/2_OurProjects/1_GR-RAR_MicaSil/2_ATACseq/3_Analysis/Cruza_con_RNAseq/ChIPpeakAnno/A+C/RADandRAvsCTRL_loss/Motifs/RADandRAvsCTRL_loss_AC_distal_withPeaksInPromoter_sequence.fasta"), format = "fasta")


#### CONTROL for STREME
# distal peaks that do not change between RADandRAvsCTRL
load('/media/usuario/Datos/eve/ifibyne/1_Genomic/2_OurProjects/1_GR-RAR_MicaSil/2_ATACseq/3_Analysis/Annotation/ChIPpeakAnno/AnnotationChIPpeakAnno_ALL.RData')
rm(AnnotationChIPpeakAnnoDEX)

AnnotationChIPpeakAnnoCTRL_distal <- plyranges::filter(AnnotationChIPpeakAnnoCTRL, ((distancetoFeature>=3000)|(distancetoFeature<=(-5000))))
AnnotationChIPpeakAnnoRA_distal <- plyranges::filter(AnnotationChIPpeakAnnoRA, ((distancetoFeature>=3000)|(distancetoFeature<=(-5000))))
AnnotationChIPpeakAnnoRAD_distal <- plyranges::filter(AnnotationChIPpeakAnnoRAD, ((distancetoFeature>=3000)|(distancetoFeature<=(-5000))))

AnnotationChIPpeakAnnoCTRL_distal_filterAC = plyranges::filter(AnnotationChIPpeakAnnoCTRL_distal, symbol %in% RADandRAvsCTRL_loss_AC_distal_withPeaksInPromoter$symbol)
AnnotationChIPpeakAnnoRA_distal_filterAC = plyranges::filter(AnnotationChIPpeakAnnoRA_distal, symbol %in% RADandRAvsCTRL_loss_AC_distal_withPeaksInPromoter$symbol)
AnnotationChIPpeakAnnoRAD_distal_filterAC = plyranges::filter(AnnotationChIPpeakAnnoRAD_distal, symbol %in% RADandRAvsCTRL_loss_AC_distal_withPeaksInPromoter$symbol)

AnnotationChIPpeakAnnoCTRLRARAD_distal_filterAC <- c(AnnotationChIPpeakAnnoCTRL_distal_filterAC, AnnotationChIPpeakAnnoRA_distal_filterAC, AnnotationChIPpeakAnnoRAD_distal_filterAC)

diff = setdiff(AnnotationChIPpeakAnnoCTRLRARAD_distal_filterAC, RADandRAvsCTRL_loss_AC_distal_withPeaksInPromoter)

RADandRAvsCTRL_gain_AC_distal_withPeaksInPromoter = read.csv('/media/usuario/Datos/eve/ifibyne/1_Genomic/2_OurProjects/1_GR-RAR_MicaSil/2_ATACseq/3_Analysis/Cruza_con_RNAseq/ChIPpeakAnno/A+C/RADandRAvsCTRL_gain/RADandRAvsCTRL_gain_AC_distal_withPeaksInPromoter.csv')
RADandRAvsCTRL_gain_AC_distal_withPeaksInPromoter = makeGRangesFromDataFrame(RADandRAvsCTRL_gain_AC_distal_withPeaksInPromoter, keep.extra.columns=TRUE)
length(unique(RADandRAvsCTRL_gain_AC_distal_withPeaksInPromoter$symbol))
newStyle <- mapSeqlevels(seqlevels(RADandRAvsCTRL_gain_AC_distal_withPeaksInPromoter), "UCSC")
RADandRAvsCTRL_gain_AC_distal_withPeaksInPromoter <- renameSeqlevels(RADandRAvsCTRL_gain_AC_distal_withPeaksInPromoter, newStyle)

diff = setdiff(diff, RADandRAvsCTRL_gain_AC_distal_withPeaksInPromoter)

#chequear que no haya solapamientos
findOverlaps(RADandRAvsCTRL_gain_AC_distal_withPeaksInPromoter, diff)
findOverlaps(RADandRAvsCTRL_loss_AC_distal_withPeaksInPromoter, diff)

diff_sequence <- diff %>%
  # Get a list of chip peaks 
  # look up the DNA sequence of each peak within each group
  get_sequence(HS.genome)

writeXStringSet(unique(diff_sequence), filepath = paste0("/media/usuario/Datos/eve/ifibyne/1_Genomic/2_OurProjects/1_GR-RAR_MicaSil/2_ATACseq/3_Analysis/Cruza_con_RNAseq/ChIPpeakAnno/A+C/RADandRAvsCTRL_loss/Motifs/CONTROL.fasta"), format = "fasta")

# Run STREME -> no lo pude hacer por un error, lo corro en la web
# streme = runStreme(unique(by_binding), control = unique(diff))




#############################################
#           PICOS EN EL PROMOTOR            #
#############################################

# Quedarme con los picos en promotores
Promoters_CTRL_filterAC = plyranges::filter(AnnotationChIPpeakAnnoCTRL, (distancetoFeature > -2000 & distancetoFeature < 500))
Promoters_RA_filterAC = plyranges::filter(AnnotationChIPpeakAnnoRA, (distancetoFeature > -2000 & distancetoFeature < 500))
Promoters_RAD_filterAC = plyranges::filter(AnnotationChIPpeakAnnoRAD, (distancetoFeature > -2000 & distancetoFeature < 500))

# Filtrar las anotaciones de todo el genoma por los genes que tienen picos diferenciales entre los tratamientos a analizar en los genes AC
Promoters_CTRL_filterAC = plyranges::filter(Promoters_CTRL_filterAC, symbol %in% RADandRAvsCTRL_gain_AC_distal_withPeaksInPromoter$symbol)
Promoters_RA_filterAC = plyranges::filter(Promoters_RA_filterAC, symbol %in% RADandRAvsCTRL_gain_AC_distal_withPeaksInPromoter$symbol)
Promoters_RAD_filterAC = plyranges::filter(Promoters_RAD_filterAC, symbol %in% RADandRAvsCTRL_gain_AC_distal_withPeaksInPromoter$symbol)
