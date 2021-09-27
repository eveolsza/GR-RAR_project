format_granges <- function(file) {
  archivo <- file
  names <- paste0(deparse(substitute(file)),'_',seq(length(file)))
  archivo <- setNames(archivo, names)
  score(archivo) <- archivo$fc
  newStyle <- mapSeqlevels(seqlevels(archivo), "UCSC")
  archivo <- renameSeqlevels(archivo, newStyle)
  return(archivo)
}

GREAT_ANNOTATION <- function(file){
  
  GREAT <- submitGreatJob(file,  
                          species = "hg19",
                          rule='basalPlusExt',
                          adv_upstream=5,
                          adv_downstream=3,
                          adv_span=1000,
                          includeCuratedRegDoms=TRUE)
  
  GREAT_tables = getEnrichmentTables(GREAT, download_by = "tsv")
  GREAT_association = plotRegionGeneAssociationGraphs(GREAT)
  
  return(GREAT_association)
} 

format_granges2 <- function(file) {
  archivo <- file
  names <- paste0(deparse(substitute(file)),'_',seq(length(file)))
  archivo <- setNames(archivo, names)
  newStyle <- mapSeqlevels(seqlevels(archivo), "NCBI")
  archivo <- renameSeqlevels(archivo, newStyle)
  return(archivo)
}

GRangestoBED <- function(GRangesfile){
  df <- data.frame(seqnames=seqnames(GRangesfile), 
                   starts=start(GRangesfile), ends=end(GRangesfile), 
                   names=paste0(deparse(substitute(GRangesfile)),'_',seq(length(GRangesfile))))
  name <- paste0(path,deparse(substitute(GRangesfile)),".bed")
  write.table(df, name, row.names=FALSE, quote=FALSE, sep = "    ")
}

GRangestoCSV <- function(GRangesfile){
  df <- data.frame(seqnames=seqnames(GRangesfile), 
                   starts=start(GRangesfile), ends=end(GRangesfile), 
                   names=paste0(deparse(substitute(GRangesfile)),'_',seq(length(GRangesfile))),
                   gene=GRangesfile$gene,
                   distTSS=GRangesfile$distTSS
                   )
  name <- paste0(path,deparse(substitute(GRangesfile)),".csv")
  write.table(df, name, row.names=FALSE, quote=FALSE, sep = ",")
}


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

plotRanges <- function(x, xlim = x, main = deparse(substitute(x)),
                       col = "black", sep = 0.5, ...)
{
  height <- 1
  if (is(xlim, "Ranges"))
    xlim <- c(min(start(xlim)), max(end(xlim)))
  bins <- disjointBins(IRanges(start(x), end(x) + 1))
  plot.new()
  plot.window(xlim, c(0, max(bins)*(height + sep)))
  ybottom <- bins * (sep + height) - height
  rect(start(x)-0.5, ybottom, end(x)+0.5, ybottom + height, col = col, ...)
  title(main)
  axis(1)
}
