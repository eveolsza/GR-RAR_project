PROMOTERS_ALL_TXDB = promoters(txdb, upstream=2000, downstream=500)
annoData[10929]

keys = c("CCDC91", "SRSF8")

myGeneSymbols <- select(org.Hs.eg.db,
                        keys = keys,
                        columns = c("SYMBOL","ENTREZID"),
                        keytype = "SYMBOL")

myGeneSymbols

myGeneSymbolsTx <- select(TxDb.Hsapiens.UCSC.hg19.knownGene,
                          keys = myGeneSymbols$ENTREZID,
                          columns = c("TXID", "TXCHROM", "TXSTART", "TXEND"),
                          keytype = "GENEID")
myGeneSymbolsTx

plyranges::filter(PROMOTERS_ALL_TXDB, tx_id == 45990)
