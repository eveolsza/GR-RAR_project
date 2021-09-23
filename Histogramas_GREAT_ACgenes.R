library(tidyverse) 

path = '/media/usuario/Datos/eve/ifibyne/1_Genomic/2_OurProjects/1_GR-RAR_MicaSil/2_ATACseq/3_Analysis/Cruza_con_RNAseq/ChIPpeakAnno/A+C/'

################################################

histogram = function(folder, name){
file <- read.csv(paste0(path,folder,name,'.csv'))

breaks=c('< -5000',
         '-5000 to 500',
         '-500 to 0',
         '0 to 500',
         '500 to 3000',
         '> 3000')

conteo = c(
sum(file$distancetoFeature < -5000),
sum(file$distancetoFeature > -5000 & file$distancetoFeature < -500),
sum(file$distancetoFeature > -500 & file$distancetoFeature < 0),
sum(file$distancetoFeature > 0 & file$distancetoFeature < 500),
sum(file$distancetoFeature > 500 & file$distancetoFeature < 3000),
sum(file$distancetoFeature > 3000)
)


pdf(paste0(path,folder,name,'_hist.pdf'))

xx = barplot(height=conteo, names.arg=breaks, 
        xlab="Distance to TSS", 
        ylab="Regions count", 
        main=name, 
        col='white',
        xaxt="n",
        ylim=c(0,max(conteo)+2))
text(x = xx, y = conteo, label = conteo, pos = 3, cex = 0.8)
text(cex=1, x=xx, y=-0.25, breaks, xpd=TRUE, srt=0)

dev.off()
}

# para estos dos usé ylim=c(0,max(conteo)+20)) y text(cex=1, x=xx, y=-2, breaks, xpd=TRUE, srt=0)
histogram('RADandRAvsCTRL_gain/', 'RADandRAvsCTRL_gain_AnnotationChIPpeakAnno_FilterAC')
histogram('RADandRAvsCTRL_loss/', 'RADandRAvsCTRL_loss_AnnotationChIPpeakAnno_FilterAC')
# para estos dos usé ylim=c(0,max(conteo)+2)) y text(cex=1, x=xx, y=-0.25, breaks, xpd=TRUE, srt=0)
histogram('RADvsRA_gain/', 'RADvsRA_gain_AnnotationChIPpeakAnno_FilterAC')
histogram('RADvsRA_loss/', 'RADvsRA_loss_AnnotationChIPpeakAnno_FilterAC')


#####################################################

RegionsPerGeneHist = function(folder, name){
        file <- read.csv(paste0(path,folder,name,'.csv'))
        
        pdf(paste0(path,folder,name,'_RegionsPerGeneHist.pdf'))
        
        df  = file %>% group_by(symbol) %>% tally()
        df2 = df %>% group_by(n) %>% count()
        
        xx = barplot(height=df2$nn, names.arg=df2$n,
                xlab="Regions per Gene",
                ylab="Genes",
                col='white',
                main=name,
           #     breaks=seq(1,max(df2$n),1),
                ylim=c(0,max(df2$nn)+20)) 
        text(x = xx, y = df2$nn, label = df2$nn, pos = 3, cex = 0.8)
       
        dev.off()
}


RegionsPerGeneHist('RADandRAvsCTRL_gain/', 'RADandRAvsCTRL_gain_AnnotationChIPpeakAnno_FilterAC')
RegionsPerGeneHist('RADandRAvsCTRL_loss/', 'RADandRAvsCTRL_loss_AnnotationChIPpeakAnno_FilterAC')
RegionsPerGeneHist('RADvsRA_gain/', 'RADvsRA_gain_AnnotationChIPpeakAnno_FilterAC')
RegionsPerGeneHist('RADvsRA_loss/', 'RADvsRA_loss_AnnotationChIPpeakAnno_FilterAC')

