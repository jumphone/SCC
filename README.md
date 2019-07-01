# SCC
Single-cell Communication

<img src="https://github.com/jumphone/Bioinformatics/raw/master/scRNAseq/try_20190424/src/SCC_LOGO.png" width="200">

### Single-Cell Communication (SCC) Toolkit 

Author: Feng Zhang
 
* New feature is only for Seurat3

### Prepare:

    hash, Seurat, Rtsne


### Usage:

    library('Seurat')
    library(dplyr)
    library(Matrix)
    
    source('https://raw.githubusercontent.com/jumphone/SCC/master/SCC.R')
    
    pbmc=load('pbmc.RDS') # load Seurat Object
    
    
    ORITAG=as.character(pbmc@active.ident)   
    VEC=pbmc@reductions$umap@cell.embeddings
    pbmc.data=as.matrix(pbmc@assays$RNA@scale.data)
    
    SAVE_DIR="./"
    
    used_cell=
    
    used_gene=
    
    USEDC=which(colnames(pbmc.data) %in% used_cell)
    USEDG=which(rownames(pbmc.data) %in% used_gene)
    
    pbmc.data=pbmc.data[USEDG,USEDC]
    VEC=VEC[USEDC,]
    ORITAG=ORITAG[USEDC]
    
    #########################################
    
    
    
    library('Rtsne')
    T1=Rtsne(VEC,dims=1)
    ONE=T1$Y  
    saveRDS(ONE,file=paste0(SAVE_DIR,'/','ONE.RDS'))
    
    NUM=100
    
    OUT=getBIN(ONE,NUM=NUM)
    BIN=OUT$BIN
    BINTAG=OUT$TAG
    saveRDS(BIN,file=paste0(SAVE_DIR,'/','BIN.RDS'))
    saveRDS(BINTAG,file=paste0(SAVE_DIR,'/','BINTAG.RDS'))
    
    pbmc@meta.data$bin=rep(NA, ncol(pbmc))
    pbmc@meta.data$bin[USEDC]=BINTAG
    pdf(paste0(SAVE_DIR,'/','1ID.pdf'),width=12,height=10)
    DimPlot(pbmc,group.by='bin',reduction.use='umap', label=T)
    dev.off()
    
<img src="https://github.com/jumphone/Bioinformatics/raw/master/scRNAseq/try_20190424/src/ID.png" width="300">

    LR=read.table('ReceptorLigand.txt.mouse',header=T,sep='\t')
    #https://github.com/jumphone/Bioinformatics/tree/master/scRNAseq/RecLig/
    
    EXP=pbmc.data
    
    MEAN=getMEAN(EXP, LR, NUM=NUM)
    saveRDS(MEAN,file=paste0(SAVE_DIR,'/','MEAN.RDS'))
        
    PMAT=getPMAT(EXP, LR, BIN, MEAN)
    saveRDS(PMAT,file=paste0(SAVE_DIR,'/','PMAT.RDS'))
    
    pdf(paste0(SAVE_DIR,'/','GCOR.pdf'),width=20,height=20)
    OUT=getPmatHEAT(PMAT,SHOW=T)
    dev.off()
    HEAT=OUT$HEAT
    DIST=OUT$DIST
    ORDER=HEAT$colInd
    
    pdf(paste0(SAVE_DIR,'/','2CLUST.pdf'),width=20,height=20)
    CLUST=getCLUST(ORDER, DIST, CCUT=0.7, SHOW=T)
    dev.off()
    

<img src="https://github.com/jumphone/Bioinformatics/raw/master/scRNAseq/try_20190424/src/CLUST.png" width="300">

    MLR=getMLR(CLUST, LR, PMAT)
    LR=MLR[,c(1:2)]

    CMAT=getCMAT(EXP,LR,PMAT,BI=TRUE)
    saveRDS(CMAT,file=paste0(SAVE_DIR,'/','CMAT.RDS'))
    
    pdf(paste0(SAVE_DIR,'/','3CMAT.pdf'),width=15,height=13)
    library('gplots')
    heatmap.2(log(CMAT+1,10),scale=c("none"),dendrogram='both',Colv=T,Rowv=T,trace='none',
      col=colorRampPalette(c('blue3','grey95','red3')) ,margins=c(10,15))
    heatmap.2(CMAT,scale=c("none"),dendrogram='none',Colv=F,Rowv=F,trace='none',
      col=colorRampPalette(c('blue3','grey95','red3')) ,margins=c(10,15))
    dev.off()

<img src="https://github.com/jumphone/Bioinformatics/raw/master/scRNAseq/try_20190424/src/CMAT.png" width="300">

    OUT=getPAIR(CMAT)
    PAIR=OUT$PAIR
    SCORE=OUT$SCORE
    RANK=OUT$RANK
    saveRDS(PAIR,file=paste0(SAVE_DIR,'/','PAIR.RDS'))
    
    pdf(paste0(SAVE_DIR,'/','4CPlot_TOP200.pdf'),width=12,height=10)
    CPlot(VEC,PAIR[1:200,],BINTAG)
    dev.off()

<img src="https://github.com/jumphone/Bioinformatics/raw/master/scRNAseq/try_20190424/src/CPlot.png" width="300">
    
    
    
    NET=getNET(PAIR, BINTAG,ORITAG )
    write.table(NET,file=paste0(SAVE_DIR,'/','NET.txt'),sep='\t',row.names=F,col.names=T,quote=F)
       
    CN=getCN(NET)
    pdf(paste0(SAVE_DIR,'/','5DPlot.pdf'),width=100,height=100)
    DP=DPlot(NET, CN, COL=3)
    dev.off()
    
    ADP=p.adjust(DP,method='fdr')
    
    IDP=ADP
    IDP[which(IDP==0)]=min(IDP[which(IDP>0)])/2
     
    DD=sort(-log(IDP,10),decreasing=T)
    CC=rep('grey',length(DD))
    CC[which(DD> -log(0.05,10))]='red'
    pdf(paste0(SAVE_DIR,'/','PVALUE.pdf'),width=100,height=100)
    par(mar=c(20,5,5,5))
    barplot(DD,las=2,ylab='-log10(adjusted p-value)',col=CC)
    dev.off()
       
<img src="https://github.com/jumphone/Bioinformatics/raw/master/scRNAseq/try_20190424/src/DPlot.png" width="300">
    
    
    
    CNET=c()
    CNETSCORE=c()
    TMP=DD[which(CC=='red')]
    i=1
    while(i <=length(TMP)){
        this_source=unlist(strsplit(names(TMP)[i],split='_to_'))[1]
        this_target=unlist(strsplit(names(TMP)[i],split='_to_'))[2]
        CNET=cbind(CNET, c(this_source,this_target))
        CNETSCORE=c(CNETSCORE, TMP[i])
        
        i=i+1}
    
    OUT=cbind(t(CNET),CNETSCORE)
    colnames(OUT)=c('source','target','score_neg_log_adp')
    write.table(OUT,file=paste0(SAVE_DIR,'/','CNET.txt'),quote=F,sep='\t',row.names=F,col.names=F)
    



    #SIG_INDEX=which(DP<0.05)    
    SIG_INDEX=which(ADP<0.05)
    SIG_PAIR=names(SIG_INDEX)
    TOP_NET=NET
    #TOP_NET=getNET(PAIR[1:500,], BINTAG,ORITAG )
    
    pdf(paste0(SAVE_DIR,'/','6LPlot.pdf'),width=100,height=100)
    OUT=c()
    #OUT_TYPE=c()
    RCN=trunc(sqrt(length(SIG_PAIR))+1)
    par(mfrow=c(RCN,RCN))
    i=1
    while(i<= length(SIG_PAIR) ){
        this_pair=SIG_PAIR[i]
        LT=unlist(strsplit(this_pair, "_to_"))[1]
        RT=unlist(strsplit(this_pair, "_to_"))[2]
        try({
        LP=LPlot(LT, RT, TOP_NET, PMAT,LR, MAIN=paste0(as.character(SIG_INDEX[i]),' ',SIG_PAIR[i]),SEED=12345,PCUT=0.05)    
        #########################
        this_out_index=which(LP[,1]>-log(0.05,10) & LP[,2]>-log(0.05,10))
        this_out=t(LP)[,c(this_out_index,this_out_index)]
        
        if(length(this_out_index)>0){ 
            colnames(this_out)=paste0(SIG_PAIR[i],'_|_',colnames(this_out))
            OUT=cbind(OUT,this_out)}
        #########################
        colnames(LP)=paste0(c('Lexp','Rexp'),'_',c(LT,RT))
        write.table(LP,file=paste0(SAVE_DIR,'/',as.character(SIG_INDEX[i]),'.tsv'),row.names=T,col.names=T,sep='\t',quote=F)
        })
        print(i)
        i=i+1}
    dev.off()
    
<img src="https://github.com/jumphone/Bioinformatics/raw/master/scRNAseq/try_20190424/src/LPlot.png" width="300">    
    
       
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    ############################################
    OUT=t(OUT)
    OUT=unique(cbind(OUT,rownames(OUT)))
    
    get_LT<-function(X){
        X=unlist(strsplit(X, "_\\|_"))[1]
        X=unlist(strsplit(X, "_to_"))[1]
        return(X)
        }
    get_RT<-function(X){
        X=unlist(strsplit(X, "_\\|_"))[1]
        X=unlist(strsplit(X, "_to_"))[2]
        return(X)
        }
        
    OUT_LT=apply(matrix(OUT[,3],ncol=1),1,get_LT)
    OUT_RT=apply(matrix(OUT[,3],ncol=1),1,get_RT)
    
    OUT=cbind(OUT,OUT_LT,OUT_RT)
    
 
    
    
    
    
    library(plotly)
    
    COL=rep('rgb(230,230,230)',nrow(OUT))
    COL[which(OUT_LT %in% c('Tumor Cells'))]='green'
    COL[which(OUT_RT %in% c('Tumor Cells'))]='blue'
    p <- plot_ly(type = 'scatter', mode = 'markers') %>%
    add_trace(
    x = OUT[,2], 
    y = OUT[,1],
    text = OUT[,3],
    hoverinfo = 'text',
    marker = list(color=COL),
    showlegend = F
    )
    htmlwidgets::saveWidget(as_widget(p), paste0(SAVE_DIR,'/',"TumorLR.html"))
    
    library(plotly)
    
    COL=rep('rgb(230,230,230)',nrow(OUT))
    COL[which(OUT_LT %in% c('Proliferating Cells'))]='green'
    COL[which(OUT_RT %in% c('Proliferating Cells'))]='blue'
    p <- plot_ly(type = 'scatter', mode = 'markers') %>%
    add_trace(
    x = OUT[,2], 
    y = OUT[,1],
    text = OUT[,3],
    hoverinfo = 'text',
    marker = list(color=COL),
    showlegend = F
    )
    htmlwidgets::saveWidget(as_widget(p), paste0(SAVE_DIR,'/',"ProliferatingLR.html"))
    
    
    
    


    TAG=groupTAG(BINTAG,LT="LGroup",RT='RGroup',LC=c(1,2,3),RC=c(4,5,6))
    TMP=getNET(PAIR, BINTAG, TAG )
    LPlot(LT='LGroup', RT='RGroup', TMP, PMAT,SEED=123)    
   
    

    
    
    
    
    
    
    
    
    




