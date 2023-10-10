.extraMitoGenes=function(organism){
  #organism: Human, Mouse
  
  if(tolower(organism)=="human"){
    gns=.extraHumanGeneAnnoAdderFn()
  } else if (tolower(organism)=="mouse"){
    gns=.extraMouseGeneAnnoAdderFn()
  } else if(tolower(organism)=="macaque"){
    gns=.extraMacaqueGeneAnnoAdderFn()
  } else {
    stop("Wrong organism name!")
  }
  
  gns=gns[gns$seqnames=="MT",]
  gns2=unlist(gns$entrezid)
  gns$entrezgene_id=gns2
  gns$ensembl_gene_id=gns$gene_id
  return(gns)
}

.extrarRNAGenes=function(organism){
  #organism: Human, Mouse
  
  if(tolower(organism)=="human"){
    gns=.extraHumanGeneAnnoAdderFn()
  } else if (tolower(organism)=="mouse"){
    gns=.extraMouseGeneAnnoAdderFn()
  } else if(tolower(organism)=="macaque"){
    gns=.extraMacaqueGeneAnnoAdderFn()
  } else {
    stop("Wrong organism name!")
  }
  
  gns=gns[which(gns$gene_biotype=="rRNA"),]
  gns2=unlist(gns$entrezid)
  gns$entrezgene_id=gns2
  gns$ensembl_gene_id=gns$gene_id
  return(gns)
}

.extraRibosomalProteinGenes=function(organism){
  #organism: Human, Mouse
  
  if(tolower(organism)=="human"){
    gns=.extraHumanGeneAnnoAdderFn()
  } else if (tolower(organism)=="mouse"){
    gns=.extraMouseGeneAnnoAdderFn()
  } else if(tolower(organism)=="macaque"){
    gns=.extraMacaqueGeneAnnoAdderFn()
  } else {
    stop("Wrong organism name!")
  }
  
  gns=gns[which(grepl("^RP[LS]",gns$gene_short_name,ignore.case = T)),]
  gns2=unlist(lapply(gns$entrezid,function(x) x[1]))
  gns$entrezgene_id=gns2
  gns$ensembl_gene_id=gns$gene_id
  return(gns)
}

.extraOXPHOSGenes=function(organism,inputGeneset=NULL){
  #organism: Human, Mouse
  
  if(tolower(organism)=="human"){
    gns=.extraHumanGeneAnnoAdderFn()
  } else if (tolower(organism)=="mouse"){
    gns=.extraMouseGeneAnnoAdderFn()
  } else if(tolower(organism)=="macaque"){
    gns=.extraMacaqueGeneAnnoAdderFn()
  } else {
    stop("Wrong organism name!")
  }
  
  if(is.null(inputGeneset)){
      if(!file.exists("~/serverFiles/OXPHOS_gene_symbols.txt")){
        stop("Error! file ~/serverFiles/OXPHOS_gene_symbols.txt is missing")
      }
      
      oxphosGenes=read.table("~/serverFiles/OXPHOS_gene_symbols.txt",sep="\t",stringsAsFactors = F)
    
  }
  
  if(sum(colnames(gns)=="gene_short_name")>0){
    gns=gns[which(toupper(gns$gene_short_name) %in% toupper(oxphosGenes[,1])),]
  } else {
    gns=gns[which(toupper(gns$symbol) %in% toupper(oxphosGenes[,1])),]
  }
  
  
  return(gns)
}

.extraIEGGenes=function(organism){
  #organism: Human, Mouse
  
  gns=.extraMouseGeneAnnoAdderFn()
  
  if(!file.exists("~/serverFiles/IEG_gene_symbols")){
    #library(googleCloudStorageR)
    stop("Error! file serverFiles/IEG_gene_symbols is missing")
  }
  IEGenes=read.table("~/serverFiles/IEG_gene_symbols",sep="\t",stringsAsFactors = F)
  
  gns=gns[tolower(gns$gene_name) %in% tolower(IEGenes$V1),]
  
  if(tolower(organism)=="human"){
    
    {
      #library(googleCloudStorageR)
      if(!file.exists("~/serverFiles/ortholog_mapping.rda")){
        stop("Error! file ~/serverFiles/ortholog_mapping.rda is missing")
        #gcs_get_object("vgazesta/serverFiles/orthologsFeb3/ortholog_mapping.rda", 
        #                 saveToDisk = "~/serverFiles/ortholog_mapping.rda",overwrite=T)
      }
      
      load("~/serverFiles/ortholog_mapping.rda")
    }
    
    
    
    mapping=mapping[mapping$mouse %in% gns$ensembl_gene_id,]
    gnsHuman=.extraHumanGeneAnnoAdderFn()
    gns=gnsHuman[gnsHuman$ensembl_gene_id %in% mapping$human,]
  } else if (tolower(organism)=="mouse"){
    
  } else if(tolower(organism)=="macaque"){
    gns=.extraMacaqueGeneAnnoAdderFn()
    gns=gns[which(tolower(gns$symbol) %in% tolower(IEGenes$V1)|tolower(gns$gene_name) %in% tolower(IEGenes$V1)),]
  } else {
    stop("Wrong organism name!")
  }
  
  gns$ensembl_gene_id=gns$gene_id
  return(gns)
}

.extraCellCycleGenes=function(organism){
  if(tolower(organism)=="human"){
    gns=.extraHumanGeneAnnoAdderFn()
  } else if (tolower(organism)=="mouse"){
    gns=.extraMouseGeneAnnoAdderFn()
  } else if(tolower(organism)=="macaque"){
    gns=.extraMacaqueGeneAnnoAdderFn()
  } else {
    stop("Wrong organism name!")
  }
  
  cc.s=Seurat::cc.genes.updated.2019$s.genes
  cc.g2m=Seurat::cc.genes.updated.2019$g2m.genes
  
  cc.g2m[cc.g2m=="PIMREG"]="FAM64A"
  cc.g2m[cc.g2m=="JPT1"]="HN1"
  
  cc.s=gns[toupper(gns$symbol) %in% toupper(cc.s),]
  cc.g2m=gns[toupper(gns$symbol) %in% toupper(cc.g2m),]
  
  cc.s$ensemble_gene_id=row.names(cc.s)
  cc.g2m$ensemble_gene_id=row.names(cc.g2m)
  
  return(list(cc.s=cc.s,cc.g2m=cc.g2m))
}

.extraHumanGeneAnnoAdderFn=function(inputGeneNames=NULL){
  #require(EnsDb.Hsapiens.v75)
  require(EnsDb.Hsapiens.v86)
  
  if(!dir.exists("~/serverFiles")){
    dir.create("~/serverFiles",recursive = T)
  }
  
  gns <- as.data.frame(genes(EnsDb.Hsapiens.v86))
  gns$gene_short_name=gns$gene_name
  gns$symbol=toupper(gns$symbol)
  gns$ensembl_gene_id=row.names(gns)
  tst=unique(as.character(gns$gene_biotype))
  tst=setdiff(tst,c("protein_coding","lincRNA"))
  tst=factor(as.character(gns$gene_biotype),levels=c("protein_coding","lincRNA",tst))
  gns=gns[order(tst,decreasing = F),]
  if(!is.null(inputGeneNames)){
    rwNames=toupper(inputGeneNames)
    psCols=c("gene_short_name","ensembl_gene_id")
    slCounts=0
    slCol=""
    if(sum(grepl("\\.",rwNames)&grepl("^ENS",rwNames))>0){
      rwNames=strsplit(rwNames,"\\.")
      rwNames=unlist(lapply(rwNames,function(x)x[1]))
    }
    {
      
      if(!file.exists("~/serverFiles/human_map_to_ensembl.rda")){
        stop("Error! file ~/serverFiles/human_map_to_ensembl.rda is missing")
      }
      
      load("~/serverFiles/human_map_to_ensembl.rda")
    }
    
    map_to_ensmbl$source=toupper(map_to_ensmbl$source)
    
    if(!file.exists("~/serverFiles/human_mapping_hg19.rda")){
      stop("Error! file ~/serverFiles/human_mapping_hg19.rda is missing")
    }
    
    load("~/serverFiles/human_mapping_hg19.rda")
    human_hg19$source=toupper(human_hg19$source)
    
    if(sum(toupper(rwNames) %in% human_hg19$source) > sum(toupper(rwNames) %in% map_to_ensmbl$source)){
      map_to_ensmbl=rbind(human_hg19,map_to_ensmbl)
      map_to_ensmbl=map_to_ensmbl[!duplicated(map_to_ensmbl$source),]
    } else {
      map_to_ensmbl=rbind(map_to_ensmbl,human_hg19)
      map_to_ensmbl=map_to_ensmbl[!duplicated(map_to_ensmbl$source),]
    }
    
    c_map=data.frame(source=gns$gene_short_name,target=gns$ensembl_gene_id,stringsAsFactors = F)
    c_map=rbind(c_map,data.frame(source=gns$ensembl_gene_id,target=gns$ensembl_gene_id,stringsAsFactors = F))
    c_map=c_map[!is.na(c_map$ensembl_gene_id),]
    map_to_ensmbl=rbind(map_to_ensmbl,c_map)
    map_to_ensmbl=map_to_ensmbl[!duplicated(map_to_ensmbl$source),]
    map_to_ensmbl=merge(map_to_ensmbl,data.frame(source=toupper(rwNames),stringsAsFactors = F),by="source",all.y=T)
    gns=merge(gns,map_to_ensmbl,by.x="ensembl_gene_id",by.y="target",all.y=T)
    
    gns=gns[match(rwNames,gns$source),]
    row.names(gns)=inputGeneNames
    gns$gene_id=inputGeneNames
    gns=gns[,-which(colnames(gns) %in% c("source","target"))]
  }
  
  return(gns)
}

.extraMouseGeneAnnoAdderFn=function(inputGeneNames=NULL){
  require(EnsDb.Mmusculus.v79)
  
  if(!dir.exists("~/serverFiles")){
    dir.create("~/serverFiles",recursive = T)
  }
  
  gns <- as.data.frame(genes(EnsDb.Mmusculus.v79))
  gns$symbol=toupper(gns$symbol)
  gns$gene_short_name=gns$gene_name
  gns$ensembl_gene_id=row.names(gns)
  tst=unique(as.character(gns$gene_biotype))
  tst=setdiff(tst,c("protein_coding","lincRNA"))
  tst=factor(as.character(gns$gene_biotype),levels=c("protein_coding","lincRNA",tst))
  gns=gns[order(tst,decreasing = F),]
  
  if(!is.null(inputGeneNames)){
    rwNames=toupper(inputGeneNames)
    
    if(sum(grepl("\\.",rwNames)&grepl("^ENS",rwNames))>0){
      rwNames=strsplit(rwNames,"\\.")
      rwNames=unlist(lapply(rwNames,function(x)x[1]))
    }
    {
      #library(googleCloudStorageR)
      if(!file.exists("~/serverFiles/mouse_map_to_ensembl.rda")){
        stop("Error! file ~/serverFiles/mouse_map_to_ensembl.rda is missing")
      }
      
      load("~/serverFiles/mouse_map_to_ensembl.rda")
      
    }
    map_to_ensmbl$source=toupper(map_to_ensmbl$source)
    
    
    c_map=data.frame(source=gns$gene_short_name,target=gns$ensembl_gene_id,stringsAsFactors = F)
    c_map=rbind(c_map,data.frame(source=gns$ensembl_gene_id,target=gns$ensembl_gene_id,stringsAsFactors = F))
    c_map=c_map[!is.na(c_map$ensembl_gene_id),]
    map_to_ensmbl=rbind(map_to_ensmbl,c_map)
    map_to_ensmbl=map_to_ensmbl[!duplicated(map_to_ensmbl$source),]
    map_to_ensmbl=merge(map_to_ensmbl,data.frame(source=toupper(rwNames),stringsAsFactors = F),by="source",all.y=T)
    
    
    gns=merge(gns,map_to_ensmbl,by.x="ensembl_gene_id",by.y="target",all.y=T)
    gns=gns[match(rwNames,gns$source),]
    row.names(gns)=inputGeneNames
    gns=gns[,-which(colnames(gns) %in% c("source","target"))]
  }
  
  
  return(gns)
}

.extraMacaqueGeneAnnoAdderFn=function(inputGeneNames=NULL){
  library(AnnotationHub)
  hub <- AnnotationHub()
  avail.resources=query(hub, c("fascicularis"))
  avail.resources=data.frame(id=avail.resources$ah_id,
                              desc=avail.resources$description,
                              genome=avail.resources$genome,
                              title=avail.resources$title,
                              stringsAsFactors = F)
  avail.resources=avail.resources[grepl("Ensembl",avail.resources$desc),]
  ensdb <- hub[[avail.resources$id[nrow(avail.resources)]]]
  
  if(!dir.exists("~/serverFiles")){
    dir.create("~/serverFiles",recursive = T)
  }
  
  gns <- as.data.frame(genes(ensdb))
  gns$symbol=toupper(gns$symbol)
  gns$gene_short_name=gns$gene_name
  gns$ensembl_gene_id=row.names(gns)
  
  if(!is.null(inputGeneNames)){
    rwNames=toupper(inputGeneNames)
    
    if(sum(grepl("\\.",rwNames)&grepl("^ENS",rwNames))>0){
      rwNames=strsplit(rwNames,"\\.")
      rwNames=unlist(lapply(rwNames,function(x)x[1]))
    }
    
    sl_cols=c('gene_id','gene_name',"symbol")
    map_to_ensmble=NULL
    for(icol in sl_cols){
      map_to_ensmble=rbind(map_to_ensmble,data.frame(source=gns[,icol],ensembl_gene_id=gns$gene_id,stringsAsFactors = F))
    }
    map_to_ensmble=map_to_ensmble[!duplicated(map_to_ensmble$source),]
    map_to_ensmble$source=toupper(map_to_ensmble$source)
    
    map_to_ensmble=merge(map_to_ensmble,data.frame(source=toupper(rwNames),stringsAsFactors = F),by="source",all.y=T)
    
    
    gns=merge(map_to_ensmble,gns,by.x="ensembl_gene_id",by.y="gene_id",all.y=T)
    gns=gns[match(rwNames,gns$source),]
    row.names(gns)=inputGeneNames
    if(sum(colnames(gns) %in% c("source","target"))>0){
      gns=gns[,-which(colnames(gns) %in% c("source","target"))]
    }
    
  }
  
  gns$gene_short_name=gns$symbol
  return(gns)
}

.extraQCAnnoAdderFn=function(inputExpSet,organism,ncores=5){
  
  {
    x=counts(inputExpSet)
    
    if(ncol(x)>10000){
      
      if(ncores==1){
        nGene=c()
        nUMI=c()
        for(i in seq(1,ncol(x),10000)){
          nGene=c(nGene,apply(x[,i:min(i+10000-1,ncol(x))],2,function(x) sum(x>0)))
          nUMI=c(nUMI,apply(x[,i:min(i+10000-1,ncol(x))],2,function(x) sum(x)))
        }
      } else {
        nGene=c()
        nUMI=c()
        nGeneUMI=parallel::mclapply(seq(1,ncol(x),10000),function(i,x){
          tmp_nGene=apply(x[,i:min(i+10000-1,ncol(x))],2,function(x) sum(x>0))
          tmp_nUMI=apply(x[,i:min(i+10000-1,ncol(x))],2,function(x) sum(x))
          return(list(nGene=tmp_nGene,nUMI=tmp_nUMI,i=i))
        },x=x,mc.cores = ncores)
        ilist=c()
        for(i in 1:length(nGeneUMI)){
          nGene=c(nGene,nGeneUMI[[i]]$nGene)
          nUMI=c(nUMI,nGeneUMI[[i]]$nUMI)
          ilist=c(ilist,nGeneUMI[[i]]$i)
        }
        
        if(any(ilist!=seq(1,ncol(x),10000),na.rm = F)){
          stop("Error in calculation of nGene/nUMI!")
        }
      }
      
    } else {
      nGene=apply(x,2,function(x) sum(x>0))
      nUMI=apply(x,2,function(x) sum(x))
    }
    
    inputExpSet$QC_Gene_total_count=nUMI
    inputExpSet$QC_Gene_unique_count=nGene
  }
  
  tmpExpData=counts(inputExpSet)
  
  inputExpSet$QC_MT.pct=.extraMitoPctFn(inputData=inputExpSet,
                                        organism = organism,
                                        recalculate_nUMI = F)
  inputExpSet$QC_IEG.pct=.extraIEGPctFn(inputData=inputExpSet,
                                        organism = organism,
                                        recalculate_nUMI = F)
  inputExpSet$QC_OXPHOS.pct=.extraOXPHOSPctFn(inputData = inputExpSet,
                                              organism = organism,
                                              recalculate_nUMI = F)
  inputExpSet$QC_rRNA.pct=.extrarRNAPctFn(inputData = inputExpSet,
                                          organism = organism,
                                          recalculate_nUMI = F)
  inputExpSet$QC_RibosomalProtein.pct=.extraRPctFn(inputData=inputExpSet,
                                                    organism = organism,
                                                    recalculate_nUMI = F)
                                                    
  tmpMedian=c()
  #tmpMean=c()
  tmpDetectionRate=c()
  #tmpVar=c()
  if(ncores==1|T){
    for(i in seq(1,nrow(tmpExpData),3000)){
      tmpMedian=c(tmpMedian,apply(tmpExpData[i:min(i+3000-1,nrow(tmpExpData)),],1,median))
      #tmpMean=c(tmpMean,apply(tmpExpData[i:min(i+3000-1,nrow(tmpExpData)),],1,mean))
      #tmpVar=c(tmpVar,apply(tmpExpData[i:min(i+3000-1,nrow(tmpExpData)),],1,var))
      tmpDetectionRate=c(tmpDetectionRate,apply(tmpExpData[i:min(i+3000-1,nrow(tmpExpData)),],1,function(x) sum(x>0)/length(x)))
    }
    
  } else {
    MedianDetectionRate=parallel::mclapply(seq(1,nrow(tmpExpData),3000),function(i,tmpExpData){
      tmpMedian=apply(tmpExpData[i:min(i+3000-1,nrow(tmpExpData)),],1,median)
      #tmpMean=c(tmpMean,apply(tmpExpData[i:min(i+3000-1,nrow(tmpExpData)),],1,mean))
      #tmpVar=c(tmpVar,apply(tmpExpData[i:min(i+3000-1,nrow(tmpExpData)),],1,var))
      tmpDetectionRate=apply(tmpExpData[i:min(i+3000-1,nrow(tmpExpData)),],1,function(x) sum(x>0)/length(x))
      return(list(median=tmpMedian,detectionRate=tmpDetectionRate,i=i))
    },tmpExpData=tmpExpData,mc.cores =ncores)
  }
  
  
  names(tmpMedian)=row.names(tmpExpData)
  rowData(inputExpSet)$QC_medianExp=round(tmpMedian,3)
  tmp=data.frame(gene=names(tmpMedian),exp=tmpMedian,stringsAsFactors = F)
  tmp=tmp[order(tmp$exp,decreasing = T),]
  tmp=tmp[which(tmp$exp>=tmp$exp[50]),]
  rowData(inputExpSet)$QC_top50_expressed="No"
  rowData(inputExpSet)$QC_top50_expressed[row.names(inputExpSet) %in% tmp$gene]="Yes"
  tmp=inputExpSet[row.names(inputExpSet) %in% tmp$gene,]
  inputExpSet$QC_top50_pct=colSums(assays(tmp)[["counts"]])/inputExpSet$QC_Gene_total_count*100#  colSums(tmpExpData)*100
  
  
  
  if(sum(colnames(rowData(inputExpSet))=="gene_biotype")==1){
    tmp=inputExpSet[which(rowData(inputExpSet)$gene_biotype=="lincRNA"),]
    
    tmp=colSums(counts(tmp))/inputExpSet$QC_Gene_total_count
    inputExpSet$QC_lncRNA_pct=tmp*100
  }
  
  #rowData(inputExpSet)$QC_meanExp=round(tmpMean,3)
  
  rowData(inputExpSet)$QC_detectionRate=round(tmpDetectionRate,3)
  
  #rowData(inputExpSet)$QC_varExp=round(tmpVar,3)
  
  #tst=cor(tmpExpData,method = "spearman")
  
  #colData(inputExpSet)$QC_raw_data_mean_cor=apply(tst,1,mean)
  #colData(inputExpSet)$QC_raw_data_median_cor=apply(tst,1,median)
  
  
  return(inputExpSet)
}

.extraMitoPctFn=function(inputData,organism,inputMTgenes=NULL,inputGeneName=NULL,recalculate_nUMI=T){
  
  rwNames=tolower(row.names(inputData))
  x=.extraMitoGenes(organism = organism)
  if(sum(grepl("\\.",rwNames))>0){
    rwNames=strsplit(rwNames,"\\.")
    rwNames=unlist(lapply(rwNames,function(x)x[1]))
  }
  if(is.null(inputGeneName)){
    tmpCols=c("ensembl_gene_id","entrezgene_id",'gene_name')
    
    slCounts=0
    slCol=""
    
    for(i in tmpCols){
      tmp=sum(rwNames %in% tolower(x[,i]))
      if(tmp>slCounts){
        slCounts=tmp
        slCol=i
      }
    }
    if(slCol!=""){
      inputMTgenes=tolower(x[,slCol])
    } else {
      inputMTgenes=""
    }
    
  } else {
    if(inputGeneName=="ensembl_gene_id"){
      inputMTgenes=tolower(x$ensembl_gene_id)
    } else {
      if(inputGeneName=="entrezgene_id"){
        x=.extraMitochondrialGenes()
        inputMTgenes=tolower(x$entrezgene_id)
      }
    }
  }
  
  print(paste("Number of MT genes in the dataset:",length(which(rwNames %in% inputMTgenes)),"/",sum(x$gene_biotype=="protein_coding")))
  
  inputMTgenes=inputMTgenes[inputMTgenes!=""]
  if(length(inputMTgenes)>0){
    tmpColSums=c()
    if(sum(colnames(colData(inputData))=='QC_Gene_total_count')==0|recalculate_nUMI){
      if(ncol(inputData)>10000){
        for(i in seq(1,ncol(inputData),10000)){
          tmpColSums=c(tmpColSums,apply(counts(inputData)[,i:min(i+10000-1,ncol(inputData))],2,sum))
        }
      } else {
        tmpColSums=apply(counts(inputData),2,sum)
      }
    } else {
      tmpColSums=inputData$QC_Gene_total_count
    }
    
    
    
    res=colSums(x = counts(inputData)[which(rwNames %in% inputMTgenes), , drop = FALSE])/tmpColSums
    res=res*100
  } else {
    res=rep(NA,ncol(inputData))
  }
  
  return(res)
}

.extraIEGPctFn=function(inputData,organism,inputIEGgenes=NULL,inputGeneName=NULL,recalculate_nUMI=T){
  
  rwNames=tolower(row.names(inputData))
  x=.extraIEGGenes(organism = organism)
  if(sum(colnames(x)=="entrezid")>0){
    colnames(x)[colnames(x)=="entrezid"]="entrezgene_id"
  }
  if(sum(grepl("\\.",rwNames))>0){
    rwNames=strsplit(rwNames,"\\.")
    rwNames=unlist(lapply(rwNames,function(x)x[1]))
  }
  if(is.null(inputGeneName)){
    tmpCols=c("ensembl_gene_id","entrezgene_id",'gene_name')
    
    slCounts=0
    slCol=""
    
    for(i in tmpCols){
      tmp=sum(rwNames %in% tolower(x[,i]))
      if(tmp>slCounts){
        slCounts=tmp
        slCol=i
      }
    }
    if(slCol!=""){
      inputIEGgenes=tolower(x[,slCol])
    } else {
      inputIEGgenes=""
    }
    
  } else {
    if(inputGeneName=="ensembl_gene_id"){
      inputIEGgenes=tolower(x$ensembl_gene_id)
    } else {
      if(inputGeneName=="entrezgene_id"){
        x=.extraMitochondrialGenes()
        inputIEGgenes=tolower(x$entrezgene_id)
      }
    }
  }
  
  print(paste("Number of IEG genes in the dataset:",length(which(rwNames %in% inputIEGgenes)),"/",sum(x$gene_biotype=="protein_coding")))
  
  inputIEGgenes=inputIEGgenes[inputIEGgenes!=""]
  if(length(inputIEGgenes)>0){
    if(sum(colnames(colData(inputData))=="QC_Gene_total_count")==0|recalculate_nUMI){
      tmpColSums=c()
      if(ncol(inputData)>10000){
        for(i in seq(1,ncol(inputData),10000)){
          tmpColSums=c(tmpColSums,apply(counts(inputData)[,i:min(i+10000-1,ncol(inputData))],2,sum))
        }
      } else {
        tmpColSums=apply(counts(inputData),2,sum)
      }
    } else {
      tmpColSums=inputData$QC_Gene_total_count
    }
    
    
    
    res=colSums(x = counts(inputData)[which(rwNames %in% inputIEGgenes), , drop = FALSE])/tmpColSums
    res=res*100
  } else {
    rep(NA,ncol(inputData))
  }
  
  return(res)
}

.extraOXPHOSPctFn=function(inputData,organism,inputGeneset=NULL,inputGeneName_col=NULL,recalculate_nUMI=T,standardize_gene_names=T){
  
  rwNames=tolower(row.names(inputData))
  x=.extraOXPHOSGenes(organism = organism,inputGeneset=inputGeneset)
  if(sum(grepl("\\.",rwNames))>0){
    rwNames=strsplit(rwNames,"\\.")
    rwNames=unlist(lapply(rwNames,function(x)x[1]))
  }
  if(is.null(inputGeneName_col)){
    tmpCols=c("ensembl_gene_id",'gene_short_name',"entrezid","gene_name")
    
    slCounts=0
    slCol=""
    
    for(i in tmpCols){
      tmp=sum(rwNames %in% tolower(x[,i]))
      if(tmp>slCounts){
        slCounts=tmp
        slCol=i
      }
    }
    if(slCol!=""){
      inputOXPHOSgenes=tolower(x[,slCol])
    } else {
      inputOXPHOSgenes=""
    }
    
  } else {
    if(inputGeneName=="ensembl_gene_id"){
      inputOXPHOSgenes=tolower(x$ensembl_gene_id)
    } else {
      if(inputGeneName=="entrezgene_id"){
        x=.extraMitochondrialGenes()
        inputOXPHOSgenes=tolower(x$entrezgene_id)
      }
    }
  }
  
  
  if(!standardize_gene_names){
    inputOXPHOSgenes=unique(tolower(inputGeneset))
  }
  inputOXPHOSgenes=inputOXPHOSgenes[inputOXPHOSgenes!=""]
  if(is.null(inputGeneset)){
    print(paste("Number of OXPHOS genes in the dataset:",length(which(rwNames %in% inputOXPHOSgenes)),"/",length(inputOXPHOSgenes)))
  } else {
    print(paste("Number of genes in the dataset:",length(which(rwNames %in% inputOXPHOSgenes)),"/",length(inputOXPHOSgenes)))
  }
  
  if(length(inputOXPHOSgenes)>0){
    if(sum(colnames(colData(inputData))=="QC_Gene_total_count")==0|recalculate_nUMI){
      tmpColSums=c()
      if(ncol(inputData)>10000){
        for(i in seq(1,ncol(inputData),10000)){
          tmpColSums=c(tmpColSums,apply(counts(inputData)[,i:min(i+10000-1,ncol(inputData))],2,sum))
        }
      } else {
        tmpColSums=apply(counts(inputData),2,sum)
      }
    } else {
      tmpColSums=inputData$QC_Gene_total_count
    }
    
    
    
    res=colSums(x = counts(inputData)[which(rwNames %in% inputOXPHOSgenes), , drop = FALSE])/tmpColSums
    res=res*100
  } else {
    res=rep(NA,ncol(inputData))
  }
  
  return(res)
}

.extrarRNAPctFn=function(inputData,organism,inputrRNAgenes=NULL,inputGeneName=NULL,recalculate_nUMI=T){
  
  rwNames=tolower(row.names(inputData))
  x=.extrarRNAGenes(organism = organism)
  if(sum(grepl("\\.",rwNames))>0){
    rwNames=strsplit(rwNames,"\\.")
    rwNames=unlist(lapply(rwNames,function(x)x[1]))
  }
  if(is.null(inputGeneName)){
    tmpCols=c("ensembl_gene_id","entrezgene_id",'gene_name')
    
    slCounts=0
    slCol=""
    
    for(i in tmpCols){
      tmp=sum(rwNames %in% tolower(x[,i]))
      if(tmp>slCounts){
        slCounts=tmp
        slCol=i
      }
    }
    if(slCol!=""){
      inputrRNAgenes=tolower(x[,slCol])
    } else {
      inputrRNAgenes=""
    }
    
  } else {
    if(inputGeneName=="ensembl_gene_id"){
      inputrRNAgenes=tolower(x$ensembl_gene_id)
    } else {
      if(inputGeneName=="entrezgene_id"){
        x=.extrarRNAGenes()
        inputrRNAgenes=tolower(x$entrezgene_id)
      }
    }
  }
  
  print(paste("Number of rRNA genes in the dataset:",length(which(rwNames %in% inputrRNAgenes)),"/",sum(x$gene_biotype=="rRNA")))
  
  inputrRNAgenes=inputrRNAgenes[inputrRNAgenes!=""]
  if(length(inputrRNAgenes)>0){
    tmpColSums=c()
    if(sum(colnames(colData(inputData))=='QC_Gene_total_count')==0|recalculate_nUMI){
      if(ncol(inputData)>10000){
        for(i in seq(1,ncol(inputData),10000)){
          tmpColSums=c(tmpColSums,apply(counts(inputData)[,i:min(i+10000-1,ncol(inputData))],2,sum))
        }
      } else {
        tmpColSums=apply(counts(inputData),2,sum)
      }
    } else {
      tmpColSums=inputData$QC_Gene_total_count
    }
    
    
    
    res=colSums(x = counts(inputData)[which(toupper(rwNames) %in% toupper(inputrRNAgenes)), , drop = FALSE])/tmpColSums
    res=res*100
  } else {
    res=rep(NA,ncol(inputData))
  }
  
  return(res)
}

# for every cell, calculate the pct rRNA protein
.extraRPctFn=function(inputData,organism,inputRPgenes=NULL,inputGeneName=NULL,recalculate_nUMI=T){
  
  # Convert row names of inputData to lowercase
  rwNames=tolower(row.names(inputData))

  # Get ribosomal protein genes if not provided
  x=.extraRibosomalProteinGenes(organism = organism)
  
  # If row names contain dot, remove the dot and everything after it
  if(sum(grepl("\\.",rwNames))>0){
    rwNames=strsplit(rwNames,"\\.")
    rwNames=unlist(lapply(rwNames,function(x)x[1]))
  }
  
  # Define inputRPgenes, the list of ribosomal protein genes
  # If inputGeneName is not provided, try to find the column name of the gene name
  if(is.null(inputGeneName)){
    tmpCols=c("ensembl_gene_id","entrezgene_id",'gene_name')
    
    slCounts=0
    slCol=""
    
    # Loop through the columns and find the column with the most matches
    for(i in tmpCols){
      tmp=sum(rwNames %in% tolower(x[,i]))
      if(tmp>slCounts){
        slCounts=tmp
        slCol=i
      }
    }
    if(slCol!=""){
      inputRPgenes=tolower(x[,slCol])
    } else {
      inputRPgenes=""
    }
    
  } 
  # If inputGeneName is provided, use the column name to get the gene names
  else {
    if(inputGeneName=="ensembl_gene_id"){
      inputRPgenes=tolower(x$ensembl_gene_id)
    } else {
      if(inputGeneName=="entrezgene_id"){
        x=.extraRibosomalProteinGenes()
        inputRPgenes=tolower(x$entrezgene_id)
      }
    }
  }
  
  print(paste("Number of ribosomal protein genes in the dataset:",length(which(rwNames %in% inputRPgenes)),"/",nrow(x)))
  
  inputRPgenes=inputRPgenes[inputRPgenes!=""]

  
  if(length(inputRPgenes)>0){
    tmpColSums=c()
    if(sum(colnames(colData(inputData))=='QC_Gene_total_count')==0|recalculate_nUMI){
      if(ncol(inputData)>10000){
        for(i in seq(1,ncol(inputData),10000)){
          tmpColSums=c(tmpColSums,apply(counts(inputData)[,i:min(i+10000-1,ncol(inputData))],2,sum))
        }
      } else {
        tmpColSums=apply(counts(inputData),2,sum)
      }
    } else {
      tmpColSums=inputData$QC_Gene_total_count
    }
    
    res=colSums(x = counts(inputData)[which(toupper(rwNames) %in% toupper(inputRPgenes)), , drop = FALSE])/tmpColSums
    res=res*100
  } else {
    res=rep(NA,ncol(inputData))
  }
  
  return(res)
}
