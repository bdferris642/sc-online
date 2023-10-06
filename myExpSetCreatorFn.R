.myExpSetCreatorFn=function(inputExpData,
                            organism,minExpCells=0,
                            inputPdata=NULL,
                            inputFdata=NULL,
                            addExtraAnno=T,
                            server=T,
                            redownload_files=T,
                            ncores=5){
  # inputExpData: Ngene x Ncell dgCMatrix of RNA counts
  # organism: str, "human", "macaque", or "mouse"
  # minExpCells: int, minimum number of cells that a gene must be expressed in to be included in the output

  require(scran)
  require(DropletUtils)
  require(Matrix)
  if(is.null(colnames(inputExpData))){
    colnames(inputExpData)=1:ncol(inputExpData)
  }
  
  if(tolower(class(inputExpData))=="seurat"){
    inputExpData=.sconline.convertSeuratToExpressionSet(object=inputExpData)
  }
  
  if(is(inputExpData,"SingleCellExperiment")){
    tmp=as.data.frame(rowData(inputExpData))
    if(is.null(inputFdata)){
      inputFdata=tmp
    } else {
      inputFdata=cbind(inputFdata,tmp)
    }
    
    tmpPdata=as.data.frame(colData(inputExpData))
    if(is.null(inputPdata)){
      inputPdata=tmpPdata
    } else {
      inputPdata=cbind(inputPdata,tmpPdata)
    }
  }
  
  if(is(inputExpData,"Seurat")){
    tmpPdata=as.data.frame(inputPdata[[]])
    if(is.null(inputPdata)){
      inputPdata=tmpPdata
    } else {
      inputPdata=cbind(inputPdata,tmpPdata)
    }
  }
  
  
  if(!is.null(inputPdata)){
    pd <- inputPdata
  } else {
    pd=data.frame(sample=colnames(inputExpData),stringsAsFactors = F)
    row.names(pd)=colnames(inputExpData)
  }
  
  if(tolower(organism)=="human"){
    tmpName=.extraHumanGeneAnnoAdderFn(row.names(inputExpData),server = server)
  } else if (tolower(organism)=="mouse"){
    tmpName=.extraMouseGeneAnnoAdderFn(row.names(inputExpData),server = server)
  } else if(tolower(organism)=="macaque"){
    tmpName=.extraMacaqueGeneAnnoAdderFn(row.names(inputExpData),server = server)
  } else {
    warning("unknown organism!")
  }
  
  fd=NULL
  if(tolower(organism) %in% c("human","mouse","macaque")){
    if(all(row.names(inputExpData)==row.names(tmpName))){
    fd=tmpName
    
    if(!is.null(inputFdata)){
      fd=cbind(inputFdata,fd)
    }
    if(sum(colnames(fd)=="ID")==0){
      fd$ID=row.names(inputExpData)
    }
    row.names(fd)=row.names(inputExpData)
    
    if(tolower(organism) %in% c("mouse","human","macaque")){
      mycc.genes=.extraCellCycleGenes(organism = organism)
      fd$QC_ssGenes="No"
      fd$QC_g2mGenes="No"
      fd$QC_ssGenes[fd$ensembl_gene_id %in% row.names(mycc.genes$cc.s)]="s_phase"
      fd$QC_g2mGenes[fd$ensembl_gene_id %in% row.names(mycc.genes$cc.g2m)]="g2m_phase"
      
      mymito.genes=.extraMitoGenes(organism=organism)
      fd$QC_mtGenes="No"
      fd$QC_mtGenes[fd$ensembl_gene_id %in% mymito.genes$ensembl_gene_id]="Yes"
      
      myIEG.genes=.extraIEGGenes(organism=organism,server = server)
      fd$QC_IEG_Genes="No"
      fd$QC_IEG_Genes[fd$ensembl_gene_id %in% myIEG.genes$ensembl_gene_id]="Yes"
    }
   
    
  } else {
      print("Error!")
    }
    if(sum(colnames(fd) %in% c("seqnames", "ranges", "strand", "start", "end", "width","element"))>0){
      colnames(fd)[colnames(fd) %in% c("seqnames", "ranges", "strand", "start", "end", "width","element")]=paste0("anno_",colnames(fd)[colnames(fd) %in% c("seqnames", "ranges", "strand", "start", "end", "width","element")])
      
    }
  } else {
    fd=inputFdata
  }
  
  
  
  if(class(inputExpData)=="SingleCellExperiment"){
    res=counts(inputExpData)
  } else {
    if(class(inputExpData)!="dgCMatrix"){
      if(class(inputExpData)=="data.frame"){
        res=Seurat::as.sparse(inputExpData)
      } else {
        res=.matrixExtraction(inputExpData)
      }
      row.names(res)=row.names(inputExpData)
      colnames(res)=colnames(inputExpData)
    } else {
      res=inputExpData
    }
    
  }
  
  colnames(res)=colnames(inputExpData)
  
  res <- as(res, "dgCMatrix")
  
  res=SingleCellExperiment(assays = list(counts = res),colData = pd,rowData=fd)
  rm(inputExpData)
  
  if(minExpCells>0){
    tmpFilter=c()
    
    if(nrow(res)>10000){
      for(i in seq(1,nrow(res),10000)){
        tmpFilter=c(tmpFilter,apply(assays(res)[["counts"]][i:min(i+10000-1,nrow(res)),],1,function(x) sum(x>0)))
      }
    } else {
      # calculates a filter for rows in the "counts" matrix, 
      # where rows with more nonzero values will have higher values in the tmpFilter vector.
      tmpFilter=apply(assays(res)[["counts"]],1,function(x) sum(x>0))
    }
    
    
    res=res[which(tmpFilter>=minExpCells),]
  }
  
  
  
  if(addExtraAnno){
    #inputExpSet = res;organism=organism;server=server;redownload_files=redownload_files
    res=.extraQCAnnoAdderFn(inputExpSet = res,organism=organism,server=server,redownload_files=redownload_files,ncores=ncores)
  }
  
  for(i in 1:ncol(colData(res))){
    if(class(colData(res)[,i])=="factor"){
      colData(res)[,i]=droplevels(colData(res)[,i])
    }
  }
  
  return(res)
}
