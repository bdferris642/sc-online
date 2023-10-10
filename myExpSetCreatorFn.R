# TODO: place in single Consts file
ANNOTATABLE_ORGANISMS = c("human", "macaque", "mouse")
ANNO_COLS_TO_PREFIX = c(
  "seqnames", "ranges", "strand", "start", "end", "width","element")

.myExpSetCreatorFn=function(inputExpData,
                            organism,
                            minExpCells=0,
                            inputPdata=NULL,
                            inputFdata=NULL,
                            addExtraAnno=T,
                            ncores=5){
  # inputExpData: Ngene x Ncell dgCMatrix of RNA counts, e.g. the @assays$RNA@counts slot of a Seurat object
  # organism: str, "human", "macaque", or "mouse"
  # minExpCells: int, minimum number of cells that a gene must be expressed in to be included in the output
  # inputPdata: dataframe or NULL
  #   if provided, df from @meta.data slot of a Seurat object with phenotype information for each cell. Must have same row dimension as inputExpData column dimension.
  #   if NULL, created from the coldata of the inputExptData
  # inputFdata: dataframe or NULL
  #   if provided, df with gene annotation information for each gene. Must have same row dimension as inputExpData row dimension.
  #   if NULL, created from the rowdata of the inputExptData
  # addExtraAnno: bool, whether to add extra columns to the output
  # ncores: int, number of cores to use

  # OUTPUT
  # SingleCellExperiment object with the following slots:
  #   where colData is specified by the input pData (if provided)
  #   and rowData is specified by the input fData (if provided)
  
  # Function does as follows:
  # 1. generate temporary fdata and pdata if it is provided (not sure this is necessary)
  # 2. extract annotations for the genome (if the organism makes this possible)
  # 3. set feature data 
  # 	a. add an `ID` column, if it doesn't already exist, based on the row names (genes) of the input data
  # 	b. add certain annotations (QC_ssGenes, QC_g2mGenes, QC_mtGenes, and QC_IEG_Genes) to the eventual row data
  # 	c. for certain columns in the feature data / future row data that come from the gene annotation, 
  # 		add an "anno_" prefix to the column name
  # 4. Coerce into a dgCMatrix of counts, and cast as a SingleCellExperiment type, with PData colData, and FData rowData
  # 5. Remove genes (rows) that are not expressed in at least `minExpCells` cells (columns)
  # 6. Optionally, additional QC columns are added to the rowData of the SingleCellExperiment
  # 7. Unused factors are removed

  require(scran)
  require(DropletUtils)
  require(Matrix)

  # if no colnames are provided, number the columns
  if(is.null(colnames(inputExpData))){
    colnames(inputExpData)=1:ncol(inputExpData)
  }
  
  # TODO: remove this and standardize input class to this function
  if(tolower(class(inputExpData))=="seurat"){
    inputExpData=.sconline.convertSeuratToExpressionSet(object=inputExpData)
  }
  
  # TODO: remove this and standardize input class to this function
  # if the input is already a SingleCellExperiment object, can use the rowData and colData to generate FData and PData for the output
  # TODO: rename tmp to something more descriptive
  if(is(inputExpData,"SingleCellExperiment")){
    tmp=as.data.frame(rowData(inputExpData))
    if(is.null(inputFdata)){
      inputFdata=tmp
    } 
    else {
      inputFdata=cbind(inputFdata,tmp)
    }
    
    tmpPdata=as.data.frame(colData(inputExpData))
    if(is.null(inputPdata)){
      inputPdata=tmpPdata
    } 
  else {
      inputPdata=cbind(inputPdata,tmpPdata)
    }
  }
  
  # TODO: remove this and standardize input class to this function
  if(is(inputExpData,"Seurat")){
    tmpPdata=as.data.frame(inputPdata[[]])
    if(is.null(inputPdata)){
      inputPdata=tmpPdata
    } else {
      inputPdata=cbind(inputPdata,tmpPdata)
    }
  }
  
  # In the base case, if no inputPdata is provided, create it from the colnames of the inputExpData (cell names)
  if(!is.null(inputPdata)){
    pd <- inputPdata
  } else {
    pd=data.frame(sample=colnames(inputExpData),stringsAsFactors = F)
    row.names(pd)=colnames(inputExpData)
  }
  
  # TODO: DRY .extra<Organism>GeneAnnoAdderFn()s into one
  if(tolower(organism)=="human"){
    tmpName=.extraHumanGeneAnnoAdderFn(row.names(inputExpData))
  } else if (tolower(organism)=="mouse"){
    tmpName=.extraMouseGeneAnnoAdderFn(row.names(inputExpData))
  } else if(tolower(organism)=="macaque"){
    tmpName=.extraMacaqueGeneAnnoAdderFn(row.names(inputExpData))
  } else {
    warning("unknown organism!")
  }
  
  # set feature data
  fd=NULL
  if(tolower(organism) %in% ANNOTATABLE_ORGANISMS){ 
    # if the row names (genes) of the input data match those of the genes, then use the gene annotation
    if(all(row.names(inputExpData)==row.names(tmpName))){
      fd=tmpName
    
    # if the inputFdata is provided, use that
    if(!is.null(inputFdata)){
      fd=cbind(inputFdata,fd)
    }

    # if the fd has no `ID` column, add one from the row names (genes) of the inputExpData
    if(sum(colnames(fd)=="ID")==0){
      fd$ID=row.names(inputExpData)
    }
    row.names(fd)=row.names(inputExpData)
    
    
    # TODO: Binarize these columns
    # Certain gene annotations are added as columns (QC_ssGenes, QC_g2mGenes, QC_mtGenes, and QC_IEG_Genes) to the eventual row data
    # These data take the form of {"No", "<Some-Other-String"}
    mycc.genes=.extraCellCycleGenes(organism = organism)
    fd$QC_ssGenes="No"
    fd$QC_g2mGenes="No"
    fd$QC_ssGenes[fd$ensembl_gene_id %in% row.names(mycc.genes$cc.s)]="s_phase"
    fd$QC_g2mGenes[fd$ensembl_gene_id %in% row.names(mycc.genes$cc.g2m)]="g2m_phase"
    
    mymito.genes=.extraMitoGenes(organism=organism)
    fd$QC_mtGenes="No"
    fd$QC_mtGenes[fd$ensembl_gene_id %in% mymito.genes$ensembl_gene_id]="Yes"
    
    myIEG.genes=.extraIEGGenes(organism=organism)
    fd$QC_IEG_Genes="No"
    fd$QC_IEG_Genes[fd$ensembl_gene_id %in% myIEG.genes$ensembl_gene_id]="Yes"
    
  # Alert the user if the genome is not annotatable 
  } else {
      warning(paste(organism, "is not in", list(ANNOTATABLE_ORGANISMS), sep=" "))
    }
  

  # for certain columns in the feature data / future row data that come from the gene annotation, 
  # add an "anno_" prefix to the column name
  if(sum(colnames(fd) %in% ANNO_COLS_TO_PREFIX) >0){
      colnames(fd)[
        colnames(fd) %in% ANNO_COLS_TO_PREFIX
        ]=paste0("anno_",colnames(fd)[colnames(fd) %in% ANNO_COLS_TO_PREFIX])
    }
  } 
  else {
    fd=inputFdata
  }
  
  # TODO: remove this and standardize input class to this function 
  # generate res, a dgCMatrix
  if(class(inputExpData)=="SingleCellExperiment"){
    res=counts(inputExpData)
  } 
  else {
    if(class(inputExpData)!="dgCMatrix"){
      if(class(inputExpData)=="data.frame"){
        res=Seurat::as.sparse(inputExpData)
      } 
      else {
        res=.matrixExtraction(inputExpData)
      }
      row.names(res)=row.names(inputExpData)
      colnames(res)=colnames(inputExpData)
    } 
    else {
      res=inputExpData
    }
  }
  
  colnames(res)=colnames(inputExpData)
  res <- as(res, "dgCMatrix")

  # phenotypic data is set as the colData, and feature data is set as the rowData
  res=SingleCellExperiment(assays=list(counts=res), colData=pd, rowData=fd)
  rm(inputExpData)
  
  # calculate a filter for rows in the "counts" matrix, 
  # where rows with more nonzero values will have higher values in the tmpFilter vector.
  # If the dataset is large, do in chunks
  if(minExpCells>0){
    tmpFilter=c()
    
    if(nrow(res)>10000){
      for(i in seq(1,nrow(res),10000)){
        tmpFilter=c(tmpFilter,
                    apply(assays(res)[["counts"]][i:min(i+10000-1,nrow(res)),],
                    1,
                    function(x) sum(x>0))
                  )
      }
    } else {
      
      tmpFilter=apply(assays(res)[["counts"]],1,function(x) sum(x>0))
    }
    
    
    res=res[which(tmpFilter>=minExpCells),]
  }
 
  if(addExtraAnno){  # Adds QC columns to the output
    res=.extraQCAnnoAdderFn(inputExpSet = res,organism=organism,ncores=ncores)
  }
  
  # remove unused factors from the colData
  for(i in 1:ncol(colData(res))){
    if(class(colData(res)[,i])=="factor"){
      colData(res)[,i]=droplevels(colData(res)[,i])
    }
  }
  
  return(res)
}