
.mycBindFn=function(inputList,batchNames=NULL,verbose_level=1){
  # Todo: rename this function to something more descriptive
  # binds multiple singleCellExpression datasets columnwise with differring number of rows.
  # inputList: the list of datasets to be merged
  # batchNames: list, nullable. the batch name to be assigned to each dataset. length(batchNames)==length(inputList)
  res_m=""
  if(!is.null(batchNames)){
    if(length(inputList)==1){
      res_m=inputList[[1]]
    } else if(length(inputList)==2){
      res_m=.extracBindDetailFn(x1=inputList[[1]],x2=inputList[[2]],batchNames = batchNames[1:2])
    } else if(length(inputList)>2){
      res_m=.extracBindDetailFn(x1=inputList[[1]],x2=inputList[[2]],batchNames = batchNames[1:2])
      for(i in 3:length(inputList)){
        res_m=.extracBindDetailFn(x1=res_m,x2=inputList[[i]],batchNames=c("",batchNames[i]))
      }
    }
  } 
  else {
    if(length(inputList)==1){
      res_m=inputList[[1]]
    } else if(length(inputList)==2){
      res_m=.extracBindDetailFn(x1=inputList[[1]],x2=inputList[[2]],batchNames = c("",""))
      # otherwise iteratively bind the datasets by running extracBindDetailFn in a loop
      # TODO: investigate whether this is the most efficient way to go about this.....
    } else if(length(inputList)>2){
      res_m=.extracBindDetailFn(x1=inputList[[1]],x2=inputList[[2]],batchNames = c("",""))
      for(i in 3:length(inputList)){
        
        res_m=.extracBindDetailFn(x1=res_m,x2=inputList[[i]],batchNames=c("",""))
        if(verbose_level>1){
          print(paste("dataset:",i,"; nrow:",nrow(res_m),"; ncol",ncol(res_m)))
        }
      }
    }
  }
  
  if(verbose_level>=1){
    print("batch information is in the anno_batch variable")
  }
  
  return(res_m)
  
}

.myOneHotFn=function (inputVector) {
  if(length(unique(inputVector))>1){
    if(T){
      #using caret package
      require(caret)
      if(sum(is.na(inputVector))>0){
        inputVector[is.na(inputVector)]="NA"
      }
      formula="~."
      data=data.frame(data=inputVector,stringsAsFactors = F)
      sep = "."
      levelsOnly = FALSE
      fullRank = FALSE
      
      formula <- as.formula(formula)
      if (!is.data.frame(data)) 
        data <- as.data.frame(data, stringsAsFactors = FALSE)
      vars <- all.vars(formula)
      if (any(vars == ".")) {
        vars <- vars[vars != "."]
        vars <- unique(c(vars, colnames(data)))
      }
      isFac <- unlist(lapply(data[, vars, drop = FALSE], is.factor))
      if (sum(isFac) > 0) {
        facVars <- vars[isFac]
        lvls <- lapply(data[, facVars, drop = FALSE], levels)
        if (levelsOnly) {
          tabs <- table(unlist(lvls))
          if (any(tabs > 1)) {
            stop(paste("You requested `levelsOnly = TRUE` but", 
                       "the following levels are not unique", "across predictors:", 
                       paste(names(tabs)[tabs > 1], collapse = ", ")))
          }
        }
      } else {
        facVars <- NULL
        lvls <- NULL
      }
      trms <- attr(model.frame(formula, data), "terms")
      out <- list(call = match.call(), form = formula, vars = vars, 
                  facVars = facVars, lvls = lvls, sep = sep, terms = trms, 
                  levelsOnly = levelsOnly, fullRank = fullRank)
      class(out) <- "dummyVars"
      trsf <- data.frame(predict(out, newdata = data))
      
      tmp=gsub("^data","",colnames(trsf))
      tmp=gsub("^\\.","",tmp)
      colnames(trsf)=tmp
    } else {
      if(length(names(inputVector))==0){
        trsf=reshape2::dcast(data=data.frame(var=1:length(inputVector),outcome=as.character(inputVector)),var ~ outcome, length)
        row.names(trsf)=trsf[,1]
        trsf=trsf[,-1]
      } else {
        inputVector=inputVector[order(inputVector)]
        inputVector=names(inputVector)
        trsf=reshape2::dcast(data=data.frame(var=1:length(inputVector),outcome=as.character(inputVector)),var ~ outcome, length)
        trsf=as.matrix(trsf)
        row.names(trsf)=trsf[,1]
        trsf=trsf[,-1]
      }
      
    }
  } else {
    trsf=matrix(1,ncol=1,nrow=length(inputVector))
    row.names(trsf)=names(inputVector)
    colnames(trsf)=unique(inputVector)
  }
  
  
  
  return(trsf)
}

.extracBindDetailFn=function(x1,x2,batchNames){
    # x1: SingleCellExperiment or Seurat object
    # x2: SingleCellExperiment or Seurat object
    # batchNames: column of strings of length 2, 
    # containing the batch names for x1 and x2 respectively
    # there is no requirement of the same rows or columns in x1 and x2

    # outputs: 
    # x_m: SingleCellExperiment object with x1 and x2 bound columnwise
    # and with the number of rows equal to the union of the rows in x1 and x2
    # the colData of this matrix will have an $anno_batch column.
    # that uses the `batchNames` argument's elements to indicate the batch of each cell
    

  require(Matrix)
  
  # tmp1 stores the row names present in x1, and not in x2
  tmp1=setdiff(row.names(x1),row.names(x2))

  # get fd_, pd_, and expr_ for x1 and x2
  if(class(x1)=="SingleCellExperiment"){
    fd1=rowData(x1)
    expr1=counts(x1)
    pd1=as.data.frame(colData(x1))
  } else if(class(x1)=="Seurat") {
    fd1=as.data.frame(x1@assays$RNA@meta.features)
    expr1=x1@assays$RNA@counts
    pd1=as.data.frame(x1@meta.data)
  } else {
    stop("Unrecognized dataset class!")
  }
  
  if(class(x2)=="SingleCellExperiment"){
    fd2=rowData(x2)
    expr2=counts(x2)
    pd2=as.data.frame(colData(x2))
  } else if(class(x2)=="Seurat") {
    fd2=as.data.frame(x2@assays$RNA@meta.features)
    expr2=x2@assays$RNA@counts
    pd2=as.data.frame(x2@meta.data)
  } else {
    stop("Unrecognized dataset class!")
  }
  
  # get the feature data for the rows in x1 that are NOT in x2
  tmpAnno1=fd1[which(row.names(x1) %in% tmp1),]
  
  # tmpMat1 is a sparse matrix with the same number of rows present in x1, and not in x2
  # and the same number of columns as x2
  tmpMat1=sparseMatrix(i=NULL,j=NULL,dims = c(length(tmp1),ncol(x2)))
  row.names(tmpMat1)=tmp1
  
  # since they have the same number of cols, you can row bind expr2 and tmpMat1
  # now x2_c has the same number of rows as union(x1, 2)
  x2_c=rbind(expr2,tmpMat1)
  x2_c_pd=as.data.frame(pd2)
  # plyr::rbind.fill rbinds a list of data frames filling missing columns with NA.
  # the missing columns will be genes that are only in x1
  x2_c_fd=as.data.frame(plyr::rbind.fill(as.data.frame(fd2),as.data.frame(tmpAnno1)))
  x2_c=SingleCellExperiment(assays = list(counts = x2_c),colData = x2_c_pd,rowData=x2_c_fd)
  
  # tmp2 is the rows in x2 that are not in x1
  # tmpAnno2 is the feature data for the rows in x2 that are not in x1
  tmp2=setdiff(row.names(x2),row.names(x1))
  tmpAnno2=fd2[which(row.names(x2) %in% tmp2),]
  # tmpMat2 is a sparse matrix with the same number of rows as genese x2, and not in x1
  # and same number of cols as x1
  tmpMat2=sparseMatrix(i=NULL,j=NULL,dims = c(length(tmp2),ncol(x1)))
  row.names(tmpMat2)=tmp2
  
  # you can now row bind expr1 and tmpMat2
  # now x1_c has the same number of rows as as union(x1, 2)
  x1_c=rbind(expr1,tmpMat2)
  
  x1_c_pd=as.data.frame(pd1)
  x1_c_fd=plyr::rbind.fill(as.data.frame(fd1),as.data.frame(tmpAnno2))
  x1_c=SingleCellExperiment(assays = list(counts = x1_c),colData = x1_c_pd,rowData=x1_c_fd)
  
  x2_c=x2_c[match(row.names(x1_c),row.names(x2_c)),]
  
  # after matching, all row names should align
  if(all(row.names(x2_c)==row.names(x1_c))){
    
    # if the batches have names, prepend them to the column names to make an `anno_batch` column
    if(batchNames[1]!=""){
      row.names(colData(x1_c))=paste0(batchNames[1],"_",colnames(x1_c))
      colnames(x1_c)=paste0(batchNames[1],"_",colnames(x1_c))
      x1_c$anno_batch=batchNames[1]
    }
    
    if(batchNames[2]!=""){
      colnames(x2_c)=paste0(batchNames[2],"_",colnames(x2_c))
      row.names(colData(x2_c))=paste0(batchNames[2],"_",colnames(x2_c))
      x2_c$anno_batch=batchNames[2]
    }
    
    # now just cbind
    x_m_exp=cbind(counts(x1_c),counts(x2_c))
    pd_m_exp=plyr::rbind.fill(as.data.frame(colData(x1_c)),as.data.frame(colData(x2_c)))
    fd=as.data.frame(rowData(x1_c))
    x_m=SingleCellExperiment(assays=list(counts=x_m_exp),colData=pd_m_exp,rowData=fd)
  } 
  else {
    stop("Error in the merging!")
  }
  return(x_m)
}

.extraDoubletMakerFn=function(
    inputData,
    label_col,
    sel_labels,
    cells_per_group=30){
# inputData: 
# label_col: str, the column name in inputData 

  print("Adding the doublet data...")
  selCells=list()
  for(i in unique(sel_labels)){
    tmp=inputData[,colData(inputData)[,label_col]==i]
    if(ncol(tmp)>0){
      tmp=tmp[,sample(ncol(tmp),min(ncol(tmp),cells_per_group))]
      selCells=c(selCells,list(new=tmp))
      names(selCells)[length(selCells)]=i
    }
  }
  
  doublet_exp=list()
  lbls=c()
  for(i in 1:(length(selCells)-1)){
    ds1=selCells[[i]]
    for(j in (i+1):length(selCells)){
      ds2=selCells[[j]]
      for(ik in 1:ncol(ds1)){
        for(jk in 1:ncol(ds2)){
          x=cbind(counts(ds1)[,ik],counts(ds2)[,jk])
          x=rowSums(x)
          x=data.frame(new=x)
          colnames(x)=paste0(names(selCells)[i],"_",names(selCells)[j])
          doublet_exp=c(doublet_exp,list(x))
          lbls=c(lbls,paste0(names(selCells)[i],"_",names(selCells)[j]))
        }
      }
    }
  }
  
  doublet_exp2=do.call("cbind",doublet_exp)
  
  colnames(doublet_exp2)=paste0(colnames(doublet_exp2),"_",LETTERS[sample(26,ncol(doublet_exp2),replace = T)])
  while(sum(duplicated(colnames(doublet_exp2)))>0){
    colnames(doublet_exp2)[duplicated(colnames(doublet_exp2))]=paste0(colnames(doublet_exp2)[duplicated(colnames(doublet_exp2))],LETTERS[sample(26,sum(duplicated(colnames(doublet_exp2))),replace = T)])
  }
  
  pd=as.data.frame(colData(inputData))
  for(i in 1:ncol(pd)){
    if(class(pd[,i])=="character"){
      pd[,i]="synthetic"
    } else{
      pd[,i]=NA
    }
  }
  
  pd=pd[1:ncol(doublet_exp2),]
  
  pd[,label_col]=lbls
  pd$is_doublet="YES"
  inputData$is_doublet="NO"
  
  doublet_exp=Matrix(as.matrix(doublet_exp2),sparse=T)
  
  fd=as.data.frame(rowData(inputData))
  if(!all(row.names(inputData)==row.names(doublet_exp))){
    stop("Error in making the doublet data!")
  }
  
  expInput=counts(inputData)
  expM=cbind(expInput,doublet_exp)
  pdM=rbind(as.data.frame(colData(inputData)),pd)
  
  res=SingleCellExperiment(
    assays=list(counts=expM),
    colData=pdM,
    rowData=fd)
  
  print("Doublet data is added.")
  return(res)
}

.extraExport2SeuratFn=function(
    inputData){
    #inputData: SingleCellExperiment

    # Outputs a Seurat object whose metadata is the pData (col data) of the input object
    # If the input object has feature data (row data), that is also added 
    # to the Seurat object's @assays$RNA@meta.features

  require(scran)
  nbt=Seurat::CreateSeuratObject(counts=counts(inputData),
                                 project = "SeuratProject",
                                 assay = "RNA",
                                 min.cells = 0,
                                 min.features = 0,
                                 names.field = 1,
                                 names.delim = "-",
                                 meta.data = .pData(inputData))
  

  if(ncol(.fData(inputData = inputData))>0){
    nbt@assays$RNA@meta.features=.fData(inputData=inputData)
  }
  
  return(nbt)
}

.mySplitObject=function(object,colName){
    # object:  the object to be split. Either a SingleCellExperiment or a Seurat object
    # colName: str, the column name in the colData / @meta.data of the object to be used for splitting

    # Outputs `res` as a list of objects, each of which is a subset of the input object
    # of the same TYPE as the input object, split columnwise. 
  res=list()
  if(class(object)=="Seurat"){
    res=Seurat::SplitObject(object,split.by = colName)
  } 
  else { 
    
    # TODO: test if the following code works and is performant 
    # res_list = split(object, object$colData[, colName])
    # names(res_list) = unique(object$colData[, colName])

    # pd is Phenotypic Data, i.e. the colData of the object
    pd=colData(object)[,colName]
    # 
    for(i in unique(pd)){
        res=c(res,list(object[,which(pd==i)]))
        names(res)[length(res)]=i
    }
    
    
  }
  return(res)
}

.pData=function(inputData){
    # takes in a singleCellExperiment or a Seurat object
    # if input arg is a singleCellExperiment, returns the colData
    # (i.e. phenotypic data, cell data) as a dataframe
    # if a seurat object is passed in, it seems to return the @meta.data, which is the colData
  pd=""
  if(class(inputData)=="Seurat"){
    pd=res=inputData[[]]
  } else {
    if(class(inputData)=="SingleCellExperiment"){
      pd=as.data.frame(colData(inputData))
    } else {
      pd=pData(inputData)
    }
    
  }
  
  return(pd)
}

.fData=function(inputData){
    # takes in a singleCellExperiment and returns the row data 
    # (i.e. featuredata, gene data) as a dataframe
  fd=as.data.frame(rowData(inputData))
  return(fd)
}

downsampleListAfterStart <- function(a_list, start_ind, step){
    # return every value of a_list from 1:start_ind, then every step-th value after that
    # a_list: list
    # stard_ind: int, the index after which to start downsampling
    # step: int, the step size for downsampling
    return(
        c(a_list[1:start_ind], 
        a_list[(start_ind+1):length(a_list)][seq(1, length(a_list)-start_ind, step)]))
}

getFracAssignableDemuxlet=function(df){
    # df: dataframe with columns `BEST` and `PRB.SNG1` 
    return(sum((df$PRB.SNG1 >= 0.9) & (substr(df$BEST, 1, 3) == 'SNG')) / nrow(df))
}

getFracAssignableVireo=function(df){

    nrows=nrow(df)
    df = df[!df$donor_id %in% c('unassigned', 'doublet'), ]
    return(sum(df$prob_max >= 0.9) / nrows)
}

orderDFByColRank=function(df, col, asc=FALSE, log_y_axis=FALSE){
    # rank order by `col`
    if(asc){
        df$rank = rank(df[[col]], ties.method = "first")     
    } else {
        df$rank = rank(-df[[col]], ties.method = "first") 
    }
    df=df[order(df$rank), ]
    if(log_y_axis){ 
        df[[col]] = log10(df[[col]]+1)
    }

    return(df)
}
