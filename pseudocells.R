source("~/sc-online/utils.R")
source("~/sc-online/extraGenes.R")


.extra_sconline.FixedSizeFn=function(
  inputExp,
  inputEmbedding=NULL,
  pseudocell_size=40,
  n.neighbors=20,
  n.trees=50,
  nPCs=30,
  k.param = 20,
  include.outlier.score=F,
  rand_pseudobulk_mod=T,
  cols_to_sum=NULL){
  #inputExp=datalist[[i]];inputEmbedding=NULL;pseudocell_size=50;nPCs=30
  #n.neighbors=20;n.trees=50;k.param = 20
  
  if(ncol(inputExp)==1){
    inputExp$pseudocell_size=1
    return(inputExp)
  }
  
  pseudocell_count=round(ncol(inputExp)/pseudocell_size)
  if(pseudocell_count<=1){
    expData=matrix(as.numeric(rowSums(counts(inputExp))),ncol=1)
    row.names(expData)=row.names(inputExp)
    colnames(expData)=colnames(inputExp)[1]
    
    pd=as.data.frame(colData(inputExp))
    pd=pd[1,]
    pd$pseudocell_outlier_score=NA
    pd$pseudocell_size=ncol(inputExp)
    
    if(!is.null(cols_to_sum)){
      for(icoltosum in cols_to_sum){
        pd[1,icoltosum]=sum(as.data.frame(colData(inputExp))[,icoltosum])
      }
    }
    
    fd=as.data.frame(rowData(inputExp))
    
    res=SingleCellExperiment(assays = list(counts = expData),colData = pd,rowData=fd)
    return(res)
  }
  
  if(rand_pseudobulk_mod){
    prop_m_hardCluster=sample(colnames(inputExp))
    prop_m_hardCluster=split(prop_m_hardCluster,floor(ecdf(seq_along(prop_m_hardCluster))(seq_along(prop_m_hardCluster))*pseudocell_count-0.001))
    prop_m_hardCluster=lapply(prop_m_hardCluster,function(x){
      x=data.frame(pseudocell=x[1],cells=x,stringsAsFactors = F)
      return(x)
    })
    prop_m_hardCluster=do.call("rbind",prop_m_hardCluster)
    
    prop_m_hardCluster$name=prop_m_hardCluster[,1]
    prop_m_hardCluster$name=gsub(" ",".",as.character(prop_m_hardCluster$name))
    
    prop_m_hardCluster$name=gsub("[[:punct:]]+", ".", prop_m_hardCluster$name)
    
    m2=.myOneHotFn(inputVector=prop_m_hardCluster$name)
    row.names(m2)=prop_m_hardCluster[,2]
    m2_colname=prop_m_hardCluster[!duplicated(prop_m_hardCluster[,1]),]
    m2_colname=m2_colname[match(colnames(m2),m2_colname$name),]
    colnames(m2)=m2_colname[,1]
    prop_m_hardCluster=t(as.matrix(m2))
    rm(m2)
  } else {
    if(is.null(inputEmbedding)){
      tmpData=suppressWarnings(.extraExport2SeuratFn(inputExp))
      tmpData = NormalizeData(tmpData,verbose =F)
      tmpData = FindVariableFeatures(tmpData, selection.method = "vst", nfeatures = 2000,verbose =F)
      tmpData <- ScaleData(tmpData,verbose =F)
      tmpData <- RunPCA(tmpData,npcs =nPCs,verbose =F)
      inputEmbedding=tmpData@reductions$pca@cell.embeddings[,1:nPCs]
    } else {
      inputEmbedding=inputEmbedding[colnames(inputExp),1:nPCs]
    }
    
    harmony_embeddings=inputEmbedding
    
    
    pseudocell_names=.extra_sconline.FixedSizekmeansFn(harmony_embeddings=harmony_embeddings,nPCs = nPCs,pseudocell_count = pseudocell_count)#,kmeansMethod=kmeans_method)
    
    
    #pca_centroid=res_clust$centers
    
    pca_centroid=harmony_embeddings[pseudocell_names,,drop=F]
    row.names(pca_centroid)=paste0("ps_",1:nrow(pca_centroid))
    sl_pseudo=data.frame(cluster=paste0("ps_",1:nrow(pca_centroid)),pseudocell=pseudocell_names,stringsAsFactors = F)
    
    
    idx=Seurat:::AnnoyBuildIndex(data = harmony_embeddings, metric = "euclidean", n.trees = n.trees)
    nn.ranked.1=Seurat:::AnnoySearch(index = idx, query = harmony_embeddings,k=k.param,include.distance = T,search.k = -1)
    
    affinities=.extra_matrix_rowNorm(input_mat = nn.ranked.1$nn.dists,rowValues = 1/(nn.ranked.1$nn.dists[,2]+0.000001))#Matrix::Diagonal(x=1/(nn.ranked.1$nn.dists[,2]+0.000001)) %*% nn.ranked.1$nn.dists
    affinities@x=-1*affinities@x^2
    affinities@x=exp(affinities@x)
    affinities[,1]=affinities[,2]
    j <- as.numeric(t(nn.ranked.1$nn.idx))
    i <- ((1:length(j)) - 1)%/%n.neighbors + 1
    x=as.numeric(t(affinities))
    adj = sparseMatrix(i = i, j = j, x = x, dims = c(nrow(inputEmbedding),nrow(inputEmbedding)))
    rownames(adj) <- row.names(inputEmbedding)
    colnames(adj)=c(row.names(inputEmbedding))
    adj= .extra_matrix_rowNorm(adj)#Matrix::Diagonal(x=1/rowSums(adj)) %*% adj
    
    
    adj_t=t(adj)
    adj_t=.extra_matrix_rowNorm(adj_t)#Matrix::Diagonal(x=1/rowSums(adj_t)) %*% adj_t
    adj=adj[sl_pseudo$pseudocell,]
    #row.names(adj)=sl_pseudo$cluster
    adj=adj %*% adj_t
    
    colMax_vals_m=qlcMatrix::colMax(adj)
    colMax_vals_m=.extra_matrix_colNorm(input_mat = adj,colValues = 1/as.numeric(colMax_vals_m))#adj %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
    itol=0.5
    tst_mat=Matrix::drop0(colMax_vals_m,tol=itol)
    tst_mat@x=rep(1,length(tst_mat@x))
    while(all(rowSums(tst_mat)>=pseudocell_size)&itol<0.95){
      itol=itol+0.05
      tst_mat=Matrix::drop0(colMax_vals_m,tol=itol)
      tst_mat@x=rep(1,length(tst_mat@x))
    }
    prop_m_hardCluster=Matrix::drop0(colMax_vals_m,tol=itol)
    prop_m_hardCluster=as.data.frame(summary(prop_m_hardCluster))
    
    o=table(prop_m_hardCluster$i)
    o=o[order(as.numeric(o),decreasing = F)]
    
    prop_m_hardCluster$groups=factor(as.character(prop_m_hardCluster$i),levels=names(o))
    prop_m_hardCluster=prop_m_hardCluster[order(prop_m_hardCluster$groups,decreasing = F),]
    prop_m_hardCluster=prop_m_hardCluster[!duplicated(prop_m_hardCluster$j),]
    
    prop_m_hardCluster=sparseMatrix(i = prop_m_hardCluster$i, j = prop_m_hardCluster$j, x = rep(1,nrow(prop_m_hardCluster)), dims = c(nrow(adj), ncol(adj)))
    row.names(prop_m_hardCluster)=row.names(adj)
    colnames(prop_m_hardCluster)=row.names(inputEmbedding)
    prop_m_hardCluster=prop_m_hardCluster[which(rowSums(prop_m_hardCluster)>0),]
    
  }
  expData=t(counts(inputExp))
  expData=t(prop_m_hardCluster %*% expData[colnames(prop_m_hardCluster),])
  
  #inputData
  #sl_pseudo=sl_pseudo[match(row.names(prop_m_hardCluster),sl_pseudo$pseudocell),]
  
  if(is.null(cols_to_sum)&sum(colnames(colData(inputExp)) %in% cols_to_sum)==0){
    pd=as.data.frame(colData(inputExp))
    pd=pd[match(row.names(prop_m_hardCluster),colnames(inputExp)),]
    pd$pseudocell_size=rowSums(prop_m_hardCluster)
    
    
  } else {
    if(length(setdiff(cols_to_sum,colnames(colData(inputExp))))>0){
      warning("some of provided cols_to_sum cols were not identified in the dataset")
    }
    
    pd=as.data.frame(colData(inputExp))[,!colnames(colData(inputExp)) %in% cols_to_sum]
    pd=pd[match(row.names(prop_m_hardCluster),colnames(inputExp)),]
    pd$pseudocell_size=rowSums(prop_m_hardCluster)
    
    
    sums_res=prop_m_hardCluster %*% as.matrix(as.data.frame(colData(inputExp)[,cols_to_sum]))
    if(any(row.names(sums_res)!=row.names(pd),na.rm = F)){
      print("Error in summing the cols!")
    }
    pd=cbind(pd,sums_res)
  }
  
  
  
  fd=as.data.frame(rowData(inputExp))
  
  pd$pseudocell_outlier_score=NA
  
  if(include.outlier.score&ncol(expData)>3){
      bkg_genes=counts(inputExp)
      bkg_genes=rowSums(bkg_genes>0)/max(ncol(bkg_genes),10)
      if(sum(bkg_genes>0.1)>100){
        bkg_genes=row.names(expData)[bkg_genes>0.1]
        
        logCPM=edgeR::cpm(expData[bkg_genes,],normalized.lib.sizes = F, log = TRUE, prior.count = 1)
        
        tocheck=require(WGCNA,quietly = T)
        if(!tocheck){
          stop("WGCNA package is missing!")
        }
        
        normadj <- (0.5+0.5*bicor(logCPM, use='pairwise.complete.obs'))^2
        netsummary <- fundamentalNetworkConcepts(normadj)
        outlier_score=scale(netsummary$Connectivity)
        pd$pseudocell_outlier_score=as.numeric(outlier_score)
      }
  }
  
  res=SingleCellExperiment(assays = list(counts = expData),colData = pd,rowData=fd)
  
  return(res)
}


.extra_sconline.PseudobulkFn=function(inputData,colName,mode="sum",min_size_limit=20,ncores=1,cols_to_sum=NULL){
  
  #colName: the column in the meta/pheno data that specifies the pseudobulk level.
  #colName: usually a column that is a combination of subjectId/libraryId + cellType
  #min_size_limit: the minimum acceptable size of the pseudobulk data.
  #min_size_limit: pseudobulks with cells less than this threshold are excluded from the analysis 
  
  #pdlist=split(as.data.frame(colData(inputData)),colData(inputData)[,colName])
  #res=matrix(0,nrow=nrow(inputData),ncol=length(pdlist))
  #colnames(res)=names(pdlist)
  
  require(Matrix)
  
  design_mat=as.matrix(.myOneHotFn(colData(inputData)[,colName]))
  if(!is.null(min_size_limit)){
    design_mat=design_mat[,colSums(design_mat)>=min_size_limit,drop=F]
  }
  
  design_mat=t(design_mat)
  design_mat=as(design_mat,"dgCMatrix")
  
  if(mode=="sum"){
    agg_mat=design_mat %*% t(counts(inputData))
    agg_mat=t(agg_mat)
    
  } else if (mode=="mean"){
    design_mat=Matrix::diag(x=1/rowSums(design_mat)) %*% design_mat
    agg_mat=design_mat %*% t(counts(inputData))
    agg_mat=t(agg_mat)
  }
  
  if(is.null(cols_to_sum)&sum(colnames(colData(inputData)) %in% cols_to_sum)==0){
    res_pd=as.data.frame(colData(inputData))
    res_pd=res_pd[!duplicated(res_pd[,colName]),]
    row.names(res_pd)=res_pd[,colName]
    res_pd=res_pd[match(colnames(agg_mat),row.names(res_pd)),]
    
  } else {
    if(length(setdiff(cols_to_sum,colnames(colData(inputData))))>0){
      warning("some of provided cols_to_sum cols were not identified in the dataset")
    }
    res_pd=as.data.frame(colData(inputData)[,!colnames(colData(inputData)) %in% cols_to_sum])
    res_pd=res_pd[!duplicated(res_pd[,colName]),]
    row.names(res_pd)=res_pd[,colName]
    res_pd=res_pd[match(colnames(agg_mat),row.names(res_pd)),]
    
    sums_res=design_mat %*% as.matrix(as.data.frame(colData(inputData)[,cols_to_sum]))
    if(any(row.names(sums_res)!=row.names(res_pd),na.rm = F)){
      print("Error in summing the cols!")
    }
    res_pd=cbind(res_pd,sums_res)
  }
  
  
  res=SingleCellExperiment(assays = list(counts = agg_mat),colData = res_pd,rowData=as.data.frame(rowData(inputData)))
  
  res$QC_Gene_total_count=apply(counts(res),2,sum)
  res$QC_Gene_unique_count=apply(counts(res),2,function(x) sum(x>0))
  lib_sizes=as.data.frame(table(colData(inputData)[,colName]))
  lib_sizes=lib_sizes[match(colData(res)[,colName],lib_sizes[,1]),]
  res$pseudocell_size=lib_sizes[,2]
  
  return(res)
  
}


.mySplitObject_v2=function(
  object,
  colName,
  min_dataset_size,
  ncores=5){

    # If the input object is a Seurat object, use Seurat's SplitObject function
    # to split the object based on the specified column name (colName).
    if(class(object)=="Seurat"){
      res=Seurat::SplitObject(object,split.by = colName)
    } else{
      # pd is the split-column data from the object 
    pd=colData(object)[,colName]
    
    # res starts as a long df representation of the sparse counts matrix
    res_sparse_long=as.data.frame(summary(counts(object)))

    # Add an anno column to the above data frame with annotations based on the previously extracted column.
    # and split the data frame into a list of data frames based on the annotation.
    res_sparse_long$anno=pd[res_sparse_long[,2]]
    res_list=split(res_sparse_long, pd[res_sparse_long[,2]])

    # Similarly, split the column data of the object based on the same annotation.
    # this creates a list of phenotypic data frames, all with the same value of `colName`
    res_pd=split(as.data.frame(colData(object)), pd)

    #subset to only include groupings present in the res_list    
    res_pd = res_pd[names(res_list)]


    output_sce_list=parallel::mclapply(1:length(res_list),function(x){

      # For each subset, create a sparse matrix with dimensions corresponding to the original object's num genes
      # and with columns corresponding to the number of rows in the corresponding element of res_pd. 
      # The values in the matrix are taken from the third column of the subset, representing ___
      y=Matrix::sparseMatrix(i = res_list[[x]][,1],
                              j = as.numeric(as.factor(as.character(res_list[[x]][,2]))),
                              x = res_list[[x]][,3],
                              dims = c(nrow(object),nrow(res_pd[[x]])))
      
      # Depending on the number of columns in the sparse matrix, create a SingleCellExperiment object.
      # This branch handles the case where the matrix has only one column.
      #res_pd[[x]]$pseudocell_size=nrow(res_pd[[x]])
      if(ncol(y)==1){#.5*min_dataset_size){
        #y=as(matrix(rowSums(y),ncol=1),"dgCMatrix")
        
        y=SingleCellExperiment(assays = list(counts = y),colData = res_pd[[x]][1,],rowData=as.data.frame(rowData(object)))
        
      } else {
        # If the matrix has more than one column, set the row names and column names according to the
        # object and the subset of res_pd, then create a SingleCellExperiment object.
        row.names(y)=row.names(object)
        colnames(y)=row.names(res_pd[[x]])
        y=SingleCellExperiment(assays = list(counts = y),colData = res_pd[[x]],rowData=as.data.frame(rowData(object)))
        
      }
      return(y)
      
    },mc.cores=ncores)
    
    if(F){
      for(i in unique(pd)){
        res=c(res,list(object[,which(pd==i)]))
        names(res)[length(res)]=i
      }
    }
    
  }
  
  # Filter the results to include only those with a number of columns (cells) greater than or equal to a minimum size threshold.  
  size_dist=unlist(lapply(output_sce_list,ncol))
  
  output_sce_list=output_sce_list[which(size_dist>=min_dataset_size)]
  
  return(output_sce_list)
}


.sconline.PseudobulkGeneration=function(
    argList=NULL,
    n_clusters=NULL, 
    parsing.col.names = c("anno_batch"), 
    use.sconline.cluster4parsing=T,
    cluster_obj=NULL,
    pseudocell.size=40,
    inputExpData=NULL,
    min_size_limit=20,
    inputPhenoData=NULL,
    inputEmbedding=NULL,
    tol_level=0.9,
    use.sconline.embeddings=F,
    nPCs=NULL,
    ncores=5,
    calculate.outlier.score=F,
    rand_pseudobulk_mod=T,
    organism,
    cols_to_sum=NULL){
  
  
  #parsing.col.names: The columns in the pheonData that will be used to parse the expression data and generate the pseudocell/pseudobulk data
  #use.sconline.cluster4parsing: use the sconline cluster as factor for the parsing
  #inputPhenoData: in case we want to run the function outside sconline space
  #min_size_limit: minimum acceptable size (ie, #cells) for each pseudobulk
  #pseudocell.size: average pseudocell size.
  #if pseudocell.size=null, function turns to a tranditional pseudobulk method, ie, all cells at each parsing level are combined together
  #nPCs: the dimension of the embedding space for the construction of pseudobulk data
  #inputEmbedding: the embedding space to be used for the generation of the pseudobulk. only needed when pseudocell.size is not null
  #use.sconline.embeddings: use the embedding space used for sconline to generate the pseudobulk samples. 
    # Only needed when the pseudocell.size is not set to null
  
  #Function adds/modifies three annotation columns: pseuodcell_size, QC_Gene_total_count, QC_Gene_unique_count
  #QC_Gene_total_count: equivalant to nUMI for the pseudobulk samples
  #QC_Gene_unique_count: equivalant to nGene for the pseudobulk samples
  #use scale(tst$QC_Gene_total_count) and scale(tst$pseudocell_size) as additional covariates for the DE analysis
  
  if(is.null(inputPhenoData)&is.null(argList)){
    stop("Both argList and inputPhenoData cannot be null!")
  }
  
  if(!is.null(pseudocell.size)){
    if(pseudocell.size<2){
      stop("pseudocell.size cannot be less than 2!")
    }
  }
  
  if(is.null(inputExpData)&!rand_pseudobulk_mod){
    warning("It's advised to provide the inputEmbedding")
  }
  
  if(is.null(nPCs)&!is.null(pseudocell.size)&!rand_pseudobulk_mod){
    if(!is.null(argList)){
      nPCs=argList$nPCs
    } else if(!is.null(inputEmbedding)){
        warning(paste0("Setting the nPCs to ", nrow(inputEmbedding)," based on the inputEmbedding"))
    } else {
      stop("nPCs argument (number of PCs for the generation of embedding) need to be provided")
    }
    
    if(nPCs>pseudocell.size){
      warning(paste0("nPCs larger than pseuodcell.size is not advised. setting nPCs to ",pseudocell.size-5))
      nPCs=pseudocell.size-5
    }
  }
  
  if(is.null(argList)){
    use.sconline.cluster4parsing=F
    use.sconline.embeddings=F
  }
  
  
  if(is.null(inputPhenoData)){
    load(.myFilePathMakerFn("UMAP_anno",argList=argList,pseudoImportant = F))
  } else {
    pd=inputPhenoData
  }
  
  
  if(use.sconline.cluster4parsing){
    parsing.col.names=c(parsing.col.names,"cluster_anno_res")
    prop_mat=qread(.myFilePathMakerFn("res_prop_mat_merged",argList=argList,uniformImportant=T,propImportant = T,qsFormat=T))
    prop_mat=.extra_matrix_rowNorm(prop_mat)#Matrix::Diagonal(x = 1 / (rowSums(prop_mat)+0.000000000001)) %*% prop_mat
    
    if(!is.null(n_clusters)){
      cat(paste0("Analysis based on ",n_clusters," clusters"))
      if(is.null(cluster_obj)){
        stop("Cluster object needs to be provided!")
      }
      diff_clust=cluster_obj$cluster_object
      d_conMat=cutree(diff_clust,k=n_clusters)
      prop_merged=t(as.matrix(.myOneHotFn(inputVector=as.factor(d_conMat))))
      prop_merged=prop_merged %*% prop_mat[colnames(prop_merged),pd$sample]
      prop_merged=.extra_matrix_rowNorm(prop_merged)#Matrix::Diagonal(x=1/rowSums(prop_merged)) %*% prop_merged
      #prop anno
      
      
    } else {
      cat(paste0("Analysis at the pseudocell level"))
      prop_merged=prop_mat
    }
    
    colMax_vals_m=qlcMatrix::colMax(prop_merged)
    #prop_merged %*% Matrix::Diagonal(x=1/as.numeric(colMax_vals_m))
    colMax_vals_m=.extra_matrix_colNorm(prop_merged,colValues = 1/as.numeric(colMax_vals_m))
    prop_m_hardCluster=colMax_vals_m=Matrix::drop0(colMax_vals_m,tol=tol_level)
    prop_m_hardCluster=as.data.frame(summary(prop_m_hardCluster))
    prop_m_hardCluster[prop_m_hardCluster$j %in% prop_m_hardCluster$j[duplicated(prop_m_hardCluster$j)],"x"]=0
    prop_m_hardCluster=prop_m_hardCluster[!duplicated(prop_m_hardCluster$j),]
    prop_m_hardCluster=prop_m_hardCluster[match(1:nrow(pd),prop_m_hardCluster$j),]
    pd$cluster_anno_res=paste0("C",as.character(prop_m_hardCluster$i))
  }
  
  # Input data handling, ensuring as ExpressionSet?
  if(is.null(inputExpData)){
    if(!file.exists(.myFilePathMakerFn("exp_merged",argList=argList,expData=T,qsFormat=T))){
      stop("Expression data is missing!")
    }
    inputExpData=qread(.myFilePathMakerFn("exp_merged",argList=argList,expData=T,qsFormat=T))
    inputExpData=.extraExport2ExpressionSetFn(
      counts=inputExpData@assays$RNA@counts,
      pd=as.data.frame(inputExpData@meta.data),
      fd=as.data.frame(inputExpData@assays$RNA@meta.features))
  } else if(tolower(class(inputExpData))=="seurat"){
    inputExpData=.sconline.convertSeuratToExpressionSet(object=inputExpData)
  } else if(class(inputExpData)!="SingleCellExperiment") {
    stop("Unrecognized inputExpression data!")
  }
  
  inputExpData=inputExpData[,colnames(inputExpData) %in% row.names(pd)]
  if(ncol(inputExpData)==0){
    stop("Expression data doesn't match with the phenoData!")
  }
  pd=pd[match(colnames(inputExpData), row.names(pd)),]
  
  for(icol in parsing.col.names){
    if(sum(colnames(colData(inputExpData))==icol)==0){
      if(sum(colnames(pd)==icol)==0){
        stop(paste0(icol," column was not identified!"))
      } else {
        if(class(pd[,icol])[1]==class(factor())){
          colData(inputExpData)[,icol]=as.character(pd[,icol])
        } else {
          colData(inputExpData)[,icol]=pd[,icol]
        }    
      }
    }
  }
  
  # set the lib_anno column, which is used for the split below
  if(length(unique(parsing.col.names))==1){
    inputExpData$lib_anno=colData(inputExpData)[,unique(parsing.col.names)]
  } else {
    inputExpData$lib_anno=apply(as.data.frame(colData(inputExpData)[,unique(parsing.col.names)]),1,function(x)paste(x,collapse="_"))
  }

  if(is.null(pseudocell.size)){
    #inputData=inputExpData;colName="lib_anno";mode="sum";cols_to_sum=cols_to_sum
    sl_data=.extra_sconline.PseudobulkFn(
      inputData=inputExpData,colName="lib_anno",mode="sum",min_size_limit=min_size_limit,cols_to_sum=cols_to_sum,ncores=ncores)
  } else {
    inputExpList=.mySplitObject_v2(inputExpData,colName="lib_anno",min_dataset_size=min_size_limit,ncores=ncores)
    print(paste("length inputExpList after .mySplitObject_v2", length(inputExpList)))
    #inputExpList2=.mySplitObject_v2(inputExpData,colName="library_anno2",min_dataset_size=min_size_limit,ncores=ncores)
    if(sum(unlist(lapply(inputExpList,class))=="SingleCellExperiment")!=length(inputExpList)){
      cat("Error: consider increasing RAM! re-trying with lower number of cores")
      inputExpList=.mySplitObject_v2(inputExpData,colName="lib_anno",min_dataset_size=min_size_limit,ncores=1)
      if(sum(unlist(lapply(inputExpList,class))=="SingleCellExperiment")!=length(inputExpList)){
        stop("Persistent error: increase RAM!")
      }
    }
    if(use.sconline.embeddings&!rand_pseudobulk_mod){
      load(.myFilePathMakerFn("harmony-embeddings",argList=argList,pseudoImportant = F))
      inputEmbedding=harmony_embeddings[,1:nPCs,drop=F]
    }
    
    data_size=unlist(lapply(inputExpList,ncol))
    inputExpList=inputExpList[data_size>=min_size_limit]
    data_size=unlist(lapply(inputExpList,ncol))
    
    if(F){
      for(icheck in inputExpList){
        #inputExp=icheck;include.outlier.score=calculate.outlier.score
        tst=.extra_sconline.FixedSizeFn(
          inputExp=icheck,inputEmbedding=inputEmbedding,pseudocell_size=pseudocell.size,nPCs=nPCs,include.outlier.score=calculate.outlier.score)
      }
    }
    
    sl_data=suppressWarnings(
      parallel::mclapply(
        inputExpList,
        .extra_sconline.FixedSizeFn,
        inputEmbedding=inputEmbedding,
        pseudocell_size=pseudocell.size,
        nPCs=nPCs,
        include.outlier.score=calculate.outlier.score,
        rand_pseudobulk_mod=rand_pseudobulk_mod,
        mc.cores = ncores,
        cols_to_sum=cols_to_sum))
    if(sum(unlist(lapply(sl_data,class))=="SingleCellExperiment")!=length(inputExpList)){
      cat("Error: consider increasing RAM! re-trying with lower number of cores")
      sl_data=suppressWarnings(
        parallel::mclapply(
          inputExpList,
          .extra_sconline.FixedSizeFn,
          inputEmbedding=inputEmbedding,
          pseudocell_size=pseudocell.size,
          nPCs=nPCs,
          include.outlier.score=calculate.outlier.score,
          cols_to_sum=cols_to_sum,
          mc.cores = 1))
      if(sum(unlist(lapply(sl_data,class))=="SingleCellExperiment")!=length(inputExpList)){
        stop("Persistent error: increase RAM!")
      }
    }
    sl_data_size=lapply(sl_data,function(x) max(x$pseudocell_size))
    sl_data=sl_data[sl_data_size>=min_size_limit]
    sl_data=lapply(sl_data,function(x) x[,which(x$pseudocell_size>=min_size_limit),drop=F])
    
    sl_data=.mycBindFn(sl_data)
    sl_data$QC_Gene_total_count=apply(counts(sl_data),2,sum)
    sl_data$QC_Gene_unique_count=apply(counts(sl_data),2,function(x) sum(x>0))
  }
  
  sl_data$QC_MT.pct=.extraMitoPctFn(inputData = sl_data,organism = organism)
  
  
  return(sl_data)
}


.sconline.Pseudobulk10=function(inputExpData,embeddings,pseudobulk_split_col,min_dataset_size=4){
  #creates pseudobulk samples of median size 10
  
  if(is.null(inputExpData)){
    stop("inputExpData should be provided!")
  } else if(tolower(class(inputExpData))=="seurat"){
    inputExpData=.sconline.convertSeuratToExpressionSet(object=inputExpData)
  } else if(class(inputExpData)!="SingleCellExperiment") {
    stop("Unrecognized inputExpression data!")
  }
  
  #data.list=.mySplitObject(inputExpData,pseudobulk_split_col)
  data.list=.mySplitObject_v2(object=inputExpData,colName=pseudobulk_split_col,min_dataset_size=min_dataset_size)
  
  
  embeddings=embeddings[match(colnames(inputExpData),row.names(embeddings)),]
  if(sum(is.na(embeddings))>0){
    stop("Error!")
  }
  
  pseudo_res=list()
  for(i in 1:length(data.list)){
    iclust=unique(colData(data.list[[i]])[,pseudobulk_split_col])
    {
      tmpExp=data.list[[i]]#[,colData(data.list[[i]])[,'pseudobulk_split_col']==iclust]
      #tmpExp=SingleCellExperiment(assays = list(counts = tmpExp),colData = as.data.frame(colData(data.list[[i]]@meta.data[data.list[[i]]@meta.data$cluster_name==iclust,]),rowData=as.data.frame(data.list[[i]]@assays$RNA@meta.features))
      if(length(setdiff(colnames(tmpExp),row.names(embeddings)))>0){
        stop("Error in input embeddings file. Embedding info was not found for some cells!")
      }
      tmp_embeddings=embeddings[match(colnames(tmpExp),row.names(embeddings)),]
      #inputExpData=tmpExp;inputPCAembeddings=tmp_embeddings; n.adaptiveKernel=5; nPropIter=3;nPCs=30;verbose=T
      icounter=0
      runCheck=T
      
      best_coverage=0
      sl_solution=NULL
      
      while(icounter<100&runCheck){
        icounter=icounter+1
        runCheck=F
        if(ncol(tmpExp)<=25&ncol(tmpExp)>13){
          s1=sample(ncol(tmpExp),ncol(tmpExp)/2)
          s2=as.matrix(counts(tmpExp))[,-s1]
          s1=as.matrix(counts(tmpExp))[,s1]
          tmp=matrix(0,nrow=nrow(s1),ncol=2)
          tmp[,1]=rowSums(s1)
          tmp[,2]=rowSums(s2)
          colnames(tmp)=paste0("group_",colnames(tmpExp)[1:2])
          sl_solution=tmp
        } else if(ncol(tmpExp)<=13){
          s1=as.matrix(counts(tmpExp))
          tmp=matrix(rowSums(s1),nrow=nrow(s1),ncol=1)
          colnames(tmp)=colnames(tmpExp)[1]
          colnames(tmp)=paste0("group_",colnames(tmp))
          sl_solution=tmp
        } else {
          
          
          tmp=tryCatch({.myPseudoCellfn(inputExpData=tmpExp,inputPCAembeddings=tmp_embeddings, n.adaptiveKernel=5, nPropIter=3,nPCs=30,verbose=F)}, error=function(e) {return(T)})
          if(class(tmp)==class(T)){
            runCheck=T
          } else {
            if(tmp$coverage>best_coverage){
              sl_solution=tmp$collapsedExpData
              best_coverage=tmp$coverage
            }
            
            if(best_coverage<0.75){
              runCheck=T
            }
          }
        }
        
      }
      
      
      
      if(runCheck&is.null(sl_solution)){
        stop("Error!")
      } else {
        tmp=sl_solution
        if(class(tmp)[1]=="numeric"){
          tmp=matrix(tmp,ncol=1)
          colnames(tmp)=colnames(tmpExp)[1]
        } else if(sum(grepl("^group",colnames(tmp)))==1){
          
          tmp2=matrix(tmp[,grepl("^group",colnames(tmp))],ncol=1)
          colnames(tmp2)=gsub("^group_","",colnames(tmp)[grepl("^group",colnames(tmp))])
          tmp=tmp2
        } else {
          tmp=tmp[,grepl("^group",colnames(tmp))]
          colnames(tmp)=gsub("^group_","",colnames(tmp))
        }
      }
      
      
      
      tmp=SingleCellExperiment(assays = list(counts = tmp),colData = as.data.frame(colData(data.list[[i]])[match(colnames(tmp),row.names(colData(data.list[[i]]))),]),rowData=as.data.frame(rowData(data.list[[i]])))
      
    }
    
    pseudo_res=c(pseudo_res,list(tmp))
    names(pseudo_res)[length(pseudo_res)]=iclust
  }
  
  pseudo_res=.mycBindFn(pseudo_res,batchNames = NULL)
  return(pseudo_res)
}
