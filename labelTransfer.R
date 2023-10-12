source("~/sc-online/utils.R")

.myLabelTransfer_harmony=function(
  dataset_source,
  dataset_target,
  source_label_col,
  indScaling,
  return_seurat_obj=T,
  inputVarGenes=NULL,
  ds_specific_hvg=T,
  target_label_col=NULL,
  source_batch_label_col=NULL,
  target_batch_label_col=NULL,
  covariates=NULL,
  prenatalDataset=F,
  calculate_depth_per_gene=T,
  doubletGroups=NULL,
  nPCs=30,
  n.adaptiveKernel=5,
  nPropIter=3,
  ncores=1){
  
  #dataset_source: dataset whose labels will be transfered
  #dataset_target: dataset for which we are predicting the cell types 
  #source_label_col: the column that specifies the cell type in the source dataset
  #target_label_col: the column that specifies the cell type in the target dataset
  #nPCs: number of PCs for the construction of knn network
  #covariates: covariates to be regressed out from the dataset
  
  if(prenatalDataset){
    covariates=unique(c(covariates,"depth_per_gene"))
  }
  
  if(is.null(target_label_col)){
    dataset_target$target_lbl="unknown"
    target_label_col="target_lbl"
  }
  
  if(class(dataset_source)!="SingleCellExperiment"){
    stop("Source dataset should be in the format of SingleCellExperiment")
  }
  
  if(class(dataset_target)!="SingleCellExperiment"){
    stop("Target dataset should be in the format of SingleCellExperiment")
  }
  
  
  if(!is.null(doubletGroups)){
    dataset_source=.extraDoubletMakerFn(
      inputData=dataset_source,
      label_col=source_label_col,
      sel_labels=doubletGroups)
  }
  
  dataset_source$madeCluster=colData(dataset_source)[,source_label_col]
  dataset_target$madeCluster=colData(dataset_target)[,target_label_col]
  
  #dsSource=dataset_source;dsTarget=dataset_target;cellType_source=colData(dataset_source)[,source_label_col];cellType_target=colData(dataset_target)[,target_label_col];covariates=covariates;calculate_depth_per_gene=calculate_depth_per_gene;source_batch_label=source_batch_label_col;target_batch_label=target_batch_label_col;indScaling=indScaling
  res_seurat=.extraHarmony_dataset_integratorFn(
    dsSource=dataset_source,
    dsTarget=dataset_target,
    cellType_source=colData(dataset_source)[,source_label_col],
    indScaling=indScaling,
    cellType_target=colData(dataset_target)[,target_label_col],
    covariates=covariates,
    calculate_depth_per_gene=calculate_depth_per_gene,
    source_batch_label=source_batch_label_col,
    target_batch_label=target_batch_label_col,
    inputVarGenes=inputVarGenes,
    ds_specific_hvg=ds_specific_hvg)
  
  if(sum(grepl("^inferred_",colnames(res_seurat@meta.data)))>0){
    colnames(res_seurat@meta.data)=gsub("^inferred_","org_inferred_",colnames(res_seurat@meta.data))
  }
  
  meta_data=as.data.frame(res_seurat@meta.data)
  if(sum(is.na(meta_data$madeCluster))>0){
    meta_data$madeCluster[is.na(meta_data$madeCluster)]=""
  }
  
  label_col="madeCluster"
  training_idx=which(res_seurat$batch_label=="source")
  
  harmony_embeddings <- harmony::HarmonyMatrix(
    res_seurat@reductions$pca@cell.embeddings[,1:nPCs], 
    meta_data,
    'batch_label2',
    do_pca = FALSE,
    verbose=FALSE)
  
  res=.myKnnLabelTransferFn(
    inputPCAembeddings=harmony_embeddings,
    meta_data=meta_data,
    training_idx=training_idx,
    label_col=label_col,
    n.adaptiveKernel=n.adaptiveKernel,
    nPropIter=nPropIter)
  
  resTest=res$test_labels
  
  slCols=colnames(resTest)[which(grepl("^inferred_",colnames(resTest)))]
  
  plotTestDf=NULL
  for(i in slCols){
    tmp=data.frame(label=resTest[,'madeCluster'],inferred=gsub("inferred_","",i),weight=resTest[,i],stringsAsFactors = F)
    plotTestDf=rbind(plotTestDf,tmp)
  }
  
  plotTestDf2=aggregate(weight~label+inferred,data=plotTestDf,sum)
  
  plotTestDf2=plotTestDf2[order(plotTestDf2$weight,decreasing = T),]
  
  plotTrainingDf=NULL
  for(i in slCols){
    tmp=data.frame(
      label=res$training_labels[,'madeCluster'],
      inferred=gsub("inferred_","",i),
      weight=res$training_labels[,i],
      stringsAsFactors=F)
    plotTrainingDf=rbind(plotTrainingDf,tmp)
  }
  
  plotTrainingDf2=aggregate(weight~label+inferred,data=plotTrainingDf,sum)
  
  res=c(res,
        list(
          dSource=plotTrainingDf2,
          dTarget=plotTestDf2,
          source_label_col=source_label_col,
          harmony_embeddings=harmony_embeddings)
      )
  
  if(return_seurat_obj){
    umap_res=.reductionUMAPFn(inputPCAembeddings = harmony_embeddings,umap.method = "uwot")
    umap_res=umap_res$embedding
    
    umap_res=umap_res[match(row.names(res$combined_labels), row.names(umap_res)),]
    if(sum(colnames(res$combined_labels) %in% c("UMAP_1","UMAP_2"))>0){
      res$combined_labels=res$combined_labels[,!grepl("UMAP",colnames(res$combined_labels))]
    }
    res$combined_labels=cbind(res$combined_labels,umap_res)
    sl_cols=colnames(res$combined_labels)[grepl("^inferred\\_",colnames(res$combined_labels))]
    res$combined_labels$anno=apply(res$combined_labels[,sl_cols],1,function(x) sl_cols[which(x==max(x,na.rm=T))[1]])
    
    arranged_data=Seurat::CreateSeuratObject(
      counts=res_seurat@assays$RNA@counts,
      project = "SeuratProject",
      assay = "RNA",
      min.cells = 0,
      min.features = 0,
      names.field = 1,
      names.delim = "-",
      meta.data = as.data.frame(res$combined_labels))
    
    if(ncol(res_seurat@assays$RNA@meta.features)>0){
      arranged_data@assays$RNA@meta.features=res_seurat@assays$RNA@meta.features
    }
    
    umapdata=as.matrix(umap_res[,c("UMAP_1","UMAP_2")])
    umap.reduction <- CreateDimReducObject(embeddings = umapdata, 
                                           key = "UMAP_", assay = "RNA", global = TRUE)
    arranged_data[["umap"]]=umap.reduction
    
    if(sum(!grepl("^PC_",colnames(harmony_embeddings)))>0){
      colnames(harmony_embeddings)=gsub("PC_","",colnames(harmony_embeddings))
    }
    
    reduction.data <- CreateDimReducObject(
      embeddings=harmony_embeddings, 
      assay="RNA", 
      key="PC_")
    arranged_data[["pca"]] <- reduction.data
    
    res=c(res,list(seurat_obj=arranged_data))
    
  }
  
  return(res)
}

.myLabelTransfer_aligned=function(pca_source,meta_source,pca_target,meta_target,source_label_col,source_data=NULL,target_data=NULL,target_label_col=NULL,nPCs=NULL,n.adaptiveKernel=5,nPropIter=3,return_seurat_obj=T){
  
  #pca_source: aligned pc space for the source dataset
  #meta_source: meta data for the source dataset
  #pca_target: aligned pc space for the target dataset
  #meta_target: meta data for the target dataset
  #source_label_col: the column that specifies the cell type in the source dataset
  #target_label_col: the column that specifies the cell type in the target dataset
  #nPCs: number of PCs for the construction of knn network; if NULL, set to ncol(pca_source)
  
  if(ncol(pca_source)!=ncol(pca_target)){
    stop("There should be identical number of PCs between pca_source and pca_target")
  }
  
  if(return_seurat_obj){
    if(is.null(source_data)|is.null(target_data)){
      stop("source_data and target_data should be provided!")
    }
  }
  
  if(nrow(pca_source)!=nrow(meta_source)){
    stop("meta data and PC space for the source dataset don't match!")
  }
  
  if(nrow(pca_target)!=nrow(meta_target)){
    stop("meta data and PC space for the target dataset don't match!")
  }
  
  
  if(!all(row.names(pca_source)==row.names(meta_source))){
    stop("row.names don't match between meta data and PC space for the source dataset!")
  }
  
  if(!all(row.names(pca_target)==row.names(meta_target))){
    stop("row.names don't match between meta data and PC space for the target dataset!")
  }
  
  if(is.null(target_label_col)){
    meta_target$target_lbl="unknown"
    target_label_col="target_lbl"
  }
  
  if(sum(colnames(meta_source)==source_label_col)==0){
    stop("source_label_col could not be found!")
  }
  
  meta_source$madeCluster=meta_source[,source_label_col]
  meta_target$madeCluster=meta_target[,target_label_col]
  
  meta_source$batch_label="source"
  meta_target$batch_label="target"
  
  meta_data=plyr::rbind.fill(meta_source,meta_target)
  row.names(meta_data)=c(row.names(meta_source),row.names(meta_target))
  
  if(sum(is.na(meta_data$madeCluster))>0){
    meta_data$madeCluster[is.na(meta_data$madeCluster)]=""
  }
  
  if(sum(grepl("^inferred_",colnames(meta_data)))>0){
    colnames(meta_data)=gsub("^inferred_","org_inferred_",colnames(meta_data))
  }
  
  label_col="madeCluster"
  training_idx=which(meta_data$batch_label=="source")
  
  harmony_embeddings <- rbind(pca_source,pca_target)
  
  row.names(harmony_embeddings)=c(row.names(pca_source),row.names(pca_target))
  
  
  res=.myKnnLabelTransferFn(inputPCAembeddings=harmony_embeddings,meta_data=meta_data,training_idx=training_idx,label_col=label_col,n.adaptiveKernel=n.adaptiveKernel,nPropIter=nPropIter)
  
  resTest=res$test_labels
  
  slCols=colnames(resTest)[which(grepl("^inferred_",colnames(resTest)))]
  
  plotTestDf=NULL
  for(i in slCols){
    tmp=data.frame(label=resTest[,'madeCluster'],inferred=gsub("inferred_","",i),weight=resTest[,i],stringsAsFactors = F)
    plotTestDf=rbind(plotTestDf,tmp)
  }
  
  plotTestDf2=aggregate(weight~label+inferred,data=plotTestDf,sum)
  
  plotTestDf2=plotTestDf2[order(plotTestDf2$weight,decreasing = T),]
  
  plotTrainingDf=NULL
  for(i in slCols){
    tmp=data.frame(label=res$training_labels[,'madeCluster'],inferred=gsub("inferred_","",i),weight=res$training_labels[,i],stringsAsFactors = F)
    plotTrainingDf=rbind(plotTrainingDf,tmp)
  }
  
  plotTrainingDf2=aggregate(weight~label+inferred,data=plotTrainingDf,sum)
  
  
  res=c(res,list(dSource=plotTrainingDf2,dTarget=plotTestDf2,source_label_col=source_label_col,harmony_embeddings=harmony_embeddings))
  
  if(return_seurat_obj){
    umap_res=.reductionUMAPFn(inputPCAembeddings = harmony_embeddings,umap.method = "uwot")
    umap_res=umap_res$embedding
    
    umap_res=umap_res[match(row.names(res$combined_labels), row.names(umap_res)),]
    if(sum(colnames(res$combined_labels) %in% c("UMAP_1","UMAP_2"))>0){
      res$combined_labels=res$combined_labels[,!grepl("UMAP",colnames(res$combined_labels))]
    }
    res$combined_labels=cbind(res$combined_labels,umap_res)
    sl_cols=colnames(res$combined_labels)[grepl("^inferred\\_",colnames(res$combined_labels))]
    res$combined_labels$anno=apply(res$combined_labels[,sl_cols],1,function(x) sl_cols[which(x==max(x,na.rm=T))[1]])
    
    res_seurat=.extraExport2SeuratFn(.mycBindFn(list(source_data,target_data)))
    res_seurat=res_seurat[,match(row.names(res$combined_labels),colnames(res_seurat))]
    arranged_data=Seurat::CreateSeuratObject(counts=res_seurat@assays$RNA@counts,
                                             project = "SeuratProject",
                                             assay = "RNA",
                                             min.cells = 0,
                                             min.features = 0,
                                             names.field = 1,
                                             names.delim = "-",
                                             meta.data = as.data.frame(res$combined_labels))
    
    if(ncol(res_seurat@assays$RNA@meta.features)>0){
      arranged_data@assays$RNA@meta.features=res_seurat@assays$RNA@meta.features
    }
    
    umapdata=as.matrix(umap_res[,c("UMAP_1","UMAP_2")])
    umap.reduction <- CreateDimReducObject(embeddings = umapdata, 
                                           key = "UMAP_", assay = "RNA", global = TRUE)
    arranged_data[["umap"]]=umap.reduction
    
    if(sum(!grepl("^PC_",colnames(harmony_embeddings)))>0){
      colnames(harmony_embeddings)=gsub("PC_","",colnames(harmony_embeddings))
    }
    
    reduction.data <- CreateDimReducObject(embeddings = harmony_embeddings, assay = "RNA", key = "PC_")
    arranged_data[["pca"]] <- reduction.data
    
    res=c(res,list(seurat_obj=arranged_data))
    
  }
  
  return(res)
}


# Plotting Functions
# TODO: move plotting functions to separate module
.mycellAssignHeatmap_binary=function(input_labelTransfer_object,confidenceLevel=0.8,target_cluster_col="Cluster"){
  require(ggplot2)
  lbls=input_labelTransfer_object$test_labels[,grepl("inferred_",colnames(input_labelTransfer_object$test_labels))&!grepl("_max",colnames(input_labelTransfer_object$test_labels))]
  lblIndx=apply(lbls,1,function(x) if(sum(!is.na(x))>0){if(max(x,na.rm = T)>confidenceLevel){which(x==max(x))[1]}else{NA}}else{NA})
  lbls=gsub("inferred_","",colnames(lbls))[lblIndx]
  
  lbls=data.frame(cluster=input_labelTransfer_object$test_labels[,target_cluster_col],inferred=lbls,stringsAsFactors = F)
  lbls=as.data.frame(table(lbls$cluster,lbls$inferred,useNA='ifany'))
  colnames(lbls)=c("cluster","inferred","count")
  lbls$proportion=0
  for(i in 1:nrow(lbls)){
    lbls$proportion[i]=lbls$count[i]/sum(lbls$count[lbls$cluster==lbls$cluster[i]],na.rm = T)
  }
  p=ggplot(lbls,aes(inferred,cluster,fill=proportion))+geom_tile(color="black")+scale_fill_gradient(low="white",high="red")+theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=1))
  return(p)
}

.mycellAssignHeatmap_prob=function(input_labelTransfer_object,confidenceLevel=0.8,target_cluster_col="Cluster"){
  require(ggplot2)
  lbls=input_labelTransfer_object$test_labels[,grepl("inferred_",colnames(input_labelTransfer_object$test_labels))&!grepl("_max",colnames(input_labelTransfer_object$test_labels))]
  lbls[lbls<confidenceLevel]=0
  lbls$cluster_name=input_labelTransfer_object$test_labels[,target_cluster_col]
  lbls=reshape2::melt(lbls,"cluster_name")
  lbls=aggregate(value~variable+cluster_name,data=lbls,sum)
  
  lbls_total=aggregate(value~cluster_name,data=lbls,sum)
  colnames(lbls_total)=c("cluster_name","total_count")
  lbls=merge(lbls,lbls_total,by="cluster_name")
  lbls$fraction=lbls$value/lbls$total_count
  colnames(lbls)[colnames(lbls)=="variable"]="Inferred"
  lbls$Inferred=gsub("inferred_","",lbls$Inferred)
  
  p=ggplot(lbls,aes(Inferred,cluster_name,fill=fraction))+geom_tile(color="black")+scale_fill_gradient(low="white",high="red")+theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=1))
  return(p)
}

# Helper Functions
.myKnnLabelTransferFn=function(inputPCAembeddings,meta_data,training_idx,label_col,n.adaptiveKernel=5,nPropIter=3,NNmethod="annoy",n.trees=50){
  
  #training_idx: the row index of training data in the inputPCAembedings
  #label_col: the column in the meta_data that specifies the cell labels (labels of test cells can be set as unknown; the labels of test sets would be excluded from the propagation step)
  
  training_data=inputPCAembeddings[training_idx,]
  training_labels=meta_data[training_idx,]
  
  test_data=inputPCAembeddings[-training_idx,]
  test_labels=meta_data[-training_idx,]
  ref_ind=NULL
  if(NNmethod=="annoy"){
    ref_ind=Seurat:::AnnoyBuildIndex(data = training_data, metric = "euclidean", 
                                 n.trees = n.trees)
    nn.ranked.1=Seurat:::AnnoySearch(index = ref_ind, query = test_data,k=n.adaptiveKernel*12,include.distance = T,search.k = -1)
    
  } else {
    nn.ranked.1 <- RANN::nn2(data = training_data,query = test_data, k = n.adaptiveKernel*12, eps = 0)
  }
  
  
  nn.ranked.1$nn.dists=cbind(rep(0,nrow(nn.ranked.1$nn.dists)),nn.ranked.1$nn.dists)
  nn.ranked.1$nn.idx=cbind(1:nrow(nn.ranked.1$nn.idx)+nrow(training_data),nn.ranked.1$nn.idx)
  nn.ranked.1$nn.dists=nn.ranked.1$nn.dists[,1:(n.adaptiveKernel*12)]
  nn.ranked.1$nn.idx=nn.ranked.1$nn.idx[,1:(n.adaptiveKernel*12)]
  
  if(NNmethod=="annoy"){
    idx=Seurat:::AnnoyBuildIndex(data = test_data, metric = "euclidean", 
                                 n.trees = n.trees)
    nn.ranked.purityIndex=Seurat:::AnnoySearch(index = idx, query = test_data,k=n.adaptiveKernel*3,include.distance = T,search.k = -1)
    
  } else {
    nn.ranked.purityIndex=RANN::nn2(data = test_data, k = n.adaptiveKernel*3, eps = 0)
  }
  
  nn.ranked.purityIndex=nn.ranked.purityIndex$nn.dists[,n.adaptiveKernel*3]
  
  if(NNmethod=="annoy"){
    nn.ranked.2=Seurat:::AnnoySearch(index = ref_ind, query = training_data,k=n.adaptiveKernel*12,include.distance = T,search.k = -1)
    
  } else {
    nn.ranked.2 <- RANN::nn2(data = training_data,query = training_data, k = n.adaptiveKernel*12, eps = 0)
    
  }
  nn.training.purityIndx=nn.ranked.2$nn.dists[,n.adaptiveKernel]
  
  test_labels$status="test_set"
  training_labels$status="training_set"
  labels=rbind(training_labels,test_labels)
  rm(test_labels,training_labels,test_data,training_data)
  
  dists=nn.ranked.1$nn.dists
  affinities=lapply(1:nrow(dists),function(x) exp((-1)*(dists[x,]/nn.ranked.purityIndex[x])^2))
  affinities=do.call("rbind",affinities)
  affinitiesThr=apply(affinities,1,function(x) x[3*n.adaptiveKernel+1])
  affinities2=affinities-affinitiesThr
  affinities[affinities2<=0]=0
  #removing the outlier matches
  toRM= which(scale(affinities[,2])<(-3))
  if(length(toRM)>0){
    for(i in 1:ncol(affinities)){
      affinities[toRM,i]=0
    }
  }
  nn.ranked.1$affinities=affinities
  rm(affinities)
  
  dists=nn.ranked.2$nn.dists
  affinities=lapply(1:nrow(dists),function(x) exp((-1)*(dists[x,]/nn.training.purityIndx[x])^2))
  affinities=do.call("rbind",affinities)
  affinitiesThr=apply(affinities,1,function(x) x[3*n.adaptiveKernel+1])
  affinities2=affinities-affinitiesThr
  affinities[affinities2<=0]=0
  nn.ranked.2$affinities=affinities
  rm(affinities)
  
  nn.ranked.org=nn.ranked.2
  nn.ranked.org$nn.idx=rbind(nn.ranked.2$nn.idx,nn.ranked.1$nn.idx)
  nn.ranked.org$nn.dists=rbind(nn.ranked.2$nn.dists,nn.ranked.1$nn.dists)
  nn.ranked.org$affinities=rbind(nn.ranked.2$affinities,nn.ranked.1$affinities)
  rm(nn.ranked.1,nn.ranked.2)
  
  
  nn.ranked=nn.ranked.org
  affinities=nn.ranked.org$affinities
  
  affCounts=apply(affinities,2,function(x) sum(x==0))
  affinities=affinities[,-which(affCounts==nrow(affinities))]
  nn.ranked$nn.idx=nn.ranked$nn.idx[,-which(affCounts==nrow(affinities))]
  j <- as.numeric(t(nn.ranked$nn.idx))
  i <- ((1:length(j)) - 1)%/%ncol(affinities) + 1
  k=as.numeric(t(affinities))
  
  graph <- sparseMatrix(i = i, j = j, x = k,dims = c(nrow(inputPCAembeddings), nrow(inputPCAembeddings)))
  rm(i,j,k)
  graph <- Matrix::Diagonal(x = 1 / (rowSums(graph)+0.0001)) %*% graph
  
  if(class(graph)!="dgCMatrix"){
    graph=as(graph,"dgCMatrix")
  }
  
  rownames(graph) <- rownames(inputPCAembeddings)
  colnames(graph) <- rownames(inputPCAembeddings)
  
  for(i in 1:nPropIter){
    print(i)
    if(i==1){
      prop_mat=graph
    } else {
      prop_mat=prop_mat%*%graph
    }
  }
  
  
  tmpLables=as.data.frame(labels)[,label_col]
  tmpLables[labels$status=="test_set"]="unknown_toexclude"
  one_hot_labels=.myOneHotFn(tmpLables)
  training_col_labels=unique(as.data.frame(labels)[training_idx,label_col])
  if(length(which(colnames(one_hot_labels)=="unknown_toexclude"))>0){
    one_hot_labels=one_hot_labels[,-which(colnames(one_hot_labels)=="unknown_toexclude")]
  }
  
  
  res=as.matrix(prop_mat %*% as.matrix(one_hot_labels))
  res=round(res,3)
  
  res2=as.data.frame(res)
  colnames(res)=paste0("inferred_",colnames(res))
  
  colnames(res2)=colnames(res)
  res=cbind(labels,res2)
  
  res_testset=res[res$status=="test_set",]
  res_training=res[res$status=="training_set",]
  
  return(list(test_labels=res_testset,training_labels=res_training,combined_labels=res))
}

.extraHarmony_dataset_integratorFn=function(dsSource,dsTarget,cellType_source,cellType_target,covariates=c("depth_per_gene"),calculate_depth_per_gene=T,source_batch_label=NULL,target_batch_label=NULL,inputVarGenes=NULL,ds_specific_hvg=T,indScaling=F){
  require(Seurat)
  require(ggplot2)
  require(ggalluvial)
  
  print("Combining the datasets...")
  
  if(any(row.names(dsSource)!=row.names(dsTarget))){
    stop("Row names should match between the two datasets")
  }
  
  ortholog_mapped=F
  if(sum(colnames(rowData(dsSource))=="mouse_ensembl_gene_id")>0|sum(colnames(rowData(dsTarget))=="mouse_ensembl_gene_id")>0){
    ortholog_mapped=T
    print("Ortholog mapping was detected")
  }
  
  
  if(sum(colnames(colData(dsSource))=="anno_batch")>0){
    dsSource$org_batch=dsSource$anno_batch
    colData(dsSource)=colData(dsSource)[,-which(colnames(colData(dsSource))=="anno_batch")]
    if(sum(colnames(colData(dsTarget))=="anno_batch")==0){
      dsTarget$org_batch="target"
    }
  }
  
  if(sum(colnames(colData(dsTarget))=="anno_batch")>0){
    dsTarget$org_batch=dsTarget$anno_batch
    colData(dsTarget)=colData(dsTarget)[,-which(colnames(colData(dsTarget))=="anno_batch")]
    if(sum(colnames(colData(dsSource))=="anno_batch")==0){
      dsSource$org_batch="source"
    }
  }
  
  data.list =list(source=dsSource,target=dsTarget)
  res=.mycBindFn(data.list,names(data.list))
  data.list =list(source=.extraExport2SeuratFn(dsSource),target=.extraExport2SeuratFn(dsTarget))
  
  
  res$batch_label=res$anno_batch
  res$batch_label2=res$anno_batch
  if(sum(colnames(colData(res))=="org_batch")>0){
    res$anno_batch=res$org_batch
  }
  
  if(!is.null(target_batch_label)){
    res$batch_label2[which(res$batch_label=="target")]=colData(res)[which(res$batch_label=="target"),target_batch_label]
  }
  
  if(!is.null(source_batch_label)){
    res$batch_label2[which(res$batch_label=="source")]=colData(res)[which(res$batch_label=="source"),source_batch_label]
  }
  
  
  res=.extraExport2SeuratFn(res)
  
  
  for (i in 1:length(data.list)) {
    data.list[[i]] <- NormalizeData(data.list[[i]], verbose = FALSE)
    data.list[[i]] <- FindVariableFeatures(data.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  }
  
  res <- NormalizeData(res, normalization.method = "LogNormalize", scale.factor = 10000)
  if(is.null(inputVarGenes)){
    if(ds_specific_hvg){
      var_genes=unique(c(data.list[[1]]@assays$RNA@var.features,data.list[[2]]@assays$RNA@var.features))
      res@assays$RNA@var.features=var_genes
    } else {
      res=FindVariableFeatures(res, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    }
  } else {
    res@assays$RNA@var.features=row.names(res)[which(toupper(row.names(res)) %in% toupper(inputVarGenes))]
  }
  
  if(ortholog_mapped){
    var_genes=as.character(res@assays$RNA@var.features)
    fd=res@assays$RNA@meta.features
    row.names(fd)=row.names(res@assays$RNA@data)
    fd=fd[match(var_genes,row.names(fd)),]
    
    #library(googleCloudStorageR)
    if(!file.exists("~/serverFiles/ortholog_mapping.rda")){
      stop("Error! file ~/serverFiles/ortholog_mapping.rda is missing")
    }
    load("~/serverFiles/ortholog_mapping.rda")
    
    var_genes=row.names(fd)[which(fd$ensembl_gene_id %in% as.character(mapping$human))]
    res@assays$RNA@var.features=var_genes
  }
  
  if(calculate_depth_per_gene){
    res$depth_per_gene=res$QC_Gene_total_count/(res$QC_Gene_unique_count+1)
  }
  
  
  
  if(indScaling){
    dataList=.mySplitObject(res,"batch_label2")
    mySeuratFn2=function(inputData,data){
      inputData@assays$RNA@var.features=data@assays$RNA@var.features
      inputData = ScaleData(inputData, verbose = FALSE,features=data@assays$RNA@var.features)
      return(inputData)
    }
    
    tmpInd=c()
    if(sum(!is.null(covariates))>0){
      for(i in covariates[!is.null(covariates)]){
        tmpInd=c(tmpInd,which(is.na(res@meta.data[,i])))
      }
      if(length(tmpInd)>0){
        print("Removing some cells with the covariate value of NA")
        res=res[,-unique(tmpInd)]
      }
    }
    
    
    dataList=parallel::mclapply(dataList,mySeuratFn2,data=res,mc.cores = length(dataList))
      
      resScaled=dataList[[1]]@assays$RNA@scale.data
      for(i in 2:length(dataList)){
        tmp=dataList[[i]]@assays$RNA@scale.data
        resScaled=cbind(resScaled,dataList[[i]]@assays$RNA@scale.data)
      }
      
      if(!all(colnames(res)==colnames(resScaled))){
        stop("Error in the matching!")
      }
      if(!all(row.names(res)[row.names(res) %in% res@assays$RNA@var.features]==row.names(resScaled))){
        stop("Error in the matching!")
      }
      
      res@assays$RNA@scale.data=resScaled
      
      res = RunPCA(res,features=res@assays$RNA@var.features,verbose = F)
      
  } else {
    res <- ScaleData(res, verbose = FALSE,features=res@assays$RNA@var.features,vars.to.regress =covariates)
    res <- RunPCA(res,verbose = F)
  }
  
  
  return(res)
  
}

.reductionUMAPFn=function (inputPCAembeddings,testPCAembeddings=NULL, umap.method = "umap-learn", n.neighbors = 30L, 
                           n.components = 2L, metric = "cosine", n.epochs = NULL, learning.rate = 1, 
                           min.dist = 0.3, spread = 1, set.op.mix.ratio = 1, local.connectivity = 1L, 
                           repulsion.strength = 1, negative.sample.rate = 5, a = NULL, 
                           b = NULL, uwot.sgd = FALSE, seed.use = 42, metric.kwds = NULL, 
                           angular.rp.forest = FALSE, reduction.key = "UMAP_", verbose = TRUE, 
                           ...) {
  #adapted from Seurat
  require(Seurat)
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  umap.model=""
  if (umap.method != "umap-learn") {
    warning("The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric", 
            "\nTo use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'", 
            call. = FALSE, immediate. = TRUE)
  }
  umap.output <- switch(EXPR = umap.method, `umap-learn` = {
    require(reticulate)
    if (! reticulate::py_module_available(module = "umap")) {
      stop("Cannot find UMAP, please install through pip (e.g. pip install umap-learn).")
    }
    if (!is.null(x = seed.use)) {
      py_set_seed(seed = seed.use)
    }
    if (typeof(x = n.epochs) == "double") {
      n.epochs <- as.integer(x = n.epochs)
    }
    umap_import <- import(module = "umap", delay_load = TRUE)
    umap <- umap_import$UMAP(n_neighbors = as.integer(x = n.neighbors), 
                             n_components = as.integer(x = n.components), metric = metric, 
                             n_epochs = n.epochs, learning_rate = learning.rate, 
                             min_dist = min.dist, spread = spread, set_op_mix_ratio = set.op.mix.ratio, 
                             local_connectivity = local.connectivity, repulsion_strength = repulsion.strength, 
                             negative_sample_rate = negative.sample.rate, a = a, 
                             b = b, metric_kwds = metric.kwds, angular_rp_forest = angular.rp.forest, 
                             verbose = verbose)
    #umap$fit_transform(as.matrix(x = inputPCAembeddings))
    
    tst=umap$fit(as.matrix(x = inputPCAembeddings))
    
    umap.output=tst$transform(inputPCAembeddings)
    
    colnames(umap.output) <- paste0(reduction.key, 1:ncol(umap.output))
    rownames(umap.output) <- rownames(inputPCAembeddings)
    
    res=c()
    if(!is.null(testPCAembeddings)){
      umap.test=tst$transform(testPCAembeddings)
      colnames(umap.test) <- paste0(reduction.key, 1:ncol(umap.test))
      rownames(umap.test) <- rownames(testPCAembeddings)
      res=list(embedding=umap.output,model=tst,test=umap.test)
    } else {
      res=list(embedding=umap.output,model=tst)
    }
    
    res
    
  }, uwot = {
    if (metric == "correlation") {
      warning("UWOT does not implement the correlation metric, using cosine instead", 
              call. = FALSE, immediate. = TRUE)
      metric <- "cosine"
    }
    resModel=uwot::umap(X = inputPCAembeddings, n_threads = 1, n_neighbors = as.integer(n.neighbors), 
               n_components = as.integer(n.components), metric = metric, 
               n_epochs = n.epochs, learning_rate = learning.rate, 
               min_dist = min.dist, spread = spread, set_op_mix_ratio = set.op.mix.ratio, 
               local_connectivity = local.connectivity, repulsion_strength = repulsion.strength, 
               negative_sample_rate = negative.sample.rate, a = a, 
               b = b, fast_sgd = uwot.sgd, verbose = verbose,ret_model = T)
    
    umap.model=resModel
    umap.output=resModel[["embedding"]]
    
    colnames(umap.output) <- paste0(reduction.key, 1:ncol(umap.output))
    rownames(umap.output) <- rownames(inputPCAembeddings)
    
    res=c()
    if(!is.null(testPCAembeddings)){
      umap.test=uwot::umap_transform(testPCAembeddings,resModel)
      res=list(embedding=umap.output,model=resModel,test=umap.test)
    } else {
      res=list(embedding=umap.output,model=umap.model)
    }
    
    res
  }, stop("Unknown umap method: ", umap.method, call. = FALSE))
  
  
  
  return(umap.output)
}

.reductionUMAPFn=function (inputPCAembeddings,testPCAembeddings=NULL, umap.method = "umap-learn", n.neighbors = 30L, 
                           n.components = 2L, metric = "cosine", n.epochs = NULL, learning.rate = 1, 
                           min.dist = 0.3, spread = 1, set.op.mix.ratio = 1, local.connectivity = 1L, 
                           repulsion.strength = 1, negative.sample.rate = 5, a = NULL, 
                           b = NULL, uwot.sgd = FALSE, seed.use = 42, metric.kwds = NULL, 
                           angular.rp.forest = FALSE, reduction.key = "UMAP_", verbose = TRUE, 
                           ...) {
  #adapted from Seurat
  require(Seurat)
  if (!is.null(x = seed.use)) {
    set.seed(seed = seed.use)
  }
  umap.model=""
  if (umap.method != "umap-learn") {
    warning("The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric", 
            "\nTo use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'", 
            call. = FALSE, immediate. = TRUE)
  }
  umap.output <- switch(EXPR = umap.method, `umap-learn` = {
    require(reticulate)
    if (! reticulate::py_module_available(module = "umap")) {
      stop("Cannot find UMAP, please install through pip (e.g. pip install umap-learn).")
    }
    if (!is.null(x = seed.use)) {
      py_set_seed(seed = seed.use)
    }
    if (typeof(x = n.epochs) == "double") {
      n.epochs <- as.integer(x = n.epochs)
    }
    umap_import <- import(module = "umap", delay_load = TRUE)
    umap <- umap_import$UMAP(n_neighbors = as.integer(x = n.neighbors), 
                             n_components = as.integer(x = n.components), metric = metric, 
                             n_epochs = n.epochs, learning_rate = learning.rate, 
                             min_dist = min.dist, spread = spread, set_op_mix_ratio = set.op.mix.ratio, 
                             local_connectivity = local.connectivity, repulsion_strength = repulsion.strength, 
                             negative_sample_rate = negative.sample.rate, a = a, 
                             b = b, metric_kwds = metric.kwds, angular_rp_forest = angular.rp.forest, 
                             verbose = verbose)
    #umap$fit_transform(as.matrix(x = inputPCAembeddings))
    
    tst=umap$fit(as.matrix(x = inputPCAembeddings))
    
    umap.output=tst$transform(inputPCAembeddings)
    
    colnames(umap.output) <- paste0(reduction.key, 1:ncol(umap.output))
    rownames(umap.output) <- rownames(inputPCAembeddings)
    
    res=c()
    if(!is.null(testPCAembeddings)){
      umap.test=tst$transform(testPCAembeddings)
      colnames(umap.test) <- paste0(reduction.key, 1:ncol(umap.test))
      rownames(umap.test) <- rownames(testPCAembeddings)
      res=list(embedding=umap.output,model=tst,test=umap.test)
    } else {
      res=list(embedding=umap.output,model=tst)
    }
    
    res
    
  }, uwot = {
    if (metric == "correlation") {
      warning("UWOT does not implement the correlation metric, using cosine instead", 
              call. = FALSE, immediate. = TRUE)
      metric <- "cosine"
    }
    resModel=uwot::umap(X = inputPCAembeddings, n_threads = 1, n_neighbors = as.integer(n.neighbors), 
               n_components = as.integer(n.components), metric = metric, 
               n_epochs = n.epochs, learning_rate = learning.rate, 
               min_dist = min.dist, spread = spread, set_op_mix_ratio = set.op.mix.ratio, 
               local_connectivity = local.connectivity, repulsion_strength = repulsion.strength, 
               negative_sample_rate = negative.sample.rate, a = a, 
               b = b, fast_sgd = uwot.sgd, verbose = verbose,ret_model = T)
    
    umap.model=resModel
    umap.output=resModel[["embedding"]]
    
    colnames(umap.output) <- paste0(reduction.key, 1:ncol(umap.output))
    rownames(umap.output) <- rownames(inputPCAembeddings)
    
    res=c()
    if(!is.null(testPCAembeddings)){
      umap.test=uwot::umap_transform(testPCAembeddings,resModel)
      res=list(embedding=umap.output,model=resModel,test=umap.test)
    } else {
      res=list(embedding=umap.output,model=umap.model)
    }
    
    res
  }, stop("Unknown umap method: ", umap.method, call. = FALSE))
  
  
  
  return(umap.output)
}

