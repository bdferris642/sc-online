.myLabelTransfer_harmony=function(dataset_source,dataset_target,source_label_col,indScaling,return_seurat_obj=T,inputVarGenes=NULL,ds_specific_hvg=T,target_label_col=NULL,source_batch_label_col=NULL,target_batch_label_col=NULL,covariates=NULL,prenatalDataset=F,calculate_depth_per_gene=T,doubletGroups=NULL,nPCs=30,n.adaptiveKernel=5,nPropIter=3){
  
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
    dataset_source=.extraDoubletMakerFn(inputData=dataset_source,label_col=source_label_col,sel_labels=doubletGroups)
  }
  
  dataset_source$madeCluster=colData(dataset_source)[,source_label_col]
  dataset_target$madeCluster=colData(dataset_target)[,target_label_col]
  
  #dsSource=dataset_source;dsTarget=dataset_target;cellType_source=colData(dataset_source)[,source_label_col];cellType_target=colData(dataset_target)[,target_label_col];covariates=covariates;calculate_depth_per_gene=calculate_depth_per_gene;source_batch_label=source_batch_label_col;target_batch_label=target_batch_label_col;indScaling=indScaling
  res_seurat=.extraHarmony_dataset_integratorFn(dsSource=dataset_source,dsTarget=dataset_target,cellType_source=colData(dataset_source)[,source_label_col],indScaling = indScaling,cellType_target=colData(dataset_target)[,target_label_col],covariates=covariates,calculate_depth_per_gene=calculate_depth_per_gene,source_batch_label=source_batch_label_col,target_batch_label=target_batch_label_col,inputVarGenes=inputVarGenes,ds_specific_hvg=ds_specific_hvg)
  
  if(sum(grepl("^inferred_",colnames(res_seurat@meta.data)))>0){
    colnames(res_seurat@meta.data)=gsub("^inferred_","org_inferred_",colnames(res_seurat@meta.data))
  }
  
  meta_data=as.data.frame(res_seurat@meta.data)
  if(sum(is.na(meta_data$madeCluster))>0){
    meta_data$madeCluster[is.na(meta_data$madeCluster)]=""
  }
  
  label_col="madeCluster"
  training_idx=which(res_seurat$batch_label=="source")
  
  harmony_embeddings <- harmony::HarmonyMatrix(res_seurat@reductions$pca@cell.embeddings[,1:nPCs], meta_data, 'batch_label2', do_pca = FALSE, verbose=FALSE)
  
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
