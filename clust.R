library(future)
library(Seurat)
library(harmony)

source("~/sc-online/utils.R")
source("~/sc-online/getData.R")

assignCellClasses = function(
    obj,
    classes,
    cluster_col="seurat_clusters",
    class_col="cell_class") {
    obj@meta.data[[class_col]] = NA  # Assign NA to the entire column in meta.data

    for (i in 1:length(classes)) {
        # Find cells that match the current cluster number
        cells = which(as.numeric(obj@meta.data[[cluster_col]]) == i)
        
        # Assign the class label to the cells in the cluster
        obj@meta.data[cells, class_col] = classes[[i]]
    }
    
    # Return the modified object
    return(obj)
}


getClusterLogFC = function(seurat_obj, cluster){
    clust_genesums = rowSums(seurat_obj[, seurat_obj$seurat_clusters==cluster])
    non_clust_genesums = rowSums(seurat_obj[, seurat_obj$seurat_clusters!=cluster])
    logfc_clust = log2((clust_genesums + 1) / (non_clust_genesums + 1))
    return(sort(logfc_clust, decreasing = TRUE))
}


normalizeScalePcaClusterUmap = function(
    sobj,
    assay = "RNA",
    var_feature_subset_col="participant_id",
    scaling_subset_col="participant_id",
    harmony_group_vars=c("participant_id"),
    n_hvgs_orig=2500, 
    n_dims_use=20, 
    resolutions=c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
    hvgs = NULL,
    regression_vars = NULL,
    run_harmony=FALSE){
    # this function takes a Seurat object, normalizes the data, scales the data, runs PCA, finds neighbors, 
    # clusters at a variety of resolutions, and runs UMAP
    # it then returns the Seurat object

    DefaultAssay(sobj) = assay
    
    if (is.null(hvgs)){
        if (is.null(var_feature_subset_col)){
            sobj = FindVariableFeatures(sobj, nfeatures=n_hvgs_orig)
            hvgs = sobj@assays$RNA@var.features
        } else {
            hvgs = getSeuratVarFeatureIntersectByCol(sobj, subset_col=var_feature_subset_col, original_nfeatures=n_hvgs_orig)
        }
    }
    sobj@assays$RNA@var.features=hvgs
    sobj = (sobj
        %>% NormalizeData() 
        %>% ScaleData(features=hvgs, split.by=scaling_subset_col, vars.to.regress=regression_vars) 
        %>% RunPCA(features=hvgs, npcs=n_dims_use) 
    )
    if (run_harmony){
        sobj = RunHarmony(sobj, group.by.vars=harmony_group_vars)
        reduction_to_use = "harmony"
    } else {
        reduction_to_use = "pca"
    }

    sobj = FindNeighbors(sobj, dims=1:n_dims_use, reduction=reduction_to_use)
    for(res in resolutions){
        sobj = sobj %>% FindClusters(resolution = res)
    }
    sobj = sobj %>% RunUMAP(dims = 1:n_dims_use, reduction = reduction_to_use)
    return(sobj)
}


printMarkersByCluster = function(marker_df, marker_tsv="~/all_markers.tsv", cluster=NULL){
    broad_markers = read.table(marker_tsv, header = TRUE, sep = "\t")
    broad_markers$Gene = toupper(broad_markers$Gene)

    broad_markers['broad_class'] = NA
    broad_markers$broad_class[grepl('microglia', tolower(broad_markers$pattern))] = 'microglia'
    broad_markers$broad_class[grepl('neuron', tolower(broad_markers$pattern))] = 'neuron' 
    broad_markers$broad_class[grepl('astrocyte', tolower(broad_markers$pattern))] = 'astrocyte' 
    broad_markers$broad_class[grepl('oligodendrocyte', tolower(broad_markers$pattern))] = 'oligo' 
    broad_markers$broad_class[grepl('endothelial', tolower(broad_markers$pattern))] = 'endo'
    broad_markers$broad_class[grepl('mural', tolower(broad_markers$pattern))] = 'endo'
    broad_markers$broad_class[grepl('fibro', tolower(broad_markers$pattern))] = 'fibroblast'
    broad_markers$broad_class[grepl('ependymal', tolower(broad_markers$pattern))] = 'ependymal'  
    broad_markers$broad_class[grepl('opc', tolower(broad_markers$pattern))] = 'opc' 
    broad_markers$broad_class[grepl('polydendro', tolower(broad_markers$pattern))] = 'opc'
    broad_markers$broad_class[grepl('b_cell', tolower(broad_markers$pattern))] = 'immune' 
    broad_markers$broad_class[grepl('t_cell', tolower(broad_markers$pattern))] = 'immune' 
    broad_markers$broad_class[grepl('neutro', tolower(broad_markers$pattern))] = 'immune'
    broad_markers$broad_class[grepl('nk_cell', tolower(broad_markers$pattern))] = 'immune' 
    broad_markers = broad_markers[!is.na(broad_markers$broad_class), c('Gene', 'broad_class')]

    broad_markers_ordered = broad_markers[order(broad_markers$broad_class),]

    # Merge broad_markers_ordered with markers
    broad_markers_ordered = merge(
        broad_markers_ordered, marker_df, by.x = "Gene", by.y="gene", all.x = TRUE)

    # Order by broad_class and cluster
    broad_markers_ordered = broad_markers_ordered[
        order(broad_markers_ordered$broad_class, 
        broad_markers_ordered$Gene, broad_markers_ordered$cluster),]

    if (!is.null(cluster)){
        broad_markers_ordered = broad_markers_ordered[broad_markers_ordered$cluster == cluster,]
        broad_markers_ordered = broad_markers_ordered[complete.cases(broad_markers_ordered),]
    }
    broad_markers_ordered$pct.1 = round(broad_markers_ordered$pct.1, 2)
    broad_markers_ordered$pct.2 = round(broad_markers_ordered$pct.2, 2)
    broad_markers_ordered$avg_log2FC = round(broad_markers_ordered$avg_log2FC, 2)
    print(broad_markers_ordered[, 
        c("Gene", "broad_class", "cluster", "avg_log2FC", "pct.1", "pct.2")], row.names=FALSE)
}

.myPCAfn=function(data, argList,projection_data=NULL,saveFiles=T,benchmark_mode=T,...){
  
  require(purrr)
  require(irlba)
  plan("multicore", workers = 12)
  plan()
  options(future.globals.maxSize = 1000 * 1024^4)
  
  UMI_cor_thr=argList$UMI_cor_thr
  
  myPrepDR=function (scaledData, features, verbose = TRUE,projection_state=F) {
    
    data.use <- scaledData
    if (nrow(x = data.use) == 0) {
      stop("Data has not been scaled. Please run ScaleData and retry")
    }
    features.keep <- unique(x = features[features %in% rownames(x = data.use)])
    if (length(x = features.keep) < length(x = features)) {
      features.exclude <- setdiff(x = features, y = features.keep)
      if (verbose) {
        warning(paste0("The following ", length(x = features.exclude),
                       " features requested have not been scaled (running reduction without them): ",
                       paste0(features.exclude, collapse = ", ")))
      }
    }
    features <- features.keep
    # TODO jonah parallize buyt make sure chunked
    features.var <- apply(X = data.use[features, ], MARGIN = 1,
                          FUN = var)
    if(!projection_state){
      features.keep <- features[features.var > 0]
    } else {
      features.keep <- features
    }
    
    if (length(x = features.keep) < length(x = features)) {
      features.exclude <- setdiff(x = features, y = features.keep)
      if (verbose) {
        warning(paste0("The following ", length(x = features.exclude),
                       " features requested have zero variance (running reduction without them): ",
                       paste0(features.exclude, collapse = ", ")))
      }
    }
    features <- features.keep
    features <- features[!is.na(x = features)]
    data.use <- data.use[features, ]
    return(data.use)
  }
  
  mySeuratFn2_org=function(inputData,varFeatures){
    plan("sequential")
    inputData@assays$RNA@var.features=varFeatures
    inputData = ScaleData(inputData, verbose = FALSE,features=varFeatures)
    return(inputData)
  }
  
  mySeuratFn2_archive=function(inputData,varFeatures,scaleData=T){
    plan("sequential")
    inputData@assays$RNA@var.features=varFeatures
    if(scaleData){
      inputData = Seurat:::ScaleData.default(inputData@assays$RNA@data, verbose = FALSE,features=varFeatures)
    } else {
      inputData=inputData@assays$RNA@data[row.names(inputData) %in% varFeatures,]
    }
    
    return(inputData)
  }
  
  internal_scaleFn=function(inputData,varFeatures,scaleData=T){
    plan("sequential")
    
    if(scaleData){
      inputData = Seurat:::ScaleData.default(Seurat::NormalizeData(counts(inputData),verbose=F), verbose = FALSE,features=varFeatures)
    } else {
      inputData=Seurat::NormalizeData(counts(inputData),verbose=F)[row.names(inputData) %in% varFeatures,]
    }
    
    return(inputData)
  }
  
  .mycBindFillFn=function(mat1,mat2){
    
    mat1c=setdiff(row.names(mat2),row.names(mat1))
    mat2c=setdiff(row.names(mat1),row.names(mat2))
    if(length(mat1c)>0){
      mat1cc=matrix(0,nrow=length(mat1c),ncol=ncol(mat1))
      row.names(mat1cc)=mat1c
      mat1=rbind(mat1,mat1cc)
    }
    if(length(mat2c)>0){
      mat2cc=matrix(0,nrow=length(mat2c),ncol=ncol(mat2))
      row.names(mat2cc)=mat2c
      mat2=rbind(mat2,mat2cc)
    }
    mat2=mat2[match(row.names(mat1),row.names(mat2)),]
    mat=cbind(mat1,mat2)
    return(mat)
  }
  
  #argList$HVG_list=unique(c(argList$HVG_list,argList$HVG_count))
  
  reRunCheck=F
  if(sum(names(argList)=="newRun")>0){
    if(argList$newRun){
      reRunCheck=T
    }
  }
  if(!saveFiles){
    reRunCheck=T
    if(length(argList$HVG_list)>1&(!benchmark_mode)){
      stop("Only one HVG count threshold can be specified")
    }
  }
  
  if(!reRunCheck){
    
    for(iHVG in argList$HVG_list){
      argList2=argList
      argList2$HVG_list=iHVG
      tmpCheck= tryCatch({load(.myFilePathMakerFn("pca_anno",argList = argList2,pseudoImportant = F));F}, error=function(e) {return(T)})
      if(!tmpCheck){
        argList$HVG_list=setdiff(argList$HVG_list,iHVG)
      }
    }
    
  }
  
  pca_final_res="Done"
  
  if(reRunCheck|length(argList$HVG_list)>0){
    tmpInd=c()
    if(sum(!is.null(argList$covariates))>0&!is.null(data$data_m)){
      for(i in argList$covariates[!is.null(argList$covariates)]){
        tmpInd=c(tmpInd,which(is.na(data$data_m@meta.data[,i])))
      }
      if(length(tmpInd)>0){
        data$data_m=data$data_m[,-unique(tmpInd)]
      }
    }
    
    if(argList$indScaling){
      dataList=data$data
      
      varFeatures=c()
      if(!is.null(argList$HVG_list)){
        for(iHVG in argList$HVG_list){
          varFeatures=c(varFeatures,data$varFeatures$Gene[which(data$varFeatures$Freq>(iHVG-1))])
        }
      } else {
        stop("HVG_list needs to be specified!")
      }
      
      varFeatures=unique(varFeatures)
      
      # browser()
      dataList2=parallel::mclapply(dataList,internal_scaleFn,varFeatures=varFeatures,mc.cores = argList$ncores)
      sl_ind=which(!unlist(lapply(dataList2,function(x) class(x)[1])) %in% c("array","matrix"))
      while(length(sl_ind)>0){
        gc()
        datalist3=parallel::mclapply(dataList[sl_ind],internal_scaleFn,varFeatures=varFeatures,mc.cores = 2)
        dataList2[sl_ind]=datalist3
        sl_ind=which(!unlist(lapply(dataList2,function(x) class(x)[1])) %in% c("array","matrix"))
        rm(datalist3)
      }
      dataList=dataList2
      rm(dataList2)
      
      if(!is.null(projection_data)){
        projectionDataList=parallel::mclapply(projection_data$data,internal_scaleFn,varFeatures=varFeatures,mc.cores = argList$ncores)
      }
      gc()
      
      resScaled=list()
      all_genes=unique(unlist(lapply(dataList,function(x) row.names(x))))
      for(i in 1:length(dataList)){
        tmp=dataList[[i]]
        tmp=tmp[match(all_genes,row.names(tmp)),,drop=F]
        if(sum(is.na(tmp))>0){
          tmp[is.na(tmp)]=0
        }
        
        row.names(tmp)=all_genes
        resScaled=c(resScaled,list(tmp))
        #dataList[[i]]=tmp
      }
      resScaled=do.call("cbind",resScaled)
      
      
      if(!is.null(projection_data)){
        projectionScaled=list()
        for(i in 1:length(projectionDataList)){
          tmp=projectionDataList[[i]]
          tmp=tmp[match(all_genes,row.names(tmp)),]
          tmp[is.na(tmp)]=0
          row.names(tmp)=all_genes
          projectionScaled=c(projectionScaled,list(tmp))
          
        }
        projectionScaled=do.call("cbind",projectionScaled)
        row.names(projectionScaled)=all_genes
      }
      
      pd=lapply(1:length(data$data),function(i){
        tmp=as.data.frame(colData(data$data[[i]]))
        tmp$sample=colnames(dataList[[i]])
        tmp
      })
      pd=do.call(eval(parse(text='plyr::rbind.fill')), pd)
      
      if(!is.null(projection_data)){
        pd_projection=lapply(1:length(projection_data$data),function(i){
          tmp=as.data.frame(colData(projection_data$data[[i]]))
          tmp$sample=colnames(projection_data$data[[i]])
          tmp
        })
        pd_projection=do.call(eval(parse(text='plyr::rbind.fill')), pd_projection)
        pd_all_projection=pd_projection
      }
      
      pd_all=pd
      gc()
      if(is.null(argList$HVG_list)){
        stop("HVG_list needs to be specified!")
      }
      
      if(!saveFiles){
        if(length(argList$HVG_list)>1){
          stop("HVG_list can't contain more than one value!")
        }
        
      }
      
      if(is.null(argList$input_highly_var_genes)){
        for(iHVG in argList$HVG_list){
          argList2=argList
          argList2$HVG_list=iHVG
          tmp_varFeatures=data$varFeatures$Gene[which(data$varFeatures$Freq>(iHVG-1))]
          
          object=myPrepDR(scaledData=resScaled,
                          features=tmp_varFeatures, verbose = TRUE)
          
          if(!is.null(projection_data)){
            projection_data=myPrepDR(scaledData=projectionScaled[row.names(object),],
                                     features=tmp_varFeatures, verbose = TRUE,projection_state = T)
          }
          
          projection_pca=NULL
          pca_res=.mySeuratRunPCA(object=object, npcs = max(50,argList2$nPCs+10),projection_data=projection_data)
          projection_pca=pca_res$projection_pca
          pca_res=pca_res$reduction.data
          
          pca_res=pca_res@cell.embeddings
          
          if(!is.null(projection_pca)){
            pca_res=rbind(pca_res,projection_pca)
            pd=plyr::rbind.fill(pd_all,pd_all_projection)
          } else {
            pd=pd_all
          }
          
          pd=pd[pd$sample %in% row.names(pca_res),]
          row.names(pd)=pd$sample
          pca_res=pca_res[match(row.names(pd),row.names(pca_res)),]
          if(saveFiles){
            save(pca_res,pd,file=.myFilePathMakerFn("pca_anno",argList=argList2,pseudoImportant = F))
          } else {
            pca_final_res=list(pd=pd,pca_res=pca_res)
          }
          
          gc()
        }
      } else {
        print("Highly variable genes are provided in argList, using it!")
        {
          
          tmp_varFeatures=as.character(argList$input_highly_var_genes)
          pca_res=Seurat:::RunPCA.default(object=myPrepDR(scaledData=resScaled,
                                                          features=tmp_varFeatures, verbose = TRUE,
                                                          ...),
                                          npcs = max(50,argList$nPCs+10))
          pca_res=pca_res@cell.embeddings
          
          pd=pd_all[pd_all$sample %in% row.names(pca_res),]
          row.names(pd)=pd$sample
          pca_res=pca_res[match(row.names(pd),row.names(pca_res)),]
          if(saveFiles){
            save(pca_res,pd,file=.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
          } else {
            pca_final_res=list(pd=pd,pca_res=pca_res)
          }
          
          gc()
        }
      }
    }
    
    if(!argList$indScaling){
      
      if(!is.null(projection_data)){
        stop("Projection is not yet implemented!")
      }
      
      dataList=data$data
      varFeatures=c()
      {
        if(is.null(argList$input_highly_var_genes)){
          if(!is.null(argList$HVG_list)){
            for(iHVG in argList$HVG_list){
              varFeatures=c(varFeatures,data$varFeatures$Gene[which(data$varFeatures$Freq>(iHVG-1))])
            }
          } else {
            stop("HVG_list needs to be specified!")
          }
          varFeatures=unique(varFeatures)
        } else {
          varFeatures=argList$input_highly_var_genes
        }
        
        
        resNorm=list()
        dataList=parallel::mclapply(dataList,internal_scaleFn,varFeatures=varFeatures,scaleData=F,mc.cores = argList$ncores)
        all_genes=unique(unlist(lapply(dataList,function(x) row.names(x))))
        for(i in 1:length(dataList)){
          tmp=as.matrix(dataList[[i]])
          tmp=tmp[match(all_genes,row.names(tmp)),,drop=F]
          
          if(sum(is.na(tmp))>0){
            if(ncol(tmp)>1){
              tmp[is.na(tmp)]=0
            } else {
              tmp[is.na(tmp)]=0
              tmp=matrix(tmp,ncol=1)
            }
            
          }
          
          row.names(tmp)=all_genes
          resNorm=c(resNorm,list(tmp))
          #dataList[[i]]=tmp
        }
        resNorm=do.call("cbind",resNorm)
        resScaled=Seurat:::ScaleData.default(as(resNorm,"dgCMatrix"),verbose=F)

        pd=lapply(1:length(dataList),function(i){
          tmp=as.data.frame(colData(data$data[[i]]))
          tmp$sample=colnames(dataList[[i]])
          tmp
        })
        pd=do.call(eval(parse(text='plyr::rbind.fill')), pd)
        
        
        pd_all=pd
        gc()
        if(is.null(argList$HVG_list)){
          stop("HVG_list needs to be specified")
        }
        
        if(!saveFiles & length(argList$HVG_list)>1){
          stop("HVG_list can't contain more than one value!")
        }
        
        if(is.null(argList$input_highly_var_genes)){
          for(iHVG in argList$HVG_list){
            argList2=argList
            argList2$HVG_list=iHVG
            tmp_varFeatures=data$varFeatures$Gene[which(data$varFeatures$Freq>(iHVG-1))]

            object=myPrepDR(scaledData=resScaled,
                            features=tmp_varFeatures, verbose = TRUE)
            
            if(!is.null(projection_data)){
              projection_data=myPrepDR(scaledData=projectionScaled[row.names(object),],
                                       features=tmp_varFeatures, verbose = TRUE,projection_state = T)
            }
            
            projection_pca=NULL
            pca_res=.mySeuratRunPCA(object=object, npcs = max(50,argList2$nPCs+10),projection_data=projection_data)
            projection_pca=pca_res$projection_pca
            pca_res=pca_res$reduction.data
          
            
            pca_res=pca_res@cell.embeddings
            
            if(!is.null(projection_pca)){
              pca_res=rbind(pca_res,projection_pca)
              pd=plyr::rbind.fill(pd_all,pd_all_projection)
            } else {
              pd=pd_all
            }
            
            pd=pd[pd$sample %in% row.names(pca_res),]
            row.names(pd)=pd$sample
            pca_res=pca_res[match(row.names(pd),row.names(pca_res)),]
            if(saveFiles){
              save(pca_res,pd,file=.myFilePathMakerFn("pca_anno",argList=argList2,pseudoImportant = F))
            } else {
              pca_final_res=list(pd=pd,pca_res=pca_res)
            }
            
            gc()
          }
        } else {
          print("Highly variable genes are provided in argList, using it!")
          {
            
            tmp_varFeatures=as.character(argList$input_highly_var_genes)
          
            pca_res=Seurat:::RunPCA.default(object=myPrepDR(scaledData=resScaled,
                                                            features=tmp_varFeatures, verbose = TRUE,
                                                            ...),
                                            npcs = max(50,argList$nPCs+10))
            pca_res=pca_res@cell.embeddings
            
            pd=pd_all[pd_all$sample %in% row.names(pca_res),]
            row.names(pd)=pd$sample
            pca_res=pca_res[match(row.names(pd),row.names(pca_res)),]
            if(saveFiles){
              argList$HVG_list=(0)
              save(pca_res,pd,file=.myFilePathMakerFn("pca_anno",argList=argList,pseudoImportant = F))
            } else {
              pca_final_res=list(pd=pd,pca_res=pca_res)
            }    
            gc()
          }
        }
      }
    }
  }
  return(pca_final_res)
}