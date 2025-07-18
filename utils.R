suppressMessages(suppressWarnings(library(biomaRt)))
suppressMessages(suppressWarnings(library(SingleCellExperiment)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(Matrix)))
suppressMessages(suppressWarnings(library(sva)))

convert_ensembl_to_symbol <- function(ensembl_ids) {
  # Connect to GRCh38 Ensembl
  mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = "108")

  # Query mapping
  mapping <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = ensembl_ids,
    mart = mart
  )

  mapping_unique <- mapping %>%
    distinct(ensembl_gene_id, .keep_all = TRUE)

  result <- data.frame(ensembl_gene_id = ensembl_ids, stringsAsFactors = FALSE) %>%
    left_join(mapping_unique, by = "ensembl_gene_id") %>%
    mutate(
      hgnc_symbol = ifelse(is.na(hgnc_symbol) | hgnc_symbol == "", ensembl_gene_id, hgnc_symbol)
    )


  # Find duplicated symbols
  dupes <- result$hgnc_symbol[duplicated(result$hgnc_symbol)]

  result <- result %>%
    mutate(
      final_symbol = ifelse(hgnc_symbol %in% dupes, paste0(hgnc_symbol, "_", ensembl_gene_id), hgnc_symbol)
    )

  return(result$final_symbol)
}

print_markers = function(df, clusters=NULL, logFC_thresh=0.5, sig_thresh=0.05){
    df = df[df$p_val_adj < sig_thresh & abs(df$avg_log2FC) >= logFC_thresh, ]
    if(!is.null(clusters)){
        df = df[df$cluster %in% clusters, ]
    }
    print(df[,c("cluster", "gene", "avg_log2FC", "pct.1", "pct.2")], row.names=F)
}

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

build_weighted_score = function(
    sobj,
    weight_df,
    md_score_col_name,
    md_norm_col = "log10_nUmi",
    weight_df_weight_col = "avg_log2FC",
    weight_df_gene_col = "gene",
    assay="SCT",
    scale_within=NULL) {

    # generate a score for each cell based on the weighted sum of gene counts
    # OUTPUTS
    # the input seurat object with new columns <md_score_col_name>, <md_score_col_name>__z, 
    # <md_score_col_name>__over__<md_norm_col>, and <md_score_col_name>__over__<md_norm_col>__z 
    

    # INPUTS
    # sobj: Seurat object
    # weight_df: data.frame with genes (in column `weight_df_gene_col`) and weights (in column `weight_df_weight_col`)
        # todo: default to weight of 1 for weight_df_weight_col = NULL
    # md_score_col_name: name of the column to be created  
    # md_norm_col: name of the column to divide the score by. Commonly nUmi or log10_nUmi
    # weight_df_weight_col: name of the column in weight_df that contains the weights
    # weight_df_gene_col: name of the column in weight_df that contains the genes
    # assay: name of the assay in sobj to use for counts
    # scale_within: if provided, z scores are calculated within the <scale_within> group within the sobj meta.data 
        # otherwise the z scores are calculated across all cells in `sobj`

    genes_orig = weight_df[[weight_df_gene_col]]
    weight_df = weight_df[genes_orig %in% rownames(sobj@assays[[assay]]@counts), ]

    weights = weight_df[[weight_df_weight_col]]
    genes = weight_df[[weight_df_gene_col]]



    filtered_counts = sobj[genes, ]

    if (!all(rownames(filtered_counts) == genes)){
        stop("Genes in weight_df out of order with respect to Counts")
    }
    
    md = sobj@meta.data 

    md[[md_score_col_name]] = colSums(filtered_counts@assays[[assay]]@counts * weights)
    if (is.null(scale_within)){
        md[[paste0(md_score_col_name, "__z")]] = scale(md[[md_score_col_name]])
    } else {
        md[[paste0(md_score_col_name, "__z")]] = ave(md[[md_score_col_name]], md[[scale_within]], FUN=scale)
    }

     if (!is.null(md_norm_col)){
         md[[paste0(md_score_col_name, "__over__", md_norm_col)]] = md[[md_score_col_name]] / md[[md_norm_col]]
         if (is.null(scale_within)){
             md[[paste0(md_score_col_name, "__over__", md_norm_col, "__z")]] = scale(md[[paste0(md_score_col_name, "__over__", md_norm_col)]])
         } else {
             md[[paste0(md_score_col_name, "__over__", md_norm_col, "__z")]] = ave(md[[paste0(md_score_col_name, "__over__", md_norm_col)]], md[[scale_within]], FUN=scale)
         }
     }

    sobj@meta.data = md

    return(sobj)
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

# get list of genes occurring in at least half of donor_sorts within a sort.
getCommonStrings <- function(list_of_lists, n) {
  # Create an empty list to store string occurrences
  string_occurrences <- list()

  # Iterate over each sublist
  for (sublist in list_of_lists) {
    # For each string in the sublist, increment its count in string_occurrences
    for (string in unique(sublist)) {
      if (!is.null(string_occurrences[[string]])) {
        string_occurrences[[string]] <- string_occurrences[[string]] + 1
      } else {
        string_occurrences[[string]] <- 1
      }
    }
  }

  # Filter the strings that occur at least n times
  common_strings <- names(string_occurrences)[string_occurrences >= n]

  return(common_strings)
}

getCountProportionDF = function(sobj, cat_col, type_col){
    # sobj: a seurat object with meta.data columns `cat_col` and `type_col` (both input as strings)
    # the counts and proportions of each `type_col` WITHIN each `cat_col` are returned as a dataframe
    # (e.g. if cat_col is donor_id and type_col is cell class, 
    # you'd get a long df of cell class counts and proportions wtihin each donor_id)

    md = sobj@meta.data
    df = (
        as.data.frame(table(md[[cat_col]], md[[type_col]])) 
        %>% group_by(Var1) %>% dplyr::mutate(total_count = sum(Freq), proportion = Freq / total_count)
    )
    colnames(df) = c(cat_col, type_col, 'count', 'total_count', 'proportion')
    return(df)
}

get_df_with_svs = function(edata, df, cols, ctr_cols=NULL, n=10){

    # edata: expression data as Matrix, genes in rows, samples in columns
    # df: dataframe with columns `cols` and `ctr_cols`, samples in rows
    # cols: list of columns in df to use as covariates
    # ctr_cols: list of columns in df to use as covariates
    # n: number of svs to calculate

    # returns a dataframe with the svs added as new columns

    model_formula = as.formula(paste0("~", paste(cols, collapse = " + ")))

    if (!is.null(ctr_cols)){
        ctr_formula = as.formula(paste0("~", paste(ctr_cols, collapse = " + ")))
    } else {
        ctr_formula = as.formula("~ 1")
    }

    cat(paste("\nRunning SVA with Formula:", deparse(model_formula)))
    cat(paste("\nAnd with Null Formula:", deparse(ctr_formula)))

    # Check if the columns in df are present in edata
    if (!all(unique(c(cols, ctr_cols)) %in% colnames(df))) {
        mensaje = paste("Error in get_df_with_svs: The following columns are not present in the dataframe:\n", 
                          paste(setdiff(unique(c(cols, ctr_cols)), colnames(df)), collapse = "\n"))
        stop(mensaje)
    }

    # throw an error if any of the sva columns contain NA values
    if (any (is.na(df[unique(c(cols, ctr_cols))]))) {
      mensaje = paste("Error in get_df_with_svs: The following columns contain NA values:\n", 
                      paste(unique(c(cols, ctr_cols))[colSums(is.na(df[unique(c(cols, ctr_cols))])) > 0], collapse = "\n"))
      stop(mensaje)
    }

    mod = model.matrix(model_formula, data = df)
    mod0 = model.matrix(ctr_formula, data = df)
    svobj = sva(edata, mod, mod0, n.sv = n)
    sv_factors = as.data.frame(svobj$sv)
    colnames(sv_factors) = paste0("SV", seq_len(ncol(sv_factors)))
    df_with_svs = cbind(df, sv_factors)

    return(df_with_svs)
}

join_df_to_sobj_metadata = function(sobj, df, metadata_join_cols, df_join_cols){
    # join a dataframe into seurat object metadata
    md = sobj@meta.data

    # save the rownames in the column or it will be lost in the merge
    md$rownames = rownames(md)

    # left join so that you don't lose any rows of meta.data
    md_new = merge(md, df, by.x=metadata_join_cols, by.y=df_join_cols, all.x=TRUE)

    # reorder new meta.data df 
    md_new = md_new[match(md$rownames, md_new$rownames),]
    rownames(md_new) = md_new$rownames

    if (! all(rownames(md_new) == rownames(md))){
        stop("Row names of new metadata do not match original metadata row names.")
    }

    md_new$rownames=NULL
    sobj@meta.data = md_new
    return(sobj)
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

getFracExpSeurat = function(
    seurat_obj,
    gene_list,
    initial_grouping_col = 'donor_id',
    final_grouping_col = 'anno'
){
    # given a seurat object and a set of genes,
    # determines the average fraction of cells in each `anno` that express each gene

    # currently averages over cells in donors, then anno
    # TODO: make grouping more flexible

    # get a col for each gene, row for each cell
    bin_counts_df = as.data.frame(t(seurat_obj@assays$RNA@counts[gene_list,])>1)

    # add grouping cols
    bin_counts_df[[initial_grouping_col]] = seurat_obj@meta.data[[initial_grouping_col]]
    bin_counts_df[[final_grouping_col]] = seurat_obj@meta.data[[final_grouping_col]]


    out_df = (bin_counts_df %>%
        group_by(!!rlang::sym(initial_grouping_col), !!rlang::sym(final_grouping_col)) %>% 
        summarize_all(mean) %>% 
        group_by(!!rlang::sym(final_grouping_col)) %>% 
        summarize_all(mean) %>% 
        as.data.frame()
    )

    rownames(out_df) = out_df[[final_grouping_col]]
    out_df = out_df[, !colnames(out_df) %in% c(initial_grouping_col, final_grouping_col)]

    return(out_df)
}

getMeanExpSeurat = function(
    seurat_obj,
    gene_list,
    slot='scale.data',
    scale_split='donor_id',
    initial_grouping_col = 'donor_id',
    final_grouping_col = 'anno',
    do_normalize=FALSE,
    do_scale=FALSE
){
    # given a seurat object and a set of genes,
    # determines the mean normalized expression within each `anno`

    # currently averages over cells in donors, then anno
    # TODO: make grouping more flexible

    if (do_normalize){
        seurat_obj = NormalizeData(seurat_obj)
    }
    if (do_scale){
        seurat_obj = ScaleData(seurat_obj, features=gene_list, split.by=scale_split)
    }

    rna = slot(seurat_obj@assays$RNA, slot)

    # get a col for each gene, row for each cell
    df = as.data.frame(t(rna[gene_list,]))
    df[[initial_grouping_col]] = seurat_obj@meta.data[[initial_grouping_col]]
    df[[final_grouping_col]] = seurat_obj@meta.data[[final_grouping_col]]

    out_df = (df %>%
        group_by(!!rlang::sym(initial_grouping_col), !!rlang::sym(final_grouping_col)) %>% 
        summarize_all(mean) %>% 
        group_by(!!rlang::sym(final_grouping_col)) %>% 
        summarize_all(mean) %>% 
        as.data.frame()
    )

    rownames(out_df) = out_df[[final_grouping_col]]
    out_df = out_df[, !colnames(out_df) %in% c(initial_grouping_col, final_grouping_col)]

    return(out_df)

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

run_fisher_exact_test = function(list_A, list_B, total_genes, alternative = "greater") {
    # Calculate the overlap and non-overlap
    a = length(intersect(list_A, list_B))
    b = length(setdiff(list_A, list_B))
    c = length(setdiff(list_B, list_A))
    d = total_genes - (a + b + c)

    # Create the contingency table
    contingency_table = matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    rownames(contingency_table) = c("In_List_B", "Not_in_List_B")
    colnames(contingency_table) = c("In_List_A", "Not_in_List_A")

    # Perform Fisher's exact test
    fisher_test_result = fisher.test(contingency_table, alternative = alternative)

    # Display the results
    return(fisher_test_result)
}


get_var_explained = function(sobj, assay = "RNA") {
    mat <- Seurat::GetAssayData(sobj, assay = assay, slot = "scale.data")
    pca <- sobj[["pca"]]
    total_variance <- sum(matrixStats::rowVars(mat))
    eigValues = (pca@stdev)^2  ## EigenValues
    varExplained = eigValues / total_variance
    return(varExplained)
}

broken_stick <- function(sobj, assay = "RNA") {
# Function to implement Broken Stick Model

    mat <- Seurat::GetAssayData(sobj, assay = assay, slot = "scale.data")
    pca <- sobj[["pca"]]

    # Get the total variance:
    total_variance <- sum(matrixStats::rowVars(mat))

    eigValues = (pca@stdev)^2  ## EigenValues
    varExplained = eigValues / total_variance

    # Extract eigenvalues (variance explained by each PC)
    n <- length(eigValues)  # Number of PCs

    # Compute Broken Stick expected values
    broken_stick_values <- sapply(1:n, function(k) sum(1 / (k:n)))

    # Create dataframe for plotting
    df <- data.frame(
        PC = 1:n,
        Variance_Explained = eigValues / sum(eigValues),  # Normalize eigenvalues
        Broken_Stick = broken_stick_values / sum(broken_stick_values)  # Normalize broken stick values
    )

    # Plot actual variance vs. broken stick model
    ggplot(df, aes(x = PC)) +
    geom_line(aes(y = Variance_Explained, color = "Variance Explained"), size = 1.2) +
    geom_line(aes(y = Broken_Stick, color = "Broken Stick Model"), size = 1.2, linetype = "dashed") +
    geom_point(aes(y = Variance_Explained), size = 2) +
    geom_point(aes(y = Broken_Stick), size = 2, shape = 17) +
    labs(title = "Broken Stick Model vs. PCA Variance Explained",
            y = "Proportion of Variance",
            x = "Principal Component") +
    scale_color_manual(values = c("Variance Explained" = "blue", "Broken Stick Model" = "red")) +
    theme_minimal()
}

softmax_cols = function(df, cols){
    # df: dataframe with columns `cols` to be transformed
    # cols: list of columns in df to be transformed
    # returns a dataframe with the softmax of the columns in `cols` added as new columns
    
    z_cols = paste0(cols, "__z")
    softmax_cols = paste0(cols, "__softmax")
    
    # Z-score normalization of the selected columns
    df[z_cols] = scale(df[cols])
    
    # Apply softmax row-wise: exp of z-scores divided by row sum of exp(z-scores)
    df[softmax_cols] =  t(apply(df[z_cols], 1, function(x) exp(x) / sum(exp(x))))
    
    return(df)
}


TopBottomGenesPerPC <- function(seurat_obj, n = 10, reduction = "pca") {
  # Ensure the reduction exists
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop(paste("Reduction", reduction, "not found in the Seurat object"))
  }
  
  # Extract loadings from the selected dimensionality reduction
  loadings <- seurat_obj[[reduction]]@feature.loadings
  
  # Iterate through each PC
  for (pc in colnames(loadings)) {
    cat("\n###", pc, "###\n")
    
    # Sort genes by loading values
    sorted_genes <- sort(loadings[, pc], decreasing = TRUE)
    
    # Get top and bottom N genes
    top_genes <- head(names(sorted_genes), n)
    bottom_genes <- tail(names(sorted_genes), n)
    
    # Print results
    cat("Top", n, "genes:\n", paste(top_genes, collapse = ", "), "\n")
    cat("Bottom", n, "genes:\n", paste(bottom_genes, collapse = ", "), "\n")
  }
}

prop = function(vector){
  if(length(vector) == 0){
    stop("Vector is empty")
  }
  return(sum(vector) / length(vector))
}


load_obj = function(path){
  if(grepl("\\.rds$", path)){
    obj = readRDS(path)
  } else if(grepl("\\.qs$", path)){
    obj = qread(path)
  } else {
    stop("Path must end with .rds or .qs")
  }
  return(obj)
}

save_obj = function(obj, path){
  # if the path ends with ".rds", save as rds. If ".qs", save as qs
  # otherwise throw an error
  if(!dir.exists(dirname(path))){
    dir.create(dirname(path), recursive = TRUE)
  }

  if(grepl("\\.rds$", path)){
    saveRDS(obj, path)
  } else if(grepl("\\.qs$", path)){
    qsave(obj, path)
  } else {
    stop("Path must end with .rds or .qs")
  }
}