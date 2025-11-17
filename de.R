library(glue)
g = glue::glue
source("~/sc-online/utils.R")
source("~/sc-online/extraGenes.R")

pr_from_df <- function(df, score_col, higher_score_is_positive = TRUE, title=NULL) {
  # df: data frame with columns is_tp, is_fp, is_tn, is_fn, and score_col
  

  if (!is.character(score_col)) score_col <- rlang::as_name(rlang::ensym(score_col))
  stopifnot(all(c("is_tp","is_fp","is_tn","is_fn", score_col) %in% names(df)))
  stopifnot(is.numeric(df[[score_col]]))

  # ground truth independent of threshold
  truth <- (df$is_tp | df$is_fn)
  scores <- df[[score_col]]
  keep <- !(is.na(truth) | is.na(scores))
  truth  <- truth[keep]
  scores <- scores[keep]

  if (!higher_score_is_positive) scores <- -scores

  pos_scores <- scores[truth]
  neg_scores <- scores[!truth]

  # PR curve and AUPRC
  pr <- PRROC::pr.curve(scores.class0 = pos_scores,
                         scores.class1 = neg_scores,
                         curve = TRUE)

  # curve columns: recall (x), precision (y), threshold
  pr_df <- data.frame(Recall = pr$curve[,1],
                      Precision = pr$curve[,2],
                      Threshold = pr$curve[,3])

  prevalence <- mean(truth)

  p <- ggplot(pr_df, aes(Recall, Precision)) +
    geom_path(size = 1) +
    geom_hline(yintercept = prevalence, linetype = "dashed") +
    coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
    theme_minimal() +
    labs(
      title = g("{ifelse(is.null(title), '', title) }\nPrecision-Recall AUPRC = {round(pr$auc.integral, 3)}"),
      subtitle = g("Score: {score_col}  |  Prevalence = {round(prevalence, 3)}"),
      x = "Recall", y = "Precision"
    )+ 
    coord_cartesian(clip = 'off') +

    theme(
        plot.title = element_text(size=16),
        plot.subtitle = element_text(size=12),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        axis.title.y = element_text(margin=margin(r=15)),
        plot.title.position = "plot"
    )

  roc = PRROC::roc.curve(scores.class0 = pos_scores,
                        scores.class1 = neg_scores,
                        curve = TRUE)
  
  roc_df = data.frame(FPR = roc$curve[,1],
                      TPR = roc$curve[,2],
                      Threshold = roc$curve[,3])
  p_roc <- ggplot(roc_df, aes(x=FPR, y=TPR)) +
    geom_path(size = 1) +
    geom_abline(slope=1, intercept=0, linetype = "dashed") +
    coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
    theme_minimal() +
    labs(
      title = g("{ifelse(is.null(title), '', title) }\nROC-AUC = {round(roc$auc, 3)}"),
      subtitle = g("Score: {score_col}  |  Baseline (random) = 0.5"),
      x = "False Positive Rate", y = "True Positive Rate"
    ) + 
    theme(
        plot.title = element_text( size=16),
        plot.subtitle = element_text(size=12),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        axis.title.y = element_text(margin=margin(r=15)),
        plot.title.position = "plot"
    )

  list(auprc = unname(pr$auc.integral), pr_curve = pr_df, pr_plot = p, auroc = unname(roc$auc), roc_curve = roc_df, roc_plot = p_roc)
}


.extra_sconline.duplicateCorrelation=function(
    object, 
    design = NULL, 
    ndups = 2, 
    spacing = 1, 
    block = NULL, 
    trim = 0.15,
    weights = NULL){

  require(limma)
  y = limma:::getEAWP(object)
  print(dim(y$design))
  M = y$exprs
  ngenes = nrow(M)
  narrays = ncol(M)
  if (is.null(design)) 
    design = y$design
  if (is.null(design)) 
    design = matrix(1, ncol(y$exprs), 1)
  else {
    design = as.matrix(design)
    if (mode(design) != "numeric") 
      stop("design must be a numeric matrix")
  }
  
  if (nrow(design) != narrays)
    stop("Number of rows of design matrix does not match number of arrays")
  ne = limma:::nonEstimable(design)
  if (!is.null(ne)) 
    cat("Coefficients not estimable:", paste(ne, collapse = " "), 
        "\n")
  nbeta = ncol(design)
  if (missing(ndups) && !is.null(y$printer$ndups)) 
    ndups = y$printer$ndups
  if (missing(spacing) && !is.null(y$printer$spacing)) 
    spacing = y$printer$spacing
  if (missing(weights) && !is.null(y$weights)) 
    weights = y$weights
  if (!is.null(weights)) {
    weights = asMatrixWeights(weights, dim(M))
    weights[weights <= 0] = NA
    M[!is.finite(weights)] = NA
  }
  if (is.null(block)) {
    if (ndups < 2) {
      warning("No duplicates: correlation between duplicates not estimable")
      return(list(cor = NA, cor.genes = rep(NA, nrow(M))))
    }
    if (is.character(spacing)) {
      if (spacing == "columns") 
        spacing = 1
      if (spacing == "rows") 
        spacing = object$printer$nspot.c
      if (spacing == "topbottom") 
        spacing = nrow(M)/2
    }
    Array = rep(1:narrays, rep(ndups, narrays))
  }
  else {
    ndups = 1
    nspacing = 1
    Array = block
  }
  if (is.null(block)) {
    M = limma:::unwrapdups(M, ndups = ndups, spacing = spacing)
    ngenes = nrow(M)
    if (!is.null(weights)) 
      weights = limma:::unwrapdups(weights, ndups = ndups, spacing = spacing)
    design = design %x% rep(1, ndups)
  }
  if (!requireNamespace("statmod", quietly = TRUE)) 
    stop("statmod package required but is not installed")
  rho = rep(NA, ngenes)
  nafun = function(e) NA
  for (i in 1:ngenes) {
    y = drop(M[i, ])
    o = is.finite(y)
    A = factor(Array[o])
    nobs = sum(o)
    nblocks = length(levels(A))
    if (nobs > (nbeta + 2) && nblocks > 1 && nblocks < nobs - 1) {
      y = y[o]
      X = design[o, , drop = FALSE]
      Z = model.matrix(~0 + A)
      if (!is.null(weights)) {
        w = drop(weights[i, ])[o]
        s = tryCatch(statmod::mixedModel2Fit(y, X, Z, 
                                              w, only.varcomp = TRUE, maxit = 20)$varcomp, error = nafun)
      } else{
        s = tryCatch(statmod::mixedModel2Fit(y, X, 
                                              Z, only.varcomp = TRUE, maxit = 20)$varcomp, 
                      error = nafun)
      } 
      if (!is.na(s[1])) 
        rho[i] = s[2]/sum(s)
    }
  }
  arho = atanh(pmax(-1, rho))
  mrho = tanh(mean(arho, trim = trim, na.rm = TRUE))
  list(consensus.correlation = mrho, cor = mrho, atanh.correlations = arho,value.list=arho)
}


.extra_sconline.Fit_LimmaDreamFn=function(
    sl_data,
    pd,
    model,
    random_effect=NULL,
    TMMnorm=F,
    VSTnorm=F,
    prior.count=1,
    quantile.norm=F,
    bkg_genes=NULL,
    no_normalization=F,
    dc.object=NULL,
    include.malat1.as.covariate=F,
    ncores=4){
  require(edgeR)
  require(limma)
  require(variancePartition)
  
  param = SnowParam(ncores, "SOCK", progressbar=TRUE)
  
  # estimate weights using linear mixed model of dream
  
  sl_data=sl_data[rowSds(as.matrix(counts(sl_data))) > 0, ]
  
  dge = DGEList(counts=counts(sl_data))
  if(is.null(bkg_genes)){
    tmpCount=apply(counts(sl_data),1,function(x) sum(x>0))
    keep=tmpCount>10
    
  } else {
    keep=row.names(sl_data) %in% bkg_genes
  }
  
  dge = dge[keep,,keep.lib.sizes=FALSE]
  
  if(TMMnorm){
    dge = calcNormFactors(dge)
  }
  
  # estimate weights using linear mixed model of dream
  #model=as.formula("~status + anno_sex +nUMI_scaled  + pseudocell_size_scale +(1 | anno_batch)")
  
  if(quantile.norm){
    vobjDream = voomWithDreamWeights( dge, model,as.data.frame(pd), BPPARAM=param, normalize.method="quantile" )
  } else {
    vobjDream = voomWithDreamWeights( dge, model,as.data.frame(pd), BPPARAM=param )
  }
  
  
  
  
  fit = dream( vobjDream, model, as.data.frame(pd) )
  
  
  return(list(fit=fit,dc=NULL,model=model,normData=NULL,blocked_analysis=NULL))
}


.extra_sconline.Fit_LimmaTrendFn=function(
    sl_data,
    model,
    random_effect=NULL,
    TMMnorm=F,
    VSTnorm=F,
    prior.count=1,
    quantile.norm=F,
    bkg_genes=NULL,
    no_normalization=F,
    dc.object=NULL,
    include.malat1.as.covariate=F){
  require(edgeR)
  require(limma)
  library(Matrix)
  
  if(VSTnorm){
    logCPM=.myRNAseqNormVSTfn(inputCountData=sl_data,fitType="parametric")
    logCPM=logCPM$originalData
    logCPM=counts(logCPM)
    if(!is.null(bkg_genes)){
      keep=row.names(sl_data) %in% bkg_genes
    } else {
      keep = rowSums(logCPM>3)>(0.05*ncol(logCPM))
    }
    
    logCPM = logCPM[keep,]
  } else if(no_normalization){
    logCPM=counts(sl_data)
    if(quantile.norm){
      logCPM=limma::normalizeQuantiles(logCPM)
    }
  } else {
    if(!is.null(bkg_genes)){
      keep=row.names(sl_data) %in% bkg_genes
    } else {
      tmpCount2=apply(counts(sl_data),1,function(x) sum(x>0))
      tmpCount=rowSums(edgeR::cpm(as.matrix(counts(sl_data))))
      keep=tmpCount>max(0.01*ncol(sl_data),min(15,ncol(sl_data)/3))
      keep=keep & tmpCount2>max(0.01*ncol(sl_data),min(10,ncol(sl_data)/3))
    }
    
    #tmpCount=apply(counts(sl_data),1,function(x) sum(x>0))
    #keep=rowData(sl_data)$exp.pct>0.03
    #keep=filterByExpr(counts(sl_data),min.count = 0.9, 
    #                  min.total.count = max(0.02*ncol(sl_data),10))
    print(paste("Number of expressed genes:",sum(keep)))
    tmpexp=counts(sl_data)[keep,]
    dge = DGEList(tmpexp)
    
    if(TMMnorm){
      dge = calcNormFactors(dge)
    }
    
    logCPM = new("EList")
    logCPM$E = edgeR::cpm(dge,normalized.lib.sizes = TMMnorm, log = TRUE, prior.count = prior.count)
    if(quantile.norm){
      logCPM$E=limma::normalizeQuantiles(logCPM$E)
    }
  }
  
  
  if(include.malat1.as.covariate){
    malat_gene=row.names(logCPM)[grepl("malat1",tolower(row.names(logCPM)))]
    if(class(logCPM)[1]=="EList"){
      model=cbind(model,malat=as.numeric(logCPM$E[malat_gene,]))
    } else {
      model=cbind(model,malat=as.numeric(logCPM[malat_gene,]))
    }
    
  }
  
  
  dc=NULL
  
  if(!is.null(random_effect)){
    if(is.null(dc.object)){

        # MAJOR TODO! WHY SHOULD THERE EVER BE A MISMATCH BETWEEN DESIGN AND LOGCPM?!
      dc = .extra_sconline.duplicateCorrelation(logCPM,design=model, block=colData(sl_data)[,random_effect])
    } else {
      dc=dc.object
    }
    
  }
  
  #
  blocked_analysis=F
  dc$consensus.correlation[dc$consensus.correlation<0.001]=0.001
  if(!is.null(dc)){
    if(!is.nan(dc$consensus.correlation)){
      if(abs(dc$consensus.correlation)<0.9){
        fit = lmFit(logCPM, model,block = as.character(colData(sl_data)[,random_effect]), correlation=dc$consensus.correlation)
        blocked_analysis=T
      } else {
        fit = lmFit(logCPM, model)
      }
    } else {
      fit = lmFit(logCPM, model)
    }
  } else {
    fit = lmFit(logCPM, model)
  }
  
  
  return(list(fit=fit,dc=dc,model=model,normData=logCPM,blocked_analysis=blocked_analysis))
}


.extra_sconline.Fit_LimmaVoomFn=function(
    sl_data,
    model,
    random_effect=NULL,
    quantile.norm=F,
    sample.weights=F,
    TMMnorm=F,
    bkg_genes=NULL,
    dc.object=NULL){
  require(edgeR)
  require(limma)
  
  
  dge = DGEList(counts=counts(sl_data))
  if(is.null(bkg_genes)){
    tmpCount=apply(counts(sl_data),1,function(x) sum(x>0))
    keep=tmpCount>10
    
  } else {
    keep=row.names(sl_data) %in% bkg_genes
  }
  
  dge = dge[keep,,keep.lib.sizes=FALSE]
  
  if(TMMnorm){
    dge = calcNormFactors(dge)
  }
  
  if(sample.weights){
    if(quantile.norm){
      logCPM = voomWithQualityWeights(dge, model, normalize.method="quantile",plot = T,save.plot=T)
    } else {
      logCPM = voomWithQualityWeights(dge, model,plot = T,save.plot=T)
    }
  } else {
    if(quantile.norm){
      logCPM = voom(dge, model, normalize.method="quantile",plot = T,save.plot=T)
    } else {
      logCPM = voom(dge, model,plot = T,save.plot=T)
    }
  }
  
  dc=NULL
  if(!is.null(random_effect)){
    if(is.null(dc.object)){
      dc = .extra_sconline.duplicateCorrelation(logCPM,design=model, block=colData(sl_data)[,random_effect])
    } else {
      dc=dc.object
    }
  }
  
  
  blocked_analysis=F
  if(!is.null(dc)){
    if(!is.nan(dc$consensus.correlation)){
      if(abs(dc$consensus.correlation)<0.9){
        fit = lmFit(logCPM, model,block = colData(sl_data)[,random_effect], correlation=dc$consensus.correlation)
        blocked_analysis=T
      } else {
        fit = lmFit(logCPM, model)
      }
    } else {
      fit = lmFit(logCPM, model)
    }
  } else {
    fit = lmFit(logCPM, model)
  }
  
  
  #logCPM=preprocessCore::normalize.quantiles(logCPM)
  
  
  return(list(fit=fit,dc=dc,model=model,blocked_analysis=blocked_analysis,normData=logCPM))
}


.sconline.fitLimmaFn=function(
    inputExpData,
    covariates,
    randomEffect,
    DEmethod="Trend",
    normalization="CPM",
    quantile.norm=F,
    bkg_genes=NULL,
    VST_fitType="parametric",
    prior.count=1,
    include.malat1.as.covariate=F,
    dc.object=NULL,
    dream_ncores=4){
  
  #for covariates use the scaled versions of nUMI (), nGene(QC_Gene_unique_count_scale), and pseudocell_size as appropriate
  #include cellType as covariate as appropriate
  
  #DEmethod="Trend";normalization="CPM";quantile.norm=F;VST_fitType="parametric";prior.count=1
  
  normalization=match.arg(normalization,c("CPM", "TMM", "VST", "rmTop50","none"))
  DEmethod=match.arg(DEmethod,c("Trend", "Voom", "VoomSampleWeights","Dream"))
  
  
  inputExpData=inputExpData[rowSums(counts(inputExpData))>5,]
  inputExpData=inputExpData[apply(counts(inputExpData),1,sd)>0.0001,]
  
  if(DEmethod=="Voom"&normalization=="VST"){
    stop("normalization method can be either Voom or VST but not both")
  }
  
  norm.tmm=normalization=="TMM"
  norm.rmTop50=normalization=="rmTop50"
  norm.vst=normalization=="VST"
  no_normalization=normalization=="none"
  
  if(is.null(inputExpData)){
    stop("Expression data is missing!")
  } else if(tolower(class(inputExpData))=="seurat"){
    inputExpData=.sconline.convertSeuratToExpressionSet(object=inputExpData)
  } else if(class(inputExpData)!="SingleCellExperiment") {
    stop("Unrecognized inputExpression data!")
  }
  
  if(sum(colnames(colData(inputExpData)) %in% covariates)<length(covariates)){
    stop(paste0("Covariates ",paste(setdiff(covariates,colnames(colData(inputExpData))),collapse = ", ")," were not identified in the inputExpData!"))
  }
  
  
  
  if(DEmethod!="Dream"){
    model_matrix="~0"
    for(icov in covariates){
      if(length(unique(colData(inputExpData)[,icov]))>1){
        if(class(colData(inputExpData)[,icov])==class(factor)){
          colData(inputExpData)[,icov]=as.character(colData(inputExpData)[,icov])
        }
        model_matrix=paste0(model_matrix,"+",icov)
      } else {
        warning(paste0("Excluding ",icov," covariate as it has only one level!"))
      }
    }
    
    model_matrix = model.matrix(as.formula(model_matrix),data=colData(inputExpData))
  } else {
    model_matrix="~"
    for(icov in covariates){
      if(length(unique(colData(inputExpData)[,icov]))>1){
        if(class(colData(inputExpData)[,icov])==class(factor)){
          colData(inputExpData)[,icov]=as.character(colData(inputExpData)[,icov])
        }
        if(model_matrix=="~"){
          model_matrix=paste0(model_matrix,icov)
        } else {
          model_matrix=paste0(model_matrix,"+",icov)
        }
        
      } else {
        warning(paste0("Excluding ",icov," covariate as it has only one level!"))
      }
    }
    model_matrix = as.formula(paste0(model_matrix," + (1|",randomEffect,")"))
    print(model_matrix)
  }
  
  if(norm.rmTop50){
    if(sum(colnames(rowData(inputExpData))=="QC_top50_expressed")>0){
      inputExpData=inputExpData[which(rowData(inputExpData)$QC_top50_expressed=="No"),]
    } else {
      warning("Column QC_top50_expressed was not identified in the gene attribute dataframe. skipping the removal of top 50 expressed genes!")
    }
  }
  print(dim(inputExpData))
  print(dim(model_matrix))
  if (ncol(inputExpData) != nrow(model_matrix)){
    print("Dimension mismatch between inputExpData and model_matrix!")
    inputExpData=inputExpData[,colnames(inputExpData) %in% row.names(model_matrix)]
  }

  
  #sl_data=inputExpData;model=model_matrix;random_effect=randomEffect;quantile.norm = quantile.norm;
  #TMMnorm=norm.tmm;VSTnorm=norm.vst;prior.count=prior.count;pd=colData(inputExpData)
  switch(DEmethod,
         Trend=.extra_sconline.Fit_LimmaTrendFn(sl_data=inputExpData,model=model_matrix,random_effect=randomEffect,
            quantile.norm = quantile.norm,TMMnorm=norm.tmm,VSTnorm=norm.vst,prior.count=prior.count,
            bkg_genes=bkg_genes,no_normalization=no_normalization,dc.object=dc.object,include.malat1.as.covariate),
         Dream=.extra_sconline.Fit_LimmaDreamFn(sl_data=inputExpData,model=model_matrix,pd=as.data.frame(colData(inputExpData)),
            random_effect=randomEffect,quantile.norm = quantile.norm,TMMnorm=norm.tmm,VSTnorm=norm.vst,prior.count=prior.count,
            bkg_genes=bkg_genes,no_normalization=no_normalization,dc.object=dc.object,include.malat1.as.covariate,ncores = dream_ncores),
         Voom=.extra_sconline.Fit_LimmaVoomFn(sl_data=inputExpData,model=model_matrix,random_effect=randomEffect,
            quantile.norm=quantile.norm,sample.weights=F,TMMnorm=norm.tmm,bkg_genes=bkg_genes,dc.object=dc.object),
         VoomSampleWeights=.extra_sconline.Fit_LimmaVoomFn(sl_data=inputExpData,model=model_matrix,random_effect=randomEffect,
            quantile.norm=quantile.norm,sample.weights=T,TMMnorm=norm.tmm,bkg_genes=bkg_genes,dc.object=dc.object))
  
}

get_de_correspondence_mats = function(de_dfs, union_min=50, mode="union"){

  # de_dfs should be a named list of either
  # (a) paths to DE csvs 
  # or (b) dfs themselves
  # containing columns gene, logFC, and adj.P.Val 

  # outputs (for the union of significant genes):
  # correlation: matrix of spearman correlations between de results
  # p: matrix of p-values for the correlations
  # consistency: matrix of sign consistency between de results
  # n: matrix of number of genes in common between de results

  # todo: add option for all genes and sig intersect, not just union

  cor_mat = matrix(NA, nrow=length(de_dfs), ncol=length(de_dfs))
  p_mat = matrix(NA, nrow=length(de_dfs), ncol=length(de_dfs))
  cons_mat = matrix(NA, nrow=length(de_dfs), ncol=length(de_dfs))
  n_mat = matrix(NA, nrow=length(de_dfs), ncol=length(de_dfs))

  i=1
  for (name1 in names(de_dfs)){
      j=1
      for(name2 in names(de_dfs)){
          if (class(de_dfs[[name1]]) == "character") {
              df1 = read.csv(de_dfs[[name1]])
          } else {
              df1 = de_dfs[[name1]]
          }

          if (class(de_dfs[[name2]]) == "character") {
              df2 = read.csv(de_dfs[[name2]])
          } else {
              df2 = de_dfs[[name2]]
          }
          
          # get genes in common, and make sure they are in the same order
          genes = intersect(df1$gene, df2$gene)

          df1 = df1[df1$gene %in% genes,]
          df2 = df2[df2$gene %in% genes,]

          if (mode == "all"){
            
            df1_sig = df1 
            df2_sig = df2

            union_sig_genes = genes
          } else {
      
            df1_sig = df1[df1$adj.P.Val < 0.05,]
            df2_sig = df2[df2$adj.P.Val < 0.05,]
            if (mode == "union"){
              union_sig_genes = sort(union(df1_sig$gene, df2_sig$gene))
            } else if (mode == "intersect"){
              union_sig_genes = sort(intersect(df1_sig$gene, df2_sig$gene))
            }
            else {
            stop("mode must be one of all, union, or intersect")
            }
          }

          n = length(union_sig_genes)
          n_mat[i, j] = n

          if (length(union_sig_genes) < union_min) {
              
              if (i == j) {
                  cor_mat[i, j] = 1
                  p_mat[i, j] = "<0.000001"
                  cons_mat[i, j] = 1
              }

              j = j + 1

              next
          }

          df1_union = df1[df1$gene %in% union_sig_genes,]
          df2_union = df2[df2$gene %in% union_sig_genes,]

          df1_union = df1_union[match(union_sig_genes, df1_union$gene),]
          df2_union = df2_union[match(union_sig_genes, df2_union$gene),]

          # get the correlation, the p-value, the number of genes in common, and the sign consistency
          cor_test = cor.test(df1_union$logFC, df2_union$logFC, method="spearman")
          

          cor_mat[i, j] = cor_test$estimate
          pval = round(cor_test$p.value, 6)
          if (pval == 0) {
              pval = "<0.000001"
          }
          p_mat[i, j] = pval
          
          cons_mat[i, j] = sum(sign(df1_union$logFC) == sign(df2_union$logFC)) / n

          j = j + 1
      }
      i = i + 1
  }

  cor_mat = as.data.frame(cor_mat)
  p_mat = as.data.frame(p_mat)
  cons_mat = as.data.frame(cons_mat)
  n_mat = as.data.frame(n_mat)

  rownames(cor_mat) = names(de_dfs)
  colnames(cor_mat) = names(de_dfs)

  rownames(p_mat) = names(de_dfs)
  colnames(p_mat) = names(de_dfs)

  rownames(cons_mat) = names(de_dfs)
  colnames(cons_mat) = names(de_dfs)

  rownames(n_mat) = names(de_dfs)
  colnames(n_mat) = names(de_dfs)

  return(list(
    "correlation" = cor_mat,
    "p" = p_mat,
    "consistency" = cons_mat,
    "n" = n_mat
  ))
}

pseudobulk_seurat = function(
    sobj, 
    grouping_cols, 
    assay="RNA", 
    min_n_cells = 10, 
    min_counts_gene = 0, 
    min_frac_gene = 0, 
    contrast_col="case_control",
    cols_to_mean=NULL,
    cols_to_weighted_avg=NULL,
    cols_to_median=NULL,
    weight_col="nCount_RNA"){
    # sobj: Seurat object containing metadata columns 
    # (commonly nUMI or nCount_RNA, pct_intronic, pct_mito, etc.)

    create_df_with_contrast_col = function() {
        setNames(data.frame(factor(levels = sort(unique(md[[contrast_col]])))), contrast_col)
    }

    df = sobj@meta.data
    if (! all(grouping_cols %in% colnames(df))){
        stop(paste0(
            "Not all grouping columns are present in the Seurat object!\n",
            "Missing Columns: ", paste(grouping_cols[!grouping_cols %in% colnames(df)], collapse=", ")
        ))
    }

    # normalize selection vectors
    if (is.null(cols_to_mean)) cols_to_mean <- character(0)
    if (is.null(cols_to_weighted_avg)) cols_to_weighted_avg <- character(0)
    if (is.null(cols_to_median)) cols_to_median <- character(0)

    # stop if there are missing columns
    stop_missing <- function(cols, label){
        missing <- setdiff(cols, colnames(df))
        if (length(missing)) {
            stop(sprintf("%s not found in metadata: %s", label, paste(missing, collapse=", ")))
        }
    }
    stop_missing(cols_to_mean, "cols_to_mean")
    stop_missing(cols_to_weighted_avg, "cols_to_weighted_avg")
    stop_missing(cols_to_median, "cols_to_median")

    # guard against overlaps across lists (ambiguous summaries)
    ov1 = intersect(cols_to_mean, cols_to_weighted_avg)
    ov2 = intersect(cols_to_mean, cols_to_median)
    ov3 = intersect(cols_to_weighted_avg, cols_to_median)
    if (length(ov1) + length(ov2) + length(ov3) > 0){
        msg <- c(
            if(length(ov1)) sprintf("mean & weighted_avg: %s", paste(ov1, collapse=", ")),
            if(length(ov2)) sprintf("mean & median: %s", paste(ov2, collapse=", ")),
            if(length(ov3)) sprintf("weighted_avg & median: %s", paste(ov3, collapse=", "))
        )
        stop("Columns cannot appear in multiple summary lists:\n  - ", paste(msg, collapse="\n  - "))
    }

    # choose a per-cell weight column
    if (!weight_col %in% colnames(df)) {
        stop(sprintf("`weight_col` '%s' not found in metadata. Provide a valid column name.", weight_col))
    }

    # get the data type of each grouping column; build a stable string "grouping" to avoid numeric weirdness
    grouping_is_numeric = sapply(grouping_cols, function(col) {is.numeric(df[[col]])})
    grouping_cols_non_numeric = grouping_cols[!grouping_is_numeric]
    df$grouping = apply(df[grouping_cols_non_numeric, drop=FALSE], 1, function(row) paste(row, collapse = "_"))
    df = sanitize_col_values(df, "grouping")
    sobj$grouping = df$grouping

    cat(paste("Grouping By Columns: ", paste(grouping_cols, collapse=", "), "\n"))
    grouping_cols = c(grouping_cols, "grouping")

    # summarise per group. Ensure that n_cells and sum_nUMI are not in any of the summary lists
    cols_to_mean <- setdiff(cols_to_mean, c("n_cells","sum_nUMI"))
    cols_to_median <- setdiff(cols_to_median, c("n_cells","sum_nUMI"))
    cols_to_weighted_avg <- setdiff(cols_to_weighted_avg, c("n_cells","sum_nUMI"))


    df_bulk = df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping_cols))) %>%
    dplyr::summarise(
        n_cells = dplyr::n(),
        sum_nUMI = sum(nCount_RNA, na.rm = TRUE),
        dplyr::across( # simple means
            dplyr::all_of(cols_to_mean),
             ~ mean(.x, na.rm = TRUE), 
             .names = "{.col}"), 
        dplyr::across( # medians
            dplyr::all_of(cols_to_median), 
            ~ stats::median(.x, na.rm = TRUE), 
            .names = "{.col}"), 
        dplyr::across( # averages weighted by `weight_col`
            dplyr::all_of(cols_to_weighted_avg),
            ~ stats::weighted.mean(.x, w = .data[[weight_col]], na.rm = TRUE),
            .names = "{.col}"),
        .groups = "drop"
    ) %>%
    dplyr::mutate(log10_sum_nUMI = log10(sum_nUMI)) %>% 
    dplyr::mutate(log10_nUMI = log10(sum_nUMI / n_cells))

    # filter on min cells
    df_bulk = df_bulk[df_bulk$n_cells >= min_n_cells, ]
    if (nrow(df_bulk) == 0){
        return(list(counts=NULL,  metadata=create_df_with_contrast_col()))
    }

    # keep only cells in retained groups
    sobj = sobj[, df$grouping %in% df_bulk$grouping]

    # aggregate counts and apply gene filters
    counts_bulk = Seurat::AggregateExpression(sobj, group.by = "grouping")[[assay]]
    counts_orig = Seurat::GetAssayData(sobj, assay, slot="counts")
    counts_bulk = counts_bulk[
        rowSums(counts_orig) >= min_counts_gene &
        rowMeans(counts_orig > 0) >= min_frac_gene, 
    ]

    # sync metadata to count columns
    colnames(counts_bulk) = gsub("-", "_", colnames(counts_bulk))
    df_bulk = df_bulk[match(colnames(counts_bulk), df_bulk$grouping), ]

    if (nrow(df_bulk) == 0){
        return(list(counts=NULL, metadata=create_df_with_contrast_col()))
    }

    return(list(counts = counts_bulk, metadata = df_bulk))
}


build_robustness_df = function( 
    master_de_csv_path, jk_dir, cov_col, 
    padj_col = "adj.P.Val", sig_thresh=0.05, neg_log_padj_ceil=30) {
    
    # master_de_csv_path is the path to the master DE csv run against all the data
    # jk_dir is a directory that should contain a set of subdirectories, each containing JK DE csvs
    # cov_col is the covariate name to use for the robustness analysis

    master_df = read.csv(master_de_csv_path)
    master_df$significant = master_df[[padj_col]] < sig_thresh
    master_df$negative_log10_padj = -log10(master_df[[padj_col]] + 1e-30)
    master_df$negative_log10_padj[master_df$negative_log10_padj > neg_log_padj_ceil] = neg_log_padj_ceil

    cnames = colnames(master_df)
    cnames[!cnames == "gene"] = paste0(cnames[!cnames == "gene"], "__master")
    colnames(master_df) = cnames


    files=list.files(jk_dir, recursive=T, full.name=T, pattern=paste0("*", cov_col, ".csv"))
    de_dfs = lapply(files, function(x) {
        df = read.csv(x)
        df$significant = df$adj.P.Val < sig_thresh
        df$negative_log10_padj = -log10(df[[padj_col]] + 1e-30)
        df$negative_log10_padj[df$negative_log10_padj > neg_log_padj_ceil] = neg_log_padj_ceil

        bndn = basename(dirname(x))
        bndn_list = strsplit(bndn, "__")[[1]]
        jk = bndn_list[[length(bndn_list)]]
        cnames = colnames(df)
        cnames[!cnames == "gene"] = paste0(cnames[!cnames == "gene"], "__jk_", jk)
        colnames(df) = cnames
        return(df)
    })

    # full outer join all of these together on the gene column
    de_df = Reduce(function(x, y) merge(x, y, by = "gene", all = T), de_dfs)
    de_df = merge(master_df, de_df, by = "gene", all = T)

    # Identify the logFC columns (excluding master)
    jk_cols = grep("^logFC__jk_", names(df), value = TRUE)

    de_df = de_df %>%
        rowwise() %>%
        mutate(
            robustness = {
                master_sign = sign(logFC__master)
                jk_signs = sign(c_across(starts_with("logFC__jk_")))
                mean(jk_signs == master_sign)
            }
        ) %>%
        ungroup()
    
    cnames_to_save = colnames(de_df)[grepl("master", colnames(de_df))]
    cnames_to_save = c("gene", "robustness", cnames_to_save)
    de_df = de_df[,cnames_to_save]

    slogan = gsub(".csv", "", basename(master_de_csv_path))
    out_path = file.path(
        dirname(master_de_csv_path),
        paste0(slogan, "_robustness.csv")
    )
    write.csv(de_df, out_path, row.names=F)
    return(de_df)
    
}

get_residual_corrected_scaled_sample_sign_agreement = function(
    res,
    contrast_col,
    case_col,
    contrast_coef_df
){
    # given a df of coefficients of interest and a limma fit object
    # residualize the expression data by removing the effect of all other covariates
    # and calculate what fraction of samples agree with the sign of the coefficients
    #   (e.g. if a coefficient is positive, 
    #   what fraction of case samples are above the median and controls below the median?)
    # scaled by the max and min values achievable given any imbalance in the number of cases vs. controls

    # res: a list, containing 
        # a `normData$EList` matrix of normalized gene expression
        # a `fit` MArrayLM, iself containing
        #   `design` matrix of covariates
        #   `coefficients` matrix of coefficients
        #   `block` vector of sample IDs for each

    # contrast_col: the name of the column in the design matrix that contains the condition of interest
    #   e.g. "condition" or "case_control"
    # case_col: the name of the column in the design matrix that contains the case/control status
    #   will usually contain the contrast_col as a substring 
    #   e.g. "case_controlPD" or "conditionXDP"
    # contrast_coef_df: a dataframe containing the coefficients of interest (logFC)
    #   must have columns logFC and gene (todo: parametrize if need be)

    fit = res[["fit"]]
    design = fit[["design"]]              # nsamp × ncoef
    coef_mat = fit[["coefficients"]]     # ngene × ncoef
    expr = res[["normData"]]$E             # ngene × nsamp
    sample_ids = fit[["block"]]           # nsamp


    # Identify the columns corresponding to contrast_col
    cond_cols = grep(contrast_col, colnames(design))  # or exact match(s)
    non_cond_cols = setdiff(seq_len(ncol(design)), cond_cols)

    # Subset design and coefficients for non-contrast_col terms
    design_nocond = design[, non_cond_cols, drop = FALSE]
    coef_nocond = coef_mat[, non_cond_cols, drop = FALSE]

    # Compute predicted values: (genes × covariates) %*% (covariates × samples)
    fitted_nocond = coef_nocond %*% t(design_nocond)  # ngene × nsamp
    residualized_expr = expr - fitted_nocond

    # We assume the presence of a random effect / sample vector stored in `fit[["block"]]`
    # and we average the residualized expression within every sample
    design_by_sample = as.data.frame(design)
    design_by_sample$block = sample_ids
    design_by_sample = (design_by_sample %>%
        group_by(block) %>%
        summarise(across(everything(), mean, na.rm = TRUE)))

    residualized_expr_by_sample = as.data.frame(t(residualized_expr))
    residualized_expr_by_sample$block = sample_ids

    residualized_expr_by_sample = (residualized_expr_by_sample %>%
        group_by(block) %>%
        summarise(across(everything(), mean, na.rm = TRUE)))

    # we want to figure out how many cases are above / controls below the median residualize expression across samples
    # from every column, subtract the median of that column
    sample_ids = residualized_expr_by_sample$block # remove `$block` for median calculation and add back in later
    residualized_expr_by_sample$block = NULL 
    residualized_expr_by_sample = sweep(
        residualized_expr_by_sample, 2, 
        apply(residualized_expr_by_sample, 2, median, na.rm = TRUE), 
        FUN = "-")
    residualized_expr_by_sample$block = sample_ids

    # Get the sign of the contrast coefficients
    contrast_coef_df$sign = sign(contrast_coef_df$logFC)
    contrast_coef_df = contrast_coef_df[match(colnames(residualized_expr_by_sample), contrast_coef_df$gene), ]
    contrast_coef_df = contrast_coef_df[!is.na(contrast_coef_df$gene), ]
    contrast_coef_df$sign = sign(contrast_coef_df$logFC)

    # Create an outer comparison matrix: samples x genes
    # 1 if (case & gene up) OR (control & gene down)
    # -1 otherwise
    genes = colnames(residualized_expr_by_sample)
    samples = residualized_expr_by_sample$block
    expected_signs = setNames(contrast_coef_df$sign, contrast_coef_df$gene)
    gene_signs = expected_signs[genes]
    sample_is_case = design_by_sample[[case_col]] == 1
    expected_direction_matrix = outer(
        sample_is_case, gene_signs,
        FUN = function(case, direction) {
            as.integer((case & direction == 1) | (!case & direction == -1))
        }
    )
    expected_direction_matrix[expected_direction_matrix == 0] = -1

    # Compute match: 1 if sign matches, 0 otherwise
    common_genes = intersect(colnames(residualized_expr_by_sample), colnames(expected_direction_matrix))
    resid_expr = residualized_expr_by_sample[, common_genes, drop = TRUE]
    expected_dir = expected_direction_matrix[, common_genes, drop = TRUE]
    match_matrix = sign(resid_expr) == sign(expected_dir)
    match_matrix = ifelse(is.na(match_matrix), 0L, as.integer(match_matrix))

    # scale match matrix to [0, 1] using min and max possible values
    case_idx = which(design_by_sample[[case_col]] == 1)
    control_idx = which(design_by_sample[[case_col]] == 0)
    n_samples = nrow(design_by_sample)
    max_possible = (n_samples / 2) + min(length(case_idx), length(control_idx))
    min_possible = (n_samples / 2) - min(length(case_idx), length(control_idx))

    gene_scores = colSums(match_matrix)
    gene_scores_scaled = (gene_scores - min_possible) / (max_possible - min_possible)
    gene_scores_scaled = gene_scores_scaled[contrast_coef_df$gene]
    contrast_coef_df$scaled_sample_sign_agreement = gene_scores_scaled

    return(contrast_coef_df)
}

purity_score_categorical = function(x, group) {
    stopifnot(length(x) == length(group))
    stopifnot(all(group %in% c(0, 1)))

    ranks = rank(x)
    case_ranks = ranks[group == 1]
    control_ranks = ranks[group == 0]

    n_case = length(case_ranks)
    n_control = length(control_ranks)

    # U-statistic: count of (case > control) pairs
    U = sum(outer(case_ranks, control_ranks, FUN = function(c, k) as.numeric(c > k)))

    score = 2 * (U / (n_case * n_control)) - 1
    return(score)
}


get_residual_corrected_purity_score = function(
    res,
    col,
    coef_df,
    categorical=TRUE,
    case_col=NULL
){
    # given a df of coefficients of interest and a limma fit object
    # residualize the expression data by removing the effect of all other covariates
    # and calculate the rank purity score for each gene
    #   (e.g. if a coefficient is positive, 
    #   what fraction of case samples are above all controls (normalized to [-1,1])

    # if a contrast col is given, this is equivalent to the Mann-Whitney U Score
    # if not, the spearman rank correlation is used

    # res: a list, containing 
        # a `normData$EList` matrix of normalized gene expression
        # a `fit` MArrayLM, iself containing
        #   `design` matrix of covariates
        #   `coefficients` matrix of coefficients
        #   `block` vector of sample IDs for each

    # col: the name of the column in the design matrix that contains the condition of interest
    #   e.g. "condition", "case_control", "age"
    # case_col: the name of the column in the design matrix that contains the case/control status
    #   will usually contain the contrast_col as a substring 
    #   e.g. "case_controlPD" or "conditionXDP"
    # coef_df: a dataframe containing the coefficients of interest
    #   must have columns logFC and gene (todo: parametrize if need be)

    fit = res[["fit"]]
    design = fit[["design"]]              # nsamp × ncoef
    coef_mat = fit[["coefficients"]]     # ngene × ncoef
    expr = res[["normData"]]$E             # ngene × nsamp
    block_ids = fit[["block"]]           # nsamp

    # what a crazy bug...if there are colnames in design by sample that are gene names, there will be an error 
    # as they are given .x and .y suffixes later on
    colnames(design)[colnames(design) %in% coef_df$gene] = paste0(
      colnames(design)[colnames(design) %in% coef_df$gene], "_design"
    )

    # Identify the columns corresponding to contrast_col
    cond_cols = grep(col, colnames(design))  # or exact match(s)
    non_cond_cols = setdiff(seq_len(ncol(design)), cond_cols)

    # Subset design and coefficients for non-contrast_col terms
    design_nocond = design[, non_cond_cols, drop = FALSE]
    coef_nocond = coef_mat[, non_cond_cols, drop = FALSE]

    # Compute predicted values: (genes × covariates) %*% (covariates × samples)
    fitted_nocond = coef_nocond %*% t(design_nocond)  # ngene × nsamp
    residualized_expr = expr - fitted_nocond

    # We assume the presence of a random effect / sample vector stored in `fit[["block"]]`
    # and we average the residualized expression within every sample
    design_by_sample = as.data.frame(design)
    design_by_sample$block = block_ids
    design_by_sample = (design_by_sample %>%
        group_by(block) %>%
        summarise(across(everything(), mean, na.rm = TRUE)))

    residualized_expr_by_sample = as.data.frame(t(residualized_expr))
    residualized_expr_by_sample$block = block_ids

    residualized_expr_by_sample = (residualized_expr_by_sample %>%
        group_by(block) %>%
        summarise(across(everything(), mean, na.rm = TRUE)))

    # we want to figure out how many cases are above / controls below the median residualize expression across samples
    # from every column, subtract the median of that column
    sample_ids = residualized_expr_by_sample$block # remove `$block` for median calculation and add back in later
    residualized_expr_by_sample$block = NULL 
    residualized_expr_by_sample = sweep(
        residualized_expr_by_sample, 2, 
        apply(residualized_expr_by_sample, 2, median, na.rm = TRUE), 
        FUN = "-")
    residualized_expr_by_sample$block = sample_ids

    residualized_expr_by_sample = merge(
        residualized_expr_by_sample, 
        design_by_sample, 
        by = "block"
    )

    coef_df = coef_df[match(colnames(residualized_expr_by_sample), coef_df$gene), ]
    coef_df = coef_df[!is.na(coef_df$gene), ]

    if (categorical){
        coef_df$rank_purity = sapply(
            colnames(residualized_expr_by_sample)[!colnames(residualized_expr_by_sample) %in% colnames(design_by_sample)],
            function(gene) {
                purity_score_categorical(residualized_expr_by_sample[[gene]], residualized_expr_by_sample[[case_col]])
            })
    } else {
        coef_df$rank_purity = sapply(
            colnames(residualized_expr_by_sample)[!colnames(residualized_expr_by_sample) %in% colnames(design_by_sample)],
            function(gene) {
                cor(residualized_expr_by_sample[[gene]], residualized_expr_by_sample[[col]], method = "spearman")
            })
    }

    # get absolute value. But first, if the sign of the rank_purity is opposite to the sign of the logFC, set it to 0
    coef_df$rank_purity[(coef_df$rank_purity * coef_df$logFC) < 0] = 0
    coef_df$rank_purity = abs(coef_df$rank_purity)

    return(coef_df)
}

run_deseq = function(counts, md, design_formula, min_frac_gene=0, min_counts_gene=0){
    # counts: a matrix of raw counts (genes x samples)
    # md: a dataframe of metadata (samples x covariates)
    # design_formula: a formula object specifying the design for DESeq2
    # min_frac_gene : minimum fraction of samples that must have non-zero counts for a gene to be kept
    # min_counts_gene : minimum total counts across all samples for a gene to be kept

    # returns a list of dataframes, one per coefficient in the design formula (excluding intercept)

    gene_sums = rowSums(counts)
    gene_nonzero_frac = rowMeans(counts > 0)
    keep_genes = (gene_sums >= min_counts_gene) & (gene_nonzero_frac >= min_frac_gene)
    if (sum(keep_genes) == 0){
        stop("No genes passed the filtering criteria!")
    }
    drop_genes = which(!keep_genes)
    if (length(drop_genes) > 0){
        cat(paste("Dropping", length(drop_genes), "genes that did not pass filtering criteria.\n"))
        counts = counts[keep_genes, ]
    }

    cov_list = all.vars(design_formula)

    # check if all columns in cov_list are in md
    missing_cols = setdiff(cov_list, colnames(md))
    if (length(missing_cols) > 0){
        stop(g("The following columns are missing from the metadata: {paste(missing_cols, collapse=', ')}"))
    }

    # for any categorical variables in cov_list, ensure they are factors
    for (col in cov_list){
        if (is.character(md[[col]]) || is.logical(md[[col]])){
            md[[col]] = as.factor(md[[col]])
        }
    }

    # check that all md columns have more than one unique value
    for (col in cov_list){
        if (length(unique(md[[col]])) < 2){
            stop(g("The column {col} has less than two unique values. DESeq2 requires at least two groups to compare."))
        }
    }
    
    cat(paste("Running DESeq2 with design formula: ", as.character(design_formula), "\n"))
    dds = DESeqDataSetFromMatrix(
        countData = counts,
        colData = md,
        design = design_formula
    )
    dds = DESeq(dds)
    cat("DESeq Complete. Coefficient Names:\n")
    cat(resultsNames(dds))

    result_list = list()
    cat("\nShrinking log fold changes using apeglm...\n")
    for (coef in resultsNames(dds)){
        if (tolower(coef) == "intercept"){
            next
        }
        cat("\n")
        cat(paste("Running lfcShrink for coef: ", coef, "\n"))
        res_shrunk = suppressMessages(suppressWarnings(lfcShrink(dds, coef=coef, type="apeglm")))
        res_shrunk = res_shrunk[!is.na(res_shrunk$padj),]
        res_shrunk = res_shrunk[order(-res_shrunk$log2FoldChange),]

        res_shrunk$gene = rownames(res_shrunk)
        res_shrunk$logFC = res_shrunk$log2FoldChange
        
        res_shrunk$log2FoldChange = NULL
        rownames(res_shrunk) = NULL

        res_shrunk = res_shrunk[,c("gene", "logFC", "padj", "baseMean", "lfcSE", "pvalue")]

        result_list[[coef]] = res_shrunk
    }
    return(result_list)
}


run_mast = function(counts, md, design_formula, min_frac_gene = 0.01, min_counts_gene = 0L, nUMI_col = NULL) {
    # counts: matrix-like (genes x cells), integer-like
    # md: data.frame with rownames matching columns of counts (cells)
    # design_formula: e.g., ~ case_control + sex + age
    # returns: named list of data.frames, one per non-intercept coefficient

    # Filter genes
    gene_sums = Matrix::rowSums(counts)
    gene_nonzero_frac = Matrix::rowMeans(counts > 0)
    keep_genes = (gene_sums >= min_counts_gene) & (gene_nonzero_frac >= min_frac_gene)
    if (sum(keep_genes) == 0) stop("No genes passed the filtering criteria!")
    if (!all(keep_genes)) {
        message(paste("Dropping", sum(!keep_genes), "genes that did not pass filtering criteria."))
        counts = counts[keep_genes, , drop = FALSE]
    }

    # Covariate checks and coercions
    cov_list = all.vars(design_formula)
    missing_cols = setdiff(cov_list, colnames(md))
    if (length(missing_cols) > 0) {
        stop(g("The following columns are missing from the metadata: {paste(missing_cols, collapse=', ')}"))
    }
    for (col in cov_list) {
        if (is.character(md[[col]]) || is.logical(md[[col]])) {
            md[[col]] = as.factor(md[[col]])
        }
        if (length(unique(md[[col]])) < 2) {
            stop(g("The column {col} has less than two unique values; MAST needs variability to estimate effects."))
        }
    }

    # Compute log2(CPM+1) for MAST
    if (!is.null(nUMI_col)) {
        if (!nUMI_col %in% colnames(md)) {
            stop(g("nUMI_col '{nUMI_col}' not found in metadata."))
        }
        libsize = md[[nUMI_col]]
    } else {
        # compute library size from counts
      libsize = Matrix::colSums(counts)
    }
    # guard against zeros
    libsize[libsize == 0] = 1
    cpm = t(t(counts) / libsize) * 1e6
    exprs = log2(cpm + 1)

    # Build SingleCellAssay
    fdata = data.frame(primerid = rownames(exprs), stringsAsFactors = FALSE, row.names = rownames(exprs))
    cdata = md
    if (is.null(rownames(cdata))) {
        stop("Metadata (md) must have rownames matching cell/sample IDs.")
    }
    # Commonly used QC covariate in MAST: detection rate (# genes expressed)
    if (!"cngeneson" %in% colnames(cdata)) {
        cdata$cngeneson = Matrix::colSums(exprs > 0)
    }

    sca = MAST::FromMatrix(
        exprsArray = as.matrix(exprs),
        cData = cdata,
        fData = fdata
    )

    # Fit hurdle model
    message(paste("Fitting MAST zlm with formula:", deparse(design_formula)))
    fit = suppressWarnings(suppressMessages({
      MAST::zlm(formula = design_formula, sca = sca, method = "glm", ebayes = TRUE)
    }))

    # Collect coefficient names (continuous part 'C')
    coefC = coef(fit, "C")  # genes x coefficients
    coef_names = colnames(coefC)
    coef_names = coef_names[coef_names != "(Intercept)"]

    result_list = list()

    # For each coefficient, perform LR test (hurdle) and collect logFC from continuous part
    for (coef_name in coef_names) {
        message(paste("Computing lrTest for:", coef_name))
        lrt = suppressWarnings(suppressMessages({
          MAST::lrTest(fit, coef_name)
        }))
        # lrt is a data.frame; 'hurdle' column has combined p-value; rownames are genes
        # Align to coefC rows
        common_genes = intersect(rownames(coefC), rownames(lrt))
        if (length(common_genes) == 0) next

        pval = lrt[common_genes, "hurdle", "Pr(>Chisq)"]
        padj = p.adjust(pval, method = "BH")
        logFC = coefC[common_genes, coef_name]

        res = data.frame(
            gene = common_genes,
            logFC = as.numeric(logFC),
            pvalue = as.numeric(pval),
            padj = as.numeric(padj),
            stringsAsFactors = FALSE
        )

        # Optional extras: add baseMean-ish proxy = mean expr (log-scale mean is not a mean of counts; include for reference only)
        res$meanExpr_log2CPM = Matrix::rowMeans(exprs[res$gene, , drop = FALSE])
        res$neg_log10_padj = -log10(res$padj + 1e-30)
        res$signed_neg_log10_padj = sign(res$logFC) * res$neg_log10_padj

        res = res[order(res$signed_neg_log10_padj), ]
        rownames(res) = NULL
        res = res[, c("gene", "logFC", "padj", "pvalue", "meanExpr_log2CPM", "neg_log10_padj", "signed_neg_log10_padj")]

        result_list[[coef_name]] = res
    }

    return(result_list)
}
