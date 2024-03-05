source("~/sc-online/utils.R")
source("~/sc-online/extraGenes.R")

.extra_sconline.duplicateCorrelation=function(
    object, 
    design = NULL, 
    ndups = 2, 
    spacing = 1, 
    block = NULL, 
    trim = 0.15,
    weights = NULL){

  require(limma)
  y <- limma:::getEAWP(object)
  print(dim(y$design))
  M <- y$exprs
  ngenes <- nrow(M)
  narrays <- ncol(M)
  if (is.null(design)) 
    design <- y$design
  if (is.null(design)) 
    design <- matrix(1, ncol(y$exprs), 1)
  else {
    design <- as.matrix(design)
    if (mode(design) != "numeric") 
      stop("design must be a numeric matrix")
  }
  
  if (nrow(design) != narrays)
    stop("Number of rows of design matrix does not match number of arrays")
  ne <- limma:::nonEstimable(design)
  if (!is.null(ne)) 
    cat("Coefficients not estimable:", paste(ne, collapse = " "), 
        "\n")
  nbeta <- ncol(design)
  if (missing(ndups) && !is.null(y$printer$ndups)) 
    ndups <- y$printer$ndups
  if (missing(spacing) && !is.null(y$printer$spacing)) 
    spacing <- y$printer$spacing
  if (missing(weights) && !is.null(y$weights)) 
    weights <- y$weights
  if (!is.null(weights)) {
    weights <- asMatrixWeights(weights, dim(M))
    weights[weights <= 0] <- NA
    M[!is.finite(weights)] <- NA
  }
  if (is.null(block)) {
    if (ndups < 2) {
      warning("No duplicates: correlation between duplicates not estimable")
      return(list(cor = NA, cor.genes = rep(NA, nrow(M))))
    }
    if (is.character(spacing)) {
      if (spacing == "columns") 
        spacing <- 1
      if (spacing == "rows") 
        spacing <- object$printer$nspot.c
      if (spacing == "topbottom") 
        spacing <- nrow(M)/2
    }
    Array <- rep(1:narrays, rep(ndups, narrays))
  }
  else {
    ndups <- 1
    nspacing <- 1
    Array <- block
  }
  if (is.null(block)) {
    M <- limma:::unwrapdups(M, ndups = ndups, spacing = spacing)
    ngenes <- nrow(M)
    if (!is.null(weights)) 
      weights <- limma:::unwrapdups(weights, ndups = ndups, spacing = spacing)
    design <- design %x% rep(1, ndups)
  }
  if (!requireNamespace("statmod", quietly = TRUE)) 
    stop("statmod package required but is not installed")
  rho <- rep(NA, ngenes)
  nafun <- function(e) NA
  for (i in 1:ngenes) {
    y <- drop(M[i, ])
    o <- is.finite(y)
    A <- factor(Array[o])
    nobs <- sum(o)
    nblocks <- length(levels(A))
    if (nobs > (nbeta + 2) && nblocks > 1 && nblocks < nobs - 1) {
      y <- y[o]
      X <- design[o, , drop = FALSE]
      Z <- model.matrix(~0 + A)
      if (!is.null(weights)) {
        w <- drop(weights[i, ])[o]
        s <- tryCatch(statmod::mixedModel2Fit(y, X, Z, 
                                              w, only.varcomp = TRUE, maxit = 20)$varcomp, error = nafun)
      } else{
        s <- tryCatch(statmod::mixedModel2Fit(y, X, 
                                              Z, only.varcomp = TRUE, maxit = 20)$varcomp, 
                      error = nafun)
      } 
      if (!is.na(s[1])) 
        rho[i] <- s[2]/sum(s)
    }
  }
  arho <- atanh(pmax(-1, rho))
  mrho <- tanh(mean(arho, trim = trim, na.rm = TRUE))
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
  
  dge <- DGEList(counts=counts(sl_data))
  if(is.null(bkg_genes)){
    tmpCount=apply(counts(sl_data),1,function(x) sum(x>0))
    keep=tmpCount>10
    
  } else {
    keep=row.names(sl_data) %in% bkg_genes
  }
  
  dge <- dge[keep,,keep.lib.sizes=FALSE]
  
  if(TMMnorm){
    dge <- calcNormFactors(dge)
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
      keep <- rowSums(logCPM>3)>(0.05*ncol(logCPM))
    }
    
    logCPM <- logCPM[keep,]
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
    dge <- DGEList(tmpexp)
    
    if(TMMnorm){
      dge <- calcNormFactors(dge)
    }
    
    logCPM <- new("EList")
    logCPM$E <- edgeR::cpm(dge,normalized.lib.sizes = TMMnorm, log = TRUE, prior.count = prior.count)
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
      dc <- .extra_sconline.duplicateCorrelation(logCPM,design=model, block=colData(sl_data)[,random_effect])
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
        fit <- lmFit(logCPM, model,block = as.character(colData(sl_data)[,random_effect]), correlation=dc$consensus.correlation)
        blocked_analysis=T
      } else {
        fit <- lmFit(logCPM, model)
      }
    } else {
      fit <- lmFit(logCPM, model)
    }
  } else {
    fit <- lmFit(logCPM, model)
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
  
  
  dge <- DGEList(counts=counts(sl_data))
  if(is.null(bkg_genes)){
    tmpCount=apply(counts(sl_data),1,function(x) sum(x>0))
    keep=tmpCount>10
    
  } else {
    keep=row.names(sl_data) %in% bkg_genes
  }
  
  dge <- dge[keep,,keep.lib.sizes=FALSE]
  
  if(TMMnorm){
    dge <- calcNormFactors(dge)
  }
  
  if(sample.weights){
    if(quantile.norm){
      logCPM <- voomWithQualityWeights(dge, model, normalize.method="quantile",plot = T,save.plot=T)
    } else {
      logCPM <- voomWithQualityWeights(dge, model,plot = T,save.plot=T)
    }
  } else {
    if(quantile.norm){
      logCPM <- voom(dge, model, normalize.method="quantile",plot = T,save.plot=T)
    } else {
      logCPM <- voom(dge, model,plot = T,save.plot=T)
    }
  }
  
  dc=NULL
  if(!is.null(random_effect)){
    if(is.null(dc.object)){
      dc <- .extra_sconline.duplicateCorrelation(logCPM,design=model, block=colData(sl_data)[,random_effect])
    } else {
      dc=dc.object
    }
  }
  
  
  blocked_analysis=F
  if(!is.null(dc)){
    if(!is.nan(dc$consensus.correlation)){
      if(abs(dc$consensus.correlation)<0.9){
        fit <- lmFit(logCPM, model,block = colData(sl_data)[,random_effect], correlation=dc$consensus.correlation)
        blocked_analysis=T
      } else {
        fit <- lmFit(logCPM, model)
      }
    } else {
      fit <- lmFit(logCPM, model)
    }
  } else {
    fit <- lmFit(logCPM, model)
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

  
  #sl_data=inputExpData;model=model_matrix;random_effect=randomEffect;quantile.norm = quantile.norm;TMMnorm=norm.tmm;VSTnorm=norm.vst;prior.count=prior.count;pd=colData(inputExpData)
  switch(DEmethod,
         Trend=.extra_sconline.Fit_LimmaTrendFn(sl_data=inputExpData,model=model_matrix,random_effect=randomEffect,quantile.norm = quantile.norm,TMMnorm=norm.tmm,VSTnorm=norm.vst,prior.count=prior.count,bkg_genes=bkg_genes,no_normalization=no_normalization,dc.object=dc.object,include.malat1.as.covariate),
         Dream=.extra_sconline.Fit_LimmaDreamFn(sl_data=inputExpData,model=model_matrix,pd=as.data.frame(colData(inputExpData)),random_effect=randomEffect,quantile.norm = quantile.norm,TMMnorm=norm.tmm,VSTnorm=norm.vst,prior.count=prior.count,bkg_genes=bkg_genes,no_normalization=no_normalization,dc.object=dc.object,include.malat1.as.covariate,ncores = dream_ncores),
         Voom=.extra_sconline.Fit_LimmaVoomFn(sl_data=inputExpData,model=model_matrix,random_effect=randomEffect,quantile.norm=quantile.norm,sample.weights=F,TMMnorm=norm.tmm,bkg_genes=bkg_genes,dc.object=dc.object),
         VoomSampleWeights=.extra_sconline.Fit_LimmaVoomFn(sl_data=inputExpData,model=model_matrix,random_effect=randomEffect,quantile.norm=quantile.norm,sample.weights=T,TMMnorm=norm.tmm,bkg_genes=bkg_genes,dc.object=dc.object))
  
}

