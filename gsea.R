library(Seurat)
library(dplyr)
library(qs)

library(Matrix)
library(dplyr)

library(fgsea)
library(qs)

source("~/sc-online/plot.R")


.extraHumanGeneAnnoAdderFn=function(inputGeneNames=NULL){
  #require(EnsDb.Hsapiens.v75)
  require(EnsDb.Hsapiens.v86)
  
  if(!dir.exists("~/serverFiles")){
    dir.create("~/serverFiles",recursive = T)
  }
  
  gns <- as.data.frame(genes(EnsDb.Hsapiens.v86))
  gns$gene_short_name=gns$gene_name
  gns$symbol=toupper(gns$symbol)
  gns$ensembl_gene_id=row.names(gns)
  
  if(!is.null(inputGeneNames)){
    rwNames=toupper(inputGeneNames)
    psCols=c("gene_short_name","ensembl_gene_id")
    slCounts=0
    slCol=""
    if(sum(grepl("\\.",rwNames)&grepl("^ENS",rwNames))>0){
      rwNames=strsplit(rwNames,"\\.")
      rwNames=unlist(lapply(rwNames,function(x)x[1]))
    }
    
    # uncomment if these files are not present
    #system(paste0("gsutil -m cp gs://fc-71ac3b81-2441-4171-8038-baf653634620/serverFiles/human_map_to_ensembl.rda ."))
    #system(paste0("gsutil -m cp gs://fc-71ac3b81-2441-4171-8038-baf653634620/serverFiles/human_mapping_hg19.rda ."))
    
    load("~/human_map_to_ensembl.rda")
    load("~/human_mapping_hg19.rda")
    
    map_to_ensmbl$source=toupper(map_to_ensmbl$source)
    
    
    human_hg19$source=toupper(human_hg19$source)
    
    if(sum(toupper(rwNames) %in% human_hg19$source) > sum(toupper(rwNames) %in% map_to_ensmbl$source)){
      map_to_ensmbl=merge(human_hg19,data.frame(source=toupper(rwNames),stringsAsFactors = F),by="source",all.y=T)
    } else {
      map_to_ensmbl=merge(map_to_ensmbl,data.frame(source=toupper(rwNames),stringsAsFactors = F),by="source",all.y=T)
    }
    
    gns=merge(gns,map_to_ensmbl,by.x="ensembl_gene_id",by.y="target",all.y=T)
    gns=gns[match(rwNames,gns$source),]
    row.names(gns)=inputGeneNames
    gns$gene_id=inputGeneNames
    gns=gns[,-which(colnames(gns) %in% c("source","target"))]
  }
  
  return(gns)
}

.sconline.GSEA.readGMT=function (file,bkg_genes=NULL,min.gs.size=NULL,max.gs.size=NULL) {
  if (!grepl("\\.gmt$", file)[1]&F) {
    stop("Pathway information must be in a .gmt file format")
  }
  geneSetDB = readLines(file)
  geneSetDB = strsplit(geneSetDB, "\t")
  names(geneSetDB) = sapply(geneSetDB, "[", 1)
  geneSetDB = lapply(geneSetDB, "[", -1:-2)
  geneSetDB = lapply(geneSetDB, function(x) {
    x[which(x != "")]
  })
  
  if(!is.null(bkg_genes)){
    for(i in 1:length(geneSetDB)){
      tmp=geneSetDB[[i]]
      tmp=bkg_genes[which(toupper(bkg_genes) %in% toupper(tmp))]
      geneSetDB[[i]]=tmp
    }
  }
  
  if(!is.null(min.gs.size)){
    size.dist=unlist(lapply(geneSetDB,length))
    geneSetDB=geneSetDB[size.dist>=min.gs.size]
  }
  
  if(!is.null(max.gs.size)){
    size.dist=unlist(lapply(geneSetDB,length))
    geneSetDB=geneSetDB[size.dist<=max.gs.size]
  }
  
  return(geneSetDB)
}

.myfGSEAfn=function(rankedVec,gs, minSize  = 15,maxSize  = 250, scoreType='std'){
  require(fgsea)
  fgseaRes <- fgsea(pathways = gs, 
                    stats    = rankedVec,
                    minSize  = minSize,
                    maxSize  = maxSize,
                    scoreType= scoreType)
  fgseaRes=fgseaRes[order(fgseaRes$pval,decreasing = F),]
  fgseaRes=as.data.frame(fgseaRes)
  fgseaRes$leadingEdge=unlist(lapply(fgseaRes$leadingEdge,function(x) paste(x,collapse = ",")))
  
  return(fgseaRes)
}

list_to_geneset = function(gene_set_list, path){
    # write a names list of gene sets to a format usable by gsea
    # gene_set_list: list of genesets. Each element is a vector of genes
    # path: path to save the text file
    txt_list = list()
    for(n in names(gene_set_list)){
      gene_list = gene_set_list[[n]]
      s = paste0(
          paste0(n, "\t\t"), paste(gene_list, collapse = "\t"))
      txt_list[[n]] = s
    }
    writeLines(unlist(txt_list), path)
}

runGSEA = function(
    de_df,
    gs_list_of_char,
    rank_col='logFC',
    gene_id_col='gene_short_name',
    desc=FALSE,
    abs=TRUE,
    scoreType='std',
    min_size=15,
    max_size=250){

    de_df=de_df[!duplicated(de_df[[gene_id_col]]),]
    
    order_vector = de_df[[rank_col]]
    if(abs){
        order_vector = abs(order_vector)
    }
    
    ranked_vec=de_df[,rank_col]
    names(ranked_vec)=de_df[[gene_id_col]]
    ranked_vec=ranked_vec[order(order_vector,decreasing = desc)]
    
    print(head(ranked_vec))

    res_fGSEA=.myfGSEAfn(rankedVec=ranked_vec,gs=gs_list_of_char, scoreType=scoreType, minSize=min_size, maxSize=max_size)
    res_fGSEA=res_fGSEA[order(res_fGSEA$padj,decreasing=F),]
    return(res_fGSEA)
}

prep_df_for_gsea = function(dataDE, protein_coding=T){
    dataDE
    anno=.extraHumanGeneAnnoAdderFn(inputGeneNames=dataDE$gene)
    anno=anno[match(dataDE$gene,anno$gene_id),]
    dataDE=cbind(dataDE,anno)
    if(protein_coding){
        dataDE=dataDE[which(dataDE$gene_biotype=="protein_coding"),]
    }
    return(dataDE)
}

prep_gsea_df_for_plotting = function(gsea_df, leading_n=10, display_n=10){

    gsea_df$leading_edge_top_n = sapply(gsea_df$leadingEdge, function(x){
        edge_list = strsplit(x, ",")[[1]]
        len = length(edge_list)
        if (len > leading_n){
            len=leading_n
        }
        return(paste(edge_list[1:len], collapse=", "))
    })

    gsea_df = gsea_df[,c("gene_set", "pathway", "padj", "NES", "size", "leading_edge_top_n")]
    for (gene_set in sort(unique(gsea_df$gene_set))){
        gsea_df_subset = gsea_df[gsea_df$gene_set == gene_set,]
        len = nrow(gsea_df_subset)
        
        gsea_df_subset_pos = gsea_df_subset[gsea_df_subset$NES > 0,]
        gsea_df_subset_neg = gsea_df_subset[gsea_df_subset$NES < 0,]

        len_pos = nrow(gsea_df_subset_pos)
        len_neg = nrow(gsea_df_subset_neg)

        if (len_pos > display_n){
            len_pos=display_n
        }
        if (len_neg > display_n){
            len_neg=display_n
        }
        
        gsea_df_subset_top_n = gsea_df_subset_pos[order(gsea_df_subset_pos$NES, decreasing=TRUE),][1:len_pos,]
        gsea_df_subset_bottom_n = gsea_df_subset_neg[order(gsea_df_subset_neg$NES, decreasing=FALSE),][1:len_neg,]   

        if (all(is.na(gsea_df_subset_top_n$leading_edge_top_n))){
            gsea_df_subset_top_n = NULL
        }
        if (all(is.na(gsea_df_subset_bottom_n$leading_edge_top_n))){
            gsea_df_subset_bottom_n = NULL
        }

        return(list(
            top_n = gsea_df_subset_top_n,
            bottom_n = gsea_df_subset_bottom_n
        ))
    }
}