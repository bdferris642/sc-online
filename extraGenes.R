.extraMitoGenes=function(organism,redownload_files=T){
  #organism: Human, Mouse
  
  if(tolower(organism)=="human"){
    gns=.extraHumanGeneAnnoAdderFn()
  } else if (tolower(organism)=="mouse"){
    gns=.extraMouseGeneAnnoAdderFn()
  } else if(tolower(organism)=="macaque"){
    gns=.extraMacaqueGeneAnnoAdderFn()
  } else {
    stop("Wrong organism name!")
  }
  
  gns=gns[gns$seqnames=="MT",]
  gns2=unlist(gns$entrezid)
  gns$entrezgene_id=gns2
  gns$ensembl_gene_id=gns$gene_id
  return(gns)
}

.extrarRNAGenes=function(organism,redownload_files=T){
  #organism: Human, Mouse
  
  if(tolower(organism)=="human"){
    gns=.extraHumanGeneAnnoAdderFn()
  } else if (tolower(organism)=="mouse"){
    gns=.extraMouseGeneAnnoAdderFn()
  } else if(tolower(organism)=="macaque"){
    gns=.extraMacaqueGeneAnnoAdderFn()
  } else {
    stop("Wrong organism name!")
  }
  
  gns=gns[which(gns$gene_biotype=="rRNA"),]
  gns2=unlist(gns$entrezid)
  gns$entrezgene_id=gns2
  gns$ensembl_gene_id=gns$gene_id
  return(gns)
}

.extraRibosomalProteinGenes=function(organism,redownload_files=T){
  #organism: Human, Mouse
  
  if(tolower(organism)=="human"){
    gns=.extraHumanGeneAnnoAdderFn()
  } else if (tolower(organism)=="mouse"){
    gns=.extraMouseGeneAnnoAdderFn()
  } else if(tolower(organism)=="macaque"){
    gns=.extraMacaqueGeneAnnoAdderFn()
  } else {
    stop("Wrong organism name!")
  }
  
  gns=gns[which(grepl("^RP[LS]",gns$gene_short_name,ignore.case = T)),]
  gns2=unlist(lapply(gns$entrezid,function(x) x[1]))
  gns$entrezgene_id=gns2
  gns$ensembl_gene_id=gns$gene_id
  return(gns)
}


.extraOXPHOSGenes=function(organism,inputGeneset=NULL,server=F){
  #organism: Human, Mouse
  
  if(tolower(organism)=="human"){
    gns=.extraHumanGeneAnnoAdderFn()
  } else if (tolower(organism)=="mouse"){
    gns=.extraMouseGeneAnnoAdderFn()
  } else if(tolower(organism)=="macaque"){
    gns=.extraMacaqueGeneAnnoAdderFn()
  } else {
    stop("Wrong organism name!")
  }
  
  if(is.null(inputGeneset)){
      if(!file.exists("~/serverFiles/OXPHOS_gene_symbols.txt")){
        stop("Error! file ~/serverFiles/OXPHOS_gene_symbols.txt is missing")
      }
      
      oxphosGenes=read.table("~/serverFiles/OXPHOS_gene_symbols.txt",sep="\t",stringsAsFactors = F)
    
  }
  
  if(sum(colnames(gns)=="gene_short_name")>0){
    gns=gns[which(toupper(gns$gene_short_name) %in% toupper(oxphosGenes[,1])),]
  } else {
    gns=gns[which(toupper(gns$symbol) %in% toupper(oxphosGenes[,1])),]
  }
  
  
  return(gns)
}

.extraIEGGenes=function(organism,server=F){
  #organism: Human, Mouse
  
  gns=.extraMouseGeneAnnoAdderFn()
  
  if(!file.exists("~/serverFiles/IEG_gene_symbols")){
    #library(googleCloudStorageR)
    stop("Error! file serverFiles/IEG_gene_symbols is missing")
  }
  IEGenes=read.table("~/serverFiles/IEG_gene_symbols",sep="\t",stringsAsFactors = F)
  
  gns=gns[tolower(gns$gene_name) %in% tolower(IEGenes$V1),]
  
  if(tolower(organism)=="human"){
    
    {
      #library(googleCloudStorageR)
      if(!file.exists("~/serverFiles/ortholog_mapping.rda")){
        stop("Error! file ~/serverFiles/ortholog_mapping.rda is missing")
        #gcs_get_object("vgazesta/serverFiles/orthologsFeb3/ortholog_mapping.rda", saveToDisk = "~/serverFiles/ortholog_mapping.rda",overwrite=T)
      }
      
      load("~/serverFiles/ortholog_mapping.rda")
    }
    
    
    
    mapping=mapping[mapping$mouse %in% gns$ensembl_gene_id,]
    gnsHuman=.extraHumanGeneAnnoAdderFn()
    gns=gnsHuman[gnsHuman$ensembl_gene_id %in% mapping$human,]
  } else if (tolower(organism)=="mouse"){
    
  } else if(tolower(organism)=="macaque"){
    gns=.extraMacaqueGeneAnnoAdderFn()
    gns=gns[which(tolower(gns$symbol) %in% tolower(IEGenes$V1)|tolower(gns$gene_name) %in% tolower(IEGenes$V1)),]
  } else {
    stop("Wrong organism name!")
  }
  
  gns$ensembl_gene_id=gns$gene_id
  return(gns)
}

.extraCellCycleGenes=function(organism){
  if(tolower(organism)=="human"){
    gns=.extraHumanGeneAnnoAdderFn()
  } else if (tolower(organism)=="mouse"){
    gns=.extraMouseGeneAnnoAdderFn()
  } else if(tolower(organism)=="macaque"){
    gns=.extraMacaqueGeneAnnoAdderFn()
  } else {
    stop("Wrong organism name!")
  }
  
  cc.s=Seurat::cc.genes.updated.2019$s.genes
  cc.g2m=Seurat::cc.genes.updated.2019$g2m.genes
  
  cc.g2m[cc.g2m=="PIMREG"]="FAM64A"
  cc.g2m[cc.g2m=="JPT1"]="HN1"
  
  cc.s=gns[toupper(gns$symbol) %in% toupper(cc.s),]
  cc.g2m=gns[toupper(gns$symbol) %in% toupper(cc.g2m),]
  
  cc.s$ensemble_gene_id=row.names(cc.s)
  cc.g2m$ensemble_gene_id=row.names(cc.g2m)
  
  return(list(cc.s=cc.s,cc.g2m=cc.g2m))
}