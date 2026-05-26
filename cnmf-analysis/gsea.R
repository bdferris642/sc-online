# gsea.R — pipeline-local copy
# Contains only the functions used by the cnmf-analysis pipeline.
# Seurat, plot.R, utils.R and their heavy dependencies (caret, ggforce, rhdf5, etc.)
# are omitted — not needed here and not installed in cnmf-analysis-env.

suppressMessages(suppressWarnings({
  library(dplyr)
  library(fgsea)
  library(qs)
  library(glue)
  g=glue::glue
}))

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

filter_redundant_gene_sets = function(
    gsea_df,
    nes_threshold = 1.0,
    fdr_threshold = 0.05,
    overlap_threshold = 0.21,
    hub_degree_threshold = 12) {

  require(igraph)

  parse_le = function(le_str) {
    if (is.na(le_str) || le_str == "") return(character(0))
    strsplit(le_str, ",")[[1]]
  }

  le_overlap = function(a, b) {
    if (length(a) == 0 || length(b) == 0) return(0)
    length(intersect(a, b)) / min(length(a), length(b))
  }

  select_winners = function(df) {
    if (nrow(df) == 0) return(character(0))
    if (nrow(df) == 1) return(df$pathway)

    le_lists = lapply(df$leadingEdge, parse_le)
    n = nrow(df)

    adj = matrix(FALSE, n, n)
    for (i in seq_len(n - 1)) {
      for (j in seq(i + 1, n)) {
        ov = le_overlap(le_lists[[i]], le_lists[[j]])
        if (ov > overlap_threshold) {
          adj[i, j] = adj[j, i] = TRUE
        }
      }
    }

    active = seq_len(n)
    repeat {
      if (length(active) == 0) break
      sub_adj = adj[active, active, drop = FALSE]
      degs = rowSums(sub_adj)
      hub_local = which(degs > hub_degree_threshold)
      if (length(hub_local) == 0) break
      active = active[-hub_local]
    }

    if (length(active) == 0) return(character(0))

    sub_adj = adj[active, active, drop = FALSE]
    g_sub = graph_from_adjacency_matrix(sub_adj, mode = "undirected")
    comps = components(g_sub)

    winners = character(0)
    for (comp_id in seq_len(comps$no)) {
      comp_local = which(comps$membership == comp_id)
      orig_idx = active[comp_local]
      comp_df = df[orig_idx, , drop = FALSE]

      comp_sub_adj = adj[orig_idx, orig_idx, drop = FALSE]
      comp_df$._comp_degree = rowSums(comp_sub_adj)

      comp_df = comp_df[order(
        -comp_df$._comp_degree,
        comp_df$padj,
        -abs(comp_df$NES)
      ), ]
      winners = c(winners, comp_df$pathway[1])
    }

    return(winners)
  }

  sig_mask = abs(gsea_df$NES) >= nes_threshold & gsea_df$padj < fdr_threshold
  sig_df   = gsea_df[sig_mask,  , drop = FALSE]
  nonsig_df = gsea_df[!sig_mask, , drop = FALSE]

  if (nrow(sig_df) == 0) {
    return(gsea_df)
  }

  pos_df = sig_df[sig_df$NES > 0, , drop = FALSE]
  neg_df = sig_df[sig_df$NES < 0, , drop = FALSE]

  pos_winners = select_winners(pos_df)
  neg_winners = select_winners(neg_df)

  all_winners = c(pos_winners, neg_winners)
  winner_sig_df = sig_df[sig_df$pathway %in% all_winners, , drop = FALSE]

  rbind(winner_sig_df, nonsig_df)
}

plot_gsea_result_hdot = function(
    df,
    title=NULL,
    xlim=c(-4, 4),
    fig_filename = NULL,
    leading_edge_n=10,
    leading_edge_linebreak_n=5,
    pathway_col="pathway",
    top_n=10,
    leading_edge_text_size=4) {

    df_pos = df[df$NES > 0,]
    df_neg = df[df$NES < 0,]

    df_pos = df_pos[order(df_pos$NES, decreasing = TRUE),]
    df_neg = df_neg[order(df_neg$NES, decreasing = FALSE),]

    df = rbind(df_pos[1:top_n,], df_neg[1:top_n,])
    df = df[order(df$NES),]
    df = df[!is.na(df$NES),]

    insert_line_breaks = function(text, n = leading_edge_linebreak_n) {
        parts = strsplit(text, ",\\s*")[[1]]
        if (length(parts) <= n) return(text)
        new_text = ""
        for (i in seq_along(parts)) {
            new_text = paste0(new_text, parts[i])
            if (i < length(parts)) new_text = paste0(new_text, ",")
            if (i %% n == 0 && i != length(parts)) {
                new_text = paste0(new_text, "\n")
            } else if (i != length(parts)) {
                new_text = paste0(new_text, " ")
            }
        }
        return(new_text)
    }

    df$leading_edge_top_n = sapply(df$leadingEdge, function(x){
        edge_list = strsplit(x, ",")[[1]]
        len = min(length(edge_list), leading_edge_n)
        return(insert_line_breaks(paste(edge_list[1:len], collapse=", ")))
    })

    if (is.null(title)) title = "GSEA NES by Pathway"

    df$pathway = gsub("(.{30}\\s)", "\\1\n", df[[pathway_col]], perl = TRUE)
    df$NES_direction = ifelse(df$NES > 0, "NES > 0", "NES < 0")

    p = ggplot(df, aes(x = NES, y = reorder(pathway, NES))) +
        geom_point(aes(size = size, color = NES_direction), shape = 1) +
        scale_color_manual(name = "NES Direction", values = c("NES > 0" = "blue", "NES < 0" = "red")) +
        scale_size_continuous(name = "Size") +
        geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
        labs(title = title, x = "NES", y = "Pathway") +
        xlim(xlim) +
        theme(
            panel.grid.major = element_line(color = "gray", linetype = "dotted"),
            panel.grid.minor = element_blank(),
            legend.position = "left",
            plot.title = element_text(size = 14),
            plot.caption = element_text(hjust = 0.1),
            plot.title.position = "plot",
            axis.line = element_line(color = "black"),
            axis.ticks = element_line(color = "black"),
            axis.text = element_text(size = 13),
            axis.title = element_text(size = 13),
            axis.text.y.right = element_text(hjust = 1)
        )

    if (leading_edge_n > 0) {
        p = p + geom_text(aes(y = reorder(pathway, NES), x = Inf, label = leading_edge_top_n),
                    hjust = -0.1, size = leading_edge_text_size) +
        coord_cartesian(clip = 'off') +
        theme(
            plot.margin = margin(5.5, 340, 5.5, 5.5),
            axis.text.y.right = element_text(hjust = 1)
        )
    }

    if (!is.null(fig_filename)) {
        ggsave(fig_filename, plot = p, width = 14, height = 2*nrow(df)/3 + 1.5, bg="white", dpi=400)
    }

    print(p)
}
