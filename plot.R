source("~/sc-online/utils.R")

library(caret)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(grid)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(rhdf5)
library(rlang)
library(tidyr)
library(viridis)
library(viridisLite)
library(textplot)

DimPlotIntegration = function(
    seurat_obj, 
    reduction="umap",
    dataset_col="dataset",
    source_label='tk',
    target_label='me',
    source_cluster_col="madeCluster",
    target_cluster_col="seurat_clusters",
    label=FALSE
){
    seurat_obj_source = seurat_obj[, seurat_obj[[dataset_col]] == source_label]
    seurat_obj_target = seurat_obj[, seurat_obj[[dataset_col]] == target_label]
    print(
        DimPlot(seurat_obj, reduction=reduction, group.by=dataset_col, label=label) 
        + ggtitle(paste('Joint', reduction, 'plot')))
    print(
        DimPlot(seurat_obj_source, reduction=reduction, group.by=source_cluster_col, label=label)
        + ggtitle(paste("Source", reduction, 'plot')))
    print(DimPlot(seurat_obj_target, reduction=reduction, group.by=target_cluster_col, label=label)
        + ggtitle(paste("Target", reduction, 'plot')))
}

getKneePlotData=function(
    df_list,
    plot_col,
    asc=FALSE,
    log_y_axis=FALSE,
    downsample_step=1000,
    downsample_start_ind=10000){
    max_len_ind = which.max(sapply(df_list, nrow))
    max_len = nrow(df_list[[max_len_ind]])
    x = orderDFByColRank(df_list[[max_len_ind]], plot_col, asc, log_y_axis)$rank
    x = downsampleListAfterStart(x, downsample_start_ind, downsample_step)
    
    y_list = lapply(df_list, function(df){
        padded = rep(0, max_len)
        ranked_col = orderDFByColRank(df, plot_col, asc, log_y_axis)[[plot_col]]
        padded[1:length(ranked_col)] = ranked_col
        padded = downsampleListAfterStart(padded, downsample_start_ind, downsample_step)
        return(padded)
    })
    return( list(x, y_list) )
}


plotColByRank = function(
    df,
    col,
    title=NULL,
    color_col=NULL,
    log_x_axis=TRUE,
    log_y_axis=TRUE,
    asc=FALSE,
    ylim=c(0, 1),
    clim=c(0, 1),
    plot_width=8,
    plot_height=8,
    title_font_size=16,
    axis_label_font_size=14,
    tick_label_font_size=11,
    fig_filename=NULL){
    
    # Useful for single-library knee plots

    options(repr.plot.width=plot_width, repr.plot.height=plot_height)

    if (class(df) != 'data.frame'){
        df = as.data.frame(df)
    }
    
    df = orderDFByColRank(df, col, asc=asc, log_y_axis=log_y_axis)

    if (is.null(title)){
        title = paste(col, 'by rank', sep=' ')
    }

    p = ggplot(df, aes_string(x = 'rank', y = col))

    # If color_col is provided, map it to the color aesthetic
    if (!is.null(color_col)) {
        p = p + geom_line(aes_string(color = color_col)) + scale_fill_viridis_c(limits = clim)
    } else {
        p = p + geom_line()
    }

    p = p + ggtitle(title) + 
        xlab('rank') +
        coord_cartesian(ylim=ylim)

    if(log_x_axis){
        p = p + scale_x_log10()
    }

    if(log_y_axis){
        p = p + scale_y_log10() + ylab(paste0('log10(', col, ')'))
    } else {
        p = p + ylab(col)
    }
        
    p = p + theme(
            plot.title = element_text(
                hjust = 0.5, # This centers the title
                size = title_font_size),          # Change plot title font size
            axis.title.x = element_text(size = axis_label_font_size),        # Change x axis label font size
            axis.title.y = element_text(size = axis_label_font_size),        # Change y axis label font size
            axis.text.x = element_text(size = tick_label_font_size),         # Change x axis tick labels font size (optional)
            axis.text.y = element_text(size = tick_label_font_size)          # Change y axis tick labels font size (optional)
        )
        
    print(p)

    if (!is.null(fig_filename)){
        ggsave(fig_filename, width=plot_width, height=plot_height, dpi=600)
    }
}

plotColorLines = function(x, y_list, line_color_by_list=NULL, clab=NULL, cmap=viridis, 
                           ylim=NULL, clim=c(0, 1), log_x_axis=FALSE, log_y_axis=FALSE, 
                           plot_width=8, plot_height=8, xlim=c(1, 3e5),
                           xlab="", ylab="", main="", show_legend=TRUE,
                           title_font_size=18, axis_label_font_size=16, tick_label_font_size=12,
                           display_plot_repr=TRUE, fig_filename=NULL, reverse=FALSE) {

    line_color_by_list = unlist(line_color_by_list)
  
    # Combine x and y_list into a data frame
    plot_data <- do.call(rbind, lapply(seq_along(y_list), function(i) {
        data.frame(x = x, y = y_list[[i]], line_id = as.factor(i))
    }))

    # Determine and set y-axis limits if ylim is NULL
    if (is.null(ylim)) {
        ymin <- min(sapply(y_list, min, na.rm = TRUE))
        ymax <- max(sapply(y_list, max, na.rm = TRUE))
        ylim <- c(ymin, ymax)
    }

    # Determine color scale type (discrete or continuous)
    if (!is.null(line_color_by_list)) {
        if (is.numeric(line_color_by_list)) {
            plot_data$color = line_color_by_list[plot_data$line_id]
            color_scale = selectColorScale(cmap, TRUE, clim, clab, num_factors=NULL, reverse=reverse)
        } else {
            plot_data$color = factor(line_color_by_list[plot_data$line_id])
            color_scale = selectColorScale(cmap, FALSE, clim, clab, num_factors=length(levels(factor(plot_data$color))), reverse=reverse)
        }
        } else {
        plot_data$color = names(y_list)[plot_data$line_id]
        color_scale = selectColorScale(cmap, FALSE, clim, clab, num_factors=length(levels(factor(plot_data$color))), reverse=reverse) 
    }
    

    # Create the ggplot
    p <- ggplot(plot_data, aes(x = x, y = y, group=line_id, color = color)) +
        geom_line() +
        color_scale +
        xlab(xlab) +
        ylab(ylab) +
        ggtitle(main) +
        theme_minimal() +
        coord_cartesian(ylim = ylim, xlim=xlim)

    if (show_legend) {
        leg_pos = "right"
    } else {
        leg_pos = "none"
    }

    p = p+theme(
        plot.title = element_text(
                hjust = 0.5,                                        # This centers the title
                size = title_font_size),                            # Change plot title font size
        axis.title.x = element_text(size = axis_label_font_size),   # Change x axis label font size
        axis.title.y = element_text(size = axis_label_font_size),   # Change y axis label font size
        axis.text.x = element_text(size = tick_label_font_size),    # Change x axis tick labels font size (optional)
        axis.text.y = element_text(size = tick_label_font_size),    # Change y axis tick labels font size (optional)
        legend.position = leg_pos                               # Change legend position
    )

    # Apply log transformations if requested
    if (log_x_axis) { p <- p + scale_x_log10() }
    if (log_y_axis) { p <- p + scale_y_log10() }

    # Save or display the plot
    options(repr.plot.height=plot_height, repr.plot.width=plot_width)
    if (display_plot_repr) { print(p) }

    if (!is.null(fig.filename)) { ggsave(fig_filename, width=plot_width, height=plot_height, dpi=600) }
}


# Function to plot the heatmap with numbers inside
plotConfusionHeatmap = function(prop_table, xlab='', ylab='', title='', 
    plot_width=10, plot_height=12, fig_filename=NULL) {
  options(repr.plot.width=plot_width, repr.plot.height=plot_height)
  
  # Prepare matrix with rounded proportions for display
  display_matrix <- matrix(sprintf("%.2f", prop_table), nrow = nrow(prop_table))

  # Create a color palette from viridis with 100 levels mapped between 0 and 1
  color_palette <- viridis(101, alpha = 1, begin = 0, end = 1)

  labs = lapply(strsplit(rownames(prop_table), '_'), function(x) x[[length(x)]])
  # Plot the heatmap
  p <- pheatmap(
      prop_table,
      color = color_palette,
      #breaks = seq(0, 1, length.out = 101), # set breaks between a specific range
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      display_numbers = display_matrix,
      number_color='White',
      fontsize = 12,
      fontsize_number = 10,
      labels_row=labs,
      labels_col=labs,
      main=title
  )
  grid.text(xlab, x=0.43, y=0.05, rot=0)
  grid.text(ylab, x=0.95, y=0.53, rot=0)


  print(p)

  if (!is.null(fig_filename)) {
      ggsave(fig_filename, plot=p, width=plot_width, height=plot_height, dpi=600)
  }
}

plotDfSummary = function(df, fig_filename){
    width=800
    height=30*nrow(df)

    png(fig_filename, width = width, height = height)
    grid_table = tableGrob(df)
    grid.draw(grid_table)
    dev.off()
}


plot_gsea_result_hbar = function(
    gsea_df_subset, title, xlim = c(-4, 4), clim=c(-4, 4), 
    fig_filename=NULL, plot_width=10, plot_height=6){

    gsea_df_subset$rank = rank(gsea_df_subset$NES)
    len = nrow(gsea_df_subset)
    if (len == 0){
        return()
    }

    # every 30 characters, replace the next space with a \n and then start over
    gsea_df_subset$pathway = gsub("(.{30}\\s)", "\\1\n", gsea_df_subset$pathway, perl=TRUE)
    gsea_df_subset$pathway = factor(gsea_df_subset$pathway, levels = rev(gsea_df_subset$pathway))
    
    p=ggplot(gsea_df_subset, aes(x=NES, y=pathway, color=NES, fill=NES)) +
        geom_bar(stat="identity", position="identity") +
        scale_fill_gradientn(colors = RColorBrewer::brewer.pal(11, "RdBu"),limits = clim) +
        scale_y_discrete(limits = levels(gsea_df_subset$pathway)) +
        scale_x_continuous(limits = xlim) +
        ggtitle(title) +
        labs(x="NES", y="Gene Set") +
        theme(
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            legend.position = "none",
            plot.title = element_text(size=14),
            plot.caption = element_text(hjust = 0.1),
            plot.title.position = "plot",
            axis.line = element_line(color = "black"),  # Add axis lines
            axis.ticks = element_line(color = "black"),  # Add axis ticks
            axis.text = element_text(size = 13),  # Increase tick label font size
            axis.title = element_text(size = 13)  # Increase axis label font size
        )
    print(p)
    if (!is.null(fig_filename)){
        ggsave(fig_filename, plot=p, width=plot_width, height=plot_height, dpi=600, bg="white")
    }
}

plot_gsea_result_hdot = function(
    df, 
    title=NULL, 
    xlim=c(-4, 4), 
    fig_filename = NULL, 
    leading_edge_n=10,
    leading_edge_linebreak_n=5,
    top_n=10) {
    #takes in a df with the following columns: NES, pathway, size, leading_edge
    # along with a title, xlim, fig_filename (nullable)
    # leading_edge_n is the max number of leading_edge genes to be annotated
    # leading_edge_linebreak_n is the number of leading_edge genes to be displayed before a line break
    # top_n is the number of pathways to be displayed on the plot for each direction 
    # (i.e. max rows is 2*top_n for up- and down-regulated genesets)

    # the df is filtered to the top_n pathways by NES above 0 and below 0

    # Y axis: pathway. every 30 characters, replace the next space with a \n and then start over
    # X axis: NES, appearing as an open circle. If NES > 0, this circle is blue, < 0 red. The size of the circle relates to the Gene Set Size. 
    # The plot is annotated on the right with the leading_edge_top_n.
    
    # first, turn the leading edge into something we can read by taking the first `leading_edge_n` elemnets
    # separating them with `, ` and finally adding line breaks
    df_pos = df[df$NES > 0,]
    df_neg = df[df$NES < 0,]

    df_pos = df_pos[order(df_pos$NES, decreasing = TRUE),]
    df_neg = df_neg[order(df_neg$NES, decreasing = FALSE),]

    df = rbind(df_pos[1:top_n,], df_neg[1:top_n,])
    df = df[order(df$NES),]
    df = df[!is.na(df$NES),]

    # Function to insert line breaks every `n` entries
    insert_line_breaks = function(text, n = leading_edge_linebreak_n) {
        parts = strsplit(text, ",\\s*")[[1]]
        if (length(parts) <= n) {
            return(text)
        }
        new_text = ""
        for (i in seq_along(parts)) {
            new_text = paste0(new_text, parts[i])
            if (i < length(parts)) {
                new_text = paste0(new_text, ",")
            }
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
        len = length(edge_list)
        if (len > leading_edge_n){
            len=leading_edge_n
        }
        return(insert_line_breaks(paste(edge_list[1:len], collapse=", ")))
    })

    if (is.null(title)){
        title = "GSEA NES by Pathway"
    }

    # Format the pathway names to include newline every 30 characters
    df$pathway = gsub("(.{30}\\s)", "\\1\n", df$pathway, perl = TRUE)
    df$NES_direction = ifelse(df$NES > 0, "NES > 0", "NES < 0")
    
    # Create the plot
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
    
        # Create a secondary y-axis with leading_edge_top_n labels
    p = p + geom_text(aes(y = reorder(pathway, NES), x = Inf, label = leading_edge_top_n), 
                      hjust = -0.1, size = 4) +
        coord_cartesian(clip = 'off') +
        theme(
            plot.margin = margin(5.5, 340, 5.5, 5.5),  # Increase right margin to make space for secondary labels
            axis.text.y.right = element_text(hjust = 1)
        )

    
    # Save the plot to file if a filename is provided
    if (!is.null(fig_filename)) {
        ggsave(fig_filename, plot = p, width = 14, height = 2*nrow(df)/3 + 1.5, bg="white", dpi=400)
    }
    
    # Print the plot to the R console
    print(p)
}


plotHistGrid=function(
    df_list,
    col,
    breaks,
    title_list,
    normalize=TRUE,
    xlab=NULL,
    ylim=c(0, 1),
    nrows=3, 
    ncols=4, 
    plot_width=8,
    plot_height=8){
    
    ngrid=nrows*ncols
    n=length(df_list)
    
    options(repr.plot.width=plot_width, repr.plot.height=plot_height)
    par(mfrow = c(nrows, ncols))
    for (i in seq(1, n, ngrid)){  
        for (j in 0:(ngrid-1)){
            name = names(df_list)[[i+j]]
            df = df_list[[name]]
            h = hist(df[[col]], plot = FALSE, breaks=breaks)
            if (normalize){
                h$counts = h$counts / sum(h$counts)
                ylab = "Probability"
            } else {
                ylab = "Count"
            }

            if (is.null(xlab)){
                xlab = col
            }

            barplot(
                probs, 
                names.arg = h$mids, 
                col="gray", 
                main=title_list[[name]], 
                xlab=xlab,
                ylab=ylab)
        }
    }
}



plotHistColorGrid=function(
    df_list,
    plot_col,
    color_col,
    title_list=NULL,
    bin_width=0.25,
    xlim=c(0, 1),
    ylim=c(0, 1),
    clim=c(0, 1),
    nrows=3,
    ncols=3,
    plot_width=16,
    plot_height=12,
    title=NULL,
    xlab=NULL,
    ylab=NULL,
    fig_filename=NULL
){

    options(repr.plot.width=plot_width, repr.plot.height=plot_height)
    breaks_vector=seq(xlim[1], xlim[2], bin_width)
    n=length(df_list)
    ngrid=nrows*ncols

    if (is.null(xlab)){
        xlab = plot_col
    }

    for (i in seq(1, n, ngrid)) {  
        plot_list = list()
        for (j in 0:(ngrid-1)) {
            if (i+j > n) break

            row = j%/%nrows+1
            col = j%%ncols+1

            name = names(df_list)[[i+j]]
            df = df_list[[name]]

            if (is.null(title_list)){
                title = paste0(name, "\nHistogram of ", plot_col, "\ncolored by ", color_col)
            } else {
                title = title_list[[name]]
            }

            # Create bins and calculate average color_col for each bin
            df_bins = df %>%
            mutate(
                bin_ind=cut(df[[plot_col]], breaks = breaks_vector, include.lowest = TRUE, right = FALSE, labels=FALSE),
                bin = as.numeric(lapply(bin_ind, function(x) {
                    if (is.na(x)) {
                        return(0)
                    } else {
                        return(breaks_vector[x])
                    }
                }))) %>%
            group_by(bin) %>%
            summarize(
                avg_color_col=mean(!!rlang::sym(color_col), na.rm = TRUE), #!!rlang::sym converts string to variable name
                count=n())
            df_bins = as.data.frame(df_bins)
            df_bins['plot_col_prob'] = df_bins$count / sum(df_bins$count)
            
            all_bins = data.frame(
                bin = breaks_vector[1:length(breaks_vector)-1]
            )
            missing_bins = anti_join(all_bins, df_bins, by = "bin")

            df_bins = arrange(bind_rows(missing_bins, df_bins), bin)
            df_bins[is.na(df_bins)] = 0

            # shift so that left edge is on left edge of bin
            df_bins$bin = df_bins$bin + bin_width/2

            p = ggplot(df_bins, aes(x = bin, y = plot_col_prob, fill = avg_color_col)) +
                geom_bar(stat = "identity") +
                scale_fill_viridis_c(limits = clim) +
                scale_x_continuous(breaks = seq(xlim[1], xlim[2], 1)) +
                ggtitle(title) +
                xlab(xlab) +
                ylab("Probability") + 
                coord_cartesian(xlim=xlim, ylim=ylim)

            # color legend on top right
            if (row==1 & col==ncols){
                p = p + theme(legend.position = "right") + labs(fill = color_col)
            } else {
                p = p + theme(legend.position = "none")
            }

            plot_list[[j+1]] = p

        }
        # Combine all plots
        combined_plot = wrap_plots(plot_list, ncol = ncols)

        # Final plot with the legend
        final_plot = combined_plot + plot_layout(guides = "collect")

        # Draw the final plot
        print(final_plot)

        if (! is.null(fig_filename)){
            ggsave(
                fig_filename,
                plot=final_plot, 
                width=plot_width, 
                height=plot_height,
                bg = 'white',
                dpi=600)
        }
    }
}

plotHistColorSingle = function(
    df,
    name,
    plot_col,
    color_col,
    title=NULL,
    bin_width=0.25,
    xlim=c(0, 1),
    ylim=c(0, 1),
    plot_width=8,
    plot_height=8,
    xlab=NULL,
    ylab=NULL,
    fig_filename=fig_filename
){
    df_list = list()
    title_list = list()
    df_list[[name]]=df
    title_list[[name]]=title

    plotHistColorGrid(
        df_list=df_list,
        plot_col=plot_col,
        color_col=color_col,
        bin_width=bin_width,
        title_list=title_list,
        xlim=xlim,
        ylim=ylim,
        nrows=1, 
        ncols=1, 
        plot_width=plot_width,
        plot_height=plot_height,
        xlab=xlab,
        ylab=ylab,
        fig_filename=fig_filename
    )
}

plotKneeColorLines=function(
    df_list,
    plot_col, 
    line_color_by_list=NULL,
    clab=NULL, cmap="viridis", 
    ylim=NULL, xlim=NULL, clim=c(0, 1), 
    log_x_axis=FALSE, log_y_axis=FALSE, asc=FALSE,
    plot_width=12, plot_height=8, xlab="", ylab="", main="", 
    display_plot_repr=TRUE, fig_filename=NULL,
    downsample_step=1000, downsample_start_ind=10000, show_legend=TRUE,
    return_plot_data=FALSE, reverse=FALSE){

    knee_plot_data = getKneePlotData(
        df_list=df_list, 
        plot_col=plot_col, 
        asc=asc, 
        log_y_axis=log_y_axis, 
        downsample_step=downsample_step, 
        downsample_start_ind=downsample_start_ind)
    x=knee_plot_data[[1]]
    y_list=knee_plot_data[[2]]
    
    plotColorLines(
        x=x, y_list=y_list, line_color_by_list=line_color_by_list, 
        xlim=xlim, ylim=ylim, clim=clim, cmap=cmap, 
        xlab=xlab, ylab=ylab, main=main, clab=clab, 
        log_x_axis=log_x_axis, log_y_axis=log_y_axis, 
        plot_width=plot_width, plot_height=plot_height, show_legend=show_legend,
        display_plot_repr=display_plot_repr, reverse=reverse, fig_filename=fig_filename)
    
    if (return_plot_data) {
        return(knee_plot_data)
    }
}

plotKneeSingle=function(
    df,
    name,
    plot_col='nUMI',
    color_col='is_in_filtered',
    clab="Called Cell", cmap="RdBu", 
    xlim=c(1, 1e6), # oh it will absolutely throw an error if xlim[[1]]==0
    ylim=c(0.5, 5.5), 
    clim=c(.15, 0.8), # I just like these shades of Blue and Red
    log_x_axis=TRUE, 
    log_y_axis=TRUE,
    asc=FALSE,
    plot_width=12,
    plot_height=8, 
    xlab=NULL, ylab=NULL, title=NULL, 
    downsample_step=1000, downsample_start_ind=10000, show_legend=TRUE,
    return_plot_data=FALSE, reverse=FALSE,
    display_plot_repr=TRUE, fig_filename=NULL,
    title_font_size=16,
    axis_label_font_size=14,
    tick_label_font_size=11,
    discretize_color_col=TRUE
){

        if (is.null(xlab)){
            xlab = 'Rank'
        }

        if (is.null(ylab)){
            if (log_y_axis) {
                ylab_prefix = 'log10 '
            } else{
                ylab_prefix = ''
            }
            ylab = paste0(ylab_prefix, plot_col)
        }

        if (is.null(title)){
            title = paste0('Knee Plot for library ', name)
        }
        
        # TODO: update getKneePlotData to take in a list of plot_cols, not just one

        df_list = list()
        df_list[[name]]=df
        xy = getKneePlotData(
            df_list=df_list, 
            plot_col=plot_col, 
            asc=asc, 
            log_y_axis=log_y_axis, 
            downsample_step=downsample_step, 
            downsample_start_ind=downsample_start_ind
        )
        xc = getKneePlotData(
            df_list=df_list, 
            plot_col=color_col, 
            asc=asc, 
            log_y_axis=FALSE, 
            downsample_step=downsample_step, 
            downsample_start_ind=downsample_start_ind
        )
        x=xy[[1]]
        y=xy[[2]][[1]] # will be a list of a single element, whose length will be length(x)
        c=xc[[2]][[1]] # will be a list of a single element, whose length will be length(x)

        if (discretize_color_col) {c = as.factor(c)}

        plot_data = data.frame(
            x=x,
            y=y,
            color_col=c
        )

        if (is.numeric(plot_data$color_col)) {
            color_scale = selectColorScale(cmap, TRUE, clim, clab, num_factors=NULL, reverse=reverse)
        } else {
            color_scale = selectColorScale(cmap, FALSE, clim, clab, num_factors=length(levels(factor(plot_data$color_col))), reverse=reverse)
        }

        # Create the ggplot
        p <- ggplot(plot_data, aes(x = x, y = y,  color = color_col)) +
            geom_line() +
            color_scale

        if(log_x_axis){
            p = p + scale_x_log10()
        }
        p = p +
            xlab(xlab) +
            ylab(ylab) +
            ggtitle(title) +
            theme_minimal() +
            coord_cartesian(ylim = ylim, xlim=xlim)

        if (show_legend) {
            leg_pos = "right"
        } else {
            leg_pos = "none"
        }

        p = p+theme(
            plot.title = element_text(
                    hjust = 0.5,                                        # This centers the title
                    size = title_font_size),                            # Change plot title font size
            axis.title.x = element_text(size = axis_label_font_size),   # Change x axis label font size
            axis.title.y = element_text(size = axis_label_font_size),   # Change y axis label font size
            axis.text.x = element_text(size = tick_label_font_size),    # Change x axis tick labels font size (optional)
            axis.text.y = element_text(size = tick_label_font_size),    # Change y axis tick labels font size (optional)
            legend.position = leg_pos                               # Change legend position
        )

        print(p)

        if (!is.null(fig_filename)){
            ggsave(fig_filename, plot=p, width=plot_width, height=plot_height, dpi=600, bg='white')
        }
    }

plotOverlappingProbabilityHistograms = function(
    long_df,
    plot_col,
    group_col,
    breaks=seq(0, 1, by=0.05),
    title=NULL,
    xlim=c(0,1)){

  if (is.null(title)) {
      title = paste0("Probability Histogram of ", plot_col, " by ", group_col)
  }
  
  long_df['plot_col'] = long_df[[plot_col]]
  long_df['group_col'] = long_df[[group_col]]

  # Binning the data and calculating normalized counts
  bin_width = breaks[2] - breaks[1]
  plot_df = long_df %>%
    mutate(bin_str = cut(plot_col, breaks = breaks, include.lowest = TRUE), # the bins will be strings like '[0.05,0.1)'
      bin = sapply(as.character(bin_str), function(x) { # extract the number of the lower bound
          lbound_str = strsplit(x, ",")[[1]][[1]]
          return(as.numeric(substr(lbound_str, 2, nchar(lbound_str)))) # remove the '[' and convert to numeric
    })) %>%
    group_by(group_col, bin) %>%
    dplyr::summarize(count = n()) %>%
    ungroup() %>%
    group_by(group_col) %>%
    dplyr::mutate(probability = count / sum(count)) %>%
    ungroup()

  plot_df$bin = plot_df$ bin + bin_width/2
  # Plotting the normalized histograms
  p = ggplot(plot_df, aes(x = bin, y = probability, fill = group_col)) +
    geom_bar(stat = "identity", position = "identity", alpha = 0.4, width = bin_width) +
    xlim(xlim) +
    theme_minimal() +
    labs(title = title, x = "Value", y = "Probability")
  
  print(p)
}

plotScatterColorGrid=function(
    df_list,
    x_col,
    y_col,
    color_col,
    title_list,
    xlab=NULL,
    ylab=NULL,
    clab=NULL,
    xlim=c(0, 1),
    ylim=c(0, 1),
    clim=c(0, 1),
    nrows=3,
    ncols=3,
    plot_height=16,
    plot_width=12,
    fig_filename=NULL
){

    n=length(df_list)
    ngrid=nrows*ncols

    if (is.null(xlab)){
        xlab = x_col
    }

    if (is.null(ylab)){
        ylab = y_col
    }

    if (is.null(clab)){
        clab = color_col
    }

    for (i in seq(1, n, ngrid)) {  
        plot_list = list()
        for (j in 0:(ngrid-1)) {
            if (i+j > n) break

            name = names(df_list)[[i+j]]
            df = df_list[[name]]

            row = j%/%nrows+1
            col = j%%ncols+1

            if (is.null(title_list) || is.null(title_list[[name]])){
                title = paste0(name, "\nScatter of ", y_col, " vs: ", x_col, " colored by ", color_col)
            } else {
                title = title_list[[name]]
            }

            p = ggplot(df, aes_string(x=x_col, y=y_col, color=color_col)) +
                geom_point(alpha=0.5) +
                scale_color_viridis_c(limits = clim) +
                scale_x_continuous(breaks = seq(xlim[1], xlim[2], 1)) +
                ggtitle(title) +
                xlab(xlab) +
                ylab(ylab) + 
                coord_cartesian(xlim=xlim, ylim=ylim)

            # color legend on top right
            if (row==1 & col==ncols){
                p = p + theme(legend.position = "right") + labs(color = color_col)
            } else {
                p = p + theme(legend.position = "none")
            }

            plot_list[[j+1]] = p

        }
        # Combine all plots
        combined_plot = wrap_plots(plot_list, ncol = ncols)

        # Final plot with the legend
        final_plot = combined_plot + plot_layout(guides = "collect")

        # Draw the final plot
        print(final_plot)

        if (! is.null(fig_filename)){
            ggsave(
                fig_filename,
                plot=final_plot, 
                width=plot_width, 
                height=plot_height,
                dpi=600)
        }
    }
}

plotScatterColorSingle=function(
    df,
    name,
    x_col,
    y_col,
    color_col,
    title=NULL,
    xlab=NULL,
    ylab=NULL,
    clab=NULL,
    xlim=c(0, 1),
    ylim=c(0, 1),
    clim=c(0, 1),
    plot_height=8,
    plot_width=8,
    fig_filename=NULL
){
    df_list = list()
    title_list = list()
    df_list[[name]]=df
    title_list[[name]]=title

    plotScatterColorGrid(
        df_list=df_list,
        x_col=x_col,
        y_col=y_col,
        color_col=color_col,
        title_list=title_list,
        xlim=xlim,
        ylim=ylim,
        clim=clim,
        nrows=1, 
        ncols=1, 
        plot_width=plot_width,
        plot_height=plot_height,
        fig_filename=fig_filename
    )
}


plotStackedBarByGroup = function(
    df,
    x_col, 
    y_col,
    fill_col,
    plot_title="",
    x_axis_label="",
    y_axis_label="",
    x_order=NULL,
    types=NULL){

    # this function takes a long dataframe (for instance, the output of utils::getCountProportionDF), 
    # and plots a stacked bar plot of y_col, with each bar corresponding to a unique value of x_col
    # the bars are colored by the unique values of fill_col
    
    if (is.null(types)) {
        types = unique(df[[fill_col]])
    }
    ntypes = length(types)
    rainbow_colors = rainbow(ntypes)
    class_colors = setNames(rainbow_colors, sort(types))

    # If x_order is provided, reorder x_col according to x_order
    if (!is.null(x_order)) {
        df[[x_col]] = factor(df[[x_col]], levels = x_order)
    } else {
        df[[x_col]] = factor(df[[x_col]])
    }


    if (plot_title == "") {
        plot_title = paste("Stacked bar plot of", y_col, "by", x_col, "colored by", fill_col)
    }

    if (x_axis_label == "") {
        x_axis_label = x_col
    }

    if (y_axis_label == "") {
        y_axis_label = y_col
    }

    p=ggplot(df, aes_string(x = x_col, y = y_col, fill = fill_col)) +
        geom_bar(stat = "identity", position = "stack", width = 0.7) + # stacked bar
        labs(title = plot_title, x = x_axis_label, y = y_axis_label) +
        theme(
            plot.title = element_text(size=16),
            axis.line = element_line(color = "black"),  # Add axis lines
            axis.ticks = element_line(color = "black"),  # Add axis ticks
            axis.text = element_text(size = 14),  # Increase tick label font size
            axis.title = element_text(size = 15),  # Increase axis label font size
           axis.text.x = element_text(angle = 45, hjust = 1)  
        ) +
        scale_fill_manual(values=class_colors)

    print(p)
}




selectColorScale = function(cmap, is_numeric, clim = c(0, 1), clab='', num_factors=NULL, reverse=FALSE) {
    if (cmap == "viridis") {
        if (is_numeric) {
            return(scale_color_viridis_c(name = clab, limits = clim))
        } else {
            return(scale_color_viridis_d(name = clab, begin=clim[1], end=clim[2]))
        }
    } else {
        base_palette = brewer.pal(11, cmap)
        color_ramp_func = colorRampPalette(base_palette)
        full_palette = color_ramp_func(101)
        if (reverse){full_palette = rev(full_palette)}
        
        if (is_numeric) {
            return(scale_color_gradientn(name = clab, colors = full_palette))
        } else {
            # Discrete scale using a portion of continuous colormap
            discrete_colors <- full_palette[seq(from = floor(clim[1] * 100), to = ceiling(clim[2] * 100), length.out = num_factors)]
            return(scale_color_manual(name = clab, values = discrete_colors))
        }
    }
}

setwh = function(w, h){
    options(repr.plot.width = w, repr.plot.height = h)
}

violin_plot_plus_hbar = function(
    sobj, features, group, assay, pt.size=0, f=mean, nrow=2, ncol=5, barsize=12) {
        plots <- lapply(features, function(gene) {
            p <- VlnPlot(
                sobj, 
                features = gene, 
                group.by = group, 
                assay = assay,  
                pt.size = pt.size
            ) + stat_summary(fun = f, geom = "point", size = barsize, colour = "black", shape = 95)

            # if the feature is not last in the list, add NoLegend(). Otherwise, return the plot as is
            if (gene != features[length(features)]) {
                p = p + NoLegend()
            }
            
            return(p)
            })

        # Combine the plots with wrap_plots
        combined_plot = wrap_plots(plots, nrow = nrow, ncol = ncol)

        # Display the combined plot
        print(combined_plot)

    }

plot_overlapping_density_histogram = function(
    df, 
    hist_col,
    fill_col,
    colors = c("blue", "red"),
    alpha=0.5,
    breaks=seq(0, 16, 1),
    title= NULL,
    xlab = NULL,
    fig_filename = NULL,
    show_legend = TRUE

){
    # hist_col is the column you're making a histogram of (e.g. nUMI)
    # fill_col is the column you're coloring by (e.g. cell_class)
    # if fig_filename is not null, the plot will be saved to that file
     
    if (is.null(xlab)){
        xlab = hist_col
    }

    if (is.null(title)){
        title = paste0("Density Histogram of ", xlab, " by ", fill_col)
    }


    p = (
        ggplot(df, aes_string(x=hist_col, fill=fill_col)) 
        + geom_histogram(aes(y=..density..), alpha=alpha, position="identity", breaks=breaks)
        + labs(title=title, x=xlab, y="Density")    
        + theme(
                plot.title = element_text(size=16),
                axis.line = element_line(color = "black"),  # Add axis lines
                axis.ticks = element_line(color = "black"),  # Add axis ticks
                axis.text = element_text(size = 14),  # Increase tick label font size
                axis.title = element_text(size = 15)  # Increase axis label font size
            ) 
        + scale_fill_manual(values=colors)   
    )

    if (!is.null(fig_filename)){
        ggsave(fig_filename, p, width=8, height=6)
    }

    if (!show_legend){
        p = p + theme(legend.position = "none")
    }

    return(p)
}

volcano_plot = function(
    data, 
    gene_col="gene",
    x_col="logFC", 
    y_col="negative_log10_p_adjusted", 
    sig_col = "significant", 
    x_label=expression(log[2]*"(Fold Change)"),
    y_label=expression(log[10]*"(adjusted p-value)"),
    title = NULL,
    annot_genes = NULL,
    xlim=c(-3,3)){

        # todo: add a vector argument for annot_x_pos and annot_y_pos

        if (is.null(title)){
            title = paste0("Volcano Plot\nFraction of Expressed Genes Significant: ", round(sum(data[[sig_col]])/nrow(data), 2))
        }

    volcano_plot = ggplot(data, aes_string(x = x_col, y = y_col)) +
        geom_point(aes_string(color = sig_col)) +
        scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
        scale_x_continuous(limits = xlim) +
        geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "blue") +
        geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "blue") +
        labs(
            title = title,
            x = x_label,
            y = y_label
        ) +
        theme(
            plot.title = element_text(size = 16),
            axis.line = element_line(color = "black"),  # Add axis lines
            axis.ticks = element_line(color = "black"),  # Add axis ticks
            axis.text = element_text(size = 14),  # Increase tick label font size
            axis.title = element_text(size = 15),  # Increase axis label font size
            axis.text.x = element_text(angle = 45, hjust = 1) 
        )

    if (!is.null(annot_genes)){
        annot_genes$nudge_x = 2*sign(annot_genes[[x_col]])
        volcano_plot = volcano_plot + geom_label_repel(
            data = annot_genes, 
            aes_string(label = gene_col),
            #vjust = -1,
            box.padding = 0.5,  # Increase this value to push labels further out
            point.padding = 0.5,  # Increase this value to add more space around the point
            nudge_y = 0.25,  # Move labels further in the y direction
            nudge_x = annot_genes$nudge_x,  # Move labels further in the x direction
            segment.color = 'grey50',
            max.overlaps = Inf
        )
    }

    # Print the plot
    print(volcano_plot)
}