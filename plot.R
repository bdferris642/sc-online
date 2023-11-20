source("~/sc-online/utils.R")

library(caret)
library(dplyr)
library(ggplot2)
library(grid)
library(patchwork)
library(pheatmap)
library(RColorBrewer)
library(rhdf5)
library(rlang)
library(tidyr)
library(viridis)
library(viridisLite)


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
  plot_df <- long_df %>%
    mutate(
      bin_str = cut(plot_col, breaks = breaks, include.lowest = TRUE), # the bins will be strings like '[0.05,0.1)'
      bin = sapply(as.character(bin_str), function(x) { # extract the number of the lower bound
          lbound_str = strsplit(x, ",")[[1]][[1]]
          return(as.numeric(substr(lbound_str, 2, nchar(lbound_str))))
    })) %>%
    group_by(group_col, bin) %>%
    summarize(count = n()) %>%
    ungroup() %>%
    group_by(group_col) %>%
    mutate(probability = count / sum(count)) %>%
    ungroup()

  plot_df$bin = plot_df$ bin + bin_width/2
  # Plotting the normalized histograms
  p = ggplot(plot_df, aes(x = bin, y = probability, fill = group_col)) +
    geom_bar(stat = "identity", position = "identity", alpha = 0.5, width = bin_width) +
    xlim(xlim) +
    theme_minimal() +
    labs(title = title, x = "Value", y = "Probability")
  
  print(p)
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

            p = ggplot(cd, aes_string(x=x_col, y=y_col, color=color_col)) +
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

#     options(repr.plot.width=16, repr.plot.height=12)

#     xlim=c(6, 17)
#     ylim=c(0, 1)

#     # For every 9 dataframes, create a 3x3 plotting layout
#     #par(mfrow = c(3, 3))
#     for (i in seq(1, n, 9)) {  
#     plot_list = list()
#     for (j in 0:8) {
#         if (i+j > n) break
        
#         name = names(donor_ids_list)[i+j]    
#         sce = sce_list[[name]]
#         cd = colData(sce)
#         donor_ids = donor_ids_list[[name]]
#         cd$log_nRead = log(cd$nRead)

#         cd = as.data.frame(cd)
#         cd$prob_donor = donor_ids$prob_max
#         cd = cd[cd$prob_doublet < (1/3),]
            
#         p = ggplot(cd, aes(x=log_nRead, y=prob_donor, color=pct_intronic)) +
#             geom_point(alpha=0.5) +
#             scale_color_viridis() +
#             facet_grid((i+j)%%3 ~ (i+j)%/%3) + 
#             ggtitle(paste(name, '\nVireo Max P(Donor) vs log(nRead)\nOverall Frac P(Donor) > 0.9: ', round(as.numeric(frac_gte_90_list[name]), 2))) + 
#             coord_cartesian(xlim=xlim, ylim=ylim)

        
#         # # color legend on top right
#         if (j==5){
#             p = p + theme(legend.position = "right")
#         } else {
#             p = p + theme(legend.position = "none")
#         }
        
#         plot_list[[j+1]] = p

#     }
#     do.call(grid.arrange, plot_list)
#     }
# }


plotStackedBarByGroup=function(
    df, cat_col, type_col, normalize=FALSE,
    plot_width=8, plot_height=8,
     x_axis_label = '', y_axis_label='', plot_title='',
     class_colors=NULL,
     fig_filename=NULL){

    # TODO: separate the concerns of preparing the data and plotting it
    # takes in a long df with categorical cat_col and type_col
    # plots stacked bar plots, one for each type_col, with the height of each sub-bar
    # dictated by counts of each catetgory in cat_col


    # make a table
    t_list = list()
    
    for (i in unique(df[[cat_col]])) {
        this_cat_types = df[df[, cat_col] == i, ][[type_col]]
        t = table(this_cat_types)
        if(normalize){
            t = t/length(this_cat_types)
        }
        t[cat_col] = i
        t_list = append(t_list, list(t))
    }

    # Convert each table to a data frame and add a group column
    df_list = lapply(1:length(t_list), function(i) {
        data.frame(type = names(t_list[[i]]), 
                    value = as.numeric(t_list[[i]]), 
                    category = paste(cat_col, '-', i))
    })

    # Combine all data frames into one
    combined_data = bind_rows(df_list)
    df2plot=combined_data[!is.na(combined_data$value),]
    
    # Define a color palette for classes
    if (is.null(class_colors)){
        u = unique(df2plot$type) 
        ntypes = length(u)
        rainbow_colors = rainbow(ntypes)
        class_colors = setNames(rainbow_colors, sort(u))
    }

    options(repr.plot.width=plot_width, repr.plot.height=plot_height)
    p=ggplot(df2plot, aes(x = category, y = value, fill = type)) +
        geom_bar(stat = "identity", position = "stack", width = 0.7) + # stacked bar
        labs(title = plot_title, x = x_axis_label, y = y_axis_label) +
        theme_minimal() +
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
