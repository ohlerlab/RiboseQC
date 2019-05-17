# Ribo-seQC, a comprehensive Ribo-seq quality control tool
#
# Authors: 
# Lorenzo Calviello (calviello.l.bio@gmail.com)
# Dominique Sydow (dominique.sydow@posteo.de)
# Dermott Harnett (Dermot.Harnett@mdc-berlin.de)
# Uwe Ohler (Uwe.Ohler@mdc-berlin.de)
#
# This software is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this software. If not, see
# <http://www.gnu.org/licenses/>.



######################################################################################################
# Methods for RiboseQC report
######################################################################################################

#' Create the Ribo-seQC analysis report in html
#'
#' This function creates the Ribo-seQC html report based on the Ribo-seQC analysis files 
#' generated with \code{RiboseQC_analysis}.
#' 
#' @keywords Ribo-seQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#' 
#' @param input_files Character vector with full paths to data files 
#' generated with \code{RiboseQC_analysis}. 
#' Must be of same length as \code{input_sample_names}.
#' 
#' @param input_sample_names Character vector containing input names (max. 5 characters per name).
#' Must be of same length as \code{input_files}.
#' 
#' @param output_file String; full path to html report file.
#' @param extended creates a large html report including codon occupancy for each read length. Defaults to \code{FALSE}
#' @return The function saves the html report file with the file path \code{output_file} and 
#' a folder containing all figures shown in the html report as RDS object files (located in the same directory as the html report).
#' 
#' @details This function creates the html report visualizing the RiboseQC analysis data. \cr \cr
#' Input are two lists of the same length: \cr \cr
#' a) \code{input_files}: list of full paths to one or multiple input files 
#' (Ribo-seQC analysis files generated with \code{RiboseQC_analysis}) and \cr \cr
#' b) \code{input_sample_names}: list of corresponding names describing the file content in max. 5 characters 
#' (these are used as names in the report). \cr \cr
#' For the report, a RMarkdown file is rendered as html document, saved as \code{output_file}. \cr \cr
#' Additionally, all figures in the report are saved as PDF figures 
#' in an extra folder in the same directory as the report html file. \cr \cr
#' Example: \cr
#' \code{output_file <- "\\mydir\\myreport.html"} will generate 
#' the html report \code{\\mydir\\myreport.html} and
#' the folder \code{\\mydir\\myreport_plots\\} for the RDS object files to be stored in.
#' 
#' @seealso \code{\link{plot_read_biotype_dist_1}}, \code{\link{plot_read_biotype_dist_2}},
#' \code{\link{plot_read_length_dist}}, \code{\link{plot_read_length_dist_by_biotype}},
#' \code{\link{plot_read_biotype_dist_by_length}}, 
#' \code{\link{get_metagene_data}}, 
#' \code{\link{plot_metagene_hm_rmd}}, \code{\link{plot_metagene_hm}},
#' \code{\link{plot_metagene_bar_rmd}}, \code{\link{plot_metagene_bar}}, 
#' \code{\link{plot_frame_dist_boxplot_rmd}}, \code{\link{plot_frame_dist_boxplot}}, 
#' \code{\link{get_rl_and_cutoffs}}, \code{\link{get_default_rl_selection}},
#' \code{\link{get_top50_mapping}}, \code{\link{get_top50_cds_genes}}, \code{\link{get_top50_all_genes}},
#' \code{\link{get_codon_usage_data}},
#' \code{\link{plot_codon_usage_positional_rmd}}, \code{\link{plot_codon_usage_positional}},
#' \code{\link{plot_codon_usage_bulk_rmd}}, \code{\link{plot_codon_usage_bulk}}
#' 
#' @export

create_html_report <- function(input_files, input_sample_names, output_file,extended=F){
    
    # get input and output file paths
    input_files <- paste(normalizePath(dirname(input_files)),basename(input_files),sep="/")
    output_file <- paste(normalizePath(dirname(output_file)),basename(output_file),sep="/")
    
    # get path to RMarkdown file (to be rendered)
    rmd_path <- paste(system.file(package="RiboseQC"),"/rmd/riboseqc_template.Rmd",sep="")
    if(extended){
        rmd_path <- paste(system.file(package="RiboseQC"),"/rmd/riboseqc_template_full.Rmd",sep="")
    }
    
    # get folder path for pdf figures
    output_fig_path <- paste(output_file,"_plots/", sep = "")
    
    # create folder for rds ojects and pdf figures
    dir.create(paste0(output_fig_path, "rds/"), recursive=TRUE, showWarnings=FALSE)
    dir.create(paste0(output_fig_path, "pdf/"), recursive=TRUE, showWarnings=FALSE)
    sink(file = paste(output_file,"_report_text_output.txt",sep = ""))
    # render RMarkdown file > html report
    suppressWarnings(render(rmd_path, 
           params = list(input_files = input_files,
                         input_sample_names = input_sample_names,
                         output_fig_path = output_fig_path),
           output_file = output_file))
    sink()
}


#' Generate PDF files from RDS object files
#'
#' This function generates figures as PDF files from RDS object files. 
#' 
#' @keywords Ribo-seQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#' 
#' @param output_rds_path String; full path to output folder for RDS object files. Example: /my_path_to/rds/
#' 
#' @return This function creates PDF files from RDS object files.
#' 
#' @seealso \code{\link{create_html_report}}
#' 
#' @export

create_pdfs_from_rds_objects <- function(output_rds_path){
    rds_dir <- paste0(output_rds_path)
    rds_files <- list.files(rds_dir)
    
    for (rds_file in rds_files) {
        # read RDS data
        rds <- readRDS(paste0(output_rds_path, rds_file))
        
        # get plot data
        g <- rds[[1]]
        
        # get format info
        format <- rds[[2]]
        
        # make sure a pdf directory exists on same level as rds directory
        dir.create(paste0(output_rds_path, "../pdf/"), recursive=TRUE, showWarnings=FALSE)
        
        # some objects can be directly printed, others need as_ggplot (arrangeGrob-derived objects)
        if (class(g)[2] == "ggplot") {
            pdf(paste0(output_rds_path, "../pdf/", rds_file, ".pdf"),
                width=format[1], height=format[2])
            print(g)
            dev.off()
        }
        else {
            pdf(paste0(output_rds_path, "../pdf/", rds_file, ".pdf"),
                width=format[1], height=format[2])
            print(as_ggplot(g))
            dev.off()
        }
    }
}


#' Generate a list of R data objects
#'
#' This function generates a list of loaded RData files 
#' to be used during Ribo-seQC report generation.
#' 
#' @keywords Ribo-seQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#' 
#' @param input_files List of RData file paths generated by \code{RiboseQC_analysis}. \cr \cr
#' Example: \cr
#' \code{input_files <- c(sample1="//path//to//sample1", sample2="//path//to//sample2")}
#' 
#' @return This function returns a list of loaded RData objects 
#' that were generated by \code{RiboseQC_analysis}.
#' 
#' @seealso \code{\link{create_html_report}}
#' 
#' @export

generate_rdata_list <- function(input_files){
    rdata_list <- list()
    for (i in 1:length(input_files)){
        load(input_files[i])
        rdata_list[[names(input_files)[i]]] <- res_all
    }
    return(rdata_list)
}


#' Plot read location distribution by biotype (and originating compartment)
#'
#' This function plots the read location distribution by biotype (and originating compartment)
#' for one input sample. \cr \cr
#' This plot is used in the Ribo-seQC report in section 1.1. 
#' 
#' @keywords Ribo-seQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#' 
#' @param pos \code{res_all$read_stats$positions} generated by \code{RiboseQC_analysis}. \cr \cr
#' \code{data.frame} containing the number of reads per biotype (rows) and originating compartment (columns).
#' 
#' @param sample String; sample name (selected from the input names 
#' given in the \code{input_sample_names} parameter of \code{create_html_report}).
#' 
#' @param output_rds_path String; full path to output folder for RDS object files created by this function.
#' Defaults to NOT save RDS; to save RDS, provide path to destination folder.
#' 
#' @return This function returns a plot that can be integrated in the html report and
#' that can be saved as RDS object file.
#' 
#' @seealso \code{\link{create_html_report}}
#' 
#' @export

plot_read_biotype_dist_1 <- function(pos, sample, output_rds_path=""){
    
    names(pos) <- gsub("reads_", "", names(pos))
    pos <- cbind(biotype=row.names(pos), pos)
    pos.m <- suppressMessages(melt(pos))
    pos.m$biotype <- factor(pos.m$biotype, levels=row.names(pos))
    pos.m$variable <- factor(pos.m$variable, levels=names(pos))
    head(pos.m)
    
    g2 <- ggplot(data = pos.m, aes(x=biotype, y=value, fill=variable)) + 
        geom_bar(stat="identity", position=position_stack(reverse=TRUE)) + 
        labs(x="", y="read count\n") +
        scale_fill_brewer(palette = "Set2") +
        theme_minimal() +
        theme(plot.margin = margin(.5, .5, .5, .5, "cm"),
              panel.background=element_blank(), 
              axis.ticks=element_blank(),
              axis.text.x=element_text(angle=90, hjust=1, size=14), 
              axis.text.y=element_text(size=14),
              axis.title=element_text(size=14),
              legend.position="bottom",
              legend.text=element_text(size=14),
              legend.title=element_blank()) 
    
    # print to HTML file
    print(g2)
    
    cat("\n\n")  # necessary for RMarkdown rendering
    
    # save plot as RDS object
    if(!output_rds_path==""){
        file_path = paste0(output_rds_path, sample, "_1_readlocdist")
        saveRDS(list(g2, c(10,5)), file = file_path, 
                ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
    }
}


#' Plot read location distribution by originating compartment (and biotype)
#' 
#' This function plots the read location distribution by originating compartment (and biotype)
#' for all input samples. \cr \cr 
#' This plot is used in the Ribo-seQC report in section 1.2.
#'
#' @keywords Ribo-seQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#' 
#' @param rdata_list List of RiboseQC analysis RData objects generated by \code{generate_rdata_list}.
#' 
#' @param output_rds_path String; full path to output folder for RDS object files created by this function.
#' Defaults to NOT save RDS; to save RDS, provide path to destination folder.
#' 
#' @return This function returns a plot that can be integrated in the html report and
#' that can be saved as RDS object file.
#' 
#' @seealso \code{\link{create_html_report}}
#' @export

plot_read_biotype_dist_2 <- function(rdata_list, output_rds_path=""){
    levels_bio <- NULL
    levels_comp <- NULL
    levels_sample <- names(rdata_list)
    
    data <- list()
    for (i in names(rdata_list)){
        res_all <- rdata_list[[i]]
        rs <- res_all$read_stats$reads_summary
        levels_bio <- row.names(rs[[1]])
        levels_comp <- names(rs)
        data[[i]] <- do.call(cbind, lapply(rs, function(x) rowSums(as.data.frame(x))))
    }
    
    data.m <- suppressMessages(melt(data))
    data.m$Var1 <- factor(data.m$Var1, levels=levels_bio)
    data.m$Var2 <- factor(data.m$Var2, levels=levels_comp)
    data.m$L1 <- factor(data.m$L1, levels=levels_sample)
    
    # read count
    g1 <- ggplot(data=data.m, aes(x=L1, y=value, fill=Var1)) +
        geom_bar(stat="identity", position=position_stack(reverse=TRUE)) +
        scale_fill_viridis(discrete=T) +
        labs(x="\nsample", y="read count\n") +
        facet_grid(Var2~., scales="free_y") +
        theme_minimal() +
        theme(plot.margin = margin(.3, .3, .3, .3, "cm"),
              panel.background=element_blank(), 
              axis.ticks.x=element_blank(),
              axis.text.x=element_text(angle=90, hjust=1, size=14), 
              axis.text.y=element_text(size=14),
              axis.title=element_text(size=14),
              legend.text=element_text(size=14),
              legend.title=element_blank(),
              strip.text.y = element_text(size=14)) +
        guides(fill=guide_legend(nrow=1))
    
    # read count percentage
    g2 <- ggplot(data=data.m, aes(x=L1, y=value, fill=Var1)) +
        geom_bar(stat="identity", position=position_fill(reverse=TRUE)) +
        scale_fill_viridis(discrete=T) +
        labs(x="\nsample", y="read count fraction\n") +
        facet_grid(Var2~., scales="free_y") +
        theme_minimal() +
        theme(plot.margin = margin(.3, .3, .3, .3, "cm"),
              panel.background=element_blank(), 
              axis.ticks.x=element_blank(),
              axis.text.x=element_text(angle=90, hjust=1, size=14), 
              axis.text.y=element_text(size=14),
              axis.title=element_text(size=14),
              legend.text=element_text(size=14),
              legend.title=element_blank(),
              strip.text.y = element_text(size=14)) +
        guides(fill=guide_legend(nrow=1))
    
    gg <- ggarrange(g1, g2, ncol = 2, common.legend=TRUE, legend="bottom") +
        theme(plot.margin = margin(.5, .5, .5, .5, "cm"))
    
    # print to HTML file
    print(gg)
    
    cat("\n\n")  # necessary for RMarkdown rendering
    
    # save plot as RDS object
    if(!output_rds_path==""){
        file_path = paste0(output_rds_path, "all_samples_1_readlocdist")
        saveRDS(list(gg, c(10,7)), file = file_path, 
                ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
    }
}



#' Plot read length distribution 
#' 
#' This function plots the read length distribution per originating compartment
#' for one input sample. \cr \cr
#' This plot is used in the Ribo-seQC report in section 2.
#' 
#' @keywords Ribo-seQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#' 
#' @param rld \code{res_all$read_stats$rld} generated by \code{RiboseQC_analysis} \cr \cr
#' \code{data.frame} containing the number of reads per originating compartment (rows) and read lengths (columns).
#' 
#' @param sample String; sample name (selected from the input names  
#' given in the \code{input_sample_names} parameter of \code{create_html_report}).
#' 
#' @param output_rds_path String; full path to output folder for RDS object files created by this function.
#' Defaults to NOT save RDS; to save RDS, provide path to destination folder.
#' 
#' @return This function returns a plot that can be integrated in the html report and
#' that can be saved as RDS object file.
#' 
#' @seealso \code{\link{create_html_report}}
#' @export

plot_read_length_dist <- function(rld, sample, output_rds_path=""){
    lw <- 1.5 # line width
    
    rld <- as.matrix(rld)
    colnames(rld) <- as.integer(gsub("reads_", "", colnames(rld)))
    
    levels_rl <- as.integer(colnames(rld))
    levels_comp <- rownames(rld)
    
    # read count
    rld.m <- suppressMessages(melt(rld))
    rld.m$Var1 <- factor(rld.m$Var1, levels=levels_comp)
    #rld.m$Var2 <- factor(rld.m$Var2, levels=levels_rl)
    rld.m$Var2 <- as.numeric(rld.m$Var2)
    #colssi<-alpha(rainbow(length(unique(rld.m$Var1))),alpha = .7)
    colssi<-alpha(c("blue","dark red","forestgreen",rainbow(length(unique(rld.m$Var1))-3)),alpha = .7)
    
    
    
    g1 <- ggplot(data=rld.m, aes(x=Var2, y=value)) + 
        geom_line(aes(colour=Var1, group=Var1), size=lw) + 
        scale_x_continuous(breaks=scales::pretty_breaks(n=8)) +
        labs(x="\nread length", y="read count\n") +
        scale_color_manual(values = colssi) +
        theme_minimal() +
        theme(plot.margin = margin(.3, .3, .3, .3, "cm"),
              axis.text=element_text(size=14),
              axis.title=element_text(size=14),
              legend.text=element_text(size=14),
              legend.title=element_blank()) 
    
    # read count fraction
    rld_norm <- rld/rowSums(rld)
    rld_norm.m <- suppressMessages(melt(rld_norm))
    rld_norm.m$Var1 <- factor(rld_norm.m$Var1, levels=levels_comp)
    #rld_norm.m$Var2 <- factor(rld_norm.m$Var2, levels=levels_rl)
    #colssi<-alpha(rainbow(length(unique(rld.m$Var1))),alpha = .7)
    colssi<-alpha(c("blue","dark red","forestgreen",rainbow(length(unique(rld.m$Var1))-3)),alpha = .7)
    
    g2 <- ggplot(data=rld_norm.m, aes(x=Var2, y=value)) + 
        geom_line(aes(colour=Var1, group=Var1), size=lw) + 
        scale_x_continuous(breaks=scales::pretty_breaks(n=8)) +
        labs(x="\nread length", y="read count fraction\n") +
        scale_color_manual(values = colssi) +
        theme_minimal() +
        theme(plot.margin = margin(.3, .3, .3, .3, "cm"),
              axis.text=element_text(size=14),
              axis.title=element_text(size=14),
              legend.text=element_text(size=14),
              legend.title=element_blank()) 
    
    gg <- ggarrange(g1, g2, ncol = 2, common.legend=TRUE, legend="bottom") +
        theme(plot.margin = margin(.5, .5, .5, .5, "cm"))
    
    # print to HTML file
    print(gg)
    
    cat("\n\n")  # necessary for RMarkdown rendering
    
    # save plot as RDS object
    if(!output_rds_path==""){
        file_path = paste0(output_rds_path, sample, "_2_readlendist")
        saveRDS(list(gg, c(10,5)), file = file_path, 
                ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
    }
    
}



#' Plot read length and location distribution 
#' (distribution per biotype)
#' 
#' This funtion plots the read length distribution for each originating compartment as distribtion per biotype
#' for one input sample (displayed as read count and as read count fraction). \cr \cr
#' This plot is used in the Ribo-seQC report in section 3.1.
#' 
#' @keywords Ribo-seQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#' 
#' @param reads_summary \code{res_all$read_stats$reads_summary} generated by \code{RiboseQC_analysis} \cr \cr
#' List of DataFrames:
#' one DataFrame for each originating compartment,
#' each DataFrame contains read counts per biotype (rows) and read lengths (columns).
#' 
#' @param sample String; sample name (selected from the input names  
#' given in the \code{input_sample_names} parameter of \code{create_html_report}).
#' 
#' @param output_rds_path String; full path to output folder for RDS object files created by this function.
#' Defaults to NOT save RDS; to save RDS, provide path to destination folder.
#' 
#' @return This function returns a plot that can be integrated in the html report and 
#' that can be saved as RDS object file.
#' 
#' @seealso \code{\link{create_html_report}}
#' @export

plot_read_length_dist_by_biotype <- function(reads_summary, sample, output_rds_path=""){
    lw <- 1.5 # line width
    
    levels_bio <- row.names(reads_summary[[1]])
    levels_comp <- names(reads_summary)
    levels_rl <- as.integer(gsub("reads_", "", names(reads_summary[[1]])))
    
    rld <- lapply(reads_summary, function(x) as.matrix(x))
    
    # read count
    rld.m <- suppressMessages(melt(rld))
    rld.m$Var1 <- factor(rld.m$Var1, levels=levels_bio)
    rld.m$Var2 <- as.integer(gsub("reads_", "", rld.m$Var2))
    #rld.m$Var2 <- factor(rld.m$Var2, levels=levels_rl)
    rld.m$L1 <- factor(rld.m$L1, levels=levels_comp)
    
    g1 <- ggplot(data=rld.m, aes(x=Var2, y=value)) + 
        geom_line(aes(colour=Var1, group=Var1), size=lw) + 
        scale_x_continuous(breaks=scales::pretty_breaks(n=8)) +
        labs(x="\nread length", y="read count\n") +
        scale_color_viridis(discrete=T,alpha = .8) +
        facet_grid(L1~., scales="free_y") +
        theme_minimal() +
        theme(plot.margin = margin(.3, .3, .3, .3, "cm"),
              axis.text=element_text(size=14),
              axis.title=element_text(size=14),
              legend.text=element_text(size=14),
              legend.title=element_blank(),
              strip.text.y = element_text(size=14)) +
        guides(colour=guide_legend(nrow=1))
    
    # read count fraction
    rld_norm <- list()
    for (i in names(rld)){
        rld_norm[[i]] <- rld[[i]]/(rowSums(rld[[i]]))
        rld_norm[[i]][is.na(rld_norm[[i]])] = 0 # replace NaN fields with 0
    }
    rld_norm.m <- suppressMessages(melt(rld_norm))
    rld_norm.m$Var1 <- factor(rld_norm.m$Var1, levels=levels_bio)
    rld_norm.m$Var2 <- as.integer(gsub("reads_", "", rld_norm.m$Var2))
    #rld_norm.m$Var2 <- factor(rld_norm.m$Var2, levels=levels_rl)
    rld_norm.m$L1 <- factor(rld_norm.m$L1, levels=levels_comp)
    
    g2 <- ggplot(data=rld_norm.m, aes(x=Var2, y=value)) + 
        geom_line(aes(colour=Var1, group=Var1), size=lw) + 
        scale_x_continuous(breaks=scales::pretty_breaks(n=8)) +
        labs(x="\nread length", y="read count fraction\n") +
        scale_color_viridis(discrete=T,alpha = .8) +
        facet_grid(L1~., scales="free_y") +
        theme_minimal() +
        theme(plot.margin = margin(.3, .3, .3, .3, "cm"),
              axis.text=element_text(size=14),
              axis.title=element_text(size=14),
              legend.text=element_text(size=14),
              legend.title=element_blank(),
              strip.text.y = element_text(size=14)) +
        guides(colour=guide_legend(nrow=1))
    
    gg <- ggarrange(g1, g2, ncol = 2, common.legend=TRUE, legend="bottom") +
        theme(plot.margin = margin(.5, .5, .5, .5, "cm"))
    
    # print to HTML file
    plot(gg)
    
    cat("\n\n")  # necessary for RMarkdown rendering
    
    # save plot as RDS object
    if(!output_rds_path==""){
        file_path = paste0(output_rds_path, sample, "_3_readloclendist1")
        saveRDS(list(gg, c(10,8)), file = file_path, 
                ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
    }
}



#' Plot read length and location distribution 
#' (distribution per read length)
#' 
#' This funtion plots the biotype distribution for each originating compartment as distribtion per read length
#' for one input sample (displayed as read count and as read count fraction). \cr \cr
#' This plot is used in the Ribo-seQC report in section 3.2.
#' 
#' @keywords Ribo-seQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#' 
#' @param reads_summary \code{res_all$read_stats$reads_summary} generated by \code{RiboseQC_analysis} \cr \cr
#' List of DataFrames:
#' one DataFrame for each originating compartment,
#' each DataFrame contains read counts per biotype (rows) and read lengths (columns).

#' @param sample String; sample name (selected from the input names  
#' given in the \code{input_sample_names} parameter of \code{create_html_report}).
#' 
#' @param output_rds_path String; full path to output folder for RDS object files created by this function.
#' Defaults to NOT save RDS; to save RDS, provide path to destination folder.
#' 
#' @return This function returns a plot that can be integrated in the html report and
#' that can be saved as RDS object file.
#' 
#' @seealso \code{\link{create_html_report}}
#' @export

plot_read_biotype_dist_by_length <- function(reads_summary, sample, output_rds_path=""){
    levels_bio <- row.names(reads_summary[[1]])
    levels_comp <- names(reads_summary)
    levels_rl <- as.integer(gsub("reads_", "", names(reads_summary[[1]])))
    
    rld <- lapply(reads_summary, function(x) as.matrix(x))
    
    # read count
    rld.m <- suppressMessages(melt(rld))
    rld.m$Var1 <- factor(rld.m$Var1, levels=levels_bio)
    rld.m$Var2 <- as.integer(gsub("reads_", "", rld.m$Var2))
    #rld.m$Var2 <- factor(rld.m$Var2, levels=levels_rl)
    rld.m$L1 <- factor(rld.m$L1, levels=levels_comp)
    
    #ggplot(data=rld.m, aes(x=Var2, y=value, fill=Var1)) +
    #    geom_bar(stat="identity")
    
    
    # read count
    g1 <- ggplot(data=rld.m, aes(x=Var2, y=value, fill=Var1)) +
        geom_bar(stat="identity", position=position_stack(reverse=TRUE)) +
        scale_x_continuous(breaks=scales::pretty_breaks(n=8)) +
        scale_fill_viridis(discrete=T) +
        labs(x="\nread length (nt)", y="read count\n") +
        facet_grid(L1~., scales="free_y") +
        theme_minimal() +
        theme(plot.margin = margin(.3, .3, .3, .3, "cm"),
              panel.background=element_blank(), 
              #axis.ticks.x=element_blank(),
              axis.text.x=element_text(angle=90, vjust=0.5, size=14), 
              axis.text.y=element_text(size=14),
              axis.title=element_text(size=14),
              legend.text=element_text(size=14),
              legend.title=element_blank(),
              strip.text.y = element_text(size=14)) +
        guides(fill=guide_legend(nrow=1))
    
    # read count percentage
    g2 <- ggplot(data=rld.m, aes(x=Var2, y=value, fill=Var1)) +
        geom_bar(stat="identity", position=position_fill(reverse=TRUE)) +
        scale_x_continuous(breaks=scales::pretty_breaks(n=8)) +
        scale_fill_viridis(discrete=T) +
        labs(x="\nread length (nt)", y="read count fraction\n") +
        facet_grid(L1~., scales="free_y") +
        theme_minimal() +
        theme(plot.margin = margin(.3, .3, .3, .3, "cm"),
              panel.background=element_blank(), 
              #axis.ticks.x=element_blank(),
              axis.text.x=element_text(angle = 90, vjust=0.5, size=14), 
              axis.text.y=element_text(size=14),
              axis.title=element_text(size=14),
              legend.text=element_text(size=14),
              legend.title=element_blank(),
              legend.position = "bottom",
              strip.text.y = element_text(size=14)) +
        guides(fill=guide_legend(nrow=1))
    
    gg <- ggarrange(g1, g2, ncol = 2, common.legend=TRUE, legend="bottom") +
        theme(plot.margin=margin(.5, .5, .5, .5, "cm"))
    
    # print to HTML file
    print(gg)
    
    cat("\n\n")  # necessary for RMarkdown rendering
    
    # save plot as RDS object
    if(!output_rds_path==""){
        file_path = paste0(output_rds_path, sample, "_3_readloclendist2")
        saveRDS(list(gg, c(10,8)), file = file_path, 
                ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
    }
}


#' Get 5'/P-site profile data for metagene analysis
#' 
#' This function processes profile data generated by \code{RiboseQC_analysis}. \cr \cr
#' This data is used as input in \code{plot_metagene_hm} to generate plot for the Ribo-seQC report in section 4.1/4.3.
#'
#' @keywords Ribo-seQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#' 
#' @param profiles 5' or P-site profile data generated by \code{RiboseQC_analysis}: \cr \cr
#' \code{res_all$profiles_fivepr} or \cr
#' \code{res_all$profiles_P_sites} \cr \cr
#' Consists of DataFrames each 
#' containing counts of 5' or P-site profiles, 
#' calculated for different resolution types (see parameter \code{res}), 
#' originating compartments (see parameter \code{comp}), and
#' read lengths per input sample. \cr \cr
#' Example to access DataFrame: \cr 
#' \code{res_all$profiles_fivepr[[res]][[comp]][[read_length]]} or \cr
#' \code{res_all$profiles_P_sites[[res]][[comp]][[read_length]]}
#' 
#' @param comp String for originating compartment \cr \cr
#' Check for available originating compartments in the data set using:
#' \code{names(res_all$profiles_fivepr$five_prime_subcodon)} 
#'  
#' @param res Resolution (subcodon or bins)  \cr \cr
#' \bold{Subcodon} resolution: \cr
#' \code{five_prime_subcodon} \cr 
#' (in order to call \code{res_all$profiles_fivepr$five_prime_subcodon}) or \cr 
#' \code{P_sites_subcodon} \cr 
#' (in order to call \code{res_all$profiles_P_sites$P_sites_subcodon}) \cr \cr
#' Read coverage for 
#' the first 25nt after the transcription start site (TSS), 
#' 25nt before and 33nt after the start codon, 
#' 33nt from the middle of the CDS, 
#' 33nt before and 25nt after the stop codon, 
#' and the last 25nt before the transcription end site (TES). \cr \cr
#' \bold{Bins}: \cr
#' \code{five_prime_bins} \cr
#' (in order to call \code{res_all$profiles_fivepr$five_prime_bins}) or \cr
#' \code{P_sites_bins} \cr
#' (in order to call \code{res_all$profiles_P_sites$P_sites_bins}) \cr \cr
#' Read coverage for 
#' 50 bins between TSS and start codon, 100 bins for the CDS, and 50 after stop codon to TES.
#' 
#' @return This function returns data \code{profile_data} 
#' as \code{list(data_single, data_all, res)} with profile data
#' \itemize{
#' \item for all read lengths individually (\code{profile_data[1]]} is \code{data_single}) and
#' \item for all read lenghts summarized (\code{profile_data[[2]]} is \code{data_all}).
#' }
#' with different scaling: no scaling (\code{none}), log2 scaling (\code{log2}), and z scoring (\code{zscore}), 
#' accessable via e.g. \code{profile_data[[1]]$none} or \code{profile_data[[2]]$zscore}. \cr \cr
#' 
#' \code{profile_data[[3]]} saves the resolution type \code{res} for later use during plotting. \cr \cr
#' 
#' \code{profile_data[[4]]} stores information on whether data represents 5' or P site profiles (\code{profile_type}) 
#' for later use during plotting.
#' 
#' @seealso \code{\link{create_html_report}}
#' @export

get_metagene_data <- function(data, profile_type, res, comp){
    
    # test input parameters for existence
    if((profile_type != "profiles_fivepr") & (profile_type != "profiles_P_sites")){
        stop("Variable profile_type takes only: \"profiles_fivepr\" or \"profiles_P_sites\"")
    }
    if(profile_type == "profiles_fivepr") {
        if((res != "five_prime_bins") & (res != "five_prime_subcodon")){
            stop("Variable res takes only: \"five_prime_bins\" or \"five_prime_subcodon\"")
        }
    }
    if(profile_type == "profiles_P_sites") {
        if((res != "P_sites_bins") & (res != "P_sites_subcodon")){
            stop("Variable res takes only: \"P_sites_bins\" or \"P_sites_subcodon\"")
        }
    }
    if(!(comp %in% names(data[[profile_type]][[res]]))) {
        stop("Variable res takes only: ", paste(names(data[[profile_type]][[res]]), collapse=", "))
    }
    
    # get metagene data
    signal <- data[[profile_type]][[res]][[comp]]
    # get available read lengths
    rl_ok <- names(signal)
    # sort read lengths in ascending order and "all" at the beginning
    rl_ok <- c(rl_ok[1], as.character(sort(as.integer(rl_ok[-1], decreasing=F))))
    
    # prepare BARPLOTS data
    
    data_bar = list()
    
    for (rl in rl_ok) {
        # get read length data and format
        x <- signal[[rl]]
        x.m <- suppressMessages(melt(x))
        x.m$pos <- 1:length(x)
        x.m$frame <- paste0("frame ", rep(c(3,1,2),67)[1:length(x)])
        
        # add read length data to list
        data_bar[[rl]] <- x.m
    }
    
    # prepare HEATMAP data
    
    # we want to make two separate heatmaps for 
    # a) single read lenghts (single) and 
    # b) summarized read lengths (all).
    
    # single read lengths
    rl_ok_single <- rl_ok[-1]
    rl_ok_single <- c(as.character(sort(as.integer(rl_ok_single), decreasing = F)))
    data_none <- do.call(signal[rl_ok_single], what=rbind)
    rownames(data_none) <- rl_ok_single
    colnames(data_none) <- 1:dim(data_none)[2]
    data_log2 <- log2(data_none+1)
    #modified to direct scaling
    #data_zscore <- t(scale(t(data_none), center=TRUE, scale=TRUE)) # scale centers (true) and scales (true) the columns of a numeric matrix, we want to scale per read length > transpose data_none to have rl as column
    data_zscore <- scale(data_none, center=TRUE, scale=TRUE) # scale centers (true) and scales (true) the columns of a numeric matrix, we want to scale per read length > transpose data_none to have rl as column
    data <- list(none=suppressMessages(melt(t(data_none))), log2=suppressMessages(melt(t(data_log2))), zscore=suppressMessages(melt(t(data_zscore))))
    data$none$Var2 <- factor(data$none$Var2,levels=rl_ok_single)
    data$log2$Var2 <- factor(data$log2$Var2,levels=rl_ok_single)
    data$zscore$Var2 <- factor(data$zscore$Var2,levels=rl_ok_single)
    data_hm_single <- data
    
    # read length "all"
    rl_ok_all <- rl_ok[1]
    data_none <- do.call(signal[rl_ok_all], what=rbind)
    rownames(data_none) <- names(signal[rl_ok_all])
    colnames(data_none) <- 1:dim(data_none)[2]
    data_log2 <- log2(data_none+1)
    #it needs the transpose as it's a 1-row matrix
    data_zscore <- t(scale(t(data_none), center=TRUE, scale=TRUE))
    data <- list(none=suppressMessages(melt(t(data_none))), log2=suppressMessages(melt(t(data_log2))), zscore=suppressMessages(melt(t(data_zscore))))
    data$none$Var2 <- data$none$Var2
    data$log2$Var2 <- data$log2$Var2
    data$zscore$Var2 <- data$zscore$Var2
    data_hm_all <- data
    
    # save barplot, heatmap and meta data!
    metagene_data <- list(data_bar=data_bar, 
                          data_hm_single=data_hm_single, 
                          data_hm_all=data_hm_all, 
                          profile_type=profile_type,
                          res=res, 
                          comp=comp)
    
    return(metagene_data)
}



#' Plot 5'/P-site profiles as heatmaps within the RMarkdown document
#' 
#' This function generates iteratively all heatmap plots 
#' (iteration over originating compartments, resolution types, and scaling methods). \cr \cr
#' These plots are displayed in the Ribo-seQC report in section 4.1/4.3.
#' 
#' @keywords Ribo-seQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#' 
#' @param profiles 5' or P-site profile data generated by \code{RiboseQC_analysis}: \cr \cr
#' \code{res_all$profiles_fivepr} or \cr
#' \code{res_all$profiles_P_sites} \cr \cr
#' Consists of DataFrames each 
#' containing counts of 5' or P-site profiles, 
#' calculated for different resolution types (see parameter \code{res}), 
#' originating compartments (see parameter \code{comp}), and
#' read lengths per input sample. \cr \cr
#' Example to access DataFrame: \cr 
#' \code{res_all$profiles_fivepr[[res]][[comp]][[read_length]]} or \cr
#' \code{res_all$profiles_P_sites[[res]][[comp]][[read_length]]}
#' 
#' @param sample String; sample name (selected from the input names 
#' given in the \code{input_sample_names} parameter of \code{create_html_report}).
#' 
#' @param output_rds_path String; full path to output folder for RDS object files created by this function.
#' Defaults to NOT save RDS; to save RDS, provide path to destination folder.
#' 
#' @return This function returns iteratively all 5' or P-site profile plots 
#' for the html report and saves the same plots as RDS object file.
#' 
#' @seealso \code{\link{create_html_report}}, \code{\link{get_metagene_data}}, \code{\link{plot_metagene_hm}}
#' @export

plot_metagene_hm_rmd <- function(data, profile_type, sample="", output_rds_path=""){
    
    # split data by resolution (subcodon vs. bins)
    profile <- data[[profile_type]]
    ress <- rev(names(profile))
    ress <- c(ress[grepl("subcodon", ress)], ress[grepl("bins", ress)])
    
    comps <- names(profile[[1]])
    
    for (comp in comps){
        # cat commands necessary for RMarkdown rendering
        cat("\n\n")
        cat("#### ", comp, " {.tabset} \n\n")
        cat("Select a *resolution* (subcodon or bins):\n")
        cat("In case of *subcodon* resolution, read coverage is shown for
            the first 25nt after the transcription start site (TSS), 
            25nt before and 33nt after the start codon, 
            33nt from the middle of the CDS,
            33nt before and 25nt after the stop codon, 
            and the last 25nt before the transcription end site (TES).\n")
        cat("In case of *bins*, read coverage is shown for 50 bins between TSS and start codon (5'UTR), 
            100 bins for the CDS, and 50 after stop codon (3'UTR). \n\n")
        cat("Select a *scaling* method: *none* (no scaling), *log2* scaling, and *z-score* scaling.\n\n")
        cat("**Note**: Select *subcodon / none* or *bins / none* scaling to display the read coverage for each read length separately.\n\n")
        
        for(res in ress){
            
            # get heatmap data
            metagene_data <- get_metagene_data(data, profile_type, res, comp)
            
            # get different scaling factors (none, log2, z-score)
            scals <- names(metagene_data$data_hm_single)
            
            for (scal in scals){
                cat("\n\n")  # necessary for RMarkdown rendering
                
                # get resolution-scale combo, e.g. subcodon / log2
                cat("##### ", strsplit(res, split="_")[[1]][3], "/", scal, " {.tabset} \n\n")  # necessary for RMarkdown rendering
                
                # plot heatmap (displaying single read lengths and summarized read lengths (all))
                plot_metagene_hm(metagene_data, scal, sample, output_rds_path)
            }
        }
    }
}


#' Plot 5'/P-site profiles (for all read lengths) as heatmap
#' 
#' This function plots a 5' or P-site profiles for all read lengths as heatmap
#' (for a specific originating compartment, resolution type, and scaling method). \cr \cr
#' This plot is used in the Ribo-seQC report in section 4.1/4.3.
#'
#' @keywords Ribo-seQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#' 
#' @param data_profile Profile data for a specific originating compartment and resolution type,
#' generated using \code{\link{get_metagene_data}}.
#' 
#' @param scal Scaling method:
#' no scaling (\code{none}), log2 scaling (\code{log2}), or z scoring (\code{zscore}).
#' 
#' @param sample String; sample name (selected from the input names  
#' iven in the \code{input_sample_names} parameter of \code{create_html_report}).
#' 
#' @param output_rds_path String; full path to output folder for RDS object files created by this function.
#' Defaults to NOT save RDS; to save RDS, provide path to destination folder.
#' 
#' @return This function returns a plot that can be integrated in the html report and  #' that can be saved as RDS object file.
#' 
#' @seealso \code{\link{create_html_report}}
#' @export

plot_metagene_hm <- function(metagene_data, scal, sample="", output_rds_path="") {
    
    # heatmap data
    data_hm_single <- metagene_data$data_hm_single
    data_hm_all <- metagene_data$data_hm_all
    
    # metadata
    profile_type <- metagene_data$profile_type
    res <- metagene_data$res
    comp <- metagene_data$comp
    
    # set subcodon/binned resolution-specific parameters
    if (grepl("subcodon", res)) {
        res_breaks <- c(1, 26, 51, 84, 117, 150, 176, 200)
        res_labels <- c("TSS", "", "start\ncodon", "", "", "stop\ncodon", "", "TES")
        res_x_lab <- "\nposition (nucleotide resolution)"
        
    } else {
        res_breaks <- c(1, 51, 150, 200)
        res_labels <- c("TSS", "start\ncodon", "stop\ncodon", "TES")
        res_x_lab <- "\nbins"
    }
    
    # heatmap for summarized read lengths (all)
    g1 <- ggplot(data_hm_all[[scal]], aes(x=Var1, y=Var2)) + 
        geom_tile(aes(fill=value)) +
        scale_fill_gradient(low="white", high="steelblue", na.value="gray80") +
        scale_x_continuous(breaks=res_breaks, labels=res_labels) +
        scale_y_discrete(expand=c(0, 0)) +
        theme(axis.text=element_text(size=18), 
              axis.title.y=element_text(size=18, colour="white"),
              axis.ticks.y=element_blank(),
              panel.background=element_blank(), 
              axis.title.x=element_blank(), 
              axis.ticks.x=element_blank(), 
              axis.text.x=element_blank(),
              legend.position = "top", 
              legend.direction = "horizontal",
              legend.title=element_text(size=18),
              legend.text=element_text(size=14)) +
        labs(colour="r", y="lg\n", x=res_x_lab) +
        geom_vline(xintercept=c(51,150), lty=3, color="black") +
        geom_point(data=data_hm_all[[scal]], aes(size="NA"), shape=NA, colour="gray80") +
        guides(fill=guide_legend(title=paste0("read count (scale: ", scal, ")"), nrow=1),
               size=guide_legend(title=""))    
    
    data_hm_single[[scal]]$Var2<-as.numeric(as.character(data_hm_single[[scal]]$Var2))
    # heatmap for single read lenghts (single)
    g2 <- ggplot(data_hm_single[[scal]], aes(x=Var1, y=Var2)) + 
        geom_tile(aes(fill=value)) +
        scale_fill_gradient(low="white", high="steelblue", na.value="gray80") +
        scale_x_continuous(breaks=res_breaks, labels=res_labels) +
        scale_y_continuous(breaks=scales::pretty_breaks(n=10)) +
        #scale_y_discrete(expand=c(0, 0)) +
        theme(axis.text=element_text(size=18), 
              axis.title=element_text(size=18),
              #axis.ticks.y=element_blank(),
              panel.background=element_blank(),
              legend.position = "bottom", 
              legend.direction = "horizontal",
              legend.title=element_text(size=18),
              legend.text=element_text(size=14,angle = 45)) +
        labs(colour="r", y="read length (nt)\n", x=res_x_lab) +
        geom_vline(xintercept=c(51,150), lty=3, color="black") +
        geom_point(data=data_hm_single[[scal]], aes(size="NA"), shape=NA, colour="gray80") +
        labs(fill=paste0("read count (scale: ", scal, ")")) + guides(size=guide_legend(title="")) 
    
    gg <- arrangeGrob(g1, g2, ncol=1, heights = c(1, 6))
    
    suppressMessages(if(scal=="zscore"){gg<-arrangeGrob(g1 + scale_fill_gradient2(low = "blue",mid = "white",high = "firebrick",na.value = "gray80",guide = "colourbar"),
                    g2 + scale_fill_gradient2(low = "blue",mid = "white",high = "firebrick",na.value = "gray80",guide = "colourbar"),
                    ncol=1, heights = c(1, 6))})
    
    
    
    cat("\n\n")  # necessary for RMarkdown rendering
    
    # save plot as RDS object
    if(!output_rds_path==""){
        file_path = paste0(output_rds_path, 
                           sample, "_", 
                           comp, "_",
                           "4_", profile_type, "_",
                           "metagene_", unlist(strsplit(res, "_"))[3], "_", 
                           scal)
        saveRDS(list(gg, c(12,10)), file = file_path, 
                ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
    }
    
    # print to HTML file
    print(as_ggplot(gg))
    
    # generate per-read-length metagene plots if scaling is none
    if(scal=="none" & grepl("subcodon", res)){
        plot_metagene_bar_rmd(metagene_data, sample, output_rds_path)
    }
}


#' Plot 5'/P-site profile per read length as barplot in R Markdown
#' 
#' This function plots a 5' or P-site profile as barplot
#' (for a specific originating compartment, resolution type, and read length). \cr \cr
#' This plot is used in the Ribo-seQC report in section 4.1/4.3.
#'
#' @keywords Ribo-seQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#' 
#' @param metagene_data ...
#' 
#' 
#' @param sample String; sample name (selected from the input names given in the \code{input_sample_names} parameter of \code{create_html_report}).
#' 
#' @param output_rds_path String; full path to output folder for RDS object files created by this function. Defaults to NOT save RDS; to save RDS, provide path to destination folder.
#' 
#' @return ...
#' 
#' @seealso \code{\link{create_html_report}}
#' @export

plot_metagene_bar_rmd <- function(metagene_data, sample="", output_rds_path="") {
    
    data_bar <- metagene_data$data_bar
    profile_type <- metagene_data$profile_type
    rls <- names(data_bar) # read lengths
    
    # first, check which read lengths you really want to display:
    # only keep read lengths with a cutoff coverage
    rls_selected = rls
    
    # do this only for the 5' profiles 
    # (for P-site profiles, we want to see everything that is there, except empty read lengths)
    # profile_type 1 refers to the 5' profiles!
    if(profile_type=="profiles_fivepr"){ 
        # set cutoff (coverage per read length)
        cutoff <- 100
        # calculate coverage per read length
        rls_cov <- unlist(lapply(lapply(data_bar, "[", , 1), sum))[rls]
        # select read lengths according to cutoff
        rls_selected <- names(rls_cov[rls_cov > cutoff])
    }
    else{
        # set cutoff (coverage per read length)
        cutoff <- 0
        # calculate coverage per read length
        rls_cov <- unlist(lapply(lapply(data_bar, "[", , 1), sum))[rls]
        # select read lengths according to cutoff
        rls_selected <- names(rls_cov[rls_cov > cutoff])
    }
    
    for (rl in rls_selected){
        cat("\n\n")
        cat("###### ", rl, " \n\n")
        
        plot_metagene_bar(metagene_data, rl, sample, output_rds_path)
    }
}


#' Plot 5'/P-site profile per read length as barplot
#' 
#' This function plots a 5' or P-site profile as barplot
#' (for a specific originating compartment, resolution type, and read length). \cr \cr
#' This plot is used in the Ribo-seQC report in section 4.1/4.3.
#'
#' @keywords Ribo-seQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#' 
#' @param data Profile data for a resolution type, specific originating compartment, and read length. \cr \cr
#' Example: \cr
#' \code{\strong{res_all$profiles_fivepr}$five_prime_subcodon$nucl$'30'} \cr or \cr
#' \code{\strong{res_all$profiles_P_sites}$five_prime_subcodon$nucl$'30'} \cr \cr
#' Shown in bold are fixed list names, remaining list names should be adapted to 
#' the resolution type, specific originating compartment, and read length of interest.
#' 
#' 
#' @param sample String; sample name (selected from the input names  
#' iven in the \code{input_sample_names} parameter of \code{create_html_report}).
#' 
#' @param output_rds_path String; full path to output folder for RDS object files created by this function.
#' Defaults to NOT save RDS; to save RDS, provide path to destination folder.
#' 
#' @return This function returns a plot that can be integrated in the html report and that can be saved as RDS object file.
#' 
#' @seealso \code{\link{create_html_report}}
#' @export

plot_metagene_bar <- function(metagene_data, rl, sample="", output_rds_path="") {
    
    # bar data
    data_bar <- metagene_data$data_bar
    # applied to selected read length
    x.m <- data_bar[[rl]]
    
    # metadata
    profile_type <- metagene_data$profile_type
    res <- metagene_data$res
    comp <- metagene_data$comp
    
    # set subcodon/binned resolution-specific parameters
    if(grepl("subcodon", res)) {
        res_breaks <- c(1, 26, 51, 84, 117, 150, 176, 200)
        res_labels <- c("TSS", "", "start\ncodon", "", "", "stop\ncodon", "", "TES")
        res_x_lab <- "\nposition (nucleotide resolution)"
        
    } else {
        res_breaks <- c(1, 51, 150, 200)
        res_labels <- c("TSS", "start\ncodon", "stop\ncodon", "TES")
        res_x_lab <- "\nbins"
    }
    
    # we need to treat subcodon and binned data differently:
    # subcodon data has frame information, therefore frames should be colored accordingly
    # binned data has NO frame information, therefore colour everything in black, please, dear computer
    
    # plot data
    if(grepl("subcodon", res)){
        
        x.m$frame<-as.character(x.m$frame)
        x.m$frame[x.m$frame=="frame 2"]<-"Frame 0"
        x.m$frame[x.m$frame=="frame 1"]<-"Frame 2"
        x.m$frame[x.m$frame=="frame 3"]<-"Frame 1"
        x.m$frame<-gsub(x.m$frame,pattern = "F",replacement = "f")
        x.m$frame<-factor(x.m$frame,levels=c("frame 0","frame 1","frame 2"))
        
        g3 <- ggplot(data=x.m, aes(x=pos, xend=pos, y=0, yend=value)) +
            geom_segment(aes(colour=frame, group=frame)) +
            scale_color_manual(values = c("red","forestgreen","blue")) +
            scale_x_continuous(breaks=res_breaks, labels=res_labels) +
            #scale_y_continuous(expand=c(0, 0)) +
            labs(x=res_x_lab, y="read count\n") +
            theme_minimal() +
            theme(legend.position = "bottom",
                  legend.text = element_text(size = 18),
                  legend.title = element_blank(),
                  axis.text = element_text(size = 18),
                  axis.title = element_text(size = 18))
    } else {
        g3 <- ggplot(data=x.m, aes(x=pos, xend=pos, y=0, yend=value)) +
            geom_segment(aes(group=frame), colour="black") +
            scale_x_continuous(breaks=res_breaks, labels=res_labels) +
            #scale_y_continuous(expand=c(0, 0)) +
            labs(x=res_x_lab, y="read count\n") +
            theme_minimal() +
            theme(legend.position = "bottom",
                  legend.text = element_text(size = 18),
                  legend.title = element_blank(),
                  axis.text = element_text(size = 18),
                  axis.title = element_text(size = 18))
    }
    
    # print to HTML file
    print(g3)
    
    cat("\n\n")  # necessary for RMarkdown rendering
    
    # save plot as RDS object
    if(!output_rds_path==""){
        file_path = paste0(output_rds_path, sample, "_", 
                           comp, "_",
                           "4_", profile_type, "_",
                           "metagene_", 
                           unlist(strsplit(res, "_"))[3], "_", 
                           rl)
        saveRDS(list(g3, c(12,10)), file = file_path, 
                ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
    }
}


#' Plot frame coverage of 5'-site profiles within the RMarkdown document
#' 
#' This function generates iteratively all frame coverage boxplots
#' (iteration over originating compartments). \cr \cr
#' These plots are displayed in the Ribo-seQC report in section 4.2.1.
#'
#' @keywords Ribo-seQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#' 
#' @param analysis_frame_cutoff Object containing statistics on frame coverage
#' (per originating compartment and read length) generated by \code{RiboseQC_analysis}. \cr \cr
#' Example: \cr
#' \code{res_all$selection_cutoffs$analysis_frame_cutoff}
#' 
#' @param sample String; sample name (selected from the input names  
#' given in the \code{input_sample_names} parameter of \code{create_html_report}).
#' 
#' @param output_rds_path String; full path to output folder for RDS object files created by this function.
#' Defaults to NOT save RDS; to save RDS, provide path to destination folder.
#' 
#' @return This function integrates plots in the html report.
#' 
#' @seealso \code{\link{create_html_report}}
#' @export

plot_frame_dist_boxplot_rmd <- function(analysis_frame_cutoff, sample="", output_rds_path=""){
    comps <- names(analysis_frame_cutoff)
    # only compartments with data
    comps <- comps[unlist(lapply(comps, function(x) length(analysis_frame_cutoff[[x]])))!=0]
    
    for (comp in comps){
        cat("\n\n")  # necessary for RMarkdown rendering
        cat("##### ", comp, " {.tabset} \n\n")  # necessary for RMarkdown rendering
        
        plot_frame_dist_boxplot(analysis_frame_cutoff, comp, sample, output_rds_path)
    }
}



#' Plot frame coverage of 5'-site profiles
#' 
#' This function plots the frame coverage of 5'-site profiles, 
#' i.e. the fraction of reads (their 5' end) in frame 0, 1, and 2 (defined by start codon),
#' as boxplots (per-frame distributions over all transcripts)
#' for individual read lengths as well as for all read lengths summarized. 
#' 
#' This plot is used in the Ribo-seQC report in section 4.2.1.
#' 
#' @keywords Ribo-seQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#' 
#' @param analysis_frame_cutoff Object containing statistics on frame coverage
#' (per originating compartment and read length) generated by \code{RiboseQC_analysis}. \cr \cr
#' Example: \cr
#' \code{res_all$selection_cutoffs$analysis_frame_cutoff}
#' 
#' @param comp String for originating compartment \cr \cr
#' Check for available originating compartments in the data set using:
#' \code{names(res_all$profiles_fivepr$five_prime_subcodon)}
#' 
#' @param sample String; sample name (selected from the input names  
#' given in the \code{input_sample_names} parameter of \code{create_html_report}).
#' 
#' @param output_rds_path String; full path to output folder for RDS object files created by this function.
#' Defaults to NOT save RDS; to save RDS, provide path to destination folder.
#' 
#' @return This function returns a plot that can be integrated in the html report and 
#' that can be saved as RDS object file.
#' 
#' @seealso \code{\link{create_html_report}}
#' @export

plot_frame_dist_boxplot <- function(analysis_frame_cutoff, comp, sample="", output_rds_path="") {
    rl <- names(analysis_frame_cutoff[[comp]])
    rl_list <- NULL
    rl_list_mean <- NULL
    
    for (k in rl){
        data <- analysis_frame_cutoff[[comp]][[k]]$frames
        colnames(data) <- c("frame 0", "frame 1", "frame 2")
        rl_list[[k]] <- data
        rl_list_mean[[k]] <- colMeans(data)[1]
    }
    data <- suppressMessages(melt(rl_list))
    
    ## order read lengths
    data$L1 <- factor(data$L1, levels=sort(names(analysis_frame_cutoff[[comp]])))
    # or by fraction of reads in maximal frame
    # data$L1 <- factor(data$L1, levels=c(row.names(res_all$selection_cutoffs$results_choice[[comp]]$data), "all")) 
    
    #from boxplot to stat_summary
    # gg <- ggplot(data, aes(x = L1, y = value, fill=Var2)) +
    #     geom_boxplot() +
    #     #stat_summary(fun.y = mean, geom="point",colour="magenta", size=2) +
    #     scale_x_discrete(name="read length") +
    #     scale_y_continuous(name="read distribution between frames") +
    #     theme(legend.position = "bottom", panel.background=element_blank()) +
    #     labs(fill = "frame") +
    #     guides(fill=guide_legend(title=""))
    
    gg <- ggplot(data, aes(x = L1, y = value, color=Var2, group=Var2)) +
        stat_summary(fun.data = mean_se,size=1.5) +
        #stat_summary(fun.y = mean, geom="point",colour="magenta", size=2) +
        scale_x_discrete(name="read length") +
        scale_color_manual(values = alpha(c("red","forestgreen","blue"),alpha = .8)) +
        scale_y_continuous(name="read distribution between frames",limits=c(0,1)) +
        theme(legend.position = "bottom", panel.background=element_blank()) +
        labs(color = "frame") +
        guides(fill=guide_legend(title=""))
    
    
    # print to HTML file
    print(gg)
    
    cat("\n\n")  # necessary for RMarkdown rendering
    
    # save plot as RDS object
    if(!output_rds_path==""){
        file_path = paste0(output_rds_path, sample, "_", comp, "_4_selection5toP")
        saveRDS(list(gg, c(7,5)), file = file_path, 
                ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
    }
}


#' Get selected read lengths and cutoffs
#' 
#' This function retrieves data on selected read lengths (per originating compartment),
#' e.g. their cutoff values, frame preference and codon gain. \cr \cr
#' This data is used in the Ribo-seQC report in section 4.2.2 (displayed as table).
#'
#' @keywords Ribo-seQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#' 
#' @param rdata_list List of RiboseQC analysis RData objects generated by \code{generate_rdata_list}.
#' 
#' @return This function returns data.
#' 
#' @seealso \code{\link{create_html_report}}
#' @export

get_rl_and_cutoffs <- function(rdata_list) {
    datasets <- NULL
    for (i in c(1:length(rdata_list))) {
        res_all <- rdata_list[[names(rdata_list)[i]]]
        
        sample <- names(rdata_list)[[i]]
        comps <- names(res_all$selection_cutoffs$results_choice)
        
        for (comp in comps){
            d <- as.data.frame(res_all$selection_cutoffs$results_choice[[comp]]$data)[c(3, 1, 2, 4, 5)]
            names(d) <- c("read length", "cutoff", "frame preference (%)", "gain_codon?", "gain_new_codons?")
            datasets[[sample]][[comp]] <- d
        }
    }
    return(datasets)
}



#' Get default choice of read lengths 
#' 
#' This function returns selected read lengths per originating compartment. 
#' These read lengths build the basis for all P-site based calculations such as
#' metagene and codon usage analyses 
#' (e.g. \code{plot_metagene_hm} and \code{plot_codon_usage}) \cr \cr
#' This data is used in the Ribo-seQC report in section 4.2.3 (displayed as table).
#'
#' @keywords Ribo-seQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#' 
#' @param rdata_list List of RiboseQC analysis RData objects generated by \code{generate_rdata_list}.
#' 
#' @return This function returns data.
#' 
#' @seealso \code{\link{create_html_report}}
#' @export

get_default_rl_selection <- function(rdata_list){
    datasets <- NULL
    for (i in 1:length(rdata_list)){
        res_all <- rdata_list[[names(rdata_list)[i]]]
        sample <- names(rdata_list)[[i]]
        comps <- names(res_all$selection_cutoffs$results_choice)
        
        rl_choices_default <- NULL
        for (j in comps){
            data <- res_all$selection_cutoffs$results_choice[[j]]$data
            data_default <- data[data[["max_coverage"]],][["read_length"]] # choose only read lengths with max_coverage=TRUE (default choice)
            rl_choices_default <- rbind(rl_choices_default, c(compartment=j, 
                                                              `read lengths of choice (ordered by quality)`=paste(data_default, collapse = ", ")))
        }
        datasets[[sample]] <- rl_choices_default
    }
    return(datasets)
}



#' Get top 50 mapping positions
#' 
#' This function retrieves data on the top 50 positions (nucleotide resolution) are listed
#' where most reads (their 5' end) map to, revealing possibly contaminating sequences. \cr \cr
#' This data is used in the Ribo-seQC report in section 5 (displayed as table).
#'
#' @keywords Ribo-seQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#' 
#' @param rdata_list List of RiboseQC analysis RData objects generated by \code{generate_rdata_list}.
#' 
#' @return This function returns data to be displayed as table in the html report.
#' 
#' @seealso \code{\link{create_html_report}}
#' @export

get_top50_mapping <- function(rdata_list) {
    datasets <- list()
    for (i in c(1:length(rdata_list))){
        res_all <- rdata_list[[names(rdata_list)[i]]]
        data <- as.data.frame(res_all$sequence_analysis)
        data <- data[,c("score", "pct", "seqnames", "region", "gene_biotype", "gene", "seq")]
        names(data) <- c("read count", "library fraction (%)", "seqnames", "region", "gene biotype", "gene", "sequence")
        datasets[[names(rdata_list)[i]]] <- data
    }
    return(datasets)
}



#' Get top 50 abundant genes (CDS regions for protein coding genes)  
#' 
#' This function retrieves data on the top 50 abundant CDS regions for protein coding genes. \cr \cr
#' This data is used in the Ribo-seQC report in section 6 (displayed as table).
#'
#' @keywords Ribo-seQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#' 
#' @param rdata_list List of RiboseQC analysis RData objects generated by \code{generate_rdata_list}.
#' 
#' @return This function returns data to be displayed as table in the html report.
#' 
#' @seealso \code{\link{create_html_report}}
#' @export

get_top50_cds_genes <- function(rdata_list) {
    datasets <- list()
    for (i in c(1:length(rdata_list))){
        res_all <- rdata_list[[names(rdata_list)[i]]]
        rc_cds <- as.data.frame(res_all$read_stats$counts_cds_genes)
        rc_cds <- rc_cds[with(rc_cds, order(-reads, chromosome)), ][1:50,]
        datasets[[names(rdata_list)[i]]] <- rc_cds
    }
    return(datasets)
}

#' Get top 50 abundant genes (all genes) 
#'
#' This function retrieves data on the top 50 abundant genes. \cr \cr
#' This data is used in the Ribo-seQC report in section 6 (displayed as table).
#'
#' @keywords Ribo-seQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#' 
#' @param rdata_list List of RiboseQC analysis RData objects generated by \code{generate_rdata_list}.
#' 
#' @return This function returns data to be displayed as table in the html report.
#' 
#' @seealso \code{\link{create_html_report}}
#' @export

get_top50_all_genes <- function(rdata_list) {
    datasets <- list()
    for (i in c(1:length(rdata_list))){
        res_all <- rdata_list[[names(rdata_list)[i]]]
        rc_all <- as.data.frame(res_all$read_stats$counts_all_genes)
        rc_all <- rc_all[with(rc_all, order(-reads, chromosome)), ][1:100,]
        datasets[[names(rdata_list)[i]]] <- rc_all
    }
    return(datasets)
}


#' Get codon usage data (positional and bulk)
#' 
#' This function processes codon usage data generated by \code{RiboseQC_analysis}, i.e.
#' codon usage 
#' \itemize{
#' \item per nucleotide position (i.e. positional codon usage) as well as
#' \item summed up over all positions* (i.e. bulk codon usage)
#' }
#' for a specific data type, originating compartment, and read length, as well as based on a user-defined genetic code. \cr \cr
#' This data is used as input in \code{plot_codon_usage_bulk} 
#' to generate a bar plot, and in \code{plot_codon_usage_bulk_rmd} to iterativly generate plots for the Ribo-seQC report in section 7.2. \cr\cr
#' Please check \code{plot_codon_usage} for positional codon usage (instead of bulk codon usage). \cr\cr
#' * Information on positions included in analysis: \cr
#' Based on P-site corrected reads, the codon usage within CDS regions for protein coding genes 
#' is examplatory calculated 
#' \itemize{
#' \item for the first 11 codons of the CDS (referred to as \emph{start}), 
#' \item for 11 codons from the middle of the CDS (referred to as \emph{middle}), and 
#' \item for the last 11 codons of the CDS (referred to as \emph{stop}). 
#' }
#' 
#' @keywords Ribo-seQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#' 
#' @param data Object (list of lists) generated by \code{RiboseQC_analysis}: \code{res_all}.
#' 
#' @param data_type String; select one of the following: \cr
#' \itemize{
#' \item \code{Codon_counts}: codon count in defined positions*, \cr
#' \item \code{P_sites_percodon}: P-sites count in defined positions*, or \cr
#' \item \code{P_sites_percodon_ratio}: ratio of P-sites counts to codon counts in defined positions*.
#' \item \code{E_sites_percodon}: E-sites count in defined positions*, or \cr
#' \item \code{E_sites_percodon_ratio}: ratio of E-sites counts to codon counts in defined positions*.
#' \item \code{A_sites_percodon}: A-sites count in defined positions*, or \cr
#' \item \code{A_sites_percodon_ratio}: ratio of A-sites counts to codon counts in defined positions*.
#' }
#' 
#' @param comp String for originating compartment. \cr
#' Check for available originating compartments in the data set using: \cr
#' \code{names(res_all$profiles_P_sites$Codon_count)} 
#' 
#' @param rl String for read length. \cr
#' Check for available read lengths in the data set using: \cr
#' \code{names(res_all$profiles_P_sites[[data_type]][[comp]])} 
#' 
#' @return This function returns a list (e.g. called codon_usage_bulk_data) with information on bulk codon usage: \cr
#' \itemize{
#' \item \code{codon_usage_bulk_data$data} contains the data on bulk codon usage,
#' \item \code{codon_usage_bulk_data$data_type} saves the data type used (see parameter \code{data_type})
#' \item \code{codon_usage_bulk_data$comp} saves the originating compartment used (see parameter \code{comp}), and
#' \item \code{codon_usage_bulk_data$rl} saves the read length used (see parameter \code{rl}).
#' }
#' 
#' @seealso \code{\link{create_html_report}}
#' @export

get_codon_usage_data <- function(data, data_type, comp, rl) {
    
    # test input parameters for existence
    if((data_type != "Codon_counts") & (data_type != "P_sites_percodon") & (data_type != "P_sites_percodon_ratio") & (data_type != "A_sites_percodon") & (data_type != "A_sites_percodon_ratio") & (data_type != "E_sites_percodon") & (data_type != "E_sites_percodon_ratio")){
        stop("Variable data_type takes only: \"Codon_counts\" or \"P_sites_percodon\" or \"P_sites_percodon_ratio\" or \"A_sites_percodon\" or \"A_sites_percodon_ratio\" \"E_sites_percodon\" or \"E_sites_percodon_ratio\"")
    }
    if(!(comp %in% names(data$profiles_P_sites[[data_type]]))) {
        stop("Variable res takes only: ", paste(names(data$profiles_P_sites[[data_type]]), collapse=", "))
    }
    if(!(rl %in% names(data$profiles_P_sites[[data_type]][[comp]]))) {
        stop("Variable rl takes only: ", paste(names(data[["profiles_P_sites"]][[data_type]][[comp]]), collapse=", "))
    }
    
    # POSITIONAL CODON USAGE
    
    # get codon usage (per nucleotide position)
    d <- data[["profiles_P_sites"]][[data_type]][[comp]][[rl]]
    d <- as.data.frame(d)
    d[is.na(d)] <- 0 # set all NAs to zero
    names(d) <- seq(1,length(names(d)))
    
    # get row names (example row name: "GGG;G" corrsponding to codon;amino acid)
    row_names <- row.names(d)
    
    # get codons
    codons <- unlist(lapply(row_names, function(x) strsplit(x, ";")[[1]][1]))
    codons <- factor(codons, levels=rev(sort(codons)))
    
    # add amino acids
    gencode <- unlist(lapply(row_names, function(x) strsplit(x, ";")[[1]][2]))
    gencode <- factor(gencode, levels=c("M", 
                                        "G", "A", "V", "I", "L", "F", "Y", "W", 
                                        "C", "P",
                                        "S", "T", "N", "Q", 
                                        "R", "H", "K",
                                        "D", "E",
                                        "*"))
    
    # add scaling (none, log2, zscore)
    # none
    d_none <- cbind(aa=gencode, codon=codons, d)
    # log2
    d_log2 <- cbind(aa=gencode, codon=codons, log2(d+1))
    # zscore
    d_zscore <- cbind(aa=gencode, codon=codons, data.frame(scale(d, center=TRUE, scale=TRUE)))
    colnames(d_zscore) <- c("aa", "codon", seq(1,length(names(d)))) # need to set position names again
    # save all scaling types as list
    data_positional <- list(none=d_none, log2=d_log2, zscore=d_zscore)
    
    # BULK CODON USAGE
    
    # summarize codon usage over all nucleotide positions
    d_bulk <- rowSums(d)
    d_bulk <- as.data.frame(d_bulk)
    d_bulk <- cbind(aa=gencode, codon=codons, d_bulk)
    
    # save codon usage bulk data and metadata information
    codon_usage_data <- list(data_positional=data_positional,
                             data_bulk=d_bulk, 
                             data_type=data_type, 
                             comp=comp, 
                             rl=rl)
    
    return(codon_usage_data)
}


#' Plot positional codon usage heatmaps within the RMarkdown document
#' 
#' This function generates iteratively all positional codon usage heatmaps; 
#' iteration over originating compartments, read length, data type (codon count, read count, codon-read count ratio),
#' and scaling method (none, log2, zscale). \cr \cr
#' These plots are displayed in the Ribo-seQC report in section 7.
#'
#' @keywords Ribo-seQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#' 
#' @param data Object (list of lists) generated by \code{RiboseQC_analysis}: \code{res_all}.
#' 
#' @param sample String; sample name (selected from the input names 
#' given in the \code{input_sample_names} parameter of \code{create_html_report}).
#' 
#' @param output_rds_path String; full path to output folder for RDS object files created by this function.
#' Defaults to NOT save RDS; to save RDS, provide path to destination folder.
#' 
#' @return This function returns iteratively all positional codon usage plots 
#' for the html report and saves the same plots as RDS object file.
#' 
#' @seealso \code{\link{create_html_report}}
#' @export

plot_codon_usage_positional_rmd <- function(data, sample="", output_rds_path="") {
    
    # different data types for codon usage available: codon count, read count, codon-read ratio 
    data_types <- c("Codon counts"="Codon_counts", 
                    "A-sites counts"="A_sites_percodon", 
                    "P-sites counts"="P_sites_percodon", 
                    "E-sites counts"="E_sites_percodon", 
                    "A-sites per codon"="A_sites_percodon_ratio",
                    "P-sites per codon"="P_sites_percodon_ratio",
                    "E-sites per codon"="E_sites_percodon_ratio")
    
    # assumption that list structure is exactly the same for Codon_counts, P_sites_percodon, P_sites_percodon_ratio
    comps <- names(data[["profiles_P_sites"]][["Codon_counts"]]) 
    
    # per compartment
    for (comp in comps) {
        
        # cat commands necessary for RMarkdown rendering
        cat("\n\n")
        cat("### ", comp, " {.tabset} \n\n")
        cat("\n\n")
        rls <- as.vector(lapply(data[["profiles_P_sites"]][["Codon_counts"]], function(x) names(x))[[comp]])
        
        # per read length
        for (rl in rls) {
            cat("\n\n")
            cat("#### ", rl, " {.tabset} \n\n")
            cat("\n\n")
            
            # per codon count, read count, or codon-read ratio
            for(i in 1:length(data_types)){
                data_type <- data_types[i]
                cat("\n\n")
                cat("##### ", names(data_type), " {.tabset} \n\n")
                cat("\n\n")
                
                # get data
                codon_usage_data <- get_codon_usage_data(data, data_type, comp, rl)
                
                # per scaling method
                scals <- names(codon_usage_data$data_positional)
                for (scal in scals){
                    cat("\n\n")
                    cat("###### ", scal, " {.tabset} \n\n")
                    cat("\n\n")
                    
                    # plot positional codon usage
                    plot_codon_usage_positional(codon_usage_data, scal, sample, output_rds_path)
                }
            }
        }
    }
}


#' Plot positional codon usage heatmap
#' 
#' This function plots the codon usage per nucleotide position* (positional codon usage) as heatmap
#' for a specific data type, originating compartment, read length, and scaling method, as well as based on a user-defined genetic code. \cr \cr
#' * Information on positions included in analysis: \cr
#' Based on P-site corrected reads, the codon usage within CDS regions for protein coding genes 
#' is examplatory calculated 
#' \itemize{
#' \item for the first 11 codons of the CDS (referred to as \emph{start}), 
#' \item for 11 codons from the middle of the CDS (referred to as \emph{middle}), and 
#' \item for the last 11 codons of the CDS (referred to as \emph{stop}). 
#' }
#' 
#' @keywords Ribo-seQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#' 
#' @param codon_usage_bulk_data List containing codon usage bulk data and meta data, generated by \code{get_codon_usage_data}.
#' 
#' @param sample String; sample name (selected from the input names  
#' given in the \code{input_sample_names} parameter of \code{create_html_report}).
#' 
#' @param output_rds_path String; full path to output folder for RDS object files created by this function.
#' Defaults to NOT save RDS; to save RDS, provide path to destination folder.
#' 
#' @return This function returns a plot that can be integrated in the html report and  
#' that can be saved as RDS object file.
#' 
#' @seealso \code{\link{create_html_report}}
#' @export

plot_codon_usage_positional <- function(codon_usage_data, scal, sample="", output_rds_path="") {
    
    # get data and meta information
    d_positional <- codon_usage_data$data_positional[[scal]]
    data_type <- codon_usage_data$data_type
    comp <- codon_usage_data$comp
    rl <- codon_usage_data$rl
    
    # translate data_type into a string suitable for the y axis label
    data_types <- c("Codon counts"="Codon_counts", 
                    "A-sites counts"="A_sites_percodon", 
                    "P-sites counts"="P_sites_percodon", 
                    "E-sites counts"="E_sites_percodon", 
                    "A-sites per codon"="A_sites_percodon_ratio",
                    "P-sites per codon"="P_sites_percodon_ratio",
                    "E-sites per codon"="E_sites_percodon_ratio")
    data_type <- names(data_types[data_types==data_type])
    
    # since ATG contains much more data than the other codons, 
    # we will plot two heatmaps: ATG and non-ATG
    
    # heatmap for start and stop codons
    d <- d_positional[d_positional$aa=="*" | d_positional$codon=="ATG", ]  # get start and stop codons
    d.m <- suppressMessages(melt(d))
    #d.m$variable<-factor(rep(1:11,length(unique(d.m$codon))))
    d.m$position <- factor(c(rep("start", nrow(d.m)/3), 
                             rep("middle", nrow(d.m)/3), 
                             rep("stop", nrow(d.m)/3)), 
                           levels=c("start", "middle", "stop"))
    g1 <- ggplot(d.m, aes(x=variable, y=codon)) +
        geom_tile(aes(fill=value), colour="darkgrey") +
        facet_grid(aa~position, scales = "free", space = "free") +
        scale_fill_gradient(low="white", high="steelblue", na.value="gray80") +
        scale_y_discrete(expand=c(0, 0)) +
        theme(axis.text.y=element_text(size=18, family="mono"),
              axis.title.y=element_text(size=18),
              axis.text.x=element_blank(),
              axis.title.x=element_blank(),
              axis.ticks=element_blank(),
              panel.background=element_blank(),
              legend.position = "top",
              legend.direction = "horizontal",
              legend.title=element_text(size=18),
              legend.text=element_text(size=14),
              strip.text = element_text(size=18)) +
        labs(colour="r", y="", x="\nposition (codon resolution)") +
        geom_point(data=d.m, aes(size="NA"), shape=NA, colour="gray80") +
        guides(fill=guide_legend(title=paste0(data_type,"(scale: ", scal, ")"), nrow=1), 
               size=guide_legend(title=""))    
    # heatmap for NON start and stop codons
    d <- d_positional[!(d_positional$aa=="*" | d_positional$codon=="ATG"), ]  # remove start and stop codons
    d.m <- suppressMessages(melt(d))
    #d.m$variable<-factor(rep(1:11,length(unique(d.m$codon))))
    d.m$position <- factor(c(rep("start", nrow(d.m)/3), 
                             rep("middle", nrow(d.m)/3), 
                             rep("stop", nrow(d.m)/3)), 
                           levels=c("start", "middle", "stop"))
    g2 <- ggplot(d.m, aes(x=variable, y=codon)) +
        geom_tile(aes(fill=value), colour="darkgrey") +
        facet_grid(aa~position, scales="free", space="free") +
        scale_fill_gradient(low="white", high="steelblue", na.value="gray80") +
        scale_y_discrete(expand=c(0, 0)) +
        theme(axis.text.x=element_text(size=18),
              axis.text.y=element_text(size=18, family="mono"),
              axis.title=element_text(size=18),
              axis.ticks = element_blank(),
              panel.background=element_blank(),
              legend.position = "bottom",
              legend.direction = "horizontal",
              legend.title=element_text(size=18),
              legend.text=element_text(size=14,angle = 45),
              strip.text.y = element_text(size=18),
              strip.text.x = element_blank()) +
        labs(colour="r", y="", x="\nposition (codon resolution)") +
        geom_point(data=d.m, aes(size="NA"), shape=NA, colour="gray80") +
        labs(fill=paste0(data_type,"(scale: ", scal, ")")) + guides(size=guide_legend(title="")) 
    
    gg <- arrangeGrob(g1, g2, ncol=1, heights = c(1, 8))
    
    suppressMessages(if(scal=="zscore"){gg<-arrangeGrob(g1 + scale_fill_gradient2(low = "blue",mid = "white",high = "firebrick",na.value = "gray80",guide = "colourbar"),
                                       g2 + scale_fill_gradient2(low = "blue",mid = "white",high = "firebrick",na.value = "gray80",guide = "colourbar"),
                                       ncol=1, heights = c(1, 8))})
    
    cat("\n\n")  # necessary for RMarkdown rendering
    
    # save plot as RDS object
    if(!output_rds_path==""){
        file_path = paste0(output_rds_path, sample, "_", 
                           comp, "_",
                           "7_codonusage_positional_", 
                           gsub(" ", "_", data_type), "_",
                           rl, "_", 
                           scal)
        saveRDS(list(gg, c(12,16)), file = file_path, 
                ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
    }
    
    # print to HTML file
    print(as_ggplot(gg))
}


#' Plot bulk codon usage bar plots within the RMarkdown document
#' 
#' This function generates iteratively all bulk codon usage bar plots; 
#' iteration over originating compartments, read length, and data type (codon count, read count, codon-read count ratio). \cr \cr
#' These plots are displayed in the Ribo-seQC report in section 8.
#' 
#' @keywords Ribo-seQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#' 
#' @param data Object (list of lists) generated by \code{RiboseQC_analysis}: \code{res_all}.
#' 
#' @param sample String; sample name (selected from the input names 
#' given in the \code{input_sample_names} parameter of \code{create_html_report}).
#' 
#' @param output_rds_path String; full path to output folder for RDS object files created by this function.
#' Defaults to NOT save RDS; to save RDS, provide path to destination folder.
#' 
#' @return This function returns iteratively all bulk codon usage plots 
#' for the html report and saves the same plots as RDS object file.
#' 
#' @seealso \code{\link{create_html_report}}, \code{\link{get_metagene_data}}, \code{\link{plot_metagene_hm}}
#' @export

plot_codon_usage_bulk_rmd <- function(data, sample="", output_rds_path="") {
    
    # different data types for codon usage available: codon count, read count, codon-read ratio 
    data_types <- c("Codon counts"="Codon_counts", 
                    "A-sites counts"="A_sites_percodon", 
                    "P-sites counts"="P_sites_percodon", 
                    "E-sites counts"="E_sites_percodon", 
                    "A-sites per codon"="A_sites_percodon_ratio",
                    "P-sites per codon"="P_sites_percodon_ratio",
                    "E-sites per codon"="E_sites_percodon_ratio")
    
    # assumption that list structure is exactly the same for Codon_counts, P_sites_percodon, P_sites_percodon_ratio
    comps <- names(data[["profiles_P_sites"]][["Codon_counts"]]) 
    rls <- lapply(data[["profiles_P_sites"]][["Codon_counts"]], function(x) names(x))
    
    for (comp in comps) {
        
        # cat commands necessary for RMarkdown rendering
        cat("\n\n")
        cat("### ", comp, " {.tabset} \n\n")
        cat("\n\n")
        rr <- as.vector(lapply(data[["profiles_P_sites"]][["Codon_counts"]], function(x) names(x))[[comp]])
        
        for (rl in rr) {
            cat("\n\n")
            cat("#### ", rl, " {.tabset} \n\n")
            cat("\n\n")
            
            for(i in 1:length(data_types)){
                data_type <- data_types[i]
                cat("\n\n")
                cat("##### ", names(data_type), " {.tabset} \n\n")
                
                # get data
                codon_usage_data <- get_codon_usage_data(data, data_type, comp, rl)
                
                # plot data
                plot_codon_usage_bulk(codon_usage_data, sample, output_rds_path)
            }
        }
    }
}


#' Plot bulk codon usage bar plots
#' 
#' This function plots the codon usage summed up over all positions* (bulk codon usage) as bar plot
#' for a specific data type, originating compartment, and read length, as well as based on a user-defined genetic code. \cr \cr
#' * Information on positions included in analysis: \cr
#' Based on P-site corrected reads, the codon usage within CDS regions for protein coding genes 
#' is examplatory calculated 
#' \itemize{
#' \item for the first 11 codons of the CDS (referred to as \emph{start}), 
#' \item for 11 codons from the middle of the CDS (referred to as \emph{middle}), and 
#' \item for the last 11 codons of the CDS (referred to as \emph{stop}). 
#' }
#' 
#' @keywords Ribo-seQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#' 
#' @param codon_usage_bulk_data List containing codon usage bulk data and meta data, generated by \code{get_codon_usage_data}.
#' 
#' @param sample String; sample name (selected from the input names  
#' given in the \code{input_sample_names} parameter of \code{create_html_report}).
#' 
#' @param output_rds_path String; full path to output folder for RDS object files created by this function.
#' Defaults to NOT save RDS; to save RDS, provide path to destination folder.
#' 
#' @return This function returns a plot that can be integrated in the html report and  
#' that can be saved as RDS object file.
#' 
#' @seealso \code{\link{create_html_report}}
#' @export

plot_codon_usage_bulk <- function(codon_usage_data, sample="", output_rds_path="") {
    
    d_bulk <- codon_usage_data$data_bulk
    d_bulk.m <- suppressMessages(melt(d_bulk))
    
    data_type <- codon_usage_data$data_type
    comp <- codon_usage_data$comp
    rl <- codon_usage_data$rl
    
    # translate data_type into a string suitable for the y axis label
    data_types <- c("Codon counts"="Codon_counts", 
                    "A-sites counts"="A_sites_percodon", 
                    "P-sites counts"="P_sites_percodon", 
                    "E-sites counts"="E_sites_percodon", 
                    "A-sites per codon"="A_sites_percodon_ratio",
                    "P-sites per codon"="P_sites_percodon_ratio",
                    "E-sites per codon"="E_sites_percodon_ratio")
    data_type <- names(data_types[data_types==data_type])
    
    # prepare subplots (grobs): g1 and g2
    g1 <- ggplot(d_bulk.m[match("ATG", d_bulk.m$codon), ], aes(x=codon, y=value)) +
        geom_bar(aes(fill=aa), stat="identity") + 
        facet_grid(.~aa, scales="free", space="free") +
        labs(x="", y=data_type) +
        scale_y_continuous(expand=c(0, 0)) +  # remove y axis offset
        #scale_colour_hue() +  # hue colors
        scale_fill_manual(values=c("gray")) + # same color for all amino acids
        #coord_fixed(ratio=0.2) +
        theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=16, family="Courier"), 
              axis.text.y=element_text(size=16),
              axis.title.x=element_blank(),
              axis.title.y=element_text(size=16),
              axis.ticks.x=element_blank(),
              panel.background=element_rect(),  # get panel background
              panel.border = element_blank(),  # remove borders from panel
              panel.grid.major = element_blank(),  # remove grid lines from panel
              panel.grid.minor = element_blank(),  # remove grid lines from panel
              strip.text=element_text(size=16),
              legend.position="none") 
    
    g2 <- ggplot(d_bulk.m[-match("ATG", d_bulk.m$codon), ], aes(x=codon, y=value)) +
        geom_bar(aes(fill=aa), stat="identity") + 
        facet_grid(.~aa, scales="free", space="free") +
        labs(x="", y=data_type) +
        scale_y_continuous(expand=c(0, 0)) +  # remove y axis offset
        scale_colour_hue() +  # hue colors
        #scale_fill_manual(values=rep("steelblue", 21)) + # same color for all amino acids
        #coord_fixed(ratio=0.2) +
        theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=16, family="Courier"), 
              axis.text.y=element_text(size=16),
              axis.title=element_blank(),
              axis.ticks.x=element_blank(),
              panel.background=element_rect(),  # get panel background
              panel.border = element_blank(),  # remove borders from panel
              panel.grid.major = element_blank(),  # remove grid lines from panel
              panel.grid.minor = element_blank(),  # remove grid lines from panel
              strip.text=element_text(size=16),
              legend.position="none") 
    
    # arrange plots
    gg <- arrangeGrob(g1, g2, ncol=2, widths = c(1, 9))
    
    cat("\n\n")
    
    # save plot as RDS object
    if(!output_rds_path==""){
        file_path = paste0(output_rds_path, sample, "_", 
                           comp, "_",
                           "8_codonusage_bulk_", 
                           gsub(" ", "_", data_type), "_", 
                           rl)
        saveRDS(list(gg, c(12,8)), file = file_path, 
                ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
    }
    
    # print to HTML file
    print(as_ggplot(gg))
}


######################################################################################################
# Methods for annotation and RiboseQC analysis
######################################################################################################


setMethods(f = "unlist","GRanges",function(x){return(x)})

setMethods(f = "unlist","GAlignments",function(x){return(x)})


#' Filter read lengths for P-sites position calculation
#'
#' This function selects a subset of readlenghts to be used in the P-sites calculation step
#' @keywords Ribo-seQC
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param summary_data output data from the \code{calc_cutoffs_from_profiles} function
#' @param choice Method used to select readlengths, defaults to "max_coverage"
#' @param nt_signals Profiles of 5'ends around start codons
#' @details Three different methods are available to choose readlengths: the "max_coverage" method
#' selects all read lenghts with more in-frame signal compared to out-of-frame signal, on all codons;
#' the "max_inframe" method starts with the most accurate read length and
#' progressively selects read lengths which add in-frame signals in codons not covered by previous read lengths;
#' the "all" method selects all available read lengths
#' @return a list object containing different compartments. Each sub-list contains \code{final_choice}, 
#' the set of chosen read lengths with cutoffs, and \code{data}, the complete stats for each selection method 
#' @seealso \code{\link{calc_cutoffs_from_profiles}}
#' @export

choose_readlengths<-function(summary_data,choice="max_coverage",nt_signals){
    
    if(!choice%in%c("max_coverage","max_inframe","all")){
        stop("choice must be one of 'max_coverage','max_inframe' or 'all'")
    }
    
    choices_res<-list()
    for(i in names(summary_data)){
        choices_res[[i]]<-c()
        data<-DataFrame(summary_data[[i]])
        sign<-nt_signals[[i]]
        data<-data[grep(rownames(data),pattern = "all",invert = T),]
        if(dim(data)[1]==0){next}
        data<-data[order(data$frame_preference,decreasing = T),]
        mat_fr<-matrix(0,nrow = nrow(sign[[1]]),ncol = 99)
        net_all<-c()
        net_all_cov<-c()
        net_all_cov_ok<-c()
        gains_cod_all<-c()
        gains_sign_all<-c()
        gains_cod<-c()
        gains_sign<-c()
        for(cnt in 1:length(rownames(data))){
            if(cnt>length(rownames(data))){break}
            rl<-rownames(data)[cnt]
            cutf<-data$cutoff[cnt]
            sig<-as.matrix(sign[[rl]])
            #sig<-do.call(args = sig[c("5-UTR_1","5-UTR_2","CDS_1","CDS_2","3-UTR_1")],what = cbind.data.frame)
            #sig<-sig[,-c(1,2)]
            cdss<-grep(colnames(sig),pattern = "CDS")
            sig_ctf<-cbind(matrix(0,ncol = cutf,nrow = nrow(sig)),sig[,1:(length(sig[1,])-cutf)])[,cdss]
            sig_ctf<-as.matrix(sig_ctf)
            sig_ctf_masked<-sig_ctf
            sig_ctf_masked[mat_fr>0]<-NA
            
            codon_cov<-t(apply(as.matrix(sig_ctf),1,function(y){
                fr1<-y[seq(1,99,by=3)]
                fr2<-y[seq(2,99,by=3)]
                fr3<-y[seq(3,99,by=3)]
                fr1<-fr1[!is.na(fr1)]
                fr2<-fr2[!is.na(fr2)]
                fr3<-fr3[!is.na(fr3)]
                codongain<-c(sum(fr1>0),sum(fr2>0),sum(fr3>0))
                codongain<-codongain/sum(codongain)
                codongain
                
            }))
            codon_cov<-codon_cov[!apply(codon_cov,1,FUN = function(y){sum(is.na(y))>0}),]
            
            
            codon_gain<-t(apply(sig_ctf_masked,1,function(y){
                fr1<-y[seq(1,99,by=3)]
                fr2<-y[seq(2,99,by=3)]
                fr3<-y[seq(3,99,by=3)]
                fr1<-fr1[!is.na(fr1)]
                fr2<-fr2[!is.na(fr2)]
                fr3<-fr3[!is.na(fr3)]
                codongain<-c(sum(fr1>0),sum(fr2>0),sum(fr3>0))
                codongain<-codongain/sum(codongain)
                codongain
                #signalgain<-c(sum(fr1),sum(fr2),sum(fr3))
                
            }))
            codon_gain<-codon_gain[!apply(codon_gain,1,FUN = function(y){sum(is.na(y))>0}),]
            
            
            av_codoncov<-matrix(colMeans(codon_cov),ncol=3)
            av_codongain<-matrix(colMeans(codon_gain),ncol=3)
            colnames(av_codongain)<-c("gain_in_frame","gain_off_1","gain_off_2")
            rownames(av_codongain)<-rl
            colnames(av_codoncov)<-c("in_frame","off_1","off_2")
            rownames(av_codoncov)<-rl
            
            
            
            gain_cod<-av_codongain[1]-(av_codongain[2]+av_codongain[3])
            gain_sign<-av_codoncov[1]-(av_codoncov[2]+av_codoncov[3])
            gains_cod_all<-c(gains_cod_all,gain_cod)
            gains_sign_all<-c(gains_sign_all,gain_sign)
            
            if(gain_cod<=0 & gain_sign<=0){next}
            if(gain_cod>0){
                mat_fr[sig_ctf>1]<-1
                net_all<-rbind(net_all,av_codongain)
                gains_cod<-c(gains_cod,gain_cod)
                
            }
            
            if(gain_sign>0){
                net_all_cov<-rbind(net_all_cov,av_codoncov)
                gains_sign<-c(gains_sign,gain_sign)
                
            }
            
            
        }
        
        if(length(gains_sign)>0){
            ln<-length(unique(gains_sign))
            if(ln<4){
                net_all_cov_ok<-net_all_cov
            }
            if(ln>3){
                km<-kmeans(gains_sign,3)
                net_all_cov_ok<-net_all_cov[!km$cluster==km$cluster[which.min(gains_sign)][1],]
            }
        }
        
        data$read_length<-rownames(data)
        data$gain_codons<-gains_sign_all
        data$gain_new_codons<-gains_cod_all
        
        
        data$all<-TRUE
        data$max_inframe<-FALSE
        if(length(net_all)>0){
            mtc<-match(rownames(net_all),rownames(data))
            data$max_inframe[mtc]<-TRUE
        }
        
        data$max_coverage<-FALSE
        if(length(net_all_cov)>0){
            mtc<-match(rownames(net_all_cov),rownames(data))
            data$max_coverage[mtc]<-TRUE
        }
        
        # data$best_coverage<-FALSE
        # if(length(net_all_cov)>0){
        #     mtc<-match(rownames(net_all_cov_ok),rownames(data))
        #     data$best_coverage[mtc]<-TRUE
        # }
        # 
        final_choice<-data[,c("read_length","cutoff")][data[,choice],]
        if(length(final_choice)==0){final_choice<-c()}
        choices<-list(final_choice=final_choice,data=data)
        choices_res[[i]]<-choices
    }
    choices_res
}


#' Offset spliced reads on plus strand
#'
#' This function calculates P-sites positions for spliced reads on the plus strand
#' @keywords Ribo-seQC
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param x a \code{GAlignments} object with a cigar string
#' @param cutoff number representing the offset value
#' @return a \code{GRanges} object with offset reads
#' @export

get_ps_fromspliceplus<-function(x,cutoff){
    rang<-cigarRangesAlongReferenceSpace(cigar(x), pos=start(x),ops="M")
    cs<-lapply(rang,function(x){cumsum(x@width)})
    rangok<-lapply(which(IntegerList(cs)>cutoff),"[[",1)
    rangok<-unlist(rangok)
    gr<-as(x,"GRanges")
    #do the splitaslist thingy
    ones<-rangok==1
    mores<-rangok>1
    psmores<-GRanges()
    psones<-GRanges()
    
    if(sum(ones)>0){
        psones<-shift(resize(gr[ones],width=1,fix="start"),shift=cutoff)
    }
    
    if(sum(mores)>0){
        rangmore<-rang[mores]
        rangok<-rangok[mores]
        cms<-cumsum(width(rangmore))
        shft<-c()
        for(i in 1:length(rangok)){shft<-c(shft,cutoff-cms[[i]][rangok[i]-1])}
        stt<-start(rangmore)
        stok<-c()
        for(i in 1:length(shft)){stok<-c(stok,stt[[i]][rangok[i]]+shft[i])}
        psmores<-GRanges(IRanges(start=stok,width = 1),seqnames = seqnames(x[mores]),strand=strand(x[mores]),seqlengths=seqlengths(x[mores]))
        
    }
    ps<-sort(c(psones,psmores))
    return(ps)
}

#' Offset spliced reads on minus strand
#'
#' This function calculates P-sites positions for spliced reads on the minus strand
#' @keywords Ribo-seQC
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param x a \code{GAlignments} object with a cigar string
#' @param cutoff number representing the offset value
#' @return a \code{GRanges} object with offset reads
#' @export

get_ps_fromsplicemin<-function(x,cutoff){
    rang<-cigarRangesAlongReferenceSpace(cigar(x), pos=start(x),ops="M")
    rang<-endoapply(rang,rev)
    cs<-lapply(rang,function(x){cumsum(x@width)})
    rangok<-lapply(which(IntegerList(cs)>cutoff),"[[",1)
    rangok<-unlist(rangok)
    gr<-as(x,"GRanges")
    #do the splitaslist thingy
    ones<-rangok==1
    mores<-rangok>1
    psmores<-GRanges()
    psones<-GRanges()
    
    
    if(sum(ones)>0){
        psones<-shift(resize(gr[ones],width=1,fix="start"),shift=-cutoff)
    }
    
    
    if(sum(mores)>0){
        rangmore<-rang[mores]
        rangok<-rangok[mores]
        cms<-cumsum(width(rangmore))
        shft<-c()
        for(i in 1:length(rangok)){shft<-c(shft,cutoff-cms[[i]][rangok[i]-1])}
        #start?
        stt<-end(rangmore)
        stok<-c()
        for(i in 1:length(shft)){stok<-c(stok,stt[[i]][rangok[i]]-shft[i])}
        psmores<-GRanges(IRanges(start=stok,width = 1),seqnames = seqnames(x[mores]),strand=strand(x[mores]),seqlengths=seqlengths(x[mores]))
        
    }
    ps<-sort(c(psones,psmores))
    return(ps)
    
    
    if(rangok!=1){
        
        ps<-shift(resize(GRanges(rang[rangok],seqnames=seqnames(x),strand=strand(x),seqlengths=seqlengths(x)),width=1,fix="start"),shift=-(cutoff-sum(rang[1:(rangok-1)]@width)))
    }
    return(ps)
}

#' Load genomic features and genome sequence
#'
#' This function loads the annotation created by the \code{prepare_annotation_files function}
#' @keywords Ribo-seQC
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param path Full path to the *Rannot R file in the annotation directory used in the \code{prepare_annotation_files function}
#' @return introduces a \code{GTF_annotation} object and a \code{genome_seq} object in the parent environment
#' @seealso \code{\link{prepare_annotation_files}}
#' @export

load_annotation<-function(path){
    GTF_annotation<-get(load(path))
    #genome_sequence<-get(library(GTF_annotation$genome_package,character.only = T))
    
    library(GTF_annotation$genome_package,character.only = T)
    genome_sequence<-get(GTF_annotation$genome_package)
    
    GTF_annotation<<-GTF_annotation
    genome_seq<<-genome_sequence
}

#' Calculate offsets from 5'end profiles
#'
#' This function calculates cutoffs and frame resolution for Ribo-seq reads,
#'  for each read length and compartment. 
#' @details Three methods are used and combined in the final choice: 
#'  the position of maximum coverage around start codon is calculated for each transcript, and 
#'  the most frequent one is stored in the "*_tab" objects. Such frequency values are also
#'  subjected to k-means clustering (centers=3) and the first value belonging to the highest cluster
#'  is selected, output as "*km_tab" objects. Analysis of aggregate plots, instead of frequencies, is performed
#'  again using kmeans (centers=3) using the same analysis above and stored in the "*km_meta" objects,
#'  and by simply calculating the maximum value in the profile, stored in the "*meta" objects.
#'  for each method, all reads ("absolute_") or only in-frame positions ("in_frame_") are considered. The final choice takes the most 
#'  frequent cutoff chosen in all methods applied to in-frame positions. 
#' @param reads_profile Profile of 5'ends around start and stop codon, as a DataFrame object with
#'   tx_ids as rows and positions as columns
#' @param length_max Maximum cutoff to use
#' @keywords Ribo-seQC
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @return a list with a \code{final_cutoff} object, the frame analysis containing the % of transcripts 
#' displaying the max frame and the average % of frame precision per transcript in \code{frames_res}, 
#' all the calculated cutoffs in \code{cutoffs}, data used for the frame analysis in \code{frames}, 
#' and profiles around start codons in  \code{profiles_start}.
#' @seealso \code{\link{RiboseQC_analysis}}
#' @export

calc_cutoffs_from_profiles<-function(reads_profile,length_max){
    profok<-as.matrix(reads_profile[,grep("CDS",x = colnames(reads_profile))])
    #profok<-do.call(args = reads_profile[c("CDS_1","CDS_2")],what = cbind.data.frame)
    sums<-rowSums(profok)
    profok<-profok[sums>0,,drop=F]
    #what if no reads?
    if((dim(profok)[1])<5){return(c())}
    
    sums<-rowSums(profok)
    frames<-t(apply(profok,1,function(y){
        sumy<-sum(y)
        fr0<-sum(y[seq(1,length(y),by=3)])/sumy
        fr1<-sum(y[seq(2,length(y),by=3)])/sumy
        fr2<-sum(y[seq(3,length(y),by=3)])/sumy
        c(fr0,fr1,fr2)
    }))
    
    tbfra<-table(apply(frames,1,which.max))
    tbfra<-(tbfra/sum(tbfra))*100
    fra<-as.numeric(names(which.max(tbfra)))-1
    frameans<-colMeans(frames)
    frameans<-(frameans/sum(frameans))*100
    names(frameans)<-1:3
    fra2<-as.numeric(names(which.max(frameans)))-1
    
    infra<-seq(fra,-24,by=-3)
    #dists<-c(-24:123)
    #profok<-do.call(args = reads_profile[c("5-UTR_2","CDS_1","CDS_2","3-UTR_1")],what = cbind.data.frame)
    
    #just offset
    profok<-as.matrix(reads_profile)
    Atg<-which(colnames(profok)=="CDS_1")
    colnames(profok)<-as.character(seq(from=-(Atg-1),length.out = length(profok[1,]),by=1))
    
    profok<-profok[,colnames(profok)%in%as.character(-24:0)]
    profok<-profok[,abs(as.numeric(colnames(profok)))<=length_max]
    sums<-rowSums(profok)
    profok<-profok[sums>0,,drop=F]
    #what if no reads?
    if((dim(profok)[1])<5){return(c())}
    
    
    metaprof<-colSums(profok)
    names(metaprof)<-as.numeric(names(metaprof))
    metaprof_fra<-metaprof[names(metaprof)%in%infra]
    km_meta<-NA
    km_meta_fra<-NA
    if(length(unique(metaprof))>3){
        kms<-kmeans(metaprof,centers=3)$cluster
        km_meta<-abs(as.numeric(names(metaprof[which(kms==kms[which.max(metaprof)])[1]])))
    }
    if(length(unique(metaprof_fra))>3){
        kms<-kmeans(metaprof_fra,centers=3)$cluster
        km_meta_fra<-abs(as.numeric(names(metaprof_fra[which(kms==kms[which.max(metaprof_fra)])[1]])))
        
    }
    tbcut<-table(apply(profok,1,which.max))
    names(tbcut)<-abs(as.numeric(colnames(profok)[as.numeric(names(tbcut))]))
    profcut<-c(rep(0,25))
    names(profcut)<-0:24
    profcut<-profcut[as.numeric(names(profcut))<=length_max]
    profcut[names(tbcut)]<-as.numeric(tbcut)
    cutoff<-abs(as.numeric(names(tbcut)[which.max(tbcut)]))
    cutoff_km<-NA
    cutoff_km_fra<-NA
    if(length(cutoff)==0){cutoff=NA}
    tbcut_fra<-tbcut[as.numeric(names(tbcut))%in%(-infra)]
    if(length(tbcut_fra)==0){tbcut_fra=NA}
    cutoff_fra<-abs(as.numeric(names(tbcut_fra)[which.max(tbcut_fra)]))
    if(length(cutoff_fra)==0){cutoff_fra=NA}
    if(length(cutoff)>0 & length(unique(tbcut))>4){
        kms<-kmeans(tbcut,centers=3)$cluster
        cutoff_km<-abs(as.numeric(names(tbcut[which(kms==kms[which.max(tbcut)])[1]])))
        if(length(unique(tbcut_fra))>3){
            kms<-kmeans(tbcut_fra,centers=3)$cluster
            cutoff_km_fra<-abs(as.numeric(names(tbcut_fra[which(kms==kms[which.max(tbcut_fra)])[1]])))
        }
    }
    cutoff_meta<-abs(as.numeric(names(metaprof)[which.max(metaprof)]))
    cutoff_meta_fra<-abs(as.numeric(names(metaprof_fra)[which.max(metaprof_fra)]))
    if(length(cutoff_meta)==0){cutoff_meta=NA}
    if(length(cutoff_meta_fra)==0){cutoff_meta_fra=NA}
    
    cutoffs<-c(cutoff,cutoff_fra,cutoff_km,cutoff_km_fra,km_meta,km_meta_fra,cutoff_meta,cutoff_meta_fra)
    names(cutoffs)<-c("absolute_tab","in_frame_tab","absolute_km_tab","in_frame_km_tab","absolute_km_meta","in_frame_km_meta","absolute_meta","in_frame_meta")
    cutt<-cutoffs[!is.na(cutoffs)]
    final_cutoff<-as.numeric(names(sort(table(cutt[grep(names(cutt),pattern = "frame")]),decreasing = T)[1]))
    frames_res<-c(sort(tbfra,decreasing = T)[1],as.numeric(frameans[which.max(frameans)]))
    names(frames_res)<-c("max_frame_in_pctORF","mean_pct_max_frame")
    
    results<-list(final_cutoff=final_cutoff,frames_res=frames_res,cutoffs=cutoffs,frames=frames,profiles_start=profok)
    return(results)
    
}



#' Prepare comprehensive sets of annotated genomic features
#'
#' This function processes a gtf file and a twobit file (created using faToTwoBit from ucsc tools: http://hgdownload.soe.ucsc.edu/admin/exe/ ) to create a comprehensive set of genomic regions of interest in genomic and transcriptomic space (e.g. introns, UTRs, start/stop codons).
#'    In addition, by linking genome sequence and annotation, it extracts additional info, such as gene and transcript biotypes, genetic codes for different organelles, or chromosomes and transcripts lengths.
#' @keywords Ribo-seQC
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param annotation_directory The target directory which will contain the output files
#' @param twobit_file Full path to the genome file in twobit format
#' @param gtf_file Full path to the annotation file in GTF format
#' @param scientific_name A name to give to the organism studied; must be two words separated by a ".", defaults to Homo.sapiens
#' @param annotation_name A name to give to annotation used; defaults to genc25
#' @param export_bed_tables_TxDb Export coordinates and info about different genomic regions in the annotation_directory? It defaults to \code{TRUE}
#' @param forge_BSgenome Forge and install a \code{BSgenome} package? It defaults to \code{TRUE}
#' @param create_TxDb Create a \code{TxDb} object and a *Rannot object? It defaults to \code{TRUE}
#' @details This function uses the \code{makeTxDbFromGFF} function to  create a TxDb object and extract
#' genomic regions and other info to a *Rannot R file; the \code{mapToTranscripts} and \code{mapFromTranscripts} functions are used to 
#' map features to genomic or transcript-level coordinates. GTF file mist contain "exon" and "CDS" lines,
#' where each line contains "transcript_id" and "gene_id" values. Additional values such as "gene_biotype" or "gene_name" are also extracted.
#' Regarding sequences, the twobit file, together with input scientific and annotation names, is used to forge and install a 
#' BSgenome package using the \code{forgeBSgenomeDataPkg} function.\cr\cr
#' The resulting GTF_annotation object (obtained after runnning \code{load_annotation}) contains:\cr\cr
#' \code{txs}: annotated transcript boundaries.\cr
#' \code{txs_gene}: GRangesList including transcript grouped by gene.\cr
#' \code{seqinfo}: indicating chromosomes and chromosome lengths.\cr
#' \code{start_stop_codons}: the set of annotated start and stop codon, with respective transcript and gene_ids.
#' reprentative_mostcommon,reprentative_boundaries and reprentative_5len represent the most common start/stop codon,
#' the most upstream/downstream start/stop codons and the start/stop codons residing on transcripts with the longest 5'UTRs\cr
#' \code{cds_txs}: GRangesList including CDS grouped by transcript.\cr
#' \code{introns_txs}: GRangesList including introns grouped by transcript.\cr
#' \code{cds_genes}: GRangesList including CDS grouped by gene.\cr
#' \code{exons_txs}: GRangesList including exons grouped by transcript.\cr
#' \code{exons_bins}: the list of exonic bins with associated transcripts and genes.\cr
#' \code{junctions}: the list of annotated splice junctions, with associated transcripts and genes.\cr
#' \code{genes}: annotated genes coordinates.\cr
#' \code{threeutrs}: collapsed set of 3'UTR regions, with correspinding gene_ids. This set does not overlap CDS region.\cr
#' \code{fiveutrs}: collapsed set of 5'UTR regions, with correspinding gene_ids. This set does not overlap CDS region.\cr
#' \code{ncIsof}: collapsed set of exonic regions of protein_coding genes, with correspinding gene_ids. This set does not overlap CDS region.\cr
#' \code{ncRNAs}: collapsed set of exonic regions of non_coding genes, with correspinding gene_ids. This set does not overlap CDS region.\cr
#' \code{introns}: collapsed set of intronic regions, with correspinding gene_ids. This set does not overlap exonic region.\cr
#' \code{intergenicRegions}: set of intergenic regions, defined as regions with no annotated genes on either strand.\cr
#' \code{trann}: DataFrame object including (when available) the mapping between gene_id, gene_name, gene_biotypes, transcript_id and transcript_biotypes.\cr
#' \code{cds_txs_coords}: transcript-level coordinates of ORF boundaries, for each annotated coding transcript. Additional columns are the same as as for the \code{start_stop_codons} object.\cr
#' \code{genetic_codes}: an object containing the list of genetic code ids used for each chromosome/organelle. see GENETIC_CODE_TABLE for more info.\cr
#' \code{genome_package}: the name of the forged BSgenome package. Loaded with \code{load_annotation} function.\cr
#' \code{stop_in_gtf}: stop codon, as defined in the annotation.\cr
#' @return a TxDb file and a *Rannot files are created in the specified \code{annotation_directory}. 
#' In addition, a BSgenome object is forged, installed, and linked to the *Rannot object
#' @seealso \code{\link{load_annotation}}, \code{\link{forgeBSgenomeDataPkg}}, \code{\link{makeTxDbFromGFF}}.
#' @export

prepare_annotation_files<-function(annotation_directory,twobit_file,gtf_file,scientific_name="Homo.sapiens",annotation_name="genc25",export_bed_tables_TxDb=T,forge_BSgenome=T,create_TxDb=T,additional_circ_seqs=NULL){
    
   if(!is.null(additional_circ_seqs)) stopifnot(is.character(additional_circ_seqs))
    DEFAULT_CIRC_SEQS <- unique(c("chrM","MT","MtDNA","mit","Mito","mitochondrion",
                                  "dmel_mitochondrion_genome","Pltd","ChrC","Pt","chloroplast",
                                  "Chloro","2micron","2-micron","2uM",
                                  "Mt", "NC_001879.2", "NC_006581.1","ChrM",additional_circ_seqs))
    #adjust variable names (some chars not permitted)
    annotation_name<-gsub(annotation_name,pattern = "_",replacement = "")
    annotation_name<-gsub(annotation_name,pattern = "-",replacement = "")
    if(!dir.exists(annotation_directory)){dir.create(path = annotation_directory,recursive = T)}
    annotation_directory<-normalizePath(annotation_directory)
    twobit_file<-normalizePath(twobit_file)
    gtf_file<-normalizePath(gtf_file)
    
    for (f in c(twobit_file,gtf_file)){
        if(file.access(f, 0)==-1) {
            stop("
                 The following files don't exist:\n",
                 f, "\n")
        }
    }
    
    
    scientific_name_spl<-strsplit(scientific_name,"[.]")[[1]]
    ok<-length(scientific_name_spl)==2
    if(!ok){stop("\"scientific_name\" must be two words separated by a \".\", like \"Homo.sapiens\"")}
    
    #get circular sequences
    
    seqinfotwob<-seqinfo(TwoBitFile(twobit_file))
    circss<-seqnames(seqinfotwob)[which(seqnames(seqinfotwob)%in%DEFAULT_CIRC_SEQS)]
    seqinfotwob@is_circular[which(seqnames(seqinfotwob)%in%DEFAULT_CIRC_SEQS)]<-TRUE
    pkgnm<-paste("BSgenome",scientific_name,annotation_name,sep=".")
    
    circseed<-circss
    if(length(circseed)==0){circseed<-NULL}
    
    #Forge a BSGenome package
    
    if(forge_BSgenome){
        cat(paste("Creating the BSgenome package ... ",date(),"\n",sep = ""))
        seed_text<-paste("Package: BSgenome.",scientific_name,".",annotation_name,"\n",
                         "Title: Full genome sequences for ",scientific_name,", ",annotation_name,"\n",
                         "Description: Full genome sequences for ",scientific_name,", ",annotation_name,"\n",
                         "Version: 1.0","\n",
                         "organism: ",scientific_name,"\n",
                         "common_name: ",scientific_name,"\n",
                         "provider: NA","\n",
                         "provider_version: ",annotation_name,"\n",
                         "release_date: NA","\n",
                         "release_name: NA","\n",
                         "source_url: NA","\n",
                         "organism_biocview: ", scientific_name,"\n",
                         "BSgenomeObjname: ",scientific_name,"\n",
                         "seqs_srcdir: ",dirname(twobit_file),"\n",
                         "seqfile_name: ",basename(twobit_file),sep="")
        
        
        seed_dest<-paste(annotation_directory,"/",basename(twobit_file),"_",scientific_name,"_seed",sep = "")
        
        if(length(circseed)==1){
            seed_text<-paste(seed_text,"\n",
                             "circ_seqs: \"",circseed,"\"",sep="")
            writeLines(text = seed_text,con = seed_dest)
        }
        
        if(length(circseed)>1){
            circseed<-paste('c("',paste(circseed,collapse=","),'")',sep="")
            circseed<-gsub(circseed,pattern = ",",replacement='","')
            
            cat(seed_text,"\n","circ_seqs: ",circseed,"\n",sep="",file = seed_dest)
        }
        
        
        
        unlink(paste(annotation_directory,pkgnm,sep="/"),recursive=T)
        
        forgeBSgenomeDataPkg(x=seed_dest,destdir=annotation_directory,seqs_srcdir=dirname(twobit_file))
        cat(paste("Creating the BSgenome package --- Done! ",date(),"\n",sep = ""))
        
        cat(paste("Installing the BSgenome package ... ",date(),"\n",sep = ""))
        
        install(paste(annotation_directory,pkgnm,sep="/"))
        cat(paste("Installing the BSgenome package --- Done! ",date(),"\n",sep = ""))
        
    }
    
    #Create the TxDb from GTF and BSGenome info
    
    if(create_TxDb){
        cat(paste("Creating the TxDb object ... ",date(),"\n",sep = ""))
        
        annotation<-makeTxDbFromGFF(file=gtf_file,format="gtf",chrominfo = seqinfotwob)
        
        saveDb(annotation, file=paste(annotation_directory,"/",basename(gtf_file),"_TxDb",sep=""))
        cat(paste("Creating the TxDb object --- Done! ",date(),"\n",sep = ""))
        cat(paste("Extracting genomic regions ... ",date(),"\n",sep = ""))
        
        genes<-genes(annotation)
        exons_ge<-exonsBy(annotation,by="gene")
        exons_ge<-reduce(exons_ge)
        
        cds_gen<-cdsBy(annotation,"gene")
        cds_ge<-reduce(cds_gen)
        
        
        #define regions not overlapping CDS ( or exons when defining introns)
        
        threeutrs<-reduce(GenomicRanges::setdiff(unlist(threeUTRsByTranscript(annotation)),unlist(cds_ge),ignore.strand=FALSE))
        
        fiveutrs<-reduce(GenomicRanges::setdiff(unlist(fiveUTRsByTranscript(annotation)),unlist(cds_ge),ignore.strand=FALSE))
        
        introns<-reduce(GenomicRanges::setdiff(unlist(intronsByTranscript(annotation)),unlist(exons_ge),ignore.strand=FALSE))
        
        nc_exons<-reduce(GenomicRanges::setdiff(unlist(exons_ge),reduce(c(unlist(cds_ge),fiveutrs,threeutrs)),ignore.strand=FALSE))
        
        #assign gene ids (mutiple when overlapping multiple genes)
        ov<-findOverlaps(threeutrs,genes)
        ov<-split(subjectHits(ov),queryHits(ov))
        threeutrs$gene_id<-CharacterList(lapply(ov,FUN = function(x){names(genes)[x]}))
        ov<-findOverlaps(fiveutrs,genes)
        ov<-split(subjectHits(ov),queryHits(ov))
        fiveutrs$gene_id<-CharacterList(lapply(ov,FUN = function(x){names(genes)[x]}))
        ov<-findOverlaps(introns,genes)
        ov<-split(subjectHits(ov),queryHits(ov))
        introns$gene_id<-CharacterList(lapply(ov,FUN = function(x){names(genes)[x]}))
        ov<-findOverlaps(nc_exons,genes)
        ov<-split(subjectHits(ov),queryHits(ov))
        nc_exons$gene_id<-CharacterList(lapply(ov,FUN = function(x){names(genes)[x]}))
        
        intergenicRegions<-genes
        strand(intergenicRegions)<-"*"
        intergenicRegions <- gaps(reduce(intergenicRegions))
        intergenicRegions<-intergenicRegions[strand(intergenicRegions)=="*"]
        
        cds_tx<-cdsBy(annotation,"tx",use.names=T)
        txs_gene<-transcriptsBy(annotation,by="gene")
        genes_red<-reduce(sort(genes(annotation)))
        
        exons_tx<-exonsBy(annotation,"tx",use.names=T)
        
        transcripts_db<-transcripts(annotation)
        intron_names_tx<-intronsByTranscript(annotation,use.names=T)
        
        
        #define exonic bins, including regions overlapping multiple genes
        nsns<-disjointExons(annotation,aggregateGenes=T)
        
        
        
        #define tx_coordinates of ORF boundaries
        
        exsss_cds<-exons_tx[names(cds_tx)]
        chunks<-seq(1,length(cds_tx),by = 20000)
        if(chunks[length(chunks)]<length(cds_tx)){chunks<-c(chunks,length(cds_tx))}
        mapp<-GRangesList()
        for(i in 1:(length(chunks)-1)){
            if(i!=(length(chunks)-1)){
                mapp<-suppressWarnings(c(mapp,pmapToTranscripts(cds_tx[chunks[i]:(chunks[i+1]-1)],transcripts = exsss_cds[chunks[i]:(chunks[i+1]-1)])))
            }
            if(i==(length(chunks)-1)){
                mapp<-suppressWarnings(c(mapp,pmapToTranscripts(cds_tx[chunks[i]:(chunks[i+1])],transcripts = exsss_cds[chunks[i]:(chunks[i+1])])))
            }
        }
        cds_txscoords<-unlist(mapp)
        
        
        #extract biotypes and ids
        
        cat(paste("Extracting ids and biotypes ... ",date(),"\n",sep = ""))
        
        trann<-unique(mcols(import.gff2(gtf_file,colnames=c("gene_id","gene_biotype","gene_type","gene_name","gene_symbol","transcript_id","transcript_biotype","transcript_type"))))
        trann<-trann[!is.na(trann$transcript_id),]
        trann<-data.frame(unique(trann),stringsAsFactors=F)
        
        if(sum(!is.na(trann$transcript_biotype))==0 & sum(!is.na(trann$transcript_type))==0 ){
            trann$transcript_biotype<-"no_type"
        }
        if(sum(!is.na(trann$transcript_biotype))==0){trann$transcript_biotype<-NULL}
        if(sum(!is.na(trann$transcript_type))==0){trann$transcript_type<-NULL}
        
        
        if(sum(!is.na(trann$gene_biotype))==0 & sum(!is.na(trann$gene_type))==0 ){
            
            trann$gene_type<-"no_type"
            
        }
        if(sum(!is.na(trann$gene_name))==0 & sum(!is.na(trann$gene_symbol))==0 ){
            
            trann$gene_name<-"no_name"
            
        }
        if(sum(!is.na(trann$gene_biotype))==0){trann$gene_biotype<-NULL}
        if(sum(!is.na(trann$gene_type))==0){trann$gene_type<-NULL}
        if(sum(!is.na(trann$gene_name))==0){trann$gene_name<-NULL}
        if(sum(!is.na(trann$gene_symbol))==0){trann$gene_symbol<-NULL}
        colnames(trann)<-c("gene_id","gene_biotype","gene_name","transcript_id","transcript_biotype")
        
        trann<-DataFrame(trann)
        
        
        
        #introns and transcript_ids/gene_ids
        unq_intr<-sort(unique(unlist(intron_names_tx)))
        names(unq_intr)<-NULL
        all_intr<-unlist(intron_names_tx)
        
        ov<-findOverlaps(unq_intr,all_intr,type="equal")
        ov<-split(subjectHits(ov),queryHits(ov))
        a_nam<-CharacterList(lapply(ov,FUN = function(x){unique(names(all_intr)[x])}))
        
        unq_intr$type="J"
        unq_intr$tx_name<-a_nam
        
        
        mat_genes<-match(unq_intr$tx_name,trann$transcript_id)
        g<-unlist(apply(cbind(1:length(mat_genes),Y = elementNROWS(mat_genes)),FUN =function(x) rep(x[1],x[2]),MARGIN = 1))
        g2<-split(trann[unlist(mat_genes),"gene_id"],g)
        unq_intr$gene_id<-CharacterList(lapply(g2,unique))
        
        
        #filter ncRNA and ncIsof regions
        ncrnas<-nc_exons[!nc_exons%over%genes[trann$gene_id[trann$gene_biotype=="protein_coding"]]]
        ncisof<-nc_exons[nc_exons%over%genes[trann$gene_id[trann$gene_biotype=="protein_coding"]]]
        
        
        # define genetic codes to use
        # IMPORTANT : modify if needed (e.g. different organelles or species) check ids of GENETIC_CODE_TABLE for more info
        
        ifs<-seqinfo(annotation)
        translations<-as.data.frame(ifs)
        translations$genetic_code<-"1"
        
        #insert new codes for chromosome name
        
        #Mammalian mito
        translations$genetic_code[rownames(translations)%in%c("chrM","MT","MtDNA","mit","mitochondrion")]<-"2"
        
        #Yeast mito
        translations$genetic_code[rownames(translations)%in%c("Mito")]<-"3"
        
        #Drosophila mito
        translations$genetic_code[rownames(translations)%in%c("dmel_mitochondrion_genome")]<-"5"
        
        circs<-ifs@seqnames[which(ifs@is_circular)]
        
        
        #define start and stop codons (genome space)
        
        suppressPackageStartupMessages(library(pkgnm,character.only=TRUE))
        genome<-get(pkgnm)
        tocheck<-as.character(runValue(seqnames(cds_tx)))
        tocheck<-cds_tx[!tocheck%in%circs]
        seqcds<-extractTranscriptSeqs(genome,transcripts = tocheck)
        cd<-unique(translations$genetic_code[!rownames(translations)%in%circs])
        trsl<-suppressWarnings(translate(seqcds,genetic.code = getGeneticCode(cd),if.fuzzy.codon = "solve"))
        trslend<-as.character(narrow(trsl,end = width(trsl),width = 1))
        stop_inannot<-NA
        if(names(sort(table(trslend),decreasing = T)[1])=="*"){stop_inannot<-"*"}
        
        cds_txscoords$gene_id<-trann$gene_id[match(as.vector(seqnames(cds_txscoords)),trann$transcript_id)]
        cds_cc<-cds_txscoords
        strand(cds_cc)<-"*"
        sta_cc<-resize(cds_cc,width = 1,"start")
        sta_cc<-unlist(pmapFromTranscripts(sta_cc,exons_tx[seqnames(sta_cc)],ignore.strand=F))
        sta_cc$gene_id<-trann$gene_id[match(names(sta_cc),trann$transcript_id)]
        sta_cc<-sta_cc[sta_cc$hit]
        strand(sta_cc)<-structure(as.vector(strand(transcripts_db)),names=transcripts_db$tx_name)[names(sta_cc)]
        sta_cc$type<-"start_codon"
        mcols(sta_cc)<-mcols(sta_cc)[,c("exon_rank","type","gene_id")]
        
        sto_cc<-resize(cds_cc,width = 1,"end")
        #stop codon is the 1st nt, e.g. U of the UAA
        #To-do: update with regards to different organelles, and different annotations
        sto_cc<-shift(sto_cc,-2)
        if(is.na(stop_inannot)){sto_cc<-shift(sto_cc,3)}
        
        sto_cc<-unlist(pmapFromTranscripts(sto_cc,exons_tx[seqnames(sto_cc)],ignore.strand=F))
        sto_cc<-sto_cc[sto_cc$hit]
        sto_cc$gene_id<-trann$gene_id[match(names(sto_cc),trann$transcript_id)]
        strand(sto_cc)<-structure(as.vector(strand(transcripts_db)),names=transcripts_db$tx_name)[names(sto_cc)]
        sto_cc$type<-"stop_codon"
        mcols(sto_cc)<-mcols(sto_cc)[,c("exon_rank","type","gene_id")]
        
        
        #define most common, most upstream/downstream
        
        cat(paste("Defining most common start/stop codons ... ",date(),"\n",sep = ""))
        
        start_stop_cc<-sort(c(sta_cc,sto_cc))
        start_stop_cc$transcript_id<-names(start_stop_cc)
        start_stop_cc$most_up_downstream<-FALSE
        start_stop_cc$most_frequent<-FALSE
        
        df<-cbind.DataFrame(start(start_stop_cc),start_stop_cc$type,start_stop_cc$gene_id)
        colnames(df)<-c("start_pos","type","gene_id")
        upst<-by(df$start_pos,INDICES = df$gene_id,function(x){x==min(x) | x==max(x)})
        start_stop_cc$most_up_downstream<-unlist(upst[unique(df$gene_id)])
        
        mostfr<-by(df[,c("start_pos","type")],INDICES = df$gene_id,function(x){
            mfreq<-table(x)
            x$start_pos%in%as.numeric(names(which(mfreq[,1]==max(mfreq[,1])))) | x$start_pos%in%as.numeric(names(which(mfreq[,2]==max(mfreq[,2]))))
        })
        
        start_stop_cc$most_frequent<-unlist(mostfr[unique(df$gene_id)])
        
        names(start_stop_cc)<-NULL
        
        
        
        #define transcripts as containing frequent start/stop codons or most upstream ones, in relation with 5'UTR length
        
        mostupstr_tx<-sum(LogicalList(split(start_stop_cc$most_up_downstream,start_stop_cc$transcript_id)))[as.character(seqnames(cds_txscoords))]
        cds_txscoords$upstr_stasto<-mostupstr_tx
        mostfreq_tx<-sum(LogicalList(split(start_stop_cc$most_frequent,start_stop_cc$transcript_id)))[as.character(seqnames(cds_txscoords))]
        cds_txscoords$mostfreq_stasto<-mostfreq_tx
        cds_txscoords$lentx<-sum(width(exons_tx[as.character(seqnames(cds_txscoords))]))
        df<-cbind.DataFrame(as.character(seqnames(cds_txscoords)),width(cds_txscoords),start(cds_txscoords),cds_txscoords$mostfreq_stasto,cds_txscoords$gene_id)
        colnames(df)<-c("txid","cdslen","utr5len","var","gene_id")
        repres_freq<-by(df[,c("txid","cdslen","utr5len","var")],df$gene_id,function(x){
            x<-x[order(x$var,x$utr5len,x$cdslen,decreasing = T),]
            x<-x[x$var==max(x$var),]
            ok<-x$txid[which(x$cdslen==max(x$cdslen) & x$utr5len==max(x$utr5len) & x$var==max(x$var))][1]
            if(length(ok)==0 | is.na(ok[1])){ok<-x$txid[1]}
            ok
        })
        
        df<-cbind.DataFrame(as.character(seqnames(cds_txscoords)),width(cds_txscoords),start(cds_txscoords),cds_txscoords$upstr_stasto,cds_txscoords$gene_id)
        colnames(df)<-c("txid","cdslen","utr5len","var","gene_id")
        repres_upstr<-by(df[,c("txid","cdslen","utr5len","var")],df$gene_id,function(x){
            x<-x[order(x$var,x$utr5len,x$utr5len,decreasing = T),]
            x<-x[x$var==max(x$var),]
            ok<-x$txid[which(x$cdslen==max(x$cdslen) & x$utr5len==max(x$utr5len) & x$var==max(x$var))][1]
            if(length(ok)==0 | is.na(ok[1])){ok<-x$txid[1]}
            ok
        })
        df<-cbind.DataFrame(as.character(seqnames(cds_txscoords)),width(cds_txscoords),start(cds_txscoords),cds_txscoords$upstr_stasto,cds_txscoords$gene_id)
        colnames(df)<-c("txid","cdslen","utr5len","var","gene_id")
        repres_len5<-by(df[,c("txid","cdslen","utr5len","var")],df$gene_id,function(x){
            x<-x[order(x$utr5len,x$var,x$cdslen,decreasing = T),]
            ok<-x$txid[which(x$utr5len==max(x$utr5len) & x$var==max(x$var))][1]
            if(length(ok)==0 | is.na(ok[1])){ok<-x$txid[1]}
            ok
        })
        
        cds_txscoords$reprentative_mostcommon<-as.character(seqnames(cds_txscoords))%in%unlist(repres_freq)
        cds_txscoords$reprentative_boundaries<-as.character(seqnames(cds_txscoords))%in%unlist(repres_upstr)
        cds_txscoords$reprentative_5len<-as.character(seqnames(cds_txscoords))%in%unlist(repres_len5)
        unq_stst<-start_stop_cc
        mcols(unq_stst)<-NULL
        unq_stst<-sort(unique(unq_stst))
        ov<-findOverlaps(unq_stst,start_stop_cc,type="equal")
        ov<-split(subjectHits(ov),queryHits(ov))
        unq_stst$type<-CharacterList(lapply(ov,FUN = function(x){unique(start_stop_cc$type[x])}))
        unq_stst$transcript_id<-CharacterList(lapply(ov,FUN = function(x){start_stop_cc$transcript_id[x]}))
        unq_stst$gene_id<-CharacterList(lapply(ov,FUN = function(x){unique(start_stop_cc$gene_id[x])}))
        
        unq_stst$reprentative_mostcommon<-sum(!is.na(match(unq_stst$transcript_id,unlist(as(repres_freq,"CharacterList")))))>0
        unq_stst$reprentative_boundaries<-sum(!is.na(match(unq_stst$transcript_id,unlist(as(repres_upstr,"CharacterList")))))>0
        unq_stst$reprentative_5len<-sum(!is.na(match(unq_stst$transcript_id,unlist(as(repres_len5,"CharacterList")))))>0
        
        
        #put in a list
        GTF_annotation<-list(transcripts_db,txs_gene,ifs,unq_stst,cds_tx,intron_names_tx,cds_gen,exons_tx,nsns,unq_intr,genes,threeutrs,fiveutrs,ncisof,ncrnas,introns,intergenicRegions,trann,cds_txscoords,translations,pkgnm,stop_inannot)
        names(GTF_annotation)<-c("txs","txs_gene","seqinfo","start_stop_codons","cds_txs","introns_txs","cds_genes","exons_txs","exons_bins","junctions","genes","threeutrs","fiveutrs","ncIsof","ncRNAs","introns","intergenicRegions","trann","cds_txs_coords","genetic_codes","genome_package","stop_in_gtf")
        
        #Save as a RData object
        save(GTF_annotation,file=paste(annotation_directory,"/",basename(gtf_file),"_Rannot",sep=""))
        cat(paste("Rannot object created!   ",date(),"\n",sep = ""))
        
        
        #create tables and bed files (with colnames, so with header)
        if(export_bed_tables_TxDb==T){
            cat(paste("Exporting annotation tables ... ",date(),"\n",sep = ""))
            for(bed_file in c("fiveutrs","threeutrs","ncIsof","ncRNAs","introns","cds_txs_coords")){
                bf<-GTF_annotation[[bed_file]]
                bf_t<-data.frame(chromosome=seqnames(bf),start=start(bf),end=end(bf),name=".",score=width(bf),strand=strand(bf))
                meccole<-mcols(bf)
                for(mecc in names(meccole)){
                    if(is(meccole[,mecc],"CharacterList") | is(meccole[,mecc],"NumericList") | is(meccole[,mecc],"IntegerList")){
                        meccole[,mecc]<-paste(meccole[,mecc],collapse=";")
                    }
                }
                bf_t<-cbind.data.frame(bf_t,meccole)
                write.table(bf_t,file = paste(annotation_directory,"/",bed_file,"_similbed.bed",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
                
            }
            
            write.table(GTF_annotation$trann,file = paste(annotation_directory,"/table_gene_tx_IDs",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
            seqi<-as.data.frame(GTF_annotation$seqinfo)
            seqi$chromosome<-rownames(seqi)
            write.table(seqi,file = paste(annotation_directory,"/seqinfo",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
            
            gen_cod<-as.data.frame(GTF_annotation$genetic_codes)
            gen_cod$chromosome<-rownames(gen_cod)
            write.table(gen_cod,file = paste(annotation_directory,"/genetic_codes",sep=""),sep="\t",quote = FALSE,row.names = FALSE)
            cat(paste("Exporting annotation tables --- Done! ",date(),"\n",sep = ""))
            
        }
        
    }
    
}


#' Perform a Ribo-seQC analysis
#'
#' This function loads annotation created by the prepare_annotation_files function, and analyzes a BAM file.
#' @keywords Ribo-seQC
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param annotation_file Full path to the annotation file (*Rannot). Or, a vector with paths to one annotation file per bam file.
#' @param bam_files character vector containing the full path to the bam files
#' @param chunk_size the number of alignments to read at each iteration, defaults to 5000000, increase when more RAM is available. Must be between 10000 and 100000000
#' @param read_subset Select readlengths up to 99 percent of the reads, defaults to \code{TRUE}. Must be of length 1 or same length as bam_files.
#' @param readlength_choice_method Method used to subset relevant read lengths (see \code{choose_readlengths} function); defaults to "max_coverage". Must be of length 1 or same length as bam_files.
#' @param rescue_all_rls Set cutoff of 12 for read lengths ignored because of insufficient coverage. Defaults to \code{FALSE}. Must be of length 1 or same length as bam_files.
#' @param write_tmp_files Should output all the results (in *results_RiboseQC_all)? Defaults to \code{TRUE}. Must be of length 1 or same length as bam_files.
#' @param dest_names character vector containing the prefixes to use for the result output files. Defaults to same as \code{bam_files}
#' @param fast_mode Use only top 500 genes to build profiles? Defaults to \code{TRUE}. Must be of length 1 or same length as bam_files.
#' @param create_report Create an html report showing the RiboseQC analysis results. Defaults to \code{TRUE}
#' @param sample_names character vector containing the names for each sample analyzed (for the html report). Defaults to "sample1", "sample2" ...
#' @param report_file desired filename for for the html report file. Defaults to the first entry of \code{bam_files} followed by ".html"
#' @param extended_report creates a large html report including codon occupancy for each read length. Defaults to \code{FALSE}
#' @param pdf_plots creates a pdf file for each produced plot. Defaults to \code{TRUE}
#' @return the function saves a "results_RiboseQC_all" R file appended to the bam_files path including the complete list of outputs described here.
#' In addition, bigwig files for coverage value and P_sites position is appended to the bam_files path, including also a summary of P_sites selection statistics, 
#' a smaller "results_RiboseQC" R file used for creating a dynamic html report, and a "for_SaTAnn" R object that can be used in the SaTAnn pipeline.
#' @details This function loads different genomic regions created in the \code{prepare_annotation_files} step, 
#' separating features on different recognized organelles. The bam files is then analyzed in chunks to minimize RAM usage.\cr
#' The complete list of analysis and output is as follows:\cr\cr
#' \code{read_stats}: contains:\cr read length distribution (rld) per organelle, \code{positions} containes mapping statistics
#' on different genomic regions, \code{reads_pos1} contains 5' end mapping positions for each read, separated by read length.
#' \code{counts_cds_genes}: contains read mapping statistics on CDS regions of protein coding genes, including gene symbols, counts, RPKM and TPM values
#' \code{counts_all_genes}: is a similar object, but contains statistics on all annotated genes.
#' \code{reads_summary}: reports mapping statistics on different genomic regions and divided by read length and organelle.\cr\cr
#' \code{profiles_fivepr} contains:\cr
#' \code{five_prime_bins}: a DataFrame object (one for each read length and compartment) with signal values over 50 5'UTR bins,
#' 100 CDS bins and 50 3'UTR bins; one representative transcript (reprentative_mostcommon) is selected for each gene. \code{five_prime_subcodon} containes a similar structure, but for 25nt downstream the Transcription
#' Start Site (TSS), 25nt upstream start codons, 33nt donwstream the start codon, 33nt in the middle of the ORF, 33nt upstream the stop codon,
#' 25nt downstream the stop codon, and 25nt upstream the Transcription End Site (TES).\cr\cr
#' \code{selection_cutoffs} contains:\cr
#' \code{results_choice}: containing the calculated cutoffs and selected readlengths, together with \code{data} about the different 
#' methods. \code{results_cutoffs} has statistics about calculated cutoffs, while \code{analysis_frame_cutoff} has extensive
#' statistics concerning cutoff calculations and read length selection, see \code{calc_cutoffs_from_profiles} for more details.\cr\cr
#' \code{P_sites_stats}: contains the list of calculated P_sites, from all reads (P_sites_all), uniquely mapping reads (P_sites_all_uniq), 
#' or uniquely mapping reads with mismatches (P_sites_uniq_mm). \code{junctions} contains stastics on read mapping on annotated splice junctions.
#' coverage for entire reads (no 5'ends or P_sites-transformed) on different strands and for all and uniquely mapping reads are also calculated.\cr\cr
#' \code{profiles_P_sites} contains:\cr
#' \code{P_sites_bins}: profiles for each organelle and read length around binned transcript locations.\cr
#' \code{P_sites_subcodon}: profiles for each organelle and read length around transcript start/ends and ORF start/ends.\cr
#' \code{Codon_counts}: codon occurrences in the first 11 codons, middle 11 codons, and last 11 codons for each ORF.\cr
#' \code{P_sites_percodon}: P_sites counts on each codon, separated by ORF positions as described above. Values are separated by organelle and read length.\cr
#' \code{P_sites_percodon_ratio}: ratio of P_sites_percodon/Codon_counts, as a measure of P_site occupancy on each codon, divided again by organelle and read length, for different ORF positions.\cr\cr
#' \code{sequence_analysis}: contains a DataFrame object with the 50top mapping location in the genome, with the corresponding DNA sequence,
#' number of reads mapping (also in percentage of total n of reads), and genomic feature annotation.\cr\cr
#' \code{summary_P_sites}: contains a DataFrame object summarizing the P_sites calculation and read length selection, including statistics on percentage of total reads used.
#' @seealso \code{\link{prepare_annotation_files}}, \code{\link{calc_cutoffs_from_profiles}}, \code{\link{choose_readlengths}}, \code{\link{create_html_report}}.
#' @export

RiboseQC_analysis<-function(annotation_file,bam_files,read_subset=T,readlength_choice_method="max_coverage",
                            chunk_size=5000000L,write_tmp_files=T,dest_names=NA,rescue_all_rls=FALSE,fast_mode=T,create_report=T,sample_names=NA,report_file=NA,extended_report=F,pdf_plots=T){
    
    if(length(dest_names)==1){
        if(is.na(dest_names)){
            dest_names=bam_files
        }
    }
    
    if(sum(duplicated(dest_names))>0){
        stop(paste("duplicated 'dest_names' parameter (possibly same bam file). Please specify different 'dest_names'.",date(),"\n"))
    }
    
    if(is.na(report_file)){report_file=paste(bam_files[1],"_RiboseQC_report.html",sep="")}
    if(length(sample_names)==1){
        if(is.na(sample_names)){sample_names=paste("sample",1:length(bam_files),sep="")}
    }
    
    if(!is.logical(read_subset)){stop(paste("'read_subset' must be logical (TRUE or FALSE, no quotes).",date(),"\n"))}
    if(!is.logical(write_tmp_files)){stop(paste("'write_tmp_files' must be logical (TRUE or FALSE, no quotes).",date(),"\n"))}
    if(!is.logical(rescue_all_rls)){stop(paste("'rescue_all_rls' must be logical (TRUE or FALSE, no quotes).",date(),"\n"))}
    if(!is.logical(fast_mode)){stop(paste("'fast_mode' must be logical (TRUE or FALSE, no quotes).",date(),"\n"))}
    if(!is.logical(create_report)){stop(paste("'create_report' must be logical (TRUE or FALSE, no quotes).",date(),"\n"))}
    if(!is.integer(chunk_size) & !is.numeric(chunk_size)){stop(paste("'chunk_size' must be an integer value (e.g. 5000000).",date(),"\n"))}
    if(chunk_size<=100000 | chunk_size>100000000){stop(paste("'chunk_size' must be an integer value between 100000 and 100000000 (e.g. 5000000).",date(),"\n"))}
    if(!is.character(annotation_file)){stop(paste("'annotation_file' must be a character (e.g.'/home/myfolder/myfile').",date(),"\n"))}
    if(!is.character(bam_files)){stop(paste("'bam_files' must be a character (e.g.'/home/myfolder/myfile') or character vector (e.g. c('/home/myfolder/myfile1','/home/myfolder/myfile2').",date(),"\n"))}
    if(!is.character(dest_names)){stop(paste("'dest_names' must be a character (e.g.'/home/myfolder/myfile') or character vector (e.g. c('/home/myfolder/myfile1','/home/myfolder/myfile2')..",date(),"\n"))}
    if(!is.character(report_file)){stop(paste("'report_file' must be a character (e.g.'/home/myfolder/myfile').",date(),"\n"))}
    if(!is.character(sample_names)){stop(paste("'sample_names' must be a character (e.g.'failed_experiment'). or character vector (e.g. c('attempt_1','attempt_2').",date(),"\n"))}
    
    if(!length(read_subset)%in%c(1,length(bam_files))){stop(paste("length of 'read_subset' must be 1 or same as 'bam_files'.",date(),"\n"))}
    if(!length(readlength_choice_method)%in%c(1,length(bam_files))){stop(paste("length of 'readlength_choice_method' must be 1 or same as 'bam_files'.",date(),"\n"))}
    if(!length(rescue_all_rls)%in%c(1,length(bam_files))){stop(paste("length of 'rescue_all_rls' must be 1 or same as 'bam_files'.",date(),"\n"))}
    if(!length(fast_mode)%in%c(1,length(bam_files))){stop(paste("length of 'fast_mode' must be 1 or same as 'bam_files'.",date(),"\n"))}
    
    if(length(bam_files)>1){
        if(length(read_subset)==1){read_subset<-rep(read_subset,length(bam_files))}
        if(length(readlength_choice_method)==1){readlength_choice_method<-rep(readlength_choice_method,length(bam_files))}
        if(length(rescue_all_rls)==1){rescue_all_rls<-rep(rescue_all_rls,length(bam_files))}
        if(length(fast_mode)==1){fast_mode<-rep(fast_mode,length(bam_files))}
    }
    
    
    
    if(!length(dest_names)%in%c(length(bam_files))){stop(paste("length of 'dest_names' must be same as 'bam_files'.",date(),"\n"))}
    if(!length(sample_names)%in%c(length(bam_files))){stop(paste("length of 'sample_names' must be same as 'bam_files'.",date(),"\n"))}
    
    if(sum(!readlength_choice_method%in%c("max_coverage","max_inframe","all"))>0){stop(paste("'readlength_choice_method' must be one of 'max_coverage','max_inframe','all'.",date(),"\n"))}
    
    list_annotations<-list()
    
    for(annots in 1:length(annotation_file)){
        cat(paste("Loading annotation files in ",annotation_file[annots]," ... ", date(),"\n",sep = ""))
        
        load_annotation(annotation_file[annots])
        lst<-list(GTF_annotation,genome_seq)
        names(lst)<-c("GTF_annotation","genome_seq")
        list_annotations[[annots]]<-lst
        
    }
    
    cat(paste("Loading annotation files --- Done! ", date(),"\n",sep = ""))
    
    GTF_annotation<-list_annotations[[1]]$GTF_annotation
    genome_seq<-list_annotations[[1]]$genome_seq
    
    
    for(bammo in 1:length(bam_files)){
        
        chunk_size<-as.integer(chunk_size)
        
        GTF_annotation<-list_annotations[[1]]$GTF_annotation
        genome_seq<-list_annotations[[1]]$genome_seq
        
        if(length(annotation_file)>1 & bammo>1){
            GTF_annotation<-list_annotations[[bammo]]$GTF_annotation
            genome_seq<-list_annotations[[bammo]]$genome_seq
        }
        
        bam_file<-bam_files[bammo]
        dest_name<-dest_names[bammo]
        
        opts <- BamFile(file=bam_file, yieldSize=chunk_size) #read BAMfile in one chunk
        param <- ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),what=c("mapq"),tag = "MD")
        
        dira<-paste(dirname(dest_name),paste("tmp_RiboseQC",basename(dest_name),sep="_"),sep="/")
        suppressWarnings(dir.create(dira,recursive = T))
        
        seqs <- seqinfo(opts)
        circs_seq<-seqnames(GTF_annotation$seqinfo)[which(isCircular(GTF_annotation$seqinfo))]
        circs <- seqs@seqnames[which(seqs@seqnames%in%circs_seq)]
        
        red_cdss<-reduce(GTF_annotation$cds_genes)
        
        red_ex <- unlist(GTF_annotation$exons_txs)
        red_ex$gene_id<-GTF_annotation$trann$gene_id[match(names(red_ex),GTF_annotation$trann$transcript_id)]
        red_ex<-reduce(split(red_ex,red_ex$gene_id))
        
        regions <- list(reduce(unlist(GTF_annotation$cds_genes)), GTF_annotation$fiveutrs, GTF_annotation$threeutrs,
                        GTF_annotation$ncIsof, GTF_annotation$ncRNAs, GTF_annotation$introns, GTF_annotation$intergenicRegions)
        names(regions) <- c("cds","fiveutrs","threeutrs",
                            "ncIsof","ncRNAs","introns","intergenic")
        
        regions<-GRangesList(endoapply(regions,function(x){mcols(x)<-NULL;x}))
        
        list_locat <- list()
        list_locat[["nucl"]] <- regions
        all <- regions
        if(length(circs)>0){
            for(i in c("nucl", circs)){
                if(i=="nucl"){
                    regions<-all
                    for(w in 1:length(all)){
                        regions[[w]] <- regions[[w]][!seqnames(regions[[w]])%in%circs]
                    }
                    list_locat[["nucl"]] <- regions
                }
                if(i!="nucl"){
                    regions <- all
                    for(w in 1:length(all)){
                        regions[[w]] <- regions[[w]][seqnames(regions[[w]])%in%i]
                    }
                    list_locat[[i]] <- regions
                }
            }
        }
        
        #get readlengths
        
        reduc <- function(x,y){
            list(max(x[[1]],y[[1]]),min(x[[2]],y[[2]]))
        }
        yiel <- function(x){
            readGAlignments(x,param = param)
        }
        
        mapp <- function(x){
            x_I<-x[grep("I",cigar(x))]
            
            if(length(x_I)>0){
                x<-x[grep("I",cigar(x),invert=T)]
                
            }
            x_D<-x[grep("D",cigar(x))]
            if(length(x_D)>0){
                x<-x[grep("D",cigar(x),invert=T)]
                
            }
            
            
            clipp <- width(cigarRangesAlongQuerySpace(x@cigar, ops="S"))
            clipp[elementNROWS(clipp)==0] <- 0
            len_adj <- qwidth(x)-sum(clipp)
            mcols(x)$len_adj <- len_adj
            
            
            maxr<-max(mcols(x)$len_adj)
            minr<-min(mcols(x)$len_adj)
            list(maxr,minr)
        }
        
        cat(paste("Extracting read lengths from ",bam_file," ... ", date(),"\n",sep=""))
        maxmin <- reduceByYield(X=opts, YIELD=yiel, MAP=mapp, REDUCE=reduc)
        cat(paste("Extracting read lengths --- Done! ", date(),"\n",sep=""))
        readlengths<-seq(maxmin[[2]],maxmin[[1]])
        
        
        
        # body of the iteration for reducebyYield
        
        
        # what to do with chunks: x present chunk, y old chunks (cumulative)
        reduc <- function(x,y){
            all_ps<-GRangesList()
            rls<-unique(c(names(x[["reads_pos1"]]),names(y[["reads_pos1"]])))
            seql<-seqlevels(GTF_annotation$seqinfo)
            seqle<-seqlengths(GTF_annotation$seqinfo)
            for(rl in rls){
                reads_x<-GRanges()
                reads_y<-GRanges()
                seqlevels(reads_x)<-seql
                seqlevels(reads_y)<-seql
                
                seqlengths(reads_x)<-seqle
                seqlengths(reads_y)<-seqle
                
                if(sum(rl%in%names(x[["reads_pos1"]]))>0){reads_x<-x[["reads_pos1"]][[rl]]}
                if(sum(rl%in%names(y[["reads_pos1"]]))>0){reads_y<-y[["reads_pos1"]][[rl]]}
                plx<-reads_x[strand(reads_x)=="+"]
                mnx<-reads_x[strand(reads_x)=="-"]
                ply<-reads_y[strand(reads_y)=="+"]
                mny<-reads_y[strand(reads_y)=="-"]
                if(length(plx)>0){covv_pl<-coverage(plx,weight = plx$score)}else{covv_pl<-coverage(plx)}
                if(length(ply)>0){covv_pl<-covv_pl+coverage(ply,weight = ply$score)}
                
                covv_pl<-GRanges(covv_pl)
                covv_pl<-covv_pl[covv_pl$score>0]
                
                if(length(mnx)>0){covv_min<-coverage(mnx,weight = mnx$score)}else{covv_min<-coverage(mnx)}
                if(length(mny)>0){covv_min<-covv_min+coverage(mny,weight = mny$score)}
                
                covv_min<-GRanges(covv_min)
                covv_min<-covv_min[covv_min$score>0]
                
                strand(covv_pl)<-"+"
                strand(covv_min)<-"-"
                
                all_ps[[rl]]<-sort(c(covv_pl,covv_min))
            }
            
            
            cntsss <- y[["counts_cds_genes"]]
            cntsss$reads <- cntsss$reads+x[["counts_cds_genes"]]$reads
            
            cntsss_all <- y[["counts_all_genes"]]
            cntsss_all$reads <- cntsss_all$reads+x[["counts_all_genes"]]$reads
            
            
            cntsss_unq <- y[["counts_cds_genes_unq"]]
            cntsss_unq$reads <- cntsss_unq$reads+x[["counts_cds_genes_unq"]]$reads
            
            cntsss_all_unq <- y[["counts_all_genes_unq"]]
            cntsss_all_unq$reads <- cntsss_all_unq$reads+x[["counts_all_genes_unq"]]$reads
            
            
            reads_summary <- mapply(FUN = function(x,y){DataFrame(as.matrix(x)+as.matrix(y))}, x[["reads_summary"]], y[["reads_summary"]], SIMPLIFY = FALSE)
            reads_summary_unq <- mapply(FUN = function(x,y){DataFrame(as.matrix(x)+as.matrix(y))}, x[["reads_summary_unq"]], y[["reads_summary_unq"]], SIMPLIFY = FALSE)
            
            lis <- list((x[["rld"]]+y[["rld"]]),x[["rld_unq"]]+y[["rld_unq"]], (x[["positions"]]+y[["positions"]]),(x[["positions_unq"]]+y[["positions_unq"]]), all_ps,cntsss,cntsss_unq, cntsss_all,cntsss_all_unq, reads_summary,reads_summary_unq)
            names(lis) <- c("rld","rld_unq", "positions","positions_unq" ,"reads_pos1", "counts_cds_genes","counts_cds_genes_unq", "counts_all_genes","counts_all_genes_unq", "reads_summary", "reads_summary_unq")
            lis
        }
        
        # what to do with each chunk (read as alignment file)
        yiel <- function(x){
            readGAlignments(x,param = param)
        }
        
        # operations on the chunk (here count reads and whatnot)
        mapp <- function(x){
            
            # Remove Insertion and Deletion (TBA)
            
            x_I<-x[grep("I",cigar(x))]
            
            if(length(x_I)>0){
                x<-x[grep("I",cigar(x),invert=T)]
                
            }
            x_D<-x[grep("D",cigar(x))]
            if(length(x_D)>0){
                x<-x[grep("D",cigar(x),invert=T)]
                
            }
            
            emptyt <- rep(0,length(readlengths))
            names(emptyt) <- paste("reads",readlengths,sep="_")
            
            
            # softclipping part
            
            
            clipp <- width(cigarRangesAlongQuerySpace(x@cigar, ops="S"))
            clipp[elementNROWS(clipp)==0] <- 0
            len_adj <- qwidth(x)-sum(clipp)
            mcols(x)$len_adj <- len_adj
            
            # Remove S from Cigar (read positions/length are already adjusted)
            # it helps calculating P-sites positions for spliced reads
            
            cigg<-cigar(x)
            cigg_s<-grep(cigg,pattern = "S")
            if(length(cigg_s)>0){
                cigs<-cigg[cigg_s]
                cigs<-gsub(cigs,pattern = "^[0-9]+S",replacement = "")
                cigs<-gsub(cigs,pattern = "[0-9]+S$",replacement = "")
                cigg[cigg_s]<-cigs
                x@cigar<-cigg
            }
            mcols(x)$cigar_str<-x@cigar
            x_uniq<-x[x@elementMetadata$mapq>50]
            
            
            
            # rld_len, rld_loc
            
            
            # initialise read length vector per compartment
            list_vects <- list(nucl=emptyt)
            if(length(circs)>0){
                for(i in circs){
                    list_vects[[i]] <- emptyt            
                }
            }
            
            list_vects_unq <- list(nucl=emptyt)
            if(length(circs)>0){
                for(i in circs){
                    list_vects_unq[[i]] <- emptyt            
                }
            }
            
            # sort reads per compartment
            list_reads <- list(nucl=x)
            list_reads_unq <- list(nucl=x_uniq)
            if(length(circs)>0){
                nucl <- subset(x,!seqnames(x)%in%circs)
                nucl_unq <- subset(x_uniq,!seqnames(x_uniq)%in%circs)
                list_reads <- list(nucl=nucl)
                for(i in circs){
                    list_reads[[i]] <- subset(x,seqnames(x)==i)
                    list_reads_unq[[i]] <- subset(x_uniq,seqnames(x_uniq)==i)            
                }
            }
            
            # count reads per read length per compartment
            
            for(i in names(list_reads)){
                tab <- table(mcols(list_reads[[i]])$len_adj)
                tab_unq <- table(mcols(list_reads_unq[[i]])$len_adj)
                for(j in readlengths){
                    cont <- as.numeric(tab[which(names(tab)==j)])
                    cont_unq <- as.numeric(tab_unq[which(names(tab_unq)==j)])
                    if(length(cont)>0){list_vects[[i]][j-(min(readlengths)-1)] <- cont}
                    if(length(cont_unq)>0){list_vects_unq[[i]][j-(min(readlengths)-1)] <- cont_unq}
                    
                }
            }
            
            rld_len  <- do.call(what=rbind.data.frame, list_vects)
            names(rld_len) <- paste("reads", readlengths, sep="_")
            row.names(rld_len) <- names(list_reads)
            
            rld_len_unq  <- do.call(what=rbind.data.frame, list_vects_unq)
            names(rld_len_unq) <- paste("reads", readlengths, sep="_")
            row.names(rld_len_unq) <- names(list_reads)
            
            
            # count reads per biotype location per compartment
            list_rld_loc <- list()
            for(i in c("nucl", circs)){
                list_rld_loc[[i]] <- (assay(summarizeOverlaps(reads=x,
                                                              features=GRangesList(list_locat[[i]]),
                                                              ignore.strand=F,
                                                              mode="Union",
                                                              inter.feature=FALSE)))
            }
            rld_loc <- do.call(what=cbind.data.frame, list_rld_loc)
            names(rld_loc) <- paste("reads", names(list_rld_loc), sep="_")
            
            
            list_rld_loc_unq <- list()
            for(i in c("nucl", circs)){
                list_rld_loc_unq[[i]] <- (assay(summarizeOverlaps(reads=x_uniq,
                                                                  features=GRangesList(list_locat[[i]]),
                                                                  ignore.strand=F,
                                                                  mode="Union",
                                                                  inter.feature=FALSE)))
            }
            rld_loc_unq <- do.call(what=cbind.data.frame, list_rld_loc_unq)
            names(rld_loc_unq) <- paste("reads", names(list_rld_loc_unq), sep="_")
            
            
            # reads_summary
            
            reads_summary <- list()
            
            # iterate over genome types and compartments
            for(i in names(list_reads)){
                reads_by_length <- list()
                rdss<-list_reads[[i]]
                rdss<-split(rdss,mcols(rdss)$len_adj)
                
                rgs<-GRangesList(list_locat[[i]])
                ovs<-lapply(rdss,FUN = function(x){
                    suppressWarnings(assay(summarizeOverlaps(reads=x,features=rgs,ignore.strand=F,mode="Union",inter.feature=FALSE)))
                })
                nope<-readlengths[which(!readlengths%in%names(ovs))]
                if(length(nope)>0){
                    for(j in nope){
                        nop<-matrix(0,nrow = 7,ncol = 1)
                        rownames(nop)<-names(rgs)
                        ovs[[as.character(j)]]<-nop
                        
                    }
                }
                ovs<-ovs[names(ovs)%in%as.character(readlengths)]
                ovs<-ovs[as.character(readlengths)]
                ovs<-do.call(ovs,what = cbind)
                colnames(ovs)<-paste("reads",readlengths,sep = "_")
                reads_summary[[i]] <- DataFrame(ovs)
                
                
            }
            
            
            
            reads_summary_unq <- list()
            
            # iterate over genome types and compartments
            for(i in names(list_reads_unq)){
                reads_by_length <- list()
                rdss<-list_reads_unq[[i]]
                rdss<-split(rdss,mcols(rdss)$len_adj)
                
                rgs<-GRangesList(list_locat[[i]])
                ovs<-lapply(rdss,FUN = function(x){
                    suppressWarnings(assay(summarizeOverlaps(reads=x,features=rgs,ignore.strand=F,mode="Union",inter.feature=FALSE)))
                })
                nope<-readlengths[which(!readlengths%in%names(ovs))]
                if(length(nope)>0){
                    for(j in nope){
                        nop<-matrix(0,nrow = 7,ncol = 1)
                        rownames(nop)<-names(rgs)
                        ovs[[as.character(j)]]<-nop
                        
                    }
                }
                ovs<-ovs[names(ovs)%in%as.character(readlengths)]
                ovs<-ovs[as.character(readlengths)]
                ovs<-do.call(ovs,what = cbind)
                colnames(ovs)<-paste("reads",readlengths,sep = "_")
                reads_summary_unq[[i]] <- DataFrame(ovs)
                
                
            }
            
            
            # reads_pos1_stst
            
            
            # take first position of each read
            reads_pos1 <- unlist(resize(split(GRanges(x), f=mcols(x)$len_adj), 1))
            
            reads_pos1<-split(reads_pos1,reads_pos1$len_adj)
            
            reads_pos1<-GRangesList(lapply(reads_pos1,function(y){
                unq<-unique(y)
                mcols(unq)<-NULL
                unq$score<-countOverlaps(unq,y,type="equal")
                unq
                
            }))
            
            # cnts_cds_genes, cnts_all_genes
            
            seqgene <- cbind(as.vector(seqnames(GTF_annotation$genes)),
                             GTF_annotation$genes$gene_id)
            
            
            cnts_cds_genes <- assay(summarizeOverlaps(reads=x, features=red_cdss, 
                                                      ignore.strand=F, mode="Union", inter.feature=FALSE))
            
            cnts_all_genes <- assay(summarizeOverlaps(reads=x, features=red_ex, 
                                                      ignore.strand=F, mode="Union", inter.feature=FALSE))
            
            cnts_cds_genes_unq <- assay(summarizeOverlaps(reads=x_uniq, features=red_cdss, 
                                                          ignore.strand=F, mode="Union", inter.feature=FALSE))
            
            cnts_all_genes_unq <- assay(summarizeOverlaps(reads=x_uniq, features=red_ex, 
                                                          ignore.strand=F, mode="Union", inter.feature=FALSE))
            
            
            chrsss <- seqgene[match(rownames(cnts_cds_genes), seqgene[,2]), 1]
            
            cnts_cds_genes <- DataFrame(cnts_cds_genes)
            cnts_cds_genes$chromosome <- chrsss
            cnts_cds_genes$gene_name <- GTF_annotation$trann$gene_name[match(rownames(cnts_cds_genes), GTF_annotation$trann$gene_id)]
            cnts_cds_genes$gene_biotype <- GTF_annotation$trann$gene_biotype[match(rownames(cnts_cds_genes), GTF_annotation$trann$gene_id)]
            
            cnts_all_genes <- DataFrame(cnts_all_genes)
            cnts_all_genes$chromosome <- seqgene[match(rownames(cnts_all_genes), seqgene[,2]), 1]
            cnts_all_genes$gene_name <- GTF_annotation$trann$gene_name[match(rownames(cnts_all_genes), GTF_annotation$trann$gene_id)]
            cnts_all_genes$gene_biotype <- GTF_annotation$trann$gene_biotype[match(rownames(cnts_all_genes), GTF_annotation$trann$gene_id)]
            
            
            cnts_cds_genes_unq <- DataFrame(cnts_cds_genes_unq)
            cnts_cds_genes_unq$chromosome <- chrsss
            cnts_cds_genes_unq$gene_name <- GTF_annotation$trann$gene_name[match(rownames(cnts_cds_genes_unq), GTF_annotation$trann$gene_id)]
            cnts_cds_genes_unq$gene_biotype <- GTF_annotation$trann$gene_biotype[match(rownames(cnts_cds_genes_unq), GTF_annotation$trann$gene_id)]
            
            cnts_all_genes_unq <- DataFrame(cnts_all_genes_unq)
            cnts_all_genes_unq$chromosome <- seqgene[match(rownames(cnts_all_genes_unq), seqgene[,2]), 1]
            cnts_all_genes_unq$gene_name <- GTF_annotation$trann$gene_name[match(rownames(cnts_all_genes_unq), GTF_annotation$trann$gene_id)]
            cnts_all_genes_unq$gene_biotype <- GTF_annotation$trann$gene_biotype[match(rownames(cnts_all_genes_unq), GTF_annotation$trann$gene_id)]
            
            
            
            chunk_res <- list(rld_len,rld_len_unq, rld_loc, rld_loc_unq, reads_pos1, 
                              cnts_cds_genes,cnts_cds_genes_unq, cnts_all_genes, cnts_all_genes_unq,
                              reads_summary,reads_summary_unq)
            # export list to use in the reduc part above
            names(chunk_res) <- c("rld","rld_unq", "positions","positions_unq" ,"reads_pos1", "counts_cds_genes","counts_cds_genes_unq", "counts_all_genes","counts_all_genes_unq", "reads_summary", "reads_summary_unq")
            
            chunk_res
        }
        
        
        # main for reducebyYield
        
        cat(paste("Analyzing BAM file:",bam_file,"...", date(),"\n"))
        
        res_1 <- reduceByYield(X=opts, YIELD=yiel, MAP=mapp, REDUCE=reduc)
        names(res_1) <- c("rld","rld_unq", "positions","positions_unq" ,"reads_pos1", "counts_cds_genes","counts_cds_genes_unq", "counts_all_genes","counts_all_genes_unq", "reads_summary", "reads_summary_unq")
        
        cds_cnts<-res_1$counts_cds_genes
        cds_cnts$RPKM<-cds_cnts$reads/(sum(cds_cnts$reads)/1e06)
        cds_cnts$RPKM<-round(cds_cnts$RPKM/(sum(width(red_cdss))/1000),digits = 4)
        cds_cnts$TPM<-cds_cnts$reads/(sum(width(red_cdss))/1000)
        cds_cnts$TPM<-round(cds_cnts$TPM/(sum(cds_cnts$TPM)/1e06),digits = 4)
        res_1$counts_cds_genes<-cds_cnts
        
        cds_cnts_unq<-res_1$counts_cds_genes_unq
        cds_cnts_unq$RPKM<-cds_cnts_unq$reads/(sum(cds_cnts_unq$reads)/1e06)
        cds_cnts_unq$RPKM<-round(cds_cnts_unq$RPKM/(sum(width(red_cdss))/1000),digits = 4)
        cds_cnts_unq$TPM<-cds_cnts_unq$reads/(sum(width(red_cdss))/1000)
        cds_cnts_unq$TPM<-round(cds_cnts_unq$TPM/(sum(cds_cnts_unq$TPM)/1e06),digits = 4)
        res_1$counts_cds_genes_unq<-cds_cnts_unq
        
        ex_cnts<-res_1$counts_all_genes
        ex_cnts$RPKM<-ex_cnts$reads/(sum(ex_cnts$reads)/1e06)
        ex_cnts$RPKM<-round(ex_cnts$RPKM/(sum(width(red_ex))/1000),digits = 4)
        ex_cnts$TPM<-ex_cnts$reads/(sum(width(red_ex))/1000)
        ex_cnts$TPM<-round(ex_cnts$TPM/(sum(ex_cnts$TPM)/1e06),digits = 4)
        res_1$counts_all_genes<-ex_cnts
        
        ex_cnts_unq<-res_1$counts_all_genes_unq
        ex_cnts_unq$RPKM<-ex_cnts_unq$reads/(sum(ex_cnts_unq$reads)/1e06)
        ex_cnts_unq$RPKM<-round(ex_cnts_unq$RPKM/(sum(width(red_ex))/1000),digits = 4)
        ex_cnts_unq$TPM<-ex_cnts_unq$reads/(sum(width(red_ex))/1000)
        ex_cnts_unq$TPM<-round(ex_cnts_unq$TPM/(sum(ex_cnts_unq$TPM)/1e06),digits = 4)
        res_1$counts_all_genes_unq<-ex_cnts_unq
        
        save(res_1,file = paste(dira,"res_1",sep = "/"))
        cat(paste("Analyzing BAM file --- Done!", date(),"\n"))
        
        ps<-res_1$reads_pos1
        
        tile_cds<-GTF_annotation$cds_txs_coords
        tile_cds<-tile_cds[tile_cds$reprentative_mostcommon]
        
        list_txs_ok<-GRangesList()
        for(ci in circs){
            citxs<-GTF_annotation$cds_txs[seqnames(GTF_annotation$cds_txs)==ci]
            citxs<-citxs[elementNROWS(citxs)>0]
            citxs_txs<-GTF_annotation$cds_txs_coords[as.vector(seqnames(GTF_annotation$cds_txs_coords))%in%names(citxs)]
            list_txs_ok[[ci]]<-citxs_txs
        }
        list_txs_ok[["nucl"]]<-tile_cds[!tile_cds%in%unlist(list_txs_ok)]
        
        cat(paste("Building aggregate 5' profiles ...", date(),"\n"))
        ps_signals_win_all<-list()
        ps_signals_tiles_all<-list()
        
        #calculate profiles for bins and windows (single nt resolution)
        for(comp in c("nucl",circs)){
            tile_cds<-list_txs_ok[[comp]]
            strand(tile_cds)<-"+"
            tile_cds<-tile_cds[width(tile_cds)>100]
            checkgen<-TRUE
            if(comp=="nucl"){
                if((sum((start(tile_cds)>50))/length(start(tile_cds)))>.33){tile_cds<-tile_cds[start(tile_cds)>50];checkgen<-FALSE}
                if(sum(seqlengths(tile_cds)[as.vector(seqnames(tile_cds))]-end(tile_cds)>100)/length(tile_cds)>.33){
                    tile_cds<-tile_cds[seqlengths(tile_cds)[as.vector(seqnames(tile_cds))]-end(tile_cds)>100];checkgen<-FALSE}
                ps_comp<-GRangesList(lapply(ps,function(x){sort(x[!seqnames(x)%in%circs])}))
            } else {ps_comp<-GRangesList(lapply(ps,function(x){sort(x[seqnames(x)==comp])}))}
            
            #select only some read lenghts, otherwise too much output and computation
            
            summm<-res_1$reads_summary
            a<-unlist(summm[[comp]]["cds",])+unlist(summm[[comp]]["fiveutrs",])+unlist(summm[[comp]]["threeutrs",])
            names(a)<-gsub(names(a),pattern = "reads_",replacement = "")
            a<-sort(a/(sum(a)/100),decreasing = T)
            cuma<-cumsum(a)
            if(length(cuma)==0){
                ps_signals_tiles_all[[comp]]<-c()
                ps_signals_win_all[[comp]]<-c()
                next
            }
            #if read subset, discard low abundance
            if(read_subset[bammo]==T){
                rl_ok<-names(cuma[1:which(cuma>99)[1]])
                #rl_ok<-as.character(seq(sort(as.numeric(rl_ok),decreasing = F)[1],sort(as.numeric(rl_ok),decreasing = T)[1]))
                ps_comp<-ps_comp[names(ps_comp)%in%rl_ok]
            }
            
            signal_ps<-list()
            signal_ps_nt<-list()
            
            all_ps<-unlist(ps_comp)
            ps_comp[["all"]]<-all_ps
            mp<-mapToTranscripts(all_ps,transcripts = GTF_annotation$exons_txs[as.vector(seqnames(tile_cds))])
            mp$score<-all_ps$score[mp$xHits]
            #put seqlevels for compartment only to speed up stuff
            seqlevels(mp)<-as.vector(seqnames(tile_cds))
            seqlengths(mp)<-sum(width(GTF_annotation$exons_txs[as.vector(seqnames(tile_cds))]))
            covtx<-coverage(mp,weight = mp$score)
            ok_txs<-names(covtx[elementNROWS(runValue(covtx))>1])
            if(length(ok_txs)<5){
                ps_signals_tiles_all[[comp]]<-c()
                ps_signals_win_all[[comp]]<-c()
                next
            }
            #keep all txs if small n of txs
            
            if(length(ok_txs)>=5 & length(covtx)<250){
                ok_txs<-names(covtx)
            }
            
            if(length(ok_txs)>5000){
                cnts_txss_agg<-sort(sum(covtx),decreasing = T)
                ok_txs<-names(cnts_txss_agg)[1:5000]
            }
            if(fast_mode[bammo]==T){
                if(length(ok_txs)>500){
                    cnts_txss_agg<-sort(sum(covtx),decreasing = T)
                    ok_txs<-names(cnts_txss_agg)[1:500]
                }
            }
            
            
            
            fivs<-GRanges(seqnames = seqnames(GTF_annotation$cds_txs_coords),ranges = IRanges(start=1,end = start(GTF_annotation$cds_txs_coords)),strand="*")
            fivs<-fivs[as.character(seqnames(fivs))%in%ok_txs]
            fivs_gen<-unlist(pmapFromTranscripts(fivs,transcripts = GTF_annotation$exons_txs[as.character(seqnames(fivs))],ignore.strand=F))
            fivs_gen<-fivs_gen[fivs_gen$hit]
            fivs_gen<-split(fivs_gen,names(fivs_gen))
            
            
            threes<-GRanges(seqnames = seqnames(GTF_annotation$cds_txs_coords),ranges = IRanges(start=end(GTF_annotation$cds_txs_coords),end = GTF_annotation$cds_txs_coords$lentx),strand="*")
            threes<-threes[as.character(seqnames(threes))%in%ok_txs]
            threes_gen<-unlist(pmapFromTranscripts(threes,transcripts = GTF_annotation$exons_txs[as.character(seqnames(threes))],ignore.strand=F))
            threes_gen<-threes_gen[threes_gen$hit]
            threes_gen<-split(threes_gen,names(threes_gen))
            
            cds_gen<-GTF_annotation$cds_txs[ok_txs]
            
            
            tile_cds<-tile_cds[as.character(seqnames(tile_cds))%in%ok_txs]
            ok_txs<-unique(as.character(seqnames(tile_cds)))
            seqlevels(tile_cds)<-ok_txs
            list_covs<-list()
            list_covs[["all"]]<-covtx
            
            tile_5<-GRanges(seqnames(tile_cds),ranges = IRanges(start=1,end = start(tile_cds)))
            tile_3<-GRanges(seqnames(tile_cds),ranges = IRanges(start = end(tile_cds),end = seqlengths(tile_cds)[as.vector(seqnames(tile_cds))]))
            no5utr<-width(tile_5)<51
            no3utr<-width(tile_3)<51
            
            ex_annot<-GTF_annotation$exons_txs
            #if no utrs put bins to 1nt
            
            if(sum(no5utr)>0){
                tx_notok<-seqnames(tile_5)[no5utr]
                annot_notok<-ex_annot[tx_notok]
                annot_ok<-GRangesList(lapply(annot_notok, function(x){
                    if(length(x)==0){return(x)}
                    x[1]<-resize(x[1],width = width(x[1])+51,fix = "end")
                    x
                }))
                ex_annot[names(annot_ok)]<-annot_ok
                
                seqlengths(tile_cds)[as.vector(tx_notok)]<-sum(width(annot_ok))
                tile_cds[no5utr]<-shift(tile_cds[no5utr],+51)
                tile_5<-GRanges(seqnames(tile_cds),ranges = IRanges(start=1,end = start(tile_cds)))
                tile_3<-GRanges(seqnames(tile_cds),ranges = IRanges(start = end(tile_cds),end = seqlengths(tile_cds)[as.vector(seqnames(tile_cds))]))
                
                mp<-mapToTranscripts(all_ps,transcripts = ex_annot[as.vector(seqnames(tile_cds))])
                seqlevels(mp)<-seqlevels(tile_cds)
                seqlengths(mp)<-seqlengths(tile_cds)
                mp$score<-all_ps$score[mp$xHits]
                
                covtx<-coverage(mp,weight = mp$score)
                list_covs[["all"]]<-covtx
            }
            
            if(sum(no3utr)>0){
                tx_notok<-seqnames(tile_3)[no3utr]
                annot_notok<-ex_annot[tx_notok]
                annot_ok<-GRangesList(lapply(annot_notok, function(x){
                    if(length(x)==0){return(x)}
                    x[length(x)]<-resize(x[length(x)],width = width(x[length(x)])+51,fix = "start")
                    x
                }))
                ex_annot[names(annot_ok)]<-annot_ok
                
                seqlengths(tile_cds)[as.vector(tx_notok)]<-sum(width(annot_ok))
                
                tile_5<-GRanges(seqnames(tile_cds),ranges = IRanges(start=1,end = start(tile_cds)))
                tile_3<-GRanges(seqnames(tile_cds),ranges = IRanges(start = end(tile_cds),end = seqlengths(tile_cds)[as.vector(seqnames(tile_cds))]))
                
                mp<-mapToTranscripts(all_ps,transcripts = ex_annot[as.vector(seqnames(tile_cds))])
                seqlevels(mp)<-seqlevels(tile_cds)
                seqlengths(mp)<-seqlengths(tile_cds)
                mp$score<-all_ps$score[mp$xHits]
                
                covtx<-coverage(mp,weight = mp$score)
                list_covs[["all"]]<-covtx
                
            }
            
            ps_tiles<-DataFrameList()
            ps_win<-DataFrameList()
            
            for(len in c("all",names(ps_comp))){
                
                if((sum(no3utr) + sum(no5utr))==0){
                    mp<-mapToTranscripts(ps_comp[[len]],transcripts = fivs_gen)
                    mp$score<-ps_comp[[len]]$score[mp$xHits]
                    seqlevels(mp)<-names(fivs_gen)
                    seqlengths(mp)<-sum(width(fivs_gen))
                    cov_5<-coverage(mp,weight = mp$score)
                    
                    
                    mp<-mapToTranscripts(ps_comp[[len]],transcripts = threes_gen)
                    mp$score<-ps_comp[[len]]$score[mp$xHits]
                    seqlevels(mp)<-names(threes_gen)
                    seqlengths(mp)<-sum(width(threes_gen))
                    cov_3<-coverage(mp,weight = mp$score)
                    
                    mp<-mapToTranscripts(ps_comp[[len]],transcripts = cds_gen)
                    mp$score<-ps_comp[[len]]$score[mp$xHits]
                    seqlevels(mp)<-names(cds_gen)
                    seqlengths(mp)<-sum(width(cds_gen))
                    cov_cds<-coverage(mp,weight = mp$score)
                }
                
                if((sum(no3utr) + sum(no5utr))>0){
                    if(len!="all"){
                        mp<-mapToTranscripts(ps_comp[[len]],transcripts = ex_annot[as.vector(seqnames(tile_cds))])
                        mp$score<-ps_comp[[len]]$score[mp$xHits]
                        seqlevels(mp)<-seqlevels(tile_cds)
                        seqlengths(mp)<-seqlengths(tile_cds)
                        covtx<-coverage(mp,weight = mp$score)
                        covtx<-covtx[ok_txs]
                    }
                    
                    cov_5<-covtx[tile_5]
                    cov_3<-covtx[tile_3]
                    cov_cds<-covtx[tile_cds]
                    
                    
                }
                
                ps_tiles_5<-DataFrame(t(sapply(cov_5,function(x){
                    clos<-50*(round(length(x)/50,digits = 0)+1)
                    idx<-as.integer(seq(1,length(x),length.out = clos))
                    colMeans(matrix(x[idx],ncol = 50))
                })))
                ps_win_5<-DataFrame(t(sapply(cov_5,function(x){as.vector(x)[c(1:25,(length(x)-25):(length(x)-1))]})))
                
                ps_tiles_3<-DataFrame(t(sapply(cov_3,function(x){
                    clos<-50*(round(length(x)/50,digits = 0)+1)
                    idx<-as.integer(seq(1,length(x),length.out = clos))
                    colMeans(matrix(x[idx],ncol = 50))
                })))
                ps_win_3<-DataFrame(t(sapply(cov_3,function(x){as.vector(x)[c(2:26,(length(x)-24):length(x))]})))
                
                
                ps_tiles_cds<-DataFrame(t(sapply(cov_cds,function(x){
                    clos<-100*(round(length(x)/100,digits = 0)+1)
                    idx<-as.integer(seq(1,length(x),length.out = clos))
                    colMeans(matrix(x[idx],ncol = 100))
                })))
                ps_win_cds<-DataFrame(t(sapply(cov_cds,function(x){
                    rnd<-as.integer(length(x)/2)%%3
                    mid<-(as.integer(length(x)/2)-rnd)
                    as.vector(x)[c(1:33,(mid-17):(mid+15),(length(x)-32):length(x))]
                })))
                tls<-cbind(ps_tiles_5,ps_tiles_cds,ps_tiles_3)
                colnames(tls)<-c(paste("5_UTR",1:length(ps_tiles_5[1,]),sep="_"),paste("CDS",1:length(ps_tiles_cds[1,]),sep="_"),paste("3_UTR",1:length(ps_tiles_3[1,]),sep="_"))
                nts<-cbind(ps_win_5,ps_win_cds,ps_win_3)
                colnames(nts)<-c(paste("5_UTR",1:length(ps_win_5[1,]),sep="_"),paste("CDS",1:length(ps_win_cds[1,]),sep="_"),paste("3_UTR",1:length(ps_win_3[1,]),sep="_"))
                #DataFrameList?
                ps_tiles[[len]]<-tls
                ps_win[[len]]<-nts
                
            }
            
            
            ps_signals_tiles_all[[comp]]<-ps_tiles
            ps_signals_win_all[[comp]]<-ps_win
        }
        
        res_2<-list(ps_signals_tiles_all,ps_signals_win_all)
        names(res_2)<-c("five_prime_bins","five_prime_subcodon")
        save(res_2,file = paste(dira,"res_2",sep = "/"))
        
        cat(paste("Building aggregate 5' profiles --- Done!", date(),"\n"))
        
        cat(paste("Calculating P-sites offsets ...", date(),"\n"))
        
        
        #calculate cutoffs for each read length
        
        list_cutoff<-lapply(res_2[["five_prime_subcodon"]],function(x){
            list_res<-list()
            for(i in names(x)){
                if(i=="all"){lnmax=100}else{lnmax=as.numeric(i)}
                list_res[[i]]<-calc_cutoffs_from_profiles(x[[i]],length_max = lnmax)
                
            }
            list_res
        })
        
        summary_res<-lapply(list_cutoff,function(x) t(sapply(x,function (y) {
            res<-c(y$final_cutoff,y$frames_res[2])
            names(res)<-c("cutoff","frame_preference")
            res
        })))
        
        
        #choose readlengths
        
        res_rls<-choose_readlengths(summary_res,choice=readlength_choice_method[bammo],nt_signals=res_2[["five_prime_subcodon"]])
        
        cat(paste("Calculating P-sites offsets --- Done!", date(),"\n"))
        res_3<-list(res_rls,summary_res,list_cutoff)
        names(res_3)<-c("results_choice","results_cutoffs","analysis_frame_cutoff")
        save(res_3,file = paste(dira,"res_3",sep = "/"))
        
        
        
        #extract P-sites positions given selection above
        
        param <- ScanBamParam(flag=scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),what=c("mapq"),tag = "MD")
        seqllll<-seqlevels(unlist(unlist(res_1$reads_pos1)))
        seqleee<-seqlengths(unlist(unlist(res_1$reads_pos1)))
        rl_cutoffs_comp<-lapply(res_rls, function(x){x$final_choice})
        
        #rescue all read lengths
        if(rescue_all_rls[bammo]==T){
            all_rlsss<-as.numeric(names(res_1$reads_pos1))
            for(rilo in c("nucl",circs)){
                ralo<-rl_cutoffs_comp[[rilo]]
                if(is.null(ralo)){ralo<-DataFrame()}
                rlno2<-all_rlsss[!all_rlsss%in%ralo$read_length]
                ralo2<-DataFrame(read_length=rlno2,cutoff=12)
                rl_cutoffs_comp[[rilo]]<-rbind(ralo,ralo2)
            }
        }
        
        
        #what to do with chunks: x present chunk, y old chunks (cumulative)
        reduc<-function(x,y){
            #adjust merging by rls
            all_ps<-GRangesList()
            rls<-unique(c(names(x[["P_sites_all"]]),names(y[["P_sites_all"]])))
            
            for(rl in rls){
                reads_x<-GRanges()
                reads_y<-GRanges()
                
                seqlevels(reads_x)<-seqllll
                seqlevels(reads_y)<-seqllll
                
                seqlengths(reads_x)<-seqleee
                seqlengths(reads_y)<-seqleee
                
                if(sum(rl%in%names(x[["P_sites_all"]]))>0){reads_x<-x[["P_sites_all"]][[rl]]}
                if(sum(rl%in%names(y[["P_sites_all"]]))>0){reads_y<-y[["P_sites_all"]][[rl]]}
                
                plx<-reads_x[strand(reads_x)=="+"]
                mnx<-reads_x[strand(reads_x)=="-"]
                ply<-reads_y[strand(reads_y)=="+"]
                mny<-reads_y[strand(reads_y)=="-"]
                if(length(plx)>0){covv_pl<-coverage(plx,weight = plx$score)}else{covv_pl<-coverage(plx)}
                if(length(ply)>0){covv_pl<-covv_pl+coverage(ply,weight = ply$score)}
                
                covv_pl<-GRanges(covv_pl)
                covv_pl<-covv_pl[covv_pl$score>0]
                
                if(length(mnx)>0){covv_min<-coverage(mnx,weight = mnx$score)}else{covv_min<-coverage(mnx)}
                if(length(mny)>0){covv_min<-covv_min+coverage(mny,weight = mny$score)}
                
                covv_min<-GRanges(covv_min)
                covv_min<-covv_min[covv_min$score>0]
                
                strand(covv_pl)<-"+"
                strand(covv_min)<-"-"
                
                all_ps[[rl]]<-sort(c(covv_pl,covv_min))
                
            }
            
            uniq_ps<-GRangesList()
            rls<-unique(c(names(x[["P_sites_uniq"]]),names(y[["P_sites_uniq"]])))
            
            for(rl in rls){
                reads_x<-GRanges()
                reads_y<-GRanges()
                
                seqlevels(reads_x)<-seqllll
                seqlevels(reads_y)<-seqllll
                
                seqlengths(reads_x)<-seqleee
                seqlengths(reads_y)<-seqleee
                if(sum(rl%in%names(x[["P_sites_uniq"]]))>0){reads_x<-x[["P_sites_uniq"]][[rl]]}
                if(sum(rl%in%names(y[["P_sites_uniq"]]))>0){reads_y<-y[["P_sites_uniq"]][[rl]]}
                
                plx<-reads_x[strand(reads_x)=="+"]
                mnx<-reads_x[strand(reads_x)=="-"]
                ply<-reads_y[strand(reads_y)=="+"]
                mny<-reads_y[strand(reads_y)=="-"]
                if(length(plx)>0){covv_pl<-coverage(plx,weight = plx$score)}else{covv_pl<-coverage(plx)}
                if(length(ply)>0){covv_pl<-covv_pl+coverage(ply,weight = ply$score)}
                
                covv_pl<-GRanges(covv_pl)
                covv_pl<-covv_pl[covv_pl$score>0]
                
                if(length(mnx)>0){covv_min<-coverage(mnx,weight = mnx$score)}else{covv_min<-coverage(mnx)}
                if(length(mny)>0){covv_min<-covv_min+coverage(mny,weight = mny$score)}
                
                covv_min<-GRanges(covv_min)
                covv_min<-covv_min[covv_min$score>0]
                
                strand(covv_pl)<-"+"
                strand(covv_min)<-"-"
                
                uniq_ps[[rl]]<-sort(c(covv_pl,covv_min))
                
                
            }
            
            uniq_mm_ps<-GRangesList()
            rls<-unique(c(names(x[["P_sites_uniq_mm"]]),names(y[["P_sites_uniq_mm"]])))
            
            for(rl in rls){
                reads_x<-GRanges()
                reads_y<-GRanges()
                
                seqlevels(reads_x)<-seqllll
                seqlevels(reads_y)<-seqllll
                
                seqlengths(reads_x)<-seqleee
                seqlengths(reads_y)<-seqleee
                if(sum(rl%in%names(x[["P_sites_uniq_mm"]]))>0){reads_x<-x[["P_sites_uniq_mm"]][[rl]]}
                if(sum(rl%in%names(y[["P_sites_uniq_mm"]]))>0){reads_y<-y[["P_sites_uniq_mm"]][[rl]]}
                
                plx<-reads_x[strand(reads_x)=="+"]
                mnx<-reads_x[strand(reads_x)=="-"]
                ply<-reads_y[strand(reads_y)=="+"]
                mny<-reads_y[strand(reads_y)=="-"]
                if(length(plx)>0){covv_pl<-coverage(plx,weight = plx$score)}else{covv_pl<-coverage(plx)}
                if(length(ply)>0){covv_pl<-covv_pl+coverage(ply,weight = ply$score)}
                
                covv_pl<-GRanges(covv_pl)
                covv_pl<-covv_pl[covv_pl$score>0]
                
                if(length(mnx)>0){covv_min<-coverage(mnx,weight = mnx$score)}else{covv_min<-coverage(mnx)}
                if(length(mny)>0){covv_min<-covv_min+coverage(mny,weight = mny$score)}
                
                covv_min<-GRanges(covv_min)
                covv_min<-covv_min[covv_min$score>0]
                
                strand(covv_pl)<-"+"
                strand(covv_min)<-"-"
                
                uniq_mm_ps[[rl]]<-sort(c(covv_pl,covv_min))
                
            }
            
            covall_plus<-x$coverage_all_plus+y$coverage_all_plus
            covall_min<-x$coverage_all_min+y$coverage_all_min
            
            covuni_plus<-x$coverage_uniq_plus+y$coverage_uniq_plus
            covuni_min<-x$coverage_uniq_min+y$coverage_uniq_min
            
            
            rang_jun<-x$junctions
            rang_jun$reads<-rang_jun$reads+y$junctions$reads
            rang_jun$unique_reads<-rang_jun$unique_reads+y$junctions$unique_reads
            
            list_res<-list(all_ps,uniq_ps,uniq_mm_ps,rang_jun,covall_plus,covall_min,covuni_plus,covuni_min)
            names(list_res)<-c("P_sites_all","P_sites_uniq","P_sites_uniq_mm","junctions","coverage_all_plus","coverage_all_min","coverage_uniq_plus","coverage_uniq_min")
            
            
            return(list_res)
        }
        
        #what to do with each chunk (read as alignment file)
        
        yiel<-function(x){
            readGAlignments(x,param = param)
        }
        
        #operations on the chunk (here count reads and whatnot)
        
        mapp<-function(x){
            x_I<-x[grep("I",cigar(x))]
            
            if(length(x_I)>0){
                x<-x[grep("I",cigar(x),invert=T)]
                
            }
            x_D<-x[grep("D",cigar(x))]
            if(length(x_D)>0){
                x<-x[grep("D",cigar(x),invert=T)]
                
            }
            
            
            # softclipping
            
            
            clipp <- width(cigarRangesAlongQuerySpace(x@cigar, ops="S"))
            clipp[elementNROWS(clipp)==0] <- 0
            len_adj <- qwidth(x)-sum(clipp)
            mcols(x)$len_adj <- len_adj
            
            # Remove S from Cigar (read positions/length are already adjusted)
            # it helps calculating P-sites positions for spliced reads
            
            cigg<-cigar(x)
            cigg_s<-grep(cigg,pattern = "S")
            if(length(cigg_s)>0){
                cigs<-cigg[cigg_s]
                cigs<-gsub(cigs,pattern = "^[0-9]+S",replacement = "")
                cigs<-gsub(cigs,pattern = "[0-9]+S$",replacement = "")
                cigg[cigg_s]<-cigs
                x@cigar<-cigg
            }
            mcols(x)$cigar_str<-x@cigar
            x_uniq<-x[x@elementMetadata$mapq>50]
            
            pos<-x[strand(x)=="+"]
            neg<-x[strand(x)=="-"]
            
            uniq_pos<-x_uniq[strand(x_uniq)=="+"]
            uniq_neg<-x_uniq[strand(x_uniq)=="-"]
            
            
            # coverage
            
            
            covuni_plus<-coverage(uniq_pos)
            covuni_min<-coverage(uniq_neg)
            covall_plus<-coverage(pos)
            covall_min<-coverage(neg)
            
            
            # junctions
            
            
            juns<-summarizeJunctions(x)
            juns_pos<-juns
            juns_neg<-juns
            mcols(juns_pos)<-NULL
            mcols(juns_neg)<-NULL
            juns_pos$reads<-juns$plus_score
            juns_neg$reads<-juns$minus_score
            strand(juns_pos)<-"+"
            strand(juns_neg)<-"-"
            juns<-sort(c(juns_pos,juns_neg))
            juns<-juns[juns$reads>0]
            
            uniq_juns<-summarizeJunctions(x_uniq)
            uniq_juns_pos<-uniq_juns
            uniq_juns_neg<-uniq_juns
            mcols(uniq_juns_pos)<-NULL
            mcols(uniq_juns_neg)<-NULL
            uniq_juns_pos$reads<-uniq_juns$plus_score
            uniq_juns_neg$reads<-uniq_juns$minus_score
            strand(uniq_juns_pos)<-"+"
            strand(uniq_juns_neg)<-"-"
            uniq_juns<-sort(c(uniq_juns_pos,uniq_juns_neg))
            uniq_juns<-uniq_juns[uniq_juns$reads>0]
            if(length(juns)>0){
                juns$unique_reads<-0
                mat<-match(uniq_juns,juns)
                juns$unique_reads[mat]<-uniq_juns$reads
            }
            rang_jun<-GTF_annotation$junctions
            rang_jun$reads<-0
            rang_jun$unique_reads<-0
            if(length(juns)>0){
                mat<-match(juns,rang_jun)
                juns<-juns[!is.na(mat)]
                mat<-mat[!is.na(mat)]
                rang_jun$reads[mat]<-juns$reads
                rang_jun$unique_reads[mat]<-juns$unique_reads
            }
            
            
            
            # P-sites calculation
            
            
            list_pss<-list()
            for(comp in names(rl_cutoffs_comp)){
                all_rl_ps<-GRangesList()
                uniq_rl_ps<-GRangesList()
                uniq_rl_mm_ps<-GRangesList()
                
                seqlevels(all_rl_ps)<-seqllll
                seqlevels(uniq_rl_ps)<-seqllll
                seqlevels(uniq_rl_mm_ps)<-seqllll
                
                seqlengths(all_rl_ps)<-seqleee
                seqlengths(uniq_rl_ps)<-seqleee
                seqlengths(uniq_rl_mm_ps)<-seqleee
                
                
                chroms<-comp
                
                if(comp=="nucl"){chroms=seqlevels(x)[!seqlevels(x)%in%circs]}
                resul<-rl_cutoffs_comp[[comp]]
                
                for(i in seq_along(resul$read_length)){
                    
                    all_ps<-GRangesList()
                    uniq_ps<-GRangesList()
                    uniq_mm_ps<-GRangesList()
                    seqlevels(all_ps)<-seqllll
                    seqlevels(uniq_ps)<-seqllll
                    seqlevels(uniq_mm_ps)<-seqllll
                    
                    seqlengths(all_ps)<-seqleee
                    seqlengths(uniq_ps)<-seqleee
                    seqlengths(uniq_mm_ps)<-seqleee
                    
                    rl<-as.numeric(resul$read_length[i])
                    ct<-as.numeric(resul$cutoff[i])
                    ok_reads<-pos[mcols(pos)$len_adj%in%rl]
                    ok_reads<-ok_reads[as.vector(seqnames(ok_reads))%in%chroms]
                    
                    ps_plus<-GRanges()
                    seqlevels(ps_plus)<-seqllll
                    seqlengths(ps_plus)<-seqleee
                    ps_plus_uniq<-ps_plus
                    ps_plus_uniq_mm<-ps_plus
                    
                    if(length(ok_reads)>0){
                        unspl<-ok_reads[grep(pattern="N",x=cigar(ok_reads),invert=T)]
                        
                        ps_unspl<-shift(resize(GRanges(unspl),width=1,fix="start"),shift=ct)
                        
                        spl<-ok_reads[grep(pattern="N",x=cigar(ok_reads))]
                        firstb<-as.numeric(sapply(strsplit(cigar(spl),"M"),"[[",1))
                        lastb<-as.numeric(sapply(strsplit(cigar(spl),"M"),function(x){gsub(x[length(x)],pattern="^[^_]*N",replacement="")}))
                        firstok<-spl[firstb>ct]
                        firstok<-shift(resize(GRanges(firstok),width=1,fix="start"),shift=ct)
                        
                        lastok<-spl[lastb>=rl-ct]
                        lastok<-shift(resize(GRanges(lastok),width=1,fix="end"),shift=-(rl-ct-1))
                        
                        
                        multi<-spl[firstb<=ct & lastb<rl-ct]
                        
                        
                        ps_spl<-GRanges()
                        seqlevels(ps_spl)<-seqllll
                        seqlengths(ps_spl)<-seqleee
                        
                        if(length(multi)>0){
                            ps_spl<-get_ps_fromspliceplus(multi,cutoff=ct)
                            
                        }
                        mcols(ps_unspl)<-NULL
                        mcols(firstok)<-NULL
                        mcols(lastok)<-NULL
                        mcols(ps_spl)<-NULL
                        
                        seqlevels(firstok)<-seqllll
                        seqlevels(lastok)<-seqllll
                        seqlevels(ps_unspl)<-seqllll
                        seqlevels(ps_spl)<-seqllll
                        
                        seqlengths(firstok)<-seqleee
                        seqlengths(lastok)<-seqleee
                        seqlengths(ps_unspl)<-seqleee
                        seqlengths(ps_spl)<-seqleee
                        
                        ps_plus<-c(ps_unspl,firstok,lastok,ps_spl)
                        
                        ps_plus_uniq<-ps_plus[mcols(ok_reads)$mapq>50]
                        ps_plus_uniq_mm<-ps_plus[mcols(ok_reads)$mapq>50 & nchar(mcols(ok_reads)$MD)>3]
                        
                    }
                    ok_reads<-neg[mcols(neg)$len_adj%in%rl]
                    ok_reads<-ok_reads[as.vector(seqnames(ok_reads))%in%chroms]
                    
                    ps_neg<-GRanges()
                    seqlevels(ps_neg)<-seqllll
                    seqlengths(ps_neg)<-seqleee
                    ps_neg_uniq<-ps_neg
                    ps_neg_uniq_mm<-ps_neg
                    
                    if(length(ok_reads)>0){
                        unspl<-ok_reads[grep(pattern="N",x=cigar(ok_reads),invert=T)]
                        
                        ps_unspl<-shift(resize(GRanges(unspl),width=1,fix="start"),shift=-ct)
                        
                        spl<-ok_reads[grep(pattern="N",x=cigar(ok_reads))]
                        
                        firstb<-as.numeric(sapply(strsplit(cigar(spl),"M"),"[[",1))
                        lastb<-as.numeric(sapply(strsplit(cigar(spl),"M"),function(x){gsub(x[length(x)],pattern="^[^_]*N",replacement="")}))
                        lastok<-spl[lastb>ct]
                        lastok<-shift(resize(GRanges(lastok),width=1,fix="start"),shift=-ct)
                        
                        firstok<-spl[firstb>=rl-ct]
                        firstok<-shift(resize(GRanges(firstok),width=1,fix="end"),shift=(rl-ct-1))
                        
                        multi<-spl[firstb<rl-ct & lastb<=ct]
                        
                        
                        ps_spl<-GRanges()
                        seqlevels(ps_spl)<-seqllll
                        seqlengths(ps_spl)<-seqleee
                        
                        
                        if(length(multi)>0){
                            ps_spl<-get_ps_fromsplicemin(multi,cutoff=ct)
                        }
                        mcols(ps_unspl)<-NULL
                        mcols(firstok)<-NULL
                        mcols(lastok)<-NULL
                        mcols(ps_spl)<-NULL
                        
                        seqlevels(firstok)<-seqllll
                        seqlevels(lastok)<-seqllll
                        seqlevels(ps_unspl)<-seqllll
                        seqlevels(ps_spl)<-seqllll
                        
                        seqlengths(firstok)<-seqleee
                        seqlengths(lastok)<-seqleee
                        seqlengths(ps_unspl)<-seqleee
                        seqlengths(ps_spl)<-seqleee
                        
                        ps_neg<-c(ps_unspl,firstok,lastok,ps_spl)
                        ps_neg_uniq<-ps_neg[mcols(ok_reads)$mapq>50]
                        
                        ps_neg_uniq_mm<-ps_neg[mcols(ok_reads)$mapq>50 & nchar(mcols(ok_reads)$MD)>3]
                        
                    }
                    
                    all_ps<-sort(c(ps_plus,ps_neg))
                    uniq_ps<-sort(c(ps_plus_uniq,ps_neg_uniq))
                    uniq_mm_ps<-sort(c(ps_plus_uniq_mm,ps_neg_uniq_mm))
                    if(length(all_ps)>0){
                        ps_res<-unique(all_ps)
                        ps_res$score<-countOverlaps(ps_res,all_ps,type="equal")
                        all_ps<-ps_res
                    }
                    if(length(uniq_ps)>0){
                        ps_res<-unique(uniq_ps)
                        ps_res$score<-countOverlaps(ps_res,uniq_ps,type="equal")
                        uniq_ps<-ps_res
                        
                    }
                    if(length(uniq_mm_ps)>0){
                        ps_res<-unique(uniq_mm_ps)
                        ps_res$score<-countOverlaps(ps_res,uniq_mm_ps,type="equal")
                        uniq_mm_ps<-ps_res
                        
                    }
                    all_rl_ps[[as.character(rl)]]<-all_ps
                    uniq_rl_ps[[as.character(rl)]]<-uniq_ps
                    uniq_rl_mm_ps[[as.character(rl)]]<-uniq_mm_ps
                    
                    
                }
                #here comps
                list_rlct<-list(all_rl_ps,uniq_rl_ps,uniq_rl_mm_ps)
                names(list_rlct)<-c("P_sites_all","P_sites_uniq","P_sites_uniq_mm")
                list_pss[[comp]]<-list_rlct
            }
            #for rl, merge psites
            
            
            all_ps_comps<-GRangesList()
            seqlevels(all_ps_comps)<-seqllll
            seqlengths(all_ps_comps)<-seqleee
            rls_comps<-unique(unlist(lapply(list_pss,FUN=function(x) names(x[["P_sites_all"]]) )))
            for(rl in rls_comps){
                reads_rl_comp<-GRanges()
                seqlevels(reads_rl_comp)<-seqllll
                seqlengths(reads_rl_comp)<-seqleee
                for(comp in names(list_pss)){
                    if(sum(rl%in%names(list_pss[[comp]][["P_sites_all"]]))>0){
                        oth<-list_pss[[comp]][["P_sites_all"]][[rl]]
                        
                        if(!is.null(oth)){
                            seqlevels(oth)<-seqllll
                            seqlengths(oth)<-seqleee
                            reads_rl_comp<-c(reads_rl_comp,oth)
                        }
                    }          
                    all_ps_comps[[rl]]<-reads_rl_comp
                }
                
            }
            
            uniq_ps_comps<-GRangesList()
            seqlevels(uniq_ps_comps)<-seqllll
            seqlengths(uniq_ps_comps)<-seqleee
            rls_comps<-unique(unlist(lapply(list_pss,FUN=function(x) names(x[["P_sites_uniq"]]) )))
            for(rl in rls_comps){
                reads_rl_comp<-GRanges()
                seqlevels(reads_rl_comp)<-seqllll
                seqlengths(reads_rl_comp)<-seqleee
                for(comp in names(list_pss)){
                    if(sum(rl%in%names(list_pss[[comp]][["P_sites_uniq"]]))>0){
                        oth<-list_pss[[comp]][["P_sites_uniq"]][[rl]]
                        
                        if(!is.null(oth)){
                            seqlevels(oth)<-seqllll
                            seqlengths(oth)<-seqleee
                            reads_rl_comp<-c(reads_rl_comp,oth)
                        }
                    }          
                    uniq_ps_comps[[rl]]<-reads_rl_comp
                }
                
            }
            
            uniq_mm_ps_comps<-GRangesList()
            seqlevels(uniq_mm_ps_comps)<-seqllll
            seqlengths(uniq_mm_ps_comps)<-seqleee
            rls_comps<-unique(unlist(lapply(list_pss,FUN=function(x) names(x[["P_sites_uniq_mm"]]) )))
            for(rl in rls_comps){
                reads_rl_comp<-GRanges()
                seqlevels(reads_rl_comp)<-seqllll
                seqlengths(reads_rl_comp)<-seqleee
                for(comp in names(list_pss)){
                    if(sum(rl%in%names(list_pss[[comp]][["P_sites_uniq_mm"]]))>0){
                        oth<-list_pss[[comp]][["P_sites_uniq_mm"]][[rl]]
                        if(!is.null(oth)){
                            seqlevels(oth)<-seqllll
                            seqlengths(oth)<-seqleee
                            reads_rl_comp<-c(reads_rl_comp,oth)
                        }
                    }          
                    uniq_mm_ps_comps[[rl]]<-reads_rl_comp
                }
                
            }
            
            
            list_res<-list(all_ps_comps,uniq_ps_comps,uniq_mm_ps_comps,rang_jun,covall_plus,covall_min,covuni_plus,covuni_min)
            names(list_res)<-c("P_sites_all","P_sites_uniq","P_sites_uniq_mm","junctions","coverage_all_plus","coverage_all_minus","coverage_uniq_plus","coverage_uniq_minus")
            
            return(list_res)
            
            
        }
        cat(paste("Calculating P-sites positions and junctions ...", date(),"\n"))
        
        res_4<-reduceByYield(X=opts,YIELD=yiel,MAP=mapp,REDUCE=reduc)
        save(res_4,file = paste(dira,"res_4",sep = "/"))
        
        cat(paste("Calculating P-sites positions and junctions --- Done!", date(),"\n"))
        
        ps<-res_4$P_sites_all
        
        cat(paste("Building aggregate P-sites profiles ...", date(),"\n"))
        
        #here do meta and other plots on P-sites positions ()
        ps_signals_win_all<-list()
        ps_signals_tiles_all<-list()
        codons_win_all<-list()
        ps_codons_ratio<-list()
        ps_codons_counts<-list()
        
        as_codons_ratio<-list()
        as_codons_counts<-list()
        
        es_codons_ratio<-list()
        es_codons_counts<-list()
        
        if(length(ps)==0){
            print("Not enough signal or low frame preference, skipped P-sites & codon occupancy calculation. \n Set 'all' in the 'choose_readlengths' choice to ignore this warning and get P-sites positions in a new run")
        }
        if(length(ps)>0){
            
            for(comp in c("nucl",circs)){
                tile_cds<-list_txs_ok[[comp]]
                #put positive?
                strand(tile_cds)<-"+"
                tile_cds<-tile_cds[width(tile_cds)>100]
                checkgen<-TRUE
                if(comp=="nucl"){
                    if((sum((start(tile_cds)>50))/length(start(tile_cds)))>.33){tile_cds<-tile_cds[start(tile_cds)>50];checkgen<-FALSE}
                    if(sum(seqlengths(tile_cds)[as.vector(seqnames(tile_cds))]-end(tile_cds)>100)/length(tile_cds)>.33){
                        tile_cds<-tile_cds[seqlengths(tile_cds)[as.vector(seqnames(tile_cds))]-end(tile_cds)>100];checkgen<-FALSE}
                    ps_comp<-GRangesList(lapply(ps,function(x){sort(x[!seqnames(x)%in%circs])}))
                } else {ps_comp<-GRangesList(lapply(ps,function(x){sort(x[seqnames(x)==comp])}))}
                
                
                signal_ps<-list()
                signal_ps_nt<-list()
                
                all_ps<-unlist(ps_comp)
                ps_comp[["all"]]<-all_ps
                mp<-mapToTranscripts(all_ps,transcripts = GTF_annotation$exons_txs[as.vector(seqnames(tile_cds))])
                mp$score<-all_ps$score[mp$xHits]
                #put seqlevels for compartment only to speed up stuff
                seqlevels(mp)<-names(GTF_annotation$exons_txs)
                seqlengths(mp)<-sum(width(GTF_annotation$exons_txs))
                covtx<-coverage(mp,weight = mp$score)
                ok_txs<-names(covtx[elementNROWS(runValue(covtx))>1])
                if(length(ok_txs)<5){
                    ps_signals_tiles_all[[comp]]<-c()
                    ps_signals_win_all[[comp]]<-c()
                    next
                }
                #keep all txs if small n of txs
                
                if(length(ok_txs)>=5 & length(covtx)<250){
                    ok_txs<-names(covtx)
                }
                
                if(length(ok_txs)>5000){
                    cnts_txss_agg<-sort(sum(covtx),decreasing = T)
                    ok_txs<-names(cnts_txss_agg)[1:5000]
                }
                
                if(fast_mode[bammo]==T){
                    if(length(ok_txs)>500){
                        cnts_txss_agg<-sort(sum(covtx),decreasing = T)
                        ok_txs<-names(cnts_txss_agg)[1:500]
                    }
                }
                
                
                fivs<-GRanges(seqnames = seqnames(GTF_annotation$cds_txs_coords),ranges = IRanges(start=1,end = start(GTF_annotation$cds_txs_coords)),strand="*")
                fivs<-fivs[as.character(seqnames(fivs))%in%ok_txs]
                fivs_gen<-unlist(pmapFromTranscripts(fivs,transcripts = GTF_annotation$exons_txs[as.character(seqnames(fivs))],ignore.strand=F))
                fivs_gen<-fivs_gen[fivs_gen$hit]
                fivs_gen<-split(fivs_gen,names(fivs_gen))
                
                
                threes<-GRanges(seqnames = seqnames(GTF_annotation$cds_txs_coords),ranges = IRanges(start=end(GTF_annotation$cds_txs_coords),end = GTF_annotation$cds_txs_coords$lentx),strand="*")
                threes<-threes[as.character(seqnames(threes))%in%ok_txs]
                threes_gen<-unlist(pmapFromTranscripts(threes,transcripts = GTF_annotation$exons_txs[as.character(seqnames(threes))],ignore.strand=F))
                threes_gen<-threes_gen[threes_gen$hit]
                threes_gen<-split(threes_gen,names(threes_gen))
                
                cds_gen<-GTF_annotation$cds_txs[ok_txs]
                
                
                tile_cds<-tile_cds[as.character(seqnames(tile_cds))%in%ok_txs]
                ok_txs<-unique(as.character(seqnames(tile_cds)))
                seqlevels(tile_cds)<-ok_txs
                list_covs<-list()
                list_covs[["all"]]<-covtx
                
                tile_5<-GRanges(seqnames(tile_cds),ranges = IRanges(start=1,end = start(tile_cds)))
                tile_3<-GRanges(seqnames(tile_cds),ranges = IRanges(start = end(tile_cds),end = seqlengths(tile_cds)[as.vector(seqnames(tile_cds))]))
                no5utr<-width(tile_5)<51
                no3utr<-width(tile_3)<51
                
                ex_annot<-GTF_annotation$exons_txs
                
                if(sum(no5utr)>0){
                    tx_notok<-seqnames(tile_5)[no5utr]
                    annot_notok<-ex_annot[tx_notok]
                    annot_ok<-GRangesList(lapply(annot_notok, function(x){
                        if(length(x)==0){return(x)}
                        x[1]<-resize(x[1],width = width(x[1])+51,fix = "end")
                        x
                    }))
                    ex_annot[names(annot_ok)]<-annot_ok
                    
                    seqlengths(tile_cds)[as.vector(tx_notok)]<-sum(width(annot_ok))
                    tile_cds[no5utr]<-shift(tile_cds[no5utr],+51)
                    tile_5<-GRanges(seqnames(tile_cds),ranges = IRanges(start=1,end = start(tile_cds)))
                    tile_3<-GRanges(seqnames(tile_cds),ranges = IRanges(start = end(tile_cds),end = seqlengths(tile_cds)[as.vector(seqnames(tile_cds))]))
                    
                    mp<-mapToTranscripts(all_ps,transcripts = ex_annot[as.vector(seqnames(tile_cds))])
                    seqlevels(mp)<-seqlevels(tile_cds)
                    seqlengths(mp)<-seqlengths(tile_cds)
                    mp$score<-all_ps$score[mp$xHits]
                    
                    covtx<-coverage(mp,weight = mp$score)
                    list_covs[["all"]]<-covtx
                }
                
                if(sum(no3utr)>0){
                    tx_notok<-seqnames(tile_3)[no3utr]
                    annot_notok<-ex_annot[tx_notok]
                    annot_ok<-GRangesList(lapply(annot_notok, function(x){
                        if(length(x)==0){return(x)}
                        x[length(x)]<-resize(x[length(x)],width = width(x[length(x)])+51,fix = "start")
                        x
                    }))
                    ex_annot[names(annot_ok)]<-annot_ok
                    
                    seqlengths(tile_cds)[as.vector(tx_notok)]<-sum(width(annot_ok))
                    
                    tile_5<-GRanges(seqnames(tile_cds),ranges = IRanges(start=1,end = start(tile_cds)))
                    tile_3<-GRanges(seqnames(tile_cds),ranges = IRanges(start = end(tile_cds),end = seqlengths(tile_cds)[as.vector(seqnames(tile_cds))]))
                    
                    mp<-mapToTranscripts(all_ps,transcripts = ex_annot[as.vector(seqnames(tile_cds))])
                    seqlevels(mp)<-seqlevels(tile_cds)
                    seqlengths(mp)<-seqlengths(tile_cds)
                    mp$score<-all_ps$score[mp$xHits]
                    
                    covtx<-coverage(mp,weight = mp$score)
                    list_covs[["all"]]<-covtx
                    
                }
                
                
                ps_tiles<-DataFrameList()
                ps_win<-DataFrameList()
                cod_win<-DataFrameList()
                ps_cod<-DataFrameList()
                ps_cod_rat<-DataFrameList()
                
                as_cod<-DataFrameList()
                as_cod_rat<-DataFrameList()
                es_cod<-DataFrameList()
                es_cod_rat<-DataFrameList()
                
                for(len in c("all",names(ps_comp))){
                    
                    if((sum(no3utr) + sum(no5utr))==0){
                        mp<-mapToTranscripts(ps_comp[[len]],transcripts = fivs_gen)
                        mp$score<-ps_comp[[len]]$score[mp$xHits]
                        seqlevels(mp)<-names(fivs_gen)
                        seqlengths(mp)<-sum(width(fivs_gen))
                        cov_5<-coverage(mp,weight = mp$score)
                        
                        
                        mp<-mapToTranscripts(ps_comp[[len]],transcripts = threes_gen)
                        mp$score<-ps_comp[[len]]$score[mp$xHits]
                        seqlevels(mp)<-names(threes_gen)
                        seqlengths(mp)<-sum(width(threes_gen))
                        cov_3<-coverage(mp,weight = mp$score)
                        
                        mp<-mapToTranscripts(ps_comp[[len]],transcripts = cds_gen)
                        mp$score<-ps_comp[[len]]$score[mp$xHits]
                        seqlevels(mp)<-names(cds_gen)
                        seqlengths(mp)<-sum(width(cds_gen))
                        strand(mp)<-"+"
                        cov_cds<-coverage(mp,weight = mp$score)
                        cov_cds_a<-suppressWarnings(coverage(shift(mp,shift = 3),weight = mp$score))
                        cov_cds_e<-suppressWarnings(coverage(shift(mp,shift = -3),weight = mp$score))
                    }
                    
                    if((sum(no3utr) + sum(no5utr))>0){
                        
                        mp<-mapToTranscripts(ps_comp[[len]],transcripts = ex_annot[as.vector(seqnames(tile_cds))])
                        mp$score<-ps_comp[[len]]$score[mp$xHits]
                        seqlevels(mp)<-seqlevels(tile_cds)
                        seqlengths(mp)<-seqlengths(tile_cds)
                        strand(mp)<-"+"
                        covtx<-coverage(mp,weight = mp$score)
                        covtx<-covtx[ok_txs]
                        cov_5<-covtx[tile_5]
                        cov_3<-covtx[tile_3]
                        cov_cds<-covtx[tile_cds]
                        cov_cds_a<-suppressWarnings(coverage(shift(mp,shift = 3),weight = mp$score))[tile_cds]
                        cov_cds_e<-suppressWarnings(coverage(shift(mp,shift = -3),weight = mp$score))[tile_cds]
                        
                    }
                    
                    ps_tiles_5<-DataFrame(t(sapply(cov_5,function(x){
                        clos<-50*(round(length(x)/50,digits = 0)+1)
                        idx<-as.integer(seq(1,length(x),length.out = clos))
                        colMeans(matrix(x[idx],ncol = 50))
                    })))
                    ps_win_5<-DataFrame(t(sapply(cov_5,function(x){as.vector(x)[c(1:25,(length(x)-25):(length(x)-1))]})))
                    
                    
                    ps_tiles_3<-DataFrame(t(sapply(cov_3,function(x){
                        clos<-50*(round(length(x)/50,digits = 0)+1)
                        idx<-as.integer(seq(1,length(x),length.out = clos))
                        colMeans(matrix(x[idx],ncol = 50))
                    })))
                    ps_win_3<-DataFrame(t(sapply(cov_3,function(x){as.vector(x)[c(2:26,(length(x)-24):length(x))]})))
                    
                    
                    ps_tiles_cds<-DataFrame(t(sapply(cov_cds,function(x){
                        clos<-100*(round(length(x)/100,digits = 0)+1)
                        idx<-as.integer(seq(1,length(x),length.out = clos))
                        colMeans(matrix(x[idx],ncol = 100))
                    })))
                    ps_win_cds<-DataFrame(t(sapply(cov_cds,function(x){
                        rnd<-as.integer(length(x)/2)%%3
                        mid<-(as.integer(length(x)/2)-rnd)
                        as.vector(x)[c(1:33,(mid-17):(mid+15),(length(x)-32):length(x))]
                    })))
                    
                    as_win_cds<-DataFrame(t(sapply(cov_cds_a,function(x){
                        rnd<-as.integer(length(x)/2)%%3
                        mid<-(as.integer(length(x)/2)-rnd)
                        as.vector(x)[c(1:33,(mid-17):(mid+15),(length(x)-32):length(x))]
                    })))
                    
                    es_win_cds<-DataFrame(t(sapply(cov_cds_e,function(x){
                        rnd<-as.integer(length(x)/2)%%3
                        mid<-(as.integer(length(x)/2)-rnd)
                        as.vector(x)[c(1:33,(mid-17):(mid+15),(length(x)-32):length(x))]
                    })))
                    
                    #select txs, output codon usage
                    txs_seqq<-extractTranscriptSeqs(x = genome_seq,transcripts = GTF_annotation$cds_txs[names(cov_cds)])
                    
                    gco <- as.character(names(getGeneticCode("1")))
                    names(gco)<-as.character(getGeneticCode("1"))
                    
                    #Get genetic code for compartment
                    
                    gen_cods<-GTF_annotation$genetic_codes
                    circ_cods<-which(comp==rownames(gen_cods))
                    if(length(circ_cods)>0){
                        
                        
                        gco <- as.character(names(getGeneticCode(gen_cods$genetic_code[circ_cods])))
                        names(gco)<-as.character(getGeneticCode(gen_cods$genetic_code[circ_cods]))
                        
                    }
                    
                    rnd<-as.integer(width(txs_seqq)/2)%%3
                    st_pos<-rep(1,length(txs_seqq))
                    mid_pos<-(as.integer(width(txs_seqq)/2)-rnd)-17
                    end_pos<-width(txs_seqq)-32
                    
                    iters<-seq(0,10)
                    cod_counts<-c()
                    psit_counts<-c()
                    asit_counts<-c()
                    esit_counts<-c()
                    #Calculate codon occurrence and occupancy for different CDS sections
                    #OUTPUT AA?
                    for(i in iters){
                        a<-narrow(txs_seqq,start = st_pos+3*i,end = st_pos+(2+3*i))
                        coood<-as.character(a)
                        at<-table(a)
                        cod_cntt<-rep(0,length(gco))
                        names(cod_cntt)<-gco
                        cod_cntt[names(at)]<-as.numeric(at)
                        cod_counts<-cbind(cod_counts,cod_cntt)
                        pt<-rowSums(as.matrix(ps_win_cds[,(1+3*i):(3+3*i)]))
                        names(pt)<-coood
                        pt<-aggregate(pt,list(names(pt)),sum)
                        ps_cntt<-rep(0,length(gco))
                        names(ps_cntt)<-gco
                        ps_cntt[pt[,1]]<-pt[,2]
                        psit_counts<-cbind(psit_counts,ps_cntt)
                        
                        
                        pt<-rowSums(as.matrix(as_win_cds[,(1+3*i):(3+3*i)]))
                        names(pt)<-coood
                        pt<-aggregate(pt,list(names(pt)),sum)
                        ps_cntt<-rep(0,length(gco))
                        names(ps_cntt)<-gco
                        ps_cntt[pt[,1]]<-pt[,2]
                        asit_counts<-cbind(asit_counts,ps_cntt)
                        
                        
                        pt<-rowSums(as.matrix(es_win_cds[,(1+3*i):(3+3*i)]))
                        names(pt)<-coood
                        pt<-aggregate(pt,list(names(pt)),sum)
                        ps_cntt<-rep(0,length(gco))
                        names(ps_cntt)<-gco
                        ps_cntt[pt[,1]]<-pt[,2]
                        esit_counts<-cbind(esit_counts,ps_cntt)
                        
                    }
                    #mid
                    for(i in iters){
                        a<-narrow(txs_seqq,start = mid_pos+3*i,end = mid_pos+(2+3*i))
                        coood<-as.character(a)
                        at<-table(a)
                        cod_cntt<-rep(0,length(gco))
                        names(cod_cntt)<-gco
                        cod_cntt[names(at)]<-as.numeric(at)
                        cod_counts<-cbind(cod_counts,cod_cntt)
                        pt<-rowSums(as.matrix(ps_win_cds[,(34+3*i):(36+3*i)]))
                        names(pt)<-coood
                        pt<-aggregate(pt,list(names(pt)),sum)
                        ps_cntt<-rep(0,length(gco))
                        names(ps_cntt)<-gco
                        ps_cntt[pt[,1]]<-pt[,2]
                        psit_counts<-cbind(psit_counts,ps_cntt)
                        
                        pt<-rowSums(as.matrix(as_win_cds[,(34+3*i):(36+3*i)]))
                        names(pt)<-coood
                        pt<-aggregate(pt,list(names(pt)),sum)
                        ps_cntt<-rep(0,length(gco))
                        names(ps_cntt)<-gco
                        ps_cntt[pt[,1]]<-pt[,2]
                        asit_counts<-cbind(asit_counts,ps_cntt)
                        
                        pt<-rowSums(as.matrix(es_win_cds[,(34+3*i):(36+3*i)]))
                        names(pt)<-coood
                        pt<-aggregate(pt,list(names(pt)),sum)
                        ps_cntt<-rep(0,length(gco))
                        names(ps_cntt)<-gco
                        ps_cntt[pt[,1]]<-pt[,2]
                        esit_counts<-cbind(esit_counts,ps_cntt)
                    }
                    #end
                    for(i in iters){
                        a<-narrow(txs_seqq,start = end_pos+3*i,end = end_pos+(2+3*i))
                        coood<-as.character(a)
                        at<-table(a)
                        cod_cntt<-rep(0,length(gco))
                        names(cod_cntt)<-gco
                        cod_cntt[names(at)]<-as.numeric(at)
                        cod_counts<-cbind(cod_counts,cod_cntt)
                        pt<-rowSums(as.matrix(ps_win_cds[,(67+3*i):(69+3*i)]))
                        names(pt)<-coood
                        pt<-aggregate(pt,list(names(pt)),sum)
                        ps_cntt<-rep(0,length(gco))
                        names(ps_cntt)<-gco
                        ps_cntt[pt[,1]]<-pt[,2]
                        psit_counts<-cbind(psit_counts,ps_cntt)
                        
                        pt<-rowSums(as.matrix(as_win_cds[,(67+3*i):(69+3*i)]))
                        names(pt)<-coood
                        pt<-aggregate(pt,list(names(pt)),sum)
                        ps_cntt<-rep(0,length(gco))
                        names(ps_cntt)<-gco
                        ps_cntt[pt[,1]]<-pt[,2]
                        asit_counts<-cbind(asit_counts,ps_cntt)
                        
                        pt<-rowSums(as.matrix(es_win_cds[,(67+3*i):(69+3*i)]))
                        names(pt)<-coood
                        pt<-aggregate(pt,list(names(pt)),sum)
                        ps_cntt<-rep(0,length(gco))
                        names(ps_cntt)<-gco
                        ps_cntt[pt[,1]]<-pt[,2]
                        esit_counts<-cbind(esit_counts,ps_cntt)
                        
                    }
                    
                    colnames(psit_counts)<-c(paste("first",seq(1,11),sep="_"),paste("mid",seq(1,11),sep="_"),paste("last",seq(1,11),sep="_"))
                    colnames(esit_counts)<-c(paste("first",seq(1,11),sep="_"),paste("mid",seq(1,11),sep="_"),paste("last",seq(1,11),sep="_"))
                    colnames(asit_counts)<-c(paste("first",seq(1,11),sep="_"),paste("mid",seq(1,11),sep="_"),paste("last",seq(1,11),sep="_"))
                    colnames(cod_counts)<-c(paste("first",seq(1,11),sep="_"),paste("mid",seq(1,11),sep="_"),paste("last",seq(1,11),sep="_"))
                    ratio_psit_cod<-psit_counts/cod_counts
                    ratio_asit_cod<-asit_counts/cod_counts
                    ratio_esit_cod<-esit_counts/cod_counts
                    
                    psit_counts<-DataFrame(psit_counts)
                    asit_counts<-DataFrame(asit_counts)
                    esit_counts<-DataFrame(esit_counts)
                    cod_counts<-DataFrame(cod_counts)
                    ratio_psit_cod<-DataFrame(ratio_psit_cod)
                    ratio_asit_cod<-DataFrame(ratio_asit_cod)
                    ratio_esit_cod<-DataFrame(ratio_esit_cod)
                    
                    gcall<-paste(gco,names(gco),sep=";")
                    
                    rownames(psit_counts)<-gcall[match(rownames(psit_counts),gco)]
                    rownames(asit_counts)<-gcall[match(rownames(asit_counts),gco)]
                    rownames(esit_counts)<-gcall[match(rownames(esit_counts),gco)]
                    
                    rownames(cod_counts)<-gcall[match(rownames(cod_counts),gco)]
                    rownames(ratio_psit_cod)<-gcall[match(rownames(ratio_psit_cod),gco)]
                    rownames(ratio_asit_cod)<-gcall[match(rownames(ratio_asit_cod),gco)]
                    rownames(ratio_esit_cod)<-gcall[match(rownames(ratio_esit_cod),gco)]
                    
                    
                    cod_win[[len]]<-cod_counts
                    ps_cod[[len]]<-psit_counts
                    ps_cod_rat[[len]]<-ratio_psit_cod
                    
                    as_cod[[len]]<-asit_counts
                    as_cod_rat[[len]]<-ratio_asit_cod
                    
                    
                    es_cod[[len]]<-esit_counts
                    es_cod_rat[[len]]<-ratio_esit_cod
                    
                    
                    tls<-cbind(ps_tiles_5,ps_tiles_cds,ps_tiles_3)
                    colnames(tls)<-c(paste("5_UTR",1:length(ps_tiles_5[1,]),sep="_"),paste("CDS",1:length(ps_tiles_cds[1,]),sep="_"),paste("3_UTR",1:length(ps_tiles_3[1,]),sep="_"))
                    nts<-cbind(ps_win_5,ps_win_cds,ps_win_3)
                    colnames(nts)<-c(paste("5_UTR",1:length(ps_win_5[1,]),sep="_"),paste("CDS",1:length(ps_win_cds[1,]),sep="_"),paste("3_UTR",1:length(ps_win_3[1,]),sep="_"))
                    
                    
                    ps_tiles[[len]]<-tls
                    ps_win[[len]]<-nts
                    
                }
                
                codons_win_all[[comp]]<-cod_win
                ps_codons_ratio[[comp]]<-ps_cod_rat
                ps_codons_counts[[comp]]<-ps_cod
                
                as_codons_ratio[[comp]]<-as_cod_rat
                as_codons_counts[[comp]]<-as_cod
                
                
                es_codons_ratio[[comp]]<-es_cod_rat
                es_codons_counts[[comp]]<-es_cod
                
                ps_signals_tiles_all[[comp]]<-ps_tiles
                ps_signals_win_all[[comp]]<-ps_win
                
            }
        }
        
        res_5<-list(ps_signals_tiles_all,ps_signals_win_all,codons_win_all,ps_codons_counts,ps_codons_ratio,as_codons_counts,as_codons_ratio,es_codons_counts,es_codons_ratio)
        
        names(res_5)<-c("P_sites_bins","P_sites_subcodon","Codon_counts","P_sites_percodon","P_sites_percodon_ratio",
                        "A_sites_percodon","A_sites_percodon_ratio","E_sites_percodon","E_sites_percodon_ratio")
        save(res_5,file = paste(dira,"res_5",sep = "/"))
        #save(list = c("ps_signals_tiles","ps_signals_tiles_all","ps_signals_win","ps_signals_win_all"),file = "ps_results_preprjan15")
        cat(paste("Building aggregate P-sites profiles --- Done!", date(),"\n"))
        
        cat(paste("Exporting files ...", date(),"\n"))
        
        merged_all_ps<-unlist(res_4$P_sites_all)
        if(length(merged_all_ps)>0){
            covv_pl<-coverage(merged_all_ps[strand(merged_all_ps)=="+"],weight = merged_all_ps[strand(merged_all_ps)=="+"]$score)
            covv_pl<-GRanges(covv_pl)
            covv_pl<-covv_pl[covv_pl$score>0]
            covv_min<-coverage(merged_all_ps[strand(merged_all_ps)=="-"],weight = merged_all_ps[strand(merged_all_ps)=="-"]$score)
            covv_min<-GRanges(covv_min)
            covv_min<-covv_min[covv_min$score>0]
            strand(covv_pl)<-"+"
            strand(covv_min)<-"-"
            
            merged_all_ps<-sort(c(covv_pl,covv_min))
            
            
            export(covv_pl,con=paste(dest_name,"_P_sites_plus.bw",sep = ""))
            export(covv_min,con=paste(dest_name,"_P_sites_minus.bw",sep = ""))
        }
        
        merged_uniq_ps<-unlist(res_4$P_sites_uniq)
        if(length(merged_uniq_ps)>0){
            
            covv_pl<-coverage(merged_uniq_ps[strand(merged_uniq_ps)=="+"],weight = merged_uniq_ps[strand(merged_uniq_ps)=="+"]$score)
            covv_pl<-GRanges(covv_pl)
            covv_pl<-covv_pl[covv_pl$score>0]
            
            covv_min<-coverage(merged_uniq_ps[strand(merged_uniq_ps)=="-"],weight = merged_uniq_ps[strand(merged_uniq_ps)=="-"]$score)
            covv_min<-GRanges(covv_min)
            covv_min<-covv_min[covv_min$score>0]
            strand(covv_pl)<-"+"
            strand(covv_min)<-"-"
            merged_uniq_ps<-sort(c(covv_pl,covv_min))
            
            
            export(covv_pl,con=paste(dest_name,"_P_sites_uniq_plus.bw",sep = ""))
            export(covv_min,con=paste(dest_name,"_P_sites_uniq_minus.bw",sep = ""))
        }
        
        
        merged_uniq_mm_ps<-unlist(res_4$P_sites_uniq_mm)
        if(length(merged_uniq_mm_ps)>0){
            
            covv_pl<-coverage(merged_uniq_mm_ps[strand(merged_uniq_mm_ps)=="+"],weight = merged_uniq_mm_ps[strand(merged_uniq_mm_ps)=="+"]$score)
            covv_pl<-GRanges(covv_pl)
            covv_pl<-covv_pl[covv_pl$score>0]
            
            covv_min<-coverage(merged_uniq_mm_ps[strand(merged_uniq_mm_ps)=="-"],weight = merged_uniq_mm_ps[strand(merged_uniq_mm_ps)=="-"]$score)
            covv_min<-GRanges(covv_min)
            covv_min<-covv_min[covv_min$score>0]
            strand(covv_pl)<-"+"
            strand(covv_min)<-"-"
            merged_uniq_mm_ps<-sort(c(covv_pl,covv_min))
        }
        
        
        export(res_4$coverage_all_plus,con=paste(dest_name,"_coverage_plus.bw",sep = ""))
        export(res_4$coverage_all_min,con=paste(dest_name,"_coverage_minus.bw",sep = ""))
        export(res_4$coverage_uniq_plus,con=paste(dest_name,"_coverage_uniq_plus.bw",sep = ""))
        export(res_4$coverage_uniq_min,con=paste(dest_name,"_coverage_uniq_minus.bw",sep = ""))
        
        juns<-res_4$junctions
        save(juns,file=paste(dest_name,"_junctions",sep=""))
        
        
        for_SaTAnn<-list(merged_all_ps,merged_uniq_ps,merged_uniq_mm_ps,juns)
        names(for_SaTAnn)<-c("P_sites_all","P_sites_uniq","P_sites_uniq_mm","junctions")
        
        
        # Store top 50 occupied regions (possible contaminants)
        
        
        rds_st<-res_1$reads_pos1
        rds_st<-sort(unlist(rds_st))
        
        regions <- list(reduce(unlist(GTF_annotation$cds_genes)), GTF_annotation$fiveutrs, GTF_annotation$threeutrs,
                        GTF_annotation$ncIsof, GTF_annotation$ncRNAs, GTF_annotation$introns, GTF_annotation$intergenicRegions)
        names(regions) <- c("cds","fiveutrs","threeutrs",
                            "ncIsof","ncRNAs","introns","intergenic")
        
        regions<-GRangesList(endoapply(regions,function(x){mcols(x)<-NULL;x}))
        
        rds_st<-rds_st[seqnames(rds_st)%in%unique(seqlevels(regions))]
        rds_st$len_adj<-as.numeric(names(rds_st))
        names(rds_st)<-NULL
        
        rds_st<-rds_st[order(rds_st$score,decreasing = T)]
        rds_st$pct<-round(rds_st$score/(sum(rds_st$score)/100),digits = 4)
        top50<-head(rds_st,50)
        top50<-resize(top50,width = top50$len_adj,fix = "start")
        top50$seq<-getSeq(genome_seq,top50)
        top50$region<-CharacterList("")
        
        ovg<-findOverlaps(top50,regions)
        ovg<-IntegerList(split(subjectHits(ovg),queryHits(ovg)))
        subsr<-as.numeric(names(ovg))
        a_reg<-CharacterList(lapply(ovg,FUN = function(x){unique(names(regions)[x])}))
        top50$region[subsr]<-a_reg
        
        ovg<-findOverlaps(top50,GTF_annotation$genes)
        ovg<-IntegerList(split(subjectHits(ovg),queryHits(ovg)))
        subs<-as.numeric(names(ovg))
        a_nam<-CharacterList(lapply(ovg,FUN = function(x){unique(names(GTF_annotation$genes)[x])}))
        a_biot<-CharacterList(lapply(a_nam,FUN = function(x){GTF_annotation$trann$gene_biotype[match(x,GTF_annotation$trann$gene_id)]}))
        top50$gene<-CharacterList("")
        top50$gene_biotype<-CharacterList("")
        top50$gene[subs]<-a_nam
        top50$gene_biotype[subs]<-a_biot
        res_6<-top50
        
        # save as list_results
        
        save(res_6,file = paste(dira,"res_6",sep = "/"))
        
        
        res_all <- list(res_1,res_2,res_3,res_4,res_5,res_6)
        names(res_all) <- c("read_stats", "profiles_fivepr","selection_cutoffs","P_sites_stats","profiles_P_sites","sequence_analysis")
        if(length(merged_all_ps)>0){
            
            save(for_SaTAnn,file = paste(dest_name,"_for_SaTAnn",sep = ""))
        }
        finn<-lapply(res_all$selection_cutoffs$results_choice,function(x){x$data})
        nope<-length(finn)
        
        namfinn<-rep(names(finn),elementNROWS(finn))
        finn<-do.call(rbind,args = finn)
        finn$comp<-namfinn
        
        pct_lib<-round(t(res_1$rld/sum(res_1$rld)*100),digits = 4)
        rownames(pct_lib)<-gsub(rownames(pct_lib),pattern = "reads_",replacement = "")
        finn$pct_map<-0
        for(i in unique(namfinn)){finn$pct_map[finn$comp==i]<-pct_lib[match(finn$read_length[finn$comp==i],rownames(pct_lib)),i]}
        if(nope==0){finn<-NULL}
        res_all$summary_P_sites<-finn
        if(nope==0){res_all$summary_P_sites<-c("")}
        if(write_tmp_files){
            save(res_all, file=paste(dest_name,"_results_RiboseQC_all",sep = ""))
        }
        write.table(finn, file=paste(dest_name,"_P_sites_calcs",sep = ""),sep = "\t",quote = FALSE,row.names = FALSE,col.names = TRUE)
        unlink(dira,recursive = T)
        res_all$read_stats$reads_pos1<-""
        res_all$P_sites_stats<-""
        
        nms<-names(res_all$profiles_fivepr$five_prime_bins)
        for(i in nms){
            fvp<-res_all$profiles_fivepr$five_prime_bins[[i]]
            res_all$profiles_fivepr$five_prime_bins[[i]]<-lapply(fvp,function(x){colSums(as.matrix(x))})
        }
        nms<-names(res_all$profiles_fivepr$five_prime_subcodon)
        for(i in nms){
            fvp<-res_all$profiles_fivepr$five_prime_subcodon[[i]]
            res_all$profiles_fivepr$five_prime_subcodon[[i]]<-lapply(fvp,function(x){colSums(as.matrix(x))})
        }
        
        nms<-names(res_all$profiles_P_sites$P_sites_bins)
        for(i in nms){
            fvp<-res_all$profiles_P_sites$P_sites_bins[[i]]
            res_all$profiles_P_sites$P_sites_bins[[i]]<-lapply(fvp,function(x){colSums(as.matrix(x))})
        }
        nms<-names(res_all$profiles_P_sites$P_sites_subcodon)
        for(i in nms){
            fvp<-res_all$profiles_P_sites$P_sites_subcodon[[i]]
            res_all$profiles_P_sites$P_sites_subcodon[[i]]<-lapply(fvp,function(x){colSums(as.matrix(x))})
        }
        
        save(res_all, file=paste(dest_name,"_results_RiboseQC",sep = ""))
        
        cat(paste("Exporting files --- Done!", date(),"\n\n"))
    }
    
    if(create_report){
        cat(paste("Creating html report in ",report_file," ... ", date(),"\n",sep=""))
        
        create_html_report(input_files=paste(dest_names,"_results_RiboseQC",sep = ""), input_sample_names=sample_names,output_file =  report_file,extended = extended_report)
        cat(paste("Creating html report --- Done!", date(),"\n\n"))
        if(pdf_plots){create_pdfs_from_rds_objects(output_rds_path = paste(report_file,"_plots/rds/",sep=""))}
        
    }
    
}

