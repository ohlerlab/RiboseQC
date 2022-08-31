# RiboseQC, a comprehensive Ribo-seq quality control tool
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

#' @import GenomicFeatures
#' @import GenomicFiles
NULL

######################################################################################################
# Methods for RiboseQC report
######################################################################################################



#' Create the RiboseQC analysis report in html
#'
#' This function creates the RiboseQC html report based on the RiboseQC analysis files
#' generated with \code{RiboseQC_analysis}.
#'
#' @keywords RiboseQC
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
#' (RiboseQC analysis files generated with \code{RiboseQC_analysis}) and \cr \cr
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
#' @examples
#' \dontrun{
#' sampnames = c("root","shoot")
#' create_html_report(input_files=paste(sampnames,"_results_RiboseQC",sep = ""), input_sample_names=sampnames,output_file =  'test_root_shoots.html',extended = FALSE)
#' }
#' @export

create_html_report <- function(input_files, input_sample_names, output_file,extended=FALSE){

    # get input and output file paths
    input_files <- paste(normalizePath(dirname(input_files)),basename(input_files),sep="/")
    output_file <- paste(normalizePath(dirname(output_file)),basename(output_file),sep="/")

    # get path to RMarkdown file (to be rendered)
    rmd_path <- paste(system.file(package="RiboseQC",mustWork = TRUE),"/rmd/riboseqc_template.Rmd",sep="")
    if(extended){
        rmd_path <- paste(system.file(package="RiboseQC",mustWork = TRUE),"/rmd/riboseqc_template_full.Rmd",sep="")
    }

    # get folder path for pdf figures
    output_fig_path <- paste(output_file,"_plots/", sep = "")

    # create folder for rds ojects and pdf figures
    dir.create(paste0(output_fig_path, "rds/"), recursive=TRUE, showWarnings=FALSE)
    dir.create(paste0(output_fig_path, "pdf/"), recursive=TRUE, showWarnings=FALSE)
    sink(file = paste(output_file,"_report_text_output.txt",sep = ""))
    # render RMarkdown file > html report
    knitclean<-knitr::knit_meta(class=NULL, clean = TRUE)
    suppressWarnings(render(rmd_path,
           params = list(input_files = input_files,
                         input_sample_names = input_sample_names,
                         output_fig_path = output_fig_path),
           output_file = output_file))
    gici<-gc()
    sink()
}


#' Generate PDF files from RDS object files
#'
#' This function generates figures as PDF files from RDS object files -
#' it's essentially a workaround for  difficulty integrating pdfs into
#' html reports.
#'
#' @keywords RiboseQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#'
#' @param output_rds_path String; full path to output folder for RDS object files. Example: /my_path_to/rds/
#'
#' @return This function creates PDF files from RDS object files.
#'
#' @seealso \code{\link{create_html_report}}
#' @examples
#'
#' @examples
#' create_pdfs_from_rds_objects(paste('test_root_shoots.html',"_plots/rds/",sep=""))
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
#' to be used during RiboseQC report generation.
#'
#' @keywords RiboseQC
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
#' @examples
#'
#' @examples
#' rdata_list <- generate_rdata_list("root_results_RiboseQC")
#' @export

generate_rdata_list <- function(input_files){
    rdata_list <- list()
    for (i in seq_len(length(input_files))){
        load(input_files[i])
        rdata_list[[names(input_files)[i]]] <- res_all
    }
    return(rdata_list)
}


#' Plot read location distribution by biotype (and originating compartment)
#'
#' This function plots the read location distribution by biotype (and originating compartment)
#' for one input sample. \cr \cr
#' This plot is used in the RiboseQC report in section 1.1.
#'
#' @keywords RiboseQC
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
#' @examples
#' data(res_all)
#' plot_read_biotype_dist_1(
#'   res_all$read_stats$positions, 'shoot', paste0('myfigpath', "rds/"))
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
#' This plot is used in the RiboseQC report in section 1.2.
#'
#' @keywords RiboseQC
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
#' @examples
#' data(res_all)
#' plot_read_biotype_dist_2(list(res_all), paste0('myfigpath', "rds/"))
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
        scale_fill_viridis(discrete=TRUE) +
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
        scale_fill_viridis(discrete=TRUE) +
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
#' This plot is used in the RiboseQC report in section 2.
#'
#' @keywords RiboseQC
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
#' @examples
#' data(res_all)
#' plot_read_length_dist(res_all$read_stats$rld, 'shoot', paste0('myfigpath', "rds/"))
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
#' This plot is used in the RiboseQC report in section 3.1.
#'
#' @keywords RiboseQC
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
#' @examples
#' data(res_all)
#' plot_read_length_dist_by_biotype(res_all$read_stats$reads_summary, 'shoot', paste0('myfigpath', "rds/"))
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
        scale_color_viridis(discrete=TRUE,alpha = .8) +
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
        scale_color_viridis(discrete=TRUE,alpha = .8) +
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
#' This plot is used in the RiboseQC report in section 3.2.
#'
#' @keywords RiboseQC
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
#' @examples
#' data(res_all)
#' plot_read_biotype_dist_by_length(res_all$read_stats$reads_summary, 'shoot', paste0('myfigpath', "rds/"))
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
        scale_fill_viridis(discrete=TRUE) +
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
        scale_fill_viridis(discrete=TRUE) +
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
#' This data is used as input in \code{plot_metagene_hm} to generate plot for the RiboseQC report in section 4.1/4.3.
#'
#' @keywords RiboseQC
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
    rl_ok <- c(rl_ok[1], as.character(sort(as.integer(rl_ok[-1], decreasing=FALSE))))

    # prepare BARPLOTS data

    data_bar = list()

    for (rl in rl_ok) {
        # get read length data and format
        x <- signal[[rl]]
        x.m <- suppressMessages(melt(x))
        x.m$pos <- seq_len(length(x))
        x.m$frame <- paste0("frame ", rep(c(3,1,2),67)[seq_along(x)])

        # add read length data to list
        data_bar[[rl]] <- x.m
    }

    # prepare HEATMAP data

    # we want to make two separate heatmaps for
    # a) single read lenghts (single) and
    # b) summarized read lengths (all).

    # single read lengths
    rl_ok_single <- rl_ok[-1]
    rl_ok_single <- c(as.character(sort(as.integer(rl_ok_single), decreasing = FALSE)))
    data_none <- do.call(signal[rl_ok_single], what=rbind)
    rownames(data_none) <- rl_ok_single
    colnames(data_none) <- seq_len(dim(data_none)[2])
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
    colnames(data_none) <- seq_len(dim(data_none)[2])
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
#' These plots are displayed in the RiboseQC report in section 4.1/4.3.
#'
#' @keywords RiboseQC
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
#' @examples
#' data(res_all)
#' plot_metagene_hm_rmd(res_all, "profiles_fivepr", 'shoot', paste0('myfigpath', "rds/"))
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
#' This plot is used in the RiboseQC report in section 4.1/4.3.
#'
#' @keywords RiboseQC
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
#' This plot is used in the RiboseQC report in section 4.1/4.3.
#'
#' @keywords RiboseQC
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
#' @examples
#' data(res_all)
#' plot_metagene_hm_rmd(res_all, "profiles_fivepr", 'shoot', paste0('myfigpath', "rds/"))
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
#' This plot is used in the RiboseQC report in section 4.1/4.3.
#'
#' @keywords RiboseQC
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
#' These plots are displayed in the RiboseQC report in section 4.2.1.
#'
#' @keywords RiboseQC
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
#' @examples
#' data(res_all)
#' plot_frame_dist_boxplot_rmd(res_all$selection_cutoffs$analysis_frame_cutoff, 'shoot', paste0('myfigpath', "rds/"))
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
#' This plot is used in the RiboseQC report in section 4.2.1.
#'
#' @keywords RiboseQC
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
#' This data is used in the RiboseQC report in section 4.2.2 (displayed as table).
#'
#' @keywords RiboseQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#'
#' @param rdata_list List of RiboseQC analysis RData objects generated by \code{generate_rdata_list}.
#'
#' @return This function returns data.
#'
#' @seealso \code{\link{create_html_report}}

get_rl_and_cutoffs <- function(rdata_list) {
    datasets <- NULL
    for (i in c(seq_along(rdata_list))) {
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
#' This data is used in the RiboseQC report in section 4.2.3 (displayed as table).
#'
#' @keywords RiboseQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#'
#' @param rdata_list List of RiboseQC analysis RData objects generated by \code{generate_rdata_list}.
#'
#' @return This function returns data.
#'
#' @seealso \code{\link{create_html_report}}

get_default_rl_selection <- function(rdata_list){
    datasets <- NULL
    for (i in seq_along(rdata_list)){
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
#' This data is used in the RiboseQC report in section 5 (displayed as table).
#'
#' @keywords RiboseQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#'
#' @param rdata_list List of RiboseQC analysis RData objects generated by \code{generate_rdata_list}.
#'
#' @return This function returns data to be displayed as table in the html report.
#'
#' @seealso \code{\link{create_html_report}}

get_top50_mapping <- function(rdata_list) {
    datasets <- list()
    for (i in c(seq_along(rdata_list))){
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
#' This data is used in the RiboseQC report in section 6 (displayed as table).
#'
#' @keywords RiboseQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#'
#' @param rdata_list List of RiboseQC analysis RData objects generated by \code{generate_rdata_list}.
#'
#' @return This function returns data to be displayed as table in the html report.
#'
#' @seealso \code{\link{create_html_report}}

get_top50_cds_genes <- function(rdata_list) {
    datasets <- list()
    for (i in c(seq_along(rdata_list))){
        res_all <- rdata_list[[names(rdata_list)[i]]]
        rc_cds <- as.data.frame(res_all$read_stats$counts_cds_genes)
        rc_cds <- rc_cds[with(rc_cds, order(-reads, chromosome)), ][seq_len(50),]
        datasets[[names(rdata_list)[i]]] <- rc_cds
    }
    return(datasets)
}

#' Get top 50 abundant genes (all genes)
#'
#' This function retrieves data on the top 50 abundant genes. \cr \cr
#' This data is used in the RiboseQC report in section 6 (displayed as table).
#'
#' @keywords RiboseQC
#' @author Dominique Sydow, \email{dominique.sydow@@posteo.de}
#'
#' @param rdata_list List of RiboseQC analysis RData objects generated by \code{generate_rdata_list}.
#'
#' @return This function returns data to be displayed as table in the html report.
#'
#' @seealso \code{\link{create_html_report}}

get_top50_all_genes <- function(rdata_list) {
    datasets <- list()
    for (i in c(seq_along(rdata_list))){
        res_all <- rdata_list[[names(rdata_list)[i]]]
        rc_all <- as.data.frame(res_all$read_stats$counts_all_genes)
        rc_all <- rc_all[with(rc_all, order(-reads, chromosome)), ][seq_len(100),]
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
#' to generate a bar plot, and in \code{plot_codon_usage_bulk_rmd} to iterativly generate plots for the RiboseQC report in section 7.2. \cr\cr
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
#' @keywords RiboseQC
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
#' These plots are displayed in the RiboseQC report in section 7.
#'
#' @keywords RiboseQC
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
#' @examples
#' data(res_all)
#' plot_codon_usage_positional_rmd(res_all, 'shoot', paste0('myfigpath', "rds/"))
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
            for(i in seq_along(data_types)){
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
#' @keywords RiboseQC
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
#' @examples
#' data(res_all)
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
#' These plots are displayed in the RiboseQC report in section 8.
#'
#' @keywords RiboseQC
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
#' @examples
#' data(res_all)
#' plot_codon_usage_bulk_rmd(res_all, 'shoot', paste0('myfigpath', "rds/"))
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

            for(i in seq_along(data_types)){
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
#' @keywords RiboseQC
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
#' @examples
#' data(res_all)
#' plot_codon_usage_bulk_rmd(res_all, 'shoot', paste0('myfigpath', "rds/"))
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
#' @keywords RiboseQC
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
#' @examples
#' data(res_all)
#' choose_readlengths(res_all$results_cutoffs,choice="max_coverage",nt_signals=res_all$profiles_fivepr[["five_prime_subcodon"]])
#' @export

choose_readlengths<-function(summary_data,choice="max_coverage",nt_signals){

    if(!choice%in%c("max_coverage","max_inframe","all")){
        stop("choice must be one of 'max_coverage','max_inframe' or 'all'")
    }
    stopifnot(names(summary_data) == names(nt_signals))
    choices_res<-list()
    for(i in names(summary_data)){
        choices_res[[i]]<-c()
        data<-DataFrame(summary_data[[i]])
        sign<-nt_signals[[i]]
        data<-data[grep(rownames(data),pattern = "all",invert = TRUE),]
        if(dim(data)[1]==0){next}
        data<-data[order(data$frame_preference,decreasing = TRUE),]
        mat_fr<-matrix(0,nrow = nrow(sign[[1]]),ncol = 99)
        net_all<-c()
        net_all_cov<-c()
        net_all_cov_ok<-c()
        gains_cod_all<-c()
        gains_sign_all<-c()
        gains_cod<-c()
        gains_sign<-c()
        for(cnt in seq_along(rownames(data))){
            if(cnt>length(rownames(data))){break}
            rl<-rownames(data)[cnt]
            cutf<-data$cutoff[cnt]
            sig<-as.matrix(sign[[rl]])
            #sig<-do.call(args = sig[c("5-UTR_1","5-UTR_2","CDS_1","CDS_2","3-UTR_1")],what = cbind.data.frame)
            #sig<-sig[,-c(1,2)]
            cdss<-grep(colnames(sig),pattern = "CDS")
            sig_ctf<-cbind(matrix(0,ncol = cutf,nrow = nrow(sig)),sig[,seq_len(length(sig[1,])-cutf)])[,cdss]
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
#' @keywords RiboseQC
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param x a \code{GAlignments} object with a cigar string
#' @param cutoff number representing the offset value
#' @return a \code{GRanges} object with offset reads

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
        for(i in seq_along(rangok)){shft<-c(shft,cutoff-cms[[i]][rangok[i]-1])}
        stt<-start(rangmore)
        stok<-c()
        for(i in seq_along(shft)){stok<-c(stok,stt[[i]][rangok[i]]+shft[i])}
        psmores<-GRanges(IRanges(start=stok,width = 1),seqnames = seqnames(x[mores]),strand=strand(x[mores]),seqlengths=seqlengths(x[mores]))

    }
    ps<-sort(c(psones,psmores))
    return(ps)
}

#' Offset spliced reads on minus strand
#'
#' This function calculates P-sites positions for spliced reads on the minus strand
#' @keywords RiboseQC
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param x a \code{GAlignments} object with a cigar string
#' @param cutoff number representing the offset value
#' @return a \code{GRanges} object with offset reads

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
        for(i in seq_along(rangok)){shft<-c(shft,cutoff-cms[[i]][rangok[i]-1])}
        #start?
        stt<-end(rangmore)
        stok<-c()
        for(i in seq_along(shft)){stok<-c(stok,stt[[i]][rangok[i]]-shft[i])}
        psmores<-GRanges(IRanges(start=stok,width = 1),seqnames = seqnames(x[mores]),strand=strand(x[mores]),seqlengths=seqlengths(x[mores]))

    }
    ps<-sort(c(psones,psmores))
    return(ps)


    if(rangok!=1){

        ps<-shift(resize(GRanges(rang[rangok],seqnames=seqnames(x),strand=strand(x),seqlengths=seqlengths(x)),width=1,fix="start"),shift=-(cutoff-sum(rang[seq_len(rangok-1)]@width)))
    }
    return(ps)
}

#' Load genomic features and genome sequence
#'
#' This function loads the annotation created by the \code{prepare_annotation_files function}
#' @keywords RiboseQC
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param path Full path to the *Rannot R file in the annotation directory used in the \code{prepare_annotation_files function}
#' @return introduces a \code{GTF_annotation} object and a \code{genome_seq} object in the parent environment
#' @seealso \code{\link{prepare_annotation_files}}
#' @examples
#' data(res_all)
#' load_annotation('test_arabidopsis.gtf.gz_Rannot')
#' @export

load_annotation<-function(path){

    GTF_annotation<-get(load(path))
	  haspackage <- isTRUE(is.character(GTF_annotation$genome))
    if(haspackage){
      genome_sequence<-get(library(GTF_annotation$genome,character.only = TRUE))
      library(GTF_annotation$genome,character.only = TRUE)
      genome_sequence<-get(GTF_annotation$genome)
	    message(paste0('assigning genome package ',GTF_annotation$genome,' to the global workspace as genome_seq'))
      assign('genome_seq',genome_sequence,envir = parent.frame())
    }else{
      genome_sequence<-GTF_annotation$genome
      message(paste0('assigning FaFile_Circ object ',GTF_annotation$genome$path,' to the global workspace as genome_seq'))
      assign('genome_seq',genome_sequence,envir = parent.frame())
    }

	  message(paste0('assigning GTF_annotation object GTF_annotation to parent workspace'))
    assign('GTF_annotation',GTF_annotation,envir = parent.frame())

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
#' @keywords RiboseQC
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @return a list with a \code{final_cutoff} object, the frame analysis containing the % of transcripts
#' displaying the max frame and the average % of frame precision per transcript in \code{frames_res},
#' all the calculated cutoffs in \code{cutoffs}, data used for the frame analysis in \code{frames},
#' and profiles around start codons in  \code{profiles_start}.
#' @seealso \code{\link{RiboseQC_analysis}}
#' @examples
#' data(profiles_fivepr)
#' calc_cutoffs_from_profiles(profiles_fivepr[["five_prime_subcodon"]][[1]][[1]],length_max = 100)
#' @export

calc_cutoffs_from_profiles<-function(reads_profile,length_max){
    profok<-as.matrix(reads_profile[,grep("CDS",x = colnames(reads_profile))])
    #profok<-do.call(args = reads_profile[c("CDS_1","CDS_2")],what = cbind.data.frame)
    sums<-rowSums(profok)
    profok<-profok[sums>0,,drop=FALSE]
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
    names(frameans)<-seq_len(3)
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
    profok<-profok[sums>0,,drop=FALSE]
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
    final_cutoff<-as.numeric(names(sort(table(cutt[grep(names(cutt),pattern = "frame")]),decreasing = TRUE)[1]))
    frames_res<-c(sort(tbfra,decreasing = TRUE)[1],as.numeric(frameans[which.max(frameans)]))
    names(frames_res)<-c("max_frame_in_pctORF","mean_pct_max_frame")

    results<-list(final_cutoff=final_cutoff,frames_res=frames_res,cutoffs=cutoffs,frames=frames,profiles_start=profok)
    return(results)

}


#' Prepare comprehensive sets of annotated genomic features
#'
#' This function processes a gtf file and a twobit file (created using faToTwoBit from ucsc tools: http://hgdownload.soe.ucsc.edu/admin/exe/ ) to create a com
#' prehensive set of genomic regions of interest in genomic and transcriptomic space (e.g. introns, UTRs, start/stop codons).
#'    In addition, by linking genome sequence and annotation, it extracts additional info, such as gene and transcript biotypes, genetic codes for different organelles, or chromosomes and transcripts lengths.
#' @keywords RiboseQC
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param annotation_directory The target directory which will contain the output files
#' @param twobit_file Full path to the genome file in twobit format
#' @param gtf_file Full path to the annotation file in GTF format
#' @param scientific_name A name to give to the organism studied; must be two words separated by a ".", defaults to Homo.sapiens
#' @param annotation_name A name to give to annotation used; defaults to genc25
#' @param export_bed_tables_TxDb Export coordinates and info about different genomic regions in the annotation_directory? It defaults to \code{TRUE}
#' @param forge_BSgenome Forge and install a \code{BSgenome} package? It defaults to \code{TRUE}
#' @param create_TxDb Create a \code{TxDb} object and a *Rannot object? It defaults to \code{TRUE}
#' @param annot_file specify an exact file name for the rds file created by this function, defaults to annotation_directory/basename(gtf)_Rannot
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
#' \code{genome}: the name of the forged BSgenome package, or an FaFile_Circ object. Loaded with \code{load_annotation} function.\cr
#' \code{stop_in_gtf}: stop codon, as defined in the annotation.\cr
#' @return a TxDb file and a *Rannot files are created in the specified \code{annotation_directory}.
#' In addition, a BSgenome object is forged, installed, and linked to the *Rannot object
#' @seealso \code{\link{load_annotation}}, \code{\link{forgeBSgenomeDataPkg}}, \code{\link{makeTxDbFromGFF}}.
#' @examples
#' gtf_file <- system.file("extdata", "example.gtf",
#' package = "RiboseQC",mustWork = TRUE)
#' twobit_file <- system.file("extdata", "example.2bit",
#'                            package = "RiboseQC",mustWork = TRUE)
#' output_dir <- tempdir()
#' prepare_annotation_files(output_dir,twobit_file,gtf_file,
#'                    scientific_name = "Arabidopsis.thaliana",
#'                    annotation_name = "tair10.42" )
#'
#' @export


prepare_annotation_files<-function(annotation_directory,twobit_file=NULL,gtf_file,scientific_name="Homo.sapiens",
                                   annotation_name="genc25",export_bed_tables_TxDb=TRUE,forge_BSgenome=FALSE,genome_seq=NULL,circ_chroms=DEFAULT_CIRC_SEQS,create_TxDb=TRUE,annot_file=NULL){


    DEFAULT_CIRC_SEQS <- unique(c("chrM","MT","MtDNA","mit","Mito","mitochondrion",
                                  "dmel_mitochondrion_genome","Pltd","ChrC","Pt","chloroplast",
                                  "Chloro","2micron","2-micron","2uM",
                                  "Mt", "NC_001879.2", "NC_006581.1","ChrM","mitochondrion_genome"))
    #adjust variable names (some chars not permitted)
    annotation_name<-gsub(annotation_name,pattern = "_",replacement = "")
    annotation_name<-gsub(annotation_name,pattern = "-",replacement = "")
    if(!dir.exists(annotation_directory)){dir.create(path = annotation_directory,recursive = TRUE)}
    annotation_directory<-normalizePath(annotation_directory)
    gtf_file<-normalizePath(gtf_file)

    for (f in c(gtf_file)){
        if(file.access(f, 0)==-1) {
            stop("
                 The following files don't exist:\n",
                 f, "\n")
        }
    }


    #get circular sequences

    #Forge a BSGenome package

    if(forge_BSgenome){
        stopifnot(!is.null(twobit_file))
        scientific_name_spl<-strsplit(scientific_name,"[.]")[[1]]
        ok<-length(scientific_name_spl)==2
        if(!ok){stop("\"scientific_name\" must be two words separated by a \".\", like \"Homo.sapiens\"")}


        twobit_file<-normalizePath(twobit_file)

        for (f in c(twobit_file)){
            if(file.access(f, 0)==-1) {
                stop("
                     The following files don't exist:\n",
                     f, "\n")
            }
        }

        seqinfotwob<-seqinfo(TwoBitFile(twobit_file))
        circss<-seqnames(seqinfotwob)[which(seqnames(seqinfotwob)%in%circ_chroms)]
        seqinfotwob@is_circular[which(seqnames(seqinfotwob)%in%circ_chroms)]<-TRUE


        circseed<-circss
        if(length(circseed)==0){circseed<-NULL}

        pkgnm<-paste("BSgenome",scientific_name,annotation_name,sep=".")


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


        if(length(circseed)==0){
            writeLines(text = seed_text,con = seed_dest)
        }

        if(length(circseed)==1){
            seed_text<-paste(seed_text,"\n",
                             "circ_seqs: \"",circseed,"\"",sep="")
        }

        if(length(circseed)>1){
            circseed<-paste('c("',paste(circseed,collapse=","),'")',sep="")
            circseed<-gsub(circseed,pattern = ",",replacement='","')

            cat(seed_text,"\n","circ_seqs: ",circseed,"\n",sep="",file = seed_dest)
        }

        writeLines(text = seed_text,con = seed_dest)

        unlink(paste(annotation_directory,pkgnm,sep="/"),recursive=TRUE)

        forgeBSgenomeDataPkg(x=seed_dest,destdir=annotation_directory,seqs_srcdir=dirname(twobit_file))
        cat(paste("Creating the BSgenome package --- Done! ",date(),"\n",sep = ""))

        cat(paste("Installing the BSgenome package ... ",date(),"\n",sep = ""))

        install(paste(annotation_directory,pkgnm,sep="/"),upgrade = FALSE)
        cat(paste("Installing the BSgenome package --- Done! ",date(),"\n",sep = ""))

        seqinfo_genome <- seqinfotwob
    }else{
        if(!is(genome_seq,'FaFile')){
            genome_seq <- Rsamtools::FaFile(genome_seq)
        }
        Rsamtools::indexFa(genome_seq)
        if(!is(genome_seq,'FaFile_Circ')){
            genome_seq <- FaFile_Circ(genome_seq,circularRanges=circ_chroms)
        }
        seqinfo_genome<-seqinfo(genome_seq)
        seqinfo_genome@is_circular[which(seqnames(seqinfo_genome)%in%circ_chroms)]<-TRUE
    }
    #Create the TxDb from GTF and BSGenome info

    if(create_TxDb){
        cat(paste("Creating the TxDb object ... ",date(),"\n",sep = ""))

        annotation<-makeTxDbFromGFF(file=gtf_file,format="gtf",chrominfo = seqinfo_genome)

        saveDb(annotation, file=paste(annotation_directory,"/",basename(gtf_file),"_TxDb",sep=""))
        cat(paste("Creating the TxDb object --- Done! ",date(),"\n",sep = ""))
        cat(paste("Extracting genomic regions ... ",date(),"\n",sep = ""))


        genes<-genes(annotation)
        exons_ge<-exonsBy(annotation,by="gene")
        exons_ge<-GenomicRanges::reduce(exons_ge)

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

        cds_tx<-cdsBy(annotation,"tx",use.names=TRUE)

        #filter out abnormally short cds (I"m looking at you maize annotation)
        cds_tx <- cds_tx[sum(width(cds_tx))>=3]
        txs_gene<-transcriptsBy(annotation,by="gene")
        genes_red<-reduce(sort(genes(annotation)))

        exons_tx<-exonsBy(annotation,"tx",use.names=TRUE)

        transcripts_db<-transcripts(annotation)
        intron_names_tx<-intronsByTranscript(annotation,use.names=TRUE)


        #define exonic bins, including regions overlapping multiple genes
        # nsns<-disjointExons(annotation,aggregateGenes=TRUE)
   	nsns<-exonicParts(annotation, linked.to.single.gene.only=F)
  	mcols(nsns) <- DataFrame(mcols(nsns)[c(3,2)],exonic_part=NA) #disjointExons is deprecated, so use exonicParts instead and mimic back to the disjoint format in ORFquant



        #define tx_coordinates of ORF boundaries

        exsss_cds<-exons_tx[names(cds_tx)]
        chunks<-seq(1,length(cds_tx),by = 20000)
        if(chunks[length(chunks)]<length(cds_tx)){chunks<-c(chunks,length(cds_tx))}
        mapp<-GRangesList()
        for(i in seq_len(length(chunks)-1)){
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


         gtfdata <- import.gff2(gtf_file,colnames=c("gene_id","gene_biotype","gene_type","gene_name","gene_symbol","transcript_id","transcript_biotype","transcript_type","type"))
         n_transcripts = length(unique(gtfdata$transcript_id))
         stopifnot('transcript' %in% gtfdata$type)
         gtfdata <- subset(gtfdata, type=='transcript')
         gtfdata$type <- NULL
         stopifnot(length(unique(gtfdata$transcript_id))==n_transcripts)

         trann<-unique(mcols(gtfdata))


        trann<-trann[!is.na(trann$transcript_id),]
        trann<-data.frame(unique(trann),stringsAsFactors=FALSE)



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
        g<-unlist(apply(cbind(seq_along(mat_genes),Y = elementNROWS(mat_genes)),FUN =function(x) rep(x[1],x[2]),MARGIN = 1))
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

        if(forge_BSgenome){
            suppressPackageStartupMessages(library(pkgnm,character.only=TRUE))
            genome<-get(pkgnm)
        }else{
            stopifnot(!is.null(genome_seq))
            genome<-genome_seq
            pkgnm<-NA
        }
        tocheck<-as.character(runValue(seqnames(cds_tx)))
        tocheck<-cds_tx[!tocheck%in%circs]
        width(tocheck)%>%sum%>%.[.<3]
        seqcds<-extractTranscriptSeqs(genome,transcripts = tocheck)
        cd<-unique(translations$genetic_code[!rownames(translations)%in%circs])
        trsl<-suppressWarnings(translate(seqcds,genetic.code = getGeneticCode(cd),if.fuzzy.codon = "solve"))
        trslend<-as.character(narrow(trsl,end = width(trsl),width = 1))
        stop_inannot<-NA
        if(names(sort(table(trslend),decreasing = TRUE)[1])=="*"){stop_inannot<-"*"}

        cds_txscoords$gene_id<-trann$gene_id[match(as.vector(seqnames(cds_txscoords)),trann$transcript_id)]
        cds_cc<-cds_txscoords
        strand(cds_cc)<-"*"
        sta_cc<-resize(cds_cc,width = 1,"start")
        sta_cc<-unlist(pmapFromTranscripts(sta_cc,exons_tx[seqnames(sta_cc)],ignore.strand=FALSE))
        sta_cc$gene_id<-trann$gene_id[match(names(sta_cc),trann$transcript_id)]
        sta_cc<-sta_cc[sta_cc$hit]
        strand(sta_cc)<-structure(as.vector(strand(transcripts_db)),names=transcripts_db$tx_name)[names(sta_cc)]
        sta_cc$type<-"start_codon"
        mcols(sta_cc)<-mcols(sta_cc)[,c("exon_rank","type","gene_id")]

        sto_cc<-resize(cds_cc,width = 1,"end")
        #stop codon is the 1st nt, e.g. U of the UAA
        #To-do: update with regards to different organelles, and different annotations
        sto_cc<-shift(sto_cc,-2)

        if(is.na(stop_inannot)){sto_cc<-resize(trim(shift(sto_cc,3)),width = 1,fix = "end")}

        sto_cc<-unlist(pmapFromTranscripts(sto_cc,exons_tx[seqnames(sto_cc)],ignore.strand=FALSE))
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
            x<-x[order(x$var,x$utr5len,x$cdslen,decreasing = TRUE),]
            x<-x[x$var==max(x$var),]
            ok<-x$txid[which(x$cdslen==max(x$cdslen) & x$utr5len==max(x$utr5len) & x$var==max(x$var))][1]
            if(length(ok)==0 | is.na(ok[1])){ok<-x$txid[1]}
            ok
        })

        df<-cbind.DataFrame(as.character(seqnames(cds_txscoords)),width(cds_txscoords),start(cds_txscoords),cds_txscoords$upstr_stasto,cds_txscoords$gene_id)
        colnames(df)<-c("txid","cdslen","utr5len","var","gene_id")
        repres_upstr<-by(df[,c("txid","cdslen","utr5len","var")],df$gene_id,function(x){
            x<-x[order(x$var,x$utr5len,x$utr5len,decreasing = TRUE),]
            x<-x[x$var==max(x$var),]
            ok<-x$txid[which(x$cdslen==max(x$cdslen) & x$utr5len==max(x$utr5len) & x$var==max(x$var))][1]
            if(length(ok)==0 | is.na(ok[1])){ok<-x$txid[1]}
            ok
        })
        df<-cbind.DataFrame(as.character(seqnames(cds_txscoords)),width(cds_txscoords),start(cds_txscoords),cds_txscoords$upstr_stasto,cds_txscoords$gene_id)
        colnames(df)<-c("txid","cdslen","utr5len","var","gene_id")
        repres_len5<-by(df[,c("txid","cdslen","utr5len","var")],df$gene_id,function(x){
            x<-x[order(x$utr5len,x$var,x$cdslen,decreasing = TRUE),]
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
        pkgnm_or_faob<- if(is(genome_seq,'FaFile') ) {genome_seq} else {pkgnm}
        GTF_annotation<-list(transcripts_db,txs_gene,ifs,unq_stst,cds_tx,intron_names_tx,cds_gen,exons_tx,nsns,unq_intr,genes,threeutrs,fiveutrs,ncisof,ncrnas,introns,intergenicRegions,trann,cds_txscoords,translations,pkgnm_or_faob,stop_inannot)
        names(GTF_annotation)<-c("txs","txs_gene","seqinfo","start_stop_codons","cds_txs","introns_txs","cds_genes","exons_txs","exons_bins","junctions","genes","threeutrs","fiveutrs","ncIsof","ncRNAs","introns","intergenicRegions","trann","cds_txs_coords","genetic_codes","genome","stop_in_gtf")

        #Save as a RData object
        if(is.null(annot_file)){
            annot_file <- paste(annotation_directory,"/",basename(gtf_file),"_Rannot",sep="")
        }
        save(GTF_annotation,file=annot_file)
        cat(paste("Rannot object created!   ",date(),"\n",sep = ""))
        GTF_annotation

        #create tables and bed files (with colnames, so with header)
        if(export_bed_tables_TxDb==TRUE){
            cat(paste("Exporting annotation tables ... ",date(),"\n",sep = ""))
            for(bed_file in c("fiveutrs","threeutrs","ncIsof","ncRNAs","introns","cds_txs_coords")){
                bf<-GTF_annotation[[bed_file]]
                bf_t<-data.frame(chromosome=seqnames(bf),start=start(bf),end=end(bf),name=rep(".",length(bf)),score=width(bf),strand=strand(bf))
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

    return(annot_file)
}

