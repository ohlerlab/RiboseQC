#' @include riboseqc.R
#' @import Rsamtools
NULL


#' Perform a Ribo-seQC analysis
#'
#' This function loads annotation created by the prepare_annotation_files function, and analyzes a BAM file.
#' @keywords Ribo-seQC
#' @author Lorenzo Calviello, \email{calviello.l.bio@@gmail.com}
#' @param annotation_file Full path to the annotation file (*Rannot). Or, a vector with paths to one annotation file per bam file.
#' @param bam_files character vector containing the full path to the bam files
#' @param chunk_size the number of alignments to read at each iteration, defaults to 5000000, increase when more RAM is available. Must be between 10000 and 100000000
#' @param read_subset Select readlengths up to 99 percent of the reads, defaults to \code{TRUE}. Must be of length 1 or same length as bam_files.
#' @param readlength_choice_method Method used to subset relevant read lengths (see \code{choose_readlengths} function); defaults to 'max_coverage'. Must be of length 1 or same length as bam_files.
#' @param rescue_all_rls Set cutoff of 12 for read lengths ignored because of insufficient coverage. Defaults to \code{FALSE}. Must be of length 1 or same length as bam_files.
#' @param write_tmp_files Should output all the results (in *results_RiboseQC_all)? Defaults to \code{TRUE}. Must be of length 1 or same length as bam_files.
#' @param dest_names character vector containing the prefixes to use for the result output files. Defaults to same as \code{bam_files}
#' @param fast_mode Use only top 500 genes to build profiles? Defaults to \code{TRUE}. Must be of length 1 or same length as bam_files.
#' @param create_report Create an html report showing the RiboseQC analysis results. Defaults to \code{TRUE}
#' @param sample_names character vector containing the names for each sample analyzed (for the html report). Defaults to 'sample1', 'sample2' ...
#' @param report_file desired filename for for the html report file. Defaults to the first entry of \code{bam_files} followed by '.html'
#' @param extended_report creates a large html report including codon occupancy for each read length. Defaults to \code{FALSE}
#' @param pdf_plots creates a pdf file for each produced plot. Defaults to \code{TRUE}
#' @param stranded are the analyzed libraries strand-specific? TRUE, FALSE or 'inverse'. Defaults to \code{TRUE}
#' @param normalize_cov export normalized (sum to 1 million) bedgraph files for coverage tracks? Defaults to \code{TRUE}
#' @param offsets_df optionally input an offsets_df created externally, note that this does not yet fully integrate with stats such as 'proportion in phase' it will however alter the psites which are outputted. Must have cols 'read_length','cutoff'
#' @param genome_seq An FaFile object, to be used instead of a BSgenome package
#' @return the function saves a 'results_RiboseQC_all' R file appended to the bam_files path including the complete list of outputs described here.
#' In addition, bedgraph files for coverage value and P_sites position is appended to the bam_files path, including also a summary of P_sites selection statistics,
#' a smaller 'results_RiboseQC' R file used for creating a dynamic html report, and a 'for_ORFquant' R object that can be used in the ORFquant pipeline.
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
#' @import rmarkdown
#' @import rtracklayer
#' @import GenomicAlignments
#' @import BSgenome
#' @import GenomicFiles
#' @import devtools
#' @import reshape2
#' @import ggplot2
#' @import knitr
#' @import DT
#' @import gridExtra
#' @import ggpubr
#' @import viridis
#' @import Biostrings
#' @import GenomicFeatures
#' @import BiocGenerics
#' @import GenomicRanges
#' @export 

RiboseQC_analysis <- function(annotation_file, bam_files, read_subset = TRUE, readlength_choice_method = "max_coverage", 
    genome_seq = NULL, stranded = TRUE, normalize_cov = TRUE, chunk_size = 5000000L, 
    write_tmp_files = TRUE, dest_names = NA, rescue_all_rls = FALSE, fast_mode = TRUE, 
    create_report = TRUE, sample_names = NA, report_file = NA, extended_report = FALSE, 
    pdf_plots = TRUE, offsets_df = NULL) {
    
    if (length(dest_names) == 1) {
        if (is.na(dest_names)) {
            dest_names = bam_files
        }
    }
    
    if (sum(duplicated(dest_names)) > 0) {
        stop(paste("duplicated 'dest_names' parameter (possibly same bam file). Please specify different 'dest_names'.", 
            date(), "\n"))
    }
    
    if (is.na(report_file)) {
        report_file = paste(bam_files[1], "_RiboseQC_report.html", sep = "")
    }
    if (length(sample_names) == 1) {
        if (is.na(sample_names)) {
            sample_names = paste("sample", seq_along(bam_files), sep = "")
        }
    }
    
    if (!is.logical(read_subset)) {
        stop(paste("'read_subset' must be logical (TRUE or FALSE, no quotes).", date(), 
            "\n"))
    }
    if (!is.logical(write_tmp_files)) {
        stop(paste("'write_tmp_files' must be logical (TRUE or FALSE, no quotes).", 
            date(), "\n"))
    }
    if (!is.logical(rescue_all_rls)) {
        stop(paste("'rescue_all_rls' must be logical (TRUE or FALSE, no quotes).", 
            date(), "\n"))
    }
    if (!is.logical(fast_mode)) {
        stop(paste("'fast_mode' must be logical (TRUE or FALSE, no quotes).", date(), 
            "\n"))
    }
    if (!is.logical(create_report)) {
        stop(paste("'create_report' must be logical (TRUE or FALSE, no quotes).", 
            date(), "\n"))
    }
    if (!is.integer(chunk_size) & !is.numeric(chunk_size)) {
        stop(paste("'chunk_size' must be an integer value (e.g. 5000000).", date(), 
            "\n"))
    }
    if (chunk_size <= 1e+05 | chunk_size > 1e+08) {
        stop(paste("'chunk_size' must be an integer value between 100000 and 100000000 (e.g. 5000000).", 
            date(), "\n"))
    }
    if (!is.character(annotation_file)) {
        stop(paste("'annotation_file' must be a character (e.g.'/home/myfolder/myfile').", 
            date(), "\n"))
    }
    if (!is.character(bam_files)) {
        stop(paste("'bam_files' must be a character (e.g.'/home/myfolder/myfile') or character vector (e.g. c('/home/myfolder/myfile1','/home/myfolder/myfile2').", 
            date(), "\n"))
    }
    if (!is.character(dest_names)) {
        stop(paste("'dest_names' must be a character (e.g.'/home/myfolder/myfile') or character vector (e.g. c('/home/myfolder/myfile1','/home/myfolder/myfile2')..", 
            date(), "\n"))
    }
    if (!is.character(report_file)) {
        stop(paste("'report_file' must be a character (e.g.'/home/myfolder/myfile').", 
            date(), "\n"))
    }
    if (!is.character(sample_names)) {
        stop(paste("'sample_names' must be a character (e.g.'failed_experiment'). or character vector (e.g. c('attempt_1','attempt_2').", 
            date(), "\n"))
    }
    if (!is.logical(normalize_cov)) {
        stop(paste("'normalize_cov' must be logical (TRUE or FALSE, no quotes).", 
            date(), "\n"))
    }
    
    if (!length(read_subset) %in% c(1, length(bam_files))) {
        stop(paste("length of 'read_subset' must be 1 or same as 'bam_files'.", date(), 
            "\n"))
    }
    if (!length(readlength_choice_method) %in% c(1, length(bam_files))) {
        stop(paste("length of 'readlength_choice_method' must be 1 or same as 'bam_files'.", 
            date(), "\n"))
    }
    if (!length(rescue_all_rls) %in% c(1, length(bam_files))) {
        stop(paste("length of 'rescue_all_rls' must be 1 or same as 'bam_files'.", 
            date(), "\n"))
    }
    if (!length(fast_mode) %in% c(1, length(bam_files))) {
        stop(paste("length of 'fast_mode' must be 1 or same as 'bam_files'.", date(), 
            "\n"))
    }
    if (!length(stranded) %in% c(1, length(bam_files))) {
        stop(paste("length of 'stranded' must be 1 or same as 'bam_files'.", date(), 
            "\n"))
    }
    if (!length(normalize_cov) %in% c(1, length(bam_files))) {
        stop(paste("length of 'normalize_cov' must be 1 or same as 'bam_files'.", 
            date(), "\n"))
    }
    if (length(bam_files) > 1) {
        if (length(read_subset) == 1) {
            read_subset <- rep(read_subset, length(bam_files))
        }
        if (length(readlength_choice_method) == 1) {
            readlength_choice_method <- rep(readlength_choice_method, length(bam_files))
        }
        if (length(rescue_all_rls) == 1) {
            rescue_all_rls <- rep(rescue_all_rls, length(bam_files))
        }
        if (length(fast_mode) == 1) {
            fast_mode <- rep(fast_mode, length(bam_files))
        }
        if (length(stranded) == 1) {
            stranded <- rep(stranded, length(bam_files))
        }
        if (length(normalize_cov) == 1) {
            normalize_cov <- rep(normalize_cov, length(bam_files))
        }
    }
    
    
    
    if (!length(dest_names) %in% c(length(bam_files))) {
        stop(paste("length of 'dest_names' must be same as 'bam_files'.", date(), 
            "\n"))
    }
    if (!length(sample_names) %in% c(length(bam_files))) {
        stop(paste("length of 'sample_names' must be same as 'bam_files'.", date(), 
            "\n"))
    }
    
    if (sum(!readlength_choice_method %in% c("max_coverage", "max_inframe", "all")) > 
        0) {
        stop(paste("'readlength_choice_method' must be one of 'max_coverage','max_inframe','all'.", 
            date(), "\n"))
    }
    
    fastainput <- is(genome_seq, "FaFile")
    
    if (fastainput) {
        stopifnot(file.exists(genome_seq))
        genome_seq <- Rsamtools::FaFile(genome_seq)
        Rsamtools::indexFa(genome_seq)
        genome_seq <- FaFile_Circ(genome_seq)
    }
    
    list_annotations <- list()
    
    for (annots in seq_along(annotation_file)) {
        cat(paste("Loading annotation files in ", annotation_file[annots], " ... ", 
            date(), "\n", sep = ""))
        
        load_annotation(annotation_file[annots])
        lst <- list(GTF_annotation, genome_seq)
        names(lst) <- c("GTF_annotation", "genome_seq")
        list_annotations[[annots]] <- lst
        
    }
    
    cat(paste("Loading annotation files --- Done! ", date(), "\n", sep = ""))
    
    GTF_annotation <- list_annotations[[1]]$GTF_annotation
    genome_seq <- list_annotations[[1]]$genome_seq
    
    
    resfilelist <- c()
    
    for (bammo in seq_along(bam_files)) {
        
        chunk_size <- as.integer(chunk_size)
        
        GTF_annotation <- list_annotations[[1]]$GTF_annotation
        genome_seq <- list_annotations[[1]]$genome_seq
        
        if (length(annotation_file) > 1 & bammo > 1) {
            GTF_annotation <- list_annotations[[bammo]]$GTF_annotation
            genome_seq <- list_annotations[[bammo]]$genome_seq
        }
        
        bam_file <- bam_files[bammo]
        dest_name <- dest_names[bammo]
        
        opts <- BamFile(file = bam_file, yieldSize = chunk_size)  #read BAMfile in one chunk
        param <- ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE, isSecondaryAlignment = FALSE), 
            what = c("mapq"), tag = "MD")
        
        dira <- paste(dirname(dest_name), paste("tmp_RiboseQC", basename(dest_name), 
            sep = "_"), sep = "/")
        suppressWarnings(dir.create(dira, recursive = TRUE))
        
        seqs <- seqinfo(opts)
        circs_seq <- seqnames(GTF_annotation$seqinfo)[which(isCircular(GTF_annotation$seqinfo))]
        circs <- seqs@seqnames[which(seqs@seqnames %in% circs_seq)]
        
        red_cdss <- reduce(GTF_annotation$cds_genes)
        
        red_ex <- unlist(GTF_annotation$exons_txs)
        red_ex$gene_id <- GTF_annotation$trann$gene_id[match(names(red_ex), GTF_annotation$trann$transcript_id)]
        red_ex <- reduce(split(red_ex, red_ex$gene_id))
        
        regions <- list(reduce(unlist(GTF_annotation$cds_genes)), GTF_annotation$fiveutrs, 
            GTF_annotation$threeutrs, GTF_annotation$ncIsof, GTF_annotation$ncRNAs, 
            GTF_annotation$introns, GTF_annotation$intergenicRegions)
        names(regions) <- c("cds", "fiveutrs", "threeutrs", "ncIsof", "ncRNAs", "introns", 
            "intergenic")
        
        regions <- GRangesList(endoapply(regions, function(x) {
            mcols(x) <- NULL
            x
        }))
        
        list_locat <- list()
        list_locat[["nucl"]] <- regions
        all <- regions
        if (length(circs) > 0) {
            for (i in c("nucl", circs)) {
                if (i == "nucl") {
                  regions <- all
                  for (w in seq_along(all)) {
                    regions[[w]] <- regions[[w]][!seqnames(regions[[w]]) %in% circs]
                  }
                  list_locat[["nucl"]] <- regions
                }
                if (i != "nucl") {
                  regions <- all
                  for (w in seq_along(all)) {
                    regions[[w]] <- regions[[w]][seqnames(regions[[w]]) %in% i]
                  }
                  list_locat[[i]] <- regions
                }
            }
        }
        
        # get readlengths
        
        reduc <- function(x, y) {
            list(max(x[[1]], y[[1]]), min(x[[2]], y[[2]]))
        }
        yiel <- function(x) {
            readGAlignments(x, param = param)
        }
        
        mapp <- function(x) {
            x_I <- x[grep("I", cigar(x))]
            
            if (length(x_I) > 0) {
                x <- x[grep("I", cigar(x), invert = TRUE)]
                
            }
            x_D <- x[grep("D", cigar(x))]
            if (length(x_D) > 0) {
                x <- x[grep("D", cigar(x), invert = TRUE)]
                
            }
            
            
            clipp <- width(cigarRangesAlongQuerySpace(x@cigar, ops = "S"))
            clipp[elementNROWS(clipp) == 0] <- 0
            len_adj <- qwidth(x) - sum(clipp)
            mcols(x)$len_adj <- len_adj
            
            
            maxr <- max(mcols(x)$len_adj)
            minr <- min(mcols(x)$len_adj)
            list(maxr, minr)
        }
        
        cat(paste("Extracting read lengths from ", bam_file, " ... ", date(), "\n", 
            sep = ""))
        maxmin <- reduceByYield(X = opts, YIELD = yiel, MAP = mapp, REDUCE = reduc)
        cat(paste("Extracting read lengths --- Done! ", date(), "\n", sep = ""))
        readlengths <- seq(maxmin[[2]], maxmin[[1]])
        
        
        
        # body of the iteration for reducebyYield
        strandedness <- stranded[bammo]
        
        # what to do with chunks: x present chunk, y old chunks (cumulative)
        reduc <- function(x, y) {
            all_ps <- GRangesList()
            rls <- unique(c(names(x[["reads_pos1"]]), names(y[["reads_pos1"]])))
            seql <- seqlevels(GTF_annotation$seqinfo)
            seqle <- seqlengths(GTF_annotation$seqinfo)
            for (rl in rls) {
                reads_x <- GRanges()
                reads_y <- GRanges()
                
                if (sum(rl %in% names(x[["reads_pos1"]])) > 0) {
                  reads_x <- x[["reads_pos1"]][[rl]]
                }
                if (sum(rl %in% names(y[["reads_pos1"]])) > 0) {
                  reads_y <- y[["reads_pos1"]][[rl]]
                }

                #This needs to go after the above                
                seqlevels(reads_x) <- seql
                seqlevels(reads_y) <- seql
                
                seqlengths(reads_x) <- seqle
                seqlengths(reads_y) <- seqle
                
                
#                reads_y <- 
                plx <- reads_x[strand(reads_x) == "+"]
                mnx <- reads_x[strand(reads_x) == "-"]
                ply <- reads_y[strand(reads_y) == "+"]
                mny <- reads_y[strand(reads_y) == "-"]
                if (length(plx) > 0) {
                  covv_pl <- coverage(plx, weight = plx$score)
                } else {
                  covv_pl <- coverage(plx)
                }
                if (length(ply) > 0) {
                  covv_pl <- covv_pl + coverage(ply, weight = ply$score)
                }
                covv_pl <- GRanges(covv_pl)
                covv_pl <- covv_pl[covv_pl$score > 0]
                
                if (length(mnx) > 0) {
                  covv_min <- coverage(mnx, weight = mnx$score)
                } else {
                  covv_min <- coverage(mnx)
                }
                if (length(mny) > 0) {
                  covv_min <- covv_min + coverage(mny, weight = mny$score)
                }
                
                covv_min <- GRanges(covv_min)
                covv_min <- covv_min[covv_min$score > 0]
                
                strand(covv_pl) <- "+"
                strand(covv_min) <- "-"
                
                all_ps[[rl]] <- sort(c(covv_pl, covv_min))
            }
            
            
            cntsss <- y[["counts_cds_genes"]]
            cntsss$reads <- cntsss$reads + x[["counts_cds_genes"]]$reads
            
            cntsss_all <- y[["counts_all_genes"]]
            cntsss_all$reads <- cntsss_all$reads + x[["counts_all_genes"]]$reads
            
            
            cntsss_unq <- y[["counts_cds_genes_unq"]]
            cntsss_unq$reads <- cntsss_unq$reads + x[["counts_cds_genes_unq"]]$reads
            
            cntsss_all_unq <- y[["counts_all_genes_unq"]]
            cntsss_all_unq$reads <- cntsss_all_unq$reads + x[["counts_all_genes_unq"]]$reads
            
            
            reads_summary <- mapply(FUN = function(x, y) {
                DataFrame(as.matrix(x) + as.matrix(y))
            }, x[["reads_summary"]], y[["reads_summary"]], SIMPLIFY = FALSE)
            reads_summary_unq <- mapply(FUN = function(x, y) {
                DataFrame(as.matrix(x) + as.matrix(y))
            }, x[["reads_summary_unq"]], y[["reads_summary_unq"]], SIMPLIFY = FALSE)
            
            lis <- list((x[["rld"]] + y[["rld"]]), x[["rld_unq"]] + y[["rld_unq"]], 
                (x[["positions"]] + y[["positions"]]), (x[["positions_unq"]] + y[["positions_unq"]]), 
                all_ps, cntsss, cntsss_unq, cntsss_all, cntsss_all_unq, reads_summary, 
                reads_summary_unq)
            names(lis) <- c("rld", "rld_unq", "positions", "positions_unq", "reads_pos1", 
                "counts_cds_genes", "counts_cds_genes_unq", "counts_all_genes", "counts_all_genes_unq", 
                "reads_summary", "reads_summary_unq")
            lis
        }
        
        # what to do with each chunk (read as alignment file)
        yiel <- function(x) {
            readGAlignments(x, param = param)
        }
        
        # operations on the chunk (here count reads and whatnot)
        mapp <- function(x) {
            if (strandedness == "inverse") {
                x <- invertStrand(x)
                strandedness <- F
            }
            
            mcols(x)$MD[which(is.na(mcols(x)$MD))] <- "NO"
            
            # Remove Insertion and Deletion (TBA)
            
            x_I <- x[grep("I", cigar(x))]
            
            if (length(x_I) > 0) {
                x <- x[grep("I", cigar(x), invert = TRUE)]
                
            }
            x_D <- x[grep("D", cigar(x))]
            if (length(x_D) > 0) {
                x <- x[grep("D", cigar(x), invert = TRUE)]
                
            }
            
            emptyt <- rep(0, length(readlengths))
            names(emptyt) <- paste("reads", readlengths, sep = "_")
            
            
            # softclipping part
            
            
            clipp <- width(cigarRangesAlongQuerySpace(x@cigar, ops = "S"))
            clipp[elementNROWS(clipp) == 0] <- 0
            len_adj <- qwidth(x) - sum(clipp)
            mcols(x)$len_adj <- len_adj
            
            # Remove S from Cigar (read positions/length are already adjusted) it helps
            # calculating P-sites positions for spliced reads
            
            cigg <- cigar(x)
            cigg_s <- grep(cigg, pattern = "S")
            if (length(cigg_s) > 0) {
                cigs <- cigg[cigg_s]
                cigs <- gsub(cigs, pattern = "^[0-9]+S", replacement = "")
                cigs <- gsub(cigs, pattern = "[0-9]+S$", replacement = "")
                cigg[cigg_s] <- cigs
                x@cigar <- cigg
            }
            mcols(x)$cigar_str <- x@cigar
            x_uniq <- x[x@elementMetadata$mapq > 50]
            
            
            
            # rld_len, rld_loc
            
            
            # initialise read length vector per compartment
            list_vects <- list(nucl = emptyt)
            if (length(circs) > 0) {
                for (i in circs) {
                  list_vects[[i]] <- emptyt
                }
            }
            
            list_vects_unq <- list(nucl = emptyt)
            if (length(circs) > 0) {
                for (i in circs) {
                  list_vects_unq[[i]] <- emptyt
                }
            }
            
            # sort reads per compartment
            list_reads <- list(nucl = x)
            list_reads_unq <- list(nucl = x_uniq)
            if (length(circs) > 0) {
                nucl <- subset(x, !seqnames(x) %in% circs)
                nucl_unq <- subset(x_uniq, !seqnames(x_uniq) %in% circs)
                list_reads <- list(nucl = nucl)
                for (i in circs) {
                  list_reads[[i]] <- subset(x, seqnames(x) == i)
                  list_reads_unq[[i]] <- subset(x_uniq, seqnames(x_uniq) == i)
                }
            }
            
            # count reads per read length per compartment
            
            for (i in names(list_reads)) {
                tab <- table(mcols(list_reads[[i]])$len_adj)
                tab_unq <- table(mcols(list_reads_unq[[i]])$len_adj)
                for (j in readlengths) {
                  cont <- as.numeric(tab[which(names(tab) == j)])
                  cont_unq <- as.numeric(tab_unq[which(names(tab_unq) == j)])
                  if (length(cont) > 0) {
                    list_vects[[i]][j - (min(readlengths) - 1)] <- cont
                  }
                  if (length(cont_unq) > 0) {
                    list_vects_unq[[i]][j - (min(readlengths) - 1)] <- cont_unq
                  }
                  
                }
            }
            
            rld_len <- do.call(what = rbind.data.frame, list_vects)
            names(rld_len) <- paste("reads", readlengths, sep = "_")
            row.names(rld_len) <- names(list_reads)
            
            rld_len_unq <- do.call(what = rbind.data.frame, list_vects_unq)
            names(rld_len_unq) <- paste("reads", readlengths, sep = "_")
            row.names(rld_len_unq) <- names(list_reads)
            
            
            # count reads per biotype location per compartment
            list_rld_loc <- list()
            for (i in c("nucl", circs)) {
                list_rld_loc[[i]] <- (assay(summarizeOverlaps(reads = x, features = GRangesList(list_locat[[i]]), 
                  ignore.strand = !as.logical(strandedness), mode = "Union", inter.feature = FALSE)))
            }
            rld_loc <- do.call(what = cbind.data.frame, list_rld_loc)
            names(rld_loc) <- paste("reads", names(list_rld_loc), sep = "_")
            
            
            list_rld_loc_unq <- list()
            for (i in c("nucl", circs)) {
                list_rld_loc_unq[[i]] <- (assay(summarizeOverlaps(reads = x_uniq, 
                  features = GRangesList(list_locat[[i]]), ignore.strand = !as.logical(strandedness), 
                  mode = "Union", inter.feature = FALSE)))
            }
            rld_loc_unq <- do.call(what = cbind.data.frame, list_rld_loc_unq)
            names(rld_loc_unq) <- paste("reads", names(list_rld_loc_unq), sep = "_")
            
            
            # reads_summary
            
            reads_summary <- list()
            
            # iterate over genome types and compartments
            for (i in names(list_reads)) {
                reads_by_length <- list()
                rdss <- list_reads[[i]]
                rdss <- split(rdss, mcols(rdss)$len_adj)
                
                rgs <- GRangesList(list_locat[[i]])
                ovs <- lapply(rdss, FUN = function(x) {
                  suppressWarnings(assay(summarizeOverlaps(reads = x, features = rgs, 
                    ignore.strand = FALSE, mode = "Union", inter.feature = FALSE)))
                })
                nope <- readlengths[which(!readlengths %in% names(ovs))]
                if (length(nope) > 0) {
                  for (j in nope) {
                    nop <- matrix(0, nrow = 7, ncol = 1)
                    rownames(nop) <- names(rgs)
                    ovs[[as.character(j)]] <- nop
                    
                  }
                }
                ovs <- ovs[names(ovs) %in% as.character(readlengths)]
                ovs <- ovs[as.character(readlengths)]
                ovs <- do.call(ovs, what = cbind)
                colnames(ovs) <- paste("reads", readlengths, sep = "_")
                reads_summary[[i]] <- DataFrame(ovs)
                
                
            }
            
            
            
            reads_summary_unq <- list()
            
            # iterate over genome types and compartments
            for (i in names(list_reads_unq)) {
                reads_by_length <- list()
                rdss <- list_reads_unq[[i]]
                rdss <- split(rdss, mcols(rdss)$len_adj)
                
                rgs <- GRangesList(list_locat[[i]])
                ovs <- lapply(rdss, FUN = function(x) {
                  suppressWarnings(assay(summarizeOverlaps(reads = x, features = rgs, 
                    ignore.strand = FALSE, mode = "Union", inter.feature = FALSE)))
                })
                nope <- readlengths[which(!readlengths %in% names(ovs))]
                if (length(nope) > 0) {
                  for (j in nope) {
                    nop <- matrix(0, nrow = 7, ncol = 1)
                    rownames(nop) <- names(rgs)
                    ovs[[as.character(j)]] <- nop
                    
                  }
                }
                ovs <- ovs[names(ovs) %in% as.character(readlengths)]
                ovs <- ovs[as.character(readlengths)]
                ovs <- do.call(ovs, what = cbind)
                colnames(ovs) <- paste("reads", readlengths, sep = "_")
                reads_summary_unq[[i]] <- DataFrame(ovs)
                
                
            }
            
            
            # reads_pos1_stst
            
            
            # take first position of each read
            reads_pos1 <- unlist(resize(split(GRanges(x), f = mcols(x)$len_adj), 
                1))
            
            reads_pos1 <- split(reads_pos1, reads_pos1$len_adj)
            
            reads_pos1 <- GRangesList(lapply(reads_pos1, function(y) {
                unq <- unique(y)
                mcols(unq) <- NULL
                unq$score <- countOverlaps(unq, y, type = "equal")
                unq
                
            }))
            
            # cnts_cds_genes, cnts_all_genes
            
            seqgene <- cbind(as.vector(seqnames(GTF_annotation$genes)), GTF_annotation$genes$gene_id)
            
            
            cnts_cds_genes <- assay(summarizeOverlaps(reads = x, features = red_cdss, 
                ignore.strand = !as.logical(strandedness), mode = "Union", inter.feature = FALSE))
            
            cnts_all_genes <- assay(summarizeOverlaps(reads = x, features = red_ex, 
                ignore.strand = !as.logical(strandedness), mode = "Union", inter.feature = FALSE))
            
            cnts_cds_genes_unq <- assay(summarizeOverlaps(reads = x_uniq, features = red_cdss, 
                ignore.strand = !as.logical(strandedness), mode = "Union", inter.feature = FALSE))
            
            cnts_all_genes_unq <- assay(summarizeOverlaps(reads = x_uniq, features = red_ex, 
                ignore.strand = !as.logical(strandedness), mode = "Union", inter.feature = FALSE))
            
            
            chrsss <- seqgene[match(rownames(cnts_cds_genes), seqgene[, 2]), 1]
            
            cnts_cds_genes <- DataFrame(cnts_cds_genes)
            cnts_cds_genes$chromosome <- chrsss
            cnts_cds_genes$gene_name <- GTF_annotation$trann$gene_name[match(rownames(cnts_cds_genes), 
                GTF_annotation$trann$gene_id)]
            cnts_cds_genes$gene_biotype <- GTF_annotation$trann$gene_biotype[match(rownames(cnts_cds_genes), 
                GTF_annotation$trann$gene_id)]
            
            cnts_all_genes <- DataFrame(cnts_all_genes)
            cnts_all_genes$chromosome <- seqgene[match(rownames(cnts_all_genes), 
                seqgene[, 2]), 1]
            cnts_all_genes$gene_name <- GTF_annotation$trann$gene_name[match(rownames(cnts_all_genes), 
                GTF_annotation$trann$gene_id)]
            cnts_all_genes$gene_biotype <- GTF_annotation$trann$gene_biotype[match(rownames(cnts_all_genes), 
                GTF_annotation$trann$gene_id)]
            
            
            cnts_cds_genes_unq <- DataFrame(cnts_cds_genes_unq)
            cnts_cds_genes_unq$chromosome <- chrsss
            cnts_cds_genes_unq$gene_name <- GTF_annotation$trann$gene_name[match(rownames(cnts_cds_genes_unq), 
                GTF_annotation$trann$gene_id)]
            cnts_cds_genes_unq$gene_biotype <- GTF_annotation$trann$gene_biotype[match(rownames(cnts_cds_genes_unq), 
                GTF_annotation$trann$gene_id)]
            
            cnts_all_genes_unq <- DataFrame(cnts_all_genes_unq)
            cnts_all_genes_unq$chromosome <- seqgene[match(rownames(cnts_all_genes_unq), 
                seqgene[, 2]), 1]
            cnts_all_genes_unq$gene_name <- GTF_annotation$trann$gene_name[match(rownames(cnts_all_genes_unq), 
                GTF_annotation$trann$gene_id)]
            cnts_all_genes_unq$gene_biotype <- GTF_annotation$trann$gene_biotype[match(rownames(cnts_all_genes_unq), 
                GTF_annotation$trann$gene_id)]
            
            
            
            chunk_res <- list(rld_len, rld_len_unq, rld_loc, rld_loc_unq, reads_pos1, 
                cnts_cds_genes, cnts_cds_genes_unq, cnts_all_genes, cnts_all_genes_unq, 
                reads_summary, reads_summary_unq)
            # export list to use in the reduc part above
            names(chunk_res) <- c("rld", "rld_unq", "positions", "positions_unq", 
                "reads_pos1", "counts_cds_genes", "counts_cds_genes_unq", "counts_all_genes", 
                "counts_all_genes_unq", "reads_summary", "reads_summary_unq")
            
            chunk_res
        }
        
        
        # main for reducebyYield
        
        cat(paste("Analyzing BAM file:", bam_file, "...", date(), "\n"))
        
        read_stats <- reduceByYield(X = opts, YIELD = yiel, MAP = mapp, REDUCE = reduc)
        names(read_stats) <- c("rld", "rld_unq", "positions", "positions_unq", "reads_pos1", 
            "counts_cds_genes", "counts_cds_genes_unq", "counts_all_genes", "counts_all_genes_unq", 
            "reads_summary", "reads_summary_unq")
        
        cds_cnts <- read_stats$counts_cds_genes
        cds_cnts$RPKM <- cds_cnts$reads/(sum(cds_cnts$reads)/1e+06)
        cds_cnts$RPKM <- round(cds_cnts$RPKM/(sum(width(red_cdss))/1000), digits = 4)
        cds_cnts$TPM <- cds_cnts$reads/(sum(width(red_cdss))/1000)
        cds_cnts$TPM <- round(cds_cnts$TPM/(sum(cds_cnts$TPM)/1e+06), digits = 4)
        read_stats$counts_cds_genes <- cds_cnts
        
        cds_cnts_unq <- read_stats$counts_cds_genes_unq
        cds_cnts_unq$RPKM <- cds_cnts_unq$reads/(sum(cds_cnts_unq$reads)/1e+06)
        cds_cnts_unq$RPKM <- round(cds_cnts_unq$RPKM/(sum(width(red_cdss))/1000), 
            digits = 4)
        cds_cnts_unq$TPM <- cds_cnts_unq$reads/(sum(width(red_cdss))/1000)
        cds_cnts_unq$TPM <- round(cds_cnts_unq$TPM/(sum(cds_cnts_unq$TPM)/1e+06), 
            digits = 4)
        read_stats$counts_cds_genes_unq <- cds_cnts_unq
        
        ex_cnts <- read_stats$counts_all_genes
        ex_cnts$RPKM <- ex_cnts$reads/(sum(ex_cnts$reads)/1e+06)
        ex_cnts$RPKM <- round(ex_cnts$RPKM/(sum(width(red_ex))/1000), digits = 4)
        ex_cnts$TPM <- ex_cnts$reads/(sum(width(red_ex))/1000)
        ex_cnts$TPM <- round(ex_cnts$TPM/(sum(ex_cnts$TPM)/1e+06), digits = 4)
        read_stats$counts_all_genes <- ex_cnts
        
        ex_cnts_unq <- read_stats$counts_all_genes_unq
        ex_cnts_unq$RPKM <- ex_cnts_unq$reads/(sum(ex_cnts_unq$reads)/1e+06)
        ex_cnts_unq$RPKM <- round(ex_cnts_unq$RPKM/(sum(width(red_ex))/1000), digits = 4)
        ex_cnts_unq$TPM <- ex_cnts_unq$reads/(sum(width(red_ex))/1000)
        ex_cnts_unq$TPM <- round(ex_cnts_unq$TPM/(sum(ex_cnts_unq$TPM)/1e+06), digits = 4)
        read_stats$counts_all_genes_unq <- ex_cnts_unq
        
        save(read_stats, file = paste(dira, "read_stats", sep = "/"))
        cat(paste("Analyzing BAM file --- Done!", date(), "\n"))
        
        ps <- read_stats$reads_pos1
        
        tile_cds <- GTF_annotation$cds_txs_coords
        tile_cds <- tile_cds[tile_cds$reprentative_mostcommon]
        
        list_txs_ok <- GRangesList()
        for (ci in circs) {
            citxs <- GTF_annotation$cds_txs[seqnames(GTF_annotation$cds_txs) == ci]
            citxs <- citxs[elementNROWS(citxs) > 0]
            citxs_txs <- GTF_annotation$cds_txs_coords[as.vector(seqnames(GTF_annotation$cds_txs_coords)) %in% 
                names(citxs)]
            list_txs_ok[[ci]] <- citxs_txs
        }
        list_txs_ok[["nucl"]] <- tile_cds[!tile_cds %in% unlist(list_txs_ok)]
        
        cat(paste("Building aggregate 5' profiles ...", date(), "\n"))
        ps_signals_win_all <- list()
        ps_signals_tiles_all <- list()
        
        # calculate profiles for bins and windows (single nt resolution)
        for (comp in c("nucl", circs)) {
            tile_cds <- list_txs_ok[[comp]]
            strand(tile_cds) <- "+"
            tile_cds <- tile_cds[width(tile_cds) > 100]
            checkgen <- TRUE
            if (comp == "nucl") {
                if ((sum((start(tile_cds) > 50))/length(start(tile_cds))) > 0.33) {
                  tile_cds <- tile_cds[start(tile_cds) > 50]
                  checkgen <- FALSE
                }
                if (sum(seqlengths(tile_cds)[as.vector(seqnames(tile_cds))] - end(tile_cds) > 
                  100)/length(tile_cds) > 0.33) {
                  tile_cds <- tile_cds[seqlengths(tile_cds)[as.vector(seqnames(tile_cds))] - 
                    end(tile_cds) > 100]
                  checkgen <- FALSE
                }
                ps_comp <- GRangesList(lapply(ps, function(x) {
                  sort(x[!seqnames(x) %in% circs])
                }))
            } else {
                ps_comp <- GRangesList(lapply(ps, function(x) {
                  sort(x[seqnames(x) == comp])
                }))
            }
            
            # select only some read lenghts, otherwise too much output and computation
            
            summm <- read_stats$reads_summary
            a <- unlist(summm[[comp]]["cds", ]) + unlist(summm[[comp]]["fiveutrs", 
                ]) + unlist(summm[[comp]]["threeutrs", ])
            names(a) <- gsub(names(a), pattern = "reads_", replacement = "")
            a <- sort(a/(sum(a)/100), decreasing = TRUE)
            cuma <- cumsum(a)
            if (length(cuma) == 0) {
                ps_signals_tiles_all[[comp]] <- c()
                ps_signals_win_all[[comp]] <- c()
                next
            }
            # if read subset, discard low abundance
            if (read_subset[bammo] == TRUE) {
                rl_ok <- names(cuma[1:which(cuma > 99)[1]])
                # rl_ok<-as.character(seq(sort(as.numeric(rl_ok),decreasing =
                # FALSE)[1],sort(as.numeric(rl_ok),decreasing = TRUE)[1]))
                ps_comp <- ps_comp[names(ps_comp) %in% rl_ok]
            }
            
            signal_ps <- list()
            signal_ps_nt <- list()
            
            all_ps <- unlist(ps_comp)
            ps_comp[["all"]] <- all_ps
            mp <- mapToTranscripts(all_ps, transcripts = GTF_annotation$exons_txs[as.vector(seqnames(tile_cds))])
            mp$score <- all_ps$score[mp$xHits]
            # put seqlevels for compartment only to speed up stuff
            seqlevels(mp) <- as.vector(seqnames(tile_cds))
            seqlengths(mp) <- sum(width(GTF_annotation$exons_txs[as.vector(seqnames(tile_cds))]))
            covtx <- coverage(mp, weight = mp$score)
            ok_txs <- names(covtx[elementNROWS(runValue(covtx)) > 1])
            if (length(ok_txs) < 5) {
                ps_signals_tiles_all[[comp]] <- c()
                ps_signals_win_all[[comp]] <- c()
                next
            }
            # keep all txs if small n of txs
            
            if (length(ok_txs) >= 5 & length(covtx) < 250) {
                ok_txs <- names(covtx)
            }
            
            if (length(ok_txs) > 5000) {
                cnts_txss_agg <- sort(sum(covtx), decreasing = TRUE)
                ok_txs <- names(cnts_txss_agg)[1:5000]
            }
            if (fast_mode[bammo] == TRUE) {
                if (length(ok_txs) > 500) {
                  cnts_txss_agg <- sort(sum(covtx), decreasing = TRUE)
                  ok_txs <- names(cnts_txss_agg)[1:500]
                }
            }
            
            
            
            fivs <- GRanges(seqnames = seqnames(GTF_annotation$cds_txs_coords), ranges = IRanges(start = 1, 
                end = start(GTF_annotation$cds_txs_coords)), strand = "*")
            fivs <- fivs[as.character(seqnames(fivs)) %in% ok_txs]
            fivs_gen <- unlist(pmapFromTranscripts(fivs, transcripts = GTF_annotation$exons_txs[as.character(seqnames(fivs))], 
                ignore.strand = FALSE))
            fivs_gen <- fivs_gen[fivs_gen$hit]
            fivs_gen <- split(fivs_gen, names(fivs_gen))
            
            
            threes <- GRanges(seqnames = seqnames(GTF_annotation$cds_txs_coords), 
                ranges = IRanges(start = end(GTF_annotation$cds_txs_coords), end = GTF_annotation$cds_txs_coords$lentx), 
                strand = "*")
            threes <- threes[as.character(seqnames(threes)) %in% ok_txs]
            threes_gen <- unlist(pmapFromTranscripts(threes, transcripts = GTF_annotation$exons_txs[as.character(seqnames(threes))], 
                ignore.strand = FALSE))
            threes_gen <- threes_gen[threes_gen$hit]
            threes_gen <- split(threes_gen, names(threes_gen))
            
            cds_gen <- GTF_annotation$cds_txs[ok_txs]
            
            
            tile_cds <- tile_cds[as.character(seqnames(tile_cds)) %in% ok_txs]
            ok_txs <- unique(as.character(seqnames(tile_cds)))
            seqlevels(tile_cds) <- ok_txs
            list_covs <- list()
            list_covs[["all"]] <- covtx
            
            tile_5 <- GRanges(seqnames(tile_cds), ranges = IRanges(start = 1, end = start(tile_cds)))
            tile_3 <- GRanges(seqnames(tile_cds), ranges = IRanges(start = end(tile_cds), 
                end = seqlengths(tile_cds)[as.vector(seqnames(tile_cds))]))
            no5utr <- width(tile_5) < 51
            no3utr <- width(tile_3) < 51
            
            ex_annot <- GTF_annotation$exons_txs
            # if no utrs put bins to 1nt
            
            if (sum(no5utr) > 0) {
                tx_notok <- seqnames(tile_5)[no5utr]
                annot_notok <- ex_annot[tx_notok]
                annot_ok <- GRangesList(lapply(annot_notok, function(x) {
                  if (length(x) == 0) {
                    return(x)
                  }
                  x[1] <- resize(x[1], width = width(x[1]) + 51, fix = "end")
                  x
                }))
                ex_annot[names(annot_ok)] <- annot_ok
                
                seqlengths(tile_cds)[as.vector(tx_notok)] <- sum(width(annot_ok))
                tile_cds[no5utr] <- shift(tile_cds[no5utr], +51)
                tile_5 <- GRanges(seqnames(tile_cds), ranges = IRanges(start = 1, 
                  end = start(tile_cds)))
                tile_3 <- GRanges(seqnames(tile_cds), ranges = IRanges(start = end(tile_cds), 
                  end = seqlengths(tile_cds)[as.vector(seqnames(tile_cds))]))
                
                mp <- mapToTranscripts(all_ps, transcripts = ex_annot[as.vector(seqnames(tile_cds))])
                seqlevels(mp) <- seqlevels(tile_cds)
                seqlengths(mp) <- seqlengths(tile_cds)
                mp$score <- all_ps$score[mp$xHits]
                
                covtx <- coverage(mp, weight = mp$score)
                list_covs[["all"]] <- covtx
            }
            
            if (sum(no3utr) > 0) {
                tx_notok <- seqnames(tile_3)[no3utr]
                annot_notok <- ex_annot[tx_notok]
                annot_ok <- GRangesList(lapply(annot_notok, function(x) {
                  if (length(x) == 0) {
                    return(x)
                  }
                  x[length(x)] <- resize(x[length(x)], width = width(x[length(x)]) + 
                    51, fix = "start")
                  x
                }))
                ex_annot[names(annot_ok)] <- annot_ok
                
                seqlengths(tile_cds)[as.vector(tx_notok)] <- sum(width(annot_ok))
                
                tile_5 <- GRanges(seqnames(tile_cds), ranges = IRanges(start = 1, 
                  end = start(tile_cds)))
                tile_3 <- GRanges(seqnames(tile_cds), ranges = IRanges(start = end(tile_cds), 
                  end = seqlengths(tile_cds)[as.vector(seqnames(tile_cds))]))
                
                mp <- mapToTranscripts(all_ps, transcripts = ex_annot[as.vector(seqnames(tile_cds))])
                seqlevels(mp) <- seqlevels(tile_cds)
                seqlengths(mp) <- seqlengths(tile_cds)
                mp$score <- all_ps$score[mp$xHits]
                
                covtx <- coverage(mp, weight = mp$score)
                list_covs[["all"]] <- covtx
                
            }
            
            ps_tiles <- DataFrameList()
            ps_win <- DataFrameList()
            
            for (len in c("all", names(ps_comp))) {
                
                if ((sum(no3utr) + sum(no5utr)) == 0) {
                  mp <- mapToTranscripts(ps_comp[[len]], transcripts = fivs_gen)
                  mp$score <- ps_comp[[len]]$score[mp$xHits]
                  seqlevels(mp) <- names(fivs_gen)
                  seqlengths(mp) <- sum(width(fivs_gen))
                  cov_5 <- coverage(mp, weight = mp$score)
                  
                  
                  mp <- mapToTranscripts(ps_comp[[len]], transcripts = threes_gen)
                  mp$score <- ps_comp[[len]]$score[mp$xHits]
                  seqlevels(mp) <- names(threes_gen)
                  seqlengths(mp) <- sum(width(threes_gen))
                  cov_3 <- coverage(mp, weight = mp$score)
                  
                  mp <- mapToTranscripts(ps_comp[[len]], transcripts = cds_gen)
                  mp$score <- ps_comp[[len]]$score[mp$xHits]
                  seqlevels(mp) <- names(cds_gen)
                  seqlengths(mp) <- sum(width(cds_gen))
                  cov_cds <- coverage(mp, weight = mp$score)
                }
                
                if ((sum(no3utr) + sum(no5utr)) > 0) {
                  if (len != "all") {
                    mp <- mapToTranscripts(ps_comp[[len]], transcripts = ex_annot[as.vector(seqnames(tile_cds))])
                    mp$score <- ps_comp[[len]]$score[mp$xHits]
                    seqlevels(mp) <- seqlevels(tile_cds)
                    seqlengths(mp) <- seqlengths(tile_cds)
                    covtx <- coverage(mp, weight = mp$score)
                    covtx <- covtx[ok_txs]
                  }
                  
                  cov_5 <- covtx[tile_5]
                  cov_3 <- covtx[tile_3]
                  cov_cds <- covtx[tile_cds]
                  
                  
                }
                
                ps_tiles_5 <- DataFrame(t(sapply(cov_5, function(x) {
                  clos <- 50 * (round(length(x)/50, digits = 0) + 1)
                  idx <- as.integer(seq(1, length(x), length.out = clos))
                  colMeans(matrix(x[idx], ncol = 50))
                })))
                ps_win_5 <- DataFrame(t(sapply(cov_5, function(x) {
                  as.vector(x)[c(1:25, (length(x) - 25):(length(x) - 1))]
                })))
                
                ps_tiles_3 <- DataFrame(t(sapply(cov_3, function(x) {
                  clos <- 50 * (round(length(x)/50, digits = 0) + 1)
                  idx <- as.integer(seq(1, length(x), length.out = clos))
                  colMeans(matrix(x[idx], ncol = 50))
                })))
                ps_win_3 <- DataFrame(t(sapply(cov_3, function(x) {
                  as.vector(x)[c(2:26, (length(x) - 24):length(x))]
                })))
                
                
                ps_tiles_cds <- DataFrame(t(sapply(cov_cds, function(x) {
                  clos <- 100 * (round(length(x)/100, digits = 0) + 1)
                  idx <- as.integer(seq(1, length(x), length.out = clos))
                  colMeans(matrix(x[idx], ncol = 100))
                })))
                ps_win_cds <- DataFrame(t(sapply(cov_cds, function(x) {
                  rnd <- as.integer(length(x)/2)%%3
                  mid <- (as.integer(length(x)/2) - rnd)
                  as.vector(x)[c(1:33, (mid - 17):(mid + 15), (length(x) - 32):length(x))]
                })))
                tls <- cbind(ps_tiles_5, ps_tiles_cds, ps_tiles_3)
                colnames(tls) <- c(paste("5_UTR", 1:length(ps_tiles_5[1, ]), sep = "_"), 
                  paste("CDS", 1:length(ps_tiles_cds[1, ]), sep = "_"), paste("3_UTR", 
                    1:length(ps_tiles_3[1, ]), sep = "_"))
                nts <- cbind(ps_win_5, ps_win_cds, ps_win_3)
                colnames(nts) <- c(paste("5_UTR", 1:length(ps_win_5[1, ]), sep = "_"), 
                  paste("CDS", 1:length(ps_win_cds[1, ]), sep = "_"), paste("3_UTR", 
                    1:length(ps_win_3[1, ]), sep = "_"))
                # DataFrameList?
                ps_tiles[[len]] <- tls
                ps_win[[len]] <- nts
                
            }
            
            
            ps_signals_tiles_all[[comp]] <- ps_tiles
            ps_signals_win_all[[comp]] <- ps_win
        }
        
        profiles_fivepr <- list(ps_signals_tiles_all, ps_signals_win_all)
        names(profiles_fivepr) <- c("five_prime_bins", "five_prime_subcodon")
        save(profiles_fivepr, file = paste(dira, "profiles_fivepr", sep = "/"))
        
        cat(paste("Building aggregate 5' profiles --- Done!", date(), "\n"))
        
        cat(paste("Calculating P-sites offsets ...", date(), "\n"))
        
        
        # calculate cutoffs for each read length
        
        list_cutoff <- lapply(profiles_fivepr[["five_prime_subcodon"]], function(x) {
            list_res <- list()
            for (i in names(x)) {
                if (i == "all") {
                  lnmax = 100
                } else {
                  lnmax = as.numeric(i)
                }
                list_res[[i]] <- calc_cutoffs_from_profiles(x[[i]], length_max = lnmax)
                
            }
            list_res
        })
        
        summary_res <- lapply(list_cutoff, function(x) t(sapply(x, function(y) {
            res <- c(y$final_cutoff, y$frames_res[2])
            names(res) <- c("cutoff", "frame_preference")
            res
        })))
        
        
        # choose readlengths
        
        res_rls <- choose_readlengths(summary_res, choice = readlength_choice_method[bammo], 
            nt_signals = profiles_fivepr[["five_prime_subcodon"]])
        
        cat(paste("Calculating P-sites offsets --- Done!", date(), "\n"))
        selection_cutoffs <- list(res_rls, summary_res, list_cutoff)
        names(selection_cutoffs) <- c("results_choice", "results_cutoffs", "analysis_frame_cutoff")
        save(selection_cutoffs, file = paste(dira, "selection_cutoffs", sep = "/"))
        
        
        
        # extract P-sites positions given selection above
        
        param <- ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE, isSecondaryAlignment = FALSE), 
            what = c("mapq"), tag = "MD")
        seqllll <- seqlevels(unlist(unlist(read_stats$reads_pos1)))
        seqleee <- seqlengths(unlist(unlist(read_stats$reads_pos1)))
        
        # if we've provided our own choice of offsets, we should repeat it for each
        # compartment, should we need to
        if (!is.null(offsets_df)) {
            if (is.data.frame(offsets_df)) {
                res_rls$nucl$final_choice <- offsets_df %>% mutate(read_length = as.character(read_length)) %>% 
                  DataFrame
            } else {
                stopifnot(names(res_rls) %in% names(offsets_df))
                for (comp in names(res_rls)) res_rls[[comp]]$final_choice <- offsets_df[[comp]]
            }
        }
        
        rl_cutoffs_comp <- lapply(res_rls, function(x) {
            x$final_choice
        })
        
        # rescue all read lengths
        if (rescue_all_rls[bammo] == TRUE) {
            all_rlsss <- as.numeric(names(read_stats$reads_pos1))
            for (rilo in c("nucl", circs)) {
                ralo <- rl_cutoffs_comp[[rilo]]
                if (is.null(ralo)) {
                  ralo <- DataFrame()
                }
                rlno2 <- all_rlsss[!all_rlsss %in% ralo$read_length]
                ralo2 <- DataFrame(read_length = rlno2, cutoff = 12)
                rl_cutoffs_comp[[rilo]] <- rbind(ralo, ralo2)
            }
        }
        
        
        # what to do with chunks: x present chunk, y old chunks (cumulative)
        reduc <- function(x, y) {
            # adjust merging by rls
            all_ps <- GRangesList()
            rls <- unique(c(names(x[["P_sites_all"]]), names(y[["P_sites_all"]])))
            
            for (rl in rls) {
                reads_x <- GRanges()
                reads_y <- GRanges()
                
                seqlevels(reads_x) <- seqllll
                seqlevels(reads_y) <- seqllll
                
                seqlengths(reads_x) <- seqleee
                seqlengths(reads_y) <- seqleee
                
                if (sum(rl %in% names(x[["P_sites_all"]])) > 0) {
                  reads_x <- x[["P_sites_all"]][[rl]]
                }
                if (sum(rl %in% names(y[["P_sites_all"]])) > 0) {
                  reads_y <- y[["P_sites_all"]][[rl]]
                }
                
                plx <- reads_x[strand(reads_x) == "+"]
                mnx <- reads_x[strand(reads_x) == "-"]
                ply <- reads_y[strand(reads_y) == "+"]
                mny <- reads_y[strand(reads_y) == "-"]
                if (length(plx) > 0) {
                  covv_pl <- coverage(plx, weight = plx$score)
                } else {
                  covv_pl <- coverage(plx)
                }
                if (length(ply) > 0) {
                  covv_pl <- covv_pl + coverage(ply, weight = ply$score)
                }
                
                covv_pl <- GRanges(covv_pl)
                covv_pl <- covv_pl[covv_pl$score > 0]
                
                if (length(mnx) > 0) {
                  covv_min <- coverage(mnx, weight = mnx$score)
                } else {
                  covv_min <- coverage(mnx)
                }
                if (length(mny) > 0) {
                  covv_min <- covv_min + coverage(mny, weight = mny$score)
                }
                
                covv_min <- GRanges(covv_min)
                covv_min <- covv_min[covv_min$score > 0]
                
                strand(covv_pl) <- "+"
                strand(covv_min) <- "-"
                
                all_ps[[rl]] <- sort(c(covv_pl, covv_min))
                rm(covv_min, covv_pl, plx, mnx, ply, mny, reads_y, reads_x, rl)
                
            }
            
            uniq_ps <- GRangesList()
            rls <- unique(c(names(x[["P_sites_uniq"]]), names(y[["P_sites_uniq"]])))
            
            for (rl in rls) {
                reads_x <- GRanges()
                reads_y <- GRanges()
                
                seqlevels(reads_x) <- seqllll
                seqlevels(reads_y) <- seqllll
                
                seqlengths(reads_x) <- seqleee
                seqlengths(reads_y) <- seqleee
                if (sum(rl %in% names(x[["P_sites_uniq"]])) > 0) {
                  reads_x <- x[["P_sites_uniq"]][[rl]]
                }
                if (sum(rl %in% names(y[["P_sites_uniq"]])) > 0) {
                  reads_y <- y[["P_sites_uniq"]][[rl]]
                }
                
                plx <- reads_x[strand(reads_x) == "+"]
                mnx <- reads_x[strand(reads_x) == "-"]
                ply <- reads_y[strand(reads_y) == "+"]
                mny <- reads_y[strand(reads_y) == "-"]
                if (length(plx) > 0) {
                  covv_pl <- coverage(plx, weight = plx$score)
                } else {
                  covv_pl <- coverage(plx)
                }
                if (length(ply) > 0) {
                  covv_pl <- covv_pl + coverage(ply, weight = ply$score)
                }
                
                covv_pl <- GRanges(covv_pl)
                covv_pl <- covv_pl[covv_pl$score > 0]
                
                if (length(mnx) > 0) {
                  covv_min <- coverage(mnx, weight = mnx$score)
                } else {
                  covv_min <- coverage(mnx)
                }
                if (length(mny) > 0) {
                  covv_min <- covv_min + coverage(mny, weight = mny$score)
                }
                
                covv_min <- GRanges(covv_min)
                covv_min <- covv_min[covv_min$score > 0]
                
                strand(covv_pl) <- "+"
                strand(covv_min) <- "-"
                
                uniq_ps[[rl]] <- sort(c(covv_pl, covv_min))
                rm(covv_min, covv_pl, plx, mnx, ply, mny, reads_y, reads_x, rl)
                
                
            }
            
            uniq_mm_ps <- GRangesList()
            rls <- unique(c(names(x[["P_sites_uniq_mm"]]), names(y[["P_sites_uniq_mm"]])))
            
            for (rl in rls) {
                reads_x <- GRanges()
                reads_y <- GRanges()
                
                seqlevels(reads_x) <- seqllll
                seqlevels(reads_y) <- seqllll
                
                seqlengths(reads_x) <- seqleee
                seqlengths(reads_y) <- seqleee
                if (sum(rl %in% names(x[["P_sites_uniq_mm"]])) > 0) {
                  reads_x <- x[["P_sites_uniq_mm"]][[rl]]
                }
                if (sum(rl %in% names(y[["P_sites_uniq_mm"]])) > 0) {
                  reads_y <- y[["P_sites_uniq_mm"]][[rl]]
                }
                
                plx <- reads_x[strand(reads_x) == "+"]
                mnx <- reads_x[strand(reads_x) == "-"]
                ply <- reads_y[strand(reads_y) == "+"]
                mny <- reads_y[strand(reads_y) == "-"]
                if (length(plx) > 0) {
                  covv_pl <- coverage(plx, weight = plx$score)
                } else {
                  covv_pl <- coverage(plx)
                }
                if (length(ply) > 0) {
                  covv_pl <- covv_pl + coverage(ply, weight = ply$score)
                }
                
                covv_pl <- GRanges(covv_pl)
                covv_pl <- covv_pl[covv_pl$score > 0]
                
                if (length(mnx) > 0) {
                  covv_min <- coverage(mnx, weight = mnx$score)
                } else {
                  covv_min <- coverage(mnx)
                }
                if (length(mny) > 0) {
                  covv_min <- covv_min + coverage(mny, weight = mny$score)
                }
                
                covv_min <- GRanges(covv_min)
                covv_min <- covv_min[covv_min$score > 0]
                
                strand(covv_pl) <- "+"
                strand(covv_min) <- "-"
                
                uniq_mm_ps[[rl]] <- sort(c(covv_pl, covv_min))
                rm(covv_min, covv_pl, plx, mnx, ply, mny, reads_y, reads_x, rl)
                
            }
            
            covall_plus <- x$coverage_all_plus + y$coverage_all_plus
            covall_min <- x$coverage_all_min + y$coverage_all_min
            
            covuni_plus <- x$coverage_uniq_plus + y$coverage_uniq_plus
            covuni_min <- x$coverage_uniq_min + y$coverage_uniq_min
            
            
            rang_jun <- x$junctions
            rang_jun$reads <- rang_jun$reads + y$junctions$reads
            rang_jun$unique_reads <- rang_jun$unique_reads + y$junctions$unique_reads
            
            list_res <- list(all_ps, uniq_ps, uniq_mm_ps, rang_jun, covall_plus, 
                covall_min, covuni_plus, covuni_min)
            names(list_res) <- c("P_sites_all", "P_sites_uniq", "P_sites_uniq_mm", 
                "junctions", "coverage_all_plus", "coverage_all_min", "coverage_uniq_plus", 
                "coverage_uniq_min")
            
            
            return(list_res)
        }
        
        # what to do with each chunk (read as alignment file)
        
        yiel <- function(x) {
            readGAlignments(x, param = param)
        }
        
        # operations on the chunk (here count reads and whatnot)
        
        mapp <- function(x) {
            x_I <- x[grep("I", cigar(x))]
            
            if (length(x_I) > 0) {
                x <- x[grep("I", cigar(x), invert = TRUE)]
                
            }
            x_D <- x[grep("D", cigar(x))]
            if (length(x_D) > 0) {
                x <- x[grep("D", cigar(x), invert = TRUE)]
                
            }
            
            
            # softclipping
            
            
            clipp <- width(cigarRangesAlongQuerySpace(x@cigar, ops = "S"))
            clipp[elementNROWS(clipp) == 0] <- 0
            len_adj <- qwidth(x) - sum(clipp)
            mcols(x)$len_adj <- len_adj
            
            # Remove S from Cigar (read positions/length are already adjusted) it helps
            # calculating P-sites positions for spliced reads
            
            cigg <- cigar(x)
            cigg_s <- grep(cigg, pattern = "S")
            if (length(cigg_s) > 0) {
                cigs <- cigg[cigg_s]
                cigs <- gsub(cigs, pattern = "^[0-9]+S", replacement = "")
                cigs <- gsub(cigs, pattern = "[0-9]+S$", replacement = "")
                cigg[cigg_s] <- cigs
                x@cigar <- cigg
            }
            mcols(x)$cigar_str <- x@cigar
            x_uniq <- x[x@elementMetadata$mapq > 50]
            
            pos <- x[strand(x) == "+"]
            neg <- x[strand(x) == "-"]
            
            uniq_pos <- x_uniq[strand(x_uniq) == "+"]
            uniq_neg <- x_uniq[strand(x_uniq) == "-"]
            
            
            # coverage
            
            
            covuni_plus <- coverage(uniq_pos)
            covuni_min <- coverage(uniq_neg)
            covall_plus <- coverage(pos)
            covall_min <- coverage(neg)
            
            
            # junctions
            
            
            juns <- summarizeJunctions(x)
            juns_pos <- juns
            juns_neg <- juns
            mcols(juns_pos) <- NULL
            mcols(juns_neg) <- NULL
            juns_pos$reads <- juns$plus_score
            juns_neg$reads <- juns$minus_score
            strand(juns_pos) <- "+"
            strand(juns_neg) <- "-"
            juns <- sort(c(juns_pos, juns_neg))
            juns <- juns[juns$reads > 0]
            
            uniq_juns <- summarizeJunctions(x_uniq)
            uniq_juns_pos <- uniq_juns
            uniq_juns_neg <- uniq_juns
            mcols(uniq_juns_pos) <- NULL
            mcols(uniq_juns_neg) <- NULL
            uniq_juns_pos$reads <- uniq_juns$plus_score
            uniq_juns_neg$reads <- uniq_juns$minus_score
            strand(uniq_juns_pos) <- "+"
            strand(uniq_juns_neg) <- "-"
            uniq_juns <- sort(c(uniq_juns_pos, uniq_juns_neg))
            uniq_juns <- uniq_juns[uniq_juns$reads > 0]
            if (length(juns) > 0) {
                juns$unique_reads <- 0
                mat <- match(uniq_juns, juns)
                juns$unique_reads[mat] <- uniq_juns$reads
            }
            rang_jun <- GTF_annotation$junctions
            rang_jun$reads <- 0
            rang_jun$unique_reads <- 0
            if (length(juns) > 0) {
                mat <- match(juns, rang_jun)
                juns <- juns[!is.na(mat)]
                mat <- mat[!is.na(mat)]
                rang_jun$reads[mat] <- juns$reads
                rang_jun$unique_reads[mat] <- juns$unique_reads
            }
            
            
            
            # P-sites calculation
            
            
            list_pss <- list()
            for (comp in names(rl_cutoffs_comp)) {
                all_rl_ps <- GRangesList()
                uniq_rl_ps <- GRangesList()
                uniq_rl_mm_ps <- GRangesList()
                
                seqlevels(all_rl_ps) <- seqllll
                seqlevels(uniq_rl_ps) <- seqllll
                seqlevels(uniq_rl_mm_ps) <- seqllll
                
                seqlengths(all_rl_ps) <- seqleee
                seqlengths(uniq_rl_ps) <- seqleee
                seqlengths(uniq_rl_mm_ps) <- seqleee
                
                
                chroms <- comp
                
                if (comp == "nucl") {
                  chroms = seqlevels(x)[!seqlevels(x) %in% circs]
                }
                resul <- rl_cutoffs_comp[[comp]]
                
                for (i in seq_along(resul$read_length)) {
                  
                  all_ps <- GRangesList()
                  uniq_ps <- GRangesList()
                  uniq_mm_ps <- GRangesList()
                  seqlevels(all_ps) <- seqllll
                  seqlevels(uniq_ps) <- seqllll
                  seqlevels(uniq_mm_ps) <- seqllll
                  
                  seqlengths(all_ps) <- seqleee
                  seqlengths(uniq_ps) <- seqleee
                  seqlengths(uniq_mm_ps) <- seqleee
                  
                  rl <- as.numeric(resul$read_length[i])
                  ct <- as.numeric(resul$cutoff[i])
                  ok_reads <- pos[mcols(pos)$len_adj %in% rl]
                  ok_reads <- ok_reads[as.vector(seqnames(ok_reads)) %in% chroms]
                  
                  ps_plus <- GRanges()
                  seqlevels(ps_plus) <- seqllll
                  seqlengths(ps_plus) <- seqleee
                  ps_plus_uniq <- ps_plus
                  ps_plus_uniq_mm <- ps_plus
                  
                  if (length(ok_reads) > 0) {
                    unspl <- ok_reads[grep(pattern = "N", x = cigar(ok_reads), invert = TRUE)]
                    
                    ps_unspl <- shift(resize(GRanges(unspl), width = 1, fix = "start"), 
                      shift = ct)
                    
                    spl <- ok_reads[grep(pattern = "N", x = cigar(ok_reads))]
                    firstb <- as.numeric(sapply(strsplit(cigar(spl), "M"), "[[", 
                      1))
                    lastb <- as.numeric(sapply(strsplit(cigar(spl), "M"), function(x) {
                      gsub(x[length(x)], pattern = "^[^_]*N", replacement = "")
                    }))
                    firstok <- spl[firstb > ct]
                    firstok <- shift(resize(GRanges(firstok), width = 1, fix = "start"), 
                      shift = ct)
                    
                    lastok <- spl[lastb >= rl - ct]
                    lastok <- shift(resize(GRanges(lastok), width = 1, fix = "end"), 
                      shift = -(rl - ct - 1))
                    
                    
                    multi <- spl[firstb <= ct & lastb < rl - ct]
                    
                    
                    ps_spl <- GRanges()
                    seqlevels(ps_spl) <- seqllll
                    seqlengths(ps_spl) <- seqleee
                    
                    if (length(multi) > 0) {
                      ps_spl <- get_ps_fromspliceplus(multi, cutoff = ct)
                      
                    }
                    mcols(ps_spl) <- mcols(multi)
                    
                    seqlevels(firstok) <- seqllll
                    seqlevels(lastok) <- seqllll
                    seqlevels(ps_unspl) <- seqllll
                    seqlevels(ps_spl) <- seqllll
                    
                    seqlengths(firstok) <- seqleee
                    seqlengths(lastok) <- seqleee
                    seqlengths(ps_unspl) <- seqleee
                    seqlengths(ps_spl) <- seqleee
                    
                    ps_plus <- c(ps_unspl, firstok, lastok, ps_spl)
                    
                    ps_plus_uniq <- ps_plus[mcols(ps_plus)$mapq > 50]
                    
                    mapqokay <- mcols(ok_reads)$mapq > 50 & grepl(x = mcols(ok_reads)$MD, 
                      pattern = "\\W|.{3,}")
                    length(mapqokay)
                    ps_plus_uniq_mm <- ps_plus[mapqokay]
                    
                    mcols(ps_plus) <- NULL
                    mcols(ps_plus_uniq) <- NULL
                    mcols(ps_plus_uniq_mm) <- NULL
                    
                  }
                  ok_reads <- neg[mcols(neg)$len_adj %in% rl]
                  ok_reads <- ok_reads[as.vector(seqnames(ok_reads)) %in% chroms]
                  
                  ps_neg <- GRanges()
                  seqlevels(ps_neg) <- seqllll
                  seqlengths(ps_neg) <- seqleee
                  ps_neg_uniq <- ps_neg
                  ps_neg_uniq_mm <- ps_neg
                  
                  if (length(ok_reads) > 0) {
                    unspl <- ok_reads[grep(pattern = "N", x = cigar(ok_reads), invert = TRUE)]
                    
                    ps_unspl <- shift(resize(GRanges(unspl), width = 1, fix = "start"), 
                      shift = -ct)
                    
                    spl <- ok_reads[grep(pattern = "N", x = cigar(ok_reads))]
                    
                    firstb <- as.numeric(sapply(strsplit(cigar(spl), "M"), "[[", 
                      1))
                    lastb <- as.numeric(sapply(strsplit(cigar(spl), "M"), function(x) {
                      gsub(x[length(x)], pattern = "^[^_]*N", replacement = "")
                    }))
                    lastok <- spl[lastb > ct]
                    lastok <- shift(resize(GRanges(lastok), width = 1, fix = "start"), 
                      shift = -ct)
                    
                    firstok <- spl[firstb >= rl - ct]
                    firstok <- shift(resize(GRanges(firstok), width = 1, fix = "end"), 
                      shift = (rl - ct - 1))
                    
                    multi <- spl[firstb < rl - ct & lastb <= ct]
                    
                    
                    ps_spl <- GRanges()
                    seqlevels(ps_spl) <- seqllll
                    seqlengths(ps_spl) <- seqleee
                    
                    
                    if (length(multi) > 0) {
                      ps_spl <- get_ps_fromsplicemin(multi, cutoff = ct)
                    }
                    mcols(ps_spl) <- mcols(multi)
                    
                    seqlevels(firstok) <- seqllll
                    seqlevels(lastok) <- seqllll
                    seqlevels(ps_unspl) <- seqllll
                    seqlevels(ps_spl) <- seqllll
                    
                    seqlengths(firstok) <- seqleee
                    seqlengths(lastok) <- seqleee
                    seqlengths(ps_unspl) <- seqleee
                    seqlengths(ps_spl) <- seqleee
                    
                    ps_neg <- c(ps_unspl, firstok, lastok, ps_spl)
                    # HERE PROBLEM WITH MAPPING QUALITY
                    ps_neg_uniq <- ps_neg[mcols(ps_neg)$mapq > 50]
                    mapqokay <- mcols(ok_reads)$mapq > 50 & grepl(x = mcols(ok_reads)$MD, 
                      pattern = "\\W|.{3,}")
                    length(mapqokay)
                    
                    ps_neg_uniq_mm <- ok_reads[mapqokay]
                    
                    mcols(ps_neg) <- NULL
                    mcols(ps_neg_uniq) <- NULL
                    mcols(ps_neg_uniq_mm) <- NULL
                    
                  }
                  
                  all_ps <- sort(c(ps_plus, ps_neg))
                  uniq_ps <- sort(c(ps_plus_uniq, ps_neg_uniq))
                  uniq_mm_ps <- sort(c(ps_plus_uniq_mm, ps_neg_uniq_mm))
                  if (length(all_ps) > 0) {
                    ps_res <- unique(all_ps)
                    ps_res$score <- countOverlaps(ps_res, all_ps, type = "equal")
                    all_ps <- ps_res
                  }
                  if (length(uniq_ps) > 0) {
                    ps_res <- unique(uniq_ps)
                    ps_res$score <- countOverlaps(ps_res, uniq_ps, type = "equal")
                    uniq_ps <- ps_res
                    
                  }
                  if (length(uniq_mm_ps) > 0) {
                    ps_res <- unique(uniq_mm_ps)
                    ps_res$score <- countOverlaps(ps_res, uniq_mm_ps, type = "equal")
                    uniq_mm_ps <- ps_res
                    
                  }
                  all_rl_ps[[as.character(rl)]] <- all_ps
                  uniq_rl_ps[[as.character(rl)]] <- uniq_ps
                  uniq_rl_mm_ps[[as.character(rl)]] <- uniq_mm_ps
                  
                  
                }
                # here comps
                list_rlct <- list(all_rl_ps, uniq_rl_ps, uniq_rl_mm_ps)
                names(list_rlct) <- c("P_sites_all", "P_sites_uniq", "P_sites_uniq_mm")
                list_pss[[comp]] <- list_rlct
            }
            # for rl, merge psites
            
            
            all_ps_comps <- GRangesList()
            seqlevels(all_ps_comps) <- seqllll
            seqlengths(all_ps_comps) <- seqleee
            rls_comps <- unique(unlist(lapply(list_pss, FUN = function(x) names(x[["P_sites_all"]]))))
            for (rl in rls_comps) {
                reads_rl_comp <- GRanges()
                seqlevels(reads_rl_comp) <- seqllll
                seqlengths(reads_rl_comp) <- seqleee
                for (comp in names(list_pss)) {
                  if (sum(rl %in% names(list_pss[[comp]][["P_sites_all"]])) > 0) {
                    oth <- list_pss[[comp]][["P_sites_all"]][[rl]]
                    
                    if (!is.null(oth)) {
                      seqlevels(oth) <- seqllll
                      seqlengths(oth) <- seqleee
                      reads_rl_comp <- c(reads_rl_comp, oth)
                    }
                  }
                  all_ps_comps[[rl]] <- reads_rl_comp
                }
                
            }
            
            uniq_ps_comps <- GRangesList()
            seqlevels(uniq_ps_comps) <- seqllll
            seqlengths(uniq_ps_comps) <- seqleee
            rls_comps <- unique(unlist(lapply(list_pss, FUN = function(x) names(x[["P_sites_uniq"]]))))
            for (rl in rls_comps) {
                reads_rl_comp <- GRanges()
                seqlevels(reads_rl_comp) <- seqllll
                seqlengths(reads_rl_comp) <- seqleee
                for (comp in names(list_pss)) {
                  if (sum(rl %in% names(list_pss[[comp]][["P_sites_uniq"]])) > 0) {
                    oth <- list_pss[[comp]][["P_sites_uniq"]][[rl]]
                    
                    if (!is.null(oth)) {
                      seqlevels(oth) <- seqllll
                      seqlengths(oth) <- seqleee
                      reads_rl_comp <- c(reads_rl_comp, oth)
                    }
                  }
                  uniq_ps_comps[[rl]] <- reads_rl_comp
                }
                
            }
            
            uniq_mm_ps_comps <- GRangesList()
            seqlevels(uniq_mm_ps_comps) <- seqllll
            seqlengths(uniq_mm_ps_comps) <- seqleee
            rls_comps <- unique(unlist(lapply(list_pss, FUN = function(x) names(x[["P_sites_uniq_mm"]]))))
            for (rl in rls_comps) {
                reads_rl_comp <- GRanges()
                seqlevels(reads_rl_comp) <- seqllll
                seqlengths(reads_rl_comp) <- seqleee
                for (comp in names(list_pss)) {
                  if (sum(rl %in% names(list_pss[[comp]][["P_sites_uniq_mm"]])) > 
                    0) {
                    oth <- list_pss[[comp]][["P_sites_uniq_mm"]][[rl]]
                    if (!is.null(oth)) {
                      seqlevels(oth) <- seqllll
                      seqlengths(oth) <- seqleee
                      reads_rl_comp <- c(reads_rl_comp, oth)
                    }
                  }
                  uniq_mm_ps_comps[[rl]] <- reads_rl_comp
                }
                
            }
            
            
            list_res <- list(all_ps_comps, uniq_ps_comps, uniq_mm_ps_comps, rang_jun, 
                covall_plus, covall_min, covuni_plus, covuni_min)
            names(list_res) <- c("P_sites_all", "P_sites_uniq", "P_sites_uniq_mm", 
                "junctions", "coverage_all_plus", "coverage_all_minus", "coverage_uniq_plus", 
                "coverage_uniq_minus")
            
            return(list_res)
            
            
        }
        cat(paste("Calculating P-sites positions and junctions ...", date(), "\n"))
        
        P_sites_stats <- reduceByYield(X = opts, YIELD = yiel, MAP = mapp, REDUCE = reduc)
        
        save(P_sites_stats, file = paste(dira, "P_sites_stats", sep = "/"))
        
        cat(paste("Calculating P-sites positions and junctions --- Done!", date(), 
            "\n"))
        
        ps <- P_sites_stats$P_sites_all
        
        cat(paste("Building aggregate P-sites profiles ...", date(), "\n"))
        
        # here do meta and other plots on P-sites positions ()
        ps_signals_win_all <- list()
        ps_signals_tiles_all <- list()
        codons_win_all <- list()
        ps_codons_ratio <- list()
        ps_codons_counts <- list()
        
        as_codons_ratio <- list()
        as_codons_counts <- list()
        
        es_codons_ratio <- list()
        es_codons_counts <- list()
        
        if (length(ps) == 0) {
            print("Not enough signal or low frame preference, skipped P-sites & codon occupancy calculation. \n Set 'all' in the 'choose_readlengths' choice to ignore this warning and get P-sites positions in a new run")
        }
        if (length(ps) > 0) {
            
            for (comp in c("nucl", circs)) {
                tile_cds <- list_txs_ok[[comp]]
                # put positive?
                strand(tile_cds) <- "+"
                tile_cds <- tile_cds[width(tile_cds) > 100]
                checkgen <- TRUE
                if (comp == "nucl") {
                  if ((sum((start(tile_cds) > 50))/length(start(tile_cds))) > 0.33) {
                    tile_cds <- tile_cds[start(tile_cds) > 50]
                    checkgen <- FALSE
                  }
                  if (sum(seqlengths(tile_cds)[as.vector(seqnames(tile_cds))] - end(tile_cds) > 
                    100)/length(tile_cds) > 0.33) {
                    tile_cds <- tile_cds[seqlengths(tile_cds)[as.vector(seqnames(tile_cds))] - 
                      end(tile_cds) > 100]
                    checkgen <- FALSE
                  }
                  ps_comp <- GRangesList(lapply(ps, function(x) {
                    sort(x[!seqnames(x) %in% circs])
                  }))
                } else {
                  ps_comp <- GRangesList(lapply(ps, function(x) {
                    sort(x[seqnames(x) == comp])
                  }))
                }
                
                
                signal_ps <- list()
                signal_ps_nt <- list()
                
                all_ps <- unlist(ps_comp)
                ps_comp[["all"]] <- all_ps
                mp <- mapToTranscripts(all_ps, transcripts = GTF_annotation$exons_txs[as.vector(seqnames(tile_cds))])
                mp$score <- all_ps$score[mp$xHits]
                # put seqlevels for compartment only to speed up stuff
                seqlevels(mp) <- names(GTF_annotation$exons_txs)
                seqlengths(mp) <- sum(width(GTF_annotation$exons_txs))
                covtx <- coverage(mp, weight = mp$score)
                ok_txs <- names(covtx[elementNROWS(runValue(covtx)) > 1])
                if (length(ok_txs) < 5) {
                  ps_signals_tiles_all[[comp]] <- c()
                  ps_signals_win_all[[comp]] <- c()
                  next
                }
                # keep all txs if small n of txs
                
                if (length(ok_txs) >= 5 & length(covtx) < 250) {
                  ok_txs <- names(covtx)
                }
                
                if (length(ok_txs) > 5000) {
                  cnts_txss_agg <- sort(sum(covtx), decreasing = TRUE)
                  ok_txs <- names(cnts_txss_agg)[seq_len(5000)]
                }
                
                if (fast_mode[bammo] == TRUE) {
                  if (length(ok_txs) > 500) {
                    cnts_txss_agg <- sort(sum(covtx), decreasing = TRUE)
                    ok_txs <- names(cnts_txss_agg)[seq_len(500)]
                  }
                }
                
                
                fivs <- GRanges(seqnames = seqnames(GTF_annotation$cds_txs_coords), 
                  ranges = IRanges(start = 1, end = start(GTF_annotation$cds_txs_coords)), 
                  strand = "*")
                fivs <- fivs[as.character(seqnames(fivs)) %in% ok_txs]
                fivs_gen <- unlist(pmapFromTranscripts(fivs, transcripts = GTF_annotation$exons_txs[as.character(seqnames(fivs))], 
                  ignore.strand = FALSE))
                fivs_gen <- fivs_gen[fivs_gen$hit]
                fivs_gen <- split(fivs_gen, names(fivs_gen))
                
                
                threes <- GRanges(seqnames = seqnames(GTF_annotation$cds_txs_coords), 
                  ranges = IRanges(start = end(GTF_annotation$cds_txs_coords), end = GTF_annotation$cds_txs_coords$lentx), 
                  strand = "*")
                threes <- threes[as.character(seqnames(threes)) %in% ok_txs]
                threes_gen <- unlist(pmapFromTranscripts(threes, transcripts = GTF_annotation$exons_txs[as.character(seqnames(threes))], 
                  ignore.strand = FALSE))
                threes_gen <- threes_gen[threes_gen$hit]
                threes_gen <- split(threes_gen, names(threes_gen))
                
                cds_gen <- GTF_annotation$cds_txs[ok_txs]
                
                
                tile_cds <- tile_cds[as.character(seqnames(tile_cds)) %in% ok_txs]
                ok_txs <- unique(as.character(seqnames(tile_cds)))
                seqlevels(tile_cds) <- ok_txs
                list_covs <- list()
                list_covs[["all"]] <- covtx
                
                tile_5 <- GRanges(seqnames(tile_cds), ranges = IRanges(start = 1, 
                  end = start(tile_cds)))
                tile_3 <- GRanges(seqnames(tile_cds), ranges = IRanges(start = end(tile_cds), 
                  end = seqlengths(tile_cds)[as.vector(seqnames(tile_cds))]))
                no5utr <- width(tile_5) < 51
                no3utr <- width(tile_3) < 51
                
                ex_annot <- GTF_annotation$exons_txs
                
                if (sum(no5utr) > 0) {
                  tx_notok <- seqnames(tile_5)[no5utr]
                  annot_notok <- ex_annot[tx_notok]
                  annot_ok <- GRangesList(lapply(annot_notok, function(x) {
                    if (length(x) == 0) {
                      return(x)
                    }
                    x[1] <- resize(x[1], width = width(x[1]) + 51, fix = "end")
                    x
                  }))
                  ex_annot[names(annot_ok)] <- annot_ok
                  
                  seqlengths(tile_cds)[as.vector(tx_notok)] <- sum(width(annot_ok))
                  tile_cds[no5utr] <- shift(tile_cds[no5utr], +51)
                  tile_5 <- GRanges(seqnames(tile_cds), ranges = IRanges(start = 1, 
                    end = start(tile_cds)))
                  tile_3 <- GRanges(seqnames(tile_cds), ranges = IRanges(start = end(tile_cds), 
                    end = seqlengths(tile_cds)[as.vector(seqnames(tile_cds))]))
                  
                  mp <- mapToTranscripts(all_ps, transcripts = ex_annot[as.vector(seqnames(tile_cds))])
                  seqlevels(mp) <- seqlevels(tile_cds)
                  seqlengths(mp) <- seqlengths(tile_cds)
                  mp$score <- all_ps$score[mp$xHits]
                  
                  covtx <- coverage(mp, weight = mp$score)
                  list_covs[["all"]] <- covtx
                }
                
                if (sum(no3utr) > 0) {
                  tx_notok <- seqnames(tile_3)[no3utr]
                  annot_notok <- ex_annot[tx_notok]
                  annot_ok <- GRangesList(lapply(annot_notok, function(x) {
                    if (length(x) == 0) {
                      return(x)
                    }
                    x[length(x)] <- resize(x[length(x)], width = width(x[length(x)]) + 
                      51, fix = "start")
                    x
                  }))
                  ex_annot[names(annot_ok)] <- annot_ok
                  
                  seqlengths(tile_cds)[as.vector(tx_notok)] <- sum(width(annot_ok))
                  
                  tile_5 <- GRanges(seqnames(tile_cds), ranges = IRanges(start = 1, 
                    end = start(tile_cds)))
                  tile_3 <- GRanges(seqnames(tile_cds), ranges = IRanges(start = end(tile_cds), 
                    end = seqlengths(tile_cds)[as.vector(seqnames(tile_cds))]))
                  
                  mp <- mapToTranscripts(all_ps, transcripts = ex_annot[as.vector(seqnames(tile_cds))])
                  seqlevels(mp) <- seqlevels(tile_cds)
                  seqlengths(mp) <- seqlengths(tile_cds)
                  mp$score <- all_ps$score[mp$xHits]
                  
                  covtx <- coverage(mp, weight = mp$score)
                  list_covs[["all"]] <- covtx
                  
                }
                
                
                ps_tiles <- DataFrameList()
                ps_win <- DataFrameList()
                cod_win <- DataFrameList()
                ps_cod <- DataFrameList()
                ps_cod_rat <- DataFrameList()
                
                as_cod <- DataFrameList()
                as_cod_rat <- DataFrameList()
                es_cod <- DataFrameList()
                es_cod_rat <- DataFrameList()
                
                for (len in c("all", names(ps_comp))) {
                  
                  if ((sum(no3utr) + sum(no5utr)) == 0) {
                    mp <- mapToTranscripts(ps_comp[[len]], transcripts = fivs_gen)
                    mp$score <- ps_comp[[len]]$score[mp$xHits]
                    seqlevels(mp) <- names(fivs_gen)
                    seqlengths(mp) <- sum(width(fivs_gen))
                    cov_5 <- coverage(mp, weight = mp$score)
                    
                    
                    mp <- mapToTranscripts(ps_comp[[len]], transcripts = threes_gen)
                    mp$score <- ps_comp[[len]]$score[mp$xHits]
                    seqlevels(mp) <- names(threes_gen)
                    seqlengths(mp) <- sum(width(threes_gen))
                    cov_3 <- coverage(mp, weight = mp$score)
                    
                    mp <- mapToTranscripts(ps_comp[[len]], transcripts = cds_gen)
                    mp$score <- ps_comp[[len]]$score[mp$xHits]
                    seqlevels(mp) <- names(cds_gen)
                    seqlengths(mp) <- sum(width(cds_gen))
                    strand(mp) <- "+"
                    cov_cds <- coverage(mp, weight = mp$score)
                    cov_cds_a <- suppressWarnings(coverage(shift(mp, shift = 3), 
                      weight = mp$score))
                    cov_cds_e <- suppressWarnings(coverage(shift(mp, shift = -3), 
                      weight = mp$score))
                  }
                  
                  if ((sum(no3utr) + sum(no5utr)) > 0) {
                    
                    mp <- mapToTranscripts(ps_comp[[len]], transcripts = ex_annot[as.vector(seqnames(tile_cds))])
                    mp$score <- ps_comp[[len]]$score[mp$xHits]
                    seqlevels(mp) <- seqlevels(tile_cds)
                    seqlengths(mp) <- seqlengths(tile_cds)
                    strand(mp) <- "+"
                    covtx <- coverage(mp, weight = mp$score)
                    covtx <- covtx[ok_txs]
                    cov_5 <- covtx[tile_5]
                    cov_3 <- covtx[tile_3]
                    cov_cds <- covtx[tile_cds]
                    cov_cds_a <- suppressWarnings(coverage(shift(mp, shift = 3), 
                      weight = mp$score))[tile_cds]
                    cov_cds_e <- suppressWarnings(coverage(shift(mp, shift = -3), 
                      weight = mp$score))[tile_cds]
                    
                  }
                  
                  ps_tiles_5 <- DataFrame(t(sapply(cov_5, function(x) {
                    clos <- 50 * (round(length(x)/50, digits = 0) + 1)
                    idx <- as.integer(seq(1, length(x), length.out = clos))
                    colMeans(matrix(x[idx], ncol = 50))
                  })))
                  ps_win_5 <- DataFrame(t(sapply(cov_5, function(x) {
                    as.vector(x)[c(seq_len(25), (length(x) - 25):(length(x) - 1))]
                  })))
                  
                  
                  ps_tiles_3 <- DataFrame(t(sapply(cov_3, function(x) {
                    clos <- 50 * (round(length(x)/50, digits = 0) + 1)
                    idx <- as.integer(seq(1, length(x), length.out = clos))
                    colMeans(matrix(x[idx], ncol = 50))
                  })))
                  ps_win_3 <- DataFrame(t(sapply(cov_3, function(x) {
                    as.vector(x)[c(2:26, (length(x) - 24):length(x))]
                  })))
                  
                  
                  ps_tiles_cds <- DataFrame(t(sapply(cov_cds, function(x) {
                    clos <- 100 * (round(length(x)/100, digits = 0) + 1)
                    idx <- as.integer(seq(1, length(x), length.out = clos))
                    colMeans(matrix(x[idx], ncol = 100))
                  })))
                  ps_win_cds <- DataFrame(t(sapply(cov_cds, function(x) {
                    rnd <- as.integer(length(x)/2)%%3
                    mid <- (as.integer(length(x)/2) - rnd)
                    as.vector(x)[c(seq_len(33), (mid - 17):(mid + 15), (length(x) - 
                      32):length(x))]
                  })))
                  
                  as_win_cds <- DataFrame(t(sapply(cov_cds_a, function(x) {
                    rnd <- as.integer(length(x)/2)%%3
                    mid <- (as.integer(length(x)/2) - rnd)
                    as.vector(x)[c(seq_len(33), (mid - 17):(mid + 15), (length(x) - 
                      32):length(x))]
                  })))
                  
                  es_win_cds <- DataFrame(t(sapply(cov_cds_e, function(x) {
                    rnd <- as.integer(length(x)/2)%%3
                    mid <- (as.integer(length(x)/2) - rnd)
                    as.vector(x)[c(seq_len(33), (mid - 17):(mid + 15), (length(x) - 
                      32):length(x))]
                  })))
                  
                  # select txs, output codon usage
                  txs_seqq <- extractTranscriptSeqs(x = genome_seq, transcripts = GTF_annotation$cds_txs[names(cov_cds)])
                  
                  gco <- as.character(names(getGeneticCode("1")))
                  names(gco) <- as.character(getGeneticCode("1"))
                  
                  # Get genetic code for compartment
                  
                  gen_cods <- GTF_annotation$genetic_codes
                  circ_cods <- which(comp == rownames(gen_cods))
                  if (length(circ_cods) > 0) {
                    
                    
                    gco <- as.character(names(getGeneticCode(gen_cods$genetic_code[circ_cods])))
                    names(gco) <- as.character(getGeneticCode(gen_cods$genetic_code[circ_cods]))
                    
                  }
                  
                  rnd <- as.integer(width(txs_seqq)/2)%%3
                  st_pos <- rep(1, length(txs_seqq))
                  mid_pos <- (as.integer(width(txs_seqq)/2) - rnd) - 17
                  end_pos <- width(txs_seqq) - 32
                  
                  iters <- seq(0, 10)
                  cod_counts <- c()
                  psit_counts <- c()
                  asit_counts <- c()
                  esit_counts <- c()
                  # Calculate codon occurrence and occupancy for different CDS sections OUTPUT AA?
                  for (i in iters) {
                    a <- narrow(txs_seqq, start = st_pos + 3 * i, end = st_pos + 
                      (2 + 3 * i))
                    coood <- as.character(a)
                    at <- table(a)
                    cod_cntt <- rep(0, length(gco))
                    names(cod_cntt) <- gco
                    cod_cntt[names(at)] <- as.numeric(at)
                    cod_counts <- cbind(cod_counts, cod_cntt)
                    pt <- rowSums(as.matrix(ps_win_cds[, (1 + 3 * i):(3 + 3 * i)]))
                    names(pt) <- coood
                    pt <- aggregate(pt, list(names(pt)), sum)
                    ps_cntt <- rep(0, length(gco))
                    names(ps_cntt) <- gco
                    ps_cntt[pt[, 1]] <- pt[, 2]
                    psit_counts <- cbind(psit_counts, ps_cntt)
                    
                    
                    pt <- rowSums(as.matrix(as_win_cds[, (1 + 3 * i):(3 + 3 * i)]))
                    names(pt) <- coood
                    pt <- aggregate(pt, list(names(pt)), sum)
                    ps_cntt <- rep(0, length(gco))
                    names(ps_cntt) <- gco
                    ps_cntt[pt[, 1]] <- pt[, 2]
                    asit_counts <- cbind(asit_counts, ps_cntt)
                    
                    
                    pt <- rowSums(as.matrix(es_win_cds[, (1 + 3 * i):(3 + 3 * i)]))
                    names(pt) <- coood
                    pt <- aggregate(pt, list(names(pt)), sum)
                    ps_cntt <- rep(0, length(gco))
                    names(ps_cntt) <- gco
                    ps_cntt[pt[, 1]] <- pt[, 2]
                    esit_counts <- cbind(esit_counts, ps_cntt)
                    
                  }
                  # mid
                  for (i in iters) {
                    a <- narrow(txs_seqq, start = mid_pos + 3 * i, end = mid_pos + 
                      (2 + 3 * i))
                    coood <- as.character(a)
                    at <- table(a)
                    cod_cntt <- rep(0, length(gco))
                    names(cod_cntt) <- gco
                    cod_cntt[names(at)] <- as.numeric(at)
                    cod_counts <- cbind(cod_counts, cod_cntt)
                    pt <- rowSums(as.matrix(ps_win_cds[, (34 + 3 * i):(36 + 3 * i)]))
                    names(pt) <- coood
                    pt <- aggregate(pt, list(names(pt)), sum)
                    ps_cntt <- rep(0, length(gco))
                    names(ps_cntt) <- gco
                    ps_cntt[pt[, 1]] <- pt[, 2]
                    psit_counts <- cbind(psit_counts, ps_cntt)
                    
                    pt <- rowSums(as.matrix(as_win_cds[, (34 + 3 * i):(36 + 3 * i)]))
                    names(pt) <- coood
                    pt <- aggregate(pt, list(names(pt)), sum)
                    ps_cntt <- rep(0, length(gco))
                    names(ps_cntt) <- gco
                    ps_cntt[pt[, 1]] <- pt[, 2]
                    asit_counts <- cbind(asit_counts, ps_cntt)
                    
                    pt <- rowSums(as.matrix(es_win_cds[, (34 + 3 * i):(36 + 3 * i)]))
                    names(pt) <- coood
                    pt <- aggregate(pt, list(names(pt)), sum)
                    ps_cntt <- rep(0, length(gco))
                    names(ps_cntt) <- gco
                    ps_cntt[pt[, 1]] <- pt[, 2]
                    esit_counts <- cbind(esit_counts, ps_cntt)
                  }
                  # end
                  for (i in iters) {
                    a <- narrow(txs_seqq, start = end_pos + 3 * i, end = end_pos + 
                      (2 + 3 * i))
                    coood <- as.character(a)
                    at <- table(a)
                    cod_cntt <- rep(0, length(gco))
                    names(cod_cntt) <- gco
                    cod_cntt[names(at)] <- as.numeric(at)
                    cod_counts <- cbind(cod_counts, cod_cntt)
                    pt <- rowSums(as.matrix(ps_win_cds[, (67 + 3 * i):(69 + 3 * i)]))
                    names(pt) <- coood
                    pt <- aggregate(pt, list(names(pt)), sum)
                    ps_cntt <- rep(0, length(gco))
                    names(ps_cntt) <- gco
                    ps_cntt[pt[, 1]] <- pt[, 2]
                    psit_counts <- cbind(psit_counts, ps_cntt)
                    
                    pt <- rowSums(as.matrix(as_win_cds[, (67 + 3 * i):(69 + 3 * i)]))
                    names(pt) <- coood
                    pt <- aggregate(pt, list(names(pt)), sum)
                    ps_cntt <- rep(0, length(gco))
                    names(ps_cntt) <- gco
                    ps_cntt[pt[, 1]] <- pt[, 2]
                    asit_counts <- cbind(asit_counts, ps_cntt)
                    
                    pt <- rowSums(as.matrix(es_win_cds[, (67 + 3 * i):(69 + 3 * i)]))
                    names(pt) <- coood
                    pt <- aggregate(pt, list(names(pt)), sum)
                    ps_cntt <- rep(0, length(gco))
                    names(ps_cntt) <- gco
                    ps_cntt[pt[, 1]] <- pt[, 2]
                    esit_counts <- cbind(esit_counts, ps_cntt)
                    
                  }
                  
                  colnames(psit_counts) <- c(paste("first", seq(1, 11), sep = "_"), 
                    paste("mid", seq(1, 11), sep = "_"), paste("last", seq(1, 11), 
                      sep = "_"))
                  colnames(esit_counts) <- c(paste("first", seq(1, 11), sep = "_"), 
                    paste("mid", seq(1, 11), sep = "_"), paste("last", seq(1, 11), 
                      sep = "_"))
                  colnames(asit_counts) <- c(paste("first", seq(1, 11), sep = "_"), 
                    paste("mid", seq(1, 11), sep = "_"), paste("last", seq(1, 11), 
                      sep = "_"))
                  colnames(cod_counts) <- c(paste("first", seq(1, 11), sep = "_"), 
                    paste("mid", seq(1, 11), sep = "_"), paste("last", seq(1, 11), 
                      sep = "_"))
                  ratio_psit_cod <- psit_counts/cod_counts
                  ratio_asit_cod <- asit_counts/cod_counts
                  ratio_esit_cod <- esit_counts/cod_counts
                  
                  psit_counts <- DataFrame(psit_counts)
                  asit_counts <- DataFrame(asit_counts)
                  esit_counts <- DataFrame(esit_counts)
                  cod_counts <- DataFrame(cod_counts)
                  ratio_psit_cod <- DataFrame(ratio_psit_cod)
                  ratio_asit_cod <- DataFrame(ratio_asit_cod)
                  ratio_esit_cod <- DataFrame(ratio_esit_cod)
                  
                  gcall <- paste(gco, names(gco), sep = ";")
                  
                  #remove 'Ns' our counts
                  non_n_rows <- grep(x=rownames(psit_counts),v=TRUE,patt='N',inv=T)
                  psit_counts <- psit_counts[non_n_rows,]
                  asit_counts <- asit_counts[non_n_rows,]
                  esit_counts <- esit_counts[non_n_rows,]
                  ratio_psit_cod <- ratio_psit_cod[non_n_rows,]
                  ratio_asit_cod <- ratio_asit_cod[non_n_rows,]
                  ratio_esit_cod <- ratio_esit_cod[non_n_rows,]
                  cod_counts <- cod_counts[non_n_rows,]

                  rownames(psit_counts) <- gcall[match(rownames(psit_counts), gco)]
                  rownames(asit_counts) <- gcall[match(rownames(asit_counts), gco)]
                  rownames(esit_counts) <- gcall[match(rownames(esit_counts), gco)]
                  rownames(cod_counts) <- gcall[match(rownames(cod_counts), gco)]
                  
                  rownames(ratio_psit_cod) <- gcall[match(rownames(ratio_psit_cod), 
                    gco)]
                  rownames(ratio_asit_cod) <- gcall[match(rownames(ratio_asit_cod), 
                    gco)]
                  rownames(ratio_esit_cod) <- gcall[match(rownames(ratio_esit_cod), 
                    gco)]
                  
                  
                  cod_win[[len]] <- cod_counts
                  ps_cod[[len]] <- psit_counts
                  ps_cod_rat[[len]] <- ratio_psit_cod
                  
                  as_cod[[len]] <- asit_counts
                  as_cod_rat[[len]] <- ratio_asit_cod
                  
                  
                  es_cod[[len]] <- esit_counts
                  es_cod_rat[[len]] <- ratio_esit_cod
                  
                  
                  tls <- cbind(ps_tiles_5, ps_tiles_cds, ps_tiles_3)
                  colnames(tls) <- c(paste("5_UTR", seq_along(ps_tiles_5[1, ]), sep = "_"), 
                    paste("CDS", seq_along(ps_tiles_cds[1, ]), sep = "_"), paste("3_UTR", 
                      seq_along(ps_tiles_3[1, ]), sep = "_"))
                  nts <- cbind(ps_win_5, ps_win_cds, ps_win_3)
                  colnames(nts) <- c(paste("5_UTR", seq_along(ps_win_5[1, ]), sep = "_"), 
                    paste("CDS", seq_along(ps_win_cds[1, ]), sep = "_"), paste("3_UTR", 
                      seq_along(ps_win_3[1, ]), sep = "_"))
                  
                  
                  ps_tiles[[len]] <- tls
                  ps_win[[len]] <- nts
                  
                }
                
                codons_win_all[[comp]] <- cod_win
                ps_codons_ratio[[comp]] <- ps_cod_rat
                ps_codons_counts[[comp]] <- ps_cod
                
                as_codons_ratio[[comp]] <- as_cod_rat
                as_codons_counts[[comp]] <- as_cod
                
                
                es_codons_ratio[[comp]] <- es_cod_rat
                es_codons_counts[[comp]] <- es_cod
                
                ps_signals_tiles_all[[comp]] <- ps_tiles
                ps_signals_win_all[[comp]] <- ps_win
                
            }
        }
        
        profiles_P_sites <- list(ps_signals_tiles_all, ps_signals_win_all, codons_win_all, 
            ps_codons_counts, ps_codons_ratio, as_codons_counts, as_codons_ratio, 
            es_codons_counts, es_codons_ratio)
        
        names(profiles_P_sites) <- c("P_sites_bins", "P_sites_subcodon", "Codon_counts", 
            "P_sites_percodon", "P_sites_percodon_ratio", "A_sites_percodon", "A_sites_percodon_ratio", 
            "E_sites_percodon", "E_sites_percodon_ratio")
        save(profiles_P_sites, file = paste(dira, "profiles_P_sites", sep = "/"))
        # save(list =
        # c('ps_signals_tiles','ps_signals_tiles_all','ps_signals_win','ps_signals_win_all'),file
        # = 'ps_results_preprjan15')
        cat(paste("Building aggregate P-sites profiles --- Done!", date(), "\n"))
        
        cat(paste("Exporting files ...", date(), "\n"))
        
        merged_all_ps <- unlist(P_sites_stats$P_sites_all)
        if (length(merged_all_ps) > 0) {
            covv_pl <- coverage(merged_all_ps[strand(merged_all_ps) == "+"], weight = merged_all_ps[strand(merged_all_ps) == 
                "+"]$score)
            covv_pl <- GRanges(covv_pl)
            covv_pl <- covv_pl[covv_pl$score > 0]
            covv_min <- coverage(merged_all_ps[strand(merged_all_ps) == "-"], weight = merged_all_ps[strand(merged_all_ps) == 
                "-"]$score)
            covv_min <- GRanges(covv_min)
            covv_min <- covv_min[covv_min$score > 0]
            strand(covv_pl) <- "+"
            strand(covv_min) <- "-"
            
            merged_all_ps <- sort(c(covv_pl, covv_min))
            
            export(covv_pl, con = paste(dest_name, "_P_sites_plus.bedgraph", sep = ""))
            export(covv_min, con = paste(dest_name, "_P_sites_minus.bedgraph", sep = ""))
        }
        
        merged_uniq_ps <- unlist(P_sites_stats$P_sites_uniq)
        if (length(merged_uniq_ps) > 0) {
            
            covv_pl <- coverage(merged_uniq_ps[strand(merged_uniq_ps) == "+"], weight = merged_uniq_ps[strand(merged_uniq_ps) == 
                "+"]$score)
            covv_pl <- GRanges(covv_pl)
            covv_pl <- covv_pl[covv_pl$score > 0]
            
            covv_min <- coverage(merged_uniq_ps[strand(merged_uniq_ps) == "-"], weight = merged_uniq_ps[strand(merged_uniq_ps) == 
                "-"]$score)
            covv_min <- GRanges(covv_min)
            covv_min <- covv_min[covv_min$score > 0]
            strand(covv_pl) <- "+"
            strand(covv_min) <- "-"
            merged_uniq_ps <- sort(c(covv_pl, covv_min))
            
            
            export(covv_pl, con = paste(dest_name, "_P_sites_uniq_plus.bedgraph", 
                sep = ""))
            export(covv_min, con = paste(dest_name, "_P_sites_uniq_minus.bedgraph", 
                sep = ""))
        }
        
        
        merged_uniq_mm_ps <- unlist(P_sites_stats$P_sites_uniq_mm)
        if (length(merged_uniq_mm_ps) > 0) {
            
            covv_pl <- coverage(merged_uniq_mm_ps[strand(merged_uniq_mm_ps) == "+"], 
                weight = merged_uniq_mm_ps[strand(merged_uniq_mm_ps) == "+"]$score)
            covv_pl <- GRanges(covv_pl)
            covv_pl <- covv_pl[covv_pl$score > 0]
            
            covv_min <- coverage(merged_uniq_mm_ps[strand(merged_uniq_mm_ps) == "-"], 
                weight = merged_uniq_mm_ps[strand(merged_uniq_mm_ps) == "-"]$score)
            covv_min <- GRanges(covv_min)
            covv_min <- covv_min[covv_min$score > 0]
            strand(covv_pl) <- "+"
            strand(covv_min) <- "-"
            merged_uniq_mm_ps <- sort(c(covv_pl, covv_min))
        }
        
        if (!normalize_cov[bammo]) {
            export(P_sites_stats$coverage_all_plus, con = paste(dest_name, "_coverage_plus.bedgraph", 
                sep = ""))
            export(P_sites_stats$coverage_all_min, con = paste(dest_name, "_coverage_minus.bedgraph", 
                sep = ""))
            export(P_sites_stats$coverage_uniq_plus, con = paste(dest_name, "_coverage_uniq_plus.bedgraph", 
                sep = ""))
            export(P_sites_stats$coverage_uniq_min, con = paste(dest_name, "_coverage_uniq_minus.bedgraph", 
                sep = ""))
        }
        
        
        
        if (normalize_cov[bammo]) {
            bwpl <- P_sites_stats$coverage_all_plus
            bwmn <- P_sites_stats$coverage_all_min
            
            correction <- (sum(as.numeric(unlist(runValue(bwpl)) * unlist(runLength(bwpl)))) + 
                sum(as.numeric(unlist(runValue(bwmn)) * unlist(runLength(bwmn)))))/1e+06
            okround <- sort(unique(unlist(runValue(bwpl))))[2]
            okround <- nchar(format(okround/correction, digits = 3, nsmall = 2)) - 
                2
            
            runValue(bwpl) <- round(runValue(bwpl)/correction, digits = okround)
            runValue(bwmn) <- round(runValue(bwmn)/correction, digits = okround)
            
            export(bwpl, con = paste(dest_name, "_coverage_plus.bedgraph", sep = ""))
            export(bwmn, con = paste(dest_name, "_coverage_minus.bedgraph", sep = ""))
            
            
            bwpl <- P_sites_stats$coverage_uniq_plus
            bwmn <- P_sites_stats$coverage_uniq_min
            
            correction <- (sum(as.numeric(unlist(runValue(bwpl)) * unlist(runLength(bwpl)))) + 
                sum(as.numeric(unlist(runValue(bwmn)) * unlist(runLength(bwmn)))))/1e+06
            okround <- sort(unique(unlist(runValue(bwpl))))[2]
            okround <- nchar(format(okround/correction, digits = 3, nsmall = 2)) - 
                2
            
            runValue(bwpl) <- round(runValue(bwpl)/correction, digits = okround)
            runValue(bwmn) <- round(runValue(bwmn)/correction, digits = okround)
            
            export(bwpl, con = paste(dest_name, "_coverage_uniq_plus.bedgraph", sep = ""))
            export(bwmn, con = paste(dest_name, "_coverage_uniq_minus.bedgraph", 
                sep = ""))
            
        }
        
        
        juns <- P_sites_stats$junctions
        save(juns, file = paste(dest_name, "_junctions", sep = ""))
        
        
        for_ORFquant <- list(merged_all_ps, merged_uniq_ps, merged_uniq_mm_ps, juns)
        names(for_ORFquant) <- c("P_sites_all", "P_sites_uniq", "P_sites_uniq_mm", 
            "junctions")
        
        
        # Store top 50 occupied regions (possible contaminants)
        
        
        rds_st <- read_stats$reads_pos1
        rds_st <- sort(unlist(rds_st))
        
        regions <- list(reduce(unlist(GTF_annotation$cds_genes)), GTF_annotation$fiveutrs, 
            GTF_annotation$threeutrs, GTF_annotation$ncIsof, GTF_annotation$ncRNAs, 
            GTF_annotation$introns, GTF_annotation$intergenicRegions)
        names(regions) <- c("cds", "fiveutrs", "threeutrs", "ncIsof", "ncRNAs", "introns", 
            "intergenic")
        
        regions <- GRangesList(endoapply(regions, function(x) {
            mcols(x) <- NULL
            x
        }))
        
        rds_st <- rds_st[seqnames(rds_st) %in% unique(seqlevels(regions))]
        rds_st$len_adj <- as.numeric(names(rds_st))
        names(rds_st) <- NULL
        
        rds_st <- rds_st[order(rds_st$score, decreasing = TRUE)]
        rds_st$pct <- round(rds_st$score/(sum(rds_st$score)/100), digits = 4)
        top50 <- head(rds_st, 50)
        top50 <- resize(top50, width = top50$len_adj, fix = "start")
        top50$seq <- getSeq(genome_seq, top50)
        top50$region <- CharacterList("")
        
        ovg <- findOverlaps(top50, regions)
        ovg <- IntegerList(split(subjectHits(ovg), queryHits(ovg)))
        subsr <- as.numeric(names(ovg))
        a_reg <- CharacterList(lapply(ovg, FUN = function(x) {
            unique(names(regions)[x])
        }))
        top50$region[subsr] <- a_reg
        
        ovg <- findOverlaps(top50, GTF_annotation$genes)
        ovg <- IntegerList(split(subjectHits(ovg), queryHits(ovg)))
        subs <- as.numeric(names(ovg))
        a_nam <- CharacterList(lapply(ovg, FUN = function(x) {
            unique(names(GTF_annotation$genes)[x])
        }))
        a_biot <- CharacterList(lapply(a_nam, FUN = function(x) {
            GTF_annotation$trann$gene_biotype[match(x, GTF_annotation$trann$gene_id)]
        }))
        top50$gene <- CharacterList("")
        top50$gene_biotype <- CharacterList("")
        top50$gene[subs] <- a_nam
        top50$gene_biotype[subs] <- a_biot
        sequence_analysis <- top50
        
        # save as list_results
        
        save(sequence_analysis, file = paste(dira, "sequence_analysis", sep = "/"))
        
        
        res_all <- list(read_stats, profiles_fivepr, selection_cutoffs, P_sites_stats, 
            profiles_P_sites, sequence_analysis)
        
        names(res_all) <- c("read_stats", "profiles_fivepr", "selection_cutoffs", 
            "P_sites_stats", "profiles_P_sites", "sequence_analysis")
        if (length(merged_all_ps) > 0) {
            
            save(for_ORFquant, file = paste(dest_name, "_for_ORFquant", sep = ""))
        }
        finn <- lapply(res_all$selection_cutoffs$results_choice, function(x) {
            x$data
        })
        nope <- length(finn)
        
        namfinn <- rep(names(finn), elementNROWS(finn))
        finn <- do.call(rbind, args = finn)
        finn$comp <- namfinn
        
        pct_lib <- round(t(read_stats$rld/sum(read_stats$rld) * 100), digits = 4)
        rownames(pct_lib) <- gsub(rownames(pct_lib), pattern = "reads_", replacement = "")
        finn$pct_map <- 0
        for (i in unique(namfinn)) {
            finn$pct_map[finn$comp == i] <- pct_lib[match(finn$read_length[finn$comp == 
                i], rownames(pct_lib)), i]
        }
        if (nope == 0) {
            finn <- NULL
        }
        res_all$summary_P_sites <- finn
        if (nope == 0) {
            res_all$summary_P_sites <- c("")
        }
        if (write_tmp_files) {
            save(res_all, file = paste(dest_name, "_results_RiboseQC_all", sep = ""))
        }
        write.table(finn, file = paste(dest_name, "_P_sites_calcs", sep = ""), sep = "\t", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)
        unlink(dira, recursive = TRUE)
        res_all$read_stats$reads_pos1 <- ""
        res_all$P_sites_stats <- ""
        
        nms <- names(res_all$profiles_fivepr$five_prime_bins)
        for (i in nms) {
            fvp <- res_all$profiles_fivepr$five_prime_bins[[i]]
            res_all$profiles_fivepr$five_prime_bins[[i]] <- lapply(fvp, function(x) {
                colSums(as.matrix(x))
            })
        }
        nms <- names(res_all$profiles_fivepr$five_prime_subcodon)
        for (i in nms) {
            fvp <- res_all$profiles_fivepr$five_prime_subcodon[[i]]
            res_all$profiles_fivepr$five_prime_subcodon[[i]] <- lapply(fvp, function(x) {
                colSums(as.matrix(x))
            })
        }
        
        nms <- names(res_all$profiles_P_sites$P_sites_bins)
        for (i in nms) {
            fvp <- res_all$profiles_P_sites$P_sites_bins[[i]]
            res_all$profiles_P_sites$P_sites_bins[[i]] <- lapply(fvp, function(x) {
                colSums(as.matrix(x))
            })
        }
        nms <- names(res_all$profiles_P_sites$P_sites_subcodon)
        for (i in nms) {
            fvp <- res_all$profiles_P_sites$P_sites_subcodon[[i]]
            res_all$profiles_P_sites$P_sites_subcodon[[i]] <- lapply(fvp, function(x) {
                colSums(as.matrix(x))
            })
        }
        resfile <- paste(dest_name, "_results_RiboseQC", sep = "")
        message(paste0('Writing results to ',resfile))
        save(res_all, file = resfile)
        rm(res_all, read_stats, profiles_fivepr, selection_cutoffs, P_sites_stats, 
            profiles_P_sites, sequence_analysis)
        gici <- gc()
        cat(paste("Exporting files --- Done!", date(), "\n\n"))
        
        resfilelist <- c(resfilelist,resfile)
    }
    
      if (create_report) {
          cat(paste("Creating html report in ", report_file, " ... ", date(), "\n", 
              sep = ""))
          
          create_html_report(input_files = resfilelist, 
              input_sample_names = sample_names, output_file = report_file, extended = extended_report)
          cat(paste("Creating html report --- Done!", date(), "\n\n"))
          
          if (pdf_plots) {
              create_pdfs_from_rds_objects(output_rds_path = paste(report_file, "_plots/rds/", 
                  sep = ""))
          }
          
      }
    
    return(resfilelist)
    
}





