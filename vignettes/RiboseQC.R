## ----setup_vig, include=FALSE--------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----load_lib------------------------------------------------------------
library("RiboseQC")

tmp_dir <-  tempdir()



## ----arab_1--------------------------------------------------------------
arab_fasta <- system.file(package='RiboseQC',"ext_data","simp_arab.fa.gz")
file.copy(arab_fasta,tmp_dir)
system(paste0('gunzip ',file.path(tmp_dir,basename(arab_fasta))))
arab_fasta <- file.path(tmp_dir,gsub(x=basename(arab_fasta),'.gz',''))

stopifnot(is(getSeq(FaFile_Circ(FaFile(arab_fasta))),'DNAStringSet'))

arab_gtf <- system.file(package='RiboseQC',"ext_data","simp_arab.gtf.gz")



## ----create_annot_hum----------------------------------------------------

annot_file <- prepare_annotation_files(annotation_directory = tmp_dir,
                         genome_seq = arab_fasta,
                         gtf_file = arab_gtf,
                         scientific_name = "arabidopsis.test",
                         annotation_name = "araport11_custom",export_bed_tables_TxDb = TRUE,
                         forge_BSgenome = FALSE,
                         create_TxDb = TRUE)


## ----load_hum------------------------------------------------------------

load_annotation(file.path(tmp_dir,"simp_arab.gtf.gz_Rannot"))
genome_seq$circularRanges


## ----gen_hum-------------------------------------------------------------
getSeq(genome_seq,GRanges("Chr1:1-100"))
getSeq(genome_seq,GRanges("ChrM:1-100"))

## ----gtf_hum_general_1---------------------------------------------------
GTF_annotation$exons_txs

## ----gtf_hum_general_2---------------------------------------------------
GTF_annotation$cds_txs

## ----gen_hum_cds---------------------------------------------------------
getSeq(genome_seq,GTF_annotation$cds_txs[[4]])

## ----gtf_hum_general_3---------------------------------------------------
GTF_annotation$start_stop_codons

## ----gtf_hum_general_4---------------------------------------------------
GTF_annotation$cds_txs_coords

## ----gtf_hum_general_5---------------------------------------------------
GTF_annotation$trann

## ----gtf_hum_general_6---------------------------------------------------
GTF_annotation$genetic_codes
getGeneticCode(GTF_annotation$genetic_codes["Chr2","genetic_code"])
getGeneticCode(GTF_annotation$genetic_codes["ChrM","genetic_code"])

## ----gtf_hum_general_7---------------------------------------------------
GTF_annotation$genome

## ----report_1------------------------------------------------------------
rootbam <-   system.file(package='RiboseQC',"ext_data","simp_arab_root.bam")

## ----report_2------------------------------------------------------------

annotation=file.path(tmp_dir,"simp_arab.gtf.gz_Rannot")

load_annotation(annotation)

#tmp_dir2="/var/folders/5t/j5mzk77n11b5w4q_54rlq97m0000gn/T//RtmpOFmFiU"

bam_filepath=c(
  system.file(package='RiboseQC',"ext_data","simp_arab_root.bam",mustWork = T),
  system.file(package='RiboseQC',"ext_data","simp_arab_shoot.bam",mustWork = T)
)

## ----run_analysis,eval = FALSE-------------------------------------------
#  
#  resfile <- RiboseQC_analysis(
#    annotation_file=annotation,
#    bam_files = bam_filepath,
#    fast_mode = TRUE,
#  #  report_file = file.path(tmp_dir,"test_root_shoots.html"),
#    dest_names = file.path(tmp_dir,c("root","shoots")),
#    sample_names = c("root","shoots"),
#    genome_seq=arab_fasta,
#    write_tmp_files = FALSE,
#    create_report = FALSE
#  )
#  

## ----load_results--------------------------------------------------------
resfiles=c(
  system.file(package='RiboseQC',"ext_data","root_results_RiboseQC",mustWork = T),
  system.file(package='RiboseQC',"ext_data","shoots_results_RiboseQC",mustWork = T)
)

