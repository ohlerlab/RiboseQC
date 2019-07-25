## ----setup_vig, include=FALSE--------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----load_lib------------------------------------------------------------
library("RiboseQC")
message('start')
if(!exists('tmp_dir')) tmp_dir <-  tempdir()
download.file <- function(file,destfile){
  if(!file.exists(destfile)) utils::download.file(file,  destfile = destfile) 
}


## ----create_dirs_hum-----------------------------------------------------

download.file("https://drive.google.com/uc?export=download&id=19Am9-iMEyB-AcIsVRdIqrF-BVdQYrcaI",destfile = file.path(tmp_dir,"test_human.gtf"))

download.file('https://drive.google.com/uc?export=download&id=1Too_Hrsdd_rDTh2LA4KEccBwnL4ESimR',destfile = file.path(tmp_dir,"test_human.fa"))


## ----create_annot_hum----------------------------------------------------

genome_seq<-FaFile_Circ(FaFile(file.path(tmp_dir,"test_human.fa")))

annot_file <- file.path(tmp_dir,"test_human.gtf_Rannot")

if(! file.exists(annot_file)){
  annot_file <- prepare_annotation_files(annotation_directory = tmp_dir,
                         # twobit_file = file.path(tmp_dir,"test_human.2bit"),
                         genome_seq = file.path(tmp_dir,"test_human.fa"),
                         gtf_file = file.path(tmp_dir,"test_human.gtf"),
                         scientific_name = "Human.test",
                         annotation_name = "genc25_22M",
                         export_bed_tables_TxDb = FALSE,
                         forge_BSgenome = FALSE,
                         create_TxDb = TRUE)
}
message('foo2')

## ----load_hum------------------------------------------------------------

load_annotation(file.path(tmp_dir,"test_human.gtf_Rannot"))
genome_seq$circularRanges


## ----gen_hum-------------------------------------------------------------
getSeq(genome_seq,GRanges("chr22:1-100"))
getSeq(genome_seq,GRanges("chrM:1-100"))

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
getGeneticCode(GTF_annotation$genetic_codes["chr22","genetic_code"])
getGeneticCode(GTF_annotation$genetic_codes["chrM","genetic_code"])

## ----gtf_hum_general_7---------------------------------------------------
GTF_annotation$genome

## ----human_report_1------------------------------------------------------
if(!file.exists(file.path(tmp_dir,'test_human_hek.bam'))) download.file("https://drive.google.com/uc?export=download&id=11PP5y2QH7si81rbEBJsOB-Lt3l_JowRW",destfile = file.path(tmp_dir,"test_human_hek.bam"))

