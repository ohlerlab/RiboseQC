## ----setup_vig, include=FALSE--------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----load_lib------------------------------------------------------------
library("RiboseQC")
tmp_dir <- tempdir()
file.path()

## ----create_dirs_hum-----------------------------------------------------

download.file("https://drive.google.com/uc?export=download&id=1n6UA5cSz6djx0dY_7sI57177T-kyZhhT",destfile = file.path(tmp_dir,"test_human.2bit"))
download.file("https://drive.google.com/uc?export=download&id=19Am9-iMEyB-AcIsVRdIqrF-BVdQYrcaI",destfile = file.path(tmp_dir,"test_human.gtf"))


## ----create_annot_hum----------------------------------------------------

prepare_annotation_files(annotation_directory = tmp_dir,
                         twobit_file = file.path(tmp_dir,"test_human.2bit"),
                         gtf_file = file.path(tmp_dir,"test_human.gtf"),scientific_name = "Human.test",
                         annotation_name = "genc25_22M",
                         export_bed_tables_TxDb = FALSE,
                         forge_BSgenome = TRUE,
                         create_TxDb = TRUE)


## ----load_hum------------------------------------------------------------
load_annotation("test_human.gtf_Rannot")

## ----gen_hum-------------------------------------------------------------
genome_seq[["chr22"]]
genome_seq[["chrM"]]

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
GTF_annotation$genome_package

## ----human_report_1------------------------------------------------------
if(!file.exists('test_human_hek.bam')) download.file("https://drive.google.com/uc?export=download&id=11PP5y2QH7si81rbEBJsOB-Lt3l_JowRW",destfile = file.path(tmp_dir,"test_human_hek.bam"))

## ----human_report_2------------------------------------------------------
RiboseQC_analysis(annotation_file=file.path(tmp_dir,"test_human.gtf_Rannot"),bam_files = file.path(tmp_dir,"test_human_hek.bam"),report_file = file.path(tmp_dir,"test_human_hek.html"),write_tmp_files = F)

## ----human_report_3,warning=FALSE, results='asis', fig.width=12, fig.height=10, dpi=120----

readRDS("test_human_hek.html_plots/rds/sample1_nucl_4_profiles_P_sites_metagene_subcodon_all")[[1]]


## ----yeast_1-------------------------------------------------------------

download.file("https://drive.google.com/uc?export=download&id=1ul_mqVimaavfhJs76R3TaOyTcUgwJvS2",destfile = file.path(tmp_dir,"test_yeast.2bit"))

download.file("https://drive.google.com/uc?export=download&id=1yffqKkfohGx9y7XRZO05vhjYUXgjiKBF",destfile = file.path(tmp_dir,"test_yeast.gtf"))

prepare_annotation_files(annotation_directory = tmp_dir,
                         twobit_file = file.path(tmp_dir,"test_yeast.2bit"),
                         gtf_file = file.path(tmp_dir,"test_yeast.gtf"),scientific_name = "yeast.test",
                         annotation_name = "yeast_custom",
                         export_bed_tables_TxDb = TRUE,
                         forge_BSgenome = TRUE,
                         create_TxDb = TRUE)

## ----yeast_2-------------------------------------------------------------
load_annotation("test_yeast.gtf_Rannot")

genome_seq[["II"]]
genome_seq[["Mito"]]
GTF_annotation$cds_txs
GTF_annotation$trann

## ----yeast_3-------------------------------------------------------------
download.file("https://drive.google.com/uc?export=download&id=1Nwlf7u-hHC_fX0fUJuJu_dr2WxwGazVJ",destfile = file.path(tmp_dir,"test_yeast_TCP_input.bam"))

download.file("https://drive.google.com/uc?export=download&id=18gQF2SxxsGQqVfhUw-wFo2dBt1mafiAB",destfile = file.path(tmp_dir,"test_yeast_TCP_RS.bam"))

download.file("https://drive.google.com/uc?export=download&id=1HBb58mOt1vmf58yqz-wkZzfxZhBl_tmO",destfile = file.path(tmp_dir,"test_yeast_TCP_SSU.bam"))

bam_filepath_y<-c("test_yeast_TCP_input.bam","test_yeast_TCP_RS.bam","test_yeast_TCP_SSU.bam")

RiboseQC_analysis(annotation_file=file.path(tmp_dir,"test_yeast.gtf_Rannot"),bam_files = bam_filepath_y,fast_mode = T,report_file = file.path(tmp_dir,"test_yeast_TCP.html"),sample_names = c("input","RS","SSU"),dest_names = file.path(tmp_dir,c("input","RS","SSU")),write_tmp_files = F)

## ----yeast_4, warning=FALSE, results='asis', fig.width=12, fig.height=10, dpi=120----
as_ggplot(readRDS("test_yeast_TCP.html_plots/rds/SSU_nucl_4_profiles_fivepr_metagene_subcodon_log2")[[1]])

as_ggplot(readRDS("test_yeast_TCP.html_plots/rds/RS_nucl_4_profiles_fivepr_metagene_subcodon_log2")[[1]])

## ----arab_1--------------------------------------------------------------

download.file("https://drive.google.com/uc?export=download&id=1gTxlR66fRX0GsHsWtJB7piLelM8jhZGu",destfile = file.path(tmp_dir,"test_arabidopsis.2bit"))

download.file("https://drive.google.com/uc?export=download&id=1jYbdvaGfVnMx38X1P85Ounyq9jGwRVgP",destfile = file.path(tmp_dir,"test_arabidopsis.gtf.gz"))

prepare_annotation_files(annotation_directory = tmp_dir,
                         twobit_file = file.path(tmp_dir,"test_arabidopsis.2bit"),
                         gtf_file = file.path(tmp_dir,"test_arabidopsis.gtf.gz"),scientific_name = "arabidopsis.test",
                         annotation_name = "araport11_custom",export_bed_tables_TxDb = TRUE,
                         forge_BSgenome = TRUE,
                         create_TxDb = TRUE)

## ----arab_2--------------------------------------------------------------

annotation="test_arabidopsis.gtf.gz_Rannot"


download.file("https://drive.google.com/uc?export=download&id=144PIK15iJrDSshMWsGbl9zQvT27Z1fzt",destfile = file.path(tmp_dir,"test_arabidopsis_root.bam"))
download.file("https://drive.google.com/uc?export=download&id=1VVdMi0ho0LiRCld6-3vzibtg_CyPSTj3",destfile = file.path(tmp_dir,"test_arabidopsis_shoot.bam"))

bam_filepath=file.path(tmp_dir,c("test_arabidopsis_root.bam","test_arabidopsis_shoot.bam"))

RiboseQC_analysis(annotation_file=file.path(annotation,bam_files = bam_filepath),fast_mode = TRUE,report_file = file.path(tmp_dir,"test_root_shoots.html"),dest_names = c("root","shoots"),sample_names = c("root","shoots"),write_tmp_files = FALSE)

## ----arab_3--------------------------------------------------------------
readRDS("test_root_shoots.html_plots/rds/all_samples_1_readlocdist")[[1]]

## ----arab_4, warning=FALSE, results='asis', fig.width=12, fig.height=16, dpi=150----
as_ggplot(readRDS("./test_root_shoots.html_plots/rds/shoots_ChrC_7_codonusage_positional_A-sites_per_codon_all_zscore")[[1]])

## ----end_sessinf---------------------------------------------------------
session_info()

