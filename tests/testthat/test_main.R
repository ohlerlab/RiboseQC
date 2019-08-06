

testthat::test_that("prepare_annotation_files works as expected",{
  
  arab_fasta <- system.file(package='RiboseQC',"ext_data","simp_arab.fa.gz")
  tmp_dir = tempdir()
  file.copy(arab_fasta,tmp_dir)
  system(paste0('gunzip ',file.path(tmp_dir,basename(arab_fasta))))
  arab_fasta <- file.path(tmp_dir,basename(arab_fasta))
  
  arab_gtf <- system.file(package='RiboseQC',"ext_data","simp_arab.gtf.gz")
 
   
  annotation <- prepare_annotation_files(annotation_directory = tempdir(),
                                         genome_seq = arab_fasta,
                                         gtf_file = arab_gtf,
                                         scientific_name = "arabidopsis.test",
                                         annotation_name = "araport11_custom",export_bed_tables_TxDb = TRUE,
                                         forge_BSgenome = FALSE,
                                         create_TxDb = TRUE)
  
  expect_that(file.exists(annotation))
  
})

