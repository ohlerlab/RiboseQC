library(RiboseQC)
library(testthat)

data(res_all)
data(profiles_fivepr)

context("Analysis Working")

test_that("calc_cutoffs_from_profiles works as expected", {

  #using the arabidopsis data, we get a cutoff and high mean_pct_max_frame
  testlen <- '28'
  testdata_periodic <- profiles_fivepr$five_prime_subcodon$nucl[[testlen]]
  
  outcutoff <- calc_cutoffs_from_profiles(testdata_periodic,length_max = 100)
  
  expect_identical(outcutoff$final_cutoff,12)
  expect_gt(outcutoff$frames_res['mean_pct_max_frame'] ,80 )
  
  set.seed(1)
  
  #if we fake aperiodic data, we still get results, but low mean_pct_max_frame
  testdata_not_periodic <- testdata_periodic[,sample(colnames(testdata_periodic))]
  colnames(testdata_not_periodic) <- colnames(testdata_periodic)
  outcutoff_ap <- calc_cutoffs_from_profiles(testdata_not_periodic,length_max = 100)

  expect_lt( outcutoff_ap$frames_res['mean_pct_max_frame'] ,80 )
  
  
  testrl_sel_data <- list(matrix(ncol=2,c(
    outcutoff$final_cutoff,
    outcutoff$frames_res['mean_pct_max_frame'])
    )%>%
    set_colnames(c('cutoff','frame_preference'))%>%
    set_rownames(c(testlen))
  )%>%
  setNames(names(profiles_fivepr[[1]])[1])
  
  len_sel <- choose_readlengths(testrl_sel_data,'max_coverage',profiles_fivepr$five_prime_subcodon['nucl'])
  
  expect_gt(len_sel$nucl$data$gain_codons, 0.60)
  
  
  #Now create fake data for choose_readlengths
  testrl_sel_dataap <- list(matrix(ncol=2,c(
    outcutoff_ap$final_cutoff,
    outcutoff_ap$frames_res['mean_pct_max_frame'])
  )%>%
    set_colnames(c('cutoff','frame_preference'))%>%
    set_rownames(c(testlen))
  )%>%
    setNames(names(profiles_fivepr[[1]])[1])
  
  len_sel <- choose_readlengths(testrl_sel_dataap,'max_coverage',profiles_fivepr$five_prime_subcodon['nucl'])
  
  #test the randomized ata gives no final choice, and negative codon gain
  expect_equal(nrow(len_sel$nucl$final_choice), 0)
  
  expect_lt(len_sel$nucl$data$gain_codons, 0)
  
  
})
