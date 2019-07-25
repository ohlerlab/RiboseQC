library(testthat)


testthat::test_that("Test the getSeq methods with out FaFile_circ works properly",{
  #library(RiboseQC)

  #First veryify this functionality is missing

  # REF<-Rsamtools::FaFile('pipeline/my_hg19.fa')
  REF<-'test.fa'
  writeXStringSet(setNames(DNAStringSet('AAAAAAAAGG'),'chrM'),filepath=REF)
  indexFa(REF)
  seqinfo(FaFile(REF))

  #showMethod doesn't actually show code (just the package it's in)

  #selectMethodd does though
  #selectMethod('seqinfo',signature='FaFile')

  cREF <- FaFile_Circ(FaFile(REF))
  cREF <- FaFile_Circ(FaFile(REF),circularRanges='chrM')

  expect_true(isCircular(seqinfo(cREF)))
  grs <- GRanges(c('chrM:4-6:+','chrM:8-11:+','chrM:8-12:-'))


  #So indeed even if I change the seqinfo function - I have to change.scan fa
  getSeq(cREF,grs)

  selectMethod('getSeq','FaFile')
  showMethods(Rsamtools::scanFa,c('FaFile','GRanges'))
  selectMethod(Rsamtools::scanFa,c('FaFile','GRanges'))
  #this gets me to this C pointer... so I need to stop here


  expect_equal(getSeq(cREF,grs[1]),setNames(DNAStringSet('AAA'),'chrM'))
  expect_warning(getSeq(cREF,grs[1]),regexp=NA)

  isCircular(grs@seqinfo)<-TRUE
  seqlengths(grs@seqinfo)<-10

  writeXStringSet(setNames(DNAStringSet(c('AAAAAAAAGG','AAAAAAAAGG')),c('chrM','chr2')),filepath='test.fa')
  Rsamtools::indexFa('test.fa')

  cREF<-FaFile_Circ(Rsamtools::FaFile('test.fa'),circularRanges='chrM')

  expect_true(isCircular(seqinfo(cREF))['chrM'])
  expect_false(isCircular(seqinfo(cREF))['chr2'])


  expect_equal(as.character(getSeq(cREF)),setNames(c('AAAAAAAAGG','AAAAAAAAGG'),c('chrM','chr2')))


  grs <- GRanges(c('chrM:4-6:+','chrM:8-11:+','chrM:8-12:-','chr2:8-12'))
  isCircular(grs@seqinfo)['chrM']<-TRUE
  seqlengths(grs@seqinfo)[]<-c(10)
  expect_error(getSeq(cREF,grs),'not on circ')

  expect_equal(as.character(getSeq(cREF,(grs[1:3]))),setNames(c('AAA','AGGA','TTCCT'),c('chrM','chrM','chrM')))


  testcirctrs <- GRangesList(list(GRanges(c('chrM:4-6:+','chrM:8-11:+')),GRanges(c('chrM:8-11:-','chrM:4-6:-'))))
  isCircular(seqinfo(testcirctrs))<-TRUE
  seqlengths(testcirctrs)<-10

  #this doesn't work, but at least the seqinfo method is being accessed
  expect_error(extractTranscriptSeqs(FaFile(REF),testcirctrs),'truncated')

  #Check we can retrieve the transcript sequences properly
  expect_equal(extractTranscriptSeqs(cREF,testcirctrs),DNAStringSet(c('AAAAGGA','TCCTTTT')))

})
