library(testthat)
library(GenomicFeatures)
library(Rsamtools)





#extending the class::

# setMethod('getSeq',signature=c(e1 = 'Fafile', e2 = 'GRanges'),function(x,gr){
# 	message('foo')
# 	getSeq(x,gr)
# })


#Below is beofre I realized it was a refclass
# FaFile_Circ<-setClass("FaFile_Circ", 
# 	contains="FaFile",
# 	slots=representation(circularRanges='character'),
# 	prototype=list(circularRanges=DEFAULT_CIRC_SEQS)
# 	)

# FaFile_Circ <- function(seqinfo,...){
# 	fafile <- FaFile(...)
# }

# FaFile_Circ(REF)


FaFile_Circ<-setRefClass("FaFile_Circ", 
	contains="FaFile",
	fields=representation(circularRanges='character'),
	prototype=list(circularRanges=DEFAULT_CIRC_SEQS)
)

setMethod(seqinfo,'FaFile_Circ',function (x){
    gr <- scanFaIndex(x)
    Seqinfo(as.character(seqnames(gr)), width(gr), as.character(seqnames(gr)) %in% x$circularRanges )
})



setMethod('getSeq','FaFile_Circ', function (x, ...){
	dots<-list(...)
	if(length(dots)==1 && is(dots[[1]],'GRanges')){
		names(dots)[1] <- 'param'
	}
    .local <- function (x, param, ...){
        if (missing(param)) {
            scanFa(x, ...)
        }
        else {
            if (is(param, "GRanges")) {
            	seqinfo(param) <- seqinfo(x)
                idx <- as.logical(strand(param) == "-")
				chrends <- seqlengths(param)[as.character(param@seqnames)]
              	is_wrapping <- end(param) > chrends
              	is_wrapping <- purrr::map_lgl(is_wrapping,isTRUE)
              	circs_wrap <- isCircular(param)[as.character(seqnames(param[is_wrapping]))] 
              	if(!all( purrr::map_lgl(circs_wrap,isTRUE) )) {
              		stop('Out of bounds ranges detected which are not on circular chromosomes')
              	}
              	#
              	if(any(purrr::map_lgl(is_wrapping,isTRUE))) {
					wrapped_grs <- restrict(param[is_wrapping],chrends+1)
					#get ends for the wrapped ranges
					chrends_wrap <- seqlengths(wrapped_grs)[as.character(wrapped_grs@seqnames)]
					#shift these wrapped ranges back to 1
					wrapped_grs <- shift(wrapped_grs,-chrends_wrap)
					#make sure they don't wrap twice
					if(any(end(wrapped_grs) > chrends_wrap)) stop("Ranges wrapping twice isn't implemented yet...")			
					#also get the within bounds ranges
					nonwrapped_grs <- restrict(param,end=chrends)
					#scan seperately
					dna <- scanFa(x, nonwrapped_grs,...)
					wrap_dna <- scanFa(x, wrapped_grs,...)

					#append the positive wraps to the end of the '+' rnages
					if(any(is_wrapping[!idx])) dna[is_wrapping & (!idx)] <- xscat( dna[is_wrapping & (!idx)], wrap_dna[!idx[is_wrapping]] ) 

					#reverse our negative ranges nad append their wraps to the start
					if(any(idx)){
						dna[idx] <- reverseComplement(dna[idx])
						if(any(is_wrapping[idx])){
							wrap_dna[idx[is_wrapping]] <- reverseComplement(wrap_dna[idx[is_wrapping]])
							dna[is_wrapping & (idx)] <- xscat( wrap_dna[idx[is_wrapping]], dna[is_wrapping & (idx)] ) 
						}
					}
				}else{
					dna <- scanFa(x, param, ...)
					idx <- as.logical(strand(param) == "-")
                	if (any(idx)) dna[idx] <- reverseComplement(dna[idx])
				}
            }else{
            	dna <- scanFa(x, param, ...)
            }
            dna
        }
    }
    do.call(.local,c(list(x=x), dots))
})



testthat::test_that("Test the getSeq methods with out FaFile_circ works properly",{
	

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
	cREF<-FaFile_Circ(Rsamtools::FaFile('test.fa'))
	Rsamtools::indexFa('test.fa')
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
