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

