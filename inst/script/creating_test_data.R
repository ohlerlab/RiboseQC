library(tidyverse)
library(magrittr)
library(rtracklayer)
library(GenomicFeatures)
library(Rsamtools)
library(RiboseQC)
#split chrs into compartments
#using the counts from the riboseqc_analysis function
data(res_all)
chrs <- res_all$read_stats$counts_cds_genes$chromosome%>%unique
complist <- lapply(list('Chr[^CM]','ChrM','ChrC'),function(reg)str_subset(chrs,reg))
#then for each comp get the top 100 genes
comptopgenedfs<-lapply( complist,function(compchrs){
  cdsgenedf<-res_all$read_stats$counts_cds_genes%>%subset(chromosome %in% compchrs)
  topcdsinds <- cdsgenedf%>%.$reads%>%{which(rank(- . )<=20)}
  cdsgenedf[topcdsinds,]
})
#combine
comptopgenedfs%<>%lapply(.%>%as.data.frame%>%rownames_to_column('gene_id'))
comptopgenedfs%<>%bind_rows

#get all the gtf data
allarabanno<-import('../ex/test_arabidopsis.gtf.gz')

#filter all the gtf data with our top genes
filtarabanno<-allarabanno%>%subset(gene_id%in%comptopgenedfs$gene_id)

#get our fasta sequence
myarabfa <- '../ex/test_arabidopsis.fa'
myarabseq <- getSeq(x=FaFile(myarabfa))

#filter it and write it
simplifiedmyarabseq <- myarabseq
simplifiedmyarabseq[filtarabanno%>%add(100)%>%coverage%>%not] <- 'N'
simplifiedmyarabseq
simplifiedmyarabseq%>%writeXStringSet('../ex/simp_arab.fa')
#zip this
system('gzip ../ex/simp_arab.fa')
myarabfa <- '../ex/simp_arab.fa.gz'
#copy to ext data
system('cp ../ex/simp_arab.fa.gz inst/ext_data/')

#works
stopifnot(getSeq(x=FaFile(myarabfa))%>%is('DNAStringSet'))

#copy out the filtered gtf data
filtarabanno$score=0
export(filtarabanno,'../ex/simp_arab.gtf')

#zip our gtf
system('gzip ../ex/simp_arab.gtf')
system('cp ../ex/simp_arab.gtf inst/ext_data/')

#now filter the bams
readGAlignments('../ex/test_arabidopsis_root.bam',param = ScanBamParam(what='mapq'))%>%
    subsetByOverlaps(filtarabanno)%>%
    export('inst/ext_data/simp_arab_root.bam')

readGAlignments('../ex/test_arabidopsis_shoot.bam',param = ScanBamParam(what='mapq'))%>%
  subsetByOverlaps(filtarabanno)%>%
  export('inst/ext_data/simp_arab_shoots.bam')


