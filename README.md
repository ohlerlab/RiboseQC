# Ribo-seQC
An analysis tool for Ribo-seq and small RNA-seq data


**Ribo-seQC** (*RiboseQC*) is an R package (to be submitted to *Bioconductor*) that performs quality control analysis of small RNA-seq data, with a focus on Ribo-seq and related techniques. Thanks to syntax and functions present in *Bioconductor* packages like *GenomicFeatures*, *rtracklayer* or *BSgenome*, this package can perform comprehensive analyses on a variety of genomic regions. In addition, Ribo-seQC allows to automatically generage an html report for each analyzed sample, allowing for quick and interactive comparison of multiple samples at once.

This tools focuses on the analysis of different read lengths, taking into account the genomic regions they map to (e.g. coding sequence, UTRs, non-coding RNAs, mitochondria or chloroplasts, etc...). Other useful fetures, such as automatic P-sites position calculation or analysis the top mapping positions, are available in the Ribo-seQC package, and we encourage to donwload and have a look at the vignette (**RiboseQC.html**), our manual (**RiboseQC-manual.pdf**), and our manuscript (to be added soon...).


To install the package:
```
library("devtools")
install_github(repo = "lcalviell/Ribo-seQC")
library("RiboseQC")

```

Two simple steps are required to use Ribo-seQC on your data:
```
?prepare_annotation_files
```
parses a *.gtf* and a *.2bit* file. (once per annotation-genome combination, a .2bit file can be obtained froma fasta file using the *faToTwoBit* software from UCSC: https://genome.ucsc.edu/goldenpath/help/twoBit.html - http://hgdownload.soe.ucsc.edu/admin/exe/ )


and
```
?RiboseQC_analysis
```

the master function used to perform the entire analysis workflow.
Please check the vignette for some example workflows.


For any question, please email:

calviello.l.bio@gmail.com, dominique.sydow@posteo.de (analysis and data visualization), Dermot.Harnett@mdc-berlin.de (package mantainer), Uwe.ohler@mdc-berlin.de (project supervisor.)


Enjoy!


