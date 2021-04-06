# Why run JoinT
JoinT-seq combines the identification of a comprehensive map of both re-JOINing and Translocations events (JoinT) arising from a bait DSB. It includes the preparation of libraries using LAM-HTGTS, high-throughput sequencing, preprocessing of the fastq files using the TranslocPreprocess pipeline, and identification of Translocations with the TranslocWrapper pipeline (https://robinmeyers.github.io/transloc_pipeline/thedocs.html). JoinT.R is the module to recover rejoining events and combine them with translocations identified by TranslocWrapper.

# Installation
## Get the code
Clone our git rep
```
git clone 
```

## Dependencies
JoinT.R a version of R >=3.6 (https://cran.r-project.org), with the following packages : argparser, stringr, dplyr, Rbowtie2, ShortRead, which can be installed as follows
```
install.packages(c("argparser", "stringr","dplyr"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Rbowtie2","ShortRead"))
```

## Location of your github repository
In the JoinT.R script, you need to adjust the location of the folder JoinT on your machine

## Genome assembly
For an easy visualisation of the output junctions, JoinT generates bedgraph files to be opened with softwares like IGV. The code runs now to generates bedgraph files for mm9, mm10, hg19 and hg38 assembly. Any new assembly can be added with the addition in the JoinT folder of the proper file with adapted information about the chromosomes of this assembly (chromInfo_assembly.txt)

# JoinT usage
```
JoinT.R path_to_metadata.txt path_to_preprocess_folder path_to_LAM_HTGTS_result_folder path_to_output_directory 
```
JoinT uses as arguments a metadata file, together with the folder including the results from the TranslocPreprocess pipeline, the folder including the results from TranslocWrapper (Hu et al., 2016) and the desired folder for the output

## Metadata
The metadata file should be compatible with TranslocWrapper pipeline and include the required following columns: Library, Sequencing, Assembly, Chr, Start, End, Strand, MID, Primer, Adapter, 5nt_BaitEnd, 5nt_PreyStart, Amplicon
* Library : ID of your library
* Sequencing: name of the sequencing run
* Assembly: genome assembly used for the coordinates and for the identification of translocations
* Chr: chromosome of the DSB
* Start and End: refer to the Start and End position of the bait, the bait being defined by the sequence between the nested-primer and the DSB
* Strand: Strand of the bait
* MID: barcode of the library
* 5nt_BaitEnd: 5 last nucleotides of the bait
* 5nt_PreyStart: 5 first nucleotides of the prey
* Amplicon: sequence including the bait and the prey. The length of the prey should be adapted to extensively include all rejoining events

Optional columns :
* Description : optional
* Repeat_Start
* Repeat_End : to be used with the optional argument --RepSeq, see below

## TranslocPreprocess folder
JoinT needs access in this folder to the Library_Sequencing_R1.fq.gz file for each library to process. These fastq files should be demultiplexed and trimmed, as obtained after running the TranslocPreprocess.pl pipeline (Hu et al., 2016). 
Note : **They should have the extension .fq (or .fq.gz if gzipped)** and not .fastq. If they have the extension .fq.gz, the .gz extension should come from compression gzip and not from the name. When fq files are normalized using seqtk, this generates unzipped fq files which shouldn't be called .fq.gz.

## LAM-HTGTS results folder
This folder should contain subfolder with names of the individual libraries: Library_Sequencing, as obtained after running the TranslocWrapper.pl pipeline (Hu et al., 2016).

## Output folder
This is the path of the folder where the output will be written, one sub-folder per library will be created in this folder

# Optional arguments

## --which
When called, JoinT.R will process all libraries in the metadata file, unless the optional argument --which specifies which libraries to run. The number of the library to run from the metadata file can be specified.

## --print
When called, JoinT.R will only show the libraries from the metadata with their associated numbers to call with the –which argument if needed.

## --RepSeq
This argument is used to accommodate the presence of repeated sequences within the prey region. It requires to add 2 extra columns in the metadata file, named Repeat_Start and Repeat_End with the position of respectively the start and the end of the repeated region.

## --keep 
This argument is used to keep some of the intermediate files generated

# JoinT output
* the final output is the Library_Sequencing_JoinT.tlx file, which combines the information about all translocations and rejoining events found in main chromosomes.
* JoinT generates bedgraph files, for an easy visualization with IGV.
* A stats file is generated including information about number of junctions found along the script.
* A log file for each library to check the different steps.

when the --keep argument is used, the following output files are also generated:
* Library_Sequencing_unmod.tlx, which contains all the reads considered as unmodified (uncut or perfectly rejoined)
* Library_Sequencing_excluded.tlx which contains all the reads which were filtered out because of the quality of some the bases between the end of the bait and start of prey, during the filtering steps
* Library_Sequencing_rejoining.tlx: which contains all the reads considered as rejoining events after the Bowtie2 runs and the fitering steps. To be noted: this is not exactly the rejoining events from the final output Library_Sequencing_JoinT.tlx file, as there might be some re-assignations when combining the rejoining events with the translocations from TranslocWrapper pipeline.
* Library_Sequencing_result_all.tlx: which is similar to the final output Library_Sequencing_JoinT.tlx with all junctions in any chromosomes.

# Reference
Hu J, Meyers RM, Dong J, Panchakshari RA, Alt FW, Frock RL. Detecting DNA double-stranded breaks in mammalian genomes by linear amplification-mediated high-throughput genome-wide translocation sequencing. Nat Protoc. 2016;11(5):853-871. doi:10.1038/nprot.2016.043

# R session info
```
> sessionInfo()
R version 4.0.3 (2020-10-10)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /usr/lib64/libopenblasp-r0.3.3.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] foreach_1.5.1               futile.logger_1.4.3        
 [3] tryCatchLog_1.2.1           stringr_1.4.0              
 [5] dplyr_1.0.2                 ShortRead_1.48.0           
 [7] GenomicAlignments_1.26.0    SummarizedExperiment_1.20.0
 [9] Biobase_2.50.0              MatrixGenerics_1.2.0       
[11] matrixStats_0.57.0          Rsamtools_2.6.0            
[13] GenomicRanges_1.42.0        GenomeInfoDb_1.26.2        
[15] Biostrings_2.58.0           XVector_0.30.0             
[17] IRanges_2.24.0              S4Vectors_0.28.1           
[19] BiocParallel_1.24.1         BiocGenerics_0.36.0        
[21] Rbowtie2_1.12.0             argparser_0.6              

loaded via a namespace (and not attached):
 [1] formatR_1.7            pillar_1.4.7           compiler_4.0.3        
 [4] RColorBrewer_1.1-2     iterators_1.0.13       futile.options_1.0.1  
 [7] bitops_1.0-6           tools_4.0.3            zlibbioc_1.36.0       
[10] tibble_3.0.4           lifecycle_0.2.0        lattice_0.20-41       
[13] pkgconfig_2.0.3        png_0.1-7              rlang_0.4.9           
[16] Matrix_1.2-18          DelayedArray_0.16.0    GenomeInfoDbData_1.2.4
[19] hwriter_1.3.2          generics_0.1.0         vctrs_0.3.5           
[22] tidyselect_1.1.0       grid_4.0.3             glue_1.4.2            
[25] R6_2.5.0               jpeg_0.1-8.1           latticeExtra_0.6-29   
[28] lambda.r_1.2.4         purrr_0.3.4            magrittr_2.0.1        
[31] codetools_0.2-18       ellipsis_0.3.1         stringi_1.5.3         
[34] RCurl_1.98-1.2         crayon_1.3.4  
```
