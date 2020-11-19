# Why run JoinT
JoinT-seq combines the identification of a comprehensive map of both re-JOINing and Translocations events (JoinT) arising from a bait DSB. It includes the preparation of libraries using LAM-HTGTS, high-throughput sequencing, preprocessing of the fastq files using the TranslocPreprocess pipeline, and identification of Translocations with the TranslocWrapper pipeline (https://robinmeyers.github.io/transloc_pipeline/thedocs.html). JoinT.R is the module to recover rejoining events and combine them with translocations identified by TranslocWrapper.

# Installation
## Get the code
Clone our git rep
```
git clone https://github.com/marielebouteiller/JoinT-seq.git
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
JoinT.R path_to_metadata.txt path_to_LAM_HTGTS_result_folder path_to_preprocess_folder
```
JoinT uses as arguments a metadata file, together with the folder including the results from TranslocWrapper and the folder including the results from the TranslocPreprocess pipeline (Hu et al., 2016)

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

## LAM-HTGTS results folder
This folder should contain subfolder with names of the individual libraries: Library_Sequencing, as obtained after running the TranslocWrapper.pl pipeline (Hu et al., 2016).

## TranslocPreprocess folder
JoinT needs access in this folder to the Library_Sequencing_R1.fq.gz file for each library to process. These fastq files should be demultiplexed and trimmed, as obtained after running the TranslocPreprocess.pl pipeline (Hu et al., 2016).

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
* A stats file is generated including information about number of junctions found along the script

when the --keep argument is used, the following output files are also generated:
* Library_Sequencing_unmod.tlx, which contains all the reads considered as unmodified (uncut or perfectly rejoined)
* Library_Sequencing_excluded.tlx which contains all the reads which were filtered out because of the quality of some the bases between the end of the bait and start of prey, during the filtering steps
* Library_Sequencing_rejoining.tlx: which contains all the reads considered as rejoining events after the Bowtie2 runs and the fitering steps. To be noted: this is not exactly the rejoining events from the final output Library_Sequencing_JoinT.tlx file, as there might be some re-assignations when combining the rejoining events with the translocations from TranslocWrapper pipeline.
* Library_Sequencing_result_all.tlx: which is similar to the final output Library_Sequencing_JoinT.tlx with all junctions in any chromosomes.

# Reference
Hu J, Meyers RM, Dong J, Panchakshari RA, Alt FW, Frock RL. Detecting DNA double-stranded breaks in mammalian genomes by linear amplification-mediated high-throughput genome-wide translocation sequencing. Nat Protoc. 2016;11(5):853-871. doi:10.1038/nprot.2016.043
