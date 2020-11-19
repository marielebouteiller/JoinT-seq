# HTGTS_Rep_rejoin
One approach to get rejoining events combined to translocations found by wrapper is to run LAM-HTGTS followed by HTGTS_Rep_Rejoin. The HTGTSrep pipeline processed raw sequence reads using the bait sequence and the prey sequence corresponding to the other side of the bait DSB to identify sequence reads harboring proximal bait/prey alignments (Lin et al., 2016). The python script Rejoin is called to further filter sequence reads identified by HTGTSRep that included germline sequence alterations originating from the targeted bait DSB and used a 1nt step, 10bp sliding window along both bait/prey sequences to identify junction positions. Resulting junctions were then combined together with translocations identified separately using TranslocWrapper (Hu et al., 2016).

# Installation
## Get the code
Clone our git rep
```
git clone https://github.com/marielebouteiller/JoinT-seq.git
```

## Requirements
Rejoin requires that LAM-HTGTS and HTGTSrep was run on the libraries

# Rejoin usage
```
python Rejoin.py -id path_to_HTGTSrep_folder -od path_to_output_directory -r path_to_LAM_HTGTS_result_folder -m path_to_metadata_file -g genome
```

## metadata file
The metadata file should be compatible with the TRanslocWrapper and HTGTSrep pipelines

# Rejoin output
* The final output is the Library_Sequencing_PP_bcigar.tlx file
* Rejoin generates bedgraph files for a convenient visualization with IGV
* The file sorted_values.xls gives some statistics on the number of modified and unmodified reads in this library

# References
* Hu J, Meyers RM, Dong J, Panchakshari RA, Alt FW, Frock RL. Detecting DNA double-stranded breaks in mammalian genomes by linear amplification-mediated high-throughput genome-wide translocation sequencing. Nat Protoc. 2016;11(5):853-871. doi:10.1038/nprot.2016.043
* Lin SG, Ba Z, Du Z, Zhang Y, Hu J, Alt FW. 2016. Highly sensitive and unbiased approach for elucidating antibody repertoires. Proc Natl Acad Sci U S A 113: 7846-7851
