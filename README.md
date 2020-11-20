# JoinT-seq generation
JoinT-seq provides a comprehensive map of prey single break rejoining and genome-wide translocation to a specific bait broken DNA end.
## LAM-HTGTS 
Junction-enriched libraries are generated as previously described for LAM-HTGTS (Hu et al., 2016) but with the removal of the bait DSB rejoining blocking step, and libraries are sequenced. Pooled raw sequences are demultiplexed and adapter trimmed using the TranslocPreprocess pipeline, and translocations are identified using the TranslocWrapper pipeline (Hu et al., 2016). From the original TranslocWrapper pipeline, we only modified the TranslocDedup.R script to refine the definition of duplicate junctions : we required the same number of insertion or resection to occur right between the end of the bait and the start of the prey and the exact same bases inserted if any for two junctions to be called duplicates.
## JoinT
The JoinT module was developed to identify bait end rejoining identification and combining with identified translocations from the TranslocWrapper pipeline.
## HTGTS-Rep-Rejoin
As an alternative approach, HTGTS-Rep-Rejoin processes raw sequence reads through the HTGTSrep pipeline (Lin et al., 2016) and the script Rejoin.py identifies rejoining events and combines them with translocations

# References
* Hu J, Meyers RM, Dong J, Panchakshari RA, Alt FW, Frock RL. Detecting DNA double-stranded breaks in mammalian genomes by linear amplification-mediated high-throughput genome-wide translocation sequencing. Nat Protoc. 2016;11(5):853-871. doi:10.1038/nprot.2016.043
* Lin SG, Ba Z, Du Z, Zhang Y, Hu J, Alt FW. 2016. Highly sensitive and unbiased approach for elucidating antibody repertoires. Proc Natl Acad Sci U S A 113: 7846-7851
