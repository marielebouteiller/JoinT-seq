# JoinT-seq generation
JoinT-seq provides a comprehensive map of prey single break rejoining and genome-wide translocation to a specific bait broken DNA end.
## LAM-HTGTS 
Junction-enriched libraries are generated as previously described for LAM-HTGTS (Hu et al., 2016) but with the removal of the bait DSB rejoining blocking step, and libraries are sequenced. Pooled raw sequences are demultiplexed and adapter trimmed using the TranslocPreprocess pipeline, and translocations are identified using the TranslocWrapper pipeline (https://robinmeyers.github.io/transloc_pipeline/). We used a modified version of the TranslocWrapper pipeline TranslocDedup.R script to correct for the determination of a duplicated junction; the adjustment now additionally takes into account the nucleotide content of inserted sequence between bait/prey and restricts the
indels involved in duplicate determination to those occuring right at the break-site. The modified TranslocDedup.R script is available here and can simply replace the same named script associated with TranslocWrapper (Hu et al., 2016).
## JoinT
The JoinT module was developed to identify bait end rejoining identification and combining with identified translocations from the TranslocWrapper pipeline.
## HTGTS-Rep-Rejoin
As an alternative approach, HTGTS-Rep-Rejoin processes raw sequence reads through the HTGTSrep pipeline (Lin et al., 2016, https://bitbucket.org/adugduzhou/htgtsrep) and the script Rejoin.py identifies rejoining events and combines them with translocations

# References
* Hu J, Meyers RM, Dong J, Panchakshari RA, Alt FW, Frock RL. Detecting DNA double-stranded breaks in mammalian genomes by linear amplification-mediated high-throughput genome-wide translocation sequencing. Nat Protoc. 2016;11(5):853-871. doi:10.1038/nprot.2016.043
* Lin SG, Ba Z, Du Z, Zhang Y, Hu J, Alt FW. 2016. Highly sensitive and unbiased approach for elucidating antibody repertoires. Proc Natl Acad Sci U S A 113: 7846-7851
