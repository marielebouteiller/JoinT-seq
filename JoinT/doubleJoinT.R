#!/usr/bin/env Rscript
#repo <- "~/HTGTS_Frocklab/"
repo <- "/home/github_pipeline/HTGTS_Frocklab/"
#repo <- "/Users/marielebouteiller/Dropbox/Marie/Documents_dropbox/HTGTS_Frocklab/"
source(paste(repo, "JoinT/JoinT_functions.R", sep=""))

suppressPackageStartupMessages(library(argparser))
suppressPackageStartupMessages(library(Rbowtie2))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))

p <- arg_parser("Run and then combine RejoinBowtie output for 2 different cutting sites, for ex for ZFN or nickase -mediated breaks")

p <- add_argument(p, "meta", help="metadata file, should contain extra info for the 2nd break/nick site : End_2ndsite or Start_2ndsite, 5nt_BaitEnd_2ndsite and 5nt_PreyStart_2ndsite", type="character")
p <- add_argument(p, "result", help="result folder from wrapper", type="character")
p <- add_argument(p, "preprocess", help="preprocess folder from wrapper", type="character")
p <- add_argument(p, "--which", help="choose library to analyze", type="character")
p <- add_argument(p, "--RepSeq", flag = TRUE, help="treat one specific repeated region in the prey", type="character")
p <- add_argument(p, "--print", flag = TRUE, help="print the libraries in the metadata to choose which ones to process")

argv <- parse_args(p)

doubleJoinT <- function(meta, result_folder, preprocess_folder, ...){
  whichlib <- argv$which
  RepSeq <- ifelse(argv$RepSeq, TRUE, FALSE)
  print(paste("repeated sequence specified : ", ifelse(RepSeq, "yes", "no"), sep=""))
  toprint <- ifelse(argv$print, TRUE, FALSE)
  meta_all <- read.csv(meta, header=TRUE, sep="\t")
  lib_to_process <- seq_len(nrow(meta_all))
  wd <- getwd()
    
  # if print
  if (toprint){
    print(meta_all[,1:2])
    return()
  }
    
  # which libraries to process ?
  if (is.na(whichlib)){
    lib_to_process <- seq_len(nrow(meta_all))
  }
  if (!is.na(whichlib)){
    lib_to_process <- as.integer(unlist(strsplit(whichlib, ",")))
  }
    
  for (i in lib_to_process) {
    lib <- paste(meta_all$Library[i], "_", meta_all$Sequencing[i], sep="")
    print(paste("running doubleRejoinBowtie for ", lib, sep=""))
    metadata <- meta_all[i,]
    
    # create the metadata file for the 2nd cutting site
    metadata_2ndsite <- metadata
    metadata_2ndsite$Start <- metadata_2ndsite$Start_2ndsite
    metadata_2ndsite$End <- metadata_2ndsite$End_2ndsite
    metadata_2ndsite$X5nt_BaitEnd <- metadata_2ndsite$X5nt_BaitEnd_2ndsite
    metadata_2ndsite$X5nt_PreyStart <- metadata_2ndsite$X5nt_PreyStart_2ndsite
    
    master_tlx <- paste(result_folder, "/", lib,"/", lib, ".tlx", sep="")
    result_tlx <- paste(result_folder, "/", lib,"/", lib, "_result.tlx", sep="")
    result <- read.csv(result_tlx, header=TRUE, sep="\t")
    fastq <- paste(preprocess_folder, "/", lib, "_R1.fq", sep="")
      
    dir.create(lib)
    setwd(paste(wd, "/", lib, sep=""))
    
    # create a directory for RejoinBowtie for the 1st cutting site
    dir.create("RejoinBowtie_site1")
    setwd(paste(wd, "/", lib, "/", "RejoinBowtie_site1", sep=""))
    find_rejoining_events(lib, master_tlx, metadata, fastq, RepSeq)
    setwd(paste(wd, "/", lib, sep=""))
    
    # create a directory for RejoinBowtie for the 2nd cutting site
    dir.create("RejoinBowtie_site2")
    setwd(paste(wd, "/", lib, "/", "RejoinBowtie_site2", sep=""))
    find_rejoining_events(lib, master_tlx, metadata_2ndsite, fastq, RepSeq)
    setwd(paste(wd, "/", lib, sep=""))
    
    # create a directory to combine both
    dir.create("RejoinBowtie_combined")
    setwd(paste(wd, "/", lib, "/", "RejoinBowtie_combined", sep=""))
    tlx1 <- paste(wd, "/", lib, "/RejoinBowtie_site1/", lib, "_rejoining.tlx", sep="")
    tlx2 <- paste(wd, "/", lib, "/RejoinBowtie_site2/", lib, "_rejoining.tlx", sep="")
    
    # combine both files of rejoining events
    RB_combine(tlx1, tlx2, lib)
    rejoining_combined <- read.csv(paste(lib, "_rejoining_combined.tlx", sep =""), header=TRUE, sep="\t")
    
    # combine the combined result of rejoining events with the result from wrapper
    combine(rejoining_combined, result, metadata, RepSeq, lib)
    
    # create bedgraph
    if (metadata$Assembly == "mm9" | metadata$Assembly == "mm10"){
      geno <- "mouse"
    }
    if (metadata$Assembly == "hg38" | metadata$Assembly == "hg19"){
      geno <- "human"
    }
    cmd1 = paste("python ", repo,  "experimental/remove_chr.py --filename ", lib, "_result_all.tlx", sep="")
    print(cmd1)
    system(cmd1)
    cmd2 = paste( repo, "post_pipeline/tlx2bed.py -f ", lib, "_result_all_chr.tlx ", "-o ", lib , "_doubleRB", " -g ", geno , " --v3", sep='')
    print(cmd2)
    system(cmd2)
    
    setwd(wd)
  }
  print("all done, enjoy your analysis !")
}


# function to combine the 2 results from RejoinBowtie with 2 sites
# main principle: I keep junctions which are unique from each file. For Qnames common to both files, I keep the alignment for which Qstart - B_Qned is the smallest. If it is the same, I use the alignment from tlx1
RB_combine <- function (tlx1, tlx2, lib){
  n1 <- read.csv(tlx1, header=TRUE, sep="\t")
  n2 <- read.csv(tlx2, header=TRUE, sep="\t")
  n1_unique <- filter(n1, !Qname %in% n2$Qname) %>% arrange(Qname)
  n2_unique <- filter(n2, !Qname %in% n1$Qname) %>% arrange(Qname)
  n1_common <- filter(n1, Qname %in% n2$Qname) %>% arrange(Qname) %>% mutate(jct_space = Qstart - B_Qend)
  n2_common <- filter(n2, Qname %in% n1$Qname) %>% arrange(Qname) %>% mutate(jct_space = Qstart - B_Qend)
  i <- which (n2_common$jct_space < n1_common$jct_space)
  print(paste("Number of rejoining events found by RejoinBowtie with the 1st cutting site : ", nrow(n1), sep=""))
  print(paste("Number of rejoining events found by RejoinBowtie with the 2nd cutting site : ", nrow(n2), sep=""))
  print(paste("Number of Qnames common in both rejoining.tlx files : ",nrow(n1_common), sep=""))
  print(paste("Number of new Qnames added by the second RejoinBowtie : ",nrow(n2_unique), sep=""))
  output_tlx <- n1_common
  for (j in i){
    output_tlx[j,] <- n2_common[j,] 
  }
  output_tlx <- output_tlx %>% select(Qname, JuncID, Rname , Junction, Strand, Rstart, Rend, B_Rname, B_Rstart, B_Rend, B_Strand, B_Qstart, B_Qend, Qstart, Qend, Qlen, 
                                      Seq, filter_step, output)
  
  output_tlx <- rbind(output_tlx, n1_unique, n2_unique)
  write.table(output_tlx, file=paste(lib, "_rejoining_combined.tlx", sep =""), sep="\t", quote = FALSE, row.names = FALSE)
  return(output_tlx)
}

doubleJoinT(argv$meta, argv$result, argv$preprocess, argv$which, argv$locus, argv$print)
