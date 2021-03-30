#!/usr/bin/env Rscript
repo <- "/home/HTGTS_Frocklab/"
#repo <- "/home/marielb/HTGTS_Frocklab/"
#repo <- "/Users/marielebouteiller/Dropbox/Marie/Documents_dropbox/HTGTS_Frocklab/"
source(paste(repo, "/JoinT/JoinT_functions.R", sep=""))

suppressPackageStartupMessages(library(argparser))
suppressPackageStartupMessages(library(Rbowtie2))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tryCatchLog))
suppressPackageStartupMessages(library(futile.logger))
suppressPackageStartupMessages(library(foreach))

### Create a parser
p <- arg_parser("JoinT")

### Add command line arguments
p <- add_argument(p, "meta", help="metadata file for the library to analyse", type="character")
p <- add_argument(p, "preprocess", help="preprocess folder from wrapper", type="character")
p <- add_argument(p, "result", help="result folder from wrapper", type="character")
p <- add_argument(p, "output", help="output directory where the results will be written", type="character")
p <- add_argument(p, "--threads", help = "number of threads used to run libraries in parallel", default = 1, type = "integer")
p <- add_argument(p, "--which", help="choose library to analyze", type="character")
p <- add_argument(p, "--RepSeq", flag = TRUE, help="treat one specific repeated region in the prey")
p <- add_argument(p, "--print", flag = TRUE, help="print the libraries in the metadata to choose which ones to process")
p <- add_argument(p, "--keep", flag = TRUE, help="keep some of the intermediate files")

### Parse the command line arguments
argv <- parse_args(p)

# Prepare the work for each library
JoinT <- function(meta, preprocess_folder, result_folder, output_folder, ...){
  whichlib <- argv$which
  RepSeq <- ifelse(argv$RepSeq, TRUE, FALSE)
  keep <- ifelse(argv$keep, TRUE, FALSE)
  toprint <- ifelse(argv$print, TRUE, FALSE)
  meta_all <- read.csv(meta, header=TRUE, sep="\t")
  lib_to_process <- seq_len(nrow(meta_all))
  wd <- output_folder
  setwd(wd)
  
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
  
  # prepare the work to be done in parallel for each library
  cores <- argv$threads
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
 
  # for each library, prepare the files to use, and run the script
  foreach(i = lib_to_process) %dopar% {
    # useful packages and source files for additional functions
    suppressPackageStartupMessages(library(futile.logger))
    repo <- "/home/HTGTS_Frocklab/"
    #repo <- "/Users/marielebouteiller/Dropbox/Marie/Documents_dropbox/HTGTS_Frocklab/"
    #repo <- "/home/marielb/HTGTS_Frocklab/"
    source(paste(repo, "/JoinT/JoinT_functions.R", sep=""))
    
    lib <- paste(meta_all$Library[i], "_", meta_all$Sequencing[i], sep="")
    dir.create(lib)
    setwd(paste(wd, "/", lib, sep=""))
    
    flog.appender(appender.file(paste(lib, ".log", sep="")))
    flog.threshold(INFO)
    flog.info(paste("Running JoinT for ", lib, sep=""))
    flog.info(paste("Working directory : ", wd, "/", lib, sep=""))
    flog.info(paste("Repeated sequence specified : ", ifelse(RepSeq, "yes", "no"), sep=""))
    
    metadata <- meta_all[i,]
    master_tlx <- paste(result_folder, "/", lib,"/", lib, ".tlx", sep="")
    result_tlx <- paste(result_folder, "/", lib,"/", lib, "_result.tlx", sep="")
    result <- read.csv(result_tlx, header=TRUE, sep="\t")
    fastq <- paste(preprocess_folder, "/", lib, "_R1.fq", sep="")
    
    #run the pipeline for the library to find rejoining events
    tlx_rejoining <- find_rejoining_events(lib, master_tlx, metadata, fastq, RepSeq, keep)
    
    # combine the combined result of rejoining events with the result from wrapper
    combine(tlx_rejoining, result, metadata, RepSeq, lib, keep)
    
    # create bedgraph
    cmd1 = paste(repo, "/JoinT/tlx2bed.py -f ", lib, "_JoinT.tlx ", "-o ", lib, " -g ", metadata$Assembly , " --v3", sep='')
    flog.info("Creating bedgraph files")
    flog.info(cmd1)
    tryCatchLog(system(cmd1))
    
    flog.info(paste("JoinT finished for ", lib, sep=""))
    setwd(wd)
  }
  
  # stop the parallel cluster
  parallel::stopCluster(cl)
  warnings()
}

### Do work based on the passed arguments
JoinT(argv$meta, argv$preprocess, argv$result, argv$output, argv$threads, argv$which, argv$RepSeq, argv$print, argv$keep)
