#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparser))
suppressPackageStartupMessages(library(Rbowtie2))
suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))

### function to get rejoining events
find_rejoining_events <- function(lib, master_tlx, meta, fastq, RepSeq, keep){
  print(Sys.time())
  
  # Read the master tlx
  master <- read.csv(master_tlx, header=TRUE, colClasses=c("character", "integer", rep("NULL",13), "integer", rep("NULL",2), "character", rep("NULL",13)), sep="\t")
  bait_length <- as.integer(meta$End - meta$Start)
  Amplicon <- DNAString(meta$Amplicon)
  bait_seq <- Amplicon[1:bait_length]
  prey_seq <- Amplicon[(bait_length + 1):length(Amplicon)]
  seq_breaksite <- paste(meta$X5nt_BaitEnd, meta$X5nt_PreyStart, sep="")
  
  # Collect all sequences from unique Qname from the master tlx file and convert them into a fasta file
  master_unique <- filter(master, JuncID==1)
  rm(master)
  ds_set <- DNAStringSet(master_unique$Seq)
  bs <- BStringSet(master_unique$Qname)
  names(ds_set) <- bs
  master_fa_name <- paste(lib, ".fa", sep="")
  writeFasta(ds_set, master_fa_name, mode="w", width=20000L)
  
  # Build the bowtie index for the prey from the metadata and run bowtie against the prey sequence, version to keep several alignments
  prey_fa <- paste(">prey","\n",prey_seq, sep="")
  write(prey_fa, "prey.fa", sep="")
  bowtie2_build("prey.fa", "prey", overwrite=TRUE)
  prey_sam_name <- paste(lib, "_prey.sam", sep="")
  bowtie2("prey", prey_sam_name, seq1 = master_fa_name , "-f --local -D 20 -R 3 -N 0 -L 8 -i C,6 --ma 2 --mp 8,2 --np 1 --rdg 8,4 --rfg 8,4 --score-min C,40 -p 1 --no-unal --reorder -t -k 3", seq2 = NULL, overwrite=TRUE)
  print(paste("bowtie2 parameters for prey alignemnt : ", "-f --local -D 20 -R 3 -N 0 -L 8 -i C,6 --ma 2 --mp 8,2 --np 1 --rdg 8,4 --rfg 8,4 --score-min C,40 -p 1 --no-unal --reorder -t -k 3", sep=""))
  
  # I remove the alignement in negative orientation i.e. V2 == 16 or 256
  sam_prey <- read.csv(prey_sam_name, skip=3, sep="\t",col.names = c("Qname","Flag","", "Pos","","Cigar", rep("",5),"Score", rep("",8)), colClasses=c("character", "integer","NULL","integer","NULL","character", rep("NULL", 5), "character", rep("NULL",8)), header=FALSE) %>% filter(Flag != 16 & Flag != 272)
  if (!RepSeq){
    sam_prey <- sam_prey %>% group_by(Qname) %>% slice_min(Pos, with_ties=FALSE) # the first alignment on the sequnce is kept
  }
  if (RepSeq){
    if (is.null(meta$Repeat_Start) | is.null(meta$Repeat_End)){
      print("Check your metadata ! if you want to use the RepSeq argument, you need to add a Repeat_Start and a Repeat_End column !")
      return()
    }
    if (meta$Strand == "+"){
      rep_start <-  meta$Repeat_Start - meta$End
      rep_end <-  meta$Repeat_End - meta$End
    }
    if (meta$Strand == "-"){
      rep_start <-  meta$Start - meta$Repeat_End
      rep_end <-  meta$Start - meta$Repeat_Start
    }
    sam_prey <- filter(sam_prey, !(Flag ==256 & Pos >= (rep_start - 1) & Pos <= (rep_end - 10) )) %>% group_by(Qname) %>% slice_min(Pos, with_ties=FALSE) # I remove lines with Flag==256 and Pos between 34 and 70 i.e. lines with Position in the CACA repeat which are 2ndary alignment
  }
  ds_set_with_prey <- ds_set[which(names(ds_set) %in% sam_prey$Qname)]
  with_prey_fa_name <- paste(lib,"_with_prey.fa", sep="")
  writeFasta(ds_set_with_prey, with_prey_fa_name, mode="w", width=20000L)
  rm(ds_set, ds_set_with_prey, bs)
  
  # Build the bowtie index for the bait from the metadata and run bowtie against the bait sequence
  bait_fa <- paste(">bait","\n",bait_seq, sep="")
  write(bait_fa, "bait.fa", sep="")
  bowtie2_build("bait.fa", "bait", overwrite=TRUE)
  bait_sam_name <- paste(lib, "_bait.sam", sep="")
  bowtie2("bait", bait_sam_name, seq1 = with_prey_fa_name , "-f --local -D 20 -R 3 -N 0 -L 10 -i C,6 --ma 2 --mp 8,2 --np 1 --rdg 8,4 --rfg 8,4 --score-min C,40 -p 1 --no-unal --reorder -t ", seq2 = NULL, overwrite=TRUE)
  print(paste("bowtie2 parameters for bait alignemnt : ", "-f --local -D 20 -R 3 -N 0 -L 10 -i C,6 --ma 2 --mp 8,2 --np 1 --rdg 8,4 --rfg 8,4 --score-min C,40 -p 1 --no-unal --reorder -t ", sep=""))
  sam_bait <-  read.csv(bait_sam_name, skip=3, sep="\t", col.names = c("Qname","Flag","", "Pos","","Cigar", rep("",5),"Score", rep("",8)), colClasses=c("character", "integer","NULL","integer","NULL","character", rep("NULL", 5), "character", rep("NULL",8)), header=FALSE) %>% arrange(Qname) %>% filter(Flag != 16 & Flag != 272)
  sam_prey_filtered1 <- data.frame(Qname = sam_prey$Qname, Flag=sam_prey$Flag, Pos=sam_prey$Pos, Cigar=sam_prey$Cigar, Score=sam_prey$Score)
  sam_prey_filtered2 <- filter(sam_prey_filtered1, Qname %in% sam_bait$Qname) %>% arrange(Qname)
  
  # Collect the Seq and the Qlen from the master
  master_fil <- filter(master_unique, Qname %in% sam_bait$Qname) %>% arrange(Qname)
  rm(master_unique)
  
  # Create a df to collect the results, and calculate new info
  df <- data.frame(Qname=sam_prey_filtered2$Qname, position_prey = sam_prey_filtered2$Pos, Cigar_prey_alignement=sam_prey_filtered2$Cigar, position_bait = sam_bait$Pos, Cigar_bait_alignement = sam_bait$Cigar)
  rm(sam_prey, sam_prey_filtered1, sam_prey_filtered2, sam_bait)
  if (nrow(df) == 0){
    print("Check your metadata ! you do not have any read which align to both your bait and you prey: are the columns Chr, Start, End, Strand, Primer, 5nt_BaitEnd, 5nt_PreyStart, and Amplicon consistent ? I am sorry the script will crash :/")
  }
  print(Sys.time())
  df$Seq = master_fil$Seq
  df$Qlen = master_fil$Qlen
  rm(master_fil)
  df$Amplicon_B_Rstart = df$position_bait
  df$Amplicon_Rstart = df$position_prey + bait_length
  df_bait_info <- as.data.frame(t(sapply(df$Cigar_bait_alignement, cigar_info)))
  df$B_Qstart = df_bait_info$V1
  df$B_Qend = df_bait_info$V2
  df$Amplicon_B_Rend = df$Amplicon_B_Rstart + df_bait_info$V3 - 1
  rm(df_bait_info)
  df_prey_info <- as.data.frame(t(sapply(df$Cigar_prey_alignement, cigar_info)))
  df$Qstart = df_prey_info$V1
  df$Qend = df_prey_info$V2
  df$Amplicon_Rend = df$Amplicon_Rstart + df_prey_info$V3 - 1
  rm(df_prey_info)
    
  # remove columns I do not need any more
  df <- df %>% select(Qname, Seq, Qlen, B_Qstart, B_Qend, Qstart, Qend, Amplicon_B_Rstart, Amplicon_B_Rend, Amplicon_Rstart, Amplicon_Rend)     

  if (meta$Strand =="+") { ### I do not really need these columns right now for the filter steps, so might be possible to create them later
    df$B_Rstart = meta$Start + df$Amplicon_B_Rstart - 1
    df$B_Rend = meta$Start  + df$Amplicon_B_Rend - 1 
    df$Rstart = meta$Start + df$Amplicon_Rstart - 1
    df$Rend = meta$Start + df$Amplicon_Rend - 1
    df$Junction = df$Rstart
  }
  if (meta$Strand =="-") {
    df$B_Rend = meta$End - df$Amplicon_B_Rstart
    df$B_Rstart = meta$End - df$Amplicon_B_Rend
    df$Rstart = meta$End - df$Amplicon_Rend 
    df$Rend = meta$End - df$Amplicon_Rstart
    df$Junction = df$Rend
  }
  
  #write.table(df, file=paste(lib, "_df.tlx", sep =""), sep="\t", quote = FALSE, row.names = FALSE)
  print("Bowtie run and info collection done, starting the filter steps")
  
  # Filtering steps
  df$filter_step <- ""
  
  # cases Amplicon_Rstart - Amplicon_B_Rend == 1 & Qstart - B_Qend ==1 : perfect uncut or rejoining
  n_uncut <- nrow(filter(df, Amplicon_Rstart - Amplicon_B_Rend == 1 & Qstart - B_Qend ==1))
  df %>% filter(Amplicon_Rstart - Amplicon_B_Rend == 1 & Qstart - B_Qend ==1) %>% 
    mutate(filter_step = "perfectly uncut",
           JuncID = 1,
           Rname = meta$Chr,
           Strand = ifelse(meta$Strand == "+",1,-1),
           B_Rname = meta$Chr,
           B_Strand = ifelse(meta$Strand == "+",1,-1),
           output = "uncut") %>%
    select(Qname, JuncID, Rname , Junction, Strand, Rstart, Rend, B_Rname, B_Rstart, B_Rend, B_Strand, B_Qstart, B_Qend, Qstart, Qend, Qlen, Seq, filter_step, output) %>%
    write.table(file=paste(lib, "_unmod.tlx", sep =""), sep="\t", quote = FALSE, row.names = FALSE)
  
  # I remove the "clean" "uncut or perfectly rejoined from df to keep less variables and less memory
  df <- filter(df, !(Amplicon_Rstart - Amplicon_B_Rend == 1 & Qstart - B_Qend ==1))
  
  # cases insertion or deletion right at the breaksite and then direct junction, or 1 mismatch at the breaksite
  df <- df %>% 
    mutate(filter_step = if_else(Amplicon_Rstart - Amplicon_B_Rend == 1 & Qstart - B_Qend > 1, "insertion only", filter_step),
           filter_step = if_else(Amplicon_Rstart - Amplicon_B_Rend > 1 & Qstart - B_Qend == 1, "resection only", filter_step),
           filter_step = if_else(Amplicon_Rstart - Amplicon_B_Rend == Qstart - B_Qend & Qstart - B_Qend == 2 , "1 mismatch at the breaksite", filter_step)
    )
  
  # if Qstart <= B_Qend : MH meaning I had deletion before the junction
  df <- df %>% mutate(filter_step = if_else(Qstart <= B_Qend , "deletion leading to MH in the sequence bait/prey", filter_step))
  
  # case where I have 2 bases in between the bait and the prey on both the read and the amplicon i.e. 1 or 2 mismatches at the breaksite
  n9 <- filter(df, Amplicon_Rstart - Amplicon_B_Rend == Qstart - B_Qend & Qstart - B_Qend == 3)
  if (nrow(n9) > 0){
    for (l in (1:nrow(n9))){
      mm_ratio <- pp_mismatch(n9[l,]$B_Qend, n9[l,]$Qstart, n9[l,]$Amplicon_B_Rend, n9[l,]$Amplicon_Rstart, Amplicon, n9[l,]$Seq)
      n9[l,]$filter_step <- if_else(mm_ratio == 1, "2 mismatches at the breaksite", "1 mismatch out of 2 at the breaksite")
    }
  }
  n10 <- filter(n9, filter_step == "2 mismatches at the breaksite")
  
  # read the fastq file
  cmd_gunzip <- paste("gunzip ",fastq, ".gz", sep="")
  print(cmd_gunzip)
  system(cmd_gunzip)
  fastq_R1 <- readFastq(fastq)
  cmd_gzip <- paste("gzip ", fastq, sep="")
  print(cmd_gzip)
  system(cmd_gzip)
  
  if (nrow(n10) >0){
    n10$mm1 <- NA
    n10$mm2 <- NA
    n10$mm_min <- NA
    for (l in (1:nrow(n10))){
      i <- grep(paste(n10[l,]$Qname, " ", sep=""), ShortRead::id(fastq_R1), ignore.case = TRUE)
      qual <- as(quality(fastq_R1)[i], "matrix")
      if ((n10[l,]$B_Qend + 2) > length(qual)){
        n10[l,]$filter_step <- "2 mismatches at the breaksite with poor quality"
      }
      else {
        n10[l,]$mm1 <- qual[1,(n10[l,]$B_Qend + 1)]
        n10[l,]$mm2 <- qual[1,(n10[l,]$B_Qend + 2)]
        n10[l,]$mm_min <- min(n10[l,]$mm1, n10[l,]$mm2)
        n10[l,]$filter_step <- if_else(n10[l,]$mm_min >20 , "2 mismatches at the breaksite with good quality", "2 mismatches at the breaksite with poor quality")
      }
    }
  }
  
  print("n9 and n10 done")
  
  # case where I have the same number of bases (>2) in the read and in the amplicon but with some mismatches,
  # if the number of bases between B_Qend and Qstart is not too high, then I count the ratio of mismatches. if it is >0.7, then I consider it as a rejoining event (n1)
  # i.e 3# in 3 bases 3# (or more) in 4 bases, 4# (or more) in 5, 5# (or more) in 6 bases, 5# (or more) in 7 bases and 6# (or more) in 8 bases
  # if the number of bases between B_Qend and Qstart is high, then I want to check right at the breaksite, if TGATAGTAGG in the sequence, this will be treated later (n4)
  # (I have some cases where I will have an insertion/deletion at the breaksite, and then a deletion/insertion elsewhere, so I cannot find the TGATAGTAGG -> it is considered as a rejoining event, even if the deletion/insertion is up to 4 bases away from the breaksite)
  
  n1 <- filter(df, Amplicon_Rstart - Amplicon_B_Rend == Qstart - B_Qend & Qstart - B_Qend > 3 & Qstart - B_Qend < 9)
  if (nrow(n1) > 0){
    n1$mm <- NA
    for (l in (1:nrow(n1))){
      mm_ratio <- pp_mismatch(n1[l,]$B_Qend, n1[l,]$Qstart, n1[l,]$Amplicon_B_Rend, n1[l,]$Amplicon_Rstart, Amplicon, n1[l,]$Seq)
      n1[l,]$mm <- mm_ratio
    }
    n1 <- n1 %>% mutate(filter_step = if_else(mm <= 0.70, "few mismatches at the breaksite", "many mismatches at the breaksite"))
  }
  
  n11 <- filter(n1, filter_step == "many mismatches at the breaksite")
  if (nrow(n11) >0){
    n11 <- qual_average(n11, fastq_R1)
  }
  n11$filter_step <- if_else(n11$mean < 25, "many mismatches at the breaksite with poor quality", "many mismatches at the breaksite with good quality")
  
  # cases Amplicon_Rstart - Amplicon_B_Rend != Qstart - B_Qend  & Qstart - B_Qend !=1 & Amplicon_Rstart - Amplicon_B_Rend != 1 and Qstart > B_Qend
  # it means I have insertion or deletion, but it can be far away from the breaksite
  # I keep junctions if I have deletion/insertions less or = 2 bases away from the breaksite
  df <- df %>% 
    mutate(filter_step = if_else(Qstart - B_Qend == 2 & Amplicon_Rstart - Amplicon_B_Rend != 1 & Amplicon_Rstart - Amplicon_B_Rend != 2, "resection (1 or more) 1 base away from the break or mismatch + resection", filter_step),
           filter_step = if_else(Amplicon_Rstart - Amplicon_B_Rend == 2 & Qstart - B_Qend > 2, "insertion (1 or more) 1 base away from the break or mismatch + insertion", filter_step))
  
  # if I have deletion/insertion 3 bases away from the breaksite, i.e. Amplicon_Rstart - Amplicon_B_Rend == 4 and Qstart - B_Qend > 4 or vice versa
  # then I keep the junction if I have more than 1 deletion or insertion, or I have only 1 deletion/insertion but also at least 1 mismatch i.e. if the sequence inserted in the query and ref matches
  df <- df %>% 
    mutate(filter_step = if_else(Amplicon_Rstart - Amplicon_B_Rend == 3 & Qstart - B_Qend > 4, "more (strictly) than 1 insertion 2 bases or less away from the breaksite", filter_step),
           filter_step = if_else(Amplicon_Rstart - Amplicon_B_Rend > 4 & Qstart - B_Qend == 3, "more (strictly) than 1 deletion 2 bases or less away from the breaksite", filter_step))
  
  n7 <- filter(df, Amplicon_Rstart - Amplicon_B_Rend == 3 & Qstart - B_Qend == 4)
  if (nrow(n7) > 0){
    for (l in (1:nrow(n7))){
      if (length(matchPattern(ins_amplicon(n7[l,], Amplicon),ins(n7[l,]))) == 0){ # ref sequences is the amplicon, as I have insertion in this case
        n7[l,]$filter_step <- "1 insertion 2 bases max away from the breaksite and at least one mismatch"
      }
      else {
        n7[l,]$filter_step <- "only 1 insertion 2 bases away from the breaksite"
      }
    }
  }
  
  n8 <- filter(df, Amplicon_Rstart - Amplicon_B_Rend == 4 & Qstart - B_Qend == 3)
  if (nrow(n8) > 0){
    for (l in (1:nrow(n8))){
      if (length(matchPattern(ins(n8[l,]),ins_amplicon(n8[l,], Amplicon))) == 0){ # ref sequences is the seq and not the amplicon, as I have deletion in this case
        n8[l,]$filter_step <- "1 deletion 2 bases max away from the breaksite and at least one mismatch"
      }
      else {
        n8[l,]$filter_step <- "only 1 deletion 2 bases away from the breaksite"
      }
    }
  }
  print("n1, n7, n8 done")
  
  # if I have deletion/insertion 3 bases away from the breaksite, i.e. Amplicon_Rstart - Amplicon_B_Rend == 4 and Qstart - B_Qend > 4 or vice versa
  # then I keep the junction if I have more than 1 deletion or insertion, or I have only 1 deletion/insertion and also at least 1 mismatch i.e. if the sequence inserted in the query and ref matches
  df <- df %>% 
    mutate(filter_step = if_else(Amplicon_Rstart - Amplicon_B_Rend == 4 & Qstart - B_Qend > 5, "more (strictly) than 1 insertion 3 bases or less away from the breaksite", filter_step),
           filter_step = if_else(Amplicon_Rstart - Amplicon_B_Rend > 5 & Qstart - B_Qend == 4, "more (strictly) than 1 deletion 3 bases or less away from the breaksite", filter_step))
  
  n2 <- filter(df, Amplicon_Rstart - Amplicon_B_Rend == 4 & Qstart - B_Qend == 5)
  if (nrow(n2) >0) {
    for (l in (1:nrow(n2))){
      if (length(matchPattern(ins_amplicon(n2[l,], Amplicon),ins(n2[l,]))) == 0){ # ref sequences is the amplicon, as I have insertion in this case
        n2[l,]$filter_step <- "1 insertion 3 bases max away from the breaksite and at least one mismatch"
      }
      else {
        n2[l,]$filter_step <- "only 1 insertion 3 bases away from the breaksite"
      }
    }
  }
  
  n3 <- filter(df, Amplicon_Rstart - Amplicon_B_Rend == 5 & Qstart - B_Qend == 4)
  if (nrow(n3) > 0){
    for (l in (1:nrow(n3))){
      if (length(matchPattern(ins(n3[l,]),ins_amplicon(n3[l,], Amplicon))) == 0){ # ref sequences is the seq and not the amplicon, as I have deletion in this case
        n3[l,]$filter_step <- "1 deletion 3 bases max away from the breaksite and at least one mismatch"
      }
      else {
        n3[l,]$filter_step <- "only 1 deletion 3 bases away from the breaksite"
      }
    }
  }
  
  print("n2 and n3 done")
  
  # cases >5 for both 
  # If I have TGATAGTAGG with 1, 2 or 3 mismatches -> uncut
  n4 <- (filter(df, ! Qname %in% n9$Qname & ! Qname %in% n1$Qname & ! Qname %in% n7$Qname & ! Qname %in% n8$Qname & ! Qname %in% n2$Qname & ! Qname %in% n3$Qname & filter_step==""))
  
  if (nrow(n4) >0){
    n4$seq_to_test <- ""
    for (l in (1:nrow(n4))){
      minus <- if_else(n4[l,]$Amplicon_B_Rend < (bait_length - 4), n4[l,]$B_Qend , n4[l,]$B_Qend + (bait_length - n4[l,]$Amplicon_B_Rend - as.integer(4)))
      plus <- if_else(n4[l,]$Amplicon_Rstart > bait_length +5 , n4[l,]$Qstart, n4[l,]$Qstart + (bait_length + as.integer(5) - n4[l,]$Amplicon_Rstart))
      sub <- substr(as.character(n4[l,]$Seq), minus, plus)
      n4[l,]$seq_to_test <- sub
      if (length(matchPattern(seq_breaksite, sub, max.mismatch=3)) != 0){
        n4[l,]$filter_step <- "find 10-base-length breaksite sequence with max 3 mismatches"
      }
    }
  }
  
  # for reads remaining in n5, check the base call quality for the bases in between B_Qend and Qstart
  n6 <- filter(n4, filter_step == "")
  if (nrow(n6) >0){
    n6 <- qual_average(n6, fastq_R1)
  }
  n6$filter_step <- if_else(n6$mean < 25, "bad call quality", "complex junction, quality ok")
  
  # for junctions qualitfied as 'complex junction, qulity ok', if the alignemnt is only in the repeated region (window extended by 2 on both sides), then I remove it
  if (RepSeq & nrow(n6) >0){
    for (l in 1:nrow(n6)){
      if (n6[l,]$Rstart >= (meta$Repeat_Start - 2) & n6[l,]$Rend <= (meta$Repeat_End + 2) & n6[l,]$filter_step == "complex junction, quality ok"){
        n6[l,]$filter_step <- "complex junction, only in repeated region"
      }
    }
  }
  print("all filters done")
  
  # Collect all the results together
  tlx_all <- rbind(
    filter(df, ! Qname %in% n9$Qname & ! Qname %in% n1$Qname & ! Qname %in% n7$Qname & ! Qname %in% n8$Qname & ! Qname %in% n2$Qname & ! Qname %in% n3$Qname & ! Qname %in% n4$Qname) %>% 
      select(Qname, Junction, Rstart, Rend, B_Rstart, B_Rend, B_Qstart, B_Qend, Qstart, Qend, Qlen, Seq, filter_step),
    filter(n9, !Qname %in% n10$Qname) %>% select(Qname, Junction, Rstart, Rend, B_Rstart, B_Rend, B_Qstart, B_Qend, Qstart, Qend, Qlen, Seq, filter_step),
    n10 %>% select(Qname, Junction, Rstart, Rend, B_Rstart, B_Rend, B_Qstart, B_Qend, Qstart, Qend, Qlen, Seq, filter_step),
    n7 %>% select(Qname, Junction, Rstart, Rend, B_Rstart, B_Rend, B_Qstart, B_Qend, Qstart, Qend, Qlen, Seq, filter_step),
    n8 %>% select(Qname, Junction, Rstart, Rend, B_Rstart, B_Rend, B_Qstart, B_Qend, Qstart, Qend, Qlen, Seq, filter_step),
    filter(n1, !Qname %in% n11$Qname) %>% select(Qname, Junction, Rstart, Rend, B_Rstart, B_Rend, B_Qstart, B_Qend, Qstart, Qend, Qlen, Seq, filter_step),
    n11 %>% select(Qname, Junction, Rstart, Rend, B_Rstart, B_Rend, B_Qstart, B_Qend, Qstart, Qend, Qlen, Seq, filter_step),
    n2 %>% select(Qname, Junction, Rstart, Rend, B_Rstart, B_Rend, B_Qstart, B_Qend, Qstart, Qend, Qlen, Seq, filter_step),
    n3 %>% select(Qname, Junction, Rstart, Rend, B_Rstart, B_Rend, B_Qstart, B_Qend, Qstart, Qend, Qlen, Seq, filter_step),
    filter(n4, !Qname %in% n6$Qname) %>% select(Qname, Junction, Rstart, Rend, B_Rstart, B_Rend, B_Qstart, B_Qend, Qstart, Qend, Qlen, Seq, filter_step),
    n6 %>% select(Qname, Junction, Rstart, Rend, B_Rstart, B_Rend, B_Qstart, B_Qend, Qstart, Qend, Qlen, Seq, filter_step)
  )
    
  tlx_all$output <- ""
  tlx_all$output <- sapply(tlx_all$filter_step, output_result)
  
  tlx_all <- tlx_all %>% mutate(JuncID = 1,
                     Rname = meta$Chr,
                     Strand = ifelse(meta$Strand == "+",1,-1),
                     B_Rname = meta$Chr,
                     B_Strand = ifelse(meta$Strand == "+",1,-1)) %>% select(Qname, JuncID, Rname , Junction, Strand, Rstart, Rend, B_Rname, B_Rstart, B_Rend, B_Strand, B_Qstart, B_Qend, Qstart, Qend, Qlen, Seq, filter_step, output)

  tlx_rejoining <- filter(tlx_all, output == "rejoining")
  unmod <- filter(tlx_all, output == "uncut")
  excluded <- filter(tlx_all, output == "bad quality")
  
  print(paste("tlx_rejoining: ", nrow(tlx_rejoining), sep="")) 
  print(paste("unmod: ", nrow(unmod) + n_uncut, sep="")) 
  print(paste("excluded: ", nrow(excluded), sep="")) 
  print(paste("sum: ", nrow(excluded) + nrow(unmod) + nrow(tlx_rejoining) + n_uncut, sep=""))
  print(paste("total from bowtie: ", nrow(df) + n_uncut, sep=""))
  print(paste("CHECK: total from bowtie and sum of my files are the same: ", nrow(df) == nrow(excluded) + nrow(unmod) + nrow(tlx_rejoining), sep=""))
  print("collecting all rejoining events done")
  print(Sys.time())
  
  # write stats file remove unnecesserary files
  write(c(paste("From Bowtie2, number of reads with alignment to the bait and the prey : ", nrow(excluded) + nrow(unmod) + nrow(tlx_rejoining) + n_uncut, " including", sep=""),
               paste("-unmodified reads : ", nrow(unmod) + n_uncut, sep=""),
               paste("-reads excluded for sequencing quality issues : ", nrow(excluded), sep=""),
               paste("-reads found as rejoining junctions : ", nrow(tlx_rejoining), sep="")),
             file = paste(lib, "_stats.txt", sep=""))
  file.remove("bait.1.bt2", "bait.2.bt2", "bait.3.bt2", "bait.4.bt2", "bait.rev.1.bt2", "bait.rev.2.bt2","prey.1.bt2", "prey.2.bt2", "prey.3.bt2", "prey.4.bt2", "prey.rev.1.bt2", "prey.rev.2.bt2", paste(lib, ".fa", sep=""), paste(lib, "_with_prey.fa", sep=""))
  file.remove(paste(lib, "_bait.sam", sep=""), paste(lib, "_prey.sam", sep=""))
  if (!keep){
    file.remove(paste(lib, "_unmod.tlx", sep=""))
  }
  if (keep){
    write.table(tlx_rejoining, file=paste(lib, "_rejoining.tlx", sep =""), sep="\t", quote = FALSE, row.names = FALSE)
    write.table(unmod, file=paste(lib, "_unmod.tlx", sep =""), sep="\t", quote = FALSE, row.names = FALSE, col.names=FALSE, append=TRUE)
    write.table(excluded, file=paste(lib, "_excluded.tlx", sep =""), sep="\t", quote = FALSE, row.names = FALSE)
  }
  return(tlx_rejoining)
}

### SUBROUTINES

#function to collect and sum all the values for a specific code in the cigar (M, I or D)
#in the cigar from our bowtie, I do not have any N, H, P, X or =, I have only S, M, D, I
cigar_coll <- function(cigar, code){
  c1 <- str_extract_all(cigar, paste("\\d+", code, sep=""))
  c2 <- unlist(c1, use.names=FALSE)
  c3 <- as.numeric(str_extract(c2, "\\d+"))
  return(sum(c3))
}

# one function to get the info from the cigar :the start of the alignment, the length of the alignment on the sequencem and the length of the alignment on the reference genome
cigar_info <- function(cigar) {
  s1 <- str_extract(cigar,"^\\d+S")
  s2 <- as.numeric(str_extract(s1, "^\\d+"))
  s2 <- if_else(is.na(s2), 0, s2) 
  m1 <- cigar_coll(cigar,"M")
  i1 <- cigar_coll(cigar,"I")
  d1 <- cigar_coll(cigar, "D")
  return(c(s2+1, s2+m1+i1, m1+d1))
}

#function to get the collection of the insertions between bait and prey
ins <- function (df) {
  N <- nrow(df)
  coll <- character(N)
  for (i in (1:N)){
    coll[i] <- substr(as.character(df$Seq[i]), df$B_Qend[i]+1, df$Qstart[i]-1)
  }
  return(coll)
}
ins_amplicon <- function (df, Amplicon) {
  N <- nrow(df)
  coll <- character(N)
  for (i in (1:N)){
    coll[i] <- substr(as.character(Amplicon), df$Amplicon_B_Rend[i]+1, df$Amplicon_Rstart[i]-1)
  }
  return(coll)
}

# function to calculate the % of mismatches between 2 sequences
pp_mismatch <- function(B_Qend, Qstart, Amplicon_B_Rend, Amplicon_Rstart, Amplicon, Seq){
  ins_ampl <- substr(as.character(Amplicon), Amplicon_B_Rend+1, Amplicon_Rstart-1)
  nm <- nmismatch(ins_ampl, Views(as.character(Seq), start = B_Qend + 1, end = Qstart - 1 )) / (Qstart - B_Qend - 1)
  return(nm)
}

# function to get the average of the base call quality between B_Qend and Qstart
qual_average <- function(n, fq){
  n$mean <- NA
  for (l in (1:nrow(n))){
    i <- grep(paste(n[l,]$Qname, " ", sep=""), ShortRead::id(fq), ignore.case = TRUE)
    qual <- as(quality(fq)[i], "matrix")
    length_min <- min((n[l,]$Qstart -1), length(qual))
    if ((n[l,]$B_Qend + 1) > length_min){ #cases where the bait end is even after the fq range, I cannot have any info about the quality, I do not keep the read
      n[l,]$mean <- 0
    }
    else { 
      n[l,]$mean <- mean(qual[1,(n[l,]$B_Qend + 1):length_min]) 
    }
  }
  return(n)
}

# function to decide the output from the filter
output_result <- function(filter){
  if (filter %in% c("insertion only", "resection only", "deletion leading to MH in the sequence bait/prey", "2 mismatches at the breaksite with good quality",
                    "resection (1 or more) 1 base away from the break or mismatch + resection","insertion (1 or more) 1 base away from the break or mismatch + insertion",
                    "more (strictly) than 1 insertion 2 bases or less away from the breaksite", "more (strictly) than 1 deletion 2 bases or less away from the breaksite",
                    "more (strictly) than 1 insertion 3 bases or less away from the breaksite", "more (strictly) than 1 deletion 3 bases or less away from the breaksite",
                    "many mismatches at the breaksite with good quality",  "1 insertion 3 bases max away from the breaksite and at least one mismatch",
                    "1 deletion 3 bases max away from the breaksite and at least one mismatch", "1 deletion 2 bases max away from the breaksite and at least one mismatch",
                    "1 insertion 2 bases max away from the breaksite and at least one mismatch",
                    "complex junction, quality ok")){
    return("rejoining")
  }
  if (filter %in% c("perfectly uncut","1 mismatch at the breaksite", "1 mismatch out of 2 at the breaksite", "2 mismatches at the breaksite with poor quality",
                    "few mismatches at the breaksite", "only 1 insertion 2 bases away from the breaksite",
                    "only 1 deletion 2 bases away from the breaksite", "only 1 insertion 3 bases away from the breaksite",
                    "only 1 deletion 3 bases away from the breaksite", "find 10-base-length breaksite sequence with max 3 mismatches", "find 10-base-length breaksite sequence with 4 mismatches")){
    return("uncut")
  }
  if (filter %in% c("bad call quality","many mismatches at the breaksite with poor quality","complex junction, only in repeated region")){
    return("bad quality")
  }
}

# function to combine the rejoining events and the results from Wrapper: 
# for Qnames found both by wrapper, and as rejoining events here, the choice is based on the order of the sequential junctions
# If the alignment from Wrapper are in the prey region, then take the one from JoinT, if it is in a different location, then take the first one on the sequence
# If there is a repeated region, then I add one step : if the alignment from JoinT is strictly in the repeated region only : if the alignments from wrapper is longer (longer than the length of the repeated region +8, adjusted for Myc CACA repeat), then I keep the alignment from wrapper, else I keep the rejoining event.
combine <- function(tlx_rejoining, result, meta, RepSeq, lib, keep){
  print(paste("Combining results from wrapper with rejoining events, RepSeq is ", RepSeq, sep=""))
  result <- result %>% mutate(filter_step = NA, output = NA)
  tlx_rejoining <- tlx_rejoining %>% mutate(B_Cigar = NA,
                                                      Cigar = NA,
                                                      J_Seq = NA,
                                                      Barcode = NA,
                                                      unaligned = NA,
                                                      baitonly = NA,
                                                      uncut = NA,
                                                      misprimed = NA,
                                                      freqcut = NA,
                                                      largegap = NA,
                                                      mapqual = NA,
                                                      breaksite = NA,
                                                      sequential = NA,
                                                      repeatseq = NA,
                                                      duplicate = NA)
  common_RB <- filter(tlx_rejoining, Qname %in% result$Qname) %>% arrange(Qname)
  common_W <- filter(result, Qname %in% tlx_rejoining$Qname) %>% arrange(Qname)
  to_keep <- data.frame()
  # cases where I do not have any common Qnames
  if (nrow(common_RB) == 0){
    result_all <- rbind(tlx_rejoining, result) %>% select(Qname, JuncID, Rname , Junction, Strand, Rstart, Rend, B_Rname, B_Rstart, B_Rend, B_Strand, B_Qstart, B_Qend, Qstart, Qend, Qlen, 
                                                   B_Cigar, Cigar, Seq, J_Seq, Barcode, unaligned, baitonly, uncut, misprimed, freqcut, largegap, mapqual, breaksite, sequential, repeatseq, duplicate, filter_step)
    print(paste("CHECK: all Qnames are unique in the result_all file: ", nrow(result_all) == length(unique(result_all$Qname)), sep=""))
    print(paste("CHECK: I have included all Qnames from result and JoinT in result_all: ", nrow(result_all) == nrow(result) + nrow(tlx_rejoining) - nrow(filter(result, Qname %in% tlx_rejoining$Qname)), sep=""))
    if (keep){
      write.table(result_all, file=paste(lib, "_result_all.tlx", sep =""), sep="\t", quote = FALSE, row.names = FALSE)
    }
    
    # get the final output file with only junctions in the main chromosomes
    JoinToutput <- filter(result_all, Rname %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", 
                                                   "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"))
    rejoin <- filter(JoinToutput, as.character(Rname) == as.character(meta$Chr) & Junction >= start_of_prey & Junction <= end_of_prey)
    write.table(JoinToutput, file=paste(lib, "_JoinT.tlx", sep =""), sep="\t", quote = FALSE, row.names = FALSE)
    
    # add infos to stats file
    write(c("Combining the rejoining events to the results from wrapper : ",
            paste("-total number of junctions : ", nrow(result_all), sep=""),
            paste("-total number of junctions in main chromosomes : ", nrow(JoinToutput), sep=""),
            paste("-number of rejoining junctions, defined as junction with the start of prey within the amplicon : ", nrow(rejoin), sep=""),
            paste("-number of translocations : ", nrow(JoinToutput) - nrow(rejoin), sep="")),
          file = paste(lib, "_stats.txt", sep=""), append = TRUE)
    print("combining with results done")
    return() 
  }
  # treat the cases of Qnames common in Wrapper and JoinT
  Amplicon <- DNAString(meta$Amplicon)
  lA <- length(Amplicon)
  start_of_prey <- ifelse(meta$Strand == "+", meta$End, meta$End - lA)
  end_of_prey <- ifelse(meta$Strand == "+", meta$Start + lA, meta$Start)
  start_of_repeat <- ifelse(RepSeq, meta$Repeat_Start, 0)
  end_of_repeat <- ifelse(RepSeq, meta$Repeat_End, 0)
  length_repeat <- end_of_repeat - start_of_repeat
  for (i in (1:nrow(common_W))){
    # if the result from wrapper is in the prey region, I keep the alignment from RB
    if (common_W$Rname[i] == meta$Chr & common_W$Junction[i] > start_of_prey & common_W$Junction[i] < end_of_prey){
      to_keep <- rbind(to_keep, common_RB[i,])
      next
    }
    # for cases where I have a repeated regions, if the alignment from JoinT falls only in the repeated region, then I check the length of the translocation found by wrapper to choose which one I keep
    if (RepSeq & common_RB$Rstart[i] > (start_of_repeat - 10) & common_RB$Rend[i] < end_of_repeat + 10){
      if (common_W$Qend[i] - common_W$Qstart[i] > (length_repeat + 8) ){
        to_keep <- rbind(to_keep, common_W[i,])
      }
      else {
        to_keep <- rbind(to_keep, common_RB[i,])
      }
      next
    }
    # for all other cases: I keep the first alignment on the sequence
    else {
      if (common_RB$Qstart[i] <= common_W$Qstart[i]){
        to_keep <- rbind(to_keep, common_RB[i,])
      }
      else {
        to_keep <- rbind(to_keep, common_W[i,])
      }
    }
  }
  result_all <- rbind(filter(tlx_rejoining, !Qname %in% common_RB$Qname), filter(result, !Qname %in% common_W$Qname), to_keep) %>% 
    select(Qname, JuncID, Rname , Junction, Strand, Rstart, Rend, B_Rname, B_Rstart, B_Rend, B_Strand, B_Qstart, B_Qend, Qstart, Qend, Qlen, 
           B_Cigar, Cigar, Seq, J_Seq, Barcode, unaligned, baitonly, uncut, misprimed, freqcut, largegap, mapqual, breaksite, sequential, repeatseq, duplicate, filter_step)
  print(paste("CHECK: all Qnames are unique in the result_all file: ", nrow(result_all) == length(unique(result_all$Qname)), sep=""))
  print(paste("CHECK: I have included all Qnames from result and JoinT in result_all: ", nrow(result_all) == nrow(result) + nrow(tlx_rejoining) - nrow(filter(result, Qname %in% tlx_rejoining$Qname)), sep=""))
  if (keep){
    write.table(result_all, file=paste(lib, "_result_all.tlx", sep =""), sep="\t", quote = FALSE, row.names = FALSE)
  }
  
  # get the final output file with only junctions in the main chromosomes
  JoinToutput <- filter(result_all, Rname %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", 
                                                 "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"))
  rejoin <- filter(JoinToutput, as.character(Rname) == as.character(meta$Chr) & Junction >= start_of_prey & Junction <= end_of_prey)
  write.table(JoinToutput, file=paste(lib, "_JoinT.tlx", sep =""), sep="\t", quote = FALSE, row.names = FALSE)
  
  # add infos to stats file
  write(c("Combining the rejoining events to the results from wrapper : ",
               paste("-total number of junctions : ", nrow(result_all), sep=""),
               paste("-total number of junctions in main chromosomes : ", nrow(JoinToutput), sep=""),
               paste("-number of rejoining junctions, defined as junction with the start of prey within the amplicon : ", nrow(rejoin), sep=""),
               paste("-number of translocations : ", nrow(JoinToutput) - nrow(rejoin), sep="")),
             file = paste(lib, "_stats.txt", sep=""), append = TRUE)
  print("combining with results done")
  return()
}
