#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(argparser))
suppressPackageStartupMessages(library(magrittr))

parser <- arg_parser("mark duplicate juntions",
                     name="TranslocDedup.R") %>%
    add_argument("tlxfile",
                 "",
                 type="character") %>%
    add_argument("output",
                 "",
                 type="character") %>%
    add_argument("--offset.dist",
                 "",
                 default=0,
                 type="integer") %>%
    add_argument("--break.dist",
                 "",
                 default=0,
                 type="integer") %>%
    add_argument("--barcode",
                 "",
                 default=NA,
                 type="character") %>%
    add_argument("--cores",
                 "",
                 default=1,
                 type="integer")

if (interactive()) {
    argv <- parse_args(parser,argv=c("~/AltLab/Data/Alt055/preresults-sample/RF204_Alt055/RF204_Alt055.tlx",
                                     "~/AltLab/Data/Alt055/preresults-sample/RF204_Alt055/RF204_Alt055_duplicate1.txt",
                                     "--cores", 4))
} else {
    argv <- parse_args(parser)
}

suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(dplyr))

l_ply(names(argv), function(i) {
  if (i != "help" && i != "opts" && i != "") {
      assign(i, argv[[i]], envir=.GlobalEnv)
  }
})

if (break.dist > 0 || offset.dist > 0) {
    stop("Error: still working on implementation for non-zero break.dist and offset.dist")
}

if (cores > 1) {
    suppressPackageStartupMessages(library(doMC))
    registerDoMC(cores=cores)
    do_parallel <- TRUE
} else {
    do_parallel <- FALSE
}

con  <- file(tlxfile, open = "r")
tlxheader <- unlist(strsplit(readLines(con, n = 1),"\t"))
close(con)

col_collectors <- rep(list(col_skip()), length(tlxheader)) %>%
    set_names(tlxheader)

col_collectors[["Qname"]] <-
    col_collectors[["Rname"]] <-
    col_collectors[["Barcode"]] <-
    col_collectors[["Seq"]] <- col_character()
col_collectors[["Junction"]] <-
    col_collectors[["Strand"]] <-
    col_collectors[["B_Qend"]] <-
    col_collectors[["Qstart"]] <-
    col_collectors[["Qend"]] <-
    col_collectors[["B_Strand"]] <-
    col_collectors[["B_Rstart"]] <-
    col_collectors[["B_Rend"]]  <-
    col_collectors[["mapqual"]] <- col_integer()

tlxs <- read_tsv(tlxfile, col_types=col_collectors, progress=F)

if (nrow(tlxs) == 0) {
    cat("Qname\tDups",file=output)
    quit()
}

tlx.juncs <- filter(tlxs, !is.na(Junction) & Rname != "Adapter") %>%
    arrange(Qname) %>%
    mutate(Offset = ifelse(Strand==1, Junction - Qstart, Junction + Qstart),
           B_Junction = ifelse(B_Strand==1, B_Rend, B_Rstart),
           MH_insertion = Qstart - B_Qend-1,
           insertion = ifelse(MH_insertion > 0, substr(Seq,B_Qend+1,Qstart-1) ,NA))

dup.rows <- tlx.juncs %>% group_by(Rname, Strand, B_Junction,  Junction, MH_insertion, insertion)

out.df <- group_split(dup.rows) %>%
    ldply(function(d) {
        if (nrow(d) > 1) {
            return(as.data.frame(cbind(d$Qname[2:nrow(d)] ,d$Qname[1])))
        } else {
            return(NULL)
        }
    }, .id=NULL, .parallel=do_parallel)

write.table(out.df, output, sep="\t", quote=F, row.names=F, col.names=F)
