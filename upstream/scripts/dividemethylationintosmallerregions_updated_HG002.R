#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

list.of.packages <- c( "dplyr", "stringr","data.table","ggplot2", "foreach", "doParallel" )
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')
library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)
seqlast <- function (from, to, by) {
  vec <- do.call(what = seq, args = list(from, to, by))
  if ( tail(vec, 1) != to ) {
    return(c(vec, to))
  } else {
    return(vec)
  }
}
 
OUTPUT_FOLDER=args[2]
#fullset = data.frame(read.table(paste0("HG002_",args[1],"_reads.tsv"), header = T, stringsAsFactors = F))

#fullset = data.frame(read.table("newnewmethylation_calls_readstatusonlychr20complete_week20.tsv", header = T, stringsAsFactors = F))
fullset = data.frame(read.table(paste0(OUTPUT_FOLDER, "/read_format/",args[1],"_reads.tsv"), header = T, stringsAsFactors = F))
#fullset = distinct(fullset)


CpGcount = 10000
windowsize = 30
# Optional 3rd arg: shared CpG coords file (one coord per line). When provided, we use this
# as the canonical CpG coordinate set so region boundaries match across samples/time points.
# Without it, we fall back to per-sample coords (original BIC ASM behavior).
if (length(args) >= 3 && nchar(args[3]) > 0 && file.exists(args[3])) {
  allcoords = as.numeric(readLines(args[3]))
} else {
  allcoords = as.numeric(unique(unlist(lapply(fullset$startcoord, str_split, ","))))
}
allcoords = allcoords[order(allcoords)]
startend_reads = matrix(NA, ncol = 2, nrow = dim(fullset)[1])
for(read in 1:dim(fullset)[1]){
  case = unlist(str_split(fullset[read,"startcoord"],","))
  startend_reads[read,] =   as.numeric(case[c(1, length(case))])
}
orderjustincase = order(startend_reads[,1])
startend_reads = startend_reads[orderjustincase ,]
fullset = fullset[orderjustincase ,]
endposition = seqlast(1, length(allcoords), by = CpGcount)
regionends = allcoords[endposition[-1]]
regionstarts = allcoords[c(1, endposition[-c(1,length(endposition))] - (windowsize-1))]


dir.create(paste0(OUTPUT_FOLDER, "/read_format_split/", args[1], "/"), recursive = TRUE)

for(cont in 1:length(regionstarts)){
  fraction = fullset[startend_reads[,1] <= regionends[cont] & startend_reads[,2] >= regionstarts[cont],]
  write.table(fraction, paste0(OUTPUT_FOLDER, "/read_format_split/", args[1], "/methylationfraction_",regionstarts[cont],"_",regionends[cont],"_.tsv" ), quote = F, row.names = F, sep = "\t")
}
#write.table(data.frame(regionstarts[-(length(regionstarts))],regionends[-1]),"coordsofeachfractionfordgefiltering.txt", quote = F,col.names = F, row.names = F, sep = "\t")
