#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
list.of.packages <- c( "dplyr", "stringr","data.table","ggplot2", "foreach", "doParallel" )
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')
library(ggplot2)
library(dplyr)
library(data.table)
library(stringr)
OUTPUT_FOLDER=args[2]
#readsdf = data.frame(readlabel = character(), chrom =character(), startcoord= character(), status = character())

con = file(paste0(OUTPUT_FOLDER,"/modkit_referenced_splitbam_mergedbychr/modkit_extract_",args[1],"_merged_filtered_coordmod.tsv"), "r", )
# NOTE: stale file.remove of HG002_read_format/HG002_*_reads.tsv disabled — path is not
# used by this pipeline and triggers a warning. Output target is ${OUTPUT_FOLDER}/read_format/.
#file.remove(paste0("HG002_read_format/HG002_", args[1],"_reads.tsv"))
write.table( data.frame(readlabel = "readlabel", 
                        haplotype = "haplotype",
                        chrom ="chrom", 
                        #string = "string",  no string in current format
                        startcoord= "startcoord",
                        status = "status"),
             paste0(OUTPUT_FOLDER,"/read_format/", args[1],"_reads.tsv"), sep = "\t", row.names = F, quote=F, append = T, col.names  = F)
#header
readCGcoordscollection = c()
readMethylation_statuscollection = c()
line = readLines(con, n = 1)
if ( length(line) == 0 ) {
  quit(save = "no", status = 0)
}
line = readLines(con, n = 1)
if ( length(line) == 0 ) {
  quit(save = "no", status = 0)
}
case = c(str_split(line, "\t", simplify = T))
#string = case[2]  no string in current format
readlabel = case[1]
CGcoords = case[2]
chrom = case[3]
Methylation_status = 1*(case[4] > 0.5)
haplotype = case[5]
set = cbind(CGcoords, Methylation_status)

readCGcoordscollection = c(readCGcoordscollection, set[,1])
readMethylation_statuscollection = c(readMethylation_statuscollection,  rep(Methylation_status, length(set[,1])))

while(TRUE){
  line = readLines(con, n = 1)
  if ( length(line) == 0 ) {
    break
  }
  case = c(str_split(line, "\t", simplify = T))
  if(case[1] != readlabel){
    write.table( data.frame(readlabel = readlabel, 
                            haplotype = haplotype,
                            chrom =chrom, 
                            #string = string, no string in current format
                            startcoord= paste0(readCGcoordscollection, collapse=","),
                            status = paste0(readMethylation_statuscollection, collapse=","), 
                            stringsAsFactors= F),
                 paste0(OUTPUT_FOLDER,"/read_format/", args[1],"_reads.tsv"), sep = "\t", row.names = F, quote=F, append = T, col.names  = F)
    
    readCGcoordscollection = c()
    readMethylation_statuscollection = c()
  }
  #string = case[2]
  readlabel = case[1]
  CGcoords = case[2]
  chrom = case[3]
  Methylation_status = 1*(case[4] > 0.5)
  haplotype = case[5]
  set = cbind(CGcoords, Methylation_status)
  readCGcoordscollection = c(readCGcoordscollection, set[,1])
  readMethylation_statuscollection = c(readMethylation_statuscollection, rep(Methylation_status, length(set[,1])))
  
  
}

write.table( data.frame(readlabel = readlabel, 
                        haplotype = haplotype,
                        chrom =chrom, 
                        #string = string,
                        startcoord= paste0(readCGcoordscollection, collapse=","),
                        status = paste0(readMethylation_statuscollection, collapse=","), 
                        stringsAsFactors= F),
             paste0(OUTPUT_FOLDER,"/read_format/", args[1],"_reads.tsv"), sep = "\t", row.names = F, quote=F, append = T, col.names  = F)
