#Ian Rambo
#Oyster Rocks Methylation Project
#Last updated: October 15, 2015
#BLAST analysis
#==============================================================================
#GOAL: Analyze BLAST-X output from tblastn alignment of bacterial DNA
#methylase queries against target genomes. Target genomes are complete
#genomes from 91 genera comprising the 10 most abundant classes from
#PhymmBL results with class confidence >= 0.700 and order confidence >= 0.650.
#==============================================================================
# -----------               ------------
#  --------   \___________/   ---------
#    -------    O       O   ---------
#      -------\     V     /--------
#              \  _____  /        
#               \/     \/ 
#==============================================================================
library(dplyr)
library(fdrtool)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(car)
library(moments)
#==============================================================================
setwd("/Users/imrambo/Documents/BINF868/")
blast <- read.table(file = "BLAST/targetGenomes_mtase_tblastnResultsFMT.txt",
                    header = FALSE, sep = "\t")

giMethyl <- read.table("Methyltransferases/giMethylase.txt", header = FALSE,
                       sep = "\t")
giGenome <- read.table("TargetGenomes/giGenome.txt", header = FALSE,
                       sep = "\t")

colnames(blast) <- c("qgi", "sgi", "pident", "length", "mismatch",
                   "gapopen", "qstart", "qend", "sstart", "send", "evalue",
                   "bitscore")

colnames(giMethyl) <- c("qgi", "qdatabase", "qname")
colnames(giGenome) <- c("sgi", "sdatabase", "sname")

blast <- mutate(blast, qgi = as.character(qgi),
              sgi = as.character(sgi))
giMethyl <- mutate(giMethyl, qgi = as.character(qgi),
                   qdatabase = as.character(qdatabase),
                   qname = as.character(qname))
giGenome <- mutate(giGenome, sgi = as.character(sgi),
                   sdatabase = as.character(sdatabase),
                   sname = as.character(sname))

bsel <- blast %>% left_join(giGenome) %>%
  left_join(giMethyl) %>%
  mutate(slen = abs(sstart - send)) %>%
  select(qgi,qdatabase,qname,sgi,sdatabase,sname,slen,evalue,bitscore) 

hist(bsel$evalue)

blastFDR <- fdrtool(bsel$evalue, statistic = "pvalue", verbose = TRUE)
blastFDR$qval
