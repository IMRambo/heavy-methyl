#Ian Rambo
#Created: February 11, 2016
#Last updated: February 29, 2016

#PURPOSE: analyze and plot DNA curvature data stored in HDF5 format.
#==============================================================================
# -----------               ------------
#  --------   \___________/   ---------
#    -------    O       O   ---------
#      -------\     V     /--------
#              \  _____  /        
#               \/     \/ 
#==============================================================================
require("dplyr")
require("ggplot2")
#require("multidplyr")
require("tidyr")
require("grid")
require("gridExtra")
require("rhdf5")
require("data.table")
#==============================================================================
setwd("/Users/imrambo/R/ORMP/data/")

#Coefficient of variation
covar <- function(x){
  (sd(x)/mean(x))
}
#==============================================================================
#-----MAPPING FILES-----
#Taxonomic lineages
lineage <- read.table("./allBacteriaTaxonomy_uniq.txt",
                      sep = ";", header = TRUE, stringsAsFactors = FALSE) %>%
  filter(Class != "NA")
#Table of accession numbers and respective taxa
bacteriaSummary <- read.table("./bacterial_genomes_summary.txt", sep = "\t",
                              header = TRUE)
#Select only the accession number, not the version
bacteriaSummary$Accession <- substr(bacteriaSummary$Accession, 1, 9)
#Extract the genera - create a list of species name components
#and select the first element
taxnames <- strsplit(as.character(bacteriaSummary$TaxName), split = " ")
Genus <- sapply(taxnames, "[[", 1)

bacteriaSummary <- cbind(bacteriaSummary, Genus) %>%
  filter(grepl("chromosome*", Replicon)) %>%
  rename(accession = Accession) %>%
  select(accession, Genus) %>%
  mutate(Genus = as.character(Genus)) %>%
  left_join(lineage)

#Taxonomic lineage paths for creating HDF5 groups in Python
# taxonomyPaths <- bacteriaSummary[ , c(6:1)] %>% unite(fullPath, Phylum:accession, sep = "/")
# write.table(taxonomyPaths, file = "genomePaths.txt", row.names = FALSE,
#             col.names = FALSE, quote = FALSE)
#==============================================================================
#LOAD IN AND WRANGLE DATA FILES

#-----HMMER RESULTS DATA LOAD-----
#Pfam Clan mapping file
pfamClanMap <- read.table("PfamClanMapping.txt", header = TRUE) %>%
  mutate(PfamClan = as.character(PfamClan), PfamFamily = as.character(PfamFamily))

pfamHMMFiles <- list.files("02-targetGenomes_PfamHMM_bestHit")

#Combine pfam best hit results files into a data frame
pfamData <- rbindlist(lapply(pfamHMMFiles, function(x){
  x <- paste("02-targetGenomes_PfamHMM_bestHit", x, sep = "/")
  fread(x, header = FALSE, skip = 3, col.names = c("accession","cds","PfamFamily"))}))

pfamData <- pfamData %>%
  as.data.frame() %>%
  distinct() %>% 
  left_join(pfamClanMap) %>% 
  na.omit()

#   
#   #CDS Pfam clan assignments for data selection in HDF5
#   cdsClans <- pfamData %>% as.data.frame() %>% filter(header = grepl("cds[0-9]+;", header)) %>%
#   separate(header, sep = ";", into = c("cds","accession","start","stop","sense","protein","proteinName")) %>%
#   select(accession, cds, PfamFamily, evalue) %>%
#   mutate(PfamFamily = as.character(PfamFamily),
#          evalue = as.character(evalue),
#          accession = substr(accession, 1, 9)) %>%
#   left_join(pfamBestHit) %>%
#   distinct() %>% 
#   left_join(pfamClanMap) %>% 
#   na.omit() %>%
#   mutate(evalue = as.numeric(evalue)) %>%
#   #select(accession, cds, PfamFamily, PfamClan) %>%
#   as.data.table() 
#------------------------------------------------------------------------------
#-----CURVATURE RESULTS HDF5 LOAD----- 
#------------------------------------------------------------------------------
curveHDF5 <- "/Users/imrambo/R/ORMP/data/02-CurveDataHDF5/test.hdf5"
curveStructure <- h5ls(curveHDF5)
#Paths in the HDF5 structure that contain data frames of interest
curvePaths <- curveStructure %>% filter(grepl("DATASET", otype) & dclass == "COMPOUND") %>%
  select(group, name) %>% unite(path, c(group, name), sep = "/", remove = TRUE)
#Load the data sets into a data frame
curveData <- rbindlist(lapply(curvePaths$path, function(x){ h5read(curveHDF5, x) %>%
    mutate(accession = strsplit(x,"/")[[1]][7],
             cds = strsplit(x,"/")[[1]][8],
             ntRelPosition = seq(-200, 299))})) 
#Close the HDF5 file
H5close()  
#Summary of curvature results
curveSummary <- curveData %>% group_by(accession, ntRelPosition) %>%
  summarise(meanCurve = mean(curvature))

curveLine <- ggplot(curveSummary, aes(x = ntRelPosition, y = meanCurve, color = accession)) + geom_point()


