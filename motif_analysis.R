#Ian Rambo
#Created: January 5, 2016
#Last updated: January 30, 2016
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
require("multidplyr")
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

pfamData <- rbindlist(lapply(pfamHMMFiles, function(x){
  x <- paste("02-targetGenomes_PfamHMM_bestHit", x, sep = "/")
  fread(x, header = FALSE, skip = 3, col.names = c("accession","cds","PfamFamily"))}))

pfamData %>% filter(grepl("PF01498", PfamFamily)) %>% dim()

pfamData %>% as.data.frame() %>% left_join(pfamClanMap) %>% na.omit() %>% left_join(bacteriaSummary) %>% group_by(Class) %>%
  count(Class, PfamClan) %>%
  filter(n > median(n))  %>% print.data.frame() 
# pfamData <- fread("/Users/imrambo/R/ORMP/data/targetGenomes_PfamHMM/targetGenomes_PfamHMM_select/NC_013416_pfamA_select.txt", header = FALSE, skip = 3,
#                       col.names = c("header","PfamFamily", "evalue")) 
#Best hit per CDS
pfamBestHit <- pfamData %>% filter(header = grepl("cds[0-9]+;", header)) %>%
  separate(header, sep = ";", into = c("cds","accession","start","stop","sense","protein","proteinName")) %>%
  select(accession, cds, PfamFamily, evalue) %>%
    mutate(PfamFamily = as.character(PfamFamily),
         accession = substr(accession, 1, 9)) %>%
  na.omit() %>% 
  distinct() %>% #Here is where I left off; let's figure this out. Keep the best evalue, 
  filter(accession == "NC_013416", cds == "cds1163")
  group_by(accession, cds) %>%
  summarise(evalue = min(evalue)) %>%
  as.data.frame() %>%
  mutate(evalue = as.character(evalue))

pfamBestHit %>% filter(cds == "cds1003", accession == "NC_013416") %>%

#CDS Pfam clan assignments for data selection in HDF5
cdsClans <- pfamData %>% as.data.frame() %>% filter(header = grepl("cds[0-9]+;", header)) %>%
  separate(header, sep = ";", into = c("cds","accession","start","stop","sense","protein","proteinName")) %>%
  select(accession, cds, PfamFamily, evalue) %>%
  mutate(PfamFamily = as.character(PfamFamily),
         evalue = as.character(evalue),
         accession = substr(accession, 1, 9)) %>%
  left_join(pfamBestHit) %>%
  distinct() %>% 
  left_join(pfamClanMap) %>% 
  na.omit() %>%
  mutate(evalue = as.numeric(evalue)) %>%
  #select(accession, cds, PfamFamily, PfamClan) %>%
  as.data.table() 

cdsClans %>% filter(accession == "NC_013416" & cds == "cds1003") %>%
  group_by(accession, cds) %>%
  summarise(evalue = min(evalue), PfamClan = PfamClan)
#-----CURVATURE DATA LOAD-----
# curveFiles <- list.files("./01-targetGenomeCurveData", recursive = TRUE)
# curveData <- rbindlist(lapply(curveFiles, function(file){
#   dt <- fread(file)
#   }))

#--------MOTIF PIPELINE DATA LOAD--------------------

#Create a cluster. Note that there is no support for data.table with multidplyr
#cluster <- create_cluster(3)
#Set default cluster size
#set_default_cluster(cluster)


#Load HDF5 file
motifHDF5 <- "/Users/imrambo/R/ORMP/data/test.hdf5"
motifStructure <- h5ls(motifHDF5)
newStruct <- motifStructure %>% filter(grepl("cds", name))

#Paths for CDS datasets
cdsPaths <- paste(newStruct$group, newStruct$name, sep = "/")
#Combine the datasets
motifData <- rbindlist(lapply(cdsPaths, function(x){h5read(motifHDF5, x) %>%
    mutate(accession = strsplit(x,"/")[[1]][7],
           cds = strsplit(x,"/")[[1]][8],
           ntRelPosition = seq(-175, 324))})) 

motifData %>% left_join(cdsClans) %>% head()






# compClustFiles <- list.files("02-CompiledClusterData/")
# CpGnorm <- as.data.frame(fread(paste(getwd(),"/02-CompiledClusterData/CpGnorm_compiled.txt", sep = ""), header = TRUE)) %>%
#   left_join(bacteriaSummary) %>%
#   partition(accession) 
# 
# relpos <- rep.int(seq(-175, 324), nrow(CpGnorm)/500)
# as.data.frame(CpGnorm) %>%
#   mutate(ntRelPosition = relpos) %>%
#   na.omit() %>%
#   group_by(Class, ntRelPosition) %>%
#   summarise(meanScore = mean(score), cvScore = covar(score)) %>%
#   ggplot(aes(x = ntRelPosition, y = meanScore, color = Class)) +
#   geom_line() +
#   ylab("Mean frequency") +
#   xlab("Nucleotide position") +
#   ggtitle("CpG normalized") +
#   xlim(-180, 325) +
#   scale_color_manual(values = cbPalette) +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(color = "black"),
#         axis.title=element_text(size=9),
#         axis.text.x=element_text(size=9),
#         axis.text.y=element_text(size=9,face="bold",colour="black"),
#         legend.text=element_text(size=6.8),
#         plot.title=element_text(size=10),
#         legend.key.size=unit(0.5, "cm"),
#         legend.title=element_text(size=7))
# 
# #CpG normalized, Coefficient of Variation
# as.data.frame(CpGnorm) %>%
#   mutate(ntRelPosition = relpos) %>%
#   na.omit() %>%
#   group_by(Class, ntRelPosition) %>%
#   summarise(meanScore = mean(score), cvScore = covar(score)) %>%
#   ggplot(aes(x = ntRelPosition, y = cvScore, color = Class)) +
#   geom_line() +
#   ylab("Coefficient of Variation") +
#   xlab("Nucleotide position") +
#   ggtitle("CpG normalized") +
#   xlim(-180, 325) +
#   scale_color_manual(values = cbPalette) +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(color = "black"),
#         axis.title=element_text(size=9),
#         axis.text.x=element_text(size=9),
#         axis.text.y=element_text(size=9,face="bold",colour="black"),
#         legend.text=element_text(size=6.8),
#         plot.title=element_text(size=10),
#         legend.key.size=unit(0.5, "cm"),
#         legend.title=element_text(size=7))
  
#------------------------------------------------------------------------------
#Per-window GATC counts
GATCmotif <- as.data.frame(fread(paste(getwd(),"/02-CompiledClusterData/GATCmotif_compiled.txt", sep = ""))) %>%
  left_join(bacteriaSummary) %>%
  mutate(ntRelPosition = relpos) %>%
  partition(accession)

as.data.frame(GATCmotif) %>% mutate(ntRelPosition = relpos) %>%
  na.omit() %>%
  group_by(Class, ntRelPosition) %>%
  summarise(meanScore = mean(score), cvScore = covar(score)) %>%
  ggplot(aes(x = ntRelPosition, y = cvScore, color = Class)) +
  geom_line() +
  ylab("Coefficient of Variation") +
  xlab("Nucleotide position") +
  ggtitle("GATC motif") +
  xlim(-180, 325) +
  scale_color_manual(values = cbPalette) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title=element_text(size=9),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9,face="bold",colour="black"),
        legend.text=element_text(size=6.8),
        plot.title=element_text(size=10),
        legend.key.size=unit(0.5, "cm"),
        legend.title=element_text(size=7))
#------------------------------------------------------------------------------
GATCmotif %>% na.omit() %>%
  group_by(Class, ntRelPosition) %>%
  summarise(meanScore = mean(score), cvScore = covar(score)) %>%
  ggplot(aes(x = ntRelPosition, y = meanScore, color = Class)) +
  geom_line() +
  ylab("Mean Frequency") +
  xlab("Nucleotide position") +
  ggtitle("GATC motif") +
  xlim(-180, 325) +
  scale_color_manual(values = cbPalette) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title=element_text(size=9),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9,face="bold",colour="black"),
        legend.text=element_text(size=6.8),
        plot.title=element_text(size=10),
        legend.key.size=unit(0.5, "cm"),
        legend.title=element_text(size=7))
#------------------------------------------------------------------------------
GATCmotif %>% right_join(GATCexpected) %>% head()
GATCexpected <- as.data.frame(fread(paste(getwd(),"/02-CompiledClusterData/GATCexpected_compiled.txt", sep = "")))
GATCexpected %>% na.omit() %>%
  left_join(bacteriaSummary) %>%
  group_by(Class, ntPosition) %>%
  summarise(meanScore = mean(score), cvScore = covar(score)) %>%
  ggplot(aes(x = ntPosition, y = meanScore, color = Class)) +
  geom_line() +
  ylab("Mean frequency") +
  xlab("Nucleotide position") +
  ggtitle("GATC expected")

AT %>% filter(ntPosition < 200 & accession == "NC_013416")
# #Directory containing subdirectories with motif data files
# motifDirs <- list.files("01-Test_ClusterData")
# ntData <- data.frame()
# ntData$cds <- vector(mode = "character")
# ntData$accession <- vector(mode = "character")
# ntData$ntPosition <- vector(mode = "numeric")
# 
# #Loop through files and combine data into one data frame
# for (dir in motifDirs[1]){
#   dirPath <-paste(getwd(),"/01-Test_ClusterData/",dir, sep = "")
#   tempDF <- data.frame()
#   tempDF$cds <- vector(mode = "character")
#   tempDF$accession <- vector(mode = "character")
#   tempDF$ntPosition <- vector(mode = "numeric")
#   for (file in list.files(dirPath)[1:2]){
#     filePath <- paste(dirPath,"/",file, sep = "")
#     motifSpread <- read.table(filePath, sep = "\t")
#     scoreMetric <- substr(file, 1, (nchar(file) - 4))
#     colnames(motifSpread) <- c("seqID", c(1:500))
#     motifSpread$seqID <- as.character(motifSpread$seqID)
#     motifTidy <- gather(data = motifSpread, key = seqID, value = "frequency")
#     names(motifTidy)[c(2,3)] <- c("ntPosition",scoreMetric)
#     print(head(motifTidy))
#     print(str(motifTidy))
# #     motifTidy <- motifTidy %>% separate(seqID, sep = ";", into = c("cds","accession","start","stop",
# #                                                                    "strand","protein","description")) %>%
# #       mutate(cds = substr(cds, 2, nchar(cds)), ntPosition = as.numeric(ntPosition))
# #     print(dim(motifTidy))
# #     motifTidy <- motifTidy[ ,c(1,2,8,9)]
# #     rm(motifSpread)
# #     rm(scoreMetric)
# #     rm(motifTidy)
#   }
# }
# 
# 
#   #If the merged data frame does not exist, create it
#   if (!exists("ntData")){
#     ntData <- read.table(file, header = TRUE, sep = " ")
#     print("creating ntData")
#   }
#   #If the merged data frame exists, append to it
#   if (exists("ntData")){
#     print(file)
#     tempData <- read.table(file, header = TRUE, sep = " ")
#     
#     #Get the accession number
#     accession <- gsub("(.*?)_CDS*.*$", "\\1", names(tempData[2]))
#     #Position relative to TSS
#     ntRelPosition <- c(-100:300)
#     ntRelPosition <- ntRelPosition[ which(ntRelPosition != 0)]
#     if (all(is.na(tempData$X))){
#       tempData$X <- NULL
#     }
#     
#     colnames(tempData) <- c("ntPosition", "ATmotif", "CpGmotif", "CpGnorm",
#                             "GATCmotif", "GCcontent", "GCmotif", "GpCnorm",
#                             "observedExpected")
#     tempData$ntRelPosition <- ntRelPosition
#     tempData$accession <- accession
#     ntData <- rbind(ntData, tempData)
#     rm(tempData)
#     rm(accession)
#   }
# }
# 
# #Save the large dataframe for later use in RData format
# save(ntData, file = "bacteriaNTData.RData")

#Load the saved dataframe
load("/Users/imrambo/R/ORMP/data/bacteriaNTData.RData")

#Create a cluster
cluster <- create_cluster(2)
#Set default cluster size
set_default_cluster(cluster)


motifDF <- as.data.frame(fread("testMotifData.txt", sep = "\t", header = TRUE))

save(motifDF, file = "motifData.RData")
#Partition data frame by accession numbers
motifDFPartition <- partition(motifDF, accession)
motifDFPartition %>% 

motifDFPartition %>% group_by(accession, ntPosition) %>%
  filter(metric == "GATCexpected") %>%
  summarise(meanScore = mean(score)) %>%
  as.data.frame() %>%
  ggplot(aes(x = ntPosition, y = meanScore, color = accession)) +
  geom_line() +
  ggtitle("GATC expected")

motifDF %>% group_by(accession, ntPosition) %>%
  filter(metric == "GATCmotif") %>%
  summarise(meanScore = mean(score)) %>%
  ggplot(aes(x = ntPosition, y = meanScore, color = accession)) +
  geom_line() +
  ggtitle("GATC motif")
  
  #facet_wrap(~accession, ncol = 3)

#==============================================================================
#-----MAPPING FILES-----
#Taxonomic lineages
lineage <- read.table("./allBacteriaTaxonomy_uniq.txt",
                      sep = ";", header = TRUE, stringsAsFactors = FALSE)
#Table of accession numbers and respective taxa
bacteriaSummary <- read.table("./bacterial_genomes_summary.txt", sep = "\t",
                              header = TRUE)
#Select only the accession number, not the version
bacteriaSummary$Accession <- substr(bacteriaSummary$Accession, 1, 9)
taxnames <- strsplit(as.character(bacteriaSummary$TaxName), split = " ")
Genus <- sapply(taxnames, "[[", 1)

bacteriaData <- cbind(bacteriaSummary, Genus) %>%
  filter(grepl("chromosome*", Replicon)) %>%
  select(Accession, Genus) %>%
  rename(accession = Accession) %>%
  mutate(Genus = as.character(Genus)) %>%
  right_join(motifDF) %>%
  left_join(lineage) %>%
  filter(Class != "NA") %>%
  rename(Index = ntPosition)



bacteriaData %>% group_by(accession, ntRelPosition) %>%
  filter(Class %in% targetClass) %>%
  summarise(ATmotifavg = mean(ATmotif),
            CpGmotifavg = mean(CpGmotif),
            CpGnormavg = mean(CpGnorm),
            GATCmotifavg = mean(GATCmotif),
            GCcontentavg = mean(GCcontent),
            GpCnormavg = mean(GpCnorm),
            OEavg = mean(observedExpected),
            ATsd = sd(ATmotif, na.rm = TRUE),
#             CpGsd = sd(CpGmotif),
#             GATCsd = sd(GATCmotif),
#             GCsd = sd(GCcontent),
#             GpCsd = sd(GpCnorm),
#             OEsd = sd(observedExpected),
            ATcv = cv(ATmotifavg, ATsd)) %>% head()
sd(bacteriaData$CpGmotif)
# bacteriaBootClass <- bacteriaData %>% group_by(Class) %>%
#   sample_n(10000, replace = TRUE)
# bacteriaBootPhylum <- bacteriaData %>% group_by(Phylum) %>%
#   sample_n(10000, replace = TRUE)


#Mean values for nucleotide targets by phylum
# bacteriaMean <- bacteriaData %>% group_by(Phylum, ntRelPosition) %>%
#   summarise(ATmotifavg = mean(ATmotif),
#             CpGmotifavg = mean(CpGmotif),
#             CpGnormavg = mean(CpGnorm),
#             GATCmotifavg = mean(GATCmotif),
#             GCcontentavg = mean(GCcontent),
#             GpCnormavg = mean(GpCnorm),
#             observedExpectedAvg = mean(observedExpected)) 

#Accession numbers for target classes
targetAcc <- bacteriaData %>% #group_by(Class, ntRelPosition) %>%
  filter(Class %in% targetClass) %>% select(Class, accession) %>%
  distinct()
write.table(targetAcc, "./targetClassAccessions.txt",
                        sep = "\t", quote = FALSE, row.names = FALSE,
            col.names = FALSE)
#==============================================================================
#FIGURES
#==============================================================================
set.seed(1)
#Colorblind-friendly palette
cbPalette <- c("deepskyblue","darkorange2","pink","aquamarine2","darkorchid2",
               "dodgerblue3","chartreuse2","firebrick2","mediumvioletred","brown","black","green")
targetPhyla <- c("Firmicutes", "Proteobacteria", "Euryarchaeota", "Chloroflexi",
                 "Deinococcus-Thermus","Chlamydiae","Spirochaetes")
targetClass <- c("Actinobacteria","Alphaproteobacteria","Bacilli",
                 "Betaproteobacteria","Clostridia","Deinococci",
                 "Deltaproteobacteria","Gammaproteobacteria")
archaeaClass <- c("Archaeoglobi", "Methanobacteria", "Thermococci", "Thermoplasmata")


archaeaClassDF <- bacteriaData %>% group_by(Class, ntRelPosition) %>%
  filter(Class %in% archaeaClass) %>%
  sample_n(10000, replace = TRUE)


atlm <- lm(GATCmotif~ATmotif, data = bacteriaBootClass)
coef(atlm)

ggplot(bacteriaBootClass, aes(x = ATmotif, y = GATCmotif)) +
  geom_point() + geom_density2d() +
  geom_abline(intercept = coef(atlm)[1], slope = coef(atlm)[2]) +
  ylim(0,0.4)
#line plot of normalized CpG per genus
bacteriaData %>% group_by(Genus, ntRelPosition) %>%
  filter(Genus == "Geobacter" | Genus == "Staphylococcus") %>%
  sample_n(10000, replace = TRUE) %>%
  summarise(ATmotifavg = mean(ATmotif),
            CpGmotifavg = mean(CpGmotif),
            CpGnormavg = mean(CpGnorm),
            GATCmotifavg = mean(GATCmotif),
            GCcontentavg = mean(GCcontent),
            GpCnormavg = mean(GpCnorm),
            observedExpectedAvg = mean(observedExpected)) %>%
  ggplot(aes(x = ntRelPosition, color = Genus)) +
  geom_line(aes(y = CpGnormavg), position = "jitter") +
  xlab("Nucleotide position relative to TSS") +
  ylab("Average weighted motif frequency") +
  ggtitle("AT motif") +
  scale_color_manual(values = cbPalette) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title=element_text(size=9),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=9,face="bold",colour="black"),
        legend.text=element_text(size=6.8),
        plot.title=element_text(size=10),
        legend.key.size=unit(0.5, "cm"),
        legend.title=element_text(size=7)) +
  ylim(0, 0.025)
#------------------------------------------------------------------------------
#AT by Phylum
phylumAT <- bacteriaData %>% group_by(Phylum, ntRelPosition) %>%
  filter(Phylum %in% targetPhyla) %>%
  #sample_n(10000, replace = TRUE) %>%
  summarise(ATmotifavg = mean(ATmotif),
            CpGmotifavg = mean(CpGmotif),
            CpGnormavg = mean(CpGnorm),
            GATCmotifavg = mean(GATCmotif),
            GCcontentavg = mean(GCcontent),
            GpCnormavg = mean(GpCnorm),
            observedExpectedAvg = mean(observedExpected)) %>%
  ggplot(aes(x = ntRelPosition, color = Phylum)) +
  geom_line(aes(y = ATmotifavg), position = "jitter") +
  xlab("Nucleotide position relative to TSS") +
  ylab("Average weighted motif frequency") +
  ggtitle("AT motif") +
  scale_color_manual(values = cbPalette) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title=element_text(size=9),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=9,face="bold",colour="black"),
        legend.text=element_text(size=6.8),
        plot.title=element_text(size=10),
        legend.key.size=unit(0.5, "cm"),
        legend.title=element_text(size=7)) 
#------------------------------------------------------------------------------
#AT by Class
classAT <- bacteriaData %>% group_by(Class, ntRelPosition) %>%
  filter(Class %in% targetClass) %>% 
  #sample_n(10000, replace = TRUE) %>%
  summarise(ATmotifavg = mean(ATmotif),
            CpGmotifavg = mean(CpGmotif),
            CpGnormavg = mean(CpGnorm),
            GATCmotifavg = mean(GATCmotif),
            GCcontentavg = mean(GCcontent),
            GpCnormavg = mean(GpCnorm),
            observedExpectedAvg = mean(observedExpected)) %>%
  ggplot(aes(x = ntRelPosition, color = Class)) +
  geom_line(aes(y = ATmotifavg), position = "jitter") +
  xlab("Nucleotide position relative to TSS") +
  ylab("Average weighted motif frequency") +
  ggtitle("AT motif, bootstrap n = 10000") +
  scale_color_manual(values = cbPalette) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title=element_text(size=9),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=9,face="bold",colour="black"),
        legend.text=element_text(size=6.8),
        plot.title=element_text(size=10),
        legend.key.size=unit(0.5, "cm"),
        legend.title=element_text(size=7)) +
  ylim(0.2, 0.8)
#------------------------------------------------------------------------------
#GATC motif for target Phyla
phylumGATC <- bacteriaData %>% group_by(Phylum, ntRelPosition) %>%
  #sample_n(10000, replace = TRUE) %>%
  filter(Phylum %in% targetPhyla) %>%
  summarise(ATmotifavg = mean(ATmotif),
            CpGmotifavg = mean(CpGmotif),
            CpGnormavg = mean(CpGnorm),
            GATCmotifavg = mean(GATCmotif),
            GCcontentavg = mean(GCcontent),
            GpCnormavg = mean(GpCnorm),
            observedExpectedAvg = mean(observedExpected)) %>%
  ggplot(aes(x = ntRelPosition, color = Phylum)) +
  #geom_line(aes(y = CpGmotifavg), position = "jitter", linetype = "dashed") +
  geom_line(aes(y = GATCmotifavg), position = "jitter", linetype = "solid") +
  #geom_line(aes(y = ATmotifavg), position = "jitter") +
  xlab("Nucleotide position relative to TSS") +
  ylab("Average weighted motif frequency") +
  ggtitle("GATC motif, bootstrap n = 10000") +
  #facet_wrap(~Class, ncol = 3, scales = "free_x") +
  scale_color_manual(values = cbPalette) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title=element_text(size=9),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=9,face="bold",colour="black"),
        legend.text=element_text(size=6.8),
        plot.title=element_text(size=10),
        legend.key.size=unit(0.5, "cm"),
        legend.title=element_text(size=7)) +
  ylim(0.0, 0.025)
#------------------------------------------------------------------------------
#GATC motifs for target bacterial classes
lineGATCgg <- bacteriaData %>% group_by(Class, ntRelPosition) %>%
  filter(Class %in% targetClass) %>%
  sample_n(10000, replace = TRUE) %>%
  summarise(ATmotifavg = mean(ATmotif),
            CpGmotifavg = mean(CpGmotif),
            CpGnormavg = mean(CpGnorm),
            GATCmotifavg = mean(GATCmotif),
            GCcontentavg = mean(GCcontent),
            GpCnormavg = mean(GpCnorm),
            observedExpectedAvg = mean(observedExpected)) %>%
  ggplot(aes(x = ntRelPosition, color = Class)) +
  #geom_line(aes(y = CpGmotifavg), position = "jitter", linetype = "dashed") +
  geom_line(aes(y = GATCmotifavg), position = "jitter", linetype = "solid") +
  xlab("Nucleotide position relative to TSS") +
  ylab("Average weighted motif frequency") +
  ggtitle("GATC motif, bootstrap n = 10000") +
  scale_color_manual(values = cbPalette) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title=element_text(size=9),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=9,face="bold",colour="black"),
        legend.text=element_text(size=6.8),
        plot.title=element_text(size=10),
        legend.key.size=unit(0.5, "cm"),
        legend.title=element_text(size=7)) +
  ylim(0.0, 0.05)
#------------------------------------------------------------------------------
#CpG motifs for target archaeal classes
archaeaCpGgg <- bacteriaData %>% group_by(Class, ntRelPosition) %>%
  filter(Class %in% archaeaClass) %>%
  sample_n(10000, replace = TRUE) %>%
  summarise(ATmotifavg = mean(ATmotif),
            CpGmotifavg = mean(CpGmotif),
            CpGnormavg = mean(CpGnorm),
            GATCmotifavg = mean(GATCmotif),
            GCcontentavg = mean(GCcontent),
            GpCnormavg = mean(GpCnorm),
            observedExpectedAvg = mean(observedExpected)) %>%
  ggplot(aes(x = ntRelPosition, color = Class)) +
  #geom_line(aes(y = CpGmotifavg), position = "jitter", linetype = "dashed") +
  geom_line(aes(y = CpGmotifavg), position = "jitter", linetype = "solid") +
  xlab("Nucleotide position relative to TSS") +
  ylab("Average weighted motif frequency") +
  ggtitle("Archaeal CpG motif, bootstrap n = 10000") +
  scale_color_manual(values = cbPalette) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title=element_text(size=9),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=9,face="bold",colour="black"),
        legend.text=element_text(size=6.8),
        plot.title=element_text(size=10),
        legend.key.size=unit(0.5, "cm"),
        legend.title=element_text(size=7)) +
  ylim(0.0, 0.025)
#------------------------------------------------------------------------------
#Normalized CpG motifs for target bacterial classes
classCpGnorm <- bacteriaData %>% group_by(Class, ntRelPosition) %>%
  filter(Class %in% targetClass) %>%
  #Bootstrap sample
  #sample_n(10000, replace = TRUE) %>%
  summarise(ATmotifavg = mean(ATmotif),
            CpGmotifavg = mean(CpGmotif),
            CpGnormavg = mean(CpGnorm),
            GATCmotifavg = mean(GATCmotif),
            GCcontentavg = mean(GCcontent),
            GpCnormavg = mean(GpCnorm),
            observedExpectedAvg = mean(observedExpected)) %>%
  ggplot(aes(x = ntRelPosition, color = Class)) +
  geom_line(aes(y = CpGnormavg), position = "jitter", linetype = "solid") +
  xlab("Nucleotide position relative to TSS") +
  ylab("Average weighted motif frequency") +
  ggtitle("CpG normalized") +
  scale_color_manual(values = cbPalette) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title=element_text(size=9),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=9,face="bold",colour="black"),
        legend.text=element_text(size=6.8),
        plot.title=element_text(size=10),
        legend.key.size=unit(0.5, "cm"),
        legend.title=element_text(size=7)) +
  ylim(0, 0.0025)
#------------------------------------------------------------------------------
#Normalized CpG motifs for target bacterial phyla
phylumCpGnorm <- bacteriaData %>% group_by(Phylum, ntRelPosition) %>%
  filter(Phylum %in% targetPhyla) %>%
  #Bootstrap sample
  sample_n(10000, replace = TRUE) %>%
  summarise(ATmotifavg = mean(ATmotif),
            CpGmotifavg = mean(CpGmotif),
            CpGnormavg = mean(CpGnorm),
            GATCmotifavg = mean(GATCmotif),
            GCcontentavg = mean(GCcontent),
            GpCnormavg = mean(GpCnorm),
            observedExpectedAvg = mean(observedExpected)) %>%
  ggplot(aes(x = ntRelPosition, color = Phylum)) +
  geom_line(aes(y = CpGnormavg), position = "jitter", linetype = "solid") +
  xlab("Nucleotide position relative to TSS") +
  ylab("Average weighted motif frequency") +
  ggtitle("CpG normalized, bootstrap n = 10000") +
  scale_color_manual(values = cbPalette) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title=element_text(size=9),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=9,face="bold",colour="black"),
        legend.text=element_text(size=6.8),
        plot.title=element_text(size=10),
        legend.key.size=unit(0.5, "cm"),
        legend.title=element_text(size=7)) +
  ylim(0, 0.005)
#------------------------------------------------------------------------------
#Observed Expected
classCpGOE <- bacteriaData %>% group_by(Class, ntRelPosition) %>%
  filter(Class %in% targetClass) %>%
  #Bootstrap sample
  sample_n(10000, replace = TRUE) %>%
  summarise(ATmotifavg = mean(ATmotif),
            CpGmotifavg = mean(CpGmotif),
            CpGnormavg = mean(CpGnorm),
            GATCmotifavg = mean(GATCmotif),
            GCcontentavg = mean(GCcontent),
            GpCnormavg = mean(GpCnorm),
            observedExpectedAvg = mean(observedExpected)) %>% 
  ggplot(aes(x = ntRelPosition, color = Class)) +
  geom_line(aes(y = observedExpectedAvg), position = "jitter", linetype = "solid") +
  xlab("Nucleotide position relative to TSS") +
  ylab("Average weighted motif frequency") +
  ggtitle("CpG Observed/Expected, bootstrap n = 10000") +
  scale_color_manual(values = cbPalette) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title=element_text(size=9),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=9,face="bold",colour="black"),
        legend.text=element_text(size=6.8),
        plot.title=element_text(size=10),
        legend.key.size=unit(0.5, "cm"),
        legend.title=element_text(size=7)) +
  ylim(0, 0.05)
#------------------------------------------------------------------------------
#CpG Observed/Expected ratio by Phylum
phylumCpGOE <- bacteriaData %>% group_by(Phylum, ntRelPosition) %>%
  filter(Phylum %in% targetPhyla) %>%
  #Bootstrap sample
  sample_n(10000, replace = TRUE) %>%
  summarise(ATmotifavg = mean(ATmotif),
            CpGmotifavg = mean(CpGmotif),
            CpGnormavg = mean(CpGnorm),
            GATCmotifavg = mean(GATCmotif),
            GCcontentavg = mean(GCcontent),
            GpCnormavg = mean(GpCnorm),
            observedExpectedAvg = mean(observedExpected)) %>% 
  ggplot(aes(x = ntRelPosition, color = Phylum)) +
  geom_line(aes(y = observedExpectedAvg), position = "jitter", linetype = "solid") +
  xlab("Nucleotide position relative to TSS") +
  ylab("Average weighted motif frequency") +
  ggtitle("CpG Observed/Expected, bootstrap n = 10000") +
  scale_color_manual(values = cbPalette) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title=element_text(size=9),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=9,face="bold",colour="black"),
        legend.text=element_text(size=6.8),
        plot.title=element_text(size=10),
        legend.key.size=unit(0.5, "cm"),
        legend.title=element_text(size=7)) +
  ylim(0, 0.05)
#------------------------------------------------------------------------------
#GC by phyla
phylumGC <- bacteriaData %>% group_by(Phylum, ntRelPosition) %>%
  filter(Phylum %in% targetPhyla) %>%
  sample_n(10000, replace = TRUE) %>%
  summarise(ATmotifavg = mean(ATmotif),
            CpGmotifavg = mean(CpGmotif),
            CpGnormavg = mean(CpGnorm),
            GATCmotifavg = mean(GATCmotif),
            GCcontentavg = mean(GCcontent),
            GpCnormavg = mean(GpCnorm),
            observedExpectedAvg = mean(observedExpected)) %>%
  ggplot(aes(x = ntRelPosition, color = Phylum)) +
  geom_line(aes(y = GCcontentavg), position = "jitter") +
  xlab("Nucleotide position relative to TSS") +
  ylab("Average weighted motif frequency") +
  ggtitle("GC content, bootstrap n = 10000") +
  scale_color_manual(values = cbPalette) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title=element_text(size=9),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=9,face="bold",colour="black"),
        legend.text=element_text(size=6.8),
        plot.title=element_text(size=10),
        legend.key.size=unit(0.5, "cm"),
        legend.title=element_text(size=7)) + ylim(0.0, 0.15)
#------------------------------------------------------------------------------
#assign ggplot grid objects
gp1 <- ggplotGrob(lineATgg)
gp2 <- ggplotGrob(lineGATCgg)
gp3 <- ggplotGrob(phylumGATC)
gp4 <- ggplotGrob(phylumAT)

#Adjust the widths
maxWidth = grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
#gp1$widths[2:7] <- as.list(maxWidth)
#gp2$widths[2:7] <- as.list(maxWidth)
#gp1$widths <- gp2$widths
#Arrange the grobs into one graphic
grid.arrange(gp1,gp2, nrow = 2,heights=c(5, 5))


maxWidth = grid::unit.pmax(gp3$widths[2:5], gp4$widths[2:5])
grid.arrange(gp3,gp4, nrow = 2,heights=c(5, 5))
#==============================================================================
####
#DNA curvature model data
####
#==============================================================================
#Initialize 2-core cluster
cluster <- create_cluster(4)
#Set default cluster size
set_default_cluster(cluster)
setwd("~/Documents/BINF868/")
curveFiles <- list.files("./testCurve")

# dnaCurve <- read.csv("./testCurve/NC_004252_cds0.csv")
# colnames(dnaCurve) <- c("Sequence","Index","Curvature","BendAngle",
#                         "CurvatureAngle","HelixAxisX","HelixAxisY","HelixAxisZ",
#                         "Phosphate1X","Phosphate1Y","Phosphate1Z","Phosphate2X",
#                         "Phosphate2Y","Phosphate2Z","BasepairNormalX","BasepairNormalY",
#                         "BasepairNormalZ","SmoothedNormalX","SmoothedNormalY","SmoothedNormalZ")

#accession <- gsub("(.*?)_cds*.csv$", "", file)
#cds <- gsub("(.*?)_CDS*.*$", "\\1", names(tempData[2]))
# 
# ntRelPosition <- c(-100:300)
# ntRelPosition <- ntRelPosition[ which(ntRelPosition != 0)]

#accession <- gsub("(NC_[0-9]{6})_cds*.csv$", "\\1", file)
# load_data <- function(path) { 
#   files <- dir(path, pattern = '\\.csv', full.names = TRUE)
#   tables <- lapply(files, read.csv)
#   do.call(rbind, tables)
# }
#dnaCurve <- load_data("/Users/imrambo/Documents/BINF868/testCurve/")
dnaCurve <- data.frame()
dnaCurve$accession <- vector(mode = "character")
dnaCurve$cds <- vector(mode = "numeric")


for (file in curveFiles){
  filePath <- paste("/Users/imrambo/Documents/BINF868/testCurve/", file, sep = "")
  #print(file)
  #print(filePath)

  #If the merged data frame does not exist, create it
  if (!exists("dnaCurve")){
    dnaCurve <- read.csv(filePath, header = TRUE)
    #print("creating data frame: dnaCurve")
  }
  #If the merged data frame exists, append to it
  if (exists("dnaCurve")){
    #print("dnaCurve exists")
    tempData <- read.csv(filePath, header = TRUE)
  
   #Get the accession number
  accession <- sub("(.*?)_cds[0-9]+.csv", "\\1", file)
  cds <- as.numeric(sub(".*?_[0-9]+_cds(.*?).csv", "\\1", file))
  #Position relative to TSS
  #ntRelPosition <- c(-100:nrow(tempData))
  #ntRelPosition <- ntRelPosition[which(ntRelPosition != 0)]
  
  colnames(tempData) <- c("Sequence","Index","Curvature","BendAngle",
                          "CurvatureAngle","HelixAxisX","HelixAxisY","HelixAxisZ",
                          "Phosphate1X","Phosphate1Y","Phosphate1Z","Phosphate2X",
                          "Phosphate2Y","Phosphate2Z","BasepairNormalX","BasepairNormalY",
                          "BasepairNormalZ","SmoothedNormalX","SmoothedNormalY","SmoothedNormalZ")
  
  #tempData$ntRelPosition <- ntRelPosition
  tempData$accession <- accession
  tempData$cds <- cds
  dnaCurve <- rbind(dnaCurve, tempData)
  rm(tempData)
  rm(accession)
  rm(cds)
  
  }
}
dnaCurveData <- save(dnaCurve, file = "dnaCurve.RData")


load("dnaCurve.RData")
#ntRelPosition <- c(-100:300)
#ntRelPosition <- ntRelPosition[which(ntRelPosition != 0)]
# 
# dnaCurve <- dnaCurve %>% filter(as.character(Sequence) != "Sequence") %>% 
#   #mutate(ntRelPosition = rep(ntRelPosition, nrow(dnaCurve)/400)) %>%
#   group_by(Index) %>%
#   sample_n(10000, replace = TRUE)
taxLineage <- cbind(bacteriaSummary, Genus) %>%
  filter(grepl("chromosome*", Replicon)) %>%
  rename(accession = Accession) %>%
  select(accession, Genus) %>%
  mutate(Genus = as.character(Genus)) %>%
  left_join(lineage) %>%
  filter(Class != "NA" | Class != " NA")

dnaCurveCount <- dnaCurve %>% count(Sequence, Index, accession) %>%
  group_by(Index, Sequence, accession) %>%
  rename(NTcount = n)

dnaCurveSum <- dnaCurveCount %>% group_by(Index, accession) %>% summarise(CDSsum = sum(NTcount))

dnaCurveMean <- dnaCurve %>% select(Sequence:CurvatureAngle, accession, cds) %>% 
  group_by(Index, accession, cds) %>%
  #sample_n(10000, replace = TRUE) %>%
  summarise(meanCurvature = mean(Curvature),
            meanBendAngle = mean(BendAngle),
            meanCurvatureAngle = mean(CurvatureAngle)) %>%
  left_join(taxLineage) %>%
  filter(Class != "NA") %>%
  left_join(dnaCurveCount) %>%
  left_join(dnaCurveSum) %>%
  mutate(NTpercent = NTcount/CDSsum) %>%
  left_join(bacteriaData)
  #sample_n(1000, replace = TRUE) %>%


#filter(dnaCurveMean, Index < 400 & Index > 175)
dnaCurveMean %>% left_join(ntData)

ggplot(dnaCurveMean, aes(x = Index)) +
  geom_line(aes(y = meanCurvatureAngle), color = "blue") +
  geom_line(aes(y = GCcontent), color = "firebrick") + ylim(0.1, 0.2) +
  ylab("Curvature") + xlab("Nucleotide position")
  #geom_point(aes(color = Sequence, fill = NTcount))
  #geom_line(aes(y = meanCurvatureAngle), color = "red")

dnaCurveMean %>% filter(as.character(Sequence) != "G" & as.character(Sequence) != "C") %>%
ggplot(aes(x = Index, color = Sequence)) +
  geom_line(aes(y = GCcontent)) +
  geom_line(aes(y = meanCurvature), color = "purple") + ylim(0.15, 0.4) +
  ylab("Percentage") + xlab("Nucleotide position")


#Partition the dnaCurve data frame
#dnaCurvePartition <- partition(dnaCurve, Index, cluster = cluster)


  #group_by(Index) %>%
  #summarise(meanCurve = mean(as.numeric(Curvature...0.0234..10.))) %>%
  #group_by(ntRelPosition) %>%
  #summarise(meanAngle = mean(as.numeric(Curvature...0.0234..10.))) %>%
# ggplot(aes(x = as.numeric(Index), y = as.numeric(meanCurve))) +
#   geom_line()
  ggplot(dnaCurve, aes(x = as.numeric(Index), y = as.numeric(Curvature...0.0234..10.))) +
  geom_line(color = "purple")
  

