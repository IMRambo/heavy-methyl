#Ian Rambo
#University of Delaware
#Oyster Rocks Methylation Project
#Last updated: October 7, 2015
#==============================================================================
library(dplyr)
library(ggdendro)
library(ggplot2)
library(tidyr)
library(grid)
library(gridExtra)
#Override the width of data frame columns printed in dplyr - print all
options(dplyr.width = Inf)
#==============================================================================
#GOAL: create heatmaps comparing conserved functions for target taxa across
#depth/age in Oyster Rocks sediment samples. Data are from an HpaII-digested
#Illumina library sequenced on an Illumina Hi-Seq 2500. CpG methylation scoring
#data are integrated to show cytosine methylation shifts in genes across depths.
#==============================================================================
# -----------               ------------
#  --------   \___________/   ---------
#    -------    O       O   ---------
#      -------\     V     /--------
#              \  _____  /        
#               \/     \/ 
#==============================================================================
setwd("/Users/imrambo/Documents/MS_thesis")

#FASTA headers for each sample
S0306 <- read.table("S03-06cm_pkhcAnnotated_01.headers.txt", header = FALSE,
                    sep = "|", quote = "")
S1215 <- read.table("S12-15cm_pkhcAnnotated_01.headers.txt", header = FALSE,
                    sep = "|", quote = "")
S2427 <- read.table("S24-27cm_pkhcAnnotated_01.headers.txt", header = FALSE,
                    sep = "|", quote = "")

S0306$Depth <- as.factor("03-06")
S1215$Depth <- as.factor("12-15")
S2427$Depth <- as.factor("24-27")

myVariables <- c("seqID","pSpecies","pScore","pGenus","genusConf","pFamily",
               "familyConf","pOrder","orderConf","pClass","classConf",
               "pPhylum","phylumConf","kPhylum","kClass","kOrder","kFamily",
               "kGenus","kSpecies","kSubSp1","kSubSp2","KO","KO_evalue",
               "KO_score","KO_gene","KO_name","KO_EC","COG","COG_evalue",
               "COG_score","COG_fcode","COG_name","COG_function", "Depth")

#Samples, All Depths data frame
sad <- rbind(S0306, S1215, S2427)
colnames(sad) <- myVariables

str(sad)

sadf <- sad %>%
  mutate(pClass = as.character(pClass),
         pOrder = as.character(pOrder),
         kClass = as.character(kClass),
         kOrder = as.character(kOrder),
         KO = as.character(KO),
         KO_gene = as.character(KO_gene),
         KO_name = as.character(KO_name),
         COG = as.character(COG),
         COG_fcode = as.character(COG_fcode)) %>%
  select(pSpecies:classConf, kClass, kOrder, kFamily, KO:Depth) %>%
  filter(classConf >= 0.700 &
           orderConf >= 0.650 & #Filter out low-scoring PhymmBL Order hits
           pClass != "<NA>" &
           pOrder != "<NA>" &
           pClass != "none" &
           pOrder != "none" &
           #Select for COG or KO
           #KO != "<NA>" &
           COG != "<NA>" &
           #KO_evalue <= 1e-4 &
           COG_evalue <= 1e-4)
#==============================================================================
#Input for Python NCBI target genomes data mining script get_genome_ncbi.py
#==============================================================================
#Abundant Classes of interest across all depths
  targetClass <- sadf %>% count(pClass, sort = TRUE) %>%
    filter(n > mean(n)) %>% print.data.frame()
#   
targetClass <- targetClass$pClass

#Target genuses - remove E. coli
targetGenus <- sadf %>% filter(pClass %in% targetClass &
                                 pGenus != "Escherichia") %>%
count(pGenus, sort = TRUE) %>%
filter(n > mean(n)) %>% as.data.frame()
targetGenusList <- as.data.frame(as.character(targetGenus$pGenus))
colnames(targetGenusList) <- "genus"
  
write.table(targetGenusList, "targetGenusRaw.txt", sep = "\t")
#sed and awk pipe to clean up the results
system("sed 's/\"//g' targetGenusRaw.txt | awk 'NR>1 {print $2}' > targetGenusClean.txt",
       intern = FALSE)
#==============================================================================
#3-6 cm
#==============================================================================
#Spread data frame - COG function code relative abundance
#for PhymmBL-called Classes, 3-6 cm
cfs0306 <- sadf %>% filter(Depth == "03-06") %>%
  group_by(pClass) %>%
  count(COG_fcode, pClass, sort = TRUE) %>%
  mutate(ra = n/sum(n)) %>%
  select(COG_fcode, pClass, ra) %>%
  spread(COG_fcode, ra, fill = 0) %>%
  as.data.frame() 

#Change rownames to taxa
rownames(cfs0306) <- cfs0306$pClass
#Remove pClass variable
cfs0306 <- select(cfs0306, BQ:X)

#Standardize the data (mean zero, unit variance) and convert to matrix
dd0306 <- as.matrix(scale(cfs0306))

ddhc0306.col <- as.dendrogram(hclust(dist(dd0306, method = "euclidean"), method = "complete"))
col0306.ord <- order.dendrogram(ddhc0306.col)

ddhc0306.row <- as.dendrogram(hclust(dist(t(dd0306), method = "euclidean")), method = "complete")
row0306.ord <- order.dendrogram(ddhc0306.row)

x0306 <- scale(cfs0306)[col0306.ord, row0306.ord]
x0306_names <- attr(x0306, "dimnames")
xdf0306 <- as.data.frame(x0306)

colnames(xdf0306) <- x0306_names[[2]]
xdf0306$taxa <- x0306_names[[1]]
xdf0306$taxa <- with(xdf0306, factor(taxa, levels = taxa, ordered = TRUE))
xdf0306t <- gather(xdf0306, taxa)
colnames(xdf0306t) <- c("taxa","fcode","value")

cleanTheme <- theme(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_blank(),
                    axis.line = element_line(color = "black"),
                    axis.title=element_text(size=8),
                    axis.text.x=element_text(size=7),
                    axis.text.y=element_text(size=9,face="bold",colour="black"),
                    legend.text=element_text(size=6.8),
                    plot.title=element_text(size=10),
                    legend.key.size=unit(0.3, "cm"),
                    legend.title=element_text(size=7))

#Heatmap fill color scale
hfillColor <- c("lavender","lavenderblush","lightcoral","indianred","purple4")

#Heatmap, target taxa, 3-6 cm
tthmp0306 <- xdf0306t %>% filter(as.character(taxa) %in% targetTaxa) %>%
  ggplot(aes(x = fcode, y = taxa, fill = value)) +
  geom_tile(colour = "white") + xlab("COG function code") +
  ylab("Class") + ggtitle("3-6 cm, COG functions") +
  scale_fill_gradientn(colours = hfillColor, name = "Abundance") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color="black"),
        axis.title=element_text(size=10),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=9, face="bold", color="black"),
        plot.title=element_text(size=10),
        legend.key.size=unit(1, "cm"),
        legend.text=element_text(size=7),
        legend.title=element_text(size=7))

#Heatmap, all taxa, 3-6 cm
athmp0306 <- ggplot(xdf0306t, aes(x = fcode, y = taxa, fill = value)) +
  geom_tile() + xlab("COG function code") +
  ylab("Class") + ggtitle("3-6 cm, COG functions") +
  scale_fill_gradientn(colours = hfillColor, name = "Value") +
  cleanTheme

#Extract dendrogram data for y-axis (taxa)
ddata_y_0306 <- dendro_data(ddhc0306.col, type = "rectangle")
#Create dendrogram of taxa
tdendro_y_0306 <- ggdendrogram(ddata_y_0306, rotate = TRUE)

#Create ggplot Grobs for the plots, which will be combined into one graphic
gp1 <- ggplotGrob(tthmp0306)
gp2 <- ggplotGrob(athmp0306)
gp3 <- ggplotGrob(tdendro_y_0306)

maxWidth = grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5], gp3$widths[2:5])
gp1$widths[2:5] <- as.list(maxWidth)
gp2$widths[2:5] <- as.list(maxWidth)
gp3$widths[2:5] <- as.list(maxWidth)

grid.arrange(gp2,gp3, ncol=2,heights=c(4/5,1/5))
#==============================================================================
#12-15 cm
#==============================================================================
cfs1215 <- sadf %>% filter(Depth == "12-15") %>%
  group_by(pClass) %>%
  count(COG_fcode, pClass, sort = TRUE) %>%
  mutate(ra = n/sum(n)) %>%
  select(COG_fcode, pClass, ra) %>%
  spread(COG_fcode, ra, fill = 0) %>%
  as.data.frame() 

#Change rownames to taxa
rownames(cfs1215) <- cfs1215$pClass
#Remove pClass variable
cfs1215 <- select(cfs1215, BQ:X)

#Standardize the data (mean zero, unit variance) and convert to matrix
dd1215 <- as.matrix(scale(cfs1215))

ddhc1215.col <- as.dendrogram(hclust(dist(dd1215, method = "euclidean"), method = "complete"))
col1215.ord <- order.dendrogram(ddhc1215.col)

ddhc1215.row <- as.dendrogram(hclust(dist(t(dd1215), method = "euclidean")), method = "complete")
row1215.ord <- order.dendrogram(ddhc1215.row)

x1215 <- scale(cfs1215)[col1215.ord, row1215.ord]
x1215_names <- attr(x1215, "dimnames")
xdf1215 <- as.data.frame(x1215)

colnames(xdf1215) <- x1215_names[[2]]
xdf1215$taxa <- x1215_names[[1]]
xdf1215$taxa <- with(xdf1215, factor(taxa, levels = taxa, ordered = TRUE))
xdf1215t <- gather(xdf1215, taxa)
colnames(xdf1215t) <- c("taxa","fcode","value")

#Extract dendrogram data for y-axis (taxa)
ddata_y_1215 <- dendro_data(ddhc1215.col, type = "rectangle")
#Create dendrogram of taxa
tdendro_y_1215 <- ggdendrogram(ddata_y_1215, rotate = TRUE)

#Heatmap, target taxa and COG function codes, 12-15 cm
tthmp1215 <- xdf1215t %>% filter(as.character(taxa) %in% targetTaxa) %>%
  ggplot(aes(x = fcode, y = taxa, fill = value)) +
  geom_tile(colour = "gray") + xlab("COG Function Code") +
  ylab("Class") + ggtitle("12-15 cm, COG functions") +
  scale_fill_gradientn(colours = hfillColor, name = "Abundance") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color="black"),
        axis.title=element_text(size=10),
        axis.text.x=element_text(size=7),
        axis.text.y=element_text(size=9, face="bold", color="black"),
        plot.title=element_text(size=10),
        legend.key.size=unit(1, "cm"),
        legend.text=element_text(size=7),
        legend.title=element_text(size=7))
#==============================================================================
#24-27 cm
#==============================================================================
cfs2427 <- sadf %>% filter(Depth == "24-27") %>%
  group_by(pClass) %>%
  count(COG_fcode, pClass, sort = TRUE) %>%
  mutate(ra = n/sum(n)) %>%
  select(COG_fcode, pClass, ra) %>%
  spread(COG_fcode, ra, fill = 0) %>%
  as.data.frame() 

#Change rownames to taxa
rownames(cfs2427) <- cfs2427$pClass
#Remove pClass variable
cfs2427 <- select(cfs2427, A:X)

#Standardize the data (mean zero, unit variance) and convert to matrix
dd2427 <- as.matrix(scale(cfs2427))

ddhc2427.col <- as.dendrogram(hclust(dist(dd2427, method = "euclidean"), method = "complete"))
col2427.ord <- order.dendrogram(ddhc2427.col)

ddhc2427.row <- as.dendrogram(hclust(dist(t(dd2427), method = "euclidean")), method = "complete")
row2427.ord <- order.dendrogram(ddhc2427.row)

x2427 <- scale(cfs2427)[col2427.ord, row2427.ord]
x2427_names <- attr(x2427, "dimnames")
xdf2427 <- as.data.frame(x2427)

colnames(xdf2427) <- x2427_names[[2]]
xdf2427$taxa <- x2427_names[[1]]
xdf2427$taxa <- with(xdf2427, factor(taxa, levels = taxa, ordered = TRUE))
xdf2427t <- gather(xdf2427, taxa)
colnames(xdf2427t) <- c("taxa","fcode","value")

#Extract dendrogram data for y-axis (taxa)
ddata_y_2427 <- dendro_data(ddhc2427.col, type = "rectangle")
#Create dendrogram of taxa
tdendro_y_2427 <- ggdendrogram(ddata_y_2427, rotate = TRUE)
#==============================================================================
#All depths
#==============================================================================
xdf0306t$depth <- as.factor("03-06 cm")
xdf1215t$depth <- as.factor("12-15 cm")
xdf2427t$depth <- as.factor("24-27 cm")
#The majority of COG functions at 3-6 cm are conserved downcore
fCodes <- levels(xdf0306t$fcode)
asdf_t <- rbind(xdf0306t, xdf1215t, xdf2427t)

#Heatmap - changes in abundance of conserved COG functions across depth/age
asdf_t %>% filter(as.character(taxa) %in% targetTaxa &
                    as.character(fcode) %in% fCodes) %>%
  ggplot(aes(x = fcode, y = taxa, fill = value)) +
  geom_tile(colour = "gray") +
  xlab("COG Function Code") +
  ylab("Class, PhymmBL") +
  scale_fill_gradientn(colours = hfillColor, name = "Abundance") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color="black"),
        axis.title=element_text(size=10),
        axis.text.x=element_text(size=8, face="bold",
                                 angle=90, color="black"),
        axis.text.y=element_text(size=9, color="black"),
        plot.title=element_text(size=10),
        legend.key.size=unit(1, "cm"),
        legend.text=element_text(size=7),
        legend.title=element_text(size=6.5)) +
  facet_wrap(~depth, nrow=3)