#Oyster Rocks Methylation Project
#Ian Rambo
#University of Delaware
#Funding: CDEBI
#Sediment porewater ion chromatography
#==============================================================================
#GOALS
#Graph sediment properties:
#Porosity
#Ion Chromatography - Mg, Ca, K, Na, Br, Cl, SO4
#==============================================================================
library(ggplot2)
library(grid)
library(tidyr)
library(plyr)
library(dplyr)
#==============================================================================
setwd("~/Documents/MS_thesis/ion_chromatography/OR_ICdata/")
IC_files <- list.files()

sampleFiles <- IC_files[grep("^S |^L ", IC_files)]
standardFiles <- IC_files[grep("Standard", IC_files)]

#==============================================================================
#Standards
#==============================================================================
#Empty data frame that will contain standards for K, NH4, Na, SO4, and Cl
ICSTD <- data.frame()
#Read in all standards csv files, rbind data frames
for (file in 1:length(standardFiles)){
  csv <- read.csv(standardFiles[file], header = TRUE, sep = ";")
  ICSTD <- rbind(csv, ICSTD)
  ICSTD <- filter(ICSTD, Determination.start != " ")   
}

ICSTD$Dilution <- integer(length = nrow(ICSTD))
ICSTD$Dilution[grep("1:2 Standard", ICSTD$Ident)] <- 2
ICSTD$Dilution[grep("1:4 Standard", ICSTD$Ident)] <- 4
ICSTD$Dilution[grep("1:20 Standard", ICSTD$Ident)] <- 20
ICSTD$Dilution[grep("1:40 Standard", ICSTD$Ident)] <- 40
ICSTD$Dilution[grep("1:200", ICSTD$Ident)] <- 200
ICSTD$Dilution[grep("1:400", ICSTD$Ident)] <- 400

ICSTD$Concentration <- c(0.05, 0.5, 5, 0.01, 1, 10)
#Standard curve linear regression models - K, NH4, Na, SO4, Cl
K.lm <- lm(data = ICSTD,  Cations.Potassium.Area ~ Concentration)
NH4.lm <- lm(data = ICSTD, Cations.Ammonium.Area ~ Concentration)
Na.lm <- lm(data = ICSTD, Cations.Sodium.Area ~ Concentration)
SO4.lm <- lm(data = ICSTD, Anions.Sulfate.Area ~ Concentration)
Cl.lm <- lm(data = ICSTD, Anions.Chloride.Area ~ Concentration)
NH4.lm <- lm(data = ICSTD, Cations.Ammonium.Area ~ Concentration)
#csv file for Mg, Ca, and Br standards
MgCaBr <- read.csv("Mg-Ca-Br_stds.csv", header = TRUE)

#Standard curve linear regression models - Mg, Ca, Br
Mg.lm <- lm(data = MgCaBr, Mg.area ~ concentration_uM)
Ca.lm <- lm(data = MgCaBr, Ca.area ~ concentration_uM)
Br.lm <- lm(data = MgCaBr, Br.area ~ concentration_uM)

#Calculate concentration from a standard curve
calcConc <- function(area, b, m, d){((area - b)/m) * d}
#==============================================================================
#Samples
#==============================================================================
#Empty data frame that will contain sample data
ICLS <- data.frame()
#Read in all sample csv files, rbind data frames
for (file in 1:length(sampleFiles)){
  csv <- read.csv(sampleFiles[file], header = TRUE, sep = ";")
  ICLS <- rbind(csv, ICLS)
  ICLS <- filter(ICLS, Determination.start != " ")   
}

ICLS <- select(ICLS, -Determination.start)

ICLS$Core <- NA
ICLS$Core[grep("^S ", ICLS$Ident)] <- "S"
ICLS$Core[grep("^L ", ICLS$Ident)] <- "L"
ICLS$Core <- as.factor(ICLS$Core)

ICLS$Dilution <- integer(length = nrow(ICLS))
ICLS$Dilution[grep("1:100", ICLS$Ident)] <- 101
ICLS$Dilution[grep("1:10$", ICLS$Ident)] <- 11
ICLS$Dilution[grep("1:5", ICLS$Ident)] <- 6

ICLS$Depth <- integer(length = nrow(ICLS))
ICLS$Depth[grep("0-3", ICLS$Ident)] <- 3
ICLS$Depth[grep("3-6", ICLS$Ident)] <- 6
ICLS$Depth[grep("6-9", ICLS$Ident)] <- 9
ICLS$Depth[grep("9-12", ICLS$Ident)] <- 12
ICLS$Depth[grep("12-15", ICLS$Ident)] <- 15
ICLS$Depth[grep("15-18", ICLS$Ident)] <- 18
ICLS$Depth[grep("18-21", ICLS$Ident)] <- 21
ICLS$Depth[grep("21-24", ICLS$Ident)] <- 24
ICLS$Depth[grep("24-27", ICLS$Ident)] <- 27

#Create a tidied, grouped data frame with ion concentrations using a dplyr pipe
ICLS_tidy <- ICLS %>% mutate(Core = as.character(Core)) %>%
  group_by(Core, Dilution, Depth) %>%
  arrange(Core, Dilution, Depth) %>%
  transmute(
    Mg = calcConc(area = Cations.Magnesium.Area,
                               m = Mg.lm$coefficients[2],
                               b = Mg.lm$coefficients[1],
                               d = Dilution),
    Ca = calcConc(area = Cations.Calcium.Area,
                               m = Ca.lm$coefficients[2],
                               b = Ca.lm$coefficients[1],
                               d = Dilution),
    K = calcConc(area = Cations.Potassium.Area,
                     m = K.lm$coefficients[2],
                     b = K.lm$coefficients[1],
                     d = Dilution),
    Br = calcConc(area = Anions.Bromide.Area,
                      m = Br.lm$coefficients[2],
                      b = Br.lm$coefficients[1],
                      d = Dilution),
    SO4 = calcConc(area = Anions.Sulfate.Area,
                       m = SO4.lm$coefficients[2],
                       b = SO4.lm$coefficients[1],
                       d = Dilution),
    Cl = calcConc(area = Anions.Chloride.Area,
                      m = Cl.lm$coefficients[2],
                      b = Cl.lm$coefficients[1],
                      d = Dilution),
    Na = calcConc(area = Cations.Sodium.Area,
                  m = Na.lm$coefficients[2],
                  b = Na.lm$coefficients[1],
                  d = Dilution),
    NH4 = calcConc(area = Cations.Ammonium.Area,
                   m = NH4.lm$coefficients[2],
                   b = NH4.lm$coefficients[1],
                   d = Dilution)
    ) %>%
  gather(Ion, Concentration, Mg:NH4) 

#Faceted scatter plots for ions, cores L and S
icgg <- ggplot(ICLS_tidy, aes(x = Concentration, y = Depth))
icgg + geom_point(aes(color = Core, shape = factor(Dilution)), size = 2) +
  theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
        axis.line = element_line(color = "black")) +
  scale_y_continuous(breaks = seq(0, 27, 3), trans = "reverse") +
  facet_wrap(~ Ion, nrow = 4, scales = "free", ) + ylab("Depth (cm)") +
  xlab(expression(Concentration~(mu~M)))

#Filter for best dilutions from core L
ICLS_tidy_filter <- rbind(filter(ICLS_tidy, Core == "L" & Ion == "Mg", Dilution == 101),
      filter(ICLS_tidy, Core == "L" & Ion == "Ca", Dilution == 6),
      filter(ICLS_tidy, Core == "L" & Ion == "K", Dilution == c(11, 6)),
      filter(ICLS_tidy, Core == "L" & Ion == "Br", Dilution == c(11, 6)),
      filter(ICLS_tidy, Core == "L" & Ion == "SO4", Dilution == 11),
      filter(ICLS_tidy, Core == "L" & Ion == "Cl", Dilution == 101)) %>%
  na.omit()

#Colorblind-friendly palette
cbPalette <- c("deepskyblue","darkorange2","pink","aquamarine2","darkorchid2",
               "dodgerblue3","chartreuse2","firebrick2","mediumvioletred","brown","black","green")
#Plot faceted scatterplots of best dilutions from core L
ICLS_tidy_filter %>% mutate(Ion = as.character(Ion)) %>%
                              ggplot(aes(x = Concentration, y = Depth)) +
  geom_point(aes(color = factor(Dilution), shape = factor(Dilution)), size = 3) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black")) +
  scale_y_continuous(breaks = seq(0, 27, 3), trans = "reverse") +
  #facet_wrap(~ Ion, nrow = 4, scales = "free", ) + ylab("Depth (cm)") +
  xlab(expression(Concentration~(mu~M))) +
  scale_color_manual(name = "Dilution", values = c("purple", "indianred", "blue")) +
  scale_shape_manual(name = "Dilution", values = c(15, 17, 19)) +
  #geom_smooth(aes(group = Ion), method = "lm", alpha = 0.2)

ass <- ICLS_tidy %>% mutate(Ion = as.character(Ion)) %>%
  filter(Core == "L" & Ion == "SO4") %>%
  filter(Dilution == 11 | Dilution == 6) %>%
  ggplot(aes(x = Concentration, y = Depth, shape = as.factor(Dilution), color = as.factor(Dilution))) +
  geom_point(size = 2) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black")) +
  scale_y_continuous(breaks = seq(0, 27, 3), trans = "reverse") +
  #facet_wrap(~ Ion, nrow = 4, scales = "free", ) + ylab("Depth (cm)") +
  xlab(expression(Concentration~(mu~M))) +
  ylab("Depth (cm)") +
  ggtitle("Sulfate concentration, Core L") +
  scale_color_manual(name = "Dilution", values = c("purple", "orange")) +
  scale_shape_manual(name = "Dilution", values = c(15, 17))
#==============================================================================
#Sediment porosity

#----POROSITY FUNCTIONS----
poro <- function(Mw, Md, S, Pds){
  #Calculate sediment porosity
  #Mw = wet mass
  #Md = dry mass
  #1 = density of pure water
  #Pds = density of dry sediment
  #S = salinity
  (Mw/1)/((Mw/1) + (((Md - S) * (Mw/1000))/Pds))}

cylinder_volume <- function(r, h){
  #Calculate the volume of a cylinder
  pi * r^2 * h
}
#-------------------------

poreMass <- read.csv("~/R/ORMP/data/porosity_sediment.csv", header = TRUE)
poreMass <- mutate(poreMass, core = as.factor(core),
                   totalVolume = vector(mode = "numeric", length = nrow(poreMass)),
                   porosityA = vector(mode = "numeric", length = nrow(poreMass)),
                   porosityB = vector(mode = "numeric", length = nrow(poreMass)))

 
poreMass$totalVolume[which(poreMass$core == "R" & poreMass$depth <= 10)] <- cylinder_volume(5.08, 1)
poreMass$totalVolume[which(poreMass$core == "R" & poreMass$depth > 10)] <- cylinder_volume(5.08, 2)
poreMass$totalVolume[which(poreMass$core == "L")] <- cylinder_volume(0.6, 2.7)
poreMass$totalVolume[which(poreMass$core == "S")] <- cylinder_volume(0.6, 2.7)


#Calculate porosity and add new vectors
poreMass <- poreMass %>% mutate(density = dry/totalVolume,
                                water_loss = wet - dry, 
                    porosityA = poro(wet, dry, 26, density),
                    porosityB = water_loss/totalVolume)

poreAL_lm <- lm(depth~porosityA, data = filter(poreMass, core == "L"))
poreAR_lm <- lm(depth~porosityA, data = filter(poreMass, core == "R"))
#Labels for linear models
lbAL <- paste("R^2 == ", round(summary(poreAL_lm)$r.squared,4))
lbAR <- paste("R^2 == ", round(summary(poreAR_lm)$r.squared,4))

poreBL_lm <- lm(depth~porosityB, data = filter(poreMass, core == "L"))
poreBR_lm <- lm(depth~porosityB, data = filter(poreMass, core == "R"))
#Labels for linear models
lbBL <- paste("R^2 == ", round(summary(poreBL_lm)$r.squared,4))
lbBR <- paste("R^2 == ", round(summary(poreBR_lm)$r.squared,4))

#Equation of linear model
eqBR <- paste("y = ", round(summary(poreBR_lm)[[4]][2],3),"x - ", round(abs(summary(poreBR_lm)[[4]][1]), 3), sep = "")
#Porosity A
poreMass %>% filter(core == "R") %>%
  ggplot(aes(x = porosityA, y = depth, color = core, shape = core)) +
  geom_point(size = 4) +
  scale_y_reverse(breaks = seq(0, 34, 2)) +
  scale_color_manual(values = "deepskyblue") +
  guides(color = FALSE) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title=element_text(size=11),
        axis.text.x=element_text(size=9),
        axis.text.y=element_text(size=9,face="bold",colour="black"),
        legend.text=element_text(size=8),
        plot.title=element_text(size=13),
        legend.key.size=unit(0.8, "cm"),
        legend.title=element_text(size=10),
        legend.position = "none") +
  ylab("Depth (cm)") +
  xlab("Porosity") +
#   geom_smooth(method = "lm", se = FALSE) +
#     annotate("text", y = 17, x = 0.95, label = lbAR, parse=TRUE,
#              size = 6.5, color = "black") +
#   annotate("text", y = 16, x = 0.95, label = eqBR, size = 6.5) +
  ggtitle("Porosity, Core R")

#Porosity B 
poreMass %>% filter(core == "S") %>%
  ggplot(aes(x = porosityB, y = depth, color = core, shape = core)) +
  geom_point(size = 4) +
  scale_y_reverse(breaks = seq(1, 34, 2)) +
  scale_color_manual(values = c("deepskyblue","mediumvioletred")) +
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
        legend.title=element_text(size=7)) +
  ylab("Depth (cm)") +
  xlab("Porosity") +
  geom_smooth(method = "lm", se = FALSE) +
  annotate("text", y = 17, x = 0.80, label = lbBR, parse=TRUE,
           size = 5, color = "mediumvioletred") 
#   annotate("text", y = 17, x = 1.12, label = lbBL, parse=TRUE,
#            size = 5, color = "deepskyblue")


