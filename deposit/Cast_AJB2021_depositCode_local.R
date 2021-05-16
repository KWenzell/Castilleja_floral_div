#Floral trait and genetic divergence in Castilleja sessiliflora and C. purpurea species complex
#Floral trait data- cleaned and averaged by individual
library(tidyverse)
library(vegan)
library(stringr)
library(ggplot2)
library(ggpubr)
library(geodist)
library(dplyr)

setwd("~/CASE_Rdocs/Cast_AJB2021/deposit")
#all floral trait data:
ca <- read.csv("CastFloral_v2-coreTraits_corrected9-7-20.csv")



#compilation of code for making figs for AJB resubmission: 
#Floral trait and genetic distance in Castilleja sessiliflora and C. purpurea

#FIGURE 2
#FLORAL NMDS AND GENETIC PCA
#A) Genetic PCA of SNPs:

# read in data
pca <- read_table2("./CaInd.eigenvec", col_names = FALSE)
eigenval <- scan("./CaInd.eigenval")

# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

pca

# sort out the individual species and pops
# spp
spp <- rep(NA, length(pca$ind))
# spp[grep("indivisa", pca$ind)] <- "Indivsa"
# spp[grep("Outgroup_2010", pca$ind)] <- "2010"
spp[grep("citrina", pca$ind)]  <- "Cit"
spp[grep("purpurea", pca$ind)]  <- "Purp"
spp[grep("lindheimeri", pca$ind)]  <- "Lind"
spp[grep("lindheimeriXpurpurea", pca$ind)]  <- "Hyb"
spp[grep("Sessiliflora-NE", pca$ind)]  <- "SessNE"
spp[grep("Sessiliflora-S", pca$ind)]  <- "SessS"
spp[grep("Sessiliflora-W", pca$ind)]  <- "SessW"



# Genus
gen <- rep(NA, length(pca$ind))
# gen[grep("Outgroup", pca$ind)] <- "Outgroup"
gen[grep("_P-", pca$ind)]  <- "Purp"
gen[grep("Sessiliflora", pca$ind)]  <- "Sess"

# combine - if you want to plot each in different colours
spp_gen <- paste0(spp, "_", gen)

# remake data.frame
pca <- as_tibble(data.frame(pca, spp, gen, spp_gen))

# Find out number of rows
pca


# first convert to percentage variance explained
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

summary(pca)

# plot pca
b <- ggplot(pca, aes(PC1, PC2, col = spp, shape = gen)) + geom_point(size = 3)
# b <- b + scale_colour_manual(values = c("red", "blue"))
b <- b + coord_equal() + theme_light()
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))


#export pca as csv 
#want to visualize specific pops, add additional population data
pca <- as.data.frame(pca)

write.csv(pca, "CaInd_hybseq_SNP_pca_cleaned.csv")


#added pop column and consistent Region, and morph columns
#read this in for further analyses
pcaX <- read.csv("CaInd_hybseq_SNP_pca_cleaned_mod2021-02-03.csv")


#color code by floral morph
snpPca <- ggplot(pcaX, aes(PC1, PC2, col = morph, shape = sp, label= reg2)) + 
  geom_point(size = 3)+
  #geom_text(aes(label= reg2), hjust= -0.5)+
  theme_classic()+
  xlab(paste0("PC1 (6.55%)")) +
  ylab(paste0("PC2 (5.55%)"))+
  scale_y_reverse()+
  #scale_x_reverse()+ #to mirror floral NMDS
  stat_ellipse(aes(PC1, PC2, col = reg2), level= 0.9)+
  scale_colour_manual(values = c("C" = "#FFCC00", "P" = "#9900FF", "L" = "#FF6633", "typ" = "#669900", 
                                 "divY" = "#e6e600", "divP" = "#FF6699", "nS" = "#003300", "wS"= "#666600", "sS"= "#CC6666"),
                      labels= c("C" = "C. citrina", "P" = "C. purpurea", "L" = "C. lindheimeri", "typ" = "C. sess- typical", 
                                "divY" = "C. sess- div-yellow", "divP" = "C. sess- div-pink", "nS" = "C. sess- North", 
                                "wS"= "C. sess-West", "sS"= "C. sess- South"))+
  scale_shape_manual(values= c("S"= 16, "P"= 15, "C"= 18, "L"= 17))

#Fig 2A:
snpPca
#####

#Floral NMDS code in file multivariate.R (lines 11-306) (didn't need to regenerate for resubmission, just recombine with new genetic PCA)
#Figure 2B- Floral trait NMDS

ca <- read.csv("CastFloral_v2-coreTraits_corrected9-7-20.csv")


View(ca)
#NMDS with all data-
#v2= pared down floral traits (removed corW1, calLbW, brL)
#removed traits were highly correlated with other remaining traits

#library(vegan)

#columns with useful traits: 6-10 (morph); 11-13 (color)

#start with all
c.mds <- metaMDS(ca[ ,6:13], distance= "gower", zerodist="add")
c.mds
#stress= 0.19

#add values to dataframe for each individual
ca$nmds1 <- c.mds$points[,1]
ca$nmds2 <- c.mds$points[,2]

stressplot(c.mds)

#Visualize NMDS with hex colors

#use hex colors from individuals
#use hexadecimals for color values

#library(stringr)
rgb <- ca[,11:13]
#function rgb_to_hex: https://rdrr.io/github/stuart-morrison/schemr/man/rgb_to_hex.html
rgb_to_hex <- function(rgb) {
  
  # Extract colour channels from first three colummns
  red <- unlist(rgb[ , 1])
  green <- unlist(rgb[ , 2])
  blue <- unlist(rgb[ , 3])
  
  # Test if any RGB values are outside of [0, 255]
  if (any(red < 0 | red > 255) |
      any(green < 0 | green > 255) |
      any(blue < 0 | blue > 255)) {
    stop("Colour channels should be between 0 and 255.")
  }
  
  # Convert colour channels to hex
  hex_value <- paste0("#",
                      # Convert red
                      str_pad(string = as.hexmode(round(red, 0)),
                              width = 2, pad = "0", side = "left"),
                      # Convert green
                      str_pad(string = as.hexmode(round(green, 0)),
                              width = 2, pad = "0", side = "left"),
                      # Convert blue
                      str_pad(string = as.hexmode(round(blue, 0)),
                              width = 2, pad = "0", side = "left")
  )
  
  return(hex_value)
}
ca$hex <- rgb_to_hex(rgb)
#it worked

#shape code points by species
ca$pch_ind[ca$sp == "P"] <- 15 #square
ca$pch_ind[ca$sp == "S"] <- 16 #circle
ca$pch_ind[ca$sp == "C"] <- 18 #diamond
ca$pch_ind[ca$sp == "L"] <- 17 #triangle

#final NMDS figure/Fig 2B:

#plotted with hex inflor colors
fig <- ordiplot(c.mds, type = "none")
points(fig, "sites", col=ca$hex, pch= ca$pch_ind,
       bg="white", cex=1.3)
legend("topleft",  inset= c(0.1,0.05),    # location and inset- still not great
       bty="n", cex= 1,
       title="Species", xpd= TRUE,
       c("C. sessiliflora", "C. purpurea", "C. lindheimeri", "C. citrina"),
       pch= c(16, 15, 17, 18),  pt.cex= 1.5,
       pt.bg= "black")

#plot st deviation for CAPU sp and CASE regions
for(i in unique(ca$reg2)) {
  ordiellipse(c.mds$point[grep(i,ca$reg2),],draw="polygon", alpha = 50, kind= "sd",
              groups=ca$reg2[ca$reg2==i],col=ca$col_index2[grep(i,ca$reg2)], label= TRUE)
} 


#envfit
#add trait loadings

#for floral traits
fit <- envfit(c.mds, ca[, 6:13], perm = 1000)
fit
plot(fit, p.max = 0.05, col = "black", cex=1.2)

#^Figure 2B
#####


#FIGURE 3: FLORAL TRAIT DISTANCE AND WITHIN GROUP GENETIC DISTANCE

#calculate trait distances within groups (taxa/species/region)
#GOWER distance method using vegdist function, followed by betadisper

#morph and color
distg <- vegdist (ca[ ,6:13], method= "gower")

#gower by taxa, 
#vegdist using all traits
dispgt <- betadisper(distg, ca$taxa, type= "median")
dispgt

#compare average distance by groups
anova(dispgt) #3.976e-06 ***, f=22
permutest(dispgt, pairwise = TRUE, permutations = 99) #.01
dispgt.HSD <- TukeyHSD(dispgt)
dispgt.HSD
plot(dispgt.HSD)


#gower by sp
#vegdist, all traits
dispg <- betadisper(distg, ca$sp, type= "median")
dispg

anova(dispg) # < 2.2e-16 ***, f=117
permutest(dispg, pairwise = TRUE, permutations = 99)
dispg.HSD <- TukeyHSD(dispg)
dispg.HSD
plot(dispg.HSD)


#gower by region
dispgr <- betadisper(distg, ca$reg2, type= "median")
dispgr

anova(dispgr) #< 2.2e-16 ***, f=30
permutest(dispgr, pairwise = TRUE, permutations = 99)
dispgr.HSD <- TukeyHSD(dispgr)
dispgr.HSD
plot(dispgr.HSD)

#betadisper color by pop- 
#Calculate color variation (distance to median) within populations

distgCol <- vegdist(ca[ ,11:13], method= "gower")
dispgColPop <- betadisper(distgCol, ca$pop, type= "median")
dispgColPop

#make dataframe with betadisper Color pop values, plus region and sp info
#copied to above w/ other betadisper figs
bdPop <- data.frame(distToMedPop= dispgColPop$distances, Pop= dispgColPop$group, Reg2= ca$reg2, Sp= ca$sp)


bdPC <- ggline(bdPop, x = "Reg2", y = "distToMedPop", plot_type= "p",
               add = c("mean_se"), error.plot = "pointrange",
               order = c("C", "L", "P","sS", "wS", "nS"),
               ylim= c(0.05,0.19), #size = 1,
               ylab = "Avg Within-Population Color Distance", xlab = "Region",
               font.y = list(size = 14, face= "plain"),
               font.x = list(size = 14, face= "plain"))+
  stat_compare_means(method= "anova", label.y= 0.18)+
  stat_compare_means(label= "p.signif", method= "t.test", ref.group = ".all.", label.y= 0.16)+
  geom_hline(yintercept = mean(bdPop$distToMedPop), linetype = 2)

bdPC

poly <- aov(formula= distToMedPop ~ Reg2, data= bdPop)
summary(poly)
poly.HSD <- TukeyHSD(poly)
poly.HSD
plot(poly.HSD)

#Plot betadisper values among groups (Fig 3D-F)

# by taxa
bdTaxa <- data.frame(distToMed= dispgt$distances, Taxa= dispgt$group)

bdT <- ggline(bdTaxa, x = "Taxa", y = "distToMed", plot_type= "p",
              add = c("mean_se"), error.plot = "pointrange",
              order = c("PC", "S"),
              ylim= c(0.0,0.3),
              ylab = "Floral trait variation (Gower distance)", xlab = "Taxa",
              #font.y = list(size = 14, face= "plain"),
              font.x = list(size = 14, face= "plain"))+
  stat_compare_means()+#label.y= 0.18)+
  geom_violin(fill= NA)


bdT

#by Species
bdSp <- data.frame(distToMed= dispg$distances, Species= dispg$group)


bdS <- ggline(bdSp, x = "Species", y = "distToMed", plot_type= "p",
              add = c("mean_se"), error.plot = "pointrange",
              order = c("P", "L", "C", "S"),
              ylim= c(0.0,0.3),
              #ylab = "Floral trait variation (Gower distance)", 
              ylab= "",
              xlab = "Species",
              #font.y = list(size = 14, face= "plain"),
              font.x = list(size = 14, face= "plain"))+
  stat_compare_means(label.y= 0.29)+
  stat_compare_means(label= "p.signif", ref.group = ".all.")+#, label.y= 0.16)+
  geom_hline(yintercept = mean(bdSp$distToMed), linetype = 2)+
  geom_violin(fill= NA)

bdS


#plot betadisper values by region
bdReg <- data.frame(distToMed= dispgr$distances, Region= dispgr$group)

bdR <- ggline(bdReg, x = "Region", y = "distToMed", plot_type= "p",
              add = c("mean_se"), error.plot = "pointrange",
              order = c("P", "L", "C","sS", "wS", "nS"),
              ylim= c(0.0,0.3),
              # ylab = "Floral trait variation (Gower distance)"
              ylab= "",
              xlab = "Region",
              #font.y = list(size = 14, face= "plain"),
              font.x = list(size = 14, face= "plain"))+
  stat_compare_means(label.y= 0.29)+
  stat_compare_means(label= "p.signif", ref.group = ".all.")+#, label.y= 0.16)+
  geom_hline(yintercept = mean(bdReg$distToMed), linetype = 2)+
  geom_violin(fill= NA)

bdR

#Figure 3 D-F: floral trait distance by taxa/species/region
flrVar <- ggarrange(bdT, bdS, bdR, ncol= 3, nrow= 1)
flrVar

#plot WITHIN POP COLOR DISTANCE BY REGION
#copied from above for more polished figure
bdPop <- data.frame(distToMedPop= dispgColPop$distances, Pop= dispgColPop$group, Reg2= ca$reg2, Sp= ca$sp)


bdPC <- ggline(bdPop, x = "Reg2", y = "distToMedPop", plot_type= "p",
               add = c("mean_se"), error.plot = "pointrange",
               order = c("P", "L", "C","sS", "wS", "nS"),
               #ylim= c(0.05,0.19), #size = 1,
               ylab = "Within-population color distance (polymorphism)", xlab = "Region",
               font.x = list(size = 14, face= "plain"))+
  stat_compare_means(label.y= 0.51)+
  stat_compare_means(label= "p.signif", ref.group = ".all.", label.y= 0.47)+
  geom_hline(yintercept = mean(bdPop$distToMedPop), linetype = 2)+
  geom_violin(fill= NA)

#Appendix S8- within pop color distance by region
bdPC

#Figure 3A-C: Evolutionary Distance/Sequence Divergence by taxa/sp/region

#from EvolutionaryDistance_rev2021-02-25.R

#BASED ON FOURFOLD DEGENERATE SITES W/IN GROUPS: 
#Estimates of average Codon-based Evolutionary Divergence over Sequence Pairs within Groups
#Calculated in MEGA; MEGA output read in as csv for comparisons among taxa/species/regions
#(variance calculation via bootstrapping not available for this method in MEGA)

## plot MEGA output for Castilleja dataset as figures
## within group Evolutionary Divergence (#base pair substitutions per 100 sites)

#Taxa level
evoTaxa <- read.csv("Cast_WithinEvoDiv_Taxa_rev2021-02-25.csv")
evoTaxa <- na.omit(evoTaxa)

#Species level
evoSp <- read.csv("Cast_WithinEvoDiv_Species_rev2021-02-25.csv")

#Region level
evoReg <- read.csv("Cast_WithinEvoDiv_Region_rev2021-02-25.csv")

#taxa
evoT <-   ggline(evoTaxa, x = "Taxa", y = "AvgCodon_based_EvoDiv_x100", plot_type= "p", size= 2,
                 ylim= c(2.9,3.5))+
  ylab("Genetic distance within groups (fourfold-degenerate codon positions)")
#geom_errorbar(aes(ymin=EvoDiv - SE, ymax= EvoDiv + SE, width= 0))

evoT


#species
evoS <- ggline(evoSp, x = "Species", y = "AvgCodon_based_EvoDiv_x100", plot_type= "p", 
               size= 2, order = c("P", "L", "C","S"), ylim= c(2.9,3.5))+
  ylab("")
#ylab("Genetic distance within groups (fourfold-degenerate codon positions)")
#geom_errorbar(aes(ymin=EvoDiv - SE, ymax= EvoDiv + SE, width= 0))

evoS


#region
evoR <-   ggline(evoReg, x = "Region", y = "AvgCodon_based_EvoDiv_x100", plot_type= "p", 
                 size= 2, order = c("P", "L", "C","sS", "wS", "nS"), ylim= c(2.9,3.5))+
  # ylab("Genetic distance within groups (fourfold-degenerate codon positions)")
  ylab("")
#geom_errorbar(aes(ymin=EvoDiv - SE, ymax= EvoDiv + SE, width= 0))

evoR


genVar <- ggarrange(evoT, evoS, evoR, ncol=3, nrow= 1) 
genVar

#FIGURE 3
fig3 <- ggarrange(genVar, flrVar, ncol=1, nrow=2)
fig3

#FIGURE 4- PAIRWISE COMPARISONS 

#calculate population mean values for pairwise-population comparisons
#read in floral trait data
cV <- read.csv("CastFloral_v2-coreTraits_corrected9-7-20.csv")
#columns with useful traits: 6-10 (morph); 11-13 (color)

#change taxa and sp to number code
#taxa: S=1, PC=2 
cV$taxa.n <- 1 #S
cV[cV$taxa == "PC", "taxa.n"] <- 2 #PC

#species: S=1, P=2, C=3, L=4 
cV$sp.n <- 1 #S
cV[cV$sp == "P", "sp.n"] <- 2 #P
cV[cV$sp == "C", "sp.n"] <- 3 #C
cV[cV$sp == "L", "sp.n"] <- 4 #L


#create dataframe with population mean traits
popTr <- group_by(cV, pop) %>%
  summarise(taxa = mean(taxa.n), sp = mean(sp.n), corL = mean(corL), corW2 = mean(corW2),
            stigma = mean(stigma), lipL = mean(lipL), 
            brLbW = mean(brLbW), red = mean(red), green = mean(green),
            blue = mean(blue), LatN = mean(LatN), LonW = mean(LonW))

#subset by taxa
spt <- popTr[12:26, ]
pcpt <- popTr[1:11, ]

#create dataframe with population mean traits- color only
colpt <- group_by(cV, pop) %>%
  summarise(taxa = mean(taxa.n), sp = mean(sp.n), red = mean(red), green = mean(green),
            blue = mean(blue), LatN = mean(LatN), LonW = mean(LonW))

#create dataframe with population mean traits- morph only
morpt <- group_by(cV, pop) %>%
  summarise(taxa = mean(taxa.n), sp = mean(sp.n), corL = mean(corL), corW2 = mean(corW2),
            stigma = mean(stigma), lipL = mean(lipL), 
            brLbW = mean(brLbW), LatN = mean(LatN), LonW = mean(LonW))

#all groups, just color first

#also calculate pop-pairwise geographic distances
#using package geodist, method "geodesic"
cdist <- as.matrix(t(dist(colpt[ ,4:6])))
gdist <- geodist(colpt[,7:8], measure= "geodesic")

#subset to pull out within spp from between spp
Scolpt <- colpt[12:26,]
PXcolpt <- colpt[1:11,]

Ccolpt <- colpt[1:4,]
Lcolpt <- colpt[5:7,]
Pcolpt <- colpt[8:11,]

#color distances among C. sessiliflora populations:
#within S color
scdist <- as.matrix(t(dist(Scolpt[ ,4:6])))
#geographic distances among C. sessiliflora populations:
sgdist <- geodist(Scolpt[,7:8], measure= "geodesic")

#export color distance as csv
sColDist <- as.data.frame(scdist)
sColDist
write.csv(sColDist, file = "CASE_pw_colDist.csv")

#export geodesic distance as csv
sGeoDist <- as.data.frame(sgdist)
sGeoDist
write.csv(sGeoDist, file = "CASE_pw_geoDist.csv")

#color distances among C. purpurea complex populations:
#CP complex color
pxcdist <- as.matrix(t(dist(PXcolpt[ ,4:6])))
#geographic distances among C. purpurea complex populations:
pxgdist <- geodist(PXcolpt[,7:8], measure= "geodesic")
#units= meters

#export CP color distance as csv
pColDist <- as.data.frame(pxcdist)
pColDist
write.csv(pColDist, file = "CAPU_pw_colDist.csv")

#export geodesic distance as csv
pxGeoDist <- as.data.frame(pxgdist)
pxGeoDist
write.csv(pxGeoDist, file = "CAPU_pw_geoDist.csv")

#calculate pairwise genetic distance measures
#v2: genetic distance based on 4fold sites, calculated in MEGA and compiled with above color/geo distance data
#v3: genetic distance based on the PCA of SNPs, calculated below:

#v3- use PCA distances to calculate pairwise genetic distance
pcaY <- read.csv("CaInd_hybseq_SNP-PCA_eigen4dist_2021-02-18.csv")

#need to average eigenvalues by population (only pops with floral data)- so for pops w/ 2 individuals, will take mean of each PC score
pcaZ <- subset(pcaY, floral== "y")

pcaV <- group_by(pcaZ, pop, spp, gen, morph, lat, long, ) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), PC3 = mean(PC3), PC4 = mean(PC4), PC5 = mean(PC5), PC6 = mean(PC6), PC7 = mean(PC7),
            PC8 = mean(PC8), PC9 = mean(PC9), PC10 = mean(PC10), PC11 = mean(PC11), PC12 = mean(PC12), PC13 = mean(PC13), PC14 = mean(PC14),
            PC15 = mean(PC15), PC16 = mean(PC16), PC17 = mean(PC17), PC18 = mean(PC18), PC19 = mean(PC19), PC20 = mean(PC20))

#visualize population-averaged PCA (mostly for fun)
x <- ggplot(pcaV, aes(PC1, PC2, col = morph, shape = gen, label= pop)) + 
  geom_point(size = 3)+
  geom_text(aes(label= pop))
x

#use just PC1 and PC2 (which explain the most variation)
pcaD <- as.matrix(t(dist(pcaV[,7:8], method= "euclidean")))
#export as csv to compile as columns with floral color dist and geo dist
write.csv(pcaD, "Cast_SNP-PCA_popAvg_distMatrix_2021-02-18.csv")

#with exported csv's, compile data for each pairwise comparison among populations
#of 1) C. sessiliflora and 2) the C. purpurea sp complex
#append with pairwise estimates of sequence divergence calculcated from fourfold degenerate sites (=GenDist_v2)
#and with euclidean distances calculated from PC1 and PC2 from PCA of SNPs by populations (GenDist_v3)


#read in pairwise matrices compiled in dataframe
#PERFORM ANALYSES WITH THIS COMPILED DATAFRAME: 
ppw <- read.csv("CAPU_PW_gen4f_col_geo_SNPpca1-2_v3_2021-02-18.csv")
spw <- read.csv("CASE_PW_gen4f_col_geo_SNPpca1-2_v3_2021-02-18.csv")


#GENDISTV2 CALCULATED USING ONLY 4-FOLD DEGENERATE SITES
#gen x geo MORPH
sggm <- ggplot(spw, aes(x=GeodesicDist_km, y=GenDistv2, color=morph, linetype= morph)) +
  geom_point(size =2) + 
  #scale_x_continuous(labels = comma) +
  geom_smooth(method="lm", alpha= 0.2, fullrange=TRUE)+
  #scale_shape_manual(values=c(3, 16, 17))+ 
  scale_color_manual(values=c('salmon','yellowgreen', 'gold'))+
  scale_linetype_manual(values=c("pink"= "dashed", "yellow"= "solid", "typical"="dashed"))+
  labs(title= "4-fold degenerate sites", y= "Genetic distance", x="Geographic distance (km)")+
  theme_classic()

#Fig 4A
sggm


#pink morph
p_sgg <- lm(formula = GenDistv2 ~ GeodesicDist_km, data = subset(spw, morph== "pink"))
summary(p_sgg)
# Multiple R-squared:  0.02309,	Adjusted R-squared:  0.0008833 
# F-statistic:  1.04 on 1 and 44 DF,  p-value: 0.3134

#typical morph
t_sgg <- lm(formula = GenDistv2 ~ GeodesicDist_km, data = subset(spw, morph== "typical"))
summary(t_sgg)
# Multiple R-squared:  0.004697,	Adjusted R-squared:  -0.01845 
# F-statistic: 0.2029 on 1 and 43 DF,  p-value: 0.6546

#yellow morph
y_sgg <- lm(formula = GenDistv2 ~ GeodesicDist_km, data = subset(spw, morph== "yellow"))
summary(y_sgg)
#p=0.02
# Multiple R-squared:  0.3786,	Adjusted R-squared:  0.3268 
# F-statistic: 7.311 on 1 and 12 DF,  p-value: 0.01917

#GenDist_v3 calculated using PCA of SNPs
#gen x geo MORPH
sgg <- ggplot(spw, aes(x=GeodesicDist_km, y=GenDist_v3_SNP_pc1_2, color=morph, linetype= morph)) +
  geom_point(size =2) + 
  #scale_x_continuous(labels = comma) +
  geom_smooth(method="lm", alpha= 0.2, fullrange=TRUE)+
  #scale_shape_manual(values=c(3, 16, 17))+ 
  scale_color_manual(values=c('salmon','yellowgreen', 'gold'))+
  scale_linetype_manual(values=c("pink"= "solid", "yellow"= "solid", "typical"="solid"))+ #all slopes significant
  labs(title= "SNP PC1 and 2", y= "Genetic distance", x="Geographic distance (km)")+
  theme_classic()

#Fig 4B
sgg



#pink morph
p_sgg <- lm(formula = GenDist_v3_SNP_pc1_2 ~ GeodesicDist_km, data = subset(spw, morph== "pink"))
summary(p_sgg)
# Multiple R-squared:  0.6389,	Adjusted R-squared:  0.6307 
# F-statistic: 77.85 on 1 and 44 DF,  p-value: 2.739e-11

#typical morph
t_sgg <- lm(formula = GenDist_v3_SNP_pc1_2 ~ GeodesicDist_km, data = subset(spw, morph== "typical"))
summary(t_sgg)
# Multiple R-squared:  0.5291,	Adjusted R-squared:  0.5181 
# F-statistic: 48.31 on 1 and 43 DF,  p-value: 1.517e-08

#yellow morph
y_sgg <- lm(formula = GenDist_v3_SNP_pc1_2 ~ GeodesicDist_km, data = subset(spw, morph== "yellow"))
summary(y_sgg)
# Multiple R-squared:  0.5219,	Adjusted R-squared:  0.482 
# F-statistic:  13.1 on 1 and 12 DF,  p-value: 0.00352

#FOURFOLD SITES
#gen x geo 
pggm <- ggplot(ppw, aes(x=GeodesicDist_km, y=GenDistv2, color=w_bSp, linetype=w_bSp)) +
  geom_point(size =2, shape =16) + 
  ylim(0.025, 0.045)+
  xlim(0,700)+
  #scale_x_continuous(labels = comma) +
  geom_smooth(method="lm", alpha=0.2, fullrange=FALSE)+
  scale_color_manual(values=c('sienna','aquamarine3'))+
  scale_linetype_manual(values=c("within"= "dashed", "between"= "solid"))+
  labs(title= "4-fold degenerate sites", y= "Genetic distance", x="Geographic distance (km)")+
  theme_classic()

#Fig 4C
pggm


#gen x geo
#within species lm
w_pgg <- lm(formula = GenDistv2 ~ GeodesicDist_km, data = subset(ppw, w_bSp== "within"))
summary(w_pgg)
# Multiple R-squared:  0.2178,	Adjusted R-squared:  0.1577 
# F-statistic:  3.62 on 1 and 13 DF,  p-value: 0.07946

#between species
b_pgg <- lm(formula = GenDistv2 ~ GeodesicDist_km, data = subset(ppw, w_bSp== "between"))
summary(b_pgg)
# Multiple R-squared:  0.1162,	Adjusted R-squared:  0.09292 
# F-statistic: 4.995 on 1 and 38 DF,  p-value: 0.03138



#PCA OF SNPS
#gen x geo COLOR #both slopes sig=solid
pgg <- ggplot(ppw, aes(x=GeodesicDist_km, y=GenDist_v3_SNP_pc1_2, color=w_bSp)) +
  geom_point(size =2, shape =16) + 
  ylim(0,0.6)+
  xlim(0,700)+
  #scale_x_continuous(labels = comma) +
  geom_smooth(method="lm", alpha=0.2, fullrange=FALSE)+
  scale_color_manual(values=c('sienna','aquamarine3'))+
  labs(title= "SNP PC1 and 2", y= "Genetic distance", x="Geographic distance (km)")+
  theme_classic()

#Fig 4D
pgg

#gen x geo
#within species lm
w_pgg <- lm(formula = GenDist_v3_SNP_pc1_2 ~ GeodesicDist_km, data = subset(ppw, w_bSp== "within"))
summary(w_pgg)
# Multiple R-squared:  0.3504,	Adjusted R-squared:  0.3004 
# F-statistic: 7.011 on 1 and 13 DF,  p-value: 0.0201

#between species
b_pgg <- lm(formula = GenDist_v3_SNP_pc1_2 ~ GeodesicDist_km, data = subset(ppw, w_bSp== "between"))
summary(b_pgg)
# Multiple R-squared:  0.0969,	Adjusted R-squared:  0.07314 
# F-statistic: 4.077 on 1 and 38 DF,  p-value: 0.05055



#1- FOURFOLD SITES
#color x gen MORPH
scgmC <- ggplot(spw, aes(y=GenDistv2, x=ColorDist, color=morph, linetype= morph)) +
  geom_point(size =2) + 
  geom_smooth(method="lm", alpha=0.2, fullrange=FALSE)+
  #scale_shape_manual(values=c(3, 16, 17))+ 
  scale_color_manual(values=c('salmon','yellowgreen', 'gold'))+
  scale_linetype_manual(values=c("pink"= "dashed", "yellow"= "dashed", "typical"="dashed"))+
  labs(title= "4-fold degenerate sites", y= "Genetic distance", x="Color distance")+
  theme_classic()

#Fig 4E
scgmC

#pink morph
p_scg <- lm(formula = GenDistv2 ~ ColorDist, data = subset(spw, morph== "pink"))
summary(p_scg)
# Multiple R-squared:  0.004094,	Adjusted R-squared:  -0.01854 
# F-statistic: 0.1809 on 1 and 44 DF,  p-value: 0.6727

#yellow morph
y_scg <- lm(formula = GenDistv2 ~ ColorDist, data = subset(spw, morph== "yellow"))
summary(y_scg)
# Multiple R-squared:  0.05886,	Adjusted R-squared:  -0.01957 
# F-statistic: 0.7505 on 1 and 12 DF,  p-value: 0.4033

#typical
t_scg <- lm(formula = GenDistv2 ~ ColorDist, data = subset(spw, morph== "typical"))
summary(t_scg)
# Multiple R-squared:  0.02588,	Adjusted R-squared:  0.003225 
# F-statistic: 1.142 on 1 and 43 DF,  p-value: 0.2911

#2-PCA OF SNPS
#color x gen MORPH
scgC <- ggplot(spw, aes(y=GenDist_v3_SNP_pc1_2, x=ColorDist, color=morph, linetype= morph)) +
  geom_point(size =2) + 
  geom_smooth(method="lm", alpha=0.2, fullrange=FALSE)+
  #scale_shape_manual(values=c(3, 16, 17))+ 
  scale_color_manual(values=c('salmon','yellowgreen', 'gold'))+
  scale_linetype_manual(values=c("pink"= "solid", "yellow"= "solid", "typical"="dashed"))+
  labs(title= "SNP PC1 and 2", y= "Genetic distance", x="Color distance")+
  theme_classic()

#fig 4F
scgC

#pink morph
p_scg <- lm(formula = GenDist_v3_SNP_pc1_2 ~ ColorDist, data = subset(spw, morph== "pink"))
summary(p_scg)
# Multiple R-squared:  0.2838,	Adjusted R-squared:  0.2675 
# F-statistic: 17.43 on 1 and 44 DF,  p-value: 0.0001383

#yellow morph
y_scg <- lm(formula = GenDist_v3_SNP_pc1_2 ~ ColorDist, data = subset(spw, morph== "yellow"))
summary(y_scg)
# Multiple R-squared:  0.3457,	Adjusted R-squared:  0.2912 
# F-statistic:  6.34 on 1 and 12 DF,  p-value: 0.02701

#typical
t_scg <- lm(formula = GenDist_v3_SNP_pc1_2 ~ ColorDist, data = subset(spw, morph== "typical"))
summary(t_scg)
# Multiple R-squared:  0.02393,	Adjusted R-squared:  0.001231 
# F-statistic: 1.054 on 1 and 43 DF,  p-value: 0.3103


#3- FOURFOLD SITES
#color x gen COLOR
pcgmC <- ggplot(ppw, aes(y=GenDistv2, x=ColorDist, color=w_bSp, linetype=w_bSp)) +
  geom_point(size =2, shape =16) + 
  geom_smooth(method="lm", alpha=0.2, fullrange=FALSE)+
  ylim(0.025, 0.045)+
  scale_color_manual(values=c('sienna','aquamarine3'))+
  scale_linetype_manual(values=c("within"= "dashed", "between"= "dashed"))+
  labs(title= "4-fold degenerate sites", y= "Genetic distance", x="Color distance")+
  theme_classic()

#Fig 4G
pcgmC


#within species lm
w_pgg <- lm(formula = GenDistv2 ~ ColorDist, data = subset(ppw, w_bSp== "within"))
summary(w_pgg)
#ns, p=.96
# Multiple R-squared:  0.0002286,	Adjusted R-squared:  -0.07668 
# F-statistic: 0.002972 on 1 and 13 DF,  p-value: 0.9574


#between species
b_pgg <- lm(formula = GenDistv2 ~ ColorDist, data = subset(ppw, w_bSp== "between"))
summary(b_pgg)
#ns, p=.26
# Multiple R-squared:  0.03238,	Adjusted R-squared:  0.006914 
# F-statistic: 1.272 on 1 and 38 DF,  p-value: 0.2666

#4- PCA OF SNPS
#color x gen COLOR
pcgC <- ggplot(ppw, aes(y=GenDist_v3_SNP_pc1_2, x=ColorDist, color=w_bSp, linetype= w_bSp)) +
  geom_point(size =2, shape =16) + 
  ylim(0,0.6)+
  geom_smooth(method="lm", alpha=0.2, fullrange=FALSE)+
  scale_color_manual(values=c('sienna','aquamarine3'))+
  scale_linetype_manual(values = c("between"= "solid","within"= "dashed"))+
  #ylim(0,0.6)+
  labs(title= "SNP PC1 and 2", y= "Genetic distance", x="Color distance")+
  theme_classic()

#Fig 4H
pcgC

#color x gen

#within species lm
w_pgg <- lm(formula = ColorDist ~ GenDist_v3_SNP_pc1_2, data = subset(ppw, w_bSp== "within"))
summary(w_pgg)
# Multiple R-squared:  0.178,	Adjusted R-squared:  0.1148 
# F-statistic: 2.815 on 1 and 13 DF,  p-value: 0.1172

#between species
b_pgg <- lm(formula = ColorDist ~ GenDist_v3_SNP_pc1_2, data = subset(ppw, w_bSp== "between"))
summary(b_pgg)
# Multiple R-squared:  0.5967,	Adjusted R-squared:  0.5861 
# F-statistic: 56.23 on 1 and 38 DF,  p-value: 5.256e-09


#plot by taxa, both MEGA v2 and SNP pca v3 gen distances- ALT GEN X COLOR
sv3 <- ggarrange(sggm, sgg, scgmC, scgC, ncol= 2, nrow=2, common.legend = TRUE)
sv3

pv3 <- ggarrange(pggm, pgg, pcgmC, pcgC, ncol= 2, nrow=2, common.legend = TRUE)
pv3

#compiled for figure
#Figure 4
fig4 <- ggarrange(sv3, pv3, ncol=2, nrow=1)
fig4



#Appendix S1C- color distance by geographic distance:

#C purpurea complex

#color x geo COLOR
pcgeo <- ggplot(ppw, aes(x=GeodesicDist_km, y=ColorDist, color=w_bSp)) +
  geom_point(size =2, shape =16) + 
  #scale_x_continuous(labels = comma) +
  xlim(0,700)+
  geom_smooth(method="lm", alpha=.2, fullrange=TRUE, linetype= "dashed")+ #both slopes nonsig
  labs(x= "Geographic distance (km)", y="Color distance")+
  scale_color_manual(values=c('sienna','aquamarine3'))+
  theme_classic()

#Appendix S1C
pcgeo

#color x geo
#within species lm
w_pcg <- lm(formula = ColorDist ~ GeodesicDist_km, data = subset(ppw, w_bSp== "within"))
summary(w_pcg)
# Multiple R-squared:  0.03978,	Adjusted R-squared:  -0.03409 
# F-statistic: 0.5385 on 1 and 13 DF,  p-value: 0.4761

#between species
b_pcg <- lm(formula = ColorDist ~ GeodesicDist_km, data = subset(ppw, w_bSp== "between"))
summary(b_pcg)
# Multiple R-squared:  0.001111,	Adjusted R-squared:  -0.02518 
# F-statistic: 0.04226 on 1 and 38 DF,  p-value: 0.8382



#C. sessiliflora: color x geographic dist by MORPH
scgeo <- ggplot(spw, aes(x=GeodesicDist_km, y=ColorDist, color=morph)) +
  geom_point(size =2) + 
  #scale_x_continuous(labels = comma) +
  geom_smooth(method="lm", alpha=0.2, fullrange=TRUE)+
  #scale_shape_manual(values=c(3, 16, 17))+ 
  scale_color_manual(values=c('salmon','yellowgreen', 'gold'))+
  theme_classic()

scgeo


#pink morph
p_scg <- lm(formula = ColorDist ~ GeodesicDist, data = subset(spw, morph== "pink"))
summary(p_scg)
#typical morph
t_sgg <- lm(formula = ColorDist ~ GeodesicDist, data = subset(spw, morph== "typical"))
summary(t_sgg)
#yellow morph
y_sgg <- lm(formula = ColorDist ~ GeodesicDist, data = subset(spw, morph== "yellow"))
summary(y_sgg)
#all significant, p< 0.05
