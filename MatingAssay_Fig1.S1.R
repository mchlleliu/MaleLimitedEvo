###################################
#
#                             Grieshop et al. 2023
#                      Author: Ryan S. Frost, (Michelle Liu)
#             DsRed experimental evolution - Transcriptomics analysis
#                          Fitness Assay GLMM & Plots
# 
# 
###################################

rm(list=ls())
setwd("~/Desktop/UofT/SSAV_RNA/")

#Load in fitness data
DsRed_Data <- read.csv("Fitness_Assay/DsRed_Data_Formatted.csv")

# Data Prep
#######
# Set up a vector full of zeroes of the appropriate length, then...
Prop_Mated <- rep.int(0, dim(DsRed_Data)[1])

# iterate over each observation in the DsRed_Data object
for (i in 1:dim(DsRed_Data)[1]){
  # if the observation is Male for that population replicate cage and generation
  if (DsRed_Data$Sex[i] == "Male"){
    # the proportion of mated males 
    Prop_Mated[i] <- DsRed_Data$Male_Prop_Mated[i]
  }
  else{
    # the proportion of mated females for that population replicate cage and generation
    Prop_Mated[i] <- DsRed_Data$Female_Prop_Mated[i]
  }
  #... fill in said vector at any given point with the appropriate value 
  # for the proportion of sampled individuals of the given sex found to be mated...
}

#...and then add that vector as another column.
DsRed_Data <- cbind(DsRed_Data, Prop_Mated) 


# We don't need these columns anymore from the .csv file
DsRed_Data$Male_Prop_Mated <- NULL
DsRed_Data$Female_Prop_Mated <- NULL
#######

# Split the data to different subsets
######
# by sex
DsRed_Female <- DsRed_Data[DsRed_Data$Sex == "Female", ]
DsRed_Male <- DsRed_Data[DsRed_Data$Sex == "Male", ]

# Red and NonRed females
DsRed_F_NR <- DsRed_Female[DsRed_Female$Marker == "Non-Red", ] 
DsRed_F_R <- DsRed_Female[DsRed_Female$Marker == "Red", ] 

# Red and NonRed males
DsRed_M_NR <- DsRed_Male[DsRed_Male$Marker == "Non-Red", ] 
DsRed_M_R <- DsRed_Male[DsRed_Male$Marker == "Red", ] 

# NonRed or Red females not mated or mated
DsRed_F_NR_NM <- DsRed_F_NR[DsRed_F_NR$Mating_Status == "Non-Mated", ]
DsRed_F_NR_M <- DsRed_F_NR[DsRed_F_NR$Mating_Status == "Mated", ]
DsRed_F_R_NM <- DsRed_F_R[DsRed_F_R$Mating_Status == "Non-Mated", ]
DsRed_F_R_M <- DsRed_F_R[DsRed_F_R$Mating_Status == "Mated", ]

# NonRed or Red males not mated or mated
DsRed_M_NR_NM <- DsRed_M_NR[DsRed_M_NR$Mating_Status == "Non-Mated", ]
DsRed_M_NR_M <- DsRed_M_NR[DsRed_M_NR$Mating_Status == "Mated", ]
DsRed_M_R_NM <- DsRed_M_R[DsRed_M_R$Mating_Status == "Non-Mated", ]
DsRed_M_R_M <- DsRed_M_R[DsRed_M_R$Mating_Status == "Mated", ]
######


# 1. Is there any effect of time period on proportion of red among maters?
# G-test of independence for effect of Red or NonRed marker
# nominal variables are Time_Block and Marker
######
# Prepare dataset. 
# create a contingency table, with counts of Mated individuals from each time block X marker category
DsRed_Input <- read.csv("Fitness_Assay/DsRed_Data_AsRecorded.csv")

#It is easier to work with the data in this format for the purposes of the G-Test.
#This is the raw data without being collapsed as in DsRed_Data
rownames(DsRed_Input) <- DsRed_Input$Cage.ID
DsRed_Input <- DsRed_Input[, 2:dim(DsRed_Input)[2]]
#Take the first column which has the (Non numeric) cage and pop and generation indicator, make it into the rownames to play nice with colsums.
DsRed_CageBlind <- data.frame(colSums(DsRed_Input[,-1]))
# Exclude the first row which has the (Non numeric) sex and mating status and genotype indicator, 
# play nice with colsums.(Essentially count the number of maters per column)

#get subsets of the data to examine overlap
#Indices which have Male Data
X <- grep("Male", rownames(DsRed_CageBlind), ignore.case = FALSE, value = FALSE) 
#Indices which have Female Data
X2 <- grep("Female", rownames(DsRed_CageBlind), ignore.case = FALSE, value = FALSE)
#B Indices which do not have Non-mated data
Y <- grep("NM", rownames(DsRed_CageBlind), ignore.case = FALSE, value = FALSE, invert = TRUE)
#C indices which have Non-Red Data
Z <- grep("Non.Red", rownames(DsRed_CageBlind), ignore.case = FALSE, value = FALSE)

# For males
# Overlap of X, Y, but NOT Z. Takes the number of maters from Red males across the different time blocks
DsRed_CageBlind_RedMatedMales <- DsRed_CageBlind[setdiff(intersect(X,Y),Z), ]
# Overlap of X, Y, and Z. Takes the number of maters from NonRed males across the different time blocks
DsRed_CageBlind_NonRedMatedMales <- DsRed_CageBlind[intersect(intersect(X,Y),Z), ]

# contingency table for number of Red and NonRed maters across time blocks 
timeblock_maletab <-as.table(rbind(c(DsRed_CageBlind_RedMatedMales), c(DsRed_CageBlind_NonRedMatedMales)))
dimnames(timeblock_maletab) <- list(Marker=c("Red","Non-Red"),
                                    Time_Block=c("Morn D1","AN D1", "E D1", "Morn D2","AN D2", "E D2", "Morn D3"))
#Having carved up the data appropriately above, I make the table which will go into the GTest
library(DescTools)
#Remember to load the library with the GTest function
(Xsq <- GTest(timeblock_maletab))
#Then run the test


# G-test for females
#Overlap of X2, Y, but NOT Z
DsRed_CageBlind_RedMatedFemales <- DsRed_CageBlind[setdiff(intersect(X2,Y),Z), ]
#Overlap of X2, Y, and Z
DsRed_CageBlind_NonRedMatedFemales <- DsRed_CageBlind[intersect(intersect(X2,Y),Z), ]

timeblock_femaletab <-as.table(rbind(c(DsRed_CageBlind_RedMatedFemales), c(DsRed_CageBlind_NonRedMatedFemales)))
dimnames(timeblock_femaletab) <- list(Marker=c("Red","Non-Red"),
                                      Time_Block=c("Morn D1","AN D1", "E D1", "Morn D2","AN D2", "E D2", "Morn D3"))
(Xsq <- GTest(timeblock_femaletab))
#From the empty hashtag to here is repeating what I did for males, but for females.

female.markprev <- DsRed_CageBlind_RedMatedFemales/DsRed_CageBlind_NonRedMatedFemales
male.markprev <- DsRed_CageBlind_RedMatedMales/DsRed_CageBlind_NonRedMatedMales
#Calculating the prevalence of the red-marked genotype among females and males, for each time block.
markerxtimeblock <- data.frame(female.markprev, male.markprev, 1:7)
#Compiling that information into a single data frame for use with ggplot
######


# Is the proportion of red related to the amount of mating activity?
# take G-Test tables (prepared above) and test for correlation between Proportion Red:Non Red to number of individuals sampled and time block. 
#######
#Normalize by number of samples taken over time block, then. 6 for Morn and E, 8 for AN
mating_activity_maletab <-as.table(rbind(c(DsRed_CageBlind_RedMatedMales/DsRed_CageBlind_NonRedMatedMales), c(DsRed_CageBlind_RedMatedMales + DsRed_CageBlind_NonRedMatedMales)))
cor.test(mating_activity_maletab[1,], mating_activity_maletab[2,]/c(6, 8, 6, 6, 8, 6, 4))

mating_activity_femaletab <-as.table(rbind(c(DsRed_CageBlind_RedMatedFemales/DsRed_CageBlind_NonRedMatedFemales), c(DsRed_CageBlind_RedMatedFemales + DsRed_CageBlind_NonRedMatedFemales)))
cor.test(mating_activity_femaletab[1,], mating_activity_femaletab[2,]/c(6, 8, 6, 6, 8, 6, 4))
#Correlation test for proportion of mated individuals which were red with total number of mated individuals in both sexes; do the different genotypes gain a competitive advantage as intensity of mating changes?

#"Cagewise" below is for when we wish to do more nuanced analyses which consider the individual cages, as opposed to only the sum of all cages for a column.
#Indices which have Male Data
Xc <- grep("Male", colnames(DsRed_Input), ignore.case = FALSE, value = FALSE) 
#Indices which have Female Data
X2c <- grep("Female", colnames(DsRed_Input), ignore.case = FALSE, value = FALSE)
#B Indices which do not have Non-mated data
Yc <- grep("NM", colnames(DsRed_Input), ignore.case = FALSE, value = FALSE, invert = TRUE)
#C indices which have Non-Red Data
Zc <- grep("Non.Red", colnames(DsRed_Input), ignore.case = FALSE, value = FALSE)

#Overlap of X, Y, but NOT Z. All columns containing number of maters for Red males
DsRed_Cagewise_RedMatedMales <- DsRed_Input[1:dim(DsRed_Input)[1], setdiff(intersect(Xc,Yc),Zc)]
#Overlap of X, Y, and Z. All columns containing number of maters for NonRed males
DsRed_Cagewise_NonRedMatedMales <- DsRed_Input[1:dim(DsRed_Input)[1], intersect(intersect(Xc,Yc),Zc), ]

#Overlap of X2, Y, but NOT Z. All columns containing number of maters for Red females
DsRed_Cagewise_RedMatedFemales <- DsRed_Input[1:dim(DsRed_Input)[1], setdiff(intersect(X2c,Yc),Zc)]
#Overlap of X2, Y, and Z. All columns containing number of maters for NonRed females
DsRed_Cagewise_NonRedMatedFemales <- DsRed_Input[1:dim(DsRed_Input)[1], intersect(intersect(X2c,Yc),Zc), ]


#Examining if there is any relationship between total mating activity and time.
Mating_by_block <- c()
Mating_by_block_errbar <- c()
correction <- c(6, 8, 6, 6, 8, 6, 4) 
#The correction factors for each time block.
# i.e., during peak mating, cage is observed every 15 mins, but other than that every 30 mins. This correction is to account for that 

# Find the mean number of maters, correcting for the difference in time block length and sampling intensity
for(i in 1:7){
  # This sum is the measure of total mating activity for the current time block
  current_maters <- DsRed_Cagewise_NonRedMatedFemales[, i] + DsRed_Cagewise_NonRedMatedMales[, i] + DsRed_Cagewise_RedMatedFemales[, i] + DsRed_Cagewise_RedMatedMales[, i]
  # Add to the accumulating list
  Mating_by_block <- c(Mating_by_block, mean(current_maters/correction[i]*mean(correction)))
  Mating_by_block_errbar <- c(Mating_by_block_errbar, SD(current_maters/correction[i]*mean(correction))/sqrt(length(current_maters)))
  #Multiplying by mean(correction) to ensure that the correction factors are normalized
}

# the clean data set with results from the test
# add additional column 1:7 for plotting the time blocks as x-axis
mating.activity.TB <- data.frame(Mating_by_block, Mating_by_block_errbar, 1:7)
#######


# Fig. S1: Mating activity by timeblock
Fig1_suppl <- ggplot(data = mating.activity.TB) + 
  geom_errorbar(aes(x=X1.7, ymin = Mating_by_block-Mating_by_block_errbar, 
                    ymax = Mating_by_block+Mating_by_block_errbar),  color = 'black', size=1, width = 0.5) +
  geom_point(aes(x=X1.7, y = Mating_by_block), color = "black", size = 7.5) +
  ylab("Average Number of \nObserved Matings") +
  scale_x_continuous(breaks = 1:7, labels = c("Morning\nDay 13", "Afternoon\nDay 13", "Evening\nDay 13", 
                                              "Morning\nDay 14", "Afternoon\nDay 14", "Evening\nDay 14", 
                                              "Morning\nDay 15")) + 
  geom_vline(xintercept = c(seq(1.5, 6.5)), color = "grey") + 
  
  # some theme settings
  theme_classic() +
  theme(plot.title.position = c("panel"),
        legend.title = element_blank(),
        legend.position = c("none"),
        legend.text = element_text(size=30, margin = margin(10,0,10,0), color = "black", vjust = -0.2),
        axis.text.x = element_text(size=20, margin = margin(5,0,0,0), color = "black"),
        axis.text.y = element_text(size=17.5, margin = margin(0,5,0,0), color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=30, margin = margin(0,10,0,0), color = "black"),
        plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
        plot.margin = margin(6,6,6,6),
        panel.border = element_rect(colour = "black", fill=NA, size=1))


# Save Fig.S1
png(file = "~/Desktop/UofT/SSAV_RNA/Plots/final_2/png_version/Fig_S1.png",   # The directory you want to save the file in
    width = 15, # 20 15 The width of the plot in inches
    height = 8, # 12 8 The height of the plot in inches
    units = "in", res = 300)
Fig1_suppl
dev.off()





#The simpler analysis for proportion of red vs time indicates that there is something worth investigating here; 
#re-do correlation tests and plots on each individual cage within each time block to see how the correlation holds up.
#This result is shown as Figure 1B
#######
# Carving up the data; proportion red maters among total maters (red + nonRed) by (time) block. 
# Instead of Red:NonRed ratio as with the G-test above, calculating the percentage is easier to comprehend for plotting
# each list in the vector is the proportion for a given time block, replicate population, cage, and generation time
## FOR FEMALES
block1 = DsRed_Cagewise_RedMatedFemales[, 1]/(DsRed_Cagewise_RedMatedFemales[, 1] + DsRed_Cagewise_NonRedMatedFemales[, 1])
block2 = DsRed_Cagewise_RedMatedFemales[, 2]/(DsRed_Cagewise_RedMatedFemales[, 2] + DsRed_Cagewise_NonRedMatedFemales[, 2])
block3 = DsRed_Cagewise_RedMatedFemales[, 3]/(DsRed_Cagewise_RedMatedFemales[, 3] + DsRed_Cagewise_NonRedMatedFemales[, 3])
block4 = DsRed_Cagewise_RedMatedFemales[, 4]/(DsRed_Cagewise_RedMatedFemales[, 4] + DsRed_Cagewise_NonRedMatedFemales[, 4])
block5 = DsRed_Cagewise_RedMatedFemales[, 5]/(DsRed_Cagewise_RedMatedFemales[, 5] + DsRed_Cagewise_NonRedMatedFemales[, 5])
block6 = DsRed_Cagewise_RedMatedFemales[, 6]/(DsRed_Cagewise_RedMatedFemales[, 6] + DsRed_Cagewise_NonRedMatedFemales[, 6])
block7 = DsRed_Cagewise_RedMatedFemales[, 7]/(DsRed_Cagewise_RedMatedFemales[, 7] + DsRed_Cagewise_NonRedMatedFemales[, 7])
proplist = c(block1, block2, block3, block4, block5, block6, block7)
blocklist = c()
for( i in 1:7){
  blocklist = c(blocklist, rep_len(i, 72))
}
#There are 72 cages; 6 populations times 4 replicate cages times 3 generations assayed. Hence I repeat the time block indicators 1 through 7, 72 times each.

blocklist <- blocklist[!is.nan(proplist)]
proplist <- proplist[!is.nan(proplist)]
block2 <- block2[!is.nan(block2)]
block3 <- block3[!is.nan(block3)]
block7 <- block7[!is.nan(block7)]
#remove all indices from both lists with a NaN (sometimes there are no Red maters)
#In principle/in general, this should be done for all blocks, it just proves to be necessary in this case only for 2, 3, and 7.

# mean and SE across populations, replicate cage, and generations
femaletimeblock_mean <- c(mean(block1), mean(block2), mean(block3), 
                          mean(block4), mean(block5), mean(block6), mean(block7))
femaletimeblock_Errbars <- c(SD(block1)/sqrt(length(block1)), 
                             SD(block2)/sqrt(length(block2)), 
                             SD(block3)/sqrt(length(block3)), 
                             SD(block4)/sqrt(length(block4)), 
                             SD(block5)/sqrt(length(block5)), 
                             SD(block6)/sqrt(length(block6)), 
                             SD(block7)/sqrt(length(block7)))
#Error bars Standard Error of the Mean (SEM): Standard Deviation / Sqrt(Num Samples)

cor.test(proplist, blocklist)
#The correlation test; also statistically significant when considering the individual data points.
#here individual data points meaning that we did not collapse the populations, replicate cage, and generations.


## -- Doing the same thing for the MALE data as for the female data above
maleblock1 = DsRed_Cagewise_RedMatedMales[, 1]/(DsRed_Cagewise_RedMatedMales[, 1] + DsRed_Cagewise_NonRedMatedMales[, 1])
maleblock2 = DsRed_Cagewise_RedMatedMales[, 2]/(DsRed_Cagewise_RedMatedMales[, 2] + DsRed_Cagewise_NonRedMatedMales[, 2])
maleblock3 = DsRed_Cagewise_RedMatedMales[, 3]/(DsRed_Cagewise_RedMatedMales[, 3] + DsRed_Cagewise_NonRedMatedMales[, 3])
maleblock4 = DsRed_Cagewise_RedMatedMales[, 4]/(DsRed_Cagewise_RedMatedMales[, 4] + DsRed_Cagewise_NonRedMatedMales[, 4])
maleblock5 = DsRed_Cagewise_RedMatedMales[, 5]/(DsRed_Cagewise_RedMatedMales[, 5] + DsRed_Cagewise_NonRedMatedMales[, 5])
maleblock6 = DsRed_Cagewise_RedMatedMales[, 6]/(DsRed_Cagewise_RedMatedMales[, 6] + DsRed_Cagewise_NonRedMatedMales[, 6])
maleblock7 = DsRed_Cagewise_RedMatedMales[, 7]/(DsRed_Cagewise_RedMatedMales[, 7] + DsRed_Cagewise_NonRedMatedMales[, 7])
maleproplist = c(maleblock1, maleblock2, maleblock3, maleblock4, maleblock5, maleblock6, maleblock7)
maleblocklist = c()
for( i in 1:7){
  maleblocklist = c(maleblocklist, rep_len(i, 72))
}

maleblocklist <- maleblocklist[!is.nan(maleproplist)]
maleproplist <- maleproplist[!is.nan(maleproplist)]
maleblock2 <- maleblock2[!is.nan(maleblock2)]
maleblock3 <- maleblock3[!is.nan(maleblock3)]
maleblock7 <- maleblock7[!is.nan(maleblock7)]

maletimeblock_mean <- c(mean(maleblock1), mean(maleblock2), mean(maleblock3), mean(maleblock4), mean(maleblock5), mean(maleblock6), mean(maleblock7))
maletimeblock_Errbars <- c(SD(maleblock1)/sqrt(length(maleblock1)), SD(maleblock2)/sqrt(length(maleblock2)), SD(maleblock3)/sqrt(length(maleblock3)), SD(maleblock4)/sqrt(length(maleblock4)), SD(maleblock5)/sqrt(length(maleblock5)), SD(maleblock6)/sqrt(length(maleblock6)), SD(maleblock7)/sqrt(length(maleblock7)))

cor.test(maleproplist, maleblocklist)
#For males, there is no correlation between the proportion of Red maters and time block.
#######



# 
# GLMM for mating probability
# Model selection was based on the best AIC score
# effects included:
#   Sex * Marker interaction term.
#   Marker * Population interaction term (Population is a random effect). 
#   Generation * Marker interaction term (Generation is a random effect). 
# We can't include time-block in the full model since final time block perfectly co-varies with non-mated status. 
# We already did the analysis of the impact of time block separately (see above)

# Random effects (Population, cage) are always treated as such, 
# same for fixed effects (marker, generation/block- these have too few levels to fit as random). 
# Sex is a fixed effect inherently; we don't draw from a pool of all possible sexes, we actually test all possible sexes in this species.

library(lme4)
library(lmerTest)
library(car)
library(ggplot2)
library(sjPlot)
library(glmmTMB)

# 
# the models below offer the best AIC and seem logically plausible.
markerXpop.bothsex.model <- glmer(Mating_Status == "Mated" ~ Marker*Sex + (Marker|Population), binomial(link='probit'), offset = Prop_Mated, data = DsRed_Data)
summary(markerXpop.bothsex.model) 
extractAIC(markerXpop.bothsex.model)

# Final model for males
#######
#Mating status as fixed effect, population as random effect nested within the marker..?
markerXpop.maleonly.model  <- glmer(Mating_Status == "Mated" ~ Marker + (Marker|Population), binomial(link='probit'), offset = Prop_Mated, data = DsRed_Male)
summary(markerXpop.maleonly.model) # marker is a significant predictor 

extractAIC(markerXpop.maleonly.model) # separating the sexes greatly improved the model
#WARNING, DO EVERYTHING ELSE YOU NEED TO DO WITH THE DATA PRIOR TO RUNNING CONFINT.MERMOD
#IT SEEMS CALCULATING THOSE ALTERS THE BASE MODEL IN SOME WAY OR SOMETHING, CHANGING ALL THE PREVIOUS RESULTS

#SJPlot "re" with transform = NULL allows us to obtain the coef values from the main GLMM subtracting the main/global effect of marker. 
demo <- plot_model(markerXpop.maleonly.model, type = "re", transform = NULL)
demo$data$estimate[7:12] # take only the estimates for the Red marker, separated by population
fixef(markerXpop.maleonly.model)[2] # fixed effects, global, not separated by population
coef(markerXpop.maleonly.model)$Population$MarkerRed

# get these coef values to plot for Fig. 1
glob.fx <- fixef(markerXpop.maleonly.model)[2] # global fixed effects
mean.pop.fx <- demo$data$estimate[7:12] + glob.fx # fixed effects by population
lower.CI.pop.fx <- demo$data$conf.low[7:12] + glob.fx
upper.CI.pop.fx <- demo$data$conf.high[7:12] + glob.fx

preConfint.Mermod.model <- markerXpop.maleonly.model # since running confint.mermod may change the model, i save the old model here
# obtain confidence interval for the model coefficients
confints <- confint.merMod(markerXpop.maleonly.model, oldNames = FALSE)
main_eff_int <- c(confints[5], confints[10]) # confidence interval for the global effect

main_result <- data.frame(mean.pop.fx, lower.CI.pop.fx, upper.CI.pop.fx, 1:6)
#######

# Final model for only females
######
markerXpop.femaleonly.model  <- glmer(Mating_Status == "Mated" ~ Marker + (Marker|Population), binomial(link='probit'), offset = Prop_Mated, data = DsRed_Female)
summary(markerXpop.femaleonly.model) # the marker is not a significant predictor here

female.demo <- plot_model(markerXpop.femaleonly.model, type = "re", transform = NULL)

female.glob.fx <- fixef(markerXpop.femaleonly.model)[2]
mean.female.pop.fx <- female.demo$data$estimate[7:12] + female.glob.fx
female.lower.CI.pop.fx <- female.demo$data$conf.low[7:12] + female.glob.fx
female.upper.CI.pop.fx <- female.demo$data$conf.high[7:12] + female.glob.fx

# again saving the previous model in case confint.merMod changes it
preConfint.Mermod.modelFem <- markerXpop.femaleonly.model
female.confints <- confint.merMod(markerXpop.femaleonly.model, oldNames = FALSE)
female.main_eff_int <- c(female.confints[5], female.confints[10])

female.main_result <- data.frame(mean.female.pop.fx, female.lower.CI.pop.fx, female.upper.CI.pop.fx, 1:6)
######





# -------------- Plotting code Fig. 1 -----------------
library(cowplot)
library(ggpubr)


## Fig. 1A: Effect of Red chromosome on mating success
########
# We want the plots above for the models of markerXmating success, but only using the global fixed effects, without separating the coefficients by population replicate
# get only the global fixed effects and store them in an object
fem.main_eff <- data.frame(list(female.glob.fx, female.main_eff_int[1], female.main_eff_int[2], "b.Female"), row.names = NULL)
colnames(fem.main_eff) <- c("glob.fx", "q0025", "q0975", "Sex")
male.main_eff <- data.frame(list(glob.fx, main_eff_int[1], main_eff_int[2], "a.Male"), row.names = NULL)
colnames(male.main_eff) <- c("glob.fx", "q0025", "q0975", "Sex")
all.main_eff <- rbind(fem.main_eff, male.main_eff)

# plot the global effects
plot.globfx <- ggplot(all.main_eff) + 
  # make horizontal line at 0
  geom_hline(yintercept = 0, size = 1) +
  # plot error bars for coefficients
  geom_errorbar(aes(Sex, glob.fx, ymin=q0025, ymax=q0975, color=Sex), size = 2, width = 0.5) +
  # plot coefs
  geom_point(aes(x=Sex, y=glob.fx, color=Sex), size = 10) +
  # axes settings
  scale_x_discrete(labels = c("Male","Female")) +
  scale_colour_manual(values = c("#0072B2", "#D55E00")) + # "Chr-2", "Chr-3", "X-Chr"
  scale_y_continuous(limits = c(-0.25,0.25)) +
  labs(y = "Change in Mating Probability \nAssociated with Red Chromosome") +
  
  # other plot theme settings
  theme_classic() +
  theme(plot.title.position = c("panel"),
        legend.title = element_blank(),
        legend.position = c("none"),
        legend.text = element_text(size=30, margin = margin(10,0,10,0), color = "black", vjust = -0.2),
        axis.text.x = element_text(size=30, margin = margin(5,0,0,0), color = "black"),
        axis.text.y = element_text(size=15, margin = margin(0,5,0,0), color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=30, margin = margin(0,10,0,0), color = "black"),
        plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
        plot.margin = margin(6,6,6,6),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
########


## Fig. 1B: Proportion of red maters by time block
#######
# make plotting data set of point estimates for figure 1
markerxtimeblock <- data.frame(femaletimeblock_mean, femaletimeblock_Errbars, 
                               maletimeblock_mean, maletimeblock_Errbars, 1:7)

# plot female data
fem_time.block <- ggplot(data = markerxtimeblock) + 
  # dashed regression line 
  geom_abline(aes(slope = lm(proplist~blocklist)$coefficients[2], 
                  intercept = lm(proplist~blocklist)$coefficients[1]),
              size = 1.5, alpha = 0.75, show.legend = F, color = "#D55E00") + 
  # plot errorbar for point estimates
  geom_errorbar(aes(x=X1.7, ymin = femaletimeblock_mean-femaletimeblock_Errbars, 
                    ymax = femaletimeblock_mean+femaletimeblock_Errbars), color = "#D55E00", size = 1.5, width = 0.5) +
  #plot point estimates
  geom_point(aes(x =X1.7,  y=femaletimeblock_mean), color = "#D55E00", size = 7.5) + 
  # set axes
  xlab("Time Block") + ylab("Proportion of Red \nFemales among Maters") + 
  ylim(0.44, 0.61) + 
  geom_vline(xintercept = c(seq(1.5, 6.5)), color = "grey") + 
  scale_x_continuous(breaks = 1:7, 
                     labels = c("Morning\nDay 13", "Afternoon\nDay 13", "Evening\nDay 13", 
                                "Morning\nDay 14", "Afternoon\nDay 14", "Evening\nDay 14", 
                                "Morning\nDay 15")) + 
  
  # other plot theme settings
  theme_classic() +
  theme(plot.title.position = c("panel"),
        legend.title = element_blank(),
        legend.position = c("none"),
        # legend.box.background = element_rect(),
        legend.text = element_text(size=30, margin = margin(10,0,10,0), color = "black", vjust = -0.2),
        axis.text.x = element_text(size=17, margin = margin(5,0,0,0), color = "black"),
        axis.text.y = element_text(size=15, margin = margin(0,5,0,0), color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=24, margin = margin(0,10,0,10), color = "black"),
        plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
        plot.margin = margin(6,6,6,6),
        panel.border = element_rect(colour = "black", fill=NA, size=1))



# plot male data using the same settings
male_time.block <- ggplot(data = markerxtimeblock) + 
  geom_abline(aes(slope = lm(maleproplist~maleblocklist)$coefficients[2], 
                  intercept = lm(maleproplist~maleblocklist)$coefficients[1]),
              size = 1.5, alpha = 0.75, show.legend = F, color = "#0072B2") + 
  geom_errorbar(aes(x=X1.7, ymin = maletimeblock_mean-maletimeblock_Errbars, 
                    ymax = maletimeblock_mean+maletimeblock_Errbars), color = "#0072B2", size = 1.5, width = 0.5) +
  geom_point(aes(x =X1.7,  y=maletimeblock_mean), color = "#0072B2", size = 7.5) +
  xlab("Time Block") + ylab("Proportion of Red \nMales among Maters") + 
  ylim(0.44, 0.61) + 
  geom_vline(xintercept = c(seq(1.5, 6.5)), color = "grey") + 
  scale_x_continuous(breaks = 1:7, 
                     labels = c("Morning\nDay 13", "Afternoon\nDay 13", "Evening\nDay 13", 
                                "Morning\nDay 14", "Afternoon\nDay 14", "Evening\nDay 14", 
                                "Morning\nDay 15")) + 
  theme_classic() +
  theme(plot.title.position = c("panel"),
        legend.title = element_blank(),
        legend.position = c("none"),
        # legend.box.background = element_rect(),
        legend.text = element_text(size=30, margin = margin(10,0,10,0), color = "black", vjust = -0.2),
        axis.text.x = element_text(size=17, margin = margin(5,0,0,0), color = "black"),
        axis.text.y = element_text(size=15, margin = margin(0,5,0,0), color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=24, margin = margin(0,10,0,10), color = "black"),
        plot.title = element_text(size=40, margin = margin(0,0,0,0), color = "black"),
        plot.margin = margin(6,6,6,6),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

#######



# Save Fig.1
png(file = "~/Desktop/UofT/SSAV_RNA/Plots/final_2/png_version/Fig1_main.png",   # The directory you want to save the file in
    width = 20, # 20 15 The width of the plot in inches
    height = 12, # 12 8 The height of the plot in inches
    units = "in", res = 300)

ggarrange(plot.globfx,
          NA,
          ggarrange(male_time.block + theme(axis.text.x = element_blank()),
                    fem_time.block,
                    nrow = 2, heights = c(0.9, 1),
                    labels = c("B)", "C)"),
                    font.label = list(size = 30)),
          ncol = 3, widths = c(0.75, 0.075, 1.5),
          labels = "A)", font.label = list(size = 30), hjust = -0.01
)

dev.off()
