#########R code for epiandrosterone validation plots - Sen et al 2022#########

###################### Setting up R  and loading packages######################
rm(list = ls())
library(tidyr)
library(ggthemes)
library(bbmle)
library(lubridate)
library(dplyr)

#function for making transparent plots - obtained from Andy Marshall
makeTransparent <- function(color.var, alpha = my.alpha){
  newcolor.var <- col2rgb(color.var)
  apply(newcolor.var, 2, function(curcoldata)			
  {rgb(red   = curcoldata[1],
       green = curcoldata[2],
       blue  = curcoldata[3],
       alpha = alpha,
       maxColorValue =  255)})
}


#read in epiandrosterone values from assays
assay_results                 <- readRDS("/Users/sharmisen/Desktop/CURRENT MANUSCRIPT PROJECTS/1.T VALIDATION PAPER/T_MS_Input_Files/Sen_et_al_2022_validation_T_Values.RData")

#monthly ILH for all juve males
demo                          <- read.csv("/Users/sharmisen/Desktop/CURRENT MANUSCRIPT PROJECTS/1.T VALIDATION PAPER/T_MS_Input_Files/Dispersal_model.csv") 

#how many males are there in the androgen dataset = 70
length(unique(assay_results$male.ID))

assay_results$age.sample <- round(assay_results$age.sample,0)
#range of samples per male  1-83
sd(apply(table(assay_results$male.ID, assay_results$age.sample),1, sum))

#mean number of samples per age = 16
sd(apply(table(assay_results$male.ID, assay_results$age.sample),2, sum))

#range of samples per age group = 1-46 
range(table(assay_results$age.sample)) 

###############################################################
#   Analytical and biological validation plots                #
###############################################################

#Linearity and parallelism plot - Fig S1
parallelism_epi               <- read.csv("/Users/sharmisen/Desktop/CURRENT MANUSCRIPT PROJECTS/1.T VALIDATION PAPER/T_MS_Input_Files/T_MS_Parallelism_Plot.csv")
par(mar=c(4,4,1,1))
jpeg("HB_Supp_Fig1.jpeg", width = 16, height = 14, res = 340, units = "cm")
plot(x = parallelism_epi$log.pg.ul.[parallelism_epi$Category == "STD"] + 0.2, y = parallelism_epi$X.B[parallelism_epi$Category == "STD"] + 2, pch = 16, col = makeTransparent("black", 200), ylim = c(0,100), xlab = "Log Dose (pg/ul)", ylab = "Percent Binding", cex = 1.5, xlim = c(0, 3.5), axes = F)  
par(new = T)
lines(ylim = c(0,100), xlab = " ", ylab = " ", xlim = c(0, 3), x = parallelism_epi$log.pg.ul.[parallelism_epi$Category == "STD"] + 0.2, y = parallelism_epi$X.B[parallelism_epi$Category == "STD"] + 2, pch = 18, cex = 1.5, col = makeTransparent("black", 250), lty = "dotted", lwd = 1.5)
par(new = T)
points(ylim = c(0,100), xlab = " ", ylab = " ", xlim = c(0, 3), x = parallelism_epi$log.pg.ul.[parallelism_epi$Category == "Dil_Set_1"], y = parallelism_epi$X.B[parallelism_epi$Category == "Dil_Set_1"], pch = 17, cex = 1.5, col = makeTransparent("#D55E00", 250))  
par(new = T)
lines(ylim = c(0,100), xlab = " ", ylab = " ", xlim = c(0,3), x = parallelism_epi$log.pg.ul.[parallelism_epi$Category == "Dil_Set_1"] , y = parallelism_epi$X.B[parallelism_epi$Category == "Dil_Set_1"], pch = 17, cex = 1.5, col = makeTransparent("#D55E00", 250), lty = "dotted", lwd = 1.5)  
axis(side = 2, at = seq(0,100, by = 20), labels = seq(0,100, by = 20),las = 2, cex = 1.2)
axis(side = 1, at = seq(0,3, by = 0.5), labels = seq(0,3, by = 0.5), cex = 2)
legend(x = 2, y = 90, legend = c("Standard", "Male pool dilutions"), bty = "n",pch = c(16,17), col = c("black", "#D55E00"))
text(x = parallelism_epi$log.pg.ul.[parallelism_epi$Category == "STD"], y = parallelism_epi$X.B[parallelism_epi$Category == "STD"], labels = parallelism_epi$Sample[parallelism_epi$Category == "STD"], pos = 2, adj = 1, offset = 2)
dev.off()

#Accuracy plot - Fig S2
accu_plot1                    <- read.csv("/Users/sharmisen/Desktop/CURRENT MANUSCRIPT PROJECTS/1.T VALIDATION PAPER/T_MS_Input_Files/Accuracy1_EpiAndro_9.19.2019.csv")
accu_plot1                    <- accu_plot1[1:21,1:4]
accu_plot1$log.pg.ul.         <- log(accu_plot1$Conc)
range(accu_plot1$log.pg.ul.)  

plot.new()
par(mar=c(4,4,1,1))
jpeg("HB_Supp_Fig2.jpeg", width = 20, height = 14, res = 340, units = "cm")
plot(x = accu_plot1$log.pg.ul.[accu_plot1$Type == "Standard"], y = accu_plot1$Binding[accu_plot1$Type == "Standard"], pch = 16, col = makeTransparent("black", 200),  cex = 1.5,  ylim = c(0,80), xlab = "Log Dose (pg/ul)", ylab = "Percent Binding", xlim = c(1,6), axes = F)  
lines(x = accu_plot1$log.pg.ul.[accu_plot1$Type == "Standard"], y = accu_plot1$Binding[accu_plot1$Type == "Standard"], lty = "dotted", col = makeTransparent("black", 200), lwd = 2)
text(x = accu_plot1$log.pg.ul.[accu_plot1$Type == "Standard"], y = accu_plot1$Binding[accu_plot1$Type == "Standard"], labels = accu_plot1$Sample[accu_plot1$Type == "Standard"], pos = 2, adj = 1, offset = 1)
par(new = T)
points(ylim = c(0,80), xlab = " ", ylab = " ", xlim = c(1,6), x = accu_plot1$log.pg.ul.[accu_plot1$Type == "Spiked1"], y = accu_plot1$Binding[accu_plot1$Type == "Spiked1"], cex = 1.5, pch = 17, col = makeTransparent("#D55E00", 250))  
lines(ylim = c(0,80), xlab = " ", ylab = " ", xlim = c(1,6), x = accu_plot1$log.pg.ul.[accu_plot1$Type == "Spiked1"], y = accu_plot1$Binding[accu_plot1$Type == "Spiked1"],col = makeTransparent("#D55E00", 250), lwd = 2, lty = "dotted")
axis(side = 2, at = seq(0,80, by = 20), labels = seq(0,80, by = 20),las = 2, cex =1.2)
axis(side = 1, at = seq(1,6, by = 1), labels = seq(1,6, by = 1), cex =1.2)
legend(x = 4, y = 70, legend = c("Standard", "Standards spiked with \n pool 1:160"), bty = "n",pch = 16, col = c(makeTransparent("black",200), "#D55E00"))
dev.off()

###############################....................###########################

###############################################################
#                   EIA vs RIA comparisons                    #
###############################################################

#comparing two different antibodies to quantify testosterone - MPBio and Epiandrosterone
old_values                   <- read.csv("/Users/sharmisen/Desktop/CURRENT MANUSCRIPT PROJECTS/1.T VALIDATION PAPER/T MS related OLD/T_ANALYSIS FILE_OLD with RIA values.csv")
ria                          <- old_values[old_values$ID %in% demo$male_ids,] 
ria                          <- ria[which(!is.na(ria$ria.ng.g)),]
#pulling out samples that have both ria and eia values - from only natal males
ria$eia.ng.g                 <- assay_results$ng.g[match(ria$Sample, assay_results$Sample)]
ria                          <- ria[which(ria$Sample %in% assay_results$Sample),]
ria                          <- ria[ria$ID != "JUD",] #removing JUDE 
cor.test(log(ria$ria.ng.g),log(ria$eia.ng.g))

ria$ng.g.RIA                 <- ria$ria.ng.g
ria$ng.g.EIA                 <- ria$eia.ng.g
#get  male ID and age for each sample
ria$age                      <- assay_results$age.sample[match(ria$Sample, assay_results$Sample)]
ria$male                     <- assay_results$male.ID[match(ria$Sample, assay_results$Sample)]
ria                          <- ria[which(!is.na(ria$age)),]
ria$age                      <- round(ria$age, 0)

#sampling distribution

#how many males are there in the androgen dataset = 49
length(unique(ria$male))

#sampling distribution per male id per age
table(ria$male, ria$age)

#range of samples per male per age 1-74
range(apply(table(ria$male, ria$age),1, sum))

#mean number of samples per male per age = 21.38, 
#standard deviation per male per age = 21.2, 
#range per male per age = 1-74
mean(sort(apply(table(ria$male, ria$age),1, sum)))
sd(sort(apply(table(ria$male, ria$age),1, sum)))
range(sort(apply(table(ria$male, ria$age),1, sum)))

#range of samples per age group 1 - 262 - age is rounded
range(table(ria$age)) 
sd(table(ria$age))
mean(table(ria$age))
length(unique(ria$male))

comp_summ1                   <- ria %>% 
                                group_by(age, male)%>% 
                                summarise(mean.male.EIA = mean(ng.g.EIA, na.rm = T), 
                                sd.T.EIA = sd(ng.g.EIA, na.rm = T), mean.male.RIA = mean(ng.g.RIA, na.rm = T), 
                                sd.T.RIA = sd(ng.g.RIA, na.rm = T)) 
comp_summ1

comp_summ2                   <- comp_summ1 %>% 
                                group_by(age) %>% 
                                summarize(mean.age.T.RIA = mean(mean.male.RIA, na.rm = T), 
                                sd.RIA = sd(mean.male.RIA, na.rm = T),
                                mean.age.T.EIA = mean(mean.male.EIA, na.rm = T), 
                                sd.EIA = sd(mean.male.EIA, na.rm = T),
                                nmale = n()) %>% 
                                mutate(se.RIA = sd.RIA/sqrt(nmale),
                                se.EIA = sd.EIA/sqrt(nmale)) %>% as.data.frame()

#Making two panel plot comparing MPBio and EpiAndro across male ages - Fig S3
jpeg("HB_Supp_Fig3.jpeg.jpeg", width = 22, height = 20, units = "cm", res = 340)
par(mfrow = c(1,2), mar = c(3,4,0,0), oma = c(2,2,0,0))
#plotting RIA values first
plot(x = comp_summ2$age, y = comp_summ2$mean.age.T.RIA, axes = F, col= makeTransparent("#999999", 250), pch = 15, ylim = c(0,16), xlim = c(1,9), ylab = " ", xlab = " ", cex = 2)
lines(x = comp_summ2$age, y = comp_summ2$mean.age.T.RIA, lty = "dotted", lwd = 1.5, xlim = c(0,10))
axis(side = 1, labels = seq(1,9,by = 1), at = seq(1,9,by = 1), cex.axis = 1.25)
axis(side = 2, labels = seq(0,16,by = 2), at = seq(0,16,by = 2), las =2, cex.axis = 1.25)
arrows(x0 = comp_summ2$age, y0 = comp_summ2$mean.age.T.RIA - 2 * comp_summ2$se.RIA, x1 = comp_summ2$age, y1 = comp_summ2$mean.age.T.RIA + 2 * comp_summ2$se.RIA, code=3, angle=90, length=0.02)
text(x = comp_summ2$age, y = 0.2, paste(comp_summ2$nmale), font = 3)
mtext(side = 1, text = "Male age (years)", line = 3,cex = 1.5)
text(5,16, "MP Biomedicals \n Testosterone RIA", col = makeTransparent("#999999", 250), font = 2, cex = 1.5)
mtext(side =2, text = "Fecal androgen metabolite levels (ng/g)", line = 4, cex = 1.5)

#plotting EIA values next
plot(x = comp_summ2$age, y = comp_summ2$mean.age.T.EIA, col= makeTransparent("#332288",200), axes = F,  pch = 15,  ylim = c(0,1400), xlim = c(1,9), ylab = " ", xlab = " ", cex = 2)
lines(x = comp_summ2$age, y = comp_summ2$mean.age.T.EIA, lty = "dotted", lwd = 1.5,xlim = c(1,8))
arrows(x0 = comp_summ2$age, y0 = comp_summ2$mean.age.T.EIA - 2 * comp_summ2$se.EIA, x1 = comp_summ2$age, y1 = comp_summ2$mean.age.T.EIA + 2 * comp_summ2$se.EIA, code=3, angle=90, length=0.02)
text(x = comp_summ2$age, y = 0.2, paste(comp_summ2$nmale), font = 3)
axis(side = 1, labels = seq(1,9,by = 1), at = seq(1,9,by = 1), cex.axis = 1.25)
axis(side = 2, labels = seq(0,1400,by = 200), at = seq(0,1400,by = 200), las =2, cex.axis = 1.25)
mtext(side = 1, text = "Male age (years)", line = 3, cex = 1.5)
text(x = 5,y = 1400,  "Epiandrosterone \n  Enzymeimmunoassay", col = "#332288", font = 2, cex = 1.5)
dev.off()


#Plotting the correlation between EpiAndrosterone and MPBioMedical T RIA - Fig S4
fit                         <- lm(log(ria$eia.ng.g) ~ log(ria$ria.ng.g))
summary(fit)
r1                          <- cor.test(log(ria$eia.ng.g),log(ria$ria.ng.g))$estimate
r                           <- summary(fit)$adj.r.squared #$Adjusted R-squared

par(mfrow = c(1,1))
jpeg("HB_Supp_Fig4.jpeg", width = 16, height = 14, units = "cm", res = 340)
plot(x = log(ria$ria.ng.g) , y = log(ria$eia.ng.g), axes = F, pch = 16, cex = 1.5, col = makeTransparent("#332288",100), ylab = "", xlab = "", xlim = c(0,5), ylim = c(0, 8))
mtext(side = 1, text = "MP Biomedicals Testosterone RIA log(ng/g)", cex = 1.5, line = 2.5)
mtext(side = 2, text = "EpiAndrosterone EIA log(ng/g)", cex = 1.5, line = 2.5)
axis(side = 2, at = seq(0,8, by = 1), labels = seq(0,8,by = 1), las = 2, cex.axis= 1.25)
axis(side = 1, at = seq(0,5, by = 1), labels = seq(0, 5, by = 1), cex.axis= 1.25)
abline(fit, lty = "dotted", lwd = 2)
mylabel = bquote(expr = paste(italic(R)^2 == .(format(r, digits = 3)) ,"\n", italic(r) == .(format(r1, digits = 3))))
legend("right", bty="n", legend = mylabel)
dev.off()

###############################....................###########################

###############################################################
#                     Biological validation                   #
###############################################################

#Summarizing T across age  - finding the average per male per age category where the sample size = the number of males in that age category
assay_results$age.sample      <- round(assay_results$age.sample,0)
assay_results <- assay_results[assay_results$age.sample<11,]
sd(table(assay_results$age.sample))

#number of males 
length(unique(assay_results$male.ID))
#sampling distribution per male 
mean(sort(apply(table(assay_results$male.ID, assay_results$age.sample),1, sum)))
sd(sort(apply(table(assay_results$male.ID, assay_results$age.sample),1, sum)))
range(sort(apply(table(assay_results$male.ID, assay_results$age.sample),1, sum)))

testo_summ                    <- assay_results %>% 
                                 group_by(age.sample, male.ID) %>% 
                                 summarise(mean.male.T = mean(ng.g), sd.T = sd(ng.g)) 

testo_summ2                   <- testo_summ %>% 
                                group_by(age.sample) %>% 
                                summarize(mean.age.T = mean(mean.male.T, na.rm = T), 
                                sd = sd(mean.male.T, na.rm = T), nmale = n()) %>% 
                                mutate(se = sd/sqrt(nmale)) %>% 
                                as.data.frame()
testo_summ2 <- testo_summ2[testo_summ2$age.sample <11,]

#Plotting T across male age - Fig 1 (Main text)
jpeg("HB_Fig1.jpeg", width = 16, height = 12, units = "cm", res = 320)
par(mar = c(4,4,2,2), oma = c(2,2,0,2))
plot(testo_summ2$age.sample, y = testo_summ2$mean.age.T, pch = 15, col = "white", ylim = c(0,1400), xlim = (c(0,10)), axes = F, xlab = " ", ylab = " ")
rect(3.7,-10,8.9,1400, col= makeTransparent("gray", 200), border = NA)
points(x = testo_summ2$age.sample, y = testo_summ2$mean.age.T, pch = 15, cex = 1.2, col = makeTransparent("#332288",200), ylim = c(0,1400), xlim = (c(0,10)))
lines(testo_summ2$age.sample, y = testo_summ2$mean.age.T, lty = "dotted", lwd = 1.5)

#Adding axis
axis(side = 1, labels = c("<1", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"), at = 0:10, cex.axis = 1.25)
axis(side = 2, labels = seq(0,1400, by = 200), at = seq(0,1400, by = 200), las = 2, cex.axis = 1.25)
mtext(side = 1, text = "Male age (years)", line = 3, cex = 1.5)
mtext(side = 2, text = "Fecal androgen metabolites (ng/g)", line = 4, cex = 1.5)

#Adding error bars
arrows(x0 = testo_summ2$age.sample, y0 = testo_summ2$mean.age.T - 2 * testo_summ2$se, x1 = testo_summ2$age.sample, y1 = testo_summ2$mean.age.T + 2 * testo_summ2$se, code=3, angle=90, length=0.02)
abline(v = 6.6, col = "black", lwd = 2)
text(x = testo_summ2$age.sample, y = 10,paste(testo_summ2$nmale), font = 3)
dev.off()

###############################....................###########################

############################################################################
#       Seasonality plots - weather and T plots across the year            #
############################################################################
df                    <- assay_results
age.model.all         <- lm(ng.g ~ age.sample, data = df)
res.fam.all           <- residuals(age.model.all) 
df                    <- cbind(df, res.fam.all)
df$sample.month       <- month(df$sample.date)

season_temp           <- df %>% 
                         group_by(sample.month, male.ID) %>% 
                         summarise(mean.T = mean(ng.g, na.rm = T), 
                         se.T = sd(mean.T, na.rm = T)/sqrt(n()),
                         mean_res = mean(res.fam.all, na.rm = T), 
                         se.res = sd(res.fam.all, na.rm = T)/sqrt(n()), 
                         nsamples = n()) %>% as.data.frame()

season_temp

season                <- df %>% 
                         group_by(sample.month) %>% 
                         summarise(avg.min.Temp = mean(min.Temp, na.rm = T), 
                         se.min.Temp = sd(min.Temp, na.rm = T)/sqrt(n()),
                         avg.max.Temp = mean(max.Temp, na.rm = T), 
                         se.max.Temp = sd(max.Temp, na.rm = T)/sqrt(n()),
                         avg.rain = mean(cum.rain, na.rm = T), 
                         se.rain = sd(cum.rain, na.rm = T)/sqrt(n())) %>% as.data.frame()
season_all            <- season

monthly_T_all         <- season_temp %>%   
                         group_by(sample.month) %>% 
                         summarise(mean.T.month = mean(mean.T, na.rm = T), 
                         se.monthly.T = sd(mean.T, na.rm = T)/sqrt(n()),
                         mean.res.month = mean(mean_res, na.rm = T), 
                         se.res = sd(mean_res, na.rm = T)/sqrt(n()), 
                         monthly_nmales = n()) %>%as.data.frame()

monthly_T_all
season_all <- cbind(monthly_T_all, season_all)
season_all <- season_all[,-7]

#Plotting season plots with actual T values from ALL natal males - Fig S5
jpeg("HB_Supp_Fig5.jpg", width = 20, height = 14, units = "cm", res = 340)
par(mar = c(2, 4, 3, 4), oma = c(2,3,2,3))       
#Create first plot with rainfall
barcenters <- barplot(season_all$avg.rain ~ season_all$sample.month, data = season_all, 
                      border = NA, col = makeTransparent("#888888", 100), axes = F, ylab = " ", 
                      xlab = " ", ylim = c(0,700), names.arg = " ")
segments(barcenters, season_all$avg.rain - season_all$se.rain, barcenters,season_all$avg.rain + season_all$se.rain, lwd = 1.5)
arrows(barcenters, season_all$avg.rain - season_all$se.rain,barcenters , y1 = season_all$avg.rain + season_all$se.rain, code = 3, angle = 90, length = 0.02,lwd = 1.5, col = makeTransparent("black", 150))

#Create second plot with min temp
par(new = TRUE)
plot(season_all$sample.month, season_all$avg.min.Temp, pch = 25, bg = makeTransparent("#D55E00",200), col = makeTransparent("#D55E00",200), axes = F, ylab = " ", xlab = " ", ylim = c(4,24), cex = 1.5)
lines(season_all$sample.month, season_all$avg.min.Temp, lty = "dotted", col = makeTransparent("#D55E00",200))

axis(side =1, at = seq(1,12,by = 1), labels = c("Jan","Feb","Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), cex.axis = 1.5)
axis(side = 4, at = seq(4,24,by = 4), labels = seq(4,24,by = 4), cex.axis = 1.5, las = 2, col = makeTransparent("#D55E00",250))
mtext(text = paste0("Monthly temperature (\u00B0C)   ", "Max      ", "Min      ") , side = 4,  line = 4, cex = 1.5)

#Create third plot with max temp
par(new = TRUE)
plot(season_all$sample.month, season_all$avg.max.Temp, pch = 24, bg = makeTransparent("#D55E00",200), col = makeTransparent("#D55E00",200),axes = F,  xlab = " ", ylab = " ", ylim = c(4,24), cex = 1.5)
lines(season_all$sample.month, season_all$avg.max.Temp, lty = "dotted", col = makeTransparent("#D55E00", 200))

#Add new plot with mean androgen values
par(new = TRUE)                            
plot(season_all$sample.month, season_all$mean.T.month, col = makeTransparent("#332288",250), pch = 15, axes = FALSE, xlab = " ", ylab = " ",ylim = c(0,700), cex = 1.5)
arrows(x0 = season_all$sample.month, y0 = season_all$mean.T.month - season_all$se.monthly.T, x1 = season_all$sample.month, y1 =  season_all$mean.T.month + season_all$se.monthly.T, code=3, angle=90, length=0.02, col = makeTransparent("black",250))
lines(season_all$sample.month, season_all$mean.T.month, col = makeTransparent("black",150), lwd = 1.5, lty = "dashed")
axis(side = 2, at = seq(0,700, by = 100), labels = seq(0,700, by = 100), cex.axis = 1.5, las = 1)      

# Add second axis
mtext("Mean monthly fecal androgen metabolites (ng/g) or \n cumulative rainfall (mm across previous 90 days)", side = 2, line = 4, cex = 1.5)
text(x = 1:12, y = 0, labels = paste(season_all$monthly_nmales), font = 3, cex = 1)
dev.off()
###############################....................###########################

#plot distribution of outcome variable  (fAMs) - Fig S6
df_models <- read.csv("/Users/sharmisen/Desktop/CURRENT MANUSCRIPT PROJECTS/1.T VALIDATION PAPER/Androgen_data_Sen_et_al_2022.csv")
jpeg("HB_Supp_Fig6.jpg", width = 18, height = 16, res = 340, units = "cm")
par(mar = c(4,3,2,0))
par(mfrow = c(1,2))
hist(df_models$ng.g, col= makeTransparent("#332288",200), main = " ", xlab = "", ylab = " ", breaks = 20, cex.axis = 1, ylim = c(0,400), border = makeTransparent("#332288",220), las = 1)
mtext(side = 1, text ="Fecal androgen metabolite \n levels (ng/g)", line = 3)
mtext(side = 2, text ="Frequency", cex = 1, line = 3)
text(x = 50, y = 400, labels ="A", cex = 1,las = 1, font = 2)
qqnorm(df_models$ng.g, col = makeTransparent("#332288",200),  pch = 16, main = "", axes = F)
axis(side = 1, at = seq(-3,3, 1), labels = seq(-3,3, 1) , las = 1)
axis(side = 2, at = seq(0,2500, 500), labels = seq(0,2500, 500), las = 1)
qqline(df_models$ng.g)
text(x = -3, y = 2420, labels ="B", cex = 1, las = 1, font = 2)
dev.off()

###############################....................###########################

############################################################################
#                       Androgens around dispersal                         #
############################################################################
androgen_data            <- readRDS("/Users/sharmisen/Desktop/CURRENT MANUSCRIPT PROJECTS/1.T VALIDATION PAPER/T_MS_Input_Files/Androgens_from_dispersed_males.RData")
androgen_data$year       <- year(androgen_data$sample.date)
androgen_data$time_dis   <- round(time_length(difftime(androgen_data$sample.date, androgen_data$dis_date), "month"), 1)
androgen_data$age_dis    <- round(time_length(difftime(androgen_data$dis_date, androgen_data$dob), "years"), 1)

#get number of unique males = for 22 males we have pre and post dispersal samples 
length(unique(androgen_data$male.ID))

#Residual androgens 
model1                   <- lm(log(ng.g) ~ age.sample, data = androgen_data)
androgen_data$dis_status <- as.factor(androgen_data$dis_status)
res_andro                <- residuals(model1)
plot(res_andro)

#add residuals to the dataframe
androgen_data <- cbind(androgen_data, res_andro)

#get the subset of males around dispersal age - so between 12 months before and after dispersal
subset_disp             <- androgen_data[androgen_data$time_dis > -3 & androgen_data$time_dis <= 3,]
subset_disp$time_dis    <- round(subset_disp$time_dis,0)

#plot mean T before and after dispersal
subset_disp <- subset_disp %>% mutate(dis_cat = case_when(age_dis <= 5.5 ~ "Early", age_dis >5.5 & age_dis <= 7.7 ~ "Mean", age_dis > 7.7 ~ "Late"))
table(subset_disp$dis_cat)
subset_disp$dis_cat <- factor(subset_disp$dis_cat, levels = c("Early", "Mean", "Late"))

#male IDs with samples at dispersal
unique(subset_disp$male.ID[subset_disp$time_dis == 0])

#Plot overall residual log androgen levels around 3 months from dispersal 
jpeg("HB_Supp_Fig10.jpg", height = 12, width = 10, units = "cm", res = 340)
ggplot(data = subset_disp, aes(x = factor(time_dis), y = res_andro)) + 
  geom_boxplot() + geom_point() + theme_pubr(base_size = 14) + xlab("Months relative dispersal") + 
  ylab("Residual log fAMs ng/g") + geom_vline(xintercept = 4, linetype = "dashed", color = "#999999", size = 1.2)  
dev.off()


#residual T levels around dispersal overall 
subset_disp1 <- subset_disp %>% 
                group_by(male.ID, time_dis, dis_cat) %>% 
                summarise(mean.ng.g = mean(ng.g), nsamples = n()) %>% 
                group_by(time_dis, dis_cat) %>% mutate(nmales = n()) 

labels          <- subset_disp1 %>% group_by(time_dis, dis_cat) %>% 
                   summarise(nmales = n()) %>% arrange(dis_cat)
labels
jpeg("HB_Fig4.jpg" , height = 15, width = 18, units = "cm", res = 340)
ggplot(data = subset_disp1, aes(x = factor(time_dis), y = log(mean.ng.g),  col = factor(dis_cat))) + 
  geom_boxplot(position = "dodge") + geom_point() + geom_hline(yintercept = 0, size = 0.5, linetype = "dashed") + 
  theme_pubr(base_size = 18) + xlab("Months relative to dispersal") + ylab("log fAMs (ng/g)") + 
  geom_vline(xintercept = 4, size = 1, linetype = "dashed", col = "#999999") + ylim(c(4,8)) + 
  scale_color_manual(" ", values = c("#E69F00","#332288","#999999")) + 
  theme(legend.background = element_rect(color = NA)) + facet_wrap(.~ dis_cat)
dev.off()

#how many males in early, mean and late?
length(unique(subset_disp$male.ID[subset_disp$dis_cat == "Early"]))
length(unique(subset_disp$male.ID[subset_disp$dis_cat == "Mean"]))
length(unique(subset_disp$male.ID[subset_disp$dis_cat == "Late"]))

#within the same male, see a paired comparison of diff in T pre and post dispersal
paired_data <- read.csv("/Users/sharmisen/Desktop/CURRENT MANUSCRIPT PROJECTS/1.T VALIDATION PAPER/paired_dispersal_data.csv")
paired_data  

paired_data$before <- log(paired_data$before)
paired_data$after  <- log(paired_data$after)

jpeg("HB_Supp_Fig11.jpg", width = 12, height = 10, units = "cm", res = 340)
ggpaired(paired_data, cond1 = "before", cond2 =  "after", 
         xlab = "Dispersal", ylab = "log fAMs (ng/g)", line.color = makeTransparent("black",200), 
         line.size = 0.4, point.size = 3, ggtheme = theme_pubr(base_size = 16)) 
dev.off()

#wilcoxon rank test 
paired_stats <- read.csv("/Users/sharmisen/Desktop/CURRENT MANUSCRIPT PROJECTS/1.T VALIDATION PAPER/T MS Figures/paired_dispersal_data_stats.csv")
wilcox.test(ng.g ~ cat, data = paired_stats, paired = TRUE, alternative = "greater")

###############################....................###########################

