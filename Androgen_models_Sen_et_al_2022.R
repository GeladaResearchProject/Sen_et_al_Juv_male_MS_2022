#########R code for running androgen models - Sen et al 2022#########

###################### Setting up R ######################
rm(list = ls())

#set working directory
setwd("/Users/sharmisen/Desktop/CURRENT MANUSCRIPT PROJECTS/T VALIDATION PAPER/T_MS_Input_Files")

#Loading packages
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(bbmle)
library(lubridate)
library(lme4)
library(VIM)
library(lmerTest)
library(sjPlot)
library(MuMIn)
library(car)
library(broom)
library(glmm)
library(jtools)
library(officer)
library(flextable)
library(gtsummary)
library(gridExtra)
library(cowplot)

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

#function to standardize variables
z_transform <- function(var){
  value <- (var - mean(var, na.rm = T))/ (2 * sd(var, na.rm = T))
  return(value)
}

######################....................###########################

#################Reading in Androgen dataframe ######################

df_models <- read.csv("/Users/sharmisen/Desktop/CURRENT MANUSCRIPT PROJECTS/T VALIDATION PAPER/Androgen_data_Sen_et_al_2022.csv")
df_models <- df_models %>% dplyr::select(male.ID, ng.g, age.sample, mrank_birth, unit.size.cont, 
                                         pat.sibs.sample, cum.rain, min.Temp, 
                                         max.Temp, year)

#correlation between variables: max Temp and min Temp = 0.65
cor.test(df_models$max.Temp, df_models$min.Temp)

#correlation between variables: max Temp and cum Rainfall = -0.42
cor.test(df_models$max.Temp, df_models$cum.rain)

#correlation between variables: min Temp and cum Rainfall = 0.1
cor.test(df_models$min.Temp, df_models$cum.rain)

#correlation between variables: unit size and cohort size = 0.22
cor.test(df_models$unit.size.cont, df_models$pat.sibs.sample)

#Percentage of missing data in this dataste = 37%  - mostly because we don't have maternal ranks at birth for certain females in units in 2006-2009 
#when the project started following units or added them later on in 2009
mean(!complete.cases(df_models)) * 100

#setting the age cutoff for splitting data for modelling androgens separately 
cutoff                           <- 2.5
######################....................###########################


###############################################################
#   Androgen models for males <= 2.5 years : pre-independence #
###############################################################

young_males                      <- df_models[df_models$age.sample <= cutoff,]
head(young_males)

#number of males in this dataset
length(unique(young_males$male.ID))

#storing age variable for plotting purposes
young_males$age                  <- young_males$age.sample

#z-transforming all continuous predictors
young_males$age.sample           <- z_transform(young_males$age.sample)
young_males$pat.sibs.sample      <- z_transform(young_males$pat.sibs.sample)
young_males$cum.rain             <- z_transform(young_males$cum.rain) 
young_males$max.Temp             <- z_transform(young_males$max.Temp) 
young_males$mrank_birth          <- z_transform(young_males$mrank_birth)
young_males$unit.size.cont       <- z_transform(young_males$unit.size.cont)
young_males$year                 <- factor(young_males$year)

#log transforming outcome variable - fAMs
young_males$log.ng.g             <- log(young_males$ng.g)

#removing missing data
young_males                      <- na.omit(young_males)

#androgen model for males <= 2.5 years using lmer() and normal distribution
young_model_birth_lmer           <- lmer(log.ng.g ~ age.sample + mrank_birth + pat.sibs.sample +
                                         unit.size.cont + max.Temp + cum.rain + (1|male.ID), 
                                         REML = F, na.action = "na.fail", data = young_males)

#check variance inflation scores - all under 1.5
vif(young_model_birth_lmer)
summary(young_model_birth_lmer)

#dredging full model 
dredged_young                    <- dredge(young_model_birth_lmer, rank = "AICc")

#performing model averaging
avg_young                        <- model.avg(get.models(dredged_young, subset = T), type = 'link', backtransform = TRUE)
summary(avg_young)

###############################################################
#  Androgen models for males > 2.5 years : post-independence  #
###############################################################
older_males                      <- df_models[df_models$age.sample > cutoff, ]
head(older_males)

#removing all missing data
older_males <- na.omit(older_males)

#number of males in this dataset
length(unique(older_males$male.ID))

#storing age variable for plotting purposes
older_males$age                  <- older_males$age.sample

#z-transforming all continuous predictors
older_males$mrank_raw            <- older_males$mrank_birth #storing this predictor to plot rank effects on fAMs
older_males$mrank_birth          <- z_transform(older_males$mrank_birth)
older_males$age.sample           <- z_transform(older_males$age.sample)
older_males$pat.sibs.sample      <- z_transform(older_males$pat.sibs.sample)
older_males$cum.rain             <- z_transform(older_males$cum.rain) 
older_males$max.Temp             <- z_transform(older_males$max.Temp) 
older_males$unit.size.cont       <- z_transform(older_males$unit.size.cont)
older_males$year                 <- factor(older_males$year)

#log transforming outcome variable - fAMs
older_males$log.ng.g             <- log(older_males$ng.g)

#androgen models for juveniles using lmers() and normal dis     

full_model_older_lmer            <- lmer(log.ng.g ~ age.sample  + mrank_birth +  pat.sibs.sample +  
                                        unit.size.cont + cum.rain + max.Temp + (1|male.ID), 
                                         REML = F, na.action = "na.fail",data = older_males)
summary(full_model_older_lmer)
vif(full_model_older_lmer)
fm1_lmer_old                     <- dredge(full_model_older_lmer, rank = "AICc")
avg_old                          <- model.avg(get.models(fm1_lmer_old, subset = T), type = 'link', backtransform = TRUE)
summary(avg_old)

######################....................###########################

######################....................###########################
#              Androgen model coefficient plots                     #
######################....................###########################
par(mar = c(2,2,0,0), omar = c(1,2,1,1))
plot1                            <- plot_model(avg_young, type = "std2", order.terms = c(1,6,2,5,4,3), 
                                      axis.labels =  c("age.sample" = "Age", "mrank_birth" = "Maternal rank",  
                                     "pat.sibs.sample" = "Cohort size", "unit.size.cont" = "Unit size",  
                                     "max.Temp" = "Max Temp (\u00B0C)", "cum.rain" = "Cum rain (mm)"),
                                      wrap.labels = F, pred.type = "fe", value.offset = 0.3, value.size = 4.5, 
                                      show.values = T, show.p = T, line.size = 0.2, title = "Androgen model \n(males < 2.5 years)", 
                                      col = "black",  vline.color = "gray", ci.style = "whisker") + theme_pubr(base_size = 14) + 
                                      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

plot2                            <- plot_model(avg_old, type = "std2", order.terms = c(1,2,6,4,3,5), 
                                      axis.labels =  c("age.sample" = "        ", "mrank_birth" = "        ",
                                     "pat.sibs.sample" = "      ","unit.size.cont" = "        ", 
                                     "max.Temp" = "       ", "cum.rain" = "       "),
                                      wrap.labels = F, pred.type = "fe", value.offset = 0.3, value.size = 4.5, 
                                      show.values = T, show.p = T, line.size = 0.2, 
                                      title = "Androgen model \n(males > 2.5 years)", col = "black",  
                                      vline.color = "gray", line.type = "dotted", ci.style = "whisker") + theme_pubr(base_size = 14) + 
                                      theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

pdf("updated_fig_lmer_combined.pdf",width = 8, height = 6)
grid.arrange(plot1, plot2, ncol = 2, widths = c(95, 75))
dev.off()

######################....................###########################

###################################################################
#               Obtaining model coefficient tables                #
###################################################################

#Model coefficients for average dredged models from fAM values obtained from young males
avg_young %>% tbl_regression(label = list(age.sample = "Age",  mrank_birth = "Maternal rank",  
                                                pat.sibs.sample = "Cohort size", unit.size.cont = "Unit size",
                                                max.Temp = "Max Temp (\u00B0C)", cum.rain = "Cum rain (mm)"), order.terms = c(1,6,2,5,3,4)) %>% bold_p() 
coefTable(avg_young)

avg_old %>% tbl_regression() %>% bold_labels() %>% bold_p() %>% italicize_levels() 
coefTable(avg_old)

#Plotting model effects of rank for older males    

#obtaining residuals from the older male androgen full model 
res.andro                     <- residuals(full_model_older_lmer)
older_males                   <- cbind(older_males, res.andro)

#changing proportional rank into an ordinal ranks for visualization
older_males$cat.rank          <- ifelse(older_males$mrank_raw > 0.5, "High", "Low")
head(older_males)

jpeg("T_rank_effect.jpeg", width = 8, height = 6, units = "in", res = 320)
ggplot(data = older_males, aes(y = res.andro, x = age, fill = cat.rank, col = cat.rank)) + 
  geom_point(size = 2, show.legend = F) + scale_color_manual(values = c(makeTransparent("#332288",200), makeTransparent("#999999", 200))) + 
  ylab("Residual log fAMs (ng/g)\n") + xlab("\nMale age in years") + theme_pubr(base_family = "Arial", base_size = 16) + theme(legend.text=element_text(size=16), axis.title = element_text(size=16)) +
  geom_smooth(alpha = 0.5, color = "black", size = 0.8) + scale_fill_manual("Maternal rank at birth", values = c(makeTransparent("#332288",100), makeTransparent("#999999", 150))) + 
  scale_x_continuous(breaks = seq(2.5,8.5,1)) + scale_y_continuous(breaks = seq(-1.5, 1.5,0.5)) 
dev.off()

######################....................###########################

