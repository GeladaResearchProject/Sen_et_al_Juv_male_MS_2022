######### Androgen models - Sen et al 2022#########

###################### Setting up R ######################
rm(list = ls())

#set working directory
setwd("/Users/sharmisen/Desktop/CURRENT MANUSCRIPT PROJECTS/1.T VALIDATION PAPER/T_MS_Input_Files")

#Loading packages
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(lme4)
library(sjPlot)
library(MuMIn)
library(car)
library(officer)
library(flextable)
library(gtsummary)
library(gridExtra)


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

#################Reading in Androgen dataframe######################

df_models <- read.csv("/Users/sharmisen/Desktop/CURRENT MANUSCRIPT PROJECTS/1.T VALIDATION PAPER/T_MS_Input_Files/revised_androgen_dataset_6.25.2022.csv")
df_models <- df_models %>% dplyr::select(male.ID, ng.g, sample.date, age.sample, mrank_birth, mrank_prop, unit.size.cont, 
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

#correlation between variables: mrank at birth and mrank each month = 0.66
cor.test(df_models$mrank_birth, df_models$mrank_prop)

#Percentage of missing data in this dataset = 37.5%  - mostly because we don't have maternal ranks at birth for certain females in units in 2006-2009 
#when the project started following units or added them later on in 2009
mean(!complete.cases(df_models)) * 100


#setting the age cutoff for splitting data for modelling androgens separately 
cutoff                           <- 2.5

######################....................###########################


#####################################################################
#   Androgen models for males <= 2.5 years : pre-independence       #
#####################################################################
young_males                      <- df_models[df_models$age.sample <= cutoff,]
head(young_males)

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
young_males$mrank_prop           <- z_transform(young_males$mrank_prop)

#log transforming outcome variable - fAMs
young_males$log.ng.g             <- log(young_males$ng.g)

#removing missing data
young_males                       <- na.omit(young_males)

#number of males in this dataset
length(unique(young_males$male.ID))

#androgen model for males <= 2.5 years using lmer() and normal distribution
young_model_birth_lmer            <- lmer(log.ng.g ~ age.sample + mrank_birth + pat.sibs.sample +
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

#####################################################################
#  Androgen models for males > 2.5 years : post-independence        #
#####################################################################
older_males                      <- df_models[df_models$age.sample > cutoff, ]

#removing all missing data
older_males <- na.omit(older_males)

#number of males in this dataset
length(unique(older_males$male.ID))

#storing age variable for plotting purposes
older_males$age                  <- older_males$age.sample

#z-transforming all continuous predictors
older_males$mrank_raw_monthly    <- older_males$mrank_prop  #storing this predictor to plot rank effects on fAMs
older_males$mrank_raw_birth      <- older_males$mrank_birth #storing this predictor to plot rank effects on fAMs
older_males$mrank_birth          <- z_transform(older_males$mrank_birth)
older_males$age.sample           <- z_transform(older_males$age.sample)
older_males$pat.sibs.sample      <- z_transform(older_males$pat.sibs.sample)
older_males$cum.rain             <- z_transform(older_males$cum.rain) 
older_males$max.Temp             <- z_transform(older_males$max.Temp) 
older_males$unit.size.cont       <- z_transform(older_males$unit.size.cont)
older_males$year                 <- factor(older_males$year)
older_males$mrank_prop           <- z_transform(older_males$mrank_prop)

#log transforming outcome variable - fAMs
older_males$log.ng.g             <- log(older_males$ng.g)

#androgen models for juveniles using lmers() and normal dis     

full_model_older_lmer             <- lmer(log.ng.g ~ age.sample  + mrank_birth +  pat.sibs.sample +  
                                        unit.size.cont + cum.rain + max.Temp + (1|male.ID), 
                                         REML = F, na.action = "na.fail",data = older_males)
summary(full_model_older_lmer)
vif(full_model_older_lmer)
fm1_lmer_old                     <- dredge(full_model_older_lmer, rank = "AICc")
avg_old                          <- model.avg(get.models(fm1_lmer_old, subset = T), type = 'link', backtransform = TRUE)
summary(avg_old)

######################....................###########################


######################....................###########################
#         Androgen model coefficient plots + rank effects          #
######################....................###########################
plot1                             <- plot_model(avg_young, type = "std2", order.terms = c(1,6,2,5,4,3), 
                                      axis.labels =  c("age.sample" = "Age", "mrank_birth" = "Maternal rank at birth",  
                                     "pat.sibs.sample" = "Cohort size", "unit.size.cont" = "Unit size",  
                                     "max.Temp" = "Max Temp (\u00B0C)", "cum.rain" = "Cum rain (mm)"),
                                      wrap.labels = F, pred.type = "fe", value.offset = 0.3, value.size = 4.5, 
                                      show.values = T, show.p = T, line.size = 0.2, title = "Androgen model (males < 2.5 years)", 
                                      col = "black",  vline.color = "gray", ci.style = "whisker") + theme_pubr(base_size = 14) + 
                                      theme(plot.title = element_text(face="bold"), plot.margin = unit(c(0.5, 0, 0, 0), "cm")) + labs(title = "A")
plot1

plot2                             <- plot_model(avg_old, type = "std2", order.terms = c(1,2,6,4,3,5), 
                                      axis.labels =  c("age.sample" = "        ", "mrank_birth" = "        ",
                                     "pat.sibs.sample" = "      ","unit.size.cont" = "        ", 
                                     "max.Temp" = "       ", "cum.rain" = "       "),
                                      wrap.labels = F, pred.type = "fe", value.offset = 0.3, value.size = 4.5, 
                                      show.values = T, show.p = T, line.size = 0.2, 
                                      title = "Androgen model (males > 2.5 years)", col = "black",  
                                      vline.color = "gray", line.type = "dotted", ci.style = "whisker") + theme_pubr(base_size = 14) + 
                                      theme(plot.title = element_text(face="bold"),plot.margin = unit(c(0.5, 0, 0, -0.5), "cm")) + labs(title = "B")

plot2

#Plotting effects of maternal rank at birth on androgens for older males    

#changing proportional rank into an ordinal ranks for visualization
older_males$cat.rank          <- ifelse(older_males$mrank_birth >= 0.5, "High", "Low")
older_males$cat.rank          <- factor(older_males$cat.rank, levels = c("High", "Low"))
rank_col <- c("High" = makeTransparent("#117733",200), "Low" = "#44AA99")

plot3 <- ggplot(data = older_males, aes(y = log(ng.g), x = age, col = cat.rank, fill = cat.rank)) + 
  geom_point(alpha = 0.75, size = 2, show.legend = F) + scale_color_manual(values = rank_col) + 
  ylab("log fAMs (ng/g)") + xlab("Male age (years)") + theme_pubr(base_family = "Arial", base_size = 16) + 
  theme(legend.text=element_text(size=16),plot.title = element_text(face="bold", vjust = -3.5), axis.title = element_text(size=16),plot.margin = margin(t = -10,r = 10,b = 5, l = 100)) +
  geom_smooth(alpha = 0.75, size = 0.5, col = "black") + scale_fill_manual("Maternal rank at birth", values = rank_col) + 
  scale_x_continuous(breaks = seq(2.5,8.5,1.5)) + ylim(c(4,8)) + labs(title = "C")
plot3

lay <- rbind(c(1,2),c(3,3))
HB_Fig2 <- grid.arrange(grobs = list(plot1, plot2, plot3), 
                        layout_matrix = lay, widths = c(95,75), 
                        ncol = 2, nrow = 2)
ggsave("HB_Fig2.jpg", HB_Fig2, width = 20, height = 18, units = "cm")

######################....................#################################


#########################################################################
#               Obtaining model coefficient tables                      #
#########################################################################

#Model coefficients for average dredged models from fAM values obtained from young males
avg_young %>% tbl_regression(label = list(age.sample = "Age",  mrank_birth = "Maternal rank",  
                                                pat.sibs.sample = "Cohort size", unit.size.cont = "Unit size",
                                                max.Temp = "Max Temp (\u00B0C)", cum.rain = "Cum rain (mm)")) %>% bold_p() 
coefTable(avg_young)

avg_old %>% tbl_regression(label = list(age.sample = "Age",  mrank_birth = "Maternal rank",  
                                        pat.sibs.sample = "Cohort size", unit.size.cont = "Unit size",
                                        max.Temp = "Max Temp (\u00B0C)", cum.rain = "Cum rain (mm)"))  %>% bold_p() %>% italicize_levels() 
coefTable(avg_old)

######################....................#################################


######################....................################################
#       Androgen models with monthly maternal ranks                     #
######################....................################################

#pre-independence males - using maternal ranks calculated at each month
young_model_birth_lmer_2          <- lmer(log.ng.g ~ age.sample + mrank_prop + pat.sibs.sample +
                                            unit.size.cont + max.Temp + cum.rain + (1|male.ID), 
                                          REML = F, na.action = "na.fail", data = young_males)

summary(young_model_birth_lmer_2)
vif(young_model_birth_lmer_2)

fm1_lmer_young                    <- dredge(young_model_birth_lmer_2, rank = "AICc")
avg_young_monthly_rank                    <- model.avg(get.models(fm1_lmer_young, subset = T), type = 'link', backtransform = TRUE)
summary(avg_young_monthly_rank) #results not different

#post-independence males - using maternal ranks calculated at each month
full_model_older_lmer_2           <- lmer(log.ng.g ~ age.sample  + mrank_prop +  pat.sibs.sample +  
                                            unit.size.cont + cum.rain + max.Temp + (1|male.ID), 
                                          REML = F, na.action = "na.fail",data = older_males)
summary(full_model_older_lmer_2)
vif(full_model_older_lmer_2)

fm1_lmer_old                       <- dredge(full_model_older_lmer_2, rank = "AICc")
avg_old_monthly_rank               <- model.avg(get.models(fm1_lmer_old, subset = T), type = 'link', backtransform = TRUE)
summary(avg_old_monthly_rank) #results not different

#re-plotting model coefficient plots and rank effect plots for older males using monthly maternal ranks

plot4                              <- plot_model(avg_young_monthly_rank, type = "std2", order.terms = c(1,6,2,5,4,3), 
                                                axis.labels =  c("age.sample" = "Age", "mrank_prop" = "Monthly maternal rank",  
                                                "pat.sibs.sample" = "Cohort size", "unit.size.cont" = "Unit size",  
                                                "max.Temp" = "Max Temp (\u00B0C)", "cum.rain" = "Cum rain (mm)"),
                                                wrap.labels = F, pred.type = "fe", value.offset = 0.3, value.size = 4.5, 
                                                show.values = T, show.p = T, line.size = 0.2, title = "Androgen model (males < 2.5 years)", 
                                                col = "black",  vline.color = "gray", ci.style = "whisker") + theme_pubr(base_size = 14) + 
                                                theme(plot.title = element_text(face="bold"), plot.margin = unit(c(0.5, 0, 0, 0), "cm")) + 
                                                labs(title = "A")
plot4

plot5                             <- plot_model(avg_old_monthly_rank, type = "std2", order.terms = c(1,3,6,4,2,5), 
                                                axis.labels =  c("age.sample" = " ", "mrank_prop" = " ",
                                                "pat.sibs.sample" = " ","unit.size.cont" = "  ", 
                                                "max.Temp" = "       ", "cum.rain" = "       "),
                                                wrap.labels = F, pred.type = "fe", value.offset = 0.3, value.size = 4.5, 
                                                show.values = T, show.p = T, line.size = 0.2, 
                                                title = "Androgen model (males > 2.5 years)", col = "black",  
                                                vline.color = "gray", line.type = "dotted", ci.style = "whisker") + 
                                                theme_pubr(base_size = 14) + theme(plot.title = element_text(face="bold"),
                                                plot.margin = unit(c(0.5, 0, 0, -0.5), "cm")) + labs(title = "B")

plot5

#Plotting effects of maternal rank at birth on androgens for older males    

#changing proportional rank into an ordinal ranks for visualization
older_males$cat.rank.monthly          <- ifelse(older_males$mrank_raw_monthly >= 0.5, "High", "Low")
older_males$cat.rank.monthly          <- factor(older_males$cat.rank.monthly, levels = c("High", "Low"))
rank_col                              <- c("High" = makeTransparent("#117733",200), "Low" = "#44AA99")

plot6                                 <- ggplot(data = older_males, aes(y = log(ng.g), x = age, col = cat.rank.monthly, fill = cat.rank.monthly)) + 
                                                geom_point(alpha = 0.75, size = 2, show.legend = F) + scale_color_manual(values = rank_col) + 
                                                ylab("log fAMs (ng/g)") + xlab("Male age (years)") + theme_pubr(base_family = "Arial", base_size = 16) + 
                                                theme(legend.text=element_text(size=16),plot.title = element_text(face="bold", vjust = -3.5), axis.title = element_text(size=16),plot.margin = margin(t = -10,r = 10,b = 5, l = 100)) +
                                                geom_smooth(alpha = 0.75, size = 0.5, col = "black") + scale_fill_manual("Maternal rank each month", values = rank_col) + 
                                                scale_x_continuous(breaks = seq(2.5,8.5,1.5)) + ylim(c(4,8)) + labs(title = "C")

plot6

lay <- rbind(c(1,2),c(3,3))
HB_Supp_Fig8 <- grid.arrange(grobs = list(plot4, plot5, plot6), 
                        layout_matrix = lay, widths = c(95,75), 
                        ncol = 2, nrow = 2)
ggsave("HB_Supp_Fig8.jpg", HB_Supp_Fig8, width = 20, height = 18, units = "cm")

###############################....................###########################


######################....................#################################################
#   predict ng/g values: comparing rank calculated at birth and  each month               #
######################....................#################################################
older_males               <- df_models[df_models$age.sample > cutoff, ]
head(older_males)

#removing all missing data
older_males               <- na.omit(older_males)

#holding all other predictors at their mean
mean_rain                 <- mean(older_males$cum.rain)
mean_min.temp             <- mean(older_males$min.Temp)
mean_max.temp             <- mean(older_males$max.Temp)
mean_pat.sibs             <- mean(older_males$pat.sibs.sample)
mean_unit.size            <- mean(older_males$unit.size.cont)

#set the age increment interval at 0.25 
int                       <- 0.25 #interval by which age increases from 2.5 years onwards 
x                         <- (((9-2.5)/int) + 1)*2 #number of fAM values needed for each age bin



#checking how maternal rank at each interacts with age for T 
int1                      <- lmer(log(ng.g) ~ age.sample + mrank_birth +  pat.sibs.sample +  
                            unit.size.cont + cum.rain + max.Temp + (1|male.ID), REML = F,data = older_males)
summary(int1)

#set the age interval at which you want to predict log ng/g using maternal rank at each month
rank_int_1                <- data.frame(mrank_birth = c(rep(0,x/2),rep(1,x/2)), 
                                age.sample = rep(seq(2.5,9, int),2 ), cum.rain = rep(mean_rain,x), max.Temp = rep(mean_max.temp,x),
                                unit.size.cont = rep(mean_unit.size,x), pat.sibs.sample = rep(mean_pat.sibs,x))

rank_int_1$log.ng.g        <- predict(int1, newdata = rank_int_1, re.form = NA)
rank_int_1$mrank_birth_cat <- ifelse(rank_int_1$mrank_birth <=0.5, "Low", "High")
rank_int_1$mrank_birth_cat <- factor(rank_int_1$mrank_birth_cat, levels = c("High", "Low"))

#boxplots with predicted log ng/g
rank_effect_birth          <- ggplot(data = rank_int_1, aes(x = mrank_birth_cat, y = log.ng.g,  fill = mrank_birth_cat))  + 
                              geom_boxplot(show.legend = F, position = "dodge") + ylab("log fAMs (ng/g)\n") + 
                              xlab("\nMaternal rank (at birth)") + theme_pubr(base_family = "Arial", base_size = 16) + 
                              theme(legend.text=element_text(size=16), axis.title = element_text(size=16)) +
                              scale_color_manual("Maternal rank", values = rank_col)  +
                              scale_fill_manual("Maternal rank", values = rank_col) + ylim (c(5,7)) + stat_compare_means(label = "p.format")
rank_effect_birth

#checking how maternal rank at each month interacts with age for T 
int2                       <- lmer(log(ng.g) ~ age.sample + mrank_prop +  pat.sibs.sample +  
                               unit.size.cont + cum.rain + max.Temp + (1|male.ID), REML = F,data = older_males)
summary(int2)


rank_int2                   <- data.frame(mrank_prop = c(rep(0,x/2),rep(1,x/2)), 
                               age.sample = rep(seq(2.5,9, int),2 ), cum.rain = rep(mean_rain,x), max.Temp = rep(mean_max.temp,x),
                               unit.size.cont = rep(mean_unit.size,x), pat.sibs.sample = rep(mean_pat.sibs,x))
               
rank_int2$log.ng.g          <- predict(int2, newdata = rank_int2, re.form = NA)

rank_int2$mrank_prop_cat    <- ifelse(rank_int2$mrank_prop <=0.5, "Low", "High")
rank_int2$mrank_prop_cat    <- factor(rank_int2$mrank_prop_cat, levels = c("High", "Low"))
rank_col                    <- c("High" = makeTransparent("#117733", 200), "Low" = makeTransparent("#44AA99", 200))

#boxplots with predicted log ng/g - maternal rank at each month of sample collection
rank_effect_monthly         <- ggplot(data = rank_int2, aes(x = mrank_prop_cat, y = log.ng.g,  fill = mrank_prop_cat))  + 
                               geom_boxplot(show.legend = F, position = "dodge") + ylab(" ") + 
                               xlab("\nMaternal rank (monthly)") + theme_pubr(base_family = "Arial", base_size = 16) + 
                               theme(legend.text=element_text(size=16), axis.title = element_text(size=16), axis.text.y=element_blank()) +
                               scale_color_manual("Maternal rank", values = rank_col)  + 
                               scale_fill_manual("Maternal rank", values = rank_col) + ylim (c(5,7)) + stat_compare_means(label = "p.format")
rank_effect_monthly

jpeg("HB_Supp_Fig9.jpg", width = 16, height = 12, units = "cm", res = 340)
ggarrange(rank_effect_birth, rank_effect_monthly,
          labels = c("A", "B"), 
          ncol = 2, nrow = 1)
dev.off()

###############################....................###########################




