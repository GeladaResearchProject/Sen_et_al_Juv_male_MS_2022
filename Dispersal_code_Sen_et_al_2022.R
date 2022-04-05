#########R code for dispersal plots - Sen et al 2022#########

###################### Setting up R ######################
rm(list = ls())
library(dplyr)
library(ggthemes)
library(lubridate)
library(ggfortify)
library(survival)
library(survminer)
library(ggplot2)
library(MuMIn)
library(VIM)
library(sjPlot)
library(tab)
library(gridExtra)

#function for making transparent plots
makeTransparent <- function(black, alpha = my.alpha){
  newColor <- col2rgb(black)
  apply(newColor, 2, function(curcoldata)			
  {rgb(red   = curcoldata[1],
       green = curcoldata[2],
       blue  = curcoldata[3],
       alpha = alpha,
       maxColorValue =  255)})
}

tt <- function(x, t, ...) {
  x=TRUE
  x*t 
}

#####################################################
#   Plotting dispersal per male age cat plot        #
#####################################################

#This code plots the number of males dispersing and number of males censored during each age - Fig 4 in main text 
plot_prop                      <- read.csv("/Users/sharmisen/Desktop/CURRENT MANUSCRIPT PROJECTS/T VALIDATION PAPER/T_MS_Input_Files/plot_dispersal_hist.csv")
plot_prop$ages_factor          <- c("0-1", "1-2", "2-3", "3-4", "4-5", "5-6", "6-7", "7-8","8-9", "9-10")
plot_prop$ages_factor          <- factor(c("0-1", "1-2", "2-3", "3-4", "4-5", "5-6", "6-7", "7-8","8-9", "9-10"), levels = c("0 - 1", "1 - 2", "2 - 3", "3 - 4", "4 - 5", "5 - 6", "6 - 7", "7 - 8","8 - 9", "9 - 10"))

#jpeg("Fig4.jpeg", width = 10, height = 8, units = "in", res = 320)
par(mar = c(5,5,2,2))
xx                             <- barplot(plot_prop$num_dispersing ~ plot_prop$ages, axes = F, col= makeTransparent("#332288",220), xlab = "", ylab = "", ylim = c(0,20), plot = T, cex.names = 1.25)
axis(side = 2, at = seq(0, 20, 5), labels = seq(0, 20, 5), las = 1, cex.axis = 1.25)
mtext(side = 2, text = "Number of males dispersing", cex = 1.5, line = 3.5)
mtext(side = 1, text = "Male ages (years)", cex = 1.5, line = 3)
arrows(x0 = 8.01, y0 = 14, x1 = 8.01, y1 = 18, length = 0.1, angle = 30,
       code = 1, col = par("fg"), lty = par("lty"),
       lwd = 2)
text(x = 7.61, y = 19, labels = "Median age of dispersal = 6.61 years", cex = 1.25)
text(x = xx, y = plot_prop$num_dispersing, label = plot_prop$num_dispersing, pos = 1, cex = 1.25, las = 2,col = "white")
text(x = xx, y = plot_prop$num_dispersing, label = plot_prop$censored_beginning, pos = 3, cex= 1.25, font = 3)
#dev.off()

##################......................#######################

#############################################################
#   Plotting number of dispersal per month across the year #
#############################################################

#This code plots the total number of males dispersing across the year during each month - Fig S7

#strict ages dispersing per month across the year with thirty days range between min and max possible dispersal dates and lose ages dispersing with sixty days range between min and max dates across the year
dispersal_monthly                               <- read.csv("/Users/sharmisen/Desktop/CURRENT MANUSCRIPT PROJECTS/T VALIDATION PAPER/T_MS_Input_Files/dispersal_by_month_plot.csv")
dispersal_monthly$last_emigration_min_date      <- mdy(dispersal_monthly$last_emigration_min_date)
dispersal_monthly$last_emigration_max_date      <- mdy(dispersal_monthly$last_emigration_max_date)
#calculate possible emigration date by assigning middle date between max and min dispersal dates
dispersal_monthly$last_emigration               <- dispersal_monthly$last_emigration_min_date + (dispersal_monthly$last_emigration_max_date -  dispersal_monthly$last_emigration_min_date)/2
range(dispersal_monthly$difference)
dispersal_monthly$month_of_dispersal            <- month(dispersal_monthly$last_emigration)

#write.csv(dispersal_monthly$month_of_dispersal, "dispersal_by_month_plot.csv")
dispersal_monthly_thirty                        <- dispersal_monthly %>% filter(difference <= 30) %>% group_by(month_of_dispersal) %>% summarize(nmales = n())
dispersal_monthly_thirty$month                  <- c("Jan","Feb","Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
dispersal_monthly_sixty                         <- dispersal_monthly %>% 
                                                   filter(difference <= 60) %>% 
                                                   group_by(month_of_dispersal) %>% 
                                                   summarize(nmales = n())

dispersal_monthly_sixty$month                   <- c("Jan","Feb","Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

sum(dispersal_monthly_thirty$nmales)
sum(dispersal_monthly_sixty$nmales)
#jpeg("monthly_dispersal.jpeg", width = 8, height = 6, res = 320, units = "in")
par(mfrow = c(1,2), mar = c(4,2,2,0), oma = c(1,2,0,0))
barplot(dispersal_monthly_thirty$nmales ~ dispersal_monthly_thirty$month, axes = F, plot = T, col= makeTransparent("#332288",220), ylab = " ", xlab = "", ylim = c(0,10), names.arg = c("Jan","Feb","Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), cex.names = 1.2)
axis(side = 2, at = seq(0, 10, 2), labels = seq(0, 10, 2), las = 1, cex.axis = 1.25)
text(x = 12, y = 9, "A", font = 2)
mtext("Number of males dispersing", side = 2, line = 2.5, cex = 1.5)

barplot(dispersal_monthly_sixty$nmales ~ dispersal_monthly_sixty$month, axes = F, plot = T, ylab = " ", xlab = "", ylim = c(0,10), names.arg = c("Jan","Feb","Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), cex.names = 1.2, col = makeTransparent("#888888", 250))
axis(side = 2, at = seq(0, 10, 2), labels = seq(0, 10, 2), las = 1, cex.axis = 1.25)
text(x = 12, y = 9, "B", font = 2)
mtext("2012 - 2020", side = 1, line = 3, cex = 1.2, adj = -0.5)
#dev.off()

##################......................#######################

###############################################################
#         Cox models for modeling dispersal age               #
###############################################################
#reading dispersal dataframe
disp_data            <- read.csv("/Users/sharmisen/Desktop/CURRENT MANUSCRIPT PROJECTS/T VALIDATION PAPER/T_MS_Input_Files/dispersal_dataset_survival_df.csv") #has updated unit size and dead mom rank
disp_data            <- disp_data %>% mutate(unit.size.cat = case_when(unit.size.cont <= 4 ~ "Small", 
                                                           unit.size.cont > 4 & unit.size.cont < 8 ~ "Medium", 
                                                           unit.size.cont >= 8 ~ "Large"))
disp_data            <- disp_data %>% dplyr::select(male_ids, unit, dispersed, no.of.peers.1, no.of.peers, 
                                       older_peers, younger_peers, unit.size.cont, age, ysib, 
                                       mom.presence, mrank_prop, mom.id, unit, orphan.cat)
df                   <- na.omit(disp_data)
cor.test(df$no.of.peers, df$unit.size.cont)


#convert categorical variables into factors
df$ysib              <- as.factor(ifelse(df$ysib < 0,2,df$ysib))
df$ysib              <- factor(df$ysib,labels = c("no sibs", "sibs present","neither"))
df$ysib              <- relevel(df$ysib,ref = "no sibs")
df$mom.presence      <- as.factor(df$mom.presence)
df$mom.presence      <- relevel(df$mom.presence,ref = "Present")
df$male_ids          <- as.factor(df$male_ids)
df$orphan.cat        <- as.factor(ifelse(df$orphan.cat <= 2, "Yes", "No"))

#multivariate model differentiated with all predictors included 
fit_all              <- coxph(Surv(time = df$age, event = df$dispersed) ~ scale(mrank_prop) +
                                   scale(no.of.peers.1) +  ysib + scale(unit.size.cont) +
                                   mom.presence , data = df)
summary(fit_all)

#check violations for proportional hazards assumptions - time independence? 
residuals            <- cox.zph(fit_all)
residuals

#no violation detected

#plotting schoenfeld residuals for each predictor
ggcoxzph(residuals[1], font.main = 14, point.size = 2,
         point.col = makeTransparent("#332288",100), lwd=1.5, ylab="Beta(t) for maternal rank", xlab="Male age (years)")
ggcoxzph(residuals[2], font.main = 14, point.size = 2,
         point.col = makeTransparent("#332288",100), lwd=1.5, ylab="Beta(t) for cohort size", xlab="Male age (years)")
ggcoxzph(residuals[3], font.main = 14, point.size = 2,
         point.col = makeTransparent("#332288",100), lwd=1.5, ylab="Beta(t) for unit size", xlab="Male age (years)")
ggcoxzph(residuals[4], font.main = 14, point.size = 2,
         point.col = makeTransparent("#332288",100), lwd=1.5, ylab="Beta(t) for early sibling arrival", xlab="Male age (years)")
ggcoxzph(residuals[5], font.main = 14, point.size = 2,
         point.col = makeTransparent("#332288",100), lwd=1.5, ylab="Beta(t) for maternal loss", xlab="Male age (years)")

#obtaining coefficients for first full model
coxtable.1           <-  tabcoxph(fit_all, sep.char = "-",`var.labels` = list(c(mrank_prop = "Maternal rank", 
                                                                                    no.of.peers.1 = "Cohort size",
                                                                                    unit.size.cont = "Unit size",
                                                                                    ysib = "Early sibling arrival", 
                                                                                    mom.presence = "Maternal loss = YES"))) 
coxtable.1 

#diagnostic plot for full model
ggcoxdiagnostics(fit_all, type = "dfbetas", linear.predictions = FALSE, ggtheme = theme_pubr(base_size = 12, base_family = "Arial"),
                 robust.se = TRUE, decimals = 2, p.decimals = c(2, 3), p.cuts = 0.01, title = "Diagnostic plot",
                 p.lowerbound = 0.001, p.leading0 = TRUE, p.avoid1 = FALSE, n = FALSE, 
                 events = FALSE)  

#multivariate model with differentiated peer categories and all other predictors included 
fit_all.2             <- coxph(Surv(time = df$age, event = df$dispersed) ~ scale(mrank_prop) + 
                                   scale(older_peers) + scale(younger_peers) + scale(unit.size.cont) + 
                                   ysib + orphan.cat + cluster(male_ids), data = df)
summary(fit_all.2) 
residuals.2 <- cox.zph(fit_all.2)
residuals.2 

#plotting schoenfeld residuals for each predictor
par(mfrow = c(3,2))
ggcoxzph(residuals.2[1], font.main = 14, point.size = 2,
         point.col = makeTransparent("#332288",100), lwd=1.5, ylab="Beta(t) for maternal rank", xlab="Male age (years)")
ggcoxzph(residuals.2[2], font.main = 14, point.size = 2,
         point.col = makeTransparent("#332288",100), lwd=1.5, ylab="Beta(t) for cohort size", xlab="Male age (years)")
ggcoxzph(residuals.2[3], font.main = 14, point.size = 2,
         point.col = makeTransparent("#332288",100), lwd=1.5, ylab="Beta(t) for unit size", xlab="Male age (years)")
ggcoxzph(residuals.2[4], font.main = 14, point.size = 2,
         point.col = makeTransparent("#332288",100), lwd=1.5, ylab="Beta(t) for early sibling arrival", xlab="Male age (years)")
ggcoxzph(residuals.2[5], font.main = 14, point.size = 2,
         point.col = makeTransparent("#332288",100), lwd=1.5, ylab="Beta(t) for maternal loss", xlab="Male age (years)")

#mrank_prop violated assumptions for this model, so adding a time varying mrank_prop predictor to the full model - this model was not reported as a main result
fit_all.2            <- coxph(Surv(time = df$age, event = df$dispersed) ~ scale(mrank_prop) + tt(scale(mrank_prop)) +
                               scale(older_peers) + scale(younger_peers) + 
                                scale(unit.size.cont) + ysib + mom.presence + cluster(male_ids), data = df)
summary(fit_all.2) 

coxtable.2           <-  tabcoxph(fit_all.2, sep.char = "-",`var.labels` = list(c(mrank_prop = "Maternal rank", 
                                                                      no.of.peers.1 = "Cohort size",
                                                                      unit.size.cont = "Unit size",
                                                                      ysib = "Early sibling arrival", 
                                                                      mom.presence = "Maternal loss = YES"))) 
coxtable.2 

##################......................#######################

###############################################################
#         Model coefficient plots for hazards models             #
###############################################################

#first full model
fit_all             <- coxph(Surv(time = df$age, event = df$dispersed) ~ scale(mrank_prop) + 
                                   scale(no.of.peers.1) + scale(unit.size.cont) + ysib + mom.presence, data = df)
summary(fit_all)
pdf("T_MS_Hazards_model1_coef_plot_April4.pdf", height = 5, width = 8)
plot_model(fit_all, wrap.labels = F, rm.terms = c("ysibneither"), 
           axis.labels = c("scale(mrank_prop)" = "Maternal rank", "scale(no.of.peers.1)" = "Cohort size",
                           "scale(unit.size.cont)" = "Unit size", "ysibsibs present" = "Early sibling arrival = YES", "mom.presenceDead" = "Maternal loss = YES"), 
           value.offset = 0.2, value.size = 4.5, transform = "exp", 
           show.values = T, show.p = T, line.size = 0.2, title = "Dispersal model", 
           col = c("black","#332288"),  vline.color = "gray", axis.title =  "Hazards ratio",line.type = "dotted", ci.style = "whisker") + theme_pubr(base_size = 14) + ylim(c(0.001,3))
dev.off()

#second full model with differentiated by peer categories
summary(fit_all.2)
pdf("T_MS_Hazards_model2_coef_plot_with_differentiated_peers_updated_April4.pdf", height = 5, width = 8)
plot_model(fit_all.2, wrap.labels = F, order.terms = c(1,5,6,2,3,4) , rm.terms = c("tt(scale(mrank_prop))", "ysibneither"), 
           axis.labels = c("scale(mrank_prop)" = "Maternal rank", "scale(older_peers)" = "No of older peers", 
                            "scale(younger_peers)" = "No of younger peers",  "scale(unit.size.cont)" = "Unit size", 
                           "ysibsibs present" = "Early sibling arrival = YES", "mom.presenceDead" = "Maternal loss = YES"), 
           value.offset = 0.2, value.size = 4.5, transform = "exp", 
           show.values = T, show.p = T, line.size = 0.2, title = "Dispersal model", 
           col = c("black","#332288"),  vline.color = "gray", axis.title =  "Hazards ratio",line.type = "dotted", ci.style = "whisker") + theme_pubr(base_size = 14) + ylim(c(0.0001,3))
dev.off()

##################......................#######################


###############################################################
#         Generate survival curve based on predictors         #
###############################################################
df1               <- disp_data
str(df1)
#KM curves for maternal rank category
df1               <- df1 %>% mutate(rank.cat = case_when(mrank_prop > 0.5 ~ "High", mrank_prop <= 0.5 ~ "Low"))

fit1              <- survfit(Surv(age, dispersed) ~ factor(rank.cat), data = df1)
summary(fit1)

surv_by_rank      <- ggsurvplot(fit1, data = df1, size = 1.5, xlim = c(4.5,9), 
                           palette = c("#999999","#332288"), pval = T, pval.coord = c(4.5,0.24), conf.int = TRUE,
                           risk.table = F, legend.labs = c("Low","High"), ylab = "Probability of dispersing", 
                           xlab ="Age (years)", break.time.by = 0.5, surv.median.line = "hv", censor = F,
                           ggtheme = theme_pubr(base_size=12, base_family = "Arial"))
surv_by_rank      <- ggpar(surv_by_rank, font.x = c(12, "bold"),
                      font.y = c(12, "bold"),
                      font.caption = c(12, "plain"), 
                      font.legend = c(10, "bold"), 
                      font.tickslab = c(12, "plain"), legend.title = "Maternal rank")
surv_by_rank 

#KM curves for maternal loss 
fit2 <- survfit(Surv(age, dispersed) ~ factor(mom.presence), data = df1)
summary(fit2)

surv_by_maternal_loss <- ggsurvplot(fit2, data = df1, size = 1.5, xlim = c(4.5,9), 
                                    palette = c("#999999","#332288"), pval = T, pval.coord = c(4.5,0.24), conf.int = TRUE,
                                    risk.table = F, legend.labs = c("Yes","No"), ylab = "Probability of dispersing", 
                                    xlab ="Age (years)", break.time.by = 0.5, surv.median.line = "hv", censor = F,
                                    ggtheme = theme_pubr(base_size=12, base_family = "Arial"))
surv_by_maternal_loss <- ggpar(surv_by_maternal_loss, font.x = c(12, "bold"),
                               font.y = c(12, "bold"),
                               font.caption = c(12, "plain"), 
                               font.legend = c(10, "bold"), 
                               font.tickslab = c(12, "plain"), legend.title = "Maternal loss")

surv_by_maternal_loss 

#KM curves for early sibling arrival

df1$ysib              <- as.factor(ifelse(df$ysib < 0,2,df1$ysib))
df1$ysib              <- factor(df1$ysib,labels = c("no sibs", "sibs present","neither"))
df1$ysib              <- relevel(df1$ysib,ref = "no sibs")

fit3 <- survfit(Surv(age, dispersed) ~ ysib, data = df1)
summary(fit3)

surv_by_early_sib     <- ggsurvplot(fit3, data = df1, size = 1.5, xlim = c(4.5,9), 
                                palette = c("#999999","#332288","#E69F00"), pval = T, pval.coord = c(4.5,0.24),conf.int = TRUE,
                                risk.table = F, legend.labs = c("No","Yes", "Neither"), ylab = "Probability of dispersing", 
                                xlab ="Age (years)", break.time.by = 0.5, surv.median.line = "hv", censor = F,
                                ggtheme = theme_pubr(base_size=12, base_family = "Arial"))
surv_by_early_sib     <- ggpar(surv_by_early_sib, font.x = c(12, "bold"),
                           font.y = c(12, "bold"),
                           font.caption = c(12, "plain"), 
                           font.legend = c(10, "bold"), 
                           font.tickslab = c(12, "plain"), legend.title = "Early sibling arrival")
surv_by_early_sib 

#KM curves for peer group size
df1 <- df1 %>% mutate(peers = case_when(no.of.peers.1 <=1 ~ "0-1 male peers", no.of.peers.1 >1 ~ ">1 male peers"))
df1$peers <- factor(df1$peers, levels = c("0-1 male peers", ">1 male peers"))
fit4 <- survfit(Surv(age, dispersed) ~ peers, data = df1)
summary(fit4)


surv_by_peers        <- ggsurvplot(fit4, data = df1, size = 1.5, xlim = c(4.5,9),palette = c("#999999","#332288"), 
                            pval = T, pval.coord = c(4.5,0.24), risk.table = F,conf.int = T, 
                            legend.labs = c("0-1 male peers", ">1 male peers"), censor = F,
                            ylab = "Probability of dispersing", xlab ="Age (years)", 
                            break.time.by = 0.5, surv.median.line = "hv", 
                            ggtheme = theme_pubr(base_size=12, base_family = "Arial")) 

surv_by_peers        <- ggpar(surv_by_peers, font.x = c(12, "bold"),
                       font.y = c(12, "bold"),
                       font.caption = c(12, "plain"), 
                       font.legend = c(10, "bold"), 
                       font.tickslab = c(12, "plain"), legend.title = "Cohort size")
surv_by_peers  

#KM curves for unit size
df1                  <- df1 %>% mutate(unit.size.cat = case_when(unit.size.cont <= 4 ~ "Small", 
                                                              unit.size.cont  > 4 & unit.size.cont < 8 ~ "Medium", 
                                                              unit.size.cont >= 8 ~ "Large"))

fit5 <- survfit(Surv(age, dispersed) ~ unit.size.cat, data = df1)
summary(fit5)

surv_by_unit.size           <- ggsurvplot(fit5, data = df1, size = 1.5, xlim = c(4.5,9), palette = c("#332288","#999999", "#E69F00"), pval = T, 
                                          conf.int = TRUE, risk.table = F, legend.labs = c("Small", "Medium", "Large"), 
                                          pval.coord = c(4.5,0.24),
                                          ylab = "Probability of dispersing", xlab ="Age (years)", 
                                          censor = F,break.time.by = 0.5,
                                           surv.median.line = "hv",  ggtheme = theme_pubr(base_size=12, base_family = "Arial")) 
surv_by_unit.size           <- ggpar(surv_by_unit.size, font.x = c(12, "bold"), 
                                       font.y = c(12, "bold"),
                                       font.caption = c(12, "plain"), 
                                       font.legend = c(10, "bold"), 
                                       font.tickslab = c(12, "plain"), legend.title = "Unit size") 
surv_by_unit.size 

#creating multipanel plots for depicting KM curves based on each kind of predictor (maternal resources and peer resources)
gc()
jpeg("combined_plots_for_maternal_resources.jpeg", width = 12, height = 5, units = "in", res = 340)
ggarrange(surv_by_rank$plot,surv_by_early_sib$plot + rremove("ylab"),
          surv_by_maternal_loss$plot + rremove("ylab") , 
          labels = c("A", "B", "C"), legend = "top",
          ncol = 3, nrow = 1)
dev.off()

gc()
jpeg("combined_plots_for_peer_resources.jpeg", width = 10, height = 5, units = "in", res = 340)
ggarrange(surv_by_peers$plot, surv_by_unit.size$plot+ rremove("ylab"),
          labels = c("A", "B"), legend = "top",
          ncol = 2, nrow = 1)
dev.off()

