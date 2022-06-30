###############     Dispersal plots - Sen et al 2022#########

###################### Setting up R ######################
rm(list = ls())
library(ggthemes)
library(ggfortify)
library(survival)
library(survminer)
library(ggplot2)
library(tab)
library(gridExtra)
library(broom)
library(tidyr)
library(magrittr)


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
#   Number of dispersal per male age               #
#####################################################

#plot the number of males dispersing and number of males censored during each age - Fig 4 in main text 
plot_prop                      <- read.csv("/Users/sharmisen/Desktop/CURRENT MANUSCRIPT PROJECTS/1.T VALIDATION PAPER/T_MS_Input_Files/plot_dispersal_hist.csv")
plot_prop$ages_factor          <- c("0-1", "1-2", "2-3", "3-4", "4-5", "5-6", "6-7", "7-8","8-9", "9-10")
plot_prop$ages_factor          <- factor(c("0-1", "1-2", "2-3", "3-4", "4-5", "5-6", "6-7", "7-8","8-9", "9-10"), levels = c("0 - 1", "1 - 2", "2 - 3", "3 - 4", "4 - 5", "5 - 6", "6 - 7", "7 - 8","8 - 9", "9 - 10"))

jpeg("HB_Fig4.jpeg", width = 16, height = 14, units = "cm", res = 320)
par(mar = c(5,5,2,2))
xx                             <- barplot(plot_prop$num_dispersing ~ plot_prop$ages, axes = F, col= makeTransparent("#332288",250), xlab = "", ylab = "", ylim = c(0,20), plot = T, cex.names = 1.25)
axis(side = 2, at = seq(0, 20, 5), labels = seq(0, 20, 5), las = 1, cex.axis = 1.25)
mtext(side = 2, text = "Number of males dispersing", cex = 1.5, line = 3.5)
mtext(side = 1, text = "Male ages (years)", cex = 1.5, line = 3)
arrows(x0 = 8.01, y0 = 14, x1 = 8.01, y1 = 18, length = 0.1, angle = 30,
       code = 1, col = par("fg"), lty = par("lty"),
       lwd = 2)
text(x = 7.61, y = 19, labels = "Median age of dispersal = 6.61 years", cex = 1.25)
text(x = xx, y = plot_prop$num_dispersing, label = plot_prop$num_dispersing, pos = 1, cex = 1.25, las = 2,col = "white")
text(x = xx, y = plot_prop$num_dispersing, label = plot_prop$censored_beginning, pos = 3, cex= 1.25, font = 3)
dev.off()

##################......................#######################

#############################################################
#   Number of dispersal per month across the year           #
#############################################################

#plot the total number of males dispersing across the year during each month - Fig S7

#strict ages dispersing per month across the year with thirty days range between min and max possible dispersal dates and lose ages dispersing with sixty days range between min and max dates across the year
dispersal_monthly                               <- read.csv("/Users/sharmisen/Desktop/CURRENT MANUSCRIPT PROJECTS/1.T VALIDATION PAPER/T_MS_Input_Files/dispersal_by_month_plot.csv")
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
#number of males with strict dispersal dates - 30 day range
sum(dispersal_monthly_thirty$nmales)
#number of males with lose dispersal dates - 60 day range
sum(dispersal_monthly_sixty$nmales)
jpeg("HB_Supp_Fig7.jpeg", width = 18, height = 16, res = 320, units = "cm")
par(mfrow = c(1,2), mar = c(4,2,2,0), oma = c(1,2,0,0))
barplot(dispersal_monthly_thirty$nmales ~ dispersal_monthly_thirty$month, axes = F, plot = T, col= makeTransparent("#332288",250), ylab = " ", xlab = "", ylim = c(0,10), names.arg = c("Jan","Feb","Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), cex.names = 1.2)
axis(side = 2, at = seq(0, 10, 2), labels = seq(0, 10, 2), las = 1, cex.axis = 1)
text(x = 12, y = 9, "A", font = 2)
mtext("Number of males dispersing", side = 2, line = 2.5, cex = 1.5)

barplot(dispersal_monthly_sixty$nmales ~ dispersal_monthly_sixty$month, axes = F, plot = T, ylab = " ", xlab = "", ylim = c(0,10), names.arg = c("Jan","Feb","Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), cex.names = 1.2, col = makeTransparent("#888888", 250))
axis(side = 2, at = seq(0, 10, 2), labels = seq(0, 10, 2), las = 1, cex.axis = 1)
text(x = 12, y = 9, "B", font = 2)
mtext("2012 - 2020", side = 1, line = 3, cex = 1.2, adj = -0.5)
dev.off()

##################......................#######################


###############################################################
#         Cox models for modeling dispersal age               #
###############################################################
#reading dispersal dataframe

disp_data            <- read.csv("/Users/sharmisen/Desktop/CURRENT MANUSCRIPT PROJECTS/1.T VALIDATION PAPER/T_MS_Input_Files/dispersal_dataset_survival_df.csv") #has updated unit size and dead mom rank
disp_data            <- disp_data %>% dplyr::select(male_ids, unit, dispersed, no.of.peers.1, no.of.peers, 
                                       older_peers, younger_peers, unit.size.cont, age, ysib, 
                                       mom.presence, mrank_prop, mom.id, unit)
df                   <- na.omit(disp_data)

#correlation between number of peers and unit size
cor.test(df$no.of.peers.1, df$unit.size.cont)


#convert categorical variables into factors
df$ysib              <- as.factor(ifelse(df$ysib < 0, 2, df$ysib))
table(df$ysib)
df$ysib              <- factor(df$ysib,labels = c("no sibs", "sibs present","neither"))
table(df$ysib)
df$ysib              <- relevel(df$ysib,ref = "no sibs")
df$mom.presence      <- as.factor(df$mom.presence)
table(df$mom.presence)
df$mom.presence      <- relevel(df$mom.presence,ref = "Present")
df$male_ids          <- as.factor(df$male_ids)


#multivariate model differentiated with all predictors included 
fit_all              <- coxph(Surv(time = df$age, event = df$dispersed) ~ scale(mrank_prop) +
                                   scale(no.of.peers.1) +  ysib + scale(unit.size.cont) +
                                   mom.presence, data = df)
summary(fit_all)

#check violations for proportional hazards assumptions - time independence? 
residuals            <- cox.zph(fit_all)
residuals #no violation detected

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
coxtable.1           <-  tabcoxph(fit_all, sep.char = "-",
                           var.labels = list(mrank_prop = "Maternal rank", no.of.peers.1 = "Cohort size",
                           unit.size.cont = "Unit size",ysib = "Early sibling arrival", mom.presence = "Maternal loss"))
coxtable.1 

#diagnostic plot for full model
ggcoxdiagnostics(fit_all, type = "dfbetas", linear.predictions = FALSE, ggtheme = theme_pubr(base_size = 12, base_family = "Arial"),
                 robust.se = TRUE, decimals = 2, p.decimals = c(2, 3), p.cuts = 0.01, title = "Diagnostic plot",
                 p.lowerbound = 0.001, p.leading0 = TRUE, p.avoid1 = FALSE, n = FALSE, 
                 events = FALSE)  

#multivariate model with differentiated peer categories and all other predictors included 
fit_all.2             <- coxph(Surv(time = df$age, event = df$dispersed) ~ scale(mrank_prop) + 
                                   scale(older_peers) + scale(younger_peers) + scale(unit.size.cont) + 
                                   ysib + mom.presence + cluster(male_ids), data = df)
summary(fit_all.2) 
residuals.2           <- cox.zph(fit_all.2)
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
jpeg("HB_Fig5.jpeg", height = 10, width = 20, units = "cm", res = 340)
plot_model(fit_all, wrap.labels = F, rm.terms = c("ysibneither"), 
           axis.labels = c("scale(mrank_prop)" = "Maternal rank", 
                           "scale(no.of.peers.1)" = "Cohort size",
                           "scale(unit.size.cont)" = "Unit size", "ysibsibs present" = "Early sibling arrival = YES", "mom.presenceDead" = "Maternal loss = YES"), 
            show.values = T, show.p = T,  title = " ", col = c("black","#332288"), value.offset = 0.2, value.size = 4.5, line.size = 0.2,line.type = "dotted", ci.style = "whisker", vline.color = "gray", axis.title =  "Hazards ratio") +
  theme_pubr(base_size = 14) 
dev.off()  

#second full model with differentiated by peer categories
summary(fit_all.2)
jpeg("HB_Supp_FigS13.jpeg", height = 10, width = 20, units = "cm", res = 340)
plot_model(fit_all.2, wrap.labels = F, order.terms = c(1,5,6,2,3,4) , rm.terms = c("tt(scale(mrank_prop))", "ysibneither"), 
           axis.labels = c("scale(mrank_prop)" = "Maternal rank", "scale(older_peers)" = "No of older peers", 
                            "scale(younger_peers)" = "No of younger peers",  "scale(unit.size.cont)" = "Unit size", 
                           "ysibsibs present" = "Early sibling arrival = YES", "mom.presenceDead" = "Maternal loss = YES"), 
           value.offset = 0.2, value.size = 4.5, transform = "exp", 
           show.values = T, show.p = T, line.size = 0.2, title = " ", 
           col = c("black","#332288"),  vline.color = "gray", axis.title =  "Hazards ratio",line.type = "dotted", ci.style = "whisker") + theme_pubr(base_size = 14)
dev.off()

##################......................#######################

###############################################################
#       Predicted median dispersal times                      #
###############################################################
fit.coxph          <- coxph(Surv(time = df$age, event = df$dispersed) ~ mrank_prop +
                     no.of.peers.1 +  ysib + unit.size.cont + mom.presence , data = df)
summary(fit.coxph)
surv_object <- Surv(time = df$age, event = df$dispersed)

#generating predicted median values for maternal rank
newdf             <- df %$% expand.grid(unit.size.cont=6, mrank_prop=c(0,1), no.of.peers.1=0, mom.presence = "Present", ysib = "no sibs")
newdf$group       <- as.character(1:nrow(newdf))

survfit(fit.coxph, newdata = newdf)

# Low rank : 8.65 years median dispersal age
# High rank: 8.89 years median dispersal age

#generating predicted median values for unit size
newdf              <- df %$% expand.grid(unit.size.cont=c(2,6,10), mrank_prop=0.5, no.of.peers.1=0, mom.presence = "Present", ysib = "no sibs")
newdf$group        <- as.character(1:nrow(newdf))

survfit(fit.coxph, newdata = newdf) 

# Small unit  : 8.89 years median dispersal age
# Medium unit :8.65 years median dispersal age
# Large unit  : 8.33 years median dispersal age

#generating predicted median values for number of peers
newdf               <- df %$% expand.grid(unit.size.cont=6, mrank_prop=0.5, no.of.peers.1= c(0,1,2), 
                                          mom.presence = "Present", ysib = "no sibs")
newdf$group         <- as.character(1:nrow(newdf))

survfit(fit.coxph, newdata = newdf) 
# 0 peers: 8.65 years median dispersal age
# 1 peer : 8.35 years median dispersal age
# 2 peers: 8.31 years median dispersal age

#generating predicted median values for maternal loss
newdf              <- df %$% expand.grid(unit.size.cont=6, mrank_prop=0.5, no.of.peers.1= 0, mom.presence = c("Present", "Dead"), ysib = "no sibs")
newdf$group        <- as.character(1:nrow(newdf))

survfit(fit.coxph, newdata = newdf) 

# Mom present: 8.65 years median dispersal age
# Mom dead   : 8.89 years median dispersal age

#generating predicted median values for younger sibling present
newdf           <- df %$% expand.grid(unit.size.cont=6, mrank_prop=0.5, no.of.peers.1= 0, mom.presence = "Present", ysib = c("no sibs", "sibs present"))
newdf$group     <- as.character(1:nrow(newdf))

survfit(fit.coxph, newdata = newdf) 

# Sibs present: 8.65 years median dispersal age
# No sibs     : 8.65 years median dispersal age

###############################....................###########################


###############################################################
#         Generate survival curve based on predictors         #
###############################################################
df1                      <- disp_data

#KM curves for maternal rank category
df1                      <- df1 %>% mutate(rank.cat = case_when(mrank_prop > 0.75 ~ "High", mrank_prop <= 0.75 ~ "Low"))
df1$rank.cat             <- factor(df1$rank.cat, levels = c("High", "Low"))
fit1                     <- survfit(Surv(age, dispersed) ~ rank.cat, data = df1)
summary(fit1)

surv_by_rank             <- ggsurvplot(fit1, data = df1, size = 1.2, xlim = c(4.5,9), legend = "top",
                             palette = c(makeTransparent("#44AA99", 200),"#117733"), 
                             risk.table = F,  ylab = "Proportion of males dispersing", 
                             legend.labs = c("Low rank","High rank"),
                             xlab ="Age (years)", break.time.by = 1.5, censor = F,
                             ggtheme = theme_pubr(base_size=14, base_family = "Arial"))
surv_by_rank             <- ggpar(surv_by_rank, font.x = c(14, "bold"),
                             font.y = c(14, "bold"),
                             font.caption = c(12, "plain"), 
                             font.legend = c(12, "plain"), 
                             font.tickslab = c(14, "plain"), legend.title = " ")
surv_by_rank 


#KM curves for maternal loss 
fit2                     <- survfit(Surv(age, dispersed) ~ factor(mom.presence), data = df1)
summary(fit2)

surv_by_maternal_loss    <- ggsurvplot(fit2, data = df1, size = 1.2, xlim = c(4.5,9), legend = "top",
                                    palette = c("356B4E9",makeTransparent("#0072B2",150)), pval = F,  
                                    risk.table = F, legend.labs = c("Mother dead","Mother alive"), ylab = "Proportion of males dispersing", 
                                    xlab ="Age (years)", break.time.by = 1.5, censor = F,
                                    ggtheme = theme_pubr(base_size = 14, base_family = "Arial"))

surv_by_maternal_loss    <- ggpar(surv_by_maternal_loss, font.x = c(14, "bold"),
                               font.y = c(14, "bold"),
                               font.caption = c(12, "plain"), 
                               font.legend = c(12, "plain"), 
                               font.tickslab = c(14, "plain"), legend.title = " ")

surv_by_maternal_loss 

#KM curves for early sibling arrival
df1$ysib                 <- as.factor(ifelse(df1$ysib < 0,2,df1$ysib))
df1$ysib                 <- factor(df1$ysib,labels = c("no sibs", "sibs present","neither"))
df1$ysib                 <- relevel(df1$ysib,ref = "neither")

fit3 <- survfit(Surv(age, dispersed) ~ ysib, data = df1)
summary(fit3)


surv_by_early_sib        <- ggsurvplot(fit3, data = df1, size = 1.2, xlim = c(4.5,9), 
                                palette = c(makeTransparent("#E69F00",250),"#D55E00","#993404"), pval = F,
                                risk.table = F, legend = "top", legend.labs = c("No younger sib","Younger sib present", "Neither"), 
                                ylab = "Proportion of males dispersing",
                                xlab ="Age (years)", break.time.by = 1.5, censor = F,
                                ggtheme = theme_pubr(base_size=14, base_family = "Arial"))
surv_by_early_sib        <- ggpar(surv_by_early_sib, font.x = c(14, "bold"),
                           font.y = c(14, "bold"),
                           font.caption = c(12, "plain"), 
                           font.legend = c(12, "plain"), 
                           font.tickslab = c(14, "plain"), legend.title = " ")
surv_by_early_sib 



#KM curves for peer group size
df1                     <- disp_data
df1                     <- df1 %>% mutate(peers = case_when(no.of.peers.1 <=1 ~ "0-1 male peers", no.of.peers.1 >1 ~ ">1 male peers"))
df1$peers               <- factor(df1$peers, levels = c("0-1 male peers", ">1 male peers"))
unique(df1$peers)
fit4                    <- survfit(Surv(age, dispersed) ~ peers, data = df1)
summary(fit4)

surv_by_peers           <- ggsurvplot(fit4, data = df1, size = 1.3, xlim = c(4.5,9), 
                                   risk.table = F, palette = c("#999999","#332288"), 
                                   legend.labs = c("0-1 male peers", ">1 male peers"), 
                                   censor = F, xlab ="Age (years)", legend = "top", ylab = "Proportion of males dispersing \n",
                                   break.time.by = 1.5, ggtheme = theme_pubr(base_size=14, base_family = "Arial"))

surv_by_peers           <- ggpar(surv_by_peers, font.x = c(14, "bold"),
                              font.y = c(14, "bold"), 
                              font.caption = c(12, "plain"), 
                              font.legend = c(13, "plain"), 
                              font.tickslab = c(14, "plain"), legend.title = " ") 

surv_by_peers

#KM curves for unit size
df1                     <- df1 %>% mutate(unit.size.cat = case_when(unit.size.cont <= 4 ~ "Small", 
                                                        unit.size.cont  >= 5 & unit.size.cont <= 7 ~ "Medium", 
                                                        unit.size.cont >= 8 ~ "Large"))
df1$unit.size.cat       <- factor(df1$unit.size.cat, levels = c("Small", "Medium","Large"))

fit5                    <- survfit(Surv(age, dispersed) ~ unit.size.cat, data = df1)
summary(fit5)


surv_by_unit.size        <- ggsurvplot(fit5, data = df1, size = 1.3, xlim = c(4.5,9),
                                   palette =  c(makeTransparent("#CC6699",150), makeTransparent("#CC6699", 250), makeTransparent("#660033", 250)), risk.table = F,
                                   censor = F, legend = "top", ylab = "Proportion of males dispersing \n", 
                                   xlab ="Age (years)", legend.labs = c("Small","Medium","Large"),
                                   break.time.by = 1.5, ggtheme = theme_pubr(base_size = 14, base_family = "Arial"))


surv_by_unit.size        <- ggpar(surv_by_unit.size, font.x = c(14, "bold"),
                              font.y = c(14, "bold"), 
                              font.caption = c(12, "plain"), 
                              font.legend = c(13, "plain"), 
                              font.tickslab = c(14, "plain"), legend.title = " ") 
surv_by_unit.size  

#creating multipanel plots for depicting KM curves based on each kind of predictor (maternal resources and peer resources)

#Peer resources and unit size curves
jpeg("HB_Fig6.jpeg", width = 24, height = 12, units = "cm", res = 340)
ggarrange(surv_by_peers$plot, surv_by_unit.size$plot + rremove("ylab") ,
                  labels = c("A", "B"), legend = "top",
                  ncol = 2, nrow = 1)
dev.off()

#Maternal resources curves
jpeg("HB_Supp_FigS12.jpg", width = 40, height = 15, units = "cm", res = 340)
supp_figure12 <- ggarrange(surv_by_rank$plot + rremove("xlab") + rremove("ylab"), 
                           surv_by_early_sib$plot  + rremove("ylab") + rremove("xlab"),
                           surv_by_maternal_loss$plot + rremove("xlab") + rremove("ylab"), hjust = -2,  vjust = 4,
                           labels = c("A", "B", "C"), legend = "top", 
                           ncol = 3 ,nrow = 1)

annotate_figure(supp_figure12, left = textGrob("Proportion of males dispersing", rot = 90, vjust = 0.5, gp = gpar(cex = 1.2, fontface = "bold")),
                bottom = textGrob("Age (years)", gp = gpar(cex = 1.2, fontface = "bold")))
dev.off()

###############################....................###########################
