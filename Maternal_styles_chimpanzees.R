install.packages("lme4")
install.packages("MuMIn")
install.packages("performance")
install.packages("loo")
install.packages("rptR")

library(rptR)
library(loo)
library(performance)
library(lme4)
library(MuMIn)
library(dplyr)
library(car)
library(factoextra)
library(devtools)
library(FactoMineR)
library(corrr)
library(ggcorrplot)
library(RColorBrewer)
library(ggpubr)
library(ggbiplot)
library(psych)
library(stringr)
library(magrittr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(tidyverse)
library(brms)
library(cowplot)
library(parallel)
library(rstan)
library(zoo)
library(scatterplot3d)
library(lme4)
library(arm)
library(MuMIn)
library(plyr)
library(broom)
library(coda)
library(grid)
library(gridExtra)
library(broom.mixed)
library(merTools)
library(tidybayes)
library(HDInterval)

save.image("Maternal_styles_clean_code_directory_final.RData")

load("Maternal_styles_directory_final.RData")



#1. Predictors of maternal styles in chimpanzee mothers ####

# Model 1: Protection - Table 1 and Figure S3 ####

protection.model<-read.csv("protection-model-data-code.csv",sep=";",dec=",")

nrow(protection.model) #number of rows

#factors for model
protection.model$presence.older.sibling  <- factor(protection.model$presence.older.sibling  , levels=c('No','Yes'))
protection.model$Group  <- factor(protection.model$Group  , levels=c('NORD','SUD','EST'))
protection.model$Sex  <- factor(protection.model$Sex  , levels=c('F','M'))
protection.model$danger.type  <- factor(protection.model$danger.type  , levels=c('Focal receiver','Fall from tree','Play'))


#set priors
mprior.protection = get_prior(approach.grab.offspring ~ Age.scale+Age.mother.round.scale+months.pregnant.scale+presence.older.sibling +danger.type+rank.scale+Sex+Group+party.size.scale+(1|Date)+ (1+Age.scale+party.size.scale|Focal), data =  protection.model,
                              family = bernoulli)

mprior.protection
mprior.protection$prior[2:12] <- "normal(0,1)"# so here change the indexing to the rows containing your fixed effects
mprior.protection

#check if model runs
make_stancode(approach.grab.offspring ~ Age.scale+Age.mother.round.scale+months.pregnant.scale+presence.older.sibling +danger.type+rank.scale+Sex+Group+party.size.scale+ (1|Date)+(1+Age.scale+party.size.scale|Focal), data =  protection.model,
              family = bernoulli, prior = mprior.protection)

#run model
b.model.protection <- brm(approach.grab.offspring ~ Age.scale+Age.mother.round.scale+months.pregnant.scale+presence.older.sibling +danger.type+rank.scale+Sex+Group+party.size.scale+(1|Date)+ (1+Age.scale+party.size.scale|Focal), data =  protection.model,
                          family = bernoulli,
                          prior = mprior.protection, 
                          chains = 12, 
                          cores = 12,
                          iter = 2000, 
                          warmup = 1000, 
                          control = list(adapt_delta = 0.88), 
                          sample_prior = "yes")

#look at the summary of the model
summary(b.model.protection)
summary(b.model.protection,prob=0.89)

#extract the random intercepts for the dyads
protection.xx = ranef(b.model.protection ,
                      summary = TRUE)
protection<- as.data.frame(protection.xx$Focal[, , "Intercept"])
protection$Focal <- rownames(protection)
colnames(protection)<-c("protection.Intercept","q1","q2","q3","Focal")

#create dataset with random intercepts
write.csv2(protection,"protection.intercepts-code.csv")

#get the R conditional
r2_bayes(b.model.protection) 

#effects of each predictor
conditional_effects(b.model.protection)



# Model 2: Rejection - Table 2 and Figure S4 ####

rejection.model<-read.csv("rejection-model-data-code.csv",sep=";",dec=",")

#number of rows
nrow(rejection.model)

#factors for model
rejection.model$presence.older.sibling  <- factor(rejection.model$presence.older.sibling  , levels=c('Yes','No'))
rejection.model$Group  <- factor(rejection.model$Group  , levels=c('NORD', 'SUD','EST'))
rejection.model$Sex  <- factor(rejection.model$Sex  , levels=c('F', 'M'))
rejection.model$ACTIVITE.MERE  <- factor(rejection.model$ACTIVITE.MERE  , levels=c('RP','deplacement','MG','Epouillage','Jeu'))

#check correlation
v2<-lm(mother.response ~ Age.scale + I(Age.scale^2)+Age.mother.round.scale+months.pregnant.scale+rank.scale+presence.older.sibling + Sex+ Group + party.size.scale + ACTIVITE.MERE,data=rejection.model)
vif(v2)


#set priors
mprior.rejection = get_prior(mother.response ~ Age.scale + I(Age.scale^2)+Age.mother.round.scale+months.pregnant.scale+rank.scale+presence.older.sibling + Sex+ Group + party.size.scale + ACTIVITE.MERE+(1+ACTIVITE.MERE+party.size.scale|Date)+ (1+ACTIVITE.MERE+party.size.scale +Age.scale|Focal), data =  rejection.model,
                             family = bernoulli)

mprior.rejection
mprior.rejection$prior[2:15] <- "normal(0,1)"# so here change the indexing to the rows containing your fixed effects
mprior.rejection

#check if model runs
make_stancode(mother.response ~ Age.scale + I(Age.scale^2)+Age.mother.round.scale+months.pregnant.scale+rank.scale+presence.older.sibling + Sex+ Group + party.size.scale + ACTIVITE.MERE+(1+ACTIVITE.MERE+party.size.scale|Date)+ (1+ACTIVITE.MERE+party.size.scale +Age.scale|Focal), data =  rejection.model,
              family = bernoulli, prior = mprior.rejection)

#run model
b.model.rejection <- brm(mother.response ~ Age.scale + I(Age.scale^2)+Age.mother.round.scale+months.pregnant.scale+rank.scale+presence.older.sibling + Sex+ Group + party.size.scale + ACTIVITE.MERE+(1+ACTIVITE.MERE+party.size.scale|Date)+ (1+ACTIVITE.MERE+party.size.scale +Age.scale|Focal), data =  rejection.model,
                         family = bernoulli,
                         prior = mprior.rejection, 
                         chains = 12, 
                         cores = 12,
                         iter = 2000, 
                         warmup = 1000, 
                         thin = 2, #prevent autocorrelation, removing 1 iteration every 3 iterations
                         control = list(adapt_delta = 0.98), 
                         sample_prior = "yes")

#get the summary of the model
summary(b.model.rejection)
summary(b.model.rejection,prob=0.89)

#get random intercepts for dyads
rejection.xx = ranef(b.model.rejection ,
                     summary = TRUE)
rejection <- as.data.frame(rejection.xx$Focal[, , "Intercept"])
rejection$Focal <- rownames(rejection)
colnames(rejection )<-c("rejection.Intercept","q1","q2","q3","Focal")

write.csv2(rejection,"rejection-intercepts-code.csv")

#check effects of predictors
conditional_effects(b.model.rejection)

#get conditional R
r2_bayes(b.model.rejection)



# Model 3: Playing - Table 3 ####

play.model<-read.csv("play.model-data-code.csv",sep=";",dec=",")


#number of rows
nrow(play.model)

#factors for model
play.model$Group  <- factor(play.model$Group  , levels=c('NORD', 'SUD','EST'))
play.model$Sex  <- factor(play.model$Sex  , levels=c('F', 'M'))
play.model$family.member  <- factor(play.model$family.member  , levels=c('Yes', 'No'))


play.model$Group.centered  <-as.numeric(as.factor(play.model$Group))-mean(as.numeric(as.factor(play.model$Group)))
play.model$Sex.centered  <-as.numeric(as.factor(play.model$Sex))-mean(as.numeric(as.factor(play.model$Sex)))
play.model$family.member.centered  <-as.numeric(as.factor(play.model$family.member))-mean(as.numeric(as.factor(play.model$family.member)))


play.model$Duration.static<-as.numeric(play.model$Duration.static)
play.model$duration.play<-as.numeric(play.model$duration.play)

#test collinearity
v4<-lm(duration.play ~ Age.scale+ Age.mother.round.scale+months.pregnant.scale+rank.scale +Sex+family.member+Group+mean.party.size.scale,data=play.model)
vif(v4)

#set priors
mprior.play = get_prior(duration.play ~ Age.scale+I(Age.scale^2)+offset(log(Duration.static))+family.member+months.pregnant.scale  +Age.mother.round.scale+  rank.scale +Group+Sex+mean.party.size.scale + (1+Age.scale+mean.party.size.scale|Focal), data =  play.model,
                        family = negbinomial())

mprior.play
mprior.play$prior[2:11] <- "normal(0,1)"# so here change the indexing to the rows containing your fixed effects
mprior.play

make_stancode(duration.play ~ Age.scale+I(Age.scale^2)+offset(log(Duration.static))+family.member+months.pregnant.scale  +Age.mother.round.scale+  rank.scale +Group+Sex+mean.party.size.scale + (1+Age.scale+mean.party.size.scale|Focal), data =  play.model,
              family = negbinomial(), prior = mprior.play)

b.model.play <- brm(duration.play ~ Age.scale+I(Age.scale^2)+offset(log(Duration.static))+family.member+months.pregnant.scale  +Age.mother.round.scale+  rank.scale +Group+Sex+mean.party.size.scale + (1+Age.scale+mean.party.size.scale|Focal), data =  play.model,
                    family = negbinomial(),
                    prior = mprior.play, 
                    chains = 12, 
                    cores = 12,
                    iter = 2000, 
                    warmup = 1000, 
                    control = list(adapt_delta = 0.99),
                    sample_prior = "yes")


#centered
mprior.play = get_prior(duration.play ~ Age.scale+I(Age.scale^2)+offset(log(Duration.static))+family.member.centered+months.pregnant.scale  +Age.mother.round.scale+  rank.scale +Group.centered+Sex.centered+mean.party.size.scale + (1+Age.scale+mean.party.size.scale|Focal), data =  play.model,
                        family = negbinomial())

mprior.play
mprior.play$prior[2:10] <- "normal(0,1)"# so here change the indexing to the rows containing your fixed effects
mprior.play

make_stancode(duration.play ~ Age.scale+I(Age.scale^2)+offset(log(Duration.static))+family.member.centered+months.pregnant.scale  +Age.mother.round.scale+  rank.scale +Group.centered+Sex.centered+mean.party.size.scale + (1+Age.scale+mean.party.size.scale|Focal), data =  play.model,
              family = negbinomial(), prior = mprior.play)

b.model.play.centered <- brm(duration.play ~ Age.scale+I(Age.scale^2)+offset(log(Duration.static))+family.member.centered+months.pregnant.scale  +Age.mother.round.scale+  rank.scale +Group.centered+Sex.centered+mean.party.size.scale + (1+Age.scale+mean.party.size.scale|Focal), data =  play.model,
                             family = negbinomial(),
                             prior = mprior.play, 
                             chains = 12, 
                             cores = 12,
                             iter = 2000, 
                             warmup = 1000, 
                             control = list(adapt_delta = 0.99),
                             sample_prior = "yes")



#get random intercepts for dyads
play.xx= ranef(b.model.play.centered ,
               summary = TRUE)
play <- as.data.frame(play.xx$Focal[, , "Intercept"])
play$Focal <- rownames(play)
colnames(play)<-c("play.Intercept","q1","q2","q3","Focal")

write.csv2(play,"play-intercepts2-centered.csv")

#summary of model
summary(b.model.play.centered)
summary(b.model.play.centered, prob=0.89)



#get random intercepts for dyads
play.xx= ranef(b.model.play ,
               summary = TRUE)
play <- as.data.frame(play.xx$Focal[, , "Intercept"])
play$Focal <- rownames(play)
colnames(play)<-c("play.Intercept","q1","q2","q3","Focal")

write.csv2(play,"play-intercepts-code.csv")

#summary of model
summary(b.model.play)
summary(b.model.play, prob=0.89)

#effects of predictors
conditional_effects(b.model.play)


#get R conditional
r2_bayes(b.model.play)



# Model 4: Offspring groomed by mother - Table 4 ####

groomed.model<-read.csv("groomed-model-data-code.csv",sep=";",dec=",")


#check number of rows
nrow(groomed.model)

#factors for model
groomed.model$Group  <- factor(groomed.model$Group  , levels=c('NORD', 'SUD','EST'))
groomed.model$Sex  <- factor(groomed.model$Sex  , levels=c('F', 'M'))
groomed.model$family.member  <- factor(groomed.model$family.member  , levels=c('Yes', 'No'))

#test collinearity
v4<-lm(duration.groomed ~ Age.scale+ Age.mother.round.scale+months.pregnant.scale+rank.scale +Sex+family.member+Group+mean.party.size.scale,data=groomed.model)
vif(v4)


#setpriors
mprior.groomed = get_prior(duration.groomed  ~ Age.scale+I(Age.scale^2)+months.pregnant.scale+ Age.mother.round.scale+offset(log(Duration.static))+rank.scale+Sex+family.member+Group+mean.party.size.scale + (1+Age.scale+mean.party.size.scale|Focal), data =  groomed.model,family = negbinomial())

mprior.groomed
mprior.groomed$prior[2:11] <- "normal(0,1)"# so here change the indexing to the rows containing your fixed effects
mprior.groomed

#check if code runs
make_stancode(duration.groomed  ~ Age.scale+I(Age.scale^2)+months.pregnant.scale+ Age.mother.round.scale+offset(log(Duration.static))+rank.scale+Sex+family.member+Group+mean.party.size.scale + (1+Age.scale+mean.party.size.scale|Focal), data =  groomed.model,
              family = negbinomial(), prior = mprior.groomed)

#run model
b.groomed.model<- brm(duration.groomed ~ Age.scale + I(Age.scale^2) + months.pregnant.scale + Age.mother.round.scale + offset(log(Duration.static)) + rank.scale + Sex + family.member + Group + mean.party.size.scale + (1 + Age.scale + mean.party.size.scale | Focal), 
                      data = groomed.model, 
                      family = negbinomial(), 
                      prior = mprior.groomed, 
                      chains = 12, 
                      cores = 12, 
                      iter = 2000, 
                      warmup = 1000, 
                      control = list(adapt_delta = 0.8), 
                      sample_prior = "yes")


#get random intercepts for dyads
groomed.xx = ranef(b.groomed.model ,
                   summary = TRUE)
groomed <- as.data.frame(groomed.xx$Focal[, , "Intercept"])
groomed $Focal <- rownames(groomed)
colnames(groomed )<-c("groomed.Intercept","q1","q2","q3","Focal")

write.csv2(groomed,"groomed-intercepts-code.csv")

#summary of the model
summary(b.groomed.model)
summary(b.groomed.model,prob=0.89)


#effects of predictors
conditional_effects(b.groomed.model)

#R conditional
r2_bayes(b.groomed.model)




# Model 5: Proximity - Table 5 ####

proximity.model<-read.csv("proximity-model-data-code.csv",sep=";",dec=",")

#number of rows
nrow(proximity.model)

#factors for model
proximity.model$Group  <- factor(proximity.model$Group  , levels=c('NORD', 'SUD','EST'))
proximity.model$Sex  <- factor(proximity.model$Sex  , levels=c('F', 'M'))
proximity.model$family.member  <- factor(proximity.model$family.member  , levels=c('Yes', 'No'))

#check collinearity
v5<-lm(contact.one.meter ~ Age.scale +Age.mother.round.scale+rank.scale+ months.pregnant.scale+Sex+Group+mean.party.size.scale, data=proximity.model)
vif(v5)

#set priors
mprior.proximity = get_prior(contact.one.meter ~ Age.scale+Sex+months.pregnant.scale+offset(log(total.scans))+family.member+Group+Age.mother.round.scale+rank.scale+mean.party.size.scale   +(1+mean.party.size.scale+Age.scale|Focal), data =  proximity.model,
                             family = negbinomial())

mprior.proximity
mprior.proximity$prior[2:10] <- "normal(0,1)"# so here change the indexing to the rows containing your fixed effects
mprior.proximity

#check if model runs
make_stancode(contact.one.meter ~ Age.scale+Sex+family.member+months.pregnant.scale+offset(log(total.scans))+Group+Age.mother.round.scale+rank.scale+mean.party.size.scale   +(1+mean.party.size.scale+Age.scale|Focal), data =  proximity.model,
              family = negbinomial(), prior = mprior.proximity)

#run model
b.model.proximity <- brm(contact.one.meter ~ Age.scale+Sex+family.member+months.pregnant.scale+offset(log(total.scans))+Group+Age.mother.round.scale+rank.scale+mean.party.size.scale  +(1+mean.party.size.scale+Age.scale|Focal), data =  proximity.model,
                         family = negbinomial(),
                         prior = mprior.proximity, 
                         chains = 12, 
                         cores = 12,
                         iter = 2000, 
                         warmup = 1000, 
                         control = list(adapt_delta = 0.98), 
                         sample_prior = "yes")


#get random intercepts for dyads
proximity.xx= ranef(b.model.proximity ,
                    summary = TRUE)
proximity <- as.data.frame(proximity.xx$Focal[, , "Intercept"])
proximity$Focal <- rownames(proximity)
colnames(proximity)<-c("contact.Intercept","q1","q2","q3","Focal")

write.csv2(proximity,"proximity-intercepts-code.csv")

#effects of predictors
conditional_effects(b.model.proximity)

#check summary
summary(b.model.proximity)
summary(b.model.proximity,prob=0.89)

#check R conditional
r2_bayes(b.model.proximity)



#Figure 1: Variation of maternal care with offspring age ####

#put together
danger.model<-read.csv("protection-model-data.csv",sep=";",dec=",")
refusal.model<-read.csv("rejection-model-data.csv",sep=";",dec=",")
play.model<-read.csv("playing-model-data.csv",sep=";",dec=",")
groomed.model<-read.csv("groomed-model-data.csv",sep=";",dec=",")
contact.model<-read.csv("proximity-model-data.csv",sep=";",dec=",")

play.model<-read.csv("play.model-data2.csv",sep=";",dec=",")

play.model$Age.new

nrow(refusal.model)
nrow(danger.model)
nrow(play.model)
nrow(groomed.model)
nrow(contact.model)

protection.intercept<-read.csv("protection.intercepts.csv",sep=";",dec=",")
rejection.intercept<-read.csv("rejection-intercepts.csv",sep=";",dec=",")
play.intercept<-read.csv("play-intercepts2.csv",sep=";",dec=",")
groomed.intercept<-read.csv("groomed-intercepts.csv",sep=";",dec=",")
proximity.intercept<-read.csv("proximity-intercepts.csv",sep=";",dec=",")


rejection.protection<-merge(rejection.intercept,protection.intercept,by="Focal",all=T)
rejection.protection.play<-merge(rejection.protection,play.intercept,by="Focal",all=T)
rejection.protection.play.groomed<-merge(rejection.protection.play,groomed.intercept,by="Focal",all=T)
behaviours.intercept<-merge(rejection.protection.play.groomed,proximity.intercept,by="Focal",all=T)


behaviours.intercept<-data.frame(behaviours.intercept$Focal,behaviours.intercept$refusal.Intercept,behaviours.intercept$protection.Intercept,behaviours.intercept$play.Intercept,behaviours.intercept$groomed.Intercept,behaviours.intercept$contact.Intercept)
colnames(behaviours.intercept)<-c("Focal","rejection.Intercept","protection.Intercept","play.Intercept","groomed.Intercept","proximity.Intercept")

write.csv2(behaviours.intercept,"behaviours-intercept-maternal-styles.csv")


behaviours.intercept<-read.csv("behaviours-intercept-maternal-styles.csv",sep=";",dec=",")


#make the figure

groomed.model$mean_groomed<- groomed.model$duration.groomed/groomed.model$Duration.static
play.model$mean_play <- play.model$duration.play/play.model$Duration.static
proximity.model$mean_contact <- proximity.model$contact.one.meter/proximity.model$total.scans


# Aggregate the mean of 'mother.response' by 'Date' and 'Focal' and keep 'Age'
mean_refusal <- aggregate(
  mother.response ~ Date + Focal, 
  data = refusal.model, 
  FUN = function(x) mean(x, na.rm = TRUE)
)

# Now, calculate the mean of 'Age.new' by 'Date' and 'Focal'
mean_age_refusal <- aggregate(
  Age.new ~ Date + Focal, 
  data = refusal.model, 
  FUN = function(x) mean(x, na.rm = TRUE)
)

# Merge the mean of 'mother.response' and 'Age.new' data
mean_refusal_final <- merge(mean_refusal, mean_age_refusal, by = c("Date", "Focal"))

# View the result
print(mean_refusal_final)



# Aggregate the mean of 'approach.grab.offspring' by 'Date' and 'Focal' and keep 'Age'
mean_danger <- aggregate(
  approach.grab.offspring ~ Date + Focal, 
  data = danger.model, 
  FUN = function(x) mean(x, na.rm = TRUE)
)

# Now, calculate the mean of 'Age.new' by 'Date' and 'Focal'
mean_age_danger <- aggregate(
  Age.new ~ Date + Focal, 
  data = danger.model, 
  FUN = function(x) mean(x, na.rm = TRUE)
)

# Merge the mean of 'approach.grab.offspring' and 'Age.new' data
mean_danger_final <- merge(mean_danger, mean_age_danger, by = c("Date", "Focal"))

# View the result
print(mean_danger_final)


# Assuming you have the means and standard deviations for scaling:
mean_age <- mean(refusal.model$Age.new)  # Replace no your actual data
sd_age <- sd(refusal.model$Age.new)

# Create a vector for unscaled age based on the scaled values
age_unscaled <- seq(from = min(refusal.model$Age.new), to = max(refusal.model$Age.new), length.out = 100)


#mean number of hours per focal per day

play.model$Duration.hours<-play.model$Duration.static/3600
groomed.model$Duration.hours<-groomed.model$Duration.static/3600

mean.hours.play <- sum(play.model$Duration.hours)/nrow(play.model)
mean.hours.groomed <- sum(groomed.model$Duration.hours)/nrow(groomed.model)
mean.scan.contact <- sum(contact.model$total.scans)/nrow(contact.model)


# Calculate mean and standard deviation for unscaling


# Add unscaled age to the conditional effects for each model
contact_effects <- conditional_effects(b.model.proximity)$`Age.scale`
contact_effects$age_unscaled <- contact_effects$Age.scale * sd_age + mean_age

play_effects <- conditional_effects(b.model.play)$`Age.scale`
play_effects$age_unscaled <- play_effects$Age.scale * sd_age + mean_age

play_effects <- conditional_effects(b.model.play.centered)$`Age.scale`
play_effects$age_unscaled <- play_effects$Age.scale * sd_age + mean_age


groomed_effects <- conditional_effects(b.groomed.model)$`Age.scale`
groomed_effects$age_unscaled <- groomed_effects$Age.scale * sd_age + mean_age


mean_age <- mean(danger.model$Age.new)  # Replace no actual data
sd_age <- sd(danger.model$Age.new)

danger_effects <- conditional_effects(b.model.protection)$`Age.scale`
danger_effects$age_unscaled <- danger_effects$Age.scale * sd_age + mean_age

mean_age <- mean(refusal.model$Age.new)  # Replace no actual data
sd_age <- sd(refusal.model$Age.new)

refusal_effects <- conditional_effects(b.model.rejection)$`Age.scale`
refusal_effects$age_unscaled <- refusal_effects$Age.scale * sd_age + mean_age

# Create individual plots

unscaled_ages <- c(0, 60, 120)
scaled_breaks <- (unscaled_ages - mean(danger.model$Age.new)) / sd(danger.model$Age.new)


danger_plot <- ggplot(danger_effects, aes(x = age_unscaled, y = estimate__)) +
  geom_line(size = 1.2, color = "#e31a1c") +  # Red color for the model line
  scale_x_continuous(
    name = "Offspring age (in months)",
    limits = c(0, 120),
    breaks = unscaled_ages,
    labels = as.character(unscaled_ages)
  ) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#e31a1c") + 
  # Add raw data points for each Date and Focal (dots in black, smaller size, no legend for Focal)
  geom_point(data = mean_danger_final, aes(x = Age.new, y = approach.grab.offspring), 
             color = "#e31a1c", size = 0.5, alpha = 0.5, shape = 16, show.legend = FALSE) + 
  labs(title = "(a) Protection", x = "Offspring age (in months)", y = "Probability of mother approaching the\n offspring during a threat") +
  theme_minimal()+
  scale_y_continuous(
    name = "Probability of mother approaching the\n offspring during a threat",
    breaks = seq(0, 1, by = 0.2),  # Set breaks at 0, 0.2, ..., 1
    limits = c(0, 1),  # Ensure y-axis limits are from 0 to 1
    labels = function(x) round(x, 2)  # Round y-axis labels to 2 decimal places
  ) 



# Refusal Plot (Refusal of solicitations)
unscaled_ages <- c(0, 60, 120)
scaled_breaks <- (unscaled_ages - mean(refusal.model$Age.new)) / sd(refusal.model$Age.new)

refusal_plot <- ggplot(refusal_effects, aes(x = age_unscaled, y = estimate__)) +
  geom_line(size = 1.2, color = "#6a3d9a") +  # Purple color for the model line
  scale_x_continuous(
    name = "Offspring age (in months)",
    limits = c(0, 120),
    breaks = unscaled_ages,
    labels = as.character(unscaled_ages)
  ) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#6a3d9a") +
  # Add raw data points for each Date and Focal (dots in black, smaller size, no legend for Focal)
  geom_point(data = mean_refusal_final, aes(x = Age.new, y = mother.response), 
             color = "#6a3d9a", size = 0.5, alpha = 0.5, shape = 16, show.legend = FALSE) +
  labs(title = "(b) Refusal of solicitations", x = "Offspring age (in months)", y = "Probability of maternal refusal of\n offspring solicitation") +
  theme_minimal()+
  scale_y_continuous(
    name = "Probability of maternal refusal of\n offspring solicitation",
    breaks = seq(0, 1, by = 0.2),  # Set breaks at 0, 0.2, ..., 1
    limits = c(0, 1),  # Ensure y-axis limits are from 0 to 1
    labels = function(x) round(x, 2)  # Round y-axis labels to 2 decimal places
  ) 




#play

# Calculate Aggression Rate per Hour per Dyad
columns_to_select <- c("Focal", "Age.scale", "duration.play", "Duration.static")
zdata <- play.model[, columns_to_select]
b = aggregate(zdata$duration.play / zdata$Duration.static, by = list(zdata$Age.scale, zdata$Focal), mean)
a = aggregate(zdata$Duration.static, by = list(zdata$Age.scale, zdata$Focal), FUN = sum)
b$number = a$x
colnames(b) = c("Age.scale", "Focal", "proba", "number")

# Plot the Aggression Rate (in terms of hourly peering rate)
d.plot <- conditional_effects(b.model.play.centered, effects = "Age.scale", re_formula = NA,
                              robust = TRUE, prob = 0.95, method = "fitted",
                              spaghetti = FALSE, surface = FALSE, resolution = 100)


# Calculate dynamic breaks for scaled x-axis
unscaled_ages <- c(0, 60, 120)
scaled_breaks <- (unscaled_ages - mean(play.model$Age.new)) / sd(play.model$Age.new)


play_plot <- plot(d.plot, plot = FALSE)[[1]] +
  theme_minimal() +
  scale_x_continuous(
    name = "Offspring age (in months)",
    limits = NULL,
    breaks = scaled_breaks,
    labels = as.character(unscaled_ages)
  ) +
  labs(
    y = "Probability of mother-offspring\n playing per hour",
    title = "(c) Playing",
    x = "Offspring age (in months)"
  ) +
  scale_y_continuous(
    name = "Probability of mother-offspring\n playing per hour",
    limits = c(0, 0.2*3600*mean.hours.play),  # Set y-axis limits
    breaks = seq(0, 0.2*3600*mean.hours.play, by = 0.1*3600*mean.hours.play),  # Custom breaks
    labels = function(x) round(x / (3600 * mean.hours.play), 1)  # Scale and round the labels
  ) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#33a02c") +  
  geom_point(data = b, mapping = aes(x = Age.scale, y = proba * mean.hours.play * 3600), 
             alpha = 0.5, inherit.aes = FALSE, color = "#33a02c", size = 0.5, shape = 16) +
  geom_line(size = 1.2, color = "#33a02c") 

# Display the plot
print(play_plot)



#groomed

# Calculate grooming Rate per Hour per Dyad
columns_to_select <- c("Focal", "Age.scale", "duration.groomed", "Duration.static")
zdata <- groomed.model[, columns_to_select]

# grooming rate per dyad
b = aggregate(zdata$duration.groomed / zdata$Duration.static, by = list(zdata$Age.scale, zdata$Focal), mean)
a = aggregate(zdata$Duration.static, by = list(zdata$Age.scale, zdata$Focal), FUN = sum)
b$number = a$x
colnames(b) = c("Age.scale", "Focal", "proba", "number")

# Plot the Aggression Rate (in terms of hourly peering rate)
d.plot.groomed <- conditional_effects(b.groomed.model.no1800, effects = "Age.scale", re_formula = NA,
                                      robust = TRUE, prob = 0.95, method = "fitted",
                                      spaghetti = FALSE, surface = FALSE, resolution = 100)

unscaled_ages <- c(0, 60, 120)
scaled_breaks <- (unscaled_ages - mean(groomed.model$Age.new)) / sd(groomed.model$Age.new)

# Creating the plot
groomed_plot <- plot(d.plot.groomed, plot = FALSE)[[1]] +
  theme_minimal() +
  scale_x_continuous(
    name = "Offspring age (in months)",
    limits = NULL,
    breaks = scaled_breaks,
    labels = as.character(unscaled_ages)
  ) +
  labs(y = "Probability of mother grooming\n her offspring per hour",title = "(d) Mother grooming \nher offspring", x = "Offspring age (in months)") +  # Title for Y-axis
  scale_y_continuous(limits = c(0, 0.2*3600*mean.hours.groomed),  # Set y-axis limits
                     breaks = seq(0, 0.2*3600*mean.hours.groomed, by = 0.1*3600*mean.hours.groomed),
                     labels = function(x) round(x / (mean.hours.groomed * 3600), 1)  # Round to 2 decimal places
  ) + 
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#ff7f00") + # Conversion factor (from seconds to hours)
  geom_point(data = b, mapping = aes(x = Age.scale, y = proba * mean.hours.groomed * 3600), 
             alpha = 0.5, inherit.aes = FALSE, color = "#ff7f00", size = 0.5) +  # Change color of data points
  geom_line(size = 1.2, color = "#ff7f00") 

print(groomed_plot)


#contact

# Calculate grooming Rate per Hour per Dyad
columns_to_select <- c("Focal", "Age.scale", "contact.one.meter", "total.scans")
zdata <- contact.model[, columns_to_select]


b = aggregate(zdata$contact.one.meter / zdata$total.scans, by = list(zdata$Age.scale, zdata$Focal), mean)
a = aggregate(zdata$total.scans, by = list(zdata$Age.scale, zdata$Focal), FUN = sum)
b$number = a$x
colnames(b) = c("Age.scale", "Focal", "proba", "number")


d.plot.contact <- conditional_effects(b.contact.no10, effects = "Age.scale", re_formula = NA,
                                      robust = TRUE, prob = 0.95, method = "fitted",
                                      spaghetti = FALSE, surface = FALSE, resolution = 100)

unscaled_ages <- c(0, 60, 120)
scaled_breaks <- (unscaled_ages - mean(contact.model$Age.new)) / sd(contact.model$Age.new)


# Creating the plot
contact_plot <- plot(d.plot.contact, plot = FALSE)[[1]] +
  theme_minimal() +
  scale_x_continuous(
    name = "Offspring age (in months)",
    limits = NULL,
    breaks = scaled_breaks,
    labels = as.character(unscaled_ages)
  ) +
  labs(y = "Probability of mother and offspring\n noin 1m per focal",title = "(e) Proximity", x = "Offspring age (in months)") +  # Title for Y-axis
  scale_y_continuous(limits = c(0, 1.1*mean.scan.contact),  # Set y-axis limits
                     breaks = seq(0, 1.1*mean.scan.contact, by = 0.2*mean.scan.contact),
                     labels = function(x) round(x / (mean.scan.contact), 1)  # Round to 2 decimal places
  ) + 
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, fill = "#1f78b4") + # Conversion factor (from seconds to hours)
  geom_point(data = b, mapping = aes(x = Age.scale, y = proba * mean.scan.contact), 
             alpha = 0.5, inherit.aes = FALSE, color = "#1f78b4", size = 0.5) +  # Change color of data points
  geom_line(size = 1.2, color = "#1f78b4") 

print(contact_plot)


# Combine all plots into one figure
combined_plot <- (danger_plot | refusal_plot) /
  (play_plot | groomed_plot | contact_plot) + 
  plot_layout(guides = 'collect') & 
  theme(legend.position = "bottom")


# Display the combined plot
print(combined_plot)






# Principal Component Analysis: Table 3 and Figure S1 ####

mother.PCA<-read.csv("behaviours-intercept-maternal-styles.csv",sep=';',dec=",")

names(mother.PCA)[names(mother.PCA) == "rejection.Intercept"] <- "Rejection"                      
names(mother.PCA)[names(mother.PCA) == "protection.Intercept"] <- "Protection" 
names(mother.PCA)[names(mother.PCA) == "groomed.Intercept"] <- "Grooming" 
names(mother.PCA)[names(mother.PCA) == "play.Intercept"] <- "Playing" 
names(mother.PCA)[names(mother.PCA) == "proximity.Intercept"] <- "Proximity" 


mother.PCA <-  mother.PCA[ mother.PCA$Focal != "BB LUC", ]
mother.PCA<-  mother.PCA[ mother.PCA$Focal != "BB KIN", ]
mother.PCA <-  mother.PCA[ mother.PCA$Focal != "KIG", ]
mother.PCA <-  mother.PCA[ mother.PCA$Focal != "AYA", ]
mother.PCA <-  mother.PCA[ mother.PCA$Focal != "JAM", ]
mother.PCA <-  mother.PCA[ mother.PCA$Focal != "PUS", ]
mother.PCA <-  mother.PCA[ mother.PCA$Focal != "LIR", ]
mother.PCA <-  mother.PCA[ mother.PCA$Focal != "LIM", ]
mother.PCA <-  mother.PCA[ mother.PCA$Focal != "XSO", ]
mother.PCA <-  mother.PCA[ mother.PCA$Focal != "PLA", ]
mother.PCA <-  mother.PCA[ mother.PCA$Focal != "NOI", ]
mother.PCA <-  mother.PCA[ mother.PCA$Focal != "NIG", ]
mother.PCA <-  mother.PCA[ mother.PCA$Focal != "BB OPA", ]
mother.PCA <-  mother.PCA[ mother.PCA$Focal != "BB POL", ]
mother.PCA <-  mother.PCA[ mother.PCA$Focal != "BB XEL", ]
mother.PCA <-  mother.PCA[ mother.PCA$Focal != "BB UJA", ]
mother.PCA <-  mother.PCA[ mother.PCA$Focal != "EDO", ]
mother.PCA <-  mother.PCA[ mother.PCA$Focal != "HAV", ]


mother.PCA <- mother.PCA[,2:7]

#correlation test because data not normal and frequency

# Filter out the character columns
numeric_data <- mother.PCA %>% select_if(is.numeric)

#install.packages("corrr")
library(corrr)
# Calculate Kendall correlation matrix
kendall_cor_matrix <- numeric_data %>% correlate(method = "kendall")

kendall_cor_matrix


#parallel analysis
library(psych)
result <- fa.parallel(numeric_data, fm = 'ml', fa = 'fa', n.iter = 1000, show.legend = FALSE)


# check for missing values

colSums(is.na(mother.PCA))
rownames(mother.PCA)<-mother.PCA$Focal

# normalizing the data

numerical.mother.PCA <- mother.PCA[,2:6]


#prcomp method
res.pca.mother <- prcomp(numerical.mother.PCA, scale = F)


# Eigenvalues
eig.val <- get_eigenvalue(res.pca.mother)
eig.val

summary(res.pca.mother)

# visualization

mother.PCA$attachment.type[mother.PCA$Focal == "AKU"] <- "Secure"
mother.PCA$attachment.type[mother.PCA$Focal == "ILA"] <- "Secure-anxious"
mother.PCA$attachment.type[mother.PCA$Focal == "IZA"] <- "Secure-anxious"
mother.PCA$attachment.type[mother.PCA$Focal == "NAK"] <- "Secure-anxious"
mother.PCA$attachment.type[mother.PCA$Focal == "POE"] <- "Secure"
mother.PCA$attachment.type[mother.PCA$Focal == "PRI"] <- "Secure"
mother.PCA$attachment.type[mother.PCA$Focal == "PAD"] <- "Insecure-avoidant"
mother.PCA$attachment.type[mother.PCA$Focal == "OFI"] <- "Secure"
mother.PCA$attachment.type[mother.PCA$Focal == "FLE"] <- "Secure"
mother.PCA$attachment.type[mother.PCA$Focal == "HAI"] <- "Secure-anxious"
mother.PCA$attachment.type[mother.PCA$Focal == "XIA"] <- "Secure"
mother.PCA$attachment.type[mother.PCA$Focal == "NEO"] <- "Secure"
mother.PCA$attachment.type[mother.PCA$Focal == "XOR"] <- "Insecure-avoidant"
mother.PCA$attachment.type[mother.PCA$Focal == "ETO"] <- "Insecure-avoidant"
mother.PCA$attachment.type[mother.PCA$Focal == "COM"] <- "Secure"
mother.PCA$attachment.type[mother.PCA$Focal == "URA"] <- "Insecure-avoidant"
mother.PCA$attachment.type[mother.PCA$Focal == "YET"] <- "Insecure-avoidant"
mother.PCA$attachment.type[mother.PCA$Focal == "RAJ"] <- "Insecure-avoidant"


fviz_pca_biplot(
  res.pca.mother,
  col.var = "black",
  palette = c( "dodgerblue2","brown", "orange"),
  col.ind = mother.PCA$attachment.type,
  pointshape = 16,
  legend.title = "",
  legend = "top",
  pointsize = 8,
  labelsize = 6,
  geom = "point"
) +
  theme(
    legend.text = element_text(size = 18),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18)
  )




score<-res.pca.mother$x
write.csv2(score,"score.pca.csv")

# 2. Offspring attachment type predicted by maternal protection and rejection ####




# Figure 2: Variation of maternal rejection (refusal of solicitations) with maternal protection (mother providing protection to the offspring during threatening events) received by each offspring ####


refusal.intercept<-read.csv("rejection-intercepts.csv",sep=";",dec=",")
protection.intercept<-read.csv("protection.intercepts.csv",sep=";",dec=",")


behaviours.intercept<-merge(refusal.intercept,protection.intercept,by=("Focal"))


behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "BB LUC", ]
behaviours.intercept<-  behaviours.intercept[ behaviours.intercept$Focal != "BB KIN", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "KIG", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "AYA", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "JAM", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "PUS", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "LIR", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "LIM", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "XSO", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "PLA", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "NOI", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "NIG", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "BB OPA", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "BB POL", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "BB XEL", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "BB UJA", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "EDO", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "HAV", ]


#make nice plot with red lines and colors no names

behaviours.intercept$attachment.type[behaviours.intercept$Focal == "AKU"] <- "Secure"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "ILA"] <- "Secure-anxious"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "IZA"] <- "Secure-anxious"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "NAK"] <- "Secure-anxious"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "POE"] <- "Secure"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "PRI"] <- "Secure"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "PAD"] <- "Insecure-avoidant"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "OFI"] <- "Secure"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "FLE"] <- "Secure"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "HAI"] <- "Secure-anxious"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "XIA"] <- "Secure"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "NEO"] <- "Secure"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "XOR"] <- "Insecure-avoidant"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "ETO"] <- "Insecure-avoidant"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "COM"] <- "Secure"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "URA"] <- "Insecure-avoidant"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "YET"] <- "Insecure-avoidant"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "RAJ"] <- "Insecure-avoidant"


attachment_colors <- c("Secure" = "brown", "Secure-anxious" = "orange", "Insecure-avoidant" = "dodgerblue2")

# Ensure that attachment.type is a factor
behaviours.intercept$attachment.type <- factor(behaviours.intercept$attachment.type, 
                                               levels = c("Insecure-avoidant","Secure", "Secure-anxious" ))



ggplot(behaviours.intercept, aes(x = rejection.Intercept, y = protection.Intercept, color = attachment.type)) +
  geom_point(size = 7) + # Increased size of points
  scale_color_manual(values = attachment_colors) + # Custom colors
  geom_vline(xintercept = 0, color = "black", linetype = "dashed") + # Vertical black line at x = 0
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") + # Horizontal black line at y = 0
  labs(x = "Deviation of refusal of solicitations", y = "Deviation of protection", 
       title = "Maternal styles") + # Updated Labels
  theme_minimal() + # Clean theme
  theme(legend.title = element_blank(), # Remove legend title
        legend.position = "top",legend.text = element_text(size = 18),
        axis.text.x = element_text(size = 18), # Adjust x-axis text size
        axis.text.y = element_text(size = 18), # Adjust legend text size
        axis.title.x = element_text(size = 18), # Adjust x-axis tick label size
        axis.title.y = element_text(size = 18), title = element_text(size = 20, face = "bold"), # Position legend at the top
        legend.justification = "center") # Justify legend to the right



# Figure 3. Individual variation in maternal protection (a) and rejection (b) with offspring’s attachment type ####

refusal.intercept<-read.csv("rejection-intercepts.csv",sep=";",dec=",")
danger.intercept<-read.csv("protection.intercepts.csv",sep=";",dec=",")

behaviours.intercept<-merge(danger.intercept,refusal.intercept,by="Focal")


behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "BB LUC", ]
behaviours.intercept<-  behaviours.intercept[ behaviours.intercept$Focal != "BB KIN", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "KIG", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "AYA", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "JAM", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "PUS", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "LIR", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "LIM", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "XSO", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "PLA", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "NOI", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "NIG", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "BB OPA", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "BB POL", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "BB XEL", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "BB UJA", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "EDO", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "HAV", ]


#make nice plot with red lines and colors no names

behaviours.intercept$attachment.type[behaviours.intercept$Focal == "AKU"] <- "Secure"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "ILA"] <- "Secure-anxious"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "IZA"] <- "Secure-anxious"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "NAK"] <- "Secure-anxious"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "POE"] <- "Secure"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "PRI"] <- "Secure"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "PAD"] <- "Insecure-avoidant"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "OFI"] <- "Secure"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "FLE"] <- "Secure"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "HAI"] <- "Secure-anxious"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "XIA"] <- "Secure"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "NEO"] <- "Secure"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "XOR"] <- "Insecure-avoidant"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "ETO"] <- "Insecure-avoidant"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "COM"] <- "Secure"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "URA"] <- "Insecure-avoidant"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "YET"] <- "Insecure-avoidant"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "RAJ"] <- "Insecure-avoidant"


attachment_colors <- c("Secure" = "brown", "Secure-anxious" = "gold", "Insecure-avoidant" = "dodgerblue2")

# Ensure that attachment.type is a factor
behaviours.intercept$attachment.type <- factor(behaviours.intercept$attachment.type, 
                                               levels = c("Insecure-avoidant","Secure", "Secure-anxious"))

behaviours.intercept$Approach.grab<-behaviours.intercept$danger.Intercept
behaviours.intercept$Refusal<-behaviours.intercept$refusal.Intercept


#plot refusal and protectiveness


install.packages("patchwork")
library(patchwork)

# First Plot: Approach.grab (Protection)
p1 <- ggplot(behaviours.intercept, aes(x = attachment.type, y = protection.Intercept, fill = attachment.type)) +
  geom_boxplot(outlier.shape = NA, color = "black") +  # Boxplot with black outlines, fill mapped to attachment.type
  geom_jitter(width = 0.2, size = 2, color = "black") +  # Jitter points in black
  scale_fill_manual(values = c("Secure" = "brown", "Secure-anxious" = "orange", "Insecure-avoidant" = "dodgerblue2")) +  # Custom colors for box fill
  labs(x = "Attachment type", 
       y = "Protection Intercept (Individual-level)",  # Updated y-axis label
       title = "(a) Attachment type and maternal protection") +
  theme_minimal() +
  theme(axis.title = element_text(size = 13),   # Larger axis titles
        axis.text = element_text(size = 13),    # Larger axis text
        plot.title = element_text(size = 13, face = "bold", hjust = 0.5),  # Larger, bold title, centered
        legend.position = "none")  # Remove legend

# Second Plot: Refusal (Rejection)
p2 <- ggplot(behaviours.intercept, aes(x = attachment.type, y = rejection.Intercept, fill = attachment.type)) +
  geom_boxplot(outlier.shape = NA, color = "black") +  # Boxplot with black outlines, fill mapped to attachment.type
  geom_jitter(width = 0.2, size = 2, color = "black") +  # Jitter points in black
  scale_fill_manual(values = c("Secure" = "brown", "Secure-anxious" = "orange", "Insecure-avoidant" = "dodgerblue2")) +  # Custom colors for box fill
  labs(x = "Attachment type", 
       y = "Rejection Intercept (Individual-level)",  # Updated y-axis label
       title = "(b) Attachment type and maternal rejection") +
  theme_minimal() +
  theme(axis.title = element_text(size = 13),   # Larger axis titles
        axis.text = element_text(size = 13),    # Larger axis text
        plot.title = element_text(size = 13, face = "bold", hjust = 0.5),  # Larger, bold title, centered
        legend.position = "none")  # Remove legend



# Combine p1 and p2 side by side
combined_plot <- plot_grid(p1, p2, ncol = 2)

# Print the combined plot
print(combined_plot)




#violin plot

p1 <- ggplot(behaviours.intercept, aes(x = attachment.type, y = protection.Intercept, fill = attachment.type)) +
  geom_violin(color = "black", trim = FALSE, alpha = 0.4) +  # Violin plot with some transparency
  geom_boxplot(width = 0.2, color = "black", outlier.shape = NA, alpha = 1) +  # Boxplot inside the violin plot
  geom_jitter(width = 0.2, size = 2, color = "black") +  # Jittered points
  scale_fill_manual(values = c("Secure" = "brown", "Secure-anxious" = "orange", "Insecure-avoidant" = "dodgerblue2")) +  # Custom fill colors
  labs(x = "Attachment type", 
       y = "Protection Intercept (Individual-level)", 
       title = "(a) Attachment type and maternal protection") +
  theme_minimal() +
  theme(axis.title = element_text(size = 13), 
        axis.text = element_text(size = 13), 
        plot.title = element_text(size = 13, face = "bold", hjust = 0.5), 
        legend.position = "none")  # Remove legend



p2 <- ggplot(behaviours.intercept, aes(x = attachment.type, y = rejection.Intercept, fill = attachment.type)) +
  geom_violin(color = "black", trim = FALSE, alpha = 0.4) +  # Violin plot with some transparency
  geom_boxplot(width = 0.2, color = "black", outlier.shape = NA, alpha = 1) +  # Boxplot inside the violin plot
  geom_jitter(width = 0.2, size = 2, color = "black") +  # Jittered points
  scale_fill_manual(values = c("Secure" = "brown", "Secure-anxious" = "orange", "Insecure-avoidant" = "dodgerblue2")) +  # Custom fill colors
  labs(x = "Attachment type", 
       y = "Rejection Intercept (Individual-level)", 
       title = "(b) Attachment type and maternal rejection") +
  theme_minimal() +
  theme(axis.title = element_text(size = 13), 
        axis.text = element_text(size = 13), 
        plot.title = element_text(size = 13, face = "bold", hjust = 0.5), 
        legend.position = "none")  # Remove legend


# Combine p1 and p2 side by side
combined_plot <- plot_grid(p1, p2, ncol = 2)

# Print the combined plot
print(combined_plot)



# Model 6 and Table 4: Attachment types predicted by maternal protection and rejection  ####

behaviours.intercept<-read.csv("score.pca.csv",sep=";",dec=",")
colnames(behaviours.intercept)<-c("Focal","PC1","PC2","PC3","PC4","PC5")



behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "BB LUC", ]
behaviours.intercept<-  behaviours.intercept[ behaviours.intercept$Focal != "BB KIN", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "KIG", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "AYA", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "JAM", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "PUS", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "LIR", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "LIM", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "XSO", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "PLA", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "NOI", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "NIG", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "BB OPA", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "BB POL", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "BB XEL", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "BB UJA", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "EDO", ]
behaviours.intercept <-  behaviours.intercept[ behaviours.intercept$Focal != "HAV", ]


#make nice plot with red lines and colors no names

behaviours.intercept$attachment.type[behaviours.intercept$Focal == "AKU"] <- "Secure"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "ILA"] <- "Secure-anxious"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "IZA"] <- "Secure-anxious"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "NAK"] <- "Secure-anxious"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "POE"] <- "Secure"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "PRI"] <- "Secure"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "PAD"] <- "Insecure-avoidant"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "OFI"] <- "Secure"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "FLE"] <- "Secure"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "HAI"] <- "Secure-anxious"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "XIA"] <- "Secure"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "NEO"] <- "Secure"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "XOR"] <- "Insecure-avoidant"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "ETO"] <- "Insecure-avoidant"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "COM"] <- "Secure"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "URA"] <- "Insecure-avoidant"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "YET"] <- "Insecure-avoidant"
behaviours.intercept$attachment.type[behaviours.intercept$Focal == "RAJ"] <- "Insecure-avoidant"


attachment_colors <- c("Secure" = "brown", "Secure-anxious" = "gold", "Insecure-avoidant" = "dodgerblue2")

# Ensure that attachment.type is a factor
behaviours.intercept$attachment.type <- factor(behaviours.intercept$attachment.type, 
                                               levels = c("Insecure-avoidant","Secure", "Secure-anxious"))

#no interaction

behaviours.intercept$attachment.type <- factor(behaviours.intercept$attachment.type, 
                                               levels = c( "Insecure-avoidant","Secure-anxious","Secure"))

mpriorj.protectiveness = get_prior(attachment.type ~ PC1+ PC2+PC3
                                   ,
                                   data = behaviours.intercept, family = categorical(link="logit"))


# you can and should add slightly stronger priors to your predictors.
# To do this, look at the mprior object, and find which rows contain your predictor fixed effects 
# These will appear three times because you have three levels of response

mpriorj.protectiveness$prior[4:6] <- "normal(0,1)"# so here change the indexing to the rows containing your fixed effects

mpriorj.protectiveness$prior[9:11] <- "normal(0,1)"# again here.
mpriorj.protectiveness
# the next function "make_stancode" makes sure that your prior and model formula will run in stan (the model estimator).

make_stancode(attachment.type ~ PC1+ PC2+PC3
              ,
              data = behaviours.intercept, family = categorical, prior = mpriorj.protectiveness)

# if you do not get an error message, you can run your model
# if you get an error message, fix the error! it is usually if you have a variable in the above formula
# that isn't in your prior formula



b.mother.protection2.no.inte = brm(attachment.type ~ PC1+ PC2+PC3
                                   ,
                                   data = behaviours.intercept,
                                   family = categorical, 
                                   prior = mpriorj.protectiveness, 
                                   chains = 12, 
                                   cores = 12, 
                                   iter = 2000, 
                                   warmup = 1000,
                                   control = list(adapt_delta = 0.75), 
                                   sample_prior = "yes")


summary(b.mother.protection2.no.inte)

summary(b.mother.protection2.no.inte, prob=0.89)

conditional_effects(b.mother.protection2.no.inte, categorical=T)

# Extract fixed effects (including the intercept)
fixed_effects <- fixef(b.mother.protection2.no.inte)

# Display the fixed effects
print(fixed_effects)

# Get the intercept for the reference category
intercept_reference <- fixed_effects[1, "Estimate"]  # First row is the intercept
print(intercept_reference)


hypothesis(b.mother.protection2.no.inte, "muSecureanxious_PC1   =muSecure_PC1  ", class = "b", group = "",
           scope = c("standard", "ranef", "coef"), alpha = 0.11, seed = NULL)


hypothesis(b.mother.protection2.no.inte, "muSecureanxious_PC2  =muSecure_PC2  ", class = "b", group = "",
           scope = c("standard", "ranef", "coef"), alpha = 0.11, seed = NULL)

hypothesis(b.mother.protection2.no.inte, "muSecureanxious_PC3   =muSecure_PC3  ", class = "b", group = "",
           scope = c("standard", "ranef", "coef"), alpha = 0.05, seed = NULL)












behaviours.intercept$attachment.type <- factor(behaviours.intercept$attachment.type, 
                                               levels = c( "Insecure-avoidant","Secure-anxious","Secure"))

mpriorj.protectiveness = get_prior(attachment.type ~ PC1
                                   ,
                                   data = behaviours.intercept, family = categorical(link="logit"))


# you can and should add slightly stronger priors to your predictors.
# To do this, look at the mprior object, and find which rows contain your predictor fixed effects 
# These will appear three times because you have three levels of response

mpriorj.protectiveness$prior[4:5] <- "normal(0,1)"# so here change the indexing to the rows containing your fixed effects
mpriorj.protectiveness
# the next function "make_stancode" makes sure that your prior and model formula will run in stan (the model estimator).

make_stancode(attachment.type ~ PC1
              ,
              data = behaviours.intercept, family = categorical, prior = mpriorj.protectiveness)

# if you do not get an error message, you can run your model
# if you get an error message, fix the error! it is usually if you have a variable in the above formula
# that isn't in your prior formula



b.PC1 = brm(attachment.type ~ PC1
            ,
            data = behaviours.intercept,
            family = categorical, 
            prior = mpriorj.protectiveness, 
            chains = 12, 
            cores = 12, 
            iter = 2000, 
            warmup = 1000,
            control = list(adapt_delta = 0.75), 
            sample_prior = "yes")


summary(b.PC1)

summary(b.PC1, prob=0.89)

hypothesis(b.PC1, "muSecureanxious_PC1   =muSecure_PC1  ", class = "b", group = "",
           scope = c("standard", "ranef", "coef"), alpha = 0.11, seed = NULL)





behaviours.intercept$attachment.type <- factor(behaviours.intercept$attachment.type, 
                                               levels = c( "Insecure-avoidant","Secure-anxious","Secure"))

mpriorj.protectiveness = get_prior(attachment.type ~ PC2
                                   ,
                                   data = behaviours.intercept, family = categorical(link="logit"))


# you can and should add slightly stronger priors to your predictors.
# To do this, look at the mprior object, and find which rows contain your predictor fixed effects 
# These will appear three times because you have three levels of response

mpriorj.protectiveness$prior[4:5] <- "normal(0,1)"# so here change the indexing to the rows containing your fixed effects
mpriorj.protectiveness
# the next function "make_stancode" makes sure that your prior and model formula will run in stan (the model estimator).

make_stancode(attachment.type ~ PC2
              ,
              data = behaviours.intercept, family = categorical, prior = mpriorj.protectiveness)

# if you do not get an error message, you can run your model
# if you get an error message, fix the error! it is usually if you have a variable in the above formula
# that isn't in your prior formula



b.PC2 = brm(attachment.type ~ PC2
            ,
            data = behaviours.intercept,
            family = categorical, 
            prior = mpriorj.protectiveness, 
            chains = 12, 
            cores = 12, 
            iter = 2000, 
            warmup = 1000,
            control = list(adapt_delta = 0.75), 
            sample_prior = "yes")


summary(b.PC2)

summary(b.PC2, prob=0.89)

hypothesis(b.PC2, "muSecureanxious_PC2   =muSecure_PC2  ", class = "b", group = "",
           scope = c("standard", "ranef", "coef"), alpha = 0.11, seed = NULL)





behaviours.intercept$attachment.type <- factor(behaviours.intercept$attachment.type, 
                                               levels = c( "Insecure-avoidant","Secure-anxious","Secure"))

mpriorj.protectiveness = get_prior(attachment.type ~ PC3
                                   ,
                                   data = behaviours.intercept, family = categorical(link="logit"))


# you can and should add slightly stronger priors to your predictors.
# To do this, look at the mprior object, and find which rows contain your predictor fixed effects 
# These will appear three times because you have three levels of response

mpriorj.protectiveness$prior[4:5] <- "normal(0,1)"# so here change the indexing to the rows containing your fixed effects
mpriorj.protectiveness
# the next function "make_stancode" makes sure that your prior and model formula will run in stan (the model estimator).

make_stancode(attachment.type ~ PC3
              ,
              data = behaviours.intercept, family = categorical, prior = mpriorj.protectiveness)

# if you do not get an error message, you can run your model
# if you get an error message, fix the error! it is usually if you have a variable in the above formula
# that isn't in your prior formula



b.PC3 = brm(attachment.type ~ PC3
            ,
            data = behaviours.intercept,
            family = categorical, 
            prior = mpriorj.protectiveness, 
            chains = 12, 
            cores = 12, 
            iter = 2000, 
            warmup = 1000,
            control = list(adapt_delta = 0.75), 
            sample_prior = "yes")


summary(b.PC3)

summary(b.PC3, prob=0.89)


hypothesis(b.PC3, "muSecureanxious_PC3   =muSecure_PC3  ", class = "b", group = "",
           scope = c("standard", "ranef", "coef"), alpha = 0.11, seed = NULL)




#Figure 3: Individual variation in maternal protection (a) and rejection (b) with offspring’s attachment type ####




# Figure S2: Posterior predictive checks for all models  ####
library(brms)
library(bayesplot)

refusal <- pp_check(b.refusal.month.pregnant2, nsamples = 1000) + ggtitle("Rejection (Model 2)")
protection <- pp_check(b.danger.month.pregnant.2, nsamples = 1000) + ggtitle("Protection (Model 1)")
Groom <- pp_check(b.groomed.model.no1800, nsamples = 1000) + ggtitle("Grooming (Model 4)")+
  xlim(0, 100)
proximity<- pp_check(b.contact.no10, nsamples = 1000) + ggtitle("Proximity (Model 5)")
play <- pp_check(b.model.playing.no1800, nsamples = 1000) + ggtitle("Playing (Model 3)")+
  xlim(0, 100)
rejection.protection.attachment <- pp_check(b.mother.protection2.no.inte, nsamples = 1000) + ggtitle("Protection and rejection predicting attachment types (Model 6)")
rejection.protection.attachment<-pp_check(b.mother.protection2.no.inte, type = "bars")+ ggtitle("Protection and rejection \npredicting attachment types (Model 6)")

# Create a list of PPC plots
plot_list <- list(protection,refusal, play, Groom, proximity,rejection.protection.attachment )

# Ensure plot_list is not empty before calling grid.arrange
if (length(plot_list) > 0) {
  # Arrange the plots in a grid
  grid.arrange(grobs = plot_list, nrow = 3)  # Adjust the number of rows as needed
} else {
  cat("No PPC plots generated.\n")
}




grid.arrange(grobs = plot_list, nrow = 3)














