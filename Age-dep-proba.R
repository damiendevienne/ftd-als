## D de Vienne & DM de Vienne
## May 6th, 20025
require(ggplot2)

##################
## FORMAT DATA ###
##################

## Get and format data from Murphy et al. 2017 : "Age related penetrance of the C9orf72 repeat expansion. Scientific reports."
data <- read.csv("C9_Data.csv", sep=";", header=TRUE)
data_clean <- data[data$Event==1,] #remove the 23 patients without symptoms
data_clean <- data_clean[!is.na(data_clean$AgeofOnset),] #remove patient for which age of onset is missing


## Count cases per age
min_age <- 0
max_age <- max(data$AgeofOnset, na.rm = TRUE) #na.rm = TRUE to avoid NA in the max
#CODE DAMIEN TROP COMPLIQUÉ : age_onset <- sapply(min_age:max_age, function(x, onsets) sum(onsets == x), onsets = data_clean$AgeofOnset)
# CODE ALTERNATIF :
age_onset <- array()
for (i in min_age:max_age) {
  age_onset[i+1] <- sum(data_clean$AgeofOnset == i) #i+1 because age 0 is at index 1
}

df_onset <- data.frame(Age = min_age:max_age, num = age_onset)


##########################
## PROBA. COMPUTATIONS ###
##########################

## Compute Penetrance by age (proba of declaring symptoms at each age if carrying the variant), taking into account a possible proportion (immunisedprop) of immunised people in the population: 

df_onset$cumsum <- cumsum(df_onset$num)
immunisedprop <- 0 #proportion of immunised people in the population (between 0 and 1)
immunisednb<-round((sum(df_onset$num)*immunisedprop)/(1-immunisedprop)) # number of immunised people in the population
df_onset$Penetrance <- df_onset$cumsum / (max(df_onset$cumsum)+immunisednb) #Penetrance

## Compute proba of being healthy by age if carrying the variant (1 - penetrance)
df_onset$pHealthy <- 1 - df_onset$Penetrance

# Compute probability of carrying the variant if being healthy at a given age (see article for justification of the formula). 
df_onset$PvariantWhenHealthy <- (1 - df_onset$Penetrance) / (2 - df_onset$Penetrance)

###########################################
## PROBA. COMPUTATIONS FOR GRANDCHILDS ####
###########################################

# Without any other information, a grandchild has a proba of 0.25 to carry the variant. So the proba of carrying the variant if being healthy at a given age is given by (see article):

df_onset$grandchild_default <- ((df_onset$pHealthy)*0.25)/(0.75+0.25*df_onset$pHealthy)

# But the probability can be corrected for a grandchild whose parent is healthy and whose carrier status is unknown. This proba depends on the age difference between the grandchild and its parent.

# Compute for grandchild of a given age the proba that its direct parent carries the variant given that he/she is healthy: 
pv_dad_for_child <-function (pvdad, age_interval) {
  res <- 0.5*c(pvdad[(age_interval+1):length(pvdad)], rep(min(pvdad), age_interval)) #multiplied by 0.5 because child has 1/2 chance of inheriting the variant
  res
}
# for diff_age years of difference between grandchild and his parent: 
diff_age <- 31 
dad_proba_healthy <- pv_dad_for_child(df_onset$PvariantWhenHealthy, diff_age)

# This value of dad_proba_healthy replaces the 0.25 probability of carrying the variant used above. 
df_onset$grandchilddiff <- ((df_onset$pHealthy)*dad_proba_healthy)/(1-dad_proba_healthy+dad_proba_healthy*df_onset$pHealthy)

# if using dominique's equation, we have: 
df_onset$grandchilddiff2 <- (1-df_onset$Penetrance)/((1/dad_proba_healthy)-df_onset$Penetrance) #I use 1/dad_proba_healthy instead of 2/dad_proba_healthy because I already diviedd by 2.


#pdf("results.pdf")
X11()
ggplot(df_onset, aes(x=Age)) + geom_step(aes(y=Penetrance), col="black", show.legend = TRUE) + geom_step(aes(y=pHealthy), col="grey") + geom_step(aes(y=PvariantWhenHealthy), col="red", lwd=1) + ylab("Probability") + geom_step(aes(y=grandchild_default), col="darkgreen", lwd=1) + geom_step(aes(y=grandchilddiff), col="darkblue", lwd=1)
#dev.off()

