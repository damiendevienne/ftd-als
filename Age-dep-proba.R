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
max_age <- 100
## Get the number of cases per age 
age_onset <- sapply(min_age:max_age, function(x, onsets) sum(onsets == x), onsets = data_clean$AgeofOnset)

## Data frame for all computed data
df_onset <- data.frame(Age = min_age:max_age, num = age_onset)

###############################
## PROBA. COMPUTATIONS (F1) ###
###############################

# 1. COMPLETE PENETRANCE OF THE VARIANT

# Penetrance
df_onset$cumsum <- cumsum(df_onset$num)
df_onset$Penetrance <- df_onset$cumsum / (max(df_onset$cumsum)) #Penetrance
# Probability of carrying the variant when healthy ()  
df_onset$PvariantWhenHealthy <- (1 -df_onset$Penetrance) / (2 - df_onset$Penetrance)
# Probability of being sick in the future
df_onset$PsickinFuture <- df_onset$PvariantWhenHealthy

# 2. INCOMPLETE PENETRANCE OF THE VARIANT
# Proportion of penetrant individuals in the population is controled by x. 
# For varying x and save all results we create a new data frame: df_onset_expanded

x_vals <- seq(0, 1, by = 0.2) # possible values of x

# Initialize empty data frame
df_onset_expanded <- data.frame()
# Loop to compute AdjustedPenetrance for each x
for (x in x_vals) {
  temp <- data.frame(Age = df_onset$Age, x = x, 
    PvariantWhenHealthy = (1-(1-x)*df_onset$Penetrance)/(2-(1-x)*df_onset$Penetrance), 
    PsickinFuture = ((1-x)*(1-df_onset$Penetrance))/(2-(1-x)*df_onset$Penetrance)) 
  df_onset_expanded <- rbind(df_onset_expanded, temp)
}

# SAVE THE DATA FOR USE IN SHINY APP
saveRDS(df_onset, file="shiny-app/df_onset.rds")



#####################
## PLOT THE DATA ####
#####################
require(ggplot2)
library(patchwork)

####  FIGURE 1  #####
## 1A: number of cases per age
f1 <- ggplot(df_onset, aes(x = Age, y = num)) +
  geom_col(fill = "steelblue") +
  scale_x_continuous(breaks = seq(0, max(df_onset$Age), by = 10)) +
  ylab("Number of cases") +
  theme(plot.tag = element_text(size = 16, face = "bold")) + 
  labs(tag = "A")

## 1B: Penetrance
f2 <- ggplot(df_onset, aes(x = Age, y = Penetrance)) +
  geom_step(color = "steelblue", lwd = 0.6) +
  scale_x_continuous(breaks = seq(0, max(df_onset$Age), by = 10)) +
  ylab("Penetrance") +
  theme(plot.tag = element_text(size = 16, face = "bold")) +
  labs(tag = "B")

## Combine and save to pdf
pdf("fig1.pdf", width=12)
f1 + f2
dev.off()


####  FIGURE 2  #####
# Set colors
colors <- c(
  "0"   = "#4682b4",    # deep blue
  "0.2" = "#00a2c8",      # medium blue-purple
  "0.4" = "#00bfc5",      # muted magenta
  "0.6" = "#39d9ab",      # warm coral red
  "0.8" = "#9ded89",      # bright orange-red
  "1"   = "#f9f871"           # vivid red
)
## 2A: proba to carry mutation for various x 
p1 <- ggplot(df_onset_expanded, aes(x = Age, y = PvariantWhenHealthy, color = factor(x))) +
  geom_step(aes(size = (x == "0"))) +
  labs(y = "Probability of carrying the c9orf72 RE variant", color = "x") +
  scale_size_manual(values = c("TRUE" = 1, "FALSE" = 0.6), guide = "none") +
  scale_color_manual(values = colors) + ylim(c(0,1)) + 
  scale_x_continuous(breaks = seq(0, max(df_onset$Age), by = 10)) +
  labs(color = expression(italic(x))) +
  theme(plot.tag = element_text(size = 16, face = "bold"), legend.position = "none") + 
  labs(tag = "A")

## 2B: proba to die of the disease
p2 <- ggplot(df_onset_expanded, aes(x = Age, y = PsickinFuture, color = factor(x))) +
  geom_step(aes(size = (x == "0"))) +
  labs(y = "Probability of developing the disease in the future", color = "x") +
  scale_size_manual(values = c("TRUE" = 1, "FALSE" = 0.6), guide = "none") +
  scale_color_manual(values = colors) + ylim(c(0,1)) + 
  scale_x_continuous(breaks = seq(0, max(df_onset$Age), by = 10)) +
  labs(color = expression(italic(x))) +
  guides(color = guide_legend(override.aes = list(size = c(1, rep(0.6, 5))))) +
  theme(plot.tag = element_text(size = 16, face = "bold")) + 
  labs(tag = "B")

# Combine and save to pdf
pdf("fig2.pdf", width=12)
(p1 | p2) + plot_layout(guides = "collect") & theme(legend.position = "right")
dev.off()











###############################
## PROBA. COMPUTATIONS (F2) ###               TODO
###############################




#### TRY FOR GRANDCHILDRENS with 30 years difference with dad
delta <- 30 #age difference between grandchild and parent
penetrance_delta <- c(df_onset$Penetrance[(delta+1):length(df_onset$Penetrance)],rep(1, delta))

penetrances <- cbind(df_onset$Penetrance, penetrance_delta)
K <- apply(penetrances, 1, function(x) pVWH_grandchild(0.2, x[1], x[2]))
K2 <- apply(penetrances, 1, function(x) pVWH_grandchild2(0.2, x[1], x[2]))


# equation damien
pVWH_grandchild <- function(x, pi_t_prime, pi_t) {
  numerator <- pi_t_prime * x + pi_t * x - pi_t_prime - pi_t + 
               pi_t_prime * pi_t + pi_t_prime * x^2 * pi_t - 
               2 * pi_t_prime * pi_t * x + 1
  
  denominator <- pi_t_prime * x + 2 * pi_t * x - pi_t_prime - 2 * pi_t + 
                 pi_t_prime * pi_t + pi_t_prime * x^2 * pi_t - 
                 2 * pi_t_prime * pi_t * x + 4
  
  result <- numerator / denominator
  return(result)
}


# equation papa 
pVWH_grandchild2 <- function(x, pi_t_prime, pi_t) {
  numerator <- (1 - (1 - x) * pi_t) * (1 - (1 - x) * pi_t_prime)
  denominator <- 2 * (2 - (1 - x) * pi_t) - (1 - (1 - x) * pi_t) * (1 - x) * pi_t_prime
  result <- numerator / denominator
  return(result)
}

#equation damien
pSICK_grandchild <- function(x, pi_t_prime, pi_t) {
  k <- (1-x) * (1-pi_t)
  u <- (1-x) * (1-pi_t_prime)
  numerator <- x * u + k * u
  denominator <- 2+x+x^2+k+x*k+x*u+u*k
  result <- numerator / denominator  
  return(result)
}

#equation papa
pSICK_grandchild2 <- function(x, pi_t_prime, pi_t) {
  print(x)
  print(pi_t_prime)
  print(pi_t)
  numerator <- (1 - x) * ( 1 - pi_t_prime)*(1 - pi_t)
  denominator <- 2 * (2 - pi_t) - pi_t_prime * (1 - x) * (1 - pi_t)
  result <- numerator / denominator
  return(result)
}

#equation damien2
pSICK_grandchild3 <- function(x, pi_t_prime, pi_t) {
  k <- (1-x) * (1-pi_t)
  u <- (1-x) * (1-pi_t_prime)
  numerator <- u
  denominator <- u + x + 1 + ((4-2*pi_t_prime+2*pi_t*x)/(1-pi_t+pi_t*x))
  result <- numerator / denominator  
  return(result)
}

#equation papa2
pSICK_grandchild4 <- function(x, pi_t_prime, pi_t) {
  pprime <- (1 - x) * pi_t
  numerator <- (1 - x) * (1 - pi_t_prime) * (1-pprime)
  denominator <- 2 * (2 - pprime) - pi_t_prime * (1 - x) * (1 - pprime) 
  result <- numerator / denominator  
  return(result)
}



X<-0.4
S <- apply(penetrances, 1, function(x) pSICK_grandchild(X, x[1], x[2]))
S2 <- apply(penetrances, 1, function(x) pSICK_grandchild2(X, x[1], x[2]))
S3 <- apply(penetrances, 1, function(x) pSICK_grandchild3(X, x[1], x[2]))
S4 <- apply(penetrances, 1, function(x) pSICK_grandchild4(X, x[1], x[2]))
plot(S)
lines(S2, col="red")
lines(S3, col="green")
lines(S4, col="purple")


lines(S, col="green")
#plot(S2)
lines(S2)

#pdf("results.pdf")
X11()
ggplot(df_onset, aes(x=Age)) + geom_step(aes(y=Penetrance), col="black", show.legend = TRUE) + geom_step(aes(y=pHealthy), col="grey") + geom_step(aes(y=PvariantWhenHealthy), col="red", lwd=1) + ylab("Probability") + geom_step(aes(y=grandchild_default), col="darkgreen", lwd=1) + geom_step(aes(y=grandchilddiff), col="darkblue", lwd=1) + geom_step(aes(y=pAffected), col="yellow", lwd=2) 
#dev.off()

browser = "google-chrome")


options(
#FOR a given age and for all possible x values (proportion of immunised people in the population), we can compute the probability of being affected at a given age.
x<-seq(0, 1, 0.01)
age <- 72
OK<-sapply(x, function(x, penet, age) {
  pen <- penet[age+1]*(1-x) 
  probcarrybuthealthy <- (1-pen)/(2-pen)
  probsick <- (1-x)*probcarrybuthealthy
  probhealthy <- 1 - probcarrybuthealthy + x * probcarrybuthealthy
  return(c(x, pen, probcarrybuthealthy, probsick, probhealthy))
}, penet = df_onset$Penetrance, age=age)
OK<-t(OK)
X11()
plot(OK[,1], OK[,3], type="l", col="black", lwd=2, ylim=c(0,1))
lines(OK[,1], OK[,4], type="l", col="black", lwd=2)
# immunisedprop <- 0 #proportion of immunised people in the population (between 0 and 1)
# immunisednb<-round((sum(df_onset$num)*immunisedprop)/(1-immunisedprop)) # number of immunised people in the population






