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

data_clean2 <- data_clean[!is.na(data_clean$Sex),] # remove patients for which sex is missing
data_clean_F <- data_clean2[data_clean2$Sex=="F",]
data_clean_M <- data_clean2[data_clean2$Sex=="M",]


## Prepara data frame 
prepare_df_onset <- function(cleandata) {
  ## Count cases per age
  min_age <- 0
  max_age <- 100
  ## Get the number of cases per age 
  age_onset <- sapply(min_age:max_age, function(x, onsets) sum(onsets == x), onsets = cleandata$AgeofOnset)

  ## Data frame for all computed data
  df_onset <- data.frame(Age = min_age:max_age, num = age_onset)

  ## Get Penetrance from the number of cases per age
  df_onset$cumsum <- cumsum(df_onset$num)
  df_onset$Penetrance <- df_onset$cumsum / (max(df_onset$cumsum)) #Penetrance
  df_onset
}

df_onset <- prepare_df_onset(data_clean)
df_onset_F <- prepare_df_onset(data_clean_F)
df_onset_M <- prepare_df_onset(data_clean_M)

# Save data frame for use in web application
#saveRDS(df_onset, file="shiny-app/df_onset.rds")
write.csv(df_onset, "ftd-als/df_onset.csv", row.names = FALSE)

###############################
## PROBA. COMPUTATIONS      ###
###############################

### Probability of carrying the variant when healthy (pVWH)
# Child
pVWH_child <- function(x, pi_t) {
  pprime_t <- (1 - x) * pi_t
  result <- (1 - pprime_t) / (2 - pprime_t)
  result
}

# Grandchild
pVWH_grandchild <- function(x, pi_t2, pi_t1) {
  pprime_t1 <- (1 - x) * pi_t1
  pprime_t2 <- (1 - x) * pi_t2
  numerator <- (1 - pprime_t2)
  denominator <- 2 * ((2 - pprime_t1) / (1 - pprime_t1)) - pprime_t2
  result <- numerator / denominator
  result
}

### Probability of being affected at age t+n when unaffected at age t (pSICKn)
# Child
pSICKn_child <- function(x, pi_t, pi_t_n) { #pi_t_n is the penetrance at age t+n
  pprime_t <- (1 - x) * pi_t
  pprime_t_n <- (1 - x) * pi_t_n
  result <- (pprime_t_n -pprime_t)/(2-pprime_t)
  result
}

# Grandchild
pSICKn_grandchild <- function(x, pi_t2,pi_t2_n, pi_t1) { #pi_t2_n is the penetrance at age t2+n
  pprime_t1 <- (1 - x) * pi_t1
  pprime_t2 <- (1 - x) * pi_t2
  pprime_t2_n <- (1 - x) * pi_t2_n
  numerator <- (pprime_t2_n - pprime_t2)
  denominator <- 2 * ((2 - pprime_t1) / (1 - pprime_t1)) - pprime_t2
  print(denominator)
  result <- numerator / denominator
  result
}

#TEMP: we compute pMX, the proba, in the presence of the mutation, of being affected at age t+n when unaffected at age t for grandchild
pMX <- function(x, pi_t2, pi_t2_n) {
  pprime_t2 <- (1 - x) * pi_t2
  print(pprime_t2)
  pprime_t2_n <- (1 - x) * pi_t2_n
  print(pprime_t2_n)
  numerator <- (pprime_t2_n - pprime_t2)
  print(numerator)
  denominator <- 1 - pprime_t2
  print(denominator)
  result <- numerator / denominator
  result
}


################################
## COMPUTE PROBABILITIES  ######
## FROM MURPHY ET AL. DATA    ##
################################

## 1. CHILDREN - PROBABILITY OF CARRYING THE VARIANT WHEN HEALTHY
getPVWH_child <- function(df, xvals=seq(0, 1, by = 0.2)) {
  # Initialize empty data frame
  df_onset_expanded <- data.frame()
  for (x in xvals) {
      temp <- data.frame(Age = df$Age, x = x, 
      PvariantWhenHealthy = pVWH_child(x, df_onset$Penetrance))
      df_onset_expanded <- rbind(df_onset_expanded, temp)
  }
  df_onset_expanded
}

## 2. GRANDCHILDREN - PROBABILITY OF CARRYING THE VARIANT WHEN HEALTHY
getPVWH_grandchild <- function(df, xvals=seq(0, 1, by = 0.2), delta) { 
  # delta is the age difference between grandchild and its parent (child of the carrier)
  penetrance_delta <- c(df$Penetrance[(delta+1):length(df$Penetrance)],rep(1, delta))
  penetrances <- cbind(df$Penetrance, penetrance_delta) # Create a matrix with two columns: Penetrance at age t and Penetrance at age t+delta (HERE, AGE OF THE PARENT) 
  # Initialize empty data frame
  df_onset_expanded <- data.frame()
  for (x in xvals) {
    temp <- data.frame(Age=df$Age, x = x, delta=delta, 
    PvariantWhenHealthy = apply(penetrances, 1, function(k,x) pVWH_grandchild(x, k[1], k[2]),x=x))
    df_onset_expanded <- rbind(df_onset_expanded, temp)
  }
df_onset_expanded
}

## 3. CHILDREN - PROBABILITY OF GETTING SICK IN THE NEXT N YEARS
getpSICK_child <- function(df, xvals=seq(0, 1, by = 0.2), n) {
  # n is the number of years in the future to compute the probability of getting sick
  penetrance_n <- c(df$Penetrance[(n+1):length(df$Penetrance)],rep(1, n))
  penetrances <- cbind(df$Penetrance, penetrance_n) # Create a matrix with two columns: Penetrance at age t and Penetrance at age t+n
  # Initialize empty data frame
  df_onset_expanded <- data.frame()
  for (x in xvals) {
      temp <- data.frame(Age = df$Age, x = x, n=n,
      pSICK = apply(penetrances, 1, function(k,x) pSICKn_child(x, k[1], k[2]), x=x))
      df_onset_expanded <- rbind(df_onset_expanded, temp)
  }
  df_onset_expanded
}

## 4. GRANDCHILDREN - PROBABILITY OF GETTING SICK IN THE NEXT N YEARS
getpSICK_grandchild <- function(df, xvals=seq(0, 1, by = 0.2), n, delta) {
  # delta is the age difference between grandchild and its parent (child of the carrier)
  # n is the number of years in the future to compute the probability of getting sick for the grandchild
  penetrance_delta <- c(df$Penetrance[(delta+1):length(df$Penetrance)],rep(1, delta))  
  penetrance_n <- c(df$Penetrance[(n+1):length(df$Penetrance)],rep(1, n))
  penetrances <- cbind(df$Penetrance, penetrance_n, penetrance_delta) # Create a matrix with two columns: Penetrance at age t and Penetrance at age t+n
  # Initialize empty data frame
  df_onset_expanded <- data.frame()
  for (x in xvals) {
      temp <- data.frame(Age = df$Age, x = x, n=n,
      pSICK = apply(penetrances, 1, function(k,x) pSICKn_grandchild(x, k[1], k[2], k[3]), x = x))
      df_onset_expanded <- rbind(df_onset_expanded, temp)
  }
  df_onset_expanded
}

getpMX <- function(df, xvals=seq(0, 1, by = 0.2), n) {
  # delta is the age difference between grandchild and its parent (child of the carrier)
  # n is the number of years in the future to compute the probability of getting sick for the grandchild
  penetrance_n <- c(df$Penetrance[(n+1):length(df$Penetrance)],rep(1, n))
  penetrances <- cbind(df$Penetrance, penetrance_n) # Create a matrix with two columns: Penetrance at age t and Penetrance at age t+n
  # Initialize empty data frame
  df_onset_expanded <- data.frame()
  for (x in xvals) {
      temp <- data.frame(Age = df$Age, x = x, n=n,
      pMX = apply(penetrances, 1, function(k,x) pMX(x, k[1], k[2]), x = x))
      df_onset_expanded <- rbind(df_onset_expanded, temp)
  }
  df_onset_expanded
}


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

##### FIGURE 2 #####
# Set colors
colors <- c(
  "0"   = "#4682b4",    # deep blue
  "0.2" = "#00a2c8",      # medium blue-purple
  "0.4" = "#00bfc5",      # muted magenta
  "0.6" = "#39d9ab",      # warm coral red
  "0.8" = "#9ded89",      # bright orange-red
  "1"   = "#f9f871"           # vivid red
)

## A. Probability of carrying the variant when healthy for various x - CHILDS
df_pvwh_child <- getPVWH_child(df_onset, xvals = c(0,0.2))
p1 <- ggplot(df_pvwh_child, aes(x = Age, y = PvariantWhenHealthy, color = factor(x))) +
  geom_step(linewidth=1) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey", linewidth=1) +
  labs(y = "Probability for child of carrying the c9orf72 RE variant", color = "x") +
  scale_size_manual(values = c("TRUE" = 1, "FALSE" = 0.6), guide = "none") +
  scale_color_manual(values = colors) + 
  scale_x_continuous(breaks = seq(0, max(df_onset$Age), by = 10)) +
  labs(color = expression(italic(x))) +
  theme(plot.tag = element_text(size = 16, face = "bold"), legend.position = "none") + 
  ylim(0, 0.5) +
  labs(tag = "A") + 
  annotate("label", x = 10, y = 0.48, label = "Child", size = 5, fontface = "bold", fill = "white", color = "darkgrey", label.size = 0)
p1
## B. Probability of carrying the variant when healthy for various x - GRANDCHILDS
df_pvwh_grandchild <- getPVWH_grandchild(df_onset, xvals = c(0, 0.2), delta = 25)
df_pvwh_grandchild$delta <- "25"

df_pvwh_grandchild2 <- getPVWH_grandchild(df_onset, xvals = c(0, 0.2), delta = 35)
df_pvwh_grandchild2$delta <- "35"

# Combine the two dataframes
df_grandchild_all <- rbind(df_pvwh_grandchild, df_pvwh_grandchild2)

p2 <- ggplot(df_grandchild_all, aes(x = Age, y = PvariantWhenHealthy, 
                                    color = factor(x), linetype = delta)) +
  geom_step(linewidth = 1) +
  geom_hline(yintercept = 0.25, linetype = "dashed", color = "grey", linewidth = 1) +
  labs(y = "Probability for grandchild of carrying the c9orf72 RE variant", 
       color = expression(italic(x)), 
       linetype = expression(italic(delta))) +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = c("25" = "solid", "35" = "dotted")) +
  scale_x_continuous(breaks = seq(0, max(df_onset$Age), by = 10)) +
  theme(plot.tag = element_text(size = 16, face = "bold")) +
  ylim(0, 0.5) +
  labs(tag = "B") + 
  annotate("label", x = 21, y = 0.48, label = "Grandchild", size = 5, fontface = "bold", fill = "white", color = "darkgrey", label.size = 0)

pdf("fig2.pdf", width=12)
(p1 | p2) + plot_layout(guides = "collect") & theme(legend.position = "right")
dev.off()



##### FIGURE 3 #####
# Set colors
colors2 <- c(
  "0"   = scales::gradient_n_pal(c("firebrick", "green"))(0),
  "0.2" = "orange",
  "0.4" = scales::gradient_n_pal(c("firebrick", "green"))(0.4),
  "0.6" = scales::gradient_n_pal(c("firebrick", "green"))(0.6),
  "0.8" = scales::gradient_n_pal(c("firebrick", "green"))(0.8),
  "1"   = scales::gradient_n_pal(c("firebrick", "green"))(1)
)

PrepareFigure3 <- function(N, tags=c("A","B")) {
  ## A. Probability of declaring the disease in next n years - CHILDS
  df_psick_child <- getpSICK_child(df_onset, xvals = c(0,0.2), n=N)
  p1 <- ggplot(df_psick_child, aes(x = Age, y = pSICK, color = factor(x))) +
    geom_step(linewidth=1) +
    labs(y =paste("Probability of declaring the disease in the next",N,"years"), color = expression(italic(x))) +
    scale_size_manual(values = c("TRUE" = 1, "FALSE" = 0.6), guide = "none") +
    scale_color_manual(values = colors2) + 
    scale_x_continuous(breaks = seq(0, max(df_onset$Age), by = 10)) +
    labs(color = expression(italic(x))) +
    theme(plot.tag = element_text(size = 16, face = "bold")) + 
    ylim(0, 0.5) +
    labs(tag = tags[1]) + 
    annotate("label", x = 10, y = 0.48, label = "Child", size = 5, fontface = "bold", fill = "white", color = "darkgrey", label.size = 0)
  ## A. Probability of declaring the disease in next n years - GRANDCHILDS

  df_psick_grandchild <- getpSICK_grandchild(df_onset, xvals = c(0,0.2), n=N, delta = 25)
  df_psick_grandchild$delta <- "25"
  df_psick_grandchild2 <- getpSICK_grandchild(df_onset, xvals = c(0,0.2), n=N, delta=35)
  df_psick_grandchild2$delta <- "35"
  df_psick_grandchild_all <- rbind(df_psick_grandchild, df_psick_grandchild2)
  p2 <- ggplot(df_psick_grandchild_all, aes(x = Age, y = pSICK, color = factor(x), linetype=delta)) +
    geom_step(linewidth=1) +
    labs(y = paste("Probability of declaring the disease in the next",N,"years"), color = expression(italic(x)), linetype = expression(italic(delta))) +
    scale_size_manual(values = c("TRUE" = 1, "FALSE" = 0.6), guide = "none") +
    scale_color_manual(values = colors2) + 
    scale_x_continuous(breaks = seq(0, max(df_onset$Age), by = 10)) +
    labs(color = expression(italic(x))) +
    theme(plot.tag = element_text(size = 16, face = "bold")) + 
    ylim(0, 0.5) +
    labs(tag = tags[2]) + 
    annotate("label", x = 21, y = 0.48, label = "Grandchild", size = 5, fontface = "bold", fill = "white", color = "darkgrey", label.size = 0)

  return(list(p1 = p1, p2 = p2))
}

row1 <- PrepareFigure3(10, tags=c("A","B"))
row2 <- PrepareFigure3(20, tags=c("C","D"))

pdf("fig3.pdf", width=12, height = 12)
((row1$p1 | row1$p2)/(row2$p1 | row2$p2)) + plot_layout(guides = "collect") & theme(legend.position = "right")
dev.off()



#######################
## OOOOOOLLD PLOTS ####
#######################




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
  scale_color_manual(values = colors) + 
  scale_x_continuous(breaks = seq(0, max(df_onset$Age), by = 10)) +
  labs(color = expression(italic(x))) +
  theme(plot.tag = element_text(size = 16, face = "bold"), legend.position = "none") + 
  labs(tag = "A")

## 2B: proba to die of the disease in the next n years
p2 <- ggplot(df_onset_expanded, aes(x = Age, y = pSICK, color = factor(x))) +
  geom_step(aes(size = (x == "0"))) +
  labs(y = "Probability of developing the disease in the future", color = "x") +
  scale_size_manual(values = c("TRUE" = 1, "FALSE" = 0.6), guide = "none") +
  scale_color_manual(values = colors) + 
  scale_x_continuous(breaks = seq(0, max(df_onset$Age), by = 10)) +
  labs(color = expression(italic(x))) +
  guides(color = guide_legend(override.aes = list(size = c(1, rep(0.6, 5))))) +
  theme(plot.tag = element_text(size = 16, face = "bold")) + 
  labs(tag = "B")

# Combine and save to pdf
pdf("fig2.pdf", width=12)
(p1 | p2) + plot_layout(guides = "collect") & theme(legend.position = "right")
dev.off()

#### TEMP ####

df_onset_expanded <- getPVWH_grandchild(df_onset, delta=25)
p1 <- ggplot(df_onset_expanded, aes(x = Age, y = PvariantWhenHealthy, color = factor(x))) +
  geom_step(aes(size = (x == "0"))) +
  labs(y = "Probability of carrying the c9orf72 RE variant", color = "x") +
  scale_size_manual(values = c("TRUE" = 1, "FALSE" = 0.6), guide = "none") +
  scale_color_manual(values = colors) + 
  scale_x_continuous(breaks = seq(0, max(df_onset$Age), by = 10)) +
  labs(color = expression(italic(x))) +
  theme(plot.tag = element_text(size = 16, face = "bold"), legend.position = "none") + 
  labs(tag = "A")
p1


  
n <- 10
df_onset_expanded <- getpSICK_grandchild(df_onset, n=n, delta=80)
p2 <- ggplot(df_onset_expanded, aes(x = Age, y = pSICK, color = factor(x))) +
  geom_step(aes(size = (x == "0"))) +
  labs(y = paste("Probability of developing the disease in the next", n, "years"), color = "x") +
  scale_size_manual(values = c("TRUE" = 1, "FALSE" = 0.6), guide = "none") +
  scale_color_manual(values = colors) + 
  scale_x_continuous(breaks = seq(0, max(df_onset$Age), by = 10)) +
  labs(color = expression(italic(x))) +
  guides(color = guide_legend(override.aes = list(size = c(1, rep(0.6, 5))))) +
  theme(plot.tag = element_text(size = 16, face = "bold")) + 
  labs(tag = "B")
p2

n <- 20
df_onset_expanded <- getpSICK_grandchild(df_onset, n=n, delta=40, xvals=seq(0, 0.1, by = 0.01))
p2 <- ggplot(df_onset_expanded, aes(x = Age, y = pSICK, color = factor(x))) +
  geom_step() + xlim(25,60)
p2




DELTA <- 30 # age difference between grandchild and parent
N <-10
df_onset_expanded1 <- getpMX(df_onset, n=N) #proba d'être malade dans 10 ans SI ON A LE VARIANT
p1 <- ggplot(df_onset_expanded1, aes(x = Age, y = pMX, color = factor(x))) +
  geom_step(aes(size = (x == "0"))) +
  labs(y = "Probability of developing the disease IF carrying the c9orf72 RE variant", color = "x") +
  scale_size_manual(values = c("TRUE" = 1, "FALSE" = 0.6), guide = "none") +
  scale_color_manual(values = colors) + 
  scale_x_continuous(breaks = seq(0, max(df_onset$Age), by = 10)) +
  labs(color = expression(italic(x))) +
  theme(plot.tag = element_text(size = 16, face = "bold"), legend.position = "none") + 
  labs(tag = "A")
p1
df_onset_expanded2 <- getPVWH_grandchild(df_onset, delta=DELTA) #proba d'avoir le variant
p2 <- ggplot(df_onset_expanded2, aes(x = Age, y = PvariantWhenHealthy, color = factor(x))) +
  geom_step(aes(size = (x == "0"))) +
  labs(y = paste("Probability of carrying the c9orf72 RE variant for grandchild ; delta = ",DELTA), color = "x") +
  scale_size_manual(values = c("TRUE" = 1, "FALSE" = 0.6), guide = "none") +
  scale_color_manual(values = colors) + 
  scale_x_continuous(breaks = seq(0, max(df_onset$Age), by = 10)) +
  labs(color = expression(italic(x))) +
  theme(plot.tag = element_text(size = 16, face = "bold"), legend.position = "none") + 
  labs(tag = "A")
p2

DFNEW <- cbind(df_onset_expanded1, PvariantWhenHealthy=df_onset_expanded2$PvariantWhenHealthy)
DFNEW$pSICK <- DFNEW$pMX * DFNEW$PvariantWhenHealthy
p3 <- ggplot(DFNEW, aes(x = Age, y = pSICK, color = factor(x))) +
  geom_step(aes(size = (x == "0"))) +
  labs(y = paste("Probability of developing the disease in the next", n, "years ; delta = ",DELTA), color = "x") +
  scale_size_manual(values = c("TRUE" = 1, "FALSE" = 0.6), guide = "none") +
  scale_color_manual(values = colors) + 
  scale_x_continuous(breaks = seq(0, max(df_onset$Age), by = 10)) +
  labs(color = expression(italic(x))) +
  guides(color = guide_legend(override.aes = list(size = c(1, rep(0.6, 5))))) +
  theme(plot.tag = element_text(size = 16, face = "bold")) + 
  labs(tag = "C")
p3
p1+p2+p3



df_onset_expanded <- getPVWH_child(df_onset)
p1 <- ggplot(df_onset_expanded, aes(x = Age, y = PvariantWhenHealthy, color = factor(x))) +
  geom_step(aes(size = (x == "0"))) +
  labs(y = "Probability of carrying the c9orf72 RE variant", color = "x") +
  scale_size_manual(values = c("TRUE" = 1, "FALSE" = 0.6), guide = "none") +
  scale_color_manual(values = colors) + 
  scale_x_continuous(breaks = seq(0, max(df_onset$Age), by = 10)) +
  labs(color = expression(italic(x))) +
  theme(plot.tag = element_text(size = 16, face = "bold"), legend.position = "none") + 
  labs(tag = "A")
p1



df_onset_expanded <- getpSICK_child(df_onset, n=100)
p2 <- ggplot(df_onset_expanded, aes(x = Age, y = pSICK, color = factor(x))) +
  geom_step(aes(size = (x == "0"))) +
  labs(y = paste("Probability of developing the disease in the next", n, "years"), color = "x") +
  scale_size_manual(values = c("TRUE" = 1, "FALSE" = 0.6), guide = "none") +
  scale_color_manual(values = colors) + 
  scale_x_continuous(breaks = seq(0, max(df_onset$Age), by = 10)) +
  labs(color = expression(italic(x))) +
  guides(color = guide_legend(override.aes = list(size = c(1, rep(0.6, 5))))) +
  theme(plot.tag = element_text(size = 16, face = "bold")) + 
  labs(tag = "B")
p2

#################  
#################
### FROM APP.R  #
#################
# Define function to compute the probability of carrying the C9orf72 RE mutation for grandchildren
# equation papa 
pVWH_grandchild <- function(x, pi_t_prime, pi_t) {
  numerator <- (1 - (1 - x) * pi_t) * (1 - (1 - x) * pi_t_prime)
  denominator <- 2 * (2 - (1 - x) * pi_t) - (1 - (1 - x) * pi_t) * (1 - x) * pi_t_prime
  result <- numerator / denominator
  return(result)
}

#equation papa2
pSICK_grandchild <- function(x, pi_t_prime, pi_t) {
  pprime <- (1 - x) * pi_t
  numerator <- (1 - x) * (1 - pi_t_prime) * (1-pprime)
  denominator <- 2 * (2 - pprime) - pi_t_prime * (1 - x) * (1 - pprime) 
  result <- numerator / denominator  
  return(result)
}
##########
#################################

pVWH_grandchild2 <- function(x, pi_t_prime, pi_t) {
  numerator <- (1 - (1-x)*pi_t_prime)
  denominator <- 2 * ((2 - (1-x)*pi_t)/(1-(1-x)*pi_t)) - (1-x)*pi_t_prime
  result <- numerator / denominator
  return(result)
}

pSICK_grandchild2 <- function(x, pi_t_prime, pi_t) {
  pprime <- (1 - x) * pi_t
  numerator <- (1 - x) * (1 - pi_t_prime)
  denominator <- 2 * ((2 - pprime)/(1-pprime))-((1-x)*pi_t_prime)
  result <- numerator / denominator  
  return(result)
}



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











K <- function() {
  x <- 3 + 4
  print(x)
  x
}




my_function <- function(pi_t1, pi_t2, pi_t2n) {
  denom_common <- 2 * ((2 - pi_t1) / (1 - pi_t1))
  
  num1 <- (1 - pi_t2n) * pi_t2n / (denom_common - pi_t2n)
  num2 <- (1 - pi_t2) * pi_t2 / (denom_common - pi_t2)
  numerator <- num1 - num2
  
  denom <- 1 - ((1 - pi_t2) * pi_t2 / (denom_common - pi_t2))
  
  result <- numerator / denom
  return(result)
}

0.0061082024
0.9397905759