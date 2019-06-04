# Script to measure size of aggregations around FADs

# 0. Description  ----------------

# There are several objectives in this study:
# (1) Measure and correct atwarth distorsion of multibeam sonar
# (2) Compare horizontal and vertical diameters of aggregations
#     Is there any consistent pattern?
# (3) Correlate abunance with chatches


# The data we are using are:
#   horizontal sonar, 
#   vertical echosounder EK60 @ 120 kHz, and:
#   catches per species

# The sonar is a Simrad SN90? with 32 beams of 3.75º each
# According to Missund (1993), when measuring school sizes with sonar
# the alongship beam correction should be:
#   Corr.beam = Uncorr.beam - (1500*Pulse*1e-6)/2 
# and the atwarthship correction should be:
#   Corr.ring = Uncorr.ring - 2*n*(2/3)*Range*tan(3.75*pi/(2*180))
# where n is the number of beams occupied by the school
# and Range is the distance to the center of the school

# 1. Load libraries and input data --------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(broom)

# read data files:
# (1) dimensiones - sonar-based school diameters 
# (2) vert - echosounder-based school heights
# (3) lances - info about catches
# (4) acous - echosounder acoustic echointegration 

# (1) dimensiones -  sonar-based school diameters 
# done with Profos by Udane 
# It measures simultaneously longitudinal and transversal diameters of 
# the same aggregation around a FAD
# it contains multiple measures per set, so 
# needs to be summarised
dimensiones <- read.table(file = "datos/Dimensions-pp.txt", sep = "\t", 
                          header = T, dec = ".")  %>% 
  rename(distance = Diancia.del.Banco..m.)
str(dimensiones)

# Relevant variables:
# Self-explanatory: set, Date (mmddyyyy), Time (hhmm)
# Trans.tilt: transducer tilt angle (º) (0 for horizontal)
# AlongBeamSize: longitudinal measure of the school diameter (m)
# AlongRingSize: transversal measure of the school diameter (m)
# Dudoso: especially difficult or compromised diameter measures (1/0)
# Pulso: pulse duration (micro s)
# Rango: maximum detection range of the sonar (m)
# Frecuency: operating frequency of the sonar (kHz)
# distance: distance to the center of the school (m)


# (2) vert - echosounder-based school heights
# output of the lines drawn in Echoview by Iñaki
# indicating the approximate initial and final depth of each school
# as recorded by the vertical echosounder onboard the logistic boat
# It contains multiple depths per set, so 
# we have to summarise befor merging to other data
vert <- read.table("datos/profundidades.txt", header = T, sep = "\t")
str(vert)

# Relevant variables:
# set (id), Zsup, SD_Zsup, Zfin, SD_Zfin (m)

# (3) lances - info about catches
lances <- read.table("datos/Lances.txt", header = T, sep = "\t")
str(lances)  

# Relevant variables:
# Lance, Fecha, Hora
# Captura.Total, skj, bet, yft (t)

# (4) acous - echosounder acoustic echointegration 
# vertical acoustic density and NASC at 120 kHz
# it is 120 instead of 38 because the 38 failed at some point of the survey
# it is echointegrated between the ranges of the Echoview lines in "vert"
acous <- read.table("datos/Ecointegra120xlance.txt", header = T, sep = "\t")
acous <- acous[c(1,5:17, 21:23)]
acous <- acous[-c(5:9)]
str(acous)

# Relevant variables:
# Set, Date_M (yyyymmdd), Time_M (hh:mm:ss.ssss; factor)
# NASC, Sv_mean

## 1. Correction of the atwarth dimension  -----------------
# we measure the diameter of the aggregation along beam and along ring
# at complete round trajectories surroundng the FAD during purse seining ops.
# then we compare the corrected along-beam diameter estimations
# with the along-ring ones
# and try to estimate correction factors for the latter

summary(dimensiones)

# we apply Misund's corrections:
dimensiones <- dimensiones %>%
  mutate(n = case_when(
          Rango == 400 ~ 15, 
          Rango == 500 ~ 14, 
          Rango == 600 ~ 13
          ),
          beam = AlongBeamSize - (1500*Pulso*1e-6)/2, 
          ring = AlongRingSize - 2*n*(2/3)*distance*tan(3.75*pi/(2*180)))
dimensiones %>%
  filter(Dudoso != 1) %>%
  ggplot(aes(x = beam, y = ring)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = "lm") +
  theme_bw() +
  facet_grid(Rango ~ Frecuency)
  
modelo <- lm(data = dimensiones, beam ~ ring)
summary(modelo)

lmod <- function(x) lm(data = x, beam ~ AlongRingSize)

dim.l <- split(x = dimensiones, f = dimensiones$Frecuency)
mod.l <- map(.x = dim.l, .f = lmod)
map(.x = mod.l, .f = summary)
# significant*** correlation between along and atwarth diameters for all freqs 
# the best R2 for 75 kHz (61%) followed by 82 kHz (52%); the rest, well below 50% 
# Estimated correction for atwarth diameters
# corrected.atwarth.diam = 0.18*atwarth.diam - 14.5 m

dim.l2 <- split(x = dimensiones, f = dimensiones$Rango)
mod.l2 <- map(.x = dim.l2, .f = lmod)
map(.x = mod.l2, .f = summary)

# Check try the following idea: correct the atwarth diameter estimation with 
# the formula by Misund (1993). Then, see if it corrects completely
# according to the along beam estimation. If not, apply the result of 
# the linear regression obtained here.

# boxpplot comparing beam vs ring school diam
dimensiones %>%  
  gather(beam:ring, key = "type", value = "diam") %>% 
  ggplot(aes(y = diam)) + 

  geom_boxplot(aes(fill = type)) 

  geom_boxplot(aes(fill = type)) +
  facet_wrap(~ set)
  
  
  # boxpplot comparing beam vs ring school diam
  dimensiones %>%  
    gather(beam:ring, key = "type", value = "diam") %>% 
    ggplot(aes(y = diam)) + 
  geom_boxplot(aes(fill = type)) +
  geom_boxplot(aes(fill = type)) +
    facet_wrap(~ set)
  

# scatterplot of delta-diam against set number
dimensiones %>%  
  mutate(delta = ring - beam) %>% 
  group_by(set) %>% 
  summarise(delta = mean(delta), pulse = mean(Pulso)) %>% 
  ggplot(aes(y = delta, x = set)) + 
  geom_point(aes(colour = pulse)) +
  geom_hline(yintercept = 0, linetype = 2) 

## 2. Vertical vs horizontal diameters -------------------------
# we compare vert vs horiz to see whether there is a stable relationship between them

## 2.1 Summarise vertical diameters per set ---------
vert <- vert %>%
  mutate(delta.z = Zfin - Zsup) 
vert %>% summary()
vert %>% group_by(set) %>%
  summarise(mean(delta.z))

vert.sum <- vert %>% 
  mutate(delta.z = Zfin - Zsup, 
         sd.delta.z = SD_Zsup + SD_Zfin) %>%
  group_by(set) %>%
  summarise(delta.z = mean(delta.z, na.rm = T), 
            sd.z = mean(sd.delta.z, na.rm = T)) %>%
  mutate(cv.z = sd.z/delta.z)
summary(vert.sum$delta.z)

## 2.2 Summarise horizontal diameters per set -------
horiz.sum <- dimensiones %>%
  group_by(set) %>%
  summarise(diam = mean(beam, na.rm = T), 
            sd.diam = sd(beam, na.rm = T)) %>%
  mutate(cv.diam = sd.diam/diam)

## 2.3 Combine and compare vertical and horizontal info  --------------
geom <- merge(vert.sum, horiz.sum, by = "set", all = F)

geom %>%
  ggplot(aes(x = delta.z, y = diam)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw()

# unfortunately (or fortunately) the diameter seems pretty constant
# independent from the vertical extension of the aggregation
# According to these data, for vertical extensions between 20 and 70 m
# in the Atlantic Ocean the average horizontal diameter would be constant = 50 m
# Unless we find a relationship with sv or nasc...


## 3. Predicting catches with only vertical data ----------------
lances.geom <- merge(geom, lances, by.x = "set", by.y = "Lance")

lances.geom <- merge(lances.geom, acous, by.x = "set", by.y = "Set")
lances.geom <- lances.geom[c(1:8, 10, 15:24, 27:32)] 


panel.cor <- function(x, y, digits=2, prefix="", cex.cor) 
{
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  r <- abs(cor(x, y)) 
  txt <- format(c(r, 0.123456789), digits = digits)[1] 
  txt <- paste(prefix, txt, sep = "") 
  if (missing(cex.cor)) cex <- 0.8/strwidth(txt) 
  
  test <- cor.test(x,y) 
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  
  text(0.5, 0.5, txt, cex = cex * r) 
  text(.8, .8, Signif, cex = cex, col = 2) 
}
x11()
pairs(lances.geom, lower.panel = panel.smooth, upper.panel = panel.cor)
# I don't see anything interesting
