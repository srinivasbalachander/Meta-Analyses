# Umesh's Network Meta-Analysis

library(googlesheets4)
library(readxl)
library(netmeta)
library(dmetar)
library(rgl)
library(colorspace)

rm(list=ls())

gs4_auth()

# Read the dataset and the labels

df.long <- read_sheet("https://docs.google.com/spreadsheets/d/1qDLFM_8M6YJCjwqCdSxT4kpZ1Lgi0iORP9lCEOU0ZP0/edit?usp=sharing", 
                      sheet = 1)

# Prepare the data for frequentist NMA
df.nm <- pairwise(treat = Treatment,
                  n = N,
                  mean = Mean, 
                  sd = SD,
                  data = df.long, 
                  studlab = Study)

df.nm <- df.nm[, 1:11]


# Model Fitting

m.netmeta <- netmeta(TE = TE,
                     seTE = seTE,
                     treat1 = treat1,
                     treat2 = treat2,
                     n1 = n1,
                     n2 = n2,
                     studlab = studlab,
                     data = df.nm,
                     sm = "MD",
                     fixed = FALSE,
                     random = TRUE,
                     reference.group = "Waitlist",
                     details.chkmultiarm = TRUE,
                     sep.trts = " vs ")

summary(m.netmeta)

data.frame(m.netmeta$k.trts)

decomp.design(m.netmeta)


# 2-D Network Graph

netgraph(m.netmeta, 
         plastic = FALSE, # Change to false and all lines will look solid
         col = "darkblue",
         thickness = "number.of.studies",  # Can also be "w.random", or "equal" or "se.random"
         rescale.thickness = FALSE,
         number.of.studies = FALSE,
         points = TRUE,  # Show all nodes
         col.points = "maroon",
         cex.points = n.trts,  # Scale node size by number of participants within each treatment
         rescale.pointsize = FALSE,
         points.max = 12,
         points.min = 4,
         # srt.labels = "orthogonal",  # In case you want to rotate the labels
         iterate = FALSE)


netgraph(m.netmeta, dim = "3d",
         points = TRUE)

# P Scores for Ranking
netrank(m.netmeta, small.values = "good")

# Forest plot of direct + indirect (total network) effect estimates
forest(m.netmeta, 
       reference.group = "Waitlist",
       sortvar = TE,
       xlim = c(-25, 15),
       smlab = paste("Psychotherapies for \n",
                     "OCD"),
       drop.reference.group = TRUE,
       label.left = "Favors Intervention",
       label.right = "Favors Care As Usual")


# Direct versus Indirect Evidence Plot
d.evidence <- direct.evidence.plot(m.netmeta, random = TRUE)

plot(d.evidence)
