# Umesh's Network Meta-Analysis

library(readxl)
library(netmeta)
library(dmetar)
library(rgl)
library(colorspace)


rm(list=ls())

# Read the dataset and the labels
df.wide <- read.csv("UMESHNMAWIDE.csv")
labels <- read_excel("THERAPYLABELNMA.xlsx")

# Prepare the data for frequentist NMA
df.nm <- pairwise(list(T.1, T.2, T.3, T.4, T.5),
               n = list(N.1, N.2, N.3, N.4, N.5),
               mean = list(Mean.1, Mean.2, Mean.3, Mean.4, Mean.5), 
               sd = list(SD.1, SD.2, SD.3, SD.4, SD.5),
               data = df.wide, 
               studlab = Study)
 
df.nm <- df.nm[, 1:11]

# Relabeling the treatment arms
df.nm[c("treat1", "treat2")] <- lapply(df.nm[c("treat1", "treat2")],
                                       function(x) labels$Label[match(x, labels$Number)])

# Look at number of studies for each intervention/node
k.nodes <- data.frame(table(c(df.nm$treat1, df.nm$treat2)))
colnames(k.nodes) <- c("Intervention", "Number of Studies")
k.nodes

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
         iterate = FALSE
         )


?netgraph.netmeta

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
