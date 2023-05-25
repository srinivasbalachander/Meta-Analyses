---
title: "OCD Comorbidity Meta-analysis: Across the Lifespan"
author: "Sharma E et al."
date: "26/07/2021"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)
```

```{r packages, include=FALSE}
library(dplyr)
library(kableExtra)
library(metafor)
library(meta)
library(weightr)
```

```{r plotsizing, include=FALSE}

# A funciton to customize image sizes within code chunks. As some forests plots are really long and some arent/

subchunkify <- function(g, fig_height=15, fig_width=10) {g_deparsed <- paste0(deparse(function() {g}), collapse = '')   
sub_chunk <- paste0("`","``{r sub_chunk_", floor(runif(1) * 10000), ", fig.height=", fig_height, ", fig.width=", fig_width, ", echo=FALSE}",
  "\n(", g_deparsed, ")()","\n`","``")
  cat(trimws(knitr::knit(text = knitr::knit_expand(text = sub_chunk), quiet = TRUE)))
}

# Another function for easy rounding 
round2 <- function(x) {round(x, digits=3)} #Make a small function to round off digits

```

```{r input, include=FALSE}

df <- read.csv("Eesha OCD Comorb.csv")

df2 <- df %>% filter(Recruitment_Type != "Community") %>% select(StudyID, Country, Design, Total_N, Age_Mean, Age_Group, AAO_Mean, YBOCST_Mean, Instrument, Gender_Male_Perc, 29:218)   # Select only the important columns

df2$Age_Group <- factor(df2$Age_Group, levels = c("Pediatric", "Adult"))

df.var <- read.csv("Var_Names.csv") # this is a file with better labels for each disorder
df.var$No <- 1:nrow(df.var)

df.var2 <- read.csv("Var_Included.csv") # similar file, but includes only the disorders we're including in this paper
df.var2$No <- 1:nrow(df.var2)

v1 <- sapply(df.var$Var, as.character) # a character vector, useful for running for loops
v2 <- sapply(df.var2$Var, as.character)

```

```{r forestfunc, echo=F}

# Function for making forests & funnels. 

forestfunc <- function(x) {temp <- df2  %>% select(StudyID, Country, Total_N, Age_Mean, AAO_Mean, YBOCST_Mean, Age_Group, x)
                          temp <- temp %>% filter(!is.na(get(x)))
                          writeLines(paste0("-----------------------------------------------------------------------------", "\n"))
                          
                          if (nrow(temp) >4) {
                            cat(paste0("<h4>" ,"Suppl Figure 2.", df.var2[df.var2$Var == x, "No"],".1 Forest Plot for ",
                                                               df.var2[df.var2$Var == x, "Proper_Name"], "</h4>"))  # Title for the forest plot
                                 ies.darc <- escalc(xi = get(x), ni = Total_N, measure = "PFT", data=temp)          # Run the initial RMA to figure out outliers
                                 pes.darc <- rma(yi, vi, data=ies.darc, method="DL", weighted=TRUE)
                                 infp <- influence(pes.darc)
                            pes.forest <- metaprop(get(x), Total_N, StudyID, data=temp[!infp$is.infl,], byvar=Age_Group,
                                                   sm="PFT", method.ci = "NAsm", method.tau = "DL",
                                                incr = 0.5,allincr = FALSE, addincr = FALSE, title="")
                            subchunkify(meta::forest(pes.forest, comb.fixed=FALSE), fig_height = (nrow(temp)*0.24+3.5), fig_width = 10)   # Forest plot
                            cat(paste0("<h4>" ,"Suppl Figure 2.", df.var2[df.var2$Var == x, "No"],".2 Funnel Plot for ",
                                                               df.var2[df.var2$Var == x, "Proper_Name"], "</h4>"))
                            
                            egg.darc <- regtest(pes.darc)     # Egger's test
                            subchunkify(metafor::funnel.rma(pes.darc, atransf = transf.ipft.hm, targs=list(ni=temp$Total_N), 
                                                            sub= paste0("Egger's Assymetry Test Z= ", round2(egg.darc$zval), "; ", 
                                                                         "p= ", round2(egg.darc$pval)),
                                                            cex.sub=0.7), 
                                        fig_height= 5, fig_width = 7)                                      # Funnel Plot
                          writeLines(paste0("</h3>------------------------------------------------------------------------------</h3>", "\n"))
                                              }
                                 else { writeLines(paste0("Less than 5 studies here, so no meta-analysis, no forest plot!", "\n"))
                                        writeLines(paste0("</h3> ------------------------------------------------------------------------- </h3>"))
                                            cat("\n")}
                            }
```


``` {r dummy, include=FALSE}
# This is just a dummy piece of code to take care of some annoying error messages!
forestfunc("Any_AxisI_Life")

```

# ------------------------------------------------------------------------------------------

## Forest plots for each disorder

Method: Random-effects meta-analysis of proportions with double-arcsine transformation

``` {r supptables, results='asis', R.options = list(width = 100)}
# Making all the tables, the easiest possible way!
for (i in v2) {forestfunc(i)}
```
