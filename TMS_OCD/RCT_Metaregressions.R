library(metafor)
library(meta)
library(dplyr)


df <- read.csv("Analysis_DepAnx.csv")

names(df)[1] <- "StudyID"

# Doing the meta-analysis with pre-post differences

df.dep <- df %>% select(StudyID, Target, 
                        Mean_Dep_Active_Baseline, SD_Dep_Active_Baseline,
                        Mean_Dep_Active_End, SD_Dep_Active_End, N_Active,
                        Mean_Dep_Sham_Baseline, SD_Dep_Sham_Baseline,
                        Mean_Dep_Sham_End, SD_Dep_Sham_End, N_Sham) %>% na.exclude

df.anx <- df %>% select(StudyID, Target, 
                        Mean_Anx_Active_Baseline, SD_Anx_Active_Baseline,
                        Mean_Anx_Active_End, SD_Anx_Active_End, N_Active,
                        Mean_Anx_Sham_Baseline, SD_Anx_Sham_Baseline,
                        Mean_Anx_Sham_End, SD_Anx_Sham_End, N_Sham) %>% na.exclude

df.YBOCS <- df %>% select(StudyID, Target,
                          Mean_YBOCS_Active_Baseline, SD_YBOCS_Active_Baseline,
                          Mean_YBOCS_Active_End, SD_YBOCS_Active_End, N_Active,
                          Mean_YBOCS_Sham_Baseline, SD_YBOCS_Sham_Baseline,
                          Mean_YBOCS_Sham_End, SD_YBOCS_Sham_End, N_Sham) %>% na.exclude

# Calculating depression mean changes

df.dep$m.dep.ch.act <- df.dep$Mean_Dep_Active_Baseline - df.dep$Mean_Dep_Active_End
df.dep$sd.dep.ch.act <- sqrt(df.dep$SD_Dep_Active_Baseline^2 + df.dep$SD_Dep_Active_End^2 + 
                               2*0.5*df.dep$SD_Dep_Active_Baseline*df.dep$SD_Dep_Active_End)
df.dep$m.dep.ch.shm <- df.dep$Mean_Dep_Sham_Baseline - df.dep$Mean_Dep_Sham_End
df.dep$sd.dep.ch.shm <-sqrt(df.dep$SD_Dep_Sham_Baseline^2 + df.dep$SD_Dep_Sham_End^2 + 
                              2*0.5*df.dep$SD_Dep_Sham_Baseline*df.dep$SD_Dep_Sham_End)

# Calculating anxiety mean changes

df.anx$m.anx.ch.act <- df.anx$Mean_Anx_Active_Baseline - df.anx$Mean_Anx_Active_End
df.anx$sd.anx.ch.act <- sqrt(df.anx$SD_Anx_Active_Baseline^2 + df.anx$SD_Anx_Active_End^2 + 
                               2*0.5*df.anx$SD_Anx_Active_Baseline*df.anx$SD_Anx_Active_End)
df.anx$m.anx.ch.shm <- df.anx$Mean_Anx_Sham_Baseline - df.anx$Mean_Anx_Sham_End
df.anx$sd.anx.ch.shm <-sqrt(df.anx$SD_Anx_Sham_Baseline^2 + df.anx$SD_Anx_Sham_End^2 + 
                              2*0.5*df.anx$SD_Anx_Sham_Baseline*df.anx$SD_Anx_Sham_End)
# Calculating YBOCS mean changes

df.YBOCS$m.YBOCS.ch.act <- df.YBOCS$Mean_YBOCS_Active_Baseline - df.YBOCS$Mean_YBOCS_Active_End
df.YBOCS$sd.YBOCS.ch.act <- sqrt(df.YBOCS$SD_YBOCS_Active_Baseline^2 + df.YBOCS$SD_YBOCS_Active_End^2 + 
                                   2*0.5*df.YBOCS$SD_YBOCS_Active_Baseline*df.YBOCS$SD_YBOCS_Active_End)

df.YBOCS$m.YBOCS.ch.shm <- df.YBOCS$Mean_YBOCS_Sham_Baseline - df.YBOCS$Mean_YBOCS_Sham_End
df.YBOCS$sd.YBOCS.ch.shm <-sqrt(df.YBOCS$SD_YBOCS_Sham_Baseline^2 + df.YBOCS$SD_YBOCS_Sham_End^2 + 
                                  2*0.5*df.YBOCS$SD_YBOCS_Sham_Baseline*df.YBOCS$SD_YBOCS_Sham_End)

# Calculating SMCRs (Pre-Post Within Group Effect Sizes)

df.dep <- escalc(measure="SMCR", m1i=Mean_Dep_Active_End, m2i=Mean_Dep_Active_Baseline, 
                 sd1i=SD_Dep_Active_Baseline, ni=N_Active,
                 ri= rep(0.5, times= nrow(df.dep)), data=df.dep, var.names = c("es.dep.ch.act", "v.dep.ch.act"))

df.dep <- escalc(measure="SMCR", m1i=Mean_Dep_Sham_End, m2i=Mean_Dep_Sham_Baseline, 
                 sd1i=SD_Dep_Sham_Baseline, ni=N_Sham,
                 ri= rep(0.5, times= nrow(df.dep)), data=df.dep, var.names = c("es.dep.ch.shm", "v.dep.ch.shm"))

df.anx <- escalc(measure="SMCR", m1i=Mean_Anx_Active_End, m2i=Mean_Anx_Active_Baseline, 
                 sd1i=SD_Anx_Active_Baseline, ni=N_Active,
                 ri= rep(0.5, times= nrow(df.anx)), data=df.anx, var.names = c("es.anx.ch.act", "v.anx.ch.act"))

df.anx <- escalc(measure="SMCR", m1i=Mean_Anx_Sham_End, m2i=Mean_Anx_Sham_Baseline, 
                 sd1i=SD_Anx_Sham_Baseline, ni=N_Sham, 
                 ri= rep(0.5, times= nrow(df.anx)), data=df.anx, var.names = c("es.anx.ch.shm", "v.anx.ch.shm"))

df.YBOCS <- escalc(measure="SMCR", m1i=Mean_YBOCS_Active_End, m2i=Mean_YBOCS_Active_Baseline, 
                   sd1i=SD_YBOCS_Active_Baseline, ni=N_Active,
                   ri= rep(0.5, times= nrow(df.YBOCS)), data=df.YBOCS, var.names = c("es.YBOCS.ch.act", "v.YBOCS.ch.act"))

df.YBOCS <- escalc(measure="SMCR", m1i=Mean_YBOCS_Sham_End, m2i=Mean_YBOCS_Sham_Baseline, 
                   sd1i=SD_YBOCS_Sham_Baseline, ni=N_Sham, 
                   ri= rep(0.5, times= nrow(df.YBOCS)), data=df.YBOCS, var.names = c("es.YBOCS.ch.shm", "v.YBOCS.ch.shm"))

# Random-effects meta-analysis of mean changes, using the Target as sub-group

rma_dep <-  metacont(n.e=N_Active, mean.e = m.dep.ch.act ,sd.e = sd.dep.ch.act , 
                      n.c = N_Sham, mean.c = m.dep.ch.shm, sd.c = sd.dep.ch.shm, 
                      studlab = StudyID, byvar = Target,
                      sm = "SMD", method.smd= "Hedges", data=df.dep, comb.fixed = FALSE)

rma_anx <-  metacont(n.e=N_Active, mean.e = m.anx.ch.act ,sd.e = sd.anx.ch.act , 
                      n.c = N_Sham, mean.c = m.anx.ch.shm, sd.c = sd.anx.ch.shm, 
                      studlab = StudyID, byvar = Target,
                      sm = "SMD", method.smd= "Hedges", data=df.anx, comb.fixed = FALSE)

rma_YBOCS <-  metacont(n.e=N_Active, mean.e = m.YBOCS.ch.act ,sd.e = sd.YBOCS.ch.act , 
                        n.c = N_Sham, mean.c = m.YBOCS.ch.shm, sd.c = sd.YBOCS.ch.shm, 
                        studlab = StudyID, byvar = Target,
                        sm = "SMD", method.smd= "Hedges", data=df.YBOCS, comb.fixed = FALSE, tau.common = FALSE)

meta::forest(rma_dep, label.right = "Favors Active", label.left= "Favors Sham", title(main = "Depression")) 
meta::forest(rma_anx, label.right = "Favors Active", label.left= "Favors Sham", title(main = "Anxiety"))
meta::forest(rma_YBOCS, label.right = "Favors Active", label.left= "Favors Sham", title(main = "YBOCS"))


meta::funnel(rma_dep)
meta::funnel(rma_anx)
meta::funnel(rma_YBOCS)

# Meta-Regression for anxeity and depression predicting YBOCS change

df.YBOCS.dep <- merge(df.YBOCS, select(df.dep, StudyID, contains("dep.ch")), by="StudyID") 
df.YBOCS.anx <- merge(df.YBOCS, select(df.anx, StudyID, contains("anx.ch")), by="StudyID") 

df.YBOCS.dep$es.ch.dep <-  df.YBOCS.dep$es.dep.ch.shm - df.YBOCS.dep$es.dep.ch.act
df.YBOCS.anx$es.ch.anx <-  df.YBOCS.anx$es.anx.ch.shm - df.YBOCS.anx$es.anx.ch.act


rma_YBOCS.anx <- metacont(n.e=N_Active, mean.e = m.YBOCS.ch.act ,sd.e = sd.YBOCS.ch.act , 
                        n.c = N_Sham, mean.c = m.YBOCS.ch.shm, sd.c = sd.YBOCS.ch.shm, 
                        studlab = StudyID, byvar = Target,
                        sm = "SMD", method.smd= "Hedges", data=df.YBOCS.anx, comb.fixed = FALSE, tau.common = FALSE)

rma_YBOCS.dep <- metacont(n.e=N_Active, mean.e = m.YBOCS.ch.act ,sd.e = sd.YBOCS.ch.act , 
                        n.c = N_Sham, mean.c = m.YBOCS.ch.shm, sd.c = sd.YBOCS.ch.shm, 
                        studlab = StudyID, byvar = Target,
                        sm = "SMD", method.smd= "Hedges", data=df.YBOCS.dep, comb.fixed = FALSE, tau.common = FALSE)


rma_YBOCS.depreg <- metareg(rma_YBOCS.dep, ~ es.ch.dep)
rma_YBOCS.anxreg <- metareg(rma_YBOCS.anx, ~ es.ch.anx)

rma_YBOCS4.anxreg
rma_YBOCS5.depreg

bubble(rma_YBOCS4.anxreg, studlab = TRUE, xlab = "Effect Size of Anxiety Change", ylab="Effect Size of YBOCS Change")
bubble(rma_YBOCS5.depreg, studlab = TRUE, xlab = "Effect Size of Depression Change", ylab="Effect Size of YBOCS Change")
