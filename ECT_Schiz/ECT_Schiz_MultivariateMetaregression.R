library(readxl)
library(dplyr)
library(metafor)
library(meta)

df <- read_excel("ECT_SOR_Dataset.xlsx", sheet = 1)
df <- data.frame(df)

# colnames(df)[1] <- 'Study'   # Problem with Windows UTF8 encoding

# Calculate within group effect sizes for each week

df.eff <- data.frame(matrix(ncol=8))  # Blank data frame to store the effect sizes
colnames(df.eff) <- c("Study", "Week", "g1", "v1", "g2", "v2", "n1", "n2")


studies <- unique(df$Study)       # Get list of studies & weeks to run a loop
weeks <- sort(unique(df$Week))
weeks <- weeks[-1]                # Remove Week "0" from this list


# Running the for loop to calculate all within group effect sizes

for(i in studies){
  for(j in weeks){
    
    if(is.na(df[df$Study == i & df$Week == j, "mean1"][1])) {next}
    
    # Get all the values for calculating g & v in Group 1
    m1 <- df[df$Study == i &  df$Week == 0, "mean1"] 
    m2 <- df[df$Study == i & df$Week == j, "mean1"]
    sd1 <- df[df$Study == i &  df$Week == 0, "SD1"]
    sd2 <- df[df$Study == i &  df$Week == j, "SD1"]
    n1 <-  df[df$Study == i &  df$Week == j, "n1"]
    
    # Calculate Hedge's g  
   
    g1 <- ((m1 - m2)/sqrt((sd1^2 + sd2^2)/2))*sqrt((n1-2)/(n1-1))
  
    # Calculate variance of Hedge's g
    sig2.d <- sd1^2 + sd2^2 - 2*0.5*sd1*sd2
    sig4 <- ((sd1^2 + sd2^2)/2)^2
    sdc.1 <- (g1^2*(sd1^4 + sd2^4 + 2*0.5*(sd1^2)*(sd2^2))/(8*(n1-1)*sig4))
    sdc.2 <- sig2.d/(sqrt(sig4)*(n1-1))
    
    v1 <- sdc.1 + sdc.2
    
    #------------------------------
    
    # Get all the values for calculating g & v in Group 2
    m1 <- df[df$Study == i &  df$Week == 0, "mean2"]
    m2 <- df[df$Study == i & df$Week == j, "mean2"]
    sd1 <- df[df$Study == i &  df$Week == 0, "SD2"]
    sd2 <- df[df$Study == i &  df$Week == j, "SD2"]
    n2 <- df[df$Study == i &  df$Week == j, "n2"]
    
    # Calculate Hedge's g  
    
    g2 <- ((m1 - m2)/sqrt((sd1^2 + sd2^2)/2))*sqrt((n2-2)/(n2-1))
    
    # Calculate variance of Hedge's g in Group1 
    sig2.d <- sd1^2 + sd2^2 - 2*0.5*sd1*sd2
    sig4 <- ((sd1^2 + sd2^2)/2)^2
    sdc.1 <- (g2^2*(sd1^4 + sd2^4 + 2*0.5*(sd1^2)*(sd2^2))/(8*(n2-1)*sig4))
    sdc.2 <- sig2.d/(sqrt(sig4)*(n2-1))
    
    v2 <- sdc.1 + sdc.2
    
    result <- c(i, j, g1, v1, g2, v2, n1, n2 )
    
    # Ask the loop to ignore weeks that are absent for that study
    if(length(result) == 2){next}       
    
    df.eff <- rbind(df.eff, result)
    
  }
  rm(g1, g2, v1, v2, n1, n2, m1, m2, sd1, sd2, i, j, result, sdc.1, sdc.2, sig2.d, sig4)
}


df.eff <- df.eff[-1, ]
df.eff[,-1] <- lapply(df.eff[,-1], function(x) as.numeric(x) )

# Plotting the results

library(ggplot2)

ggplot(df.eff, mapping = aes(x=Week, y= g1, color = Study, group = Study)) + 
  geom_point() + 
  geom_line() + 
  geom_errorbar(aes(ymin = g1 - v1, ymax = g1 + v1)) +
  scale_x_continuous(breaks = 1:28) + ggtitle("ECT+AP")

ggplot(df.eff, mapping = aes(x=Week, y= g2, color = Study, group = Study)) + 
  geom_point() + 
  geom_line() + 
  geom_errorbar(aes(ymin = g2 - v2, ymax = g2 + v2)) +
  scale_x_continuous(breaks = 1:28) + ggtitle("AP only")



# Running the between-group meta-analysis at each week

rmas <- list()

rma.wkly <- function(x) {metacont(n.e=n1, mean.e = g1, sd.e = sqrt(v1), 
                                  n.c = n2, mean.c = g2, sd.c = sqrt(v2),
                                  studlab = Study,
                                  data = df.eff[df.eff$Week == x,],
                                  sm = "SMD", method.smd= "Hedges", fixed = FALSE, tau.common = FALSE)}


rma.wkly.forest <- function(x) {meta::forest(x, label.right = "Favors ECT+AP", label.left= "Favors AP Only", 
                                       leftlabs =  c("Study", "N", "SMC", "SD",  "N", "SMC", "SD"),
                                       digits = 2, digits.se=2, hetstat= FALSE) 
                                 grid.text(paste("Week", x[["data"]][["Week"]]), .5, .9, gp=gpar(cex=2))}

for(i in 1:8) {rmas[[paste("Week", i)]] <- rma.wkly(i)}

for(i in 1:8) {rma.wkly.forest(rmas[[i]])}


# Within-groups effect sizes in both arms at each week:

df.wg <- data.frame(matrix(ncol=9))
colnames(df.wg) <- c("Week", "gp1", "sep1", "gp2", "sep2", "np1", "np2", "k1", "k2")

for(j in weeks) { x1= metamean(n = n1, mean = g1, sd = sqrt(v1), 
                               studlab = Study,data = df.eff[df.eff$Week ==j,])
                  x2= metamean(n = n2, mean = g2, sd = sqrt(v2), 
                               studlab = Study,data = df.eff[df.eff$Week ==j,])
                 np1 <- sum(x1$n) %>% round
                 np2 <- sum(x2$n) %>% round
  
                 gp1 <- x1$TE.random %>% round(digits = 3)
                 gp2 <- x2$TE.random %>% round(digits = 3)
                 
                 sep1 <- x1$seTE.random %>% round(digits = 3)
                 sep2 <- x2$seTE.random %>% round(digits = 3)
                 
                 k1 <- x1$k.study %>% round
                 k2 <- x2$k.study %>% round
                 
                out <-  c(j, gp1, sep1, gp2, sep2, np1, np2, k1, k2)
                df.wg <- rbind(df.wg, out)
                
                rm(i, j, gp1, sep1, gp2, sep2, np1, np2, k1, k2)
                
                }

df.wg <- df.wg[-1,]


# Plotting this

df.wg.long <- reshape(df.wg, varying= list(c("gp1", "gp2"), 
                                           c("sep1", "sep2"), 
                                           c("np1", "np2"), 
                                           c("k1", "k2")),
                      v.names = c("g", "se", "n", "k"),
                      timevar = "Group",
                      idvar = "Week",
                      times = c(1,2),
                      direction = "long" ) 


df.wg.long$Group <- factor(df.wg.long$Group, labels = c("ECT+AP", "AP Only"))

ggplot(df.wg.long, mapping=aes(x=Week, y=g, ymin = g-se, ymax=g+se, group = Group)) +
  geom_point() +
  geom_errorbar(width=0.3) +
  geom_line() +
  facet_wrap(~Group, ncol=1) +
  scale_x_continuous(name="Weeks", breaks = seq(0,28,2)) +
  scale_y_continuous(name = "Pooled Hedge's g (SE)", limits = c(0,6))

# Plot only till 8 weeks

ggplot(df.wg.long[df.wg.long$Week <10,], 
       mapping=aes(x=Week, y=g, ymin = g-se, ymax=g+se, group = Group)) +
  geom_point() +
  geom_errorbar(width=0.3) +
  geom_line() +
  facet_wrap(~Group, ncol=1) +
  scale_x_continuous(name="Weeks", breaks = 0:8) +
  scale_y_continuous(name = "Pooled Hedge's g (SE)", limits = c(0,6))

  
# Running the multi-variate meta-analysis

df.eff$ID <- 1:nrow(df.eff)

df.long <- reshape(df.eff, varying = list(c("g1","g2"),
                                          c("v1","v2"),
                                          c("n1","n2")),
                   v.names = c("g", "v", "n"),
                   timevar = "Group",
                   idvar = "ID",
                   times = c(1,2),
                   direction="long")

df.long$Group <- factor(df.long$Group, labels = c("ECT+AP", "AP only"))

df.long <- df.long[,-3]

# Multivariate meta-analysis running it
rma.mv(g ~ Group*Week, v,  random = ~ Week|Study,
       struct ="CAR", data=df.long)

rma.mv(g ~ Group*Week, v,  random = ~ Week|Study,
       struct ="CAR", data=df.long[df.long$Week <10,])
