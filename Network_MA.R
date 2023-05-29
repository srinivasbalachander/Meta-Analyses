library(netmeta)
library(readxl)

df <- read.csv("UMESHNMAWIDE.csv")

labels <- read_excel("THERAPYLABELNMA.xlsx")


df.nm <- pairwise(treat = list(T.1, T.2, T.3, T.4, T.5),
               n = list(N.1, N.2, N.3, N.4, N.5),
               mean = list(Mean.1, Mean.2, Mean.3, Mean.4, Mean.5), 
               sd = list(SD.1, SD.2, SD.3, SD.4, SD.5),
               data = df, studlab = Study)

df.nm <- df.nm[,1:11]


df.nm[c("treat1", "treat2")] <- lapply(df.nm[c("treat1", "treat2")],
                                       function(x) labels$Label[match(x, labels$Number)])


m.netmeta <- netmeta(TE = TE,
                     seTE = seTE,
                     treat1 = treat1,
                     treat2 = treat2,
                     studlab = studlab,
                     data = df.nm,
                     sm = "MD",
                     fixed = FALSE,
                     random = TRUE,
                     reference.group = "Waitlist",
                     details.chkmultiarm = TRUE,
                     sep.trts = " vs ")
summary(m.netmeta)



netgraph(m.netmeta)


library(rgl)
netgraph(m.netmeta, dim = "3d")


d.evidence <- direct.evidence.plot(m.netmeta)
plot(d.evidence)

netrank(m.netmeta, small.values = "good")

funnel(m.netmeta, order = labels$Label, linreg = TRUE) 

netsplit(m.netmeta) %>% forest()



forest(m.netmeta, 
       reference.group = "Waitlist",
       sortvar = TE,
       xlim = c(-21, 13),
       smlab = paste("Psychotherapies for OCD"),
       drop.reference.group = TRUE,
       label.left = "Favors Intervention",
       label.right = "Favors Care As Usual")

