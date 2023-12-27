# Author : Orson Mestanza

library(tidyverse)
library(ggpubr)
library(rstatix)


# cargando datos bacterial

tb_01 <- read.table("/media/ins-bio/DATA02/orson_user/avance_tesis_IV/Tb_01.tsv", header = FALSE)

names(tb_01) <- c("ID","INF06_01")

tb_02 <- read.table("/media/ins-bio/DATA02/orson_user/avance_tesis_IV/Tb_02.tsv", header = FALSE)

names(tb_02) <- c("ID","INF06_02")

tb_03 <- read.table("/media/ins-bio/DATA02/orson_user/avance_tesis_IV/Tb_02.tsv", header = FALSE)

names(tb_03) <- c("ID","INF06_03")


treat_tb <- merge(tb_01,tb_02, by = "ID" , all =  TRUE)

inf06 <- merge(treat_tb, tb_03, by = "ID" , all =  TRUE)

inf06[is.na(inf06)] <- 0



tb_11 <- read.table("/media/ins-bio/DATA02/orson_user/avance_tesis_IV/Tb_11.tsv", header = FALSE)

names(tb_11) <- c("ID","INF24_01")

tb_12 <- read.table("/media/ins-bio/DATA02/orson_user/avance_tesis_IV/Tb_12.tsv", header = FALSE)

names(tb_12) <- c("ID","INF24_02")

tb_13 <- read.table("/media/ins-bio/DATA02/orson_user/avance_tesis_IV/Tb_13.tsv", header = FALSE)

names(tb_13) <- c("ID","INF24_03")


treat_tb24 <- merge(tb_11,tb_12, by = "ID" , all =  TRUE)

inf24 <- merge(treat_tb24, tb_13, by = "ID" , all =  TRUE)

inf24[is.na(inf24)] <- 0




tb_21 <- read.table("/media/ins-bio/DATA02/orson_user/avance_tesis_IV/Tb_21.tsv", header = FALSE)

names(tb_21) <- c("ID","INF48_01")

tb_22 <- read.table("/media/ins-bio/DATA02/orson_user/avance_tesis_IV/Tb_22.tsv", header = FALSE)

names(tb_22) <- c("ID","INF48_02")

tb_23 <- read.table("/media/ins-bio/DATA02/orson_user/avance_tesis_IV/Tb_23.tsv", header = FALSE)

names(tb_23) <- c("ID","INF48_03")


treat_tb48 <- merge(tb_21,tb_22, by = "ID" , all =  TRUE)

inf48 <- merge(treat_tb48, tb_23, by = "ID" , all =  TRUE)

inf48[is.na(inf48)] <- 0


infec <- merge(inf06,inf24, by = "ID" , all =  TRUE)
infectados <- merge(infec,inf48, by = "ID" , all =  TRUE)

infectados[is.na(infectados)] <- 0




######## No onfectado 


Tc_01 <- read.table("/media/ins-bio/DATA02/orson_user/avance_tesis_IV/Tc_01.tsv", header = FALSE)

names(Tc_01) <- c("ID","CON06_01")

Tc_02 <- read.table("/media/ins-bio/DATA02/orson_user/avance_tesis_IV/Tc_02.tsv", header = FALSE)

names(Tc_02) <- c("ID","CON06_02")

Tc_03 <- read.table("/media/ins-bio/DATA02/orson_user/avance_tesis_IV/Tc_02.tsv", header = FALSE)

names(Tc_03) <- c("ID","CON06_03")


treat_Tc <- merge(Tc_01,Tc_02, by = "ID" , all =  TRUE)

CON06 <- merge(treat_Tc, Tc_03, by = "ID" , all =  TRUE)

CON06[is.na(CON06)] <- 0



Tc_11 <- read.table("/media/ins-bio/DATA02/orson_user/avance_tesis_IV/Tc_11.tsv", header = FALSE)

names(Tc_11) <- c("ID","CON24_01")

Tc_12 <- read.table("/media/ins-bio/DATA02/orson_user/avance_tesis_IV/Tc_12.tsv", header = FALSE)

names(Tc_12) <- c("ID","CON24_02")

Tc_13 <- read.table("/media/ins-bio/DATA02/orson_user/avance_tesis_IV/Tc_13.tsv", header = FALSE)

names(Tc_13) <- c("ID","CON24_03")


treat_Tc24 <- merge(Tc_11,Tc_12, by = "ID" , all =  TRUE)

CON24 <- merge(treat_Tc24, Tc_13, by = "ID" , all =  TRUE)

CON24[is.na(CON24)] <- 0




Tc_21 <- read.table("/media/ins-bio/DATA02/orson_user/avance_tesis_IV/Tc_21.tsv", header = FALSE)

names(Tc_21) <- c("ID","CON48_01")

Tc_22 <- read.table("/media/ins-bio/DATA02/orson_user/avance_tesis_IV/Tc_22.tsv", header = FALSE)

names(Tc_22) <- c("ID","CON48_02")

Tc_23 <- read.table("/media/ins-bio/DATA02/orson_user/avance_tesis_IV/Tc_23.tsv", header = FALSE)

names(Tc_23) <- c("ID","CON48_03")


treat_Tc48 <- merge(Tc_21,Tc_22, by = "ID" , all =  TRUE)

CON48 <- merge(treat_Tc48, Tc_23, by = "ID" , all =  TRUE)

CON48[is.na(CON48)] <- 0


contr <- merge(CON06,CON24, by = "ID" , all =  TRUE)
controles <- merge(contr,CON48, by = "ID" , all =  TRUE)

controles[is.na(controles)] <- 0




### Terminando de unir la data

rna_seq <- merge(infectados,controles, by = "ID" , all =  TRUE)

rna_seq[is.na(rna_seq)] <- 0



### cargando metadata 


metadata <- read.table("/media/ins-bio/DATA02/orson_user/avance_tesis_IV/metadata.tsv", header = TRUE)
attach(metadata)


#### terminando de set la data


rownames(rna_seq) <- rna_seq$ID

rna_seq$ID <- NULL


##### Analisis de expreison diferencial 



library(glmmSeq)


# Dispersion

disp <- apply(rna_seq, 1, function(x){
  (var(x, na.rm=TRUE)-mean(x, na.rm=TRUE))/(mean(x, na.rm=TRUE)**2)
})

head(disp)

# Size Factors

sizeFactors <- colSums(rna_seq)  
sizeFactors <- sizeFactors / mean(sizeFactors)  # normalise to mean = 1

head(sizeFactors)


# Fitting Models

# gene expression ~ fixed effects + random effects


results <- glmmSeq(~ Timepoint * EULAR_binary  + (1 | PATID),
                   countdata = rna_seq,
                   metadata = metadata,
                   dispersion = disp,
                   progress = TRUE,
                   cores = 25)


stats <- summary(results)

summary(results, gene = "K02072")

predict = data.frame(results@predict)

results <- glmmQvals(results)

library(emmeans)

oldpar <- par(mfrow=c(1, 2))


plotColours <- c("skyblue", "goldenrod1")
modColours <- c("dodgerblue3", "goldenrod3")
shapes <- c(17, 19)

ggmodelPlot(results,
            geneName = "K17584",
            x1var = "Timepoint",
            x2var="EULAR_binary",
            xlab="Time",
            colours = plotColours,
            shapes = shapes,
            lineColours = plotColours, 
            modelColours = modColours,
            modelSize = 10)


ggmodelPlot(results,
            geneName = "K04362",
            x1var = "Timepoint",
            x2var="EULAR_binary",
            xlab="Time",
            colours = plotColours,
            shapes = shapes,
            lineColours = plotColours, 
            modelColours = modColours,
            modelSize = 10)


labels = c("K17584","K04362","K02072")

fcPlot(results, x1var = "Timepoint", x2var = "EULAR_binary", graphics = "plotly",
       pCutoff = 0.05, useAdjusted = TRUE,
       labels = labels,
       colours = c('grey', 'green3', 'gold3','blue'))
