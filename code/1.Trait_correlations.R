# Trait correlations and phylogeny

setwd("~/Dropbox/Data/FIA-Traits/Rproj")

library(readr)
library(tidyverse)
library(phytools)
library(adephylo)
library(ape)
library(phylobase)
library(picante)
library(V.PhyloMaker)

### Figure 2: Phylogenetic tree figure ###
#read in tree
tree <- read.tree("data/phylogeny.tre")
tree$tip.label <- str_replace(tree$tip.label, "_", " ")
tree$tip.label[66] <- "Callitropsis nootkatensis"

traits <- read.csv("data/species_trait_data.csv") %>% as.data.frame() %>%
          select(species, p50 = p50_mpa, rd = rd_m, rd_cat = rd_deep, Resprouter = resprout, bark = bark_thickness_mm,
                  ldmc = ldmc_g_g, sla = sla_mm2_mg, wd = wd_g_cm3, height = hmax_m, seedmass = seed_mass_mg, spcode, counts) %>%
          mutate(Rooting_depth= log(rd), Bark_thickness = log(bark), Xylem_P50= log(-p50), SLA = log(sla), LDMC = ldmc, Height = log(height))

rownames(traits) <- traits$species
p.dist.mat <- cophenetic(tree)         ###to ensure correct order
traits <- traits[row.names(p.dist.mat),] ###to ensure correct order


### Figure 2 ###
trait_data <- traits[, c("Xylem_P50","Rooting_depth","Bark_thickness","Resprouter")]
colnames(trait_data) <- c("Stem P50", "Max rooting depth", "Bark thickness", "Resprouter")

pdf(file="output/Fig2.pdf",width=11, height=12)

phylo.heatmap(tree, trait_data, 
              standardize=TRUE, 
              fsize=c(0.4,0.7,0.7),
              pts=FALSE, 
              split=c(0.7,0.3), 
              colors=viridisLite::viridis(n=7, direction=-1), 
              mar=rep(0.1,4))
text(x=-0.2, y=0.76, "Gymnosperms", adj=0.5, cex=1)
text(x=-0.2, y=0.38, "Angiosperms", adj=0.5, cex=1)
text(x=0.1, y=0.9, "Pinaceae", adj=0.5, cex=0.8)
text(x=0.13, y=0.59, "Cupressaceae", adj=0.5, cex=0.8)
text(x=0.13, y=0.39, "Salicaceae", adj=0.5, cex=0.8)
text(x=0.13, y=0.23, "Fagaceae", adj=0.5, cex=0.8)
text(x=1.05, y=0.8, "No", adj=0.5, cex=1)
text(x=1.05, y=0.4, "Yes", adj=0.5, cex=1)

dev.off()
### end Fig 2 ###


### proportion of species measured for each trait
colSums(!is.na(traits[,c("Xylem_P50","Rooting_depth","Bark_thickness","Resprouter")]))/nrow(traits)





### Figure 3: Bivariate trait plot ###

traits$resprout <- traits$Resprouter
traits$resprout[traits$resprout == 0] <- "No"
traits$resprout[traits$resprout == 1] <- "Yes"
traits$phylo <- c(rep("angio",53),rep("gymno",56)) # phylogenetic group

### 3x2 Figure ####################
pdf(file="output/Fig3.pdf",width=8, height=11)
par(mfrow=c(3,2), mar=c(5,5,2,1), cex.lab=1.4)

# p50 vs resprout 
boxplot(abs(p50) ~ resprout, data=traits, log="y", cex=1.5, pch=19, ylim=c(1,14),
        xlab="Resprouter", ylab="Stem P50 (-MPa)", outline=FALSE)#, col=c("gold","blue"))
for(i in seq_along(levels(factor(traits$resprout)))) {
  subset_data <- traits[traits$resprout == levels(factor(traits$resprout))[i], ]
  y_vals <- abs(subset_data$p50)
  y_vals <- y_vals[!is.na(y_vals)]
  stripchart(y_vals, vertical=TRUE, method="jitter", 
             at=i, add=TRUE, pch=19, cex=1.5, col=rgb(0,0,0,0.5))
}
title(main="A", adj=0)
summary(m3 <- lm(log(abs(p50)) ~ factor(resprout), data=traits))
t.test(log(abs(p50)) ~ factor(resprout), data=traits)
exp(1.513178); exp(1.091794)
title(main = "N = 69 species, t = 3.8, P = 0.0003", adj = 1)
text(1.5, 14,"***", cex=2)
polygon(x=c(0.5,1.7,1.7,0.5), y=c(1,1,1,3),
        col =  adjustcolor("grey70", alpha.f = 0.30), border = NA)
text(1,1.3, "Uncommon
phenotype", cex=1.2)
polygon(x=c(1.8,2.5,2.5), y=c(14,14,5),
        col =  adjustcolor("grey70", alpha.f = 0.30), border = NA)
text(2.25,10.8, "Uncommon
phenotype", cex=1.2)

# p50 vs bark_thickness 
plot(abs(p50) ~ bark, data=traits, log="xy", pch=19, cex=1.5, 
     xlab="Bark thickness (mm at 20 cm DBH)",
     ylab="Stem P50 (- MPa)")
title(main="B", adj=0)
cor.test(log(abs(traits$p50)), log(traits$bark))
title(main = "N = 61 species, r = 0.20, P = 0.11", adj = 1)
summary(m4 <- lm(log(abs(p50)) ~ log(bark), data=traits))
bark_new <- data.frame(bark = seq(2.8,17,by=0.1))
preds2 <- predict(m4, newdata=bark_new, se=TRUE)
lines(bark_new$bark, exp(preds2$fit), lwd=2, lty=2)

# rd vs resprout
boxplot(rd ~ resprout, data=traits, log="y", cex=1.5, pch=19,ylim=c(1,55),
        xlab="Resprouter",
        ylab="Maximum rooting depth (m)", outline=FALSE)
for(i in seq_along(levels(factor(traits$resprout)))) {
  subset_data <- traits[traits$resprout == levels(factor(traits$resprout))[i], ]
  y_vals <- abs(subset_data$rd)
  y_vals <- y_vals[!is.na(y_vals)]
  stripchart(y_vals, vertical=TRUE, method="jitter", 
             at=i, add=TRUE, pch=19, cex=1.5, col=rgb(0,0,0,0.5))
}
title(main="C", adj=0)
summary(m1 <- lm(log(rd) ~ resprout, data=traits))
t.test(log(traits$rd) ~ traits$resprout)
title(main = "N = 44 species, t = -0.2, P = 0.84", adj = 1)

# rd vs bark_thickness
plot(rd ~ bark, log="xy", data=traits, xlim=c(4,17), pch=19, cex=1.5,
     xlab="Bark thickness (mm at 20 cm DBH)",
     ylab="Maximum rooting depth (m)")
title(main="D", adj=0)
title(main = "N = 37 species, r = 0.36, P = 0.025", adj = 1)
log_root <- log(traits$rd)
log_bark <- log(traits$bark)
summary(m2 <- lm(log_root ~ log_bark))
bark_new <- data.frame(log_bark = seq(log(4),log(17),length.out=16))
preds <- predict(m2, newdata=bark_new, se=TRUE)
lines(exp(bark_new$log_bark), exp(preds$fit), lwd=2)
lines(exp(bark_new$log_bark), exp(preds$fit) + exp(preds$se.fit), lty=1)
lines(exp(bark_new$log_bark), exp(preds$fit) - exp(preds$se.fit), lty=1)
polygon(x = c(exp(bark_new$log_bark), rev(exp(bark_new$log_bark))),
        y = c(exp(preds$fit) - exp(preds$se.fit), rev(exp(preds$fit) + exp(preds$se.fit))),
        col =  adjustcolor("grey80", alpha.f = 0.3), border = NA)
summary(m2 <- lm(log(rd) ~ log(bark), data=traits))
cor.test(log(traits$rd) , log(traits$bark))
polygon(x=c(4,8, 8, 4), y=c(10,58,58,58),
        col =  adjustcolor("grey70", alpha.f = 0.30), border = NA)
polygon(x=c(10,17,17,10), y=c(0.6,0.6,2,0.6),
        col =  adjustcolor("grey70", alpha.f = 0.30), border = NA)
text(5,35, "Uncommon
phenotype", cex=1.2)
text(14,0.9, "Uncommon
phenotype", cex=1.2)

# bark vs resprout
boxplot(bark ~ resprout, data=traits, log="y", cex=1.5, pch=19,
        xlab="Resprouter",
        ylab="Bark thickness (mm at 20 cm DBH)", outline=FALSE)
for(i in seq_along(levels(factor(traits$resprout)))) {
  subset_data <- traits[traits$resprout == levels(factor(traits$resprout))[i], ]
  y_vals <- abs(subset_data$bark)
  y_vals <- y_vals[!is.na(y_vals)]
  stripchart(y_vals, vertical=TRUE, method="jitter", 
             at=i, add=TRUE, pch=19, cex=1.5, col=rgb(0,0,0,0.5))
}
title(main="E", adj=0)
summary(m6 <- lm(log(bark) ~ resprout, data=traits))
t.test(log(traits$bark) ~ traits$resprout)
title(main = "N = 81 species, t = 1.7, P = 0.08", adj = 1)

# p50 vs rooting depth
plot(rd ~ abs(p50), log="xy", data=traits, xlim=c(1.3,10), pch=19, cex=1.5,
     xlab="Stem P50 (- MPa)",
     ylab="Maximum rooting depth (m)")
title(main="F", adj=0)
summary(m5 <- lm(log(rd) ~ log(abs(p50)), data=traits))
p50_new <- data.frame(p50 = seq(1.3,10,by=0.1))
preds4 <- predict(m5, newdata=p50_new, se=TRUE)
lines(p50_new$p50, exp(preds4$fit), lwd=2, lty=2)
cor.test(log(traits$rd), abs(traits$p50))
title(main = "N = 37 species, r = 0.08, P = 0.62", adj = 1)

dev.off()
### End fig 3 ###




### Table S1 ###
tableS1 <- traits %>%
  select(species, P50=p50, Rootingdepth = rd, Barkthickness = bark, Resprouter) %>%
  arrange(species) %>%
  mutate(P50=round(P50,2), Rootingdepth=round(Rootingdepth,1),Barkthickness =round(Barkthickness,1))
tableS1$Resprouter <- ifelse(tableS1$Resprouter == 1, "Yes", "No")
tableS1$P50 <- ifelse(is.na(tableS1$P50), "---", tableS1$P50)
tableS1$Rootingdepth <- ifelse(is.na(tableS1$Rootingdepth), "---", tableS1$Rootingdepth)
tableS1$Barkthickness <- ifelse(is.na(tableS1$Barkthickness), "---", tableS1$Barkthickness)
write.csv(tableS1, file="output/TableS1.csv")
### End Table S1 ###
