### Functional ecology of western US trees ###

library(tidyverse)
library(mgcv)
library(glmmTMB)
library(dplyr)
library(tidyr)
library(visreg)
library(ggplot2)
library(maps)
library(patchwork)

setwd("~/Dropbox/Data/FIA-Traits/Rproj")

# Read in data ------------------------------------------------------------

# presence-absence data
d <- readRDS("data/pres_abs_tree.rds") %>% as_tibble()

# subsample 25% of absences
d_pres <- d %>% filter(pres==1)
d_abs <- d %>% filter(pres==0)
nrow(d_abs)/nrow(d_pres)
set.seed(02201924)
d_abs_sample <- slice_sample(d_abs, prop = 0.25)
nrow(d_abs_sample)/nrow(d_pres)
d <- rbind(d_pres, d_abs_sample) %>% arrange(province, plt_cn, species)

prov <- read.csv("data/provinces.csv", header=TRUE)
d <- d %>% left_join(prov, by="province")
# provinces https://www.fs.usda.gov/land/ecosysmgmt/colorimagemap/ecoreg1_provinces.html
# unique(d$province) > "M242" "341"  "M341" "342"  "M331" "313"  "322"  "M262" "M261" "261"  "263"  "M332" "M333" "M313" "242"  "262"

# trait data
trait_data <- read.csv("data/species_trait_data.csv") %>%
  select(species, p50 = p50_mpa, rd = rd_m, resprout = resprout, bark = bark_thickness_mm)

# environment
env <- readRDS("data/plot_environment_data.rds")

# precipitation seasonality
clim_monthly <- readRDS("data/plot_climate_data.rds")
clim <- clim_monthly %>% select(plt_cn, pptcv)

# merge data
d <- d %>%
  inner_join(trait_data) %>%
  inner_join(env) %>%
  inner_join(clim) %>%
  mutate(slog_rd = scale(log(rd)),
         slog_bark = scale(log(bark)),
         slog_p50 = scale(log(-p50)),
         resprout = factor(resprout),
         slog_arid = scale(log(arid)),
         slog_pptcv = scale(log(pptcv)),
         slog_wtd = scale(log(-wtd)),
         slog_fri = scale(log(fri)),
         species = factor(species), 
         province = factor(province),
         region = factor(region),
         plt_cn = factor(plt_cn))

# load unscale function
unscale <- function(z, center = attr(z, "scaled:center"), scale = attr(z, "scaled:scale")) {
  if(!is.null(scale))  z <- sweep(z, 2, scale, `*`)
  if(!is.null(center)) z <- sweep(z, 2, center, `+`)
  structure(z,
            "scaled:center"   = NULL,
            "scaled:scale"    = NULL,
            "unscaled:center" = center,
            "unscaled:scale"  = scale
  )
} 

# Calculate quantiles of env variables
env_quantiles <- d %>%
  summarise(
    arid_0.025 = exp(quantile(log(arid), 0.025, na.rm = TRUE)),
    arid_0.5 = exp(quantile(log(arid), 0.5, na.rm = TRUE)),
    arid_0.975 = exp(quantile(log(arid), 0.975, na.rm = TRUE)),
    pptcv_0.025 = exp(quantile(log(pptcv), 0.025, na.rm = TRUE)),
    pptcv_0.5 = exp(quantile(log(pptcv), 0.5, na.rm = TRUE)),
    pptcv_0.975 = exp(quantile(log(pptcv), 0.975, na.rm = TRUE)),
    wtd_0.025 = exp(quantile(log(abs(wtd)), 0.025, na.rm = TRUE)),
    wtd_0.5 = exp(quantile(log(abs(wtd)), 0.5, na.rm = TRUE)),
    wtd_0.975 = exp(quantile(log(abs(wtd)), 0.975, na.rm = TRUE)),
    fri_0.025 = exp(quantile(log(fri), 0.025, na.rm = TRUE)),
    fri_0.5 = exp(quantile(log(fri), 0.5, na.rm = TRUE)),
    fri_0.975 = exp(quantile(log(fri), 0.975, na.rm = TRUE))
  ) %>%
  pivot_longer(everything(), names_to = "var_quantile", values_to = "value") %>%
  separate(var_quantile, into = c("variable", "quantile"), sep = "_") %>%
  pivot_wider(names_from = quantile, values_from = value) %>%
  mutate(variable = factor(variable, 
                           levels = c("arid", "pptcv", "wtd", "fri"),
                           labels = c("Aridity (P:PET)", "Precipitation Seasonality (CV)",
                                      "Water table depth (m)", "Fire return interval (yrs)"))) %>%
  arrange(variable) %>%
  mutate(across(c(`0.025`, `0.5`, `0.975`), ~round(., 2)))

# Save as CSV
#write.csv(env_quantiles, "ms/env_quantiles.csv", row.names = FALSE)


# Calculate quantiles
quantiles <- d %>%
  summarise(
    arid_0.025 = exp(quantile(log(arid), 0.025, na.rm = TRUE)),
    arid_0.5 = exp(quantile(log(arid), 0.5, na.rm = TRUE)),
    arid_0.975 = exp(quantile(log(arid), 0.975, na.rm = TRUE)),
    pptcv_0.025 = exp(quantile(log(pptcv), 0.025, na.rm = TRUE)),
    pptcv_0.5 = exp(quantile(log(pptcv), 0.5, na.rm = TRUE)),
    pptcv_0.975 = exp(quantile(log(pptcv), 0.975, na.rm = TRUE)),
    wtd_0.025 = exp(quantile(log(abs(wtd)), 0.025, na.rm = TRUE)),
    wtd_0.5 = exp(quantile(log(abs(wtd)), 0.5, na.rm = TRUE)),
    wtd_0.975 = exp(quantile(log(abs(wtd)), 0.975, na.rm = TRUE)),
    fri_0.025 = exp(quantile(log(fri), 0.025, na.rm = TRUE)),
    fri_0.5 = exp(quantile(log(fri), 0.5, na.rm = TRUE)),
    fri_0.975 = exp(quantile(log(fri), 0.975, na.rm = TRUE))
  )
# Get US state boundaries
states <- map_data("state")
western_states <- c("washington", "oregon", "california", "idaho", "nevada", 
                    "utah", "arizona", "montana", "wyoming", "colorado", 
                    "new mexico")
states_west <- states %>% filter(region %in% western_states)

lon_range <- range(d$lon, na.rm = TRUE)
lat_range <- range(d$lat, na.rm = TRUE)

# Create table text for each variable
fri_table <- paste0("2.5%: ", round(quantiles$fri_0.025, 1), "\n",
                    "50%: ", round(quantiles$fri_0.5, 1), "\n",
                    "97.5%: ", round(quantiles$fri_0.975, 1))

arid_table <- paste0("2.5%: ", round(quantiles$arid_0.025, 2), "\n",
                     "50%: ", round(quantiles$arid_0.5, 2), "\n",
                     "97.5%: ", round(quantiles$arid_0.975, 2))

pptcv_table <- paste0("2.5%: ", round(quantiles$pptcv_0.025, 2), "\n",
                      "50%: ", round(quantiles$pptcv_0.5, 2), "\n",
                      "97.5%: ", round(quantiles$pptcv_0.975, 2))

wtd_table <- paste0("2.5%: ", round(quantiles$wtd_0.025, 1), "\n",
                    "50%: ", round(quantiles$wtd_0.5, 1), "\n",
                    "97.5%: ", round(quantiles$wtd_0.975, 1))

p1 <- d %>%
  distinct(lon, lat, fri) %>%
  ggplot() +
  geom_polygon(data = states_west, aes(x = long, y = lat, group = group),
               fill = NA, color = "gray30", linewidth = 0.3) +
  geom_point(aes(x = lon, y = lat, color = fri), size = 0.5, alpha = 0.5) +
  scale_color_viridis_c(name = "Fire return\ninterval (yrs)", trans = "log10") +
  annotate("label", x = -125, y = 35, label = fri_table, 
           size = 4, hjust = 0, vjust = 1, fill = "white", alpha = 0.8) +
  coord_quickmap(xlim = lon_range, ylim = lat_range) +
  theme_bw() +
  theme(aspect.ratio = diff(lat_range)/diff(lon_range),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7))

p2 <- d %>%
  distinct(lon, lat, arid) %>%
  ggplot() +
  geom_polygon(data = states_west, aes(x = long, y = lat, group = group),
               fill = NA, color = "gray30", linewidth = 0.3) +
  geom_point(aes(x = lon, y = lat, color = arid), size = 0.5, alpha = 0.5) +
  scale_color_viridis_c(name = "P:PET", trans = "log10") +
  annotate("label", x = -125, y = 35, label = arid_table, 
           size = 4, hjust = 0, vjust = 1, fill = "white", alpha = 0.8) +
  coord_quickmap(xlim = lon_range, ylim = lat_range) +
  theme_bw() +
  theme(aspect.ratio = diff(lat_range)/diff(lon_range),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7))

p3 <- d %>%
  distinct(lon, lat, pptcv) %>%
  ggplot() +
  geom_polygon(data = states_west, aes(x = long, y = lat, group = group),
               fill = NA, color = "gray30", linewidth = 0.3) +
  geom_point(aes(x = lon, y = lat, color = pptcv), size = 0.5, alpha = 0.5) +
  scale_color_viridis_c(name = "Precipitation\nseasonality (CV)", trans = "log10") +
  annotate("label", x = -125, y = 35, label = pptcv_table, 
           size = 4, hjust = 0, vjust = 1, fill = "white", alpha = 0.8) +
  coord_quickmap(xlim = lon_range, ylim = lat_range) +
  theme_bw() +
  theme(aspect.ratio = diff(lat_range)/diff(lon_range),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7))

p4 <- d %>%
  distinct(lon, lat, wtd) %>%
  mutate(wtd_pos = abs(wtd)) %>%
  ggplot() +
  geom_polygon(data = states_west, aes(x = long, y = lat, group = group),
               fill = NA, color = "gray30", linewidth = 0.3) +
  geom_point(aes(x = lon, y = lat, color = wtd_pos), size = 0.5, alpha = 0.5) +
  scale_color_viridis_c(name = "Water table\ndepth (m)", trans = "log10") +
  annotate("label", x = -125, y = 35, label = wtd_table, 
           size = 4, hjust = 0, vjust = 1, fill = "white", alpha = 0.8) +
  coord_quickmap(xlim = lon_range, ylim = lat_range) +
  theme_bw() +
  theme(aspect.ratio = diff(lat_range)/diff(lon_range),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7))

# Combine into 2x2 panel
combined_plot <- (p1 | p2) / (p3 | p4)

# Save
ggsave("output/FigS1.pdf", combined_plot, 
       width = 12, height = 10)




#### Models ##################################

# Define variable combinations
traits <- c("slog_p50", "slog_rd", "slog_bark", "factor(resprout)")
envs <- c("slog_arid", "slog_wtd", "slog_pptcv", "slog_fri")

# Create all combinations
combinations <- expand.grid(trait = traits,env = envs,stringsAsFactors = FALSE)

# Fit all models
models <- lapply(seq_len(nrow(combinations)), function(i) {
  trait <- combinations$trait[i]
  env <- combinations$env[i]
  # Build formula
  formula_str <- sprintf("pres ~ %s * %s + s(species, bs = 're') + s(province, bs = 're')",env, trait)
  # Fit model
  model <- bam(
    formula = as.formula(formula_str),family = binomial(),data = d,discrete = TRUE)
  # Return model with metadata
  list(
    model = model,trait = trait,env = env,formula = formula_str)
})

# Name list elements for easy access
names(models) <- sprintf("%s_x_%s",combinations$env,gsub("factor\\(|\\)|slog_", "", combinations$trait)
)
# Optional: Extract just the model objects
model_objects <- lapply(models, `[[`, "model")
names(model_objects) <- names(models)
# Access specific model
models$slog_arid_x_p50$model
models$slog_wtd_x_resprout$model
# Get summary
summary(models$slog_fri_x_bark$model)
# View all model names
names(models)
# Function to extract all model information in one table
extract_model_summary <- function(model_list) {
  summary_df <- lapply(seq_along(model_list), function(i) {
    m <- model_list[[i]]
    model_name <- names(model_list)[i]
    summ <- summary(m$model)
    # Extract model data
    model_data <- m$model$model
    # Count unique species and provinces
    n_species <- length(unique(model_data$species))
    n_provinces <- length(unique(model_data$province))
    # Get coefficients for interaction term only
    coef_summary <- summ$p.table
    interaction_row <- grep(":", rownames(coef_summary))
    if (length(interaction_row) > 0) {
      interaction_est <- coef_summary[interaction_row, "Estimate"]
      interaction_se <- coef_summary[interaction_row, "Std. Error"]
      interaction_z <- coef_summary[interaction_row, "z value"]
      interaction_p <- coef_summary[interaction_row, "Pr(>|z|)"]
      interaction_sig <- case_when(
        interaction_p < 0.001 ~ "***",
        interaction_p < 0.01 ~ "**",
        interaction_p < 0.05 ~ "*",
        interaction_p < 0.1 ~ ".",
        TRUE ~ ""
      )
      interaction_term <- rownames(coef_summary)[interaction_row]
    } else {
      interaction_est <- NA
      interaction_se <- NA
      interaction_z <- NA
      interaction_p <- NA
      interaction_sig <- ""
      interaction_term <- NA
    }
    
    data.frame(
      model = model_name,
      n_obs = nrow(model_data),
      n_species = n_species,
      n_provinces = n_provinces,
      interaction_term = interaction_term,
      interaction_est = round(interaction_est, 4),
      interaction_se = round(interaction_se, 4),
      interaction_z = round(interaction_z, 3),
      interaction_p = round(interaction_p, 4),
      interaction_sig = interaction_sig,
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()
  
  return(summary_df)
}

# Extract comprehensive summary
model_summary <- extract_model_summary(models)

# Save to CSV
write.csv(model_summary, "output/model_summary.csv", row.names = FALSE)

####################################################################


### Rename all 16 models m1 through m16
m1 <- models$slog_arid_x_p50$model
m2 <- models$slog_arid_x_rd$model
m3 <- models$slog_arid_x_bark$model
m4 <- models$slog_arid_x_resprout$model
m5 <- models$slog_wtd_x_p50$model
m6 <- models$slog_wtd_x_rd$model
m7 <- models$slog_wtd_x_bark$model
m8 <- models$slog_wtd_x_resprout$model
m9 <- models$slog_pptcv_x_p50$model
m10 <- models$slog_pptcv_x_rd$model
m11 <- models$slog_pptcv_x_bark$model
m12 <- models$slog_pptcv_x_resprout$model
m13 <- models$slog_fri_x_p50$model
m14 <- models$slog_fri_x_rd$model
m15 <- models$slog_fri_x_bark$model
m16 <- models$slog_fri_x_resprout$model

# m1: arid x p50
visreg(m1, xvar = "slog_p50", by = "slog_arid", scale = "response", overlay = TRUE, breaks = c(-2, 2))
# m2: arid x rd
visreg(m2, xvar = "slog_rd", by = "slog_arid", scale = "response", overlay = TRUE, breaks = c(-2, 2))
# m3: arid x bark
visreg(m3, xvar = "slog_bark", by = "slog_arid", scale = "response", overlay = TRUE, breaks = c(-2, 2))
# m4: arid x resprout
visreg(m4, xvar = "resprout", by = "slog_arid", scale = "response", overlay = TRUE, breaks = c(-2, 2))
# m5: wtd x p50
visreg(m5, xvar = "slog_p50", by = "slog_wtd", scale = "response", overlay = TRUE, breaks = c(-2, 2))
# m6: wtd x rd
visreg(m6, xvar = "slog_rd", by = "slog_wtd", scale = "response", overlay = TRUE, breaks = c(-2, 2))
# m7: wtd x bark
visreg(m7, xvar = "slog_bark", by = "slog_wtd", scale = "response", overlay = TRUE, breaks = c(-2, 2))
# m8: wtd x resprout
visreg(m8, xvar = "resprout", by = "slog_wtd", scale = "response", overlay = TRUE, breaks = c(-2, 2))
# m9: pptcv x p50
visreg(m9, xvar = "slog_p50", by = "slog_pptcv", scale = "response", overlay = TRUE, breaks = c(-2, 2))
# m10: pptcv x rd
visreg(m10, xvar = "slog_rd", by = "slog_pptcv", scale = "response", overlay = TRUE, breaks = c(-2, 2))
# m11: pptcv x bark
visreg(m11, xvar = "slog_bark", by = "slog_pptcv", scale = "response", overlay = TRUE, breaks = c(-2, 2))
# m12: pptcv x resprout
visreg(m12, xvar = "resprout", by = "slog_pptcv", scale = "response", overlay = TRUE, breaks = c(-2, 2))
# m13: fri x p50
visreg(m13, xvar = "slog_p50", by = "slog_fri", scale = "response", overlay = TRUE, breaks = c(-2, 2))
# m14: fri x rd
visreg(m14, xvar = "slog_rd", by = "slog_fri", scale = "response", overlay = TRUE, breaks = c(-2, 2))
# m15: fri x bark
visreg(m15, xvar = "slog_bark", by = "slog_fri", scale = "response", overlay = TRUE, breaks = c(-2, 2))
# m16: fri x resprout
visreg(m16, xvar = "resprout", by = "slog_fri", scale = "response", overlay = TRUE, breaks = c(-2, 2))



###### Fig 4 ########################
pdf(file="output/Fig4.pdf",width=8, height=11)

par(mfcol=c(3,2), mar=c(6,5,2,1))

m1$model$slog_p50_unscaled <- unscale(m1$model$slog_p50)
m1$model$slog_arid_unscaled <- unscale(m1$model$slog_arid)
l1 <- lm(slog_arid_unscaled ~ as.numeric(slog_arid), data=m1$model)
l2 <- lm(slog_p50_unscaled ~ as.numeric(slog_p50), data=m1$model)
newdat <- data.frame(slog_arid=seq(-5,5,1), slog_p50=seq(-5,5,1))
xscale <- data.frame(sd=seq(-5,5,1), arid=round(exp(predict(l1,newdat)),0), p50=round(-exp(predict(l2,newdat)),1))

### Plot TxE interactions with confidence intervals
newdat <- data.frame(
  slog_p50 = rep(seq(quantile(m1$model$slog_p50,0.025),quantile(m1$model$slog_p50,0.975),length.out = 100), 2),
  slog_arid = c(rep(quantile(m1$model$slog_arid, 0.025), 100), rep(quantile(m1$model$slog_arid, 0.975), 100)),
  lab = factor(rep(paste0(c("low.", "high."), "arid"), each = 100))) %>%
  mutate(`slog_arid:slog_p50` = slog_arid * slog_p50)
mm <- model.matrix( ~ slog_arid + slog_p50 + slog_arid:slog_p50, newdat)
newdat$eta <- mm %*% coef(m1)[1:4]
newdat$ProbOccurrence <- exp(newdat$eta) / (1 + exp(newdat$eta))
newdat$VarPred <- diag(mm %*% vcov(m1)[1:4,1:4] %*% t(mm))
newdat$sePred <- sqrt(newdat$VarPred)
newdat$seLow <-
  exp(newdat$eta - 1.96 * newdat$sePred) / (1 + exp(newdat$eta - 1.96 * newdat$sePred))
newdat$seHigh <-
  exp(newdat$eta + 1.96 * newdat$sePred) / (1 + exp(newdat$eta + 1.96 * newdat$sePred))


plot(newdat$slog_p50[1:100], newdat$ProbOccurrence[1:100], type="l", col="gold",
     ylim=c(0,0.9), xlab="", ylab="Probability of Occurrence", cex.lab=1.2)
title(sub=list(expression(paste("Stem P50 (MPa)")), cex=1.2),line=2.5)
title(sub="\n(scaled to unit variance below, log scale above)",line=3.5)
title(main="A", adj=0)
polygon(x=c(seq(min(newdat$slog_p50),max(newdat$slog_p50),length.out=100),
            rev(seq(min(newdat$slog_p50),max(newdat$slog_p50),length.out=100))),
        y=c(newdat$seLow[1:100], rev(newdat$seHigh[1:100])),
        col=adjustcolor("gold",alpha.f=0.5), border=NA)
lines(newdat$slog_p50[1:100], newdat$ProbOccurrence[1:100], col="gold")
polygon(x=c(seq(min(newdat$slog_p50),max(newdat$slog_p50),length.out=100),
            rev(seq(min(newdat$slog_p50),max(newdat$slog_p50),length.out=100))),
        y=c(newdat$seLow[101:200], rev(newdat$seHigh[101:200])),
        col=adjustcolor("darkgreen",alpha.f=0.5), border=NA)
lines(newdat$slog_p50[101:200], newdat$ProbOccurrence[101:200],col="darkgreen")
axis(3, at=xscale$sd, labels=xscale$p50, line=0)
title("",line=2.5)
legend("topright", legend=c("Arid (2.5th %)","Humid (97.5th %)"), lty=1,
       col=c("gold","darkgreen"), lwd=3, bty="n", cex=1)
text(-1.5,0.25,"Humid",col="darkgreen", cex=1.4)
text(0.9,0.4,"Arid",col="gold", cex=1.4)
text(0,0.7, "TxE: P < 0.0001")

# rd x Aridity
m2$model$slog_rd_unscaled <- unscale(m2$model$slog_rd)
m2$model$slog_arid_unscaled <- unscale(m2$model$slog_arid)
l1 <- lm(slog_arid_unscaled ~ as.numeric(slog_arid), data=m2$model)
l2 <- lm(slog_rd_unscaled ~ as.numeric(slog_rd), data=m2$model)
newdat <- data.frame(slog_arid=seq(-5,5,1), slog_rd=seq(-5,5,1))
xscale <- data.frame(sd=seq(-5,5,1), arid=round(exp(predict(l1,newdat)),2), rd=round(exp(predict(l2,newdat)),1))

### Plot TxE interactions with confidence intervals
newdat <- data.frame(
  slog_rd = rep(seq(quantile(m2$model$slog_rd,0.025),quantile(m2$model$slog_rd,0.975),length.out = 100), 2),
  slog_arid = c(rep(quantile(m2$model$slog_arid, 0.025), 100), rep(quantile(m2$model$slog_arid, 0.975), 100)),
  lab = factor(rep(paste0(c("low.", "high."), "arid"), each = 100))) %>%
  mutate(`slog_arid:slog_rd` = slog_arid * slog_rd)
mm <- model.matrix( ~ slog_arid + slog_rd + slog_arid:slog_rd, newdat)
newdat$eta <- mm %*% coef(m2)[1:4]
newdat$ProbOccurrence <- exp(newdat$eta) / (1 + exp(newdat$eta))
newdat$VarPred <- diag(mm %*% vcov(m2)[1:4,1:4] %*% t(mm))
newdat$sePred <- sqrt(newdat$VarPred)
newdat$seLow <-
  exp(newdat$eta - 1.96 * newdat$sePred) / (1 + exp(newdat$eta - 1.96 * newdat$sePred))
newdat$seHigh <-
  exp(newdat$eta + 1.96 * newdat$sePred) / (1 + exp(newdat$eta + 1.96 * newdat$sePred))

plot(newdat$slog_rd[1:100], newdat$ProbOccurrence[1:100], type="l", col="gold",
     ylim=c(0,0.8), xlab="", ylab="Probability of Occurrence", cex.lab=1.2)
title(sub=list(expression(paste("Maximum rooting depth (m)")), cex=1.2),line=2.5)
title(sub="\n(scaled to unit variance below, log scale above)",line=3.5)
title(main="B", adj=0)
polygon(x=c(seq(min(newdat$slog_rd),max(newdat$slog_rd),length.out=100),
            rev(seq(min(newdat$slog_rd),max(newdat$slog_rd),length.out=100))),
        y=c(newdat$seLow[1:100], rev(newdat$seHigh[1:100])),
        col=adjustcolor("gold",alpha.f=0.5), border=NA)
lines(newdat$slog_rd[1:100], newdat$ProbOccurrence[1:100], col="gold")
polygon(x=c(seq(min(newdat$slog_rd),max(newdat$slog_rd),length.out=100),
            rev(seq(min(newdat$slog_rd),max(newdat$slog_rd),length.out=100))),
        y=c(newdat$seLow[101:200], rev(newdat$seHigh[101:200])),
        col=adjustcolor("darkgreen",alpha.f=0.5), border=NA)
lines(newdat$slog_rd[101:200], newdat$ProbOccurrence[101:200],col="darkgreen")
axis(3, at=xscale$sd, labels=xscale$rd, line=0)
title("",line=2.5)
legend("topleft", legend=c("Arid (2.5th %)","Humid (97.5th %)"), lty=1,
       col=c("gold","darkgreen"), lwd=3, bty="n", cex=1)
text(-0.5,0.4,"Humid",col="darkgreen", cex=1.4)
text(0.9,0.4,"Arid",col="gold", cex=1.4)
text(0,0.6, "TxE: P < 0.0001")

## rd and pptcv
m10$model$slog_rd_unscaled <- unscale(m10$model$slog_rd)
m10$model$slog_pptcv_unscaled <- unscale(m10$model$slog_pptcv)
l1 <- lm(slog_pptcv_unscaled ~ as.numeric(slog_pptcv), data=m10$model)
l2 <- lm(slog_rd_unscaled ~ as.numeric(slog_rd), data=m10$model)
newdat <- data.frame(slog_pptcv=seq(-5,5,1), slog_rd=seq(-5,5,1))
xscale <- data.frame(sd=seq(-5,5,1), pptcv=round(exp(predict(l1,newdat)),2), rd=round(exp(predict(l2,newdat)),1))

### Plot TxE interactions with confidence intervals
newdat <- data.frame(
  slog_rd = rep(seq(quantile(m10$model$slog_rd,0.025),quantile(m10$model$slog_rd,0.975),length.out = 100), 2),
  slog_pptcv = c(rep(quantile(m10$model$slog_pptcv, 0.025), 100), rep(quantile(m10$model$slog_pptcv, 0.975), 100)),
  lab = factor(rep(paste0(c("low.", "high."), "pptcv"), each = 100))) %>%
  mutate(`slog_pptcv:slog_rd` = slog_pptcv * slog_rd)
mm <- model.matrix( ~ slog_pptcv + slog_rd + slog_pptcv:slog_rd, newdat)
newdat$eta <- mm %*% coef(m10)[1:4]
newdat$ProbOccurrence <- exp(newdat$eta) / (1 + exp(newdat$eta))
newdat$VarPred <- diag(mm %*% vcov(m10)[1:4,1:4] %*% t(mm))
newdat$sePred <- sqrt(newdat$VarPred)
newdat$seLow <-
  exp(newdat$eta - 1.96 * newdat$sePred) / (1 + exp(newdat$eta - 1.96 * newdat$sePred))
newdat$seHigh <-
  exp(newdat$eta + 1.96 * newdat$sePred) / (1 + exp(newdat$eta + 1.96 * newdat$sePred))

plot(newdat$slog_rd[1:100], newdat$ProbOccurrence[1:100], type="l", col="darkgreen",
     ylim=c(0,0.9), xlab="", ylab="Probability of Occurrence", cex.lab=1.2)
title(sub=list(expression(paste("Maximum rooting depth (m)")), cex=1.2),line=2.5)
title(sub="\n(scaled to unit variance below, log scale above)",line=3.5)
title(main="C", adj=0)
polygon(x=c(seq(min(newdat$slog_rd),max(newdat$slog_rd),length.out=100),
            rev(seq(min(newdat$slog_rd),max(newdat$slog_rd),length.out=100))),
        y=c(newdat$seLow[1:100], rev(newdat$seHigh[1:100])),
        col=adjustcolor("darkgreen",alpha.f=0.5), border=NA)
lines(newdat$slog_rd[1:100], newdat$ProbOccurrence[1:100], col="darkgreen")
polygon(x=c(seq(min(newdat$slog_rd),max(newdat$slog_rd),length.out=100),
            rev(seq(min(newdat$slog_rd),max(newdat$slog_rd),length.out=100))),
        y=c(newdat$seLow[101:200], rev(newdat$seHigh[101:200])),
        col=adjustcolor("gold",alpha.f=0.5), border=NA)
lines(newdat$slog_rd[101:200], newdat$ProbOccurrence[101:200],col="gold")
axis(3, at=xscale$sd, labels=xscale$rd, line=0)
title("",line=2.5)
legend("topleft", legend=c("Low Precip Seasonality (2.5th %)","High Precip Seasonality (97.5th %)"), lty=1,
       col=c("darkgreen","gold"), lwd=3, bty="n", cex=1)
text(-1,0.3,"Low",col="darkgreen", cex=1.4)
text(1.5,0.6,"High",col="gold", cex=1.4)
text(0,0.6, "TxE: P < 0.0001")

## bark and pptcv
m11$model$slog_bark_unscaled <- unscale(m11$model$slog_bark)
m11$model$slog_pptcv_unscaled <- unscale(m11$model$slog_pptcv)
l1 <- lm(slog_pptcv_unscaled ~ as.numeric(slog_pptcv), data=m11$model)
l2 <- lm(slog_bark_unscaled ~ as.numeric(slog_bark), data=m11$model)
newdat <- data.frame(slog_pptcv=seq(-5,5,1), slog_bark=seq(-5,5,1))
xscale <- data.frame(sd=seq(-5,5,1), pptcv=round(exp(predict(l1,newdat)),2), bark=round(exp(predict(l2,newdat)),1))

### Plot TxE interactions with confidence intervals
newdat <- data.frame(
  slog_bark = rep(seq(quantile(m11$model$slog_bark,0.025),quantile(m11$model$slog_bark,0.975),length.out = 100), 2),
  slog_pptcv = c(rep(quantile(m11$model$slog_pptcv, 0.025), 100), rep(quantile(m11$model$slog_pptcv, 0.975), 100)),
  lab = factor(rep(paste0(c("low.", "high."), "pptcv"), each = 100))) %>%
  mutate(`slog_pptcv:slog_bark` = slog_pptcv * slog_bark)
mm <- model.matrix( ~ slog_pptcv + slog_bark + slog_pptcv:slog_bark, newdat)
newdat$eta <- mm %*% coef(m11)[1:4]
newdat$ProbOccurrence <- exp(newdat$eta) / (1 + exp(newdat$eta))
newdat$VarPred <- diag(mm %*% vcov(m11)[1:4,1:4] %*% t(mm))
newdat$sePred <- sqrt(newdat$VarPred)
newdat$seLow <-
  exp(newdat$eta - 1.96 * newdat$sePred) / (1 + exp(newdat$eta - 1.96 * newdat$sePred))
newdat$seHigh <-
  exp(newdat$eta + 1.96 * newdat$sePred) / (1 + exp(newdat$eta + 1.96 * newdat$sePred))

plot(newdat$slog_bark[1:100], newdat$ProbOccurrence[1:100], type="l", col="darkgreen",
     ylim=c(0,0.5), xlab="", ylab="Probability of Occurrence", cex.lab=1.2)
title(sub=list(expression(paste("Bark thickness (mm at 20 cm DBH)")), cex=1.2),line=2.5)
title(sub="\n(scaled to unit variance below, log scale above)",line=3.5)
title(main="D", adj=0)
polygon(x=c(seq(min(newdat$slog_bark),max(newdat$slog_bark),length.out=100),
            rev(seq(min(newdat$slog_bark),max(newdat$slog_bark),length.out=100))),
        y=c(newdat$seLow[1:100], rev(newdat$seHigh[1:100])),
        col=adjustcolor("darkgreen",alpha.f=0.5), border=NA)
lines(newdat$slog_bark[1:100], newdat$ProbOccurrence[1:100], col="darkgreen")
polygon(x=c(seq(min(newdat$slog_bark),max(newdat$slog_bark),length.out=100),
            rev(seq(min(newdat$slog_bark),max(newdat$slog_bark),length.out=100))),
        y=c(newdat$seLow[101:200], rev(newdat$seHigh[101:200])),
        col=adjustcolor("gold",alpha.f=0.5), border=NA)
lines(newdat$slog_bark[101:200], newdat$ProbOccurrence[101:200],col="gold")
axis(3, at=xscale$sd[5:7], labels=xscale$bark[5:7], line=0)
title("",line=2.5)
legend("topright", legend=c("Low Precip Seasonality (2.5th %)","High Precip Seasonality (97.5th %)"), lty=1,
       col=c("darkgreen","gold"), lwd=3, bty="n", cex=1)
text(-1.5,0.3,"Low",col="darkgreen", cex=1.4)
text(0.7,0.3,"High",col="gold", cex=1.4)
text(0,0.4,"TxE: P < 0.0001")


# FRI x Bark
m15$model$slog_bark_unscaled <- unscale(m15$model$slog_bark)
m15$model$slog_fri_unscaled <- unscale(m15$model$slog_fri)
l1 <- lm(slog_fri_unscaled ~ as.numeric(slog_fri), data=m15$model)
l2 <- lm(slog_bark_unscaled ~ as.numeric(slog_bark), data=m15$model)
newdat <- data.frame(slog_fri=seq(-5,5,1), slog_bark=seq(-5,5,1))
xscale <- data.frame(sd=seq(-5,5,1), fri=round(exp(predict(l1,newdat)),1), bark=round(exp(predict(l2,newdat)),1))

### Plot TxE interactions with confidence intervals
newdat <- data.frame(
  slog_bark = rep(seq(quantile(m15$model$slog_bark,0.01),quantile(m15$model$slog_bark,0.9999),length.out = 100), 2),
  slog_fri = c(rep(quantile(m15$model$slog_fri, 0.025), 100), rep(quantile(m15$model$slog_fri, 0.975), 100)),
  lab = factor(rep(paste0(c("low.", "high."), "fri"), each = 100))) %>%
  mutate(`slog_fri:slog_bark` = slog_fri * slog_bark)
mm <- model.matrix( ~ slog_fri + slog_bark + slog_fri:slog_bark, newdat)
newdat$eta <- mm %*% coef(m15)[1:4]
newdat$ProbOccurrence <- exp(newdat$eta) / (1 + exp(newdat$eta))
newdat$VarPred <- diag(mm %*% vcov(m15)[1:4,1:4] %*% t(mm))
newdat$sePred <- sqrt(newdat$VarPred)
newdat$seLow <-
  exp(newdat$eta - 1.96 * newdat$sePred) / (1 + exp(newdat$eta - 1.96 * newdat$sePred))
newdat$seHigh <-
  exp(newdat$eta + 1.96 * newdat$sePred) / (1 + exp(newdat$eta + 1.96 * newdat$sePred))

plot(newdat$slog_bark[1:100], newdat$ProbOccurrence[1:100], type="l", col="orangered",
     ylim=c(0,0.35), xlab="", ylab="Probability of Occurrence", cex.lab=1.2)
title(sub=list(expression(paste("Bark thickness (mm at 20 cm DBH)")), cex=1.2),line=2.5)
title(sub="\n(scaled to unit variance below, log scale above)",line=3.5)
title(main="E", adj=0)
polygon(x=c(seq(min(newdat$slog_bark),max(newdat$slog_bark),length.out=100),
            rev(seq(min(newdat$slog_bark),max(newdat$slog_bark),length.out=100))),
        y=c(newdat$seLow[1:100], rev(newdat$seHigh[1:100])),
        col=adjustcolor("orangered",alpha.f=0.5), border=NA)
lines(newdat$slog_bark[1:100], newdat$ProbOccurrence[1:100], col="orangered")
polygon(x=c(seq(min(newdat$slog_bark),max(newdat$slog_bark),length.out=100),
            rev(seq(min(newdat$slog_bark),max(newdat$slog_bark),length.out=100))),
        y=c(newdat$seLow[101:200], rev(newdat$seHigh[101:200])),
        col=adjustcolor("blue",alpha.f=0.5), border=NA)
lines(newdat$slog_bark[101:200], newdat$ProbOccurrence[101:200],col="blue")
axis(3, at=xscale$sd[5:7], labels=xscale$bark[5:7], line=0)
title("",line=2.5)
legend("topright", legend=c("Frequent fire (2.5th %)","Infrequent fire (97.5th %)"), lty=1,
       col=c("orangered","blue"), lwd=3, bty="n", cex=1)
text(-1,0.28,"Infrequent fire",col="blue", cex=1.4)
text(0.5,0.25,"Frequent fire",col="orangered", cex=1.4)
text(0,0.3,"TxE: P < 0.0001")

# resprout and FRI
newdat <- data.frame(
  resprout = c(rep(0,1),rep(1,1),rep(0,1),rep(1,1)),
  slog_fri = c(rep(quantile(m16$model$slog_fri, 0.025), 2), rep(quantile(m16$model$slog_fri, 0.975), 2)),
  lab = factor(rep(paste0(c("low.", "high."), "fri"), each = 2))) %>%
  mutate(`slog_fri:resprout` = slog_fri * resprout)
mm <- model.matrix( ~ slog_fri + resprout + slog_fri:resprout, newdat)
newdat$eta <- mm %*% coef(m16)[1:4]
newdat$ProbOccurrence <- exp(newdat$eta) / (1 + exp(newdat$eta))
newdat$VarPred <- diag(mm %*% vcov(m16)[1:4,1:4] %*% t(mm))
newdat$sePred <- sqrt(newdat$VarPred)
newdat$seLow <-
  exp(newdat$eta - 1.96 * newdat$sePred) / (1 + exp(newdat$eta - 1.96 * newdat$sePred))
newdat$seHigh <-
  exp(newdat$eta + 1.96 * newdat$sePred) / (1 + exp(newdat$eta + 1.96 * newdat$sePred))

# alter new dataframe for plotting
newdat$resprout <- c("No","Yes","No","Yes")
newdat$lab <- c("Frequent fire","Frequent fire","Infrequent fire","Infrequent fire")
# Set colors manually for groups
cols <- c("Frequent fire" = "orangered", "Infrequent fire" = "blue")
# Set dodging offset
dodge <- 0.2
# Get unique x-values
xvals <- as.numeric(as.factor(newdat$resprout))
labs_levels <- levels(as.factor(newdat$lab))
# Compute x-positions: left/right dodge for groups
xpos <- xvals + (as.numeric(as.factor(newdat$lab)) - 1.5) * dodge
# Plot empty plot with range settings
plot(NA, xlim = range(xvals) + c(-0.5, 0.5), ylim = range(c(newdat$seLow, newdat$seHigh), na.rm = TRUE),
     xlab = "Resprouter", ylab = "Probability of Occurrence", xaxt = "n", cex.lab=1.2)
title(main="F", adj=0)
axis(1, at = xvals, labels = c(levels(as.factor(newdat$resprout)),levels(as.factor(newdat$resprout)) ))
# Points
points(xpos, newdat$ProbOccurrence,
       col = cols[newdat$lab], pch = 16, cex = 1.5)
# Error bars
arrows(xpos, newdat$seLow, xpos, newdat$seHigh,
       length = 0.05, angle = 90, code = 3, col = cols[newdat$lab], lwd = 2)
# Legend
legend("topright",legend = c("Frequent fire (2.5th %)","Infrequent fire (97.5th %)"),col = cols,
       pch = 16, title = "Fire Return Interval", bty = "n")
text(2,0.15,"TxE: P < 0.0001")

dev.off()
#### End of figure ########
