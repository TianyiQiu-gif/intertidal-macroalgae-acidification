# Meta_Analysis.R
# -------------------------------------------------------------------------
# Random-effects meta-analysis of macroalgal responses to ocean
# acidification. Effect sizes were calculated as the log response ratio
# (lnRR) and cumulative effect sizes were estimated using restricted
# maximum likelihood (REML) with the "metafor" package (Viechtbauer,
# 2010). Forest plots were generated using the "forestplot" package.
#
# Input:
#   meta_analysis_input.csv — comma-separated file with columns:
#     Study           — first author and year (character)
#     Species         — macroalgal species name (character)
#     Phylum          — Rhodophyta / Chlorophyta / Phaeophyceae
#     Latitude        — study site latitude (decimal degrees)
#     Longitude       — study site longitude (decimal degrees)
#     Trait_Category   — Growth / Photosynthesis / Oxidative_stress /
#                        Carbon_fixation / Nitrogen_assimilation
#     Trait           — specific measured variable (character)
#     Mean_Control    — mean of the control group (numeric)
#     SD_Control      — standard deviation of the control group (numeric)
#     N_Control       — sample size of the control group (integer)
#     Mean_Treatment  — mean of the treatment group (numeric)
#     SD_Treatment    — standard deviation of the treatment group (numeric)
#     N_Treatment     — sample size of the treatment group (integer)
#
# Output:
#   meta_analysis_overall.csv       — overall pooled effect size
#   meta_analysis_by_trait.csv      — pooled effect sizes by trait category
#   meta_analysis_individual.csv    — per-observation effect sizes
#   forest_plot_overall.pdf         — forest plot of all observations
#   forest_plot_by_trait.pdf        — forest plots stratified by trait
#   funnel_plot.pdf                 — funnel plot for publication bias
#   study_site_map.pdf              — geographic distribution of study sites
#
# Dependencies:
#   install.packages(c("metafor", "forestplot", "ggplot2", "sf",
#                       "rnaturalearth", "rnaturalearthdata"))
# -------------------------------------------------------------------------

library(metafor)
library(forestplot)
library(ggplot2)

# =========================================================================
# 1. Read extracted data
# =========================================================================
dat <- read.csv("meta_analysis_input.csv", stringsAsFactors = FALSE)

# =========================================================================
# 2. Calculate effect sizes (log response ratio, lnRR)
# =========================================================================
es <- escalc(
  measure = "ROM",
  m1i = Mean_Treatment, sd1i = SD_Treatment, n1i = N_Treatment,
  m2i = Mean_Control,   sd2i = SD_Control,   n2i = N_Control,
  data = dat
)

# =========================================================================
# 3. Overall random-effects model (REML)
# =========================================================================
model_overall <- rma(yi, vi, data = es, method = "REML")

overall_result <- data.frame(
  Category      = "Overall",
  k             = model_overall$k,
  lnRR          = round(model_overall$beta[1], 4),
  CI_lower      = round(model_overall$ci.lb, 4),
  CI_upper      = round(model_overall$ci.ub, 4),
  p_value       = round(model_overall$pval, 6),
  I2            = round(model_overall$I2, 2),
  tau2          = round(model_overall$tau2, 6)
)

write.csv(overall_result, "meta_analysis_overall.csv", row.names = FALSE)

# =========================================================================
# 4. Subgroup analysis by trait category
# =========================================================================
trait_categories <- unique(es$Trait_Category)
trait_results    <- data.frame()

for (tc in trait_categories) {
  sub <- es[es$Trait_Category == tc, ]
  if (nrow(sub) < 2) next
  mod <- rma(yi, vi, data = sub, method = "REML")
  trait_results <- rbind(trait_results, data.frame(
    Category = tc,
    k        = mod$k,
    lnRR     = round(mod$beta[1], 4),
    CI_lower = round(mod$ci.lb, 4),
    CI_upper = round(mod$ci.ub, 4),
    p_value  = round(mod$pval, 6),
    I2       = round(mod$I2, 2),
    tau2     = round(mod$tau2, 6)
  ))
}

write.csv(trait_results, "meta_analysis_by_trait.csv", row.names = FALSE)

# =========================================================================
# 5. Export per-observation effect sizes
# =========================================================================
individual <- data.frame(
  Study          = es$Study,
  Species        = es$Species,
  Phylum         = es$Phylum,
  Trait_Category = es$Trait_Category,
  Trait          = es$Trait,
  lnRR           = round(es$yi, 4),
  Variance       = round(es$vi, 6)
)

write.csv(individual, "meta_analysis_individual.csv", row.names = FALSE)

# =========================================================================
# 6. Forest plot — overall
# =========================================================================
pdf("forest_plot_overall.pdf", width = 10, height = max(6, 0.35 * nrow(es)))
forest(
  model_overall,
  slab  = paste0(es$Study, " — ", es$Species, " (", es$Trait, ")"),
  xlab  = "Log Response Ratio (lnRR)",
  cex   = 0.7,
  header = TRUE
)
dev.off()

# =========================================================================
# 7. Forest plots — stratified by trait category
# =========================================================================
pdf("forest_plot_by_trait.pdf", width = 10, height = 8)
for (tc in trait_categories) {
  sub <- es[es$Trait_Category == tc, ]
  if (nrow(sub) < 2) next
  mod <- rma(yi, vi, data = sub, method = "REML")
  forest(
    mod,
    slab  = paste0(sub$Study, " — ", sub$Species),
    xlab  = "Log Response Ratio (lnRR)",
    main  = tc,
    cex   = 0.75,
    header = TRUE
  )
}
dev.off()

# =========================================================================
# 8. Funnel plot for publication bias assessment
# =========================================================================
pdf("funnel_plot.pdf", width = 7, height = 6)
funnel(model_overall, xlab = "Log Response Ratio (lnRR)")
dev.off()

# =========================================================================
# 9. Geographic distribution of study sites
# =========================================================================
if (requireNamespace("sf", quietly = TRUE) &
    requireNamespace("rnaturalearth", quietly = TRUE)) {

  library(sf)
  library(rnaturalearth)

  world <- ne_countries(scale = "medium", returnclass = "sf")
  sites <- unique(dat[, c("Study", "Latitude", "Longitude")])

  p <- ggplot() +
    geom_sf(data = world, fill = "grey90", colour = "grey60", linewidth = 0.2) +
    geom_point(data = sites, aes(x = Longitude, y = Latitude),
               colour = "#D32F2F", size = 2.5, alpha = 0.8) +
    coord_sf(expand = FALSE) +
    labs(x = "Longitude", y = "Latitude") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major = element_line(colour = "grey85", linewidth = 0.3))

  ggsave("study_site_map.pdf", plot = p, width = 10, height = 5.5)
}

# =========================================================================
# 10. Heterogeneity test and Egger's regression (publication bias)
# =========================================================================
sink("meta_analysis_diagnostics.txt")

cat("=== Overall model summary ===\n")
print(summary(model_overall))

cat("\n=== Heterogeneity (Q-test) ===\n")
cat(sprintf("Q = %.4f, df = %d, p = %.6f\n",
            model_overall$QE, model_overall$k - 1, model_overall$QEp))

cat(sprintf("I-squared = %.2f%%\n", model_overall$I2))
cat(sprintf("tau-squared = %.6f\n", model_overall$tau2))

cat("\n=== Egger's regression test for funnel plot asymmetry ===\n")
reg_test <- regtest(model_overall)
print(reg_test)

sink()
