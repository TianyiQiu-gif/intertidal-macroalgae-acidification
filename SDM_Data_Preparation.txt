# SDM_Data_Preparation.R
# -------------------------------------------------------------------------
# Preparation of species occurrence records and environmental layers for
# habitat-suitability modelling with MaxEnt. This script performs:
#   (1) batch download of occurrence data from GBIF via the rgbif package;
#   (2) automated cleaning (coordinate validation, duplicate removal,
#       spatial thinning) following best practices for species distribution
#       modelling (Zizka et al., 2019);
#   (3) multicollinearity screening of Bio-ORACLE environmental layers
#       using Pearson correlation and variance inflation factor (VIF)
#       analysis to produce a refined variable set.
#
# Input:
#   species_list.csv        — CSV with column "Species" listing binomial
#                             names of target macroalgae (one per row)
#   env_layers/             — directory containing Bio-ORACLE GeoTIFF
#                             rasters (current and future SSP scenarios)
#
# Output:
#   cleaned_occurrences/    — one CSV per species with cleaned coordinates
#   occurrence_summary.csv  — record counts before and after cleaning
#   correlation_matrix.pdf  — pairwise Pearson correlation heat map
#   vif_results.csv         — VIF values for all candidate variables
#   selected_variables.csv  — retained variable names after screening
#
# Dependencies:
#   install.packages(c("rgbif", "CoordinateCleaner", "spThin",
#                       "terra", "corrplot", "usdm"))
# -------------------------------------------------------------------------

library(rgbif)
library(CoordinateCleaner)
library(spThin)
library(terra)
library(corrplot)

# =========================================================================
# 1. Read species list
# =========================================================================
species_list <- read.csv("species_list.csv", stringsAsFactors = FALSE)
species      <- species_list$Species

dir.create("cleaned_occurrences", showWarnings = FALSE)

# =========================================================================
# 2. Batch download and clean GBIF occurrence data
# =========================================================================
summary_records <- data.frame(
  Species      = character(),
  Raw          = integer(),
  After_Clean  = integer(),
  After_Thin   = integer(),
  stringsAsFactors = FALSE
)

for (sp in species) {

  cat(sprintf("\n=== Processing: %s ===\n", sp))

  # --- 2a. Download from GBIF ---
  gbif_data <- occ_search(
    scientificName = sp,
    hasCoordinate  = TRUE,
    limit          = 10000
  )$data

  if (is.null(gbif_data) || nrow(gbif_data) == 0) {
    cat(sprintf("  No records found for %s. Skipped.\n", sp))
    next
  }

  raw_n <- nrow(gbif_data)
  cat(sprintf("  Raw records: %d\n", raw_n))

  # --- 2b. Standardise columns ---
  occ <- data.frame(
    species          = gbif_data$species,
    decimalLongitude = gbif_data$decimalLongitude,
    decimalLatitude  = gbif_data$decimalLatitude,
    countryCode      = gbif_data$countryCode,
    basisOfRecord    = gbif_data$basisOfRecord
  )
  occ <- occ[complete.cases(occ$decimalLongitude, occ$decimalLatitude), ]

  # --- 2c. Coordinate cleaning (CoordinateCleaner) ---
  occ_clean <- clean_coordinates(
    occ,
    lon = "decimalLongitude",
    lat = "decimalLatitude",
    species = "species",
    tests = c(
      "capitals",      # remove records at country capitals
      "centroids",     # remove records at country centroids
      "equal",         # remove records with equal lon/lat
      "gbif",          # remove GBIF headquarters
      "institutions",  # remove biodiversity institutions
      "zeros",         # remove (0,0) coordinates
      "seas"           # flag terrestrial records for marine species
    ),
    seas_ref = NULL     # uses default; set to a land polygon if needed
  )

  occ_flagged <- occ[occ_clean$.summary, ]
  clean_n     <- nrow(occ_flagged)
  cat(sprintf("  After cleaning: %d\n", clean_n))

  if (clean_n < 5) {
    cat(sprintf("  Too few records for %s after cleaning. Skipped.\n", sp))
    next
  }

  # --- 2d. Spatial thinning (10 km minimum distance) ---
  thin_result <- thin(
    loc.data      = occ_flagged,
    lat.col       = "decimalLatitude",
    long.col      = "decimalLongitude",
    spec.col      = "species",
    thin.par      = 10,     # minimum distance in km
    reps          = 5,
    locs.thinned.list.return = TRUE,
    write.files   = FALSE,
    write.log.file = FALSE
  )

  thinned <- thin_result[[which.max(sapply(thin_result, nrow))]]
  thin_n  <- nrow(thinned)
  cat(sprintf("  After thinning: %d\n", thin_n))

  # --- 2e. Export cleaned occurrence file ---
  sp_tag <- gsub(" ", "_", sp)
  write.csv(
    thinned,
    file      = file.path("cleaned_occurrences", paste0(sp_tag, ".csv")),
    row.names = FALSE
  )

  summary_records <- rbind(summary_records, data.frame(
    Species     = sp,
    Raw         = raw_n,
    After_Clean = clean_n,
    After_Thin  = thin_n
  ))
}

write.csv(summary_records, "occurrence_summary.csv", row.names = FALSE)

# =========================================================================
# 3. Environmental variable selection
# =========================================================================

# --- 3a. Load all candidate Bio-ORACLE layers (current period) ---
env_files <- list.files("env_layers", pattern = "\\.tif$", full.names = TRUE)
env_stack <- rast(env_files)

cat(sprintf("\nCandidate environmental variables: %d\n", nlyr(env_stack)))

# --- 3b. Sample values at occurrence locations (all species pooled) ---
all_occ_files <- list.files("cleaned_occurrences", pattern = "\\.csv$",
                            full.names = TRUE)
all_coords <- do.call(rbind, lapply(all_occ_files, function(f) {
  d <- read.csv(f)
  data.frame(lon = d[,1], lat = d[,2])
}))

env_values <- extract(env_stack, vect(all_coords, geom = c("lon","lat"),
                                       crs = "EPSG:4326"))
env_values <- env_values[complete.cases(env_values), -1]  # drop ID column

# --- 3c. Pearson correlation matrix ---
cor_matrix <- cor(env_values, use = "pairwise.complete.obs")

pdf("correlation_matrix.pdf", width = 10, height = 10)
corrplot(cor_matrix, method = "color", type = "lower", tl.cex = 0.7,
         addCoef.col = "black", number.cex = 0.5,
         title = "Pairwise Pearson Correlation of Environmental Variables",
         mar = c(0, 0, 2, 0))
dev.off()

# --- 3d. Remove highly correlated pairs (|r| >= 0.7) ---
cor_threshold <- 0.7
high_cor <- findCorrelation(cor_matrix, cutoff = cor_threshold, names = TRUE)
# Note: findCorrelation is from the caret package; alternative below:
# Identify pairs exceeding threshold and iteratively remove the variable
# with the highest mean absolute correlation.
if (!requireNamespace("caret", quietly = TRUE)) {
  drop_vars <- character()
  tmp <- abs(cor_matrix)
  diag(tmp) <- 0
  while (max(tmp) >= cor_threshold) {
    mean_cor <- colMeans(tmp)
    worst    <- names(which.max(mean_cor))
    drop_vars <- c(drop_vars, worst)
    tmp <- tmp[!rownames(tmp) %in% worst, !colnames(tmp) %in% worst,
               drop = FALSE]
  }
  high_cor <- drop_vars
}

retained_vars <- setdiff(names(env_values), high_cor)
cat(sprintf("Variables removed (|r| >= %.1f): %s\n",
            cor_threshold, paste(high_cor, collapse = ", ")))
cat(sprintf("Variables retained: %s\n", paste(retained_vars, collapse = ", ")))

# --- 3e. VIF analysis on retained variables ---
if (requireNamespace("usdm", quietly = TRUE)) {
  library(usdm)
  env_retained <- env_values[, retained_vars, drop = FALSE]
  vif_result   <- vif(env_retained)

  write.csv(vif_result, "vif_results.csv", row.names = FALSE)
  cat("\nVIF values:\n")
  print(vif_result)

  # Exclude variables with VIF > 10
  final_vars <- vif_result$Variables[vif_result$VIF < 10]
} else {
  final_vars <- retained_vars
}

write.csv(data.frame(Variable = final_vars), "selected_variables.csv",
          row.names = FALSE)

cat(sprintf("\nFinal selected variables (%d): %s\n",
            length(final_vars), paste(final_vars, collapse = ", ")))

# =========================================================================
# 4. Export refined environmental layers for MaxEnt
# =========================================================================
dir.create("env_selected", showWarnings = FALSE)

# Current layers
for (v in final_vars) {
  idx <- which(names(env_stack) == v)
  if (length(idx) == 1) {
    writeRaster(env_stack[[idx]],
                file.path("env_selected", paste0(v, ".tif")),
                overwrite = TRUE)
  }
}

# Future layers (repeat for each SSP scenario directory if structured as
# env_layers_ssp119/ and env_layers_ssp585/)
future_dirs <- c("env_layers_ssp119", "env_layers_ssp585")

for (fd in future_dirs) {
  if (!dir.exists(fd)) next
  out_dir <- paste0("env_selected_", basename(fd))
  dir.create(out_dir, showWarnings = FALSE)

  fut_files <- list.files(fd, pattern = "\\.tif$", full.names = TRUE)
  fut_stack <- rast(fut_files)

  for (v in final_vars) {
    idx <- which(names(fut_stack) == v)
    if (length(idx) == 1) {
      writeRaster(fut_stack[[idx]],
                  file.path(out_dir, paste0(v, ".tif")),
                  overwrite = TRUE)
    }
  }
}

cat("\nRefined environmental layers exported to env_selected/\n")
cat("Data preparation complete.\n")
