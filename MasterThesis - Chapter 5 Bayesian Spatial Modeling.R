#install.packages("leaflet")
#install.packages("ggplot2")
#install.packages(c("gstat", "sf"))
#install.packages("rgal")
#install.packages("spdep")
#install.packages("spBayes")
library(leaflet)
library(ggplot2)
library(dplyr)
library(sf)
library(sp)
library(gstat)
library(raster)
library(rgdal)
library(spdep)
library(spBayes)
library(gridExtra)
library(coda)
library(tidyr)



#Code to subset the file to 500 rows
#spatial_data <- read.csv("C:/Users/018647575SA/OneDrive - csulb/Desktop/housing_edited.csv")
#names(spatial_data)
#set.seed(84323)  
#k <- min(500, nrow(spatial_data))               
#idx <- sample.int(n = nrow(spatial_data), size = k, replace = FALSE)
#spatial_500 <- spatial_data[idx, , drop = FALSE]

# Save 
#write.csv(spatial_500, "C:/Users/018647575SA/OneDrive - csulb/Desktop/housing_edited_subset_500.csv", row.names = FALSE)

# quick checks
#nrow(spatial_500)
#head(spatial_500)


################################
s_data <- read.csv("C:/Users/garci/OneDrive/Desktop/housing_edited_subset_500.csv")
names(s_data)

###########################################################
#DESCRIPTIVE STATISTICS
library(leaflet)

#OpenStreetMap of the Iowa Housing Area
leaflet(s_data) %>%
  addTiles() %>%
  addCircleMarkers(
    lng = ~Longitude,
    lat = ~Latitude,
    radius = 4,
    color = "blue",
    stroke = FALSE,
    fillOpacity = 0.5
  )



# Relative Frequency Bar Graph of Building Type
data_pct <- s_data %>%
  group_by(Bldg_Type) %>%
  summarise(n = n()) %>%
  mutate(percent = n / sum(n) * 100)

# Plot
ggplot(data_pct, aes(x = Bldg_Type, y = percent)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    title = "Percentage of Bldg_Type Categories",
    x = "Bldg_Type",
    y = "Percentage (%)"
  ) +
  geom_text(aes(label = paste0(round(percent, 1), "%")), vjust = -0.5) +
  theme_minimal()

#########################################################################
# STEP 1: Remove rows with missing values from Longitude, Latitude, and Sale Price
data_clean <- subset(s_data, !is.na(Longitude) & !is.na(Latitude) & !is.na(Sale_Price))


# STEP 2: Convert to sf object
# Converts the cleaned dataframe into a simple features (sf) spatial object using Longitude and Latitude as coordinates
# CRS 4326 means WGS84 Geographic Coordinate System (degrees of latitude/longitude)

data_sf <- st_as_sf(
  data_clean,
  coords = c("Longitude", "Latitude"),
  crs = 4326 # WGS84
)

# STEP 3:  Project to meters (for variogram distance)
# Transforms the spatial data from geographic coordinates (Degrees) to projected coordinates in meters using EPSG:3857 (Web Mercator)
# Variogram distance computations require units in meters not degrees
data_sf_proj <- st_transform(data_sf, crs = 3857)

# STEP 4: Convert to Spatial object
# Converts the sf object to a spatial object using the sp packages
# Variogram operates on Spatial objects, not sf 
data_sp <- as(data_sf_proj, "Spatial")

# STEP 5: Create gstat object
# This creates a gstat model for the response variable, Sale_Price, using a simple constant mean model (~1)
# Setting up the data for the variogram calculation.
gs <- gstat::gstat(
  formula = Sale_Price ~ 1,
  data = data_sp
)


# STEP 6A: Compute empirical variogram
# Quantifies how the spatial dependence (semivariance) of Sale Price changes with distance
# Estimeates empirical semivariance between points at different spatial lags
emp_variogram <- variogram(gs)

# STEP 6B: View variogram
# Distance (x-axis), Semivariance (y-axis)
plot(emp_variogram, main = "Empirical Variogram of Sale Price")



# STEP 7A: Fit a theoretical variogram model
# Fits an exponential theoretical model to the empirical variogram.
# Could also try "Sph" (spherical) or "Gau"(Gaussian)
vgm_model <- fit.variogram(emp_variogram, model = vgm("Exp"))

# STEP 7B: Plot with fitted model
# Overlays the fitted model on top of the empirical variogram to compare
plot(emp_variogram, model = vgm_model, main = "Variogram of Sale Price (Fitted Model)")

vgm_model

########
#Includes Matern

# Refit variogram models (now including Matérn)
vgm_exp <- fit.variogram(emp_variogram, model = vgm("Exp"))
vgm_sph <- fit.variogram(emp_variogram, model = vgm("Sph"))
vgm_gau <- fit.variogram(emp_variogram, model = vgm("Gau"))
vgm_mat <- fit.variogram(emp_variogram, model = vgm("Mat"))  # Matérn

# Convert variogram to data frame for plotting
emp_df <- as.data.frame(emp_variogram)

# Function to generate a ggplot variogram model
plot_variogram_model <- function(model, model_name, color) {
  model_line <- variogramLine(model, maxdist = max(emp_df$dist), n = 100)
  ggplot() +
    geom_point(data = emp_df, aes(x = dist, y = gamma)) +
    geom_line(data = model_line, aes(x = dist, y = gamma), color = color, size = 1.2) +
    labs(title = model_name, x = "Distance", y = "Semivariance") +
    theme_minimal()
}

# Create plots
p1 <- plot_variogram_model(vgm_mat, "Matérn", "purple")
p2 <- plot_variogram_model(vgm_exp, "Exponential", "blue")
p3 <- plot_variogram_model(vgm_gau, "Gaussian", "red")
p4 <- plot_variogram_model(vgm_sph, "Spherical", "forestgreen")


# Display all four side by side
grid.arrange(p1, p2, p3, p4, ncol = 2)

# Print fitted models
vgm_mat
vgm_exp
vgm_gau
vgm_sph

##########################################################
# Determining best fit via RSS

# ---- Show the three fitted models side-by-side ----
grid.arrange(p1, p2, p3, ncol = 3)

# ---- Model comparison table (lower RSS is better) ----
# Helper to extract parameters cleanly from a gstat variogram model
extract_params <- function(vgm_obj) {
  # rows: "Nug", then the structural model (Exp/Sph/Gau)
  nugget_row <- vgm_obj[vgm_obj$model == "Nug", , drop = FALSE]
  struct_row <- vgm_obj[vgm_obj$model != "Nug", , drop = FALSE]
  nugget     <- if (nrow(nugget_row)) nugget_row$psill[1] else 0
  psill      <- struct_row$psill[1]
  range_raw  <- struct_row$range[1]
  model_name <- struct_row$model[1]
#############################################





# Practical range helper for Matérn
practical_range_matern <- function(vgm_obj, maxdist = max(emp_df$dist)) {
  mod_line <- variogramLine(vgm_obj, maxdist = maxdist, n = 800)
  nug  <- vgm_obj$psill[vgm_obj$model == "Nug"][1]
  psil <- vgm_obj$psill[vgm_obj$model != "Nug"][1]
  target <- nug + 0.95 * psil
  hit <- mod_line$dist[which(mod_line$gamma >= target)[1]]
  if (is.na(hit)) NA_real_ else hit
}

# Extract parameters robustly
extract_params <- function(vgm_obj) {
  nug_row <- vgm_obj[vgm_obj$model == "Nug", , drop = FALSE]
  str_row <- vgm_obj[vgm_obj$model != "Nug", , drop = FALSE]

  nugget     <- if (nrow(nug_row)) nug_row$psill[1] else 0
  psill      <- str_row$psill[1]
  range_raw  <- str_row$range[1]
  model_name <- as.character(str_row$model[1])
  kappa_val  <- if ("kappa" %in% names(str_row)) str_row$kappa[1] else NA_real_

  practical_range <- switch(
    model_name,
    "Sph" = range_raw,
    "Exp" = 3 * range_raw,
    "Gau" = sqrt(3) * range_raw,
    "Mat" = practical_range_matern(vgm_obj, max(emp_df$dist)),
    range_raw
  )

  data.frame(
    Nugget = nugget,
    Partial_Sill = psill,
    Total_Sill = nugget + psill,
    Range_Raw = range_raw,
    Kappa = kappa_val,
    Practical_Range = practical_range
  )
}

# Comparison table
cmp <- rbind(
  cbind(Model = "Matérn",      extract_params(vgm_mat), RSS = attr(vgm_mat, "SSErr")),
  cbind(Model = "Exponential", extract_params(vgm_exp), RSS = attr(vgm_exp, "SSErr")),
  cbind(Model = "Gaussian",    extract_params(vgm_gau), RSS = attr(vgm_gau, "SSErr")),
  cbind(Model = "Spherical",   extract_params(vgm_sph), RSS = attr(vgm_sph, "SSErr"))
)

cmp <- cmp[order(cmp$RSS), ]
print(cmp, row.names = FALSE)

cat("\nBest variogram model by RSS:", cmp$Model[1], "\n")







#################################################
# ORDINARY KRIGING
# Spatially interpolate Sale Price using the fitted variogram model

# Convert SpatialPointsDataFrame to data frame
data_df <- as.data.frame(data_sp)
data_df$X <- coordinates(data_sp)[,1]
data_df$Y <- coordinates(data_sp)[,2]

# Aggregate by unique coordinate pairs
data_clean <- data_df %>%
  group_by(X, Y) %>%
  summarise(
    Sale_Price = mean(Sale_Price),
    .groups = 'drop'
  )

# Convert back to spatial
coordinates(data_clean) <- ~ X + Y
proj4string(data_clean) <- proj4string(data_sp)



###
# Create prediction grid

grid <- spsample(data_clean, type = "regular", cellsize = 200)
gridded(grid) <- TRUE
proj4string(grid) <- proj4string(data_clean)
grid_sp <- grid  # Needed for kriging


# Kriging
kriging_result <- krige(
  Sale_Price ~ 1,
  locations = data_clean,
  newdata = grid_sp,
  model = vgm_model,
  nmax = 30,
  maxdist = 5000
)


st_bbox(data_sf_proj)
summary(kriging_result)

mean(data_clean$Sale_Price, na.rm = TRUE)
##################################
#test run
# Extract the predicted local means
local_means <- kriging_result$var1.pred

# Quick summary
summary(local_means)

# Average of all local means (should be close to global mean, but smoothed)
mean(local_means, na.rm = TRUE)

# Convert to raster for mapping
library(raster)
kriging_raster <- raster(kriging_result["var1.pred"])
plot(kriging_raster, main = "Ordinary Kriging Local Means")
points(data_clean, pch = 16, cex = 0.5)





####################################

#ALTERNATIVE PLOT
#install.packages("rgal")
#library(rgdal)

# Step 1: Convert kriging predictions to raster
r_pred <- rasterFromXYZ(kriging_df[, c("x", "y", "var1.pred")])

# Step 2: Assign CRS (EPSG:3857 — Web Mercator, as used earlier)
crs(r_pred) <- CRS("+init=epsg:3857")

# Step 3: Reproject raster to EPSG:4326 for leaflet (lat/lon)
r_pred_wgs84 <- projectRaster(r_pred, crs = CRS("+init=epsg:4326"))

# Step 4: Convert raster to leaflet-compatible object
pal <- colorNumeric("plasma", values(r_pred_wgs84), na.color = "transparent")

leaflet() %>%
  addTiles() %>%  # Default OpenStreetMap basemap
  addRasterImage(r_pred_wgs84, colors = pal, opacity = 0.6) %>%
  addLegend(pal = pal, values = values(r_pred_wgs84),
            title = "Predicted Sale Price") %>%
  setView(lng = mean(data$Longitude, na.rm = TRUE),
          lat = mean(data$Latitude, na.rm = TRUE),
          zoom = 12)



#MORAN'S I TEST
# Step 1: Prepare coordinates
coords_obs <- coordinates(data_clean)

# Step 2: Create neighbors list using 8 nearest neighbors
nb_obs <- knn2nb(knearneigh(coords_obs, k = 8))

# Step 3: Convert to spatial weights list
listw_obs <- nb2listw(nb_obs, style = "W")

# Step 4: Moran's I test
moran_obs <- moran.test(data_clean$Sale_Price, listw_obs)
print(moran_obs)


#Global Moran's I Monte Carlo Test
library(spdep)

# Step 1: Prepare coordinates
coords_obs <- coordinates(data_clean)

# Step 2: Create neighbors list using 8 nearest neighbors
nb_obs <- knn2nb(knearneigh(coords_obs, k = 8))

# Step 3: Convert to spatial weights list
listw_obs <- nb2listw(nb_obs, style = "W")

# Step 4: Global Moran's I Monte Carlo test
set.seed(123)  # for reproducibility
moran_mc <- moran.mc(data_clean$Sale_Price, listw_obs, nsim = 999)

print(moran_mc)

# Monte Carlo Global Moran's I plot
hist(moran_mc$res, 
     breaks = 20, 
     main = "Monte Carlo Simulation of Global Moran's I",
     xlab = "Simulated Moran's I",
     col = "lightblue", 
     border = "white")

# Add observed Moran's I as a red vertical line
abline(v = moran_mc$statistic, col = "red", lwd = 2)

# Add legend
legend("topright", legend = c("Observed Moran's I"), 
       col = "red", lwd = 2)

########################################################
#GEARY'S C

geary_result <- geary.test(
  data_clean$Sale_Price,
  listw_obs
)

print(geary_result)


library(spdep)

# Step 1: Prepare coordinates
coords_obs <- coordinates(data_clean)

# Step 2: Create neighbors list using 8 nearest neighbors
nb_obs <- knn2nb(knearneigh(coords_obs, k = 8))

# Step 3: Convert to spatial weights list
listw_obs <- nb2listw(nb_obs, style = "W")

# Step 4: Geary's C Monte Carlo test (999 permutations)
set.seed(123)
geary_mc <- geary.mc(data_clean$Sale_Price, listw_obs, nsim = 999)

print(geary_mc)



hist(geary_mc$res, 
     breaks = 20,
     main = "Monte Carlo Simulation of Geary's C",
     xlab = "Simulated Geary's C",
     col = "lightblue", 
     border = "white")

# Add observed Geary's C as a red vertical line
abline(v = geary_mc$statistic, col = "red", lwd = 2)

# Add legend
legend("topright", legend = c("Observed Geary's C"), 
       col = "red", lwd = 2)










############################################################################################################
############################################################################################################
library(spBayes)
library(sf)

# Use the projected sf from your pipeline
# data_sf_proj: POINT geometry in meters, has Sale_Price column

dat    <- sf::st_drop_geometry(data_sf_proj)         # plain data.frame
coords <- sf::st_coordinates(data_sf_proj)           # matrix [n x 2] in meters

set.seed(1)
m <- spLM(
  Sale_Price ~ 1,
  data   = dat,
  coords = coords,
  starting = list(phi = 3e-4, sigma.sq = var(dat$Sale_Price)/2, tau.sq = var(dat$Sale_Price)/10),
  tuning   = list(phi = 0.05, sigma.sq = 0.05, tau.sq = 0.05),
  priors   = list(
    phi.Unif    = c(1e-5, 1e-2),   # in m^-1 (≈ practical range 300 m to 300 km)
    sigma.sq.IG = c(2, 2),
    tau.sq.IG   = c(2, 2)
  ),
  cov.model = "exponential",
  n.samples = 40000,
  verbose   = TRUE
)

# Recover posterior samples after burn-in/thin
burn <- 10000; thin <- 40
mrec <- spRecover(m, start = burn, thin = thin, verbose = FALSE)

# Quick posterior summaries
theta_samps <- mrec$p.theta.recover.samples[, c("sigma.sq","tau.sq","phi")]
colMeans(theta_samps)
apply(theta_samps, 2, quantile, probs = c(0.025, 0.5, 0.975))

# Optional: map phi to practical range (meters)
practical_range <- 3 / theta_samps[, "phi"]
mean(practical_range)



######################################################
#using this one in final thesis

# Data prep
dat    <- sf::st_drop_geometry(data_sf_proj)
coords <- sf::st_coordinates(data_sf_proj)

# Function to run one chain
run_chain <- function(seed) {
  set.seed(seed)
  spLM(
    Sale_Price ~ 1,
    data   = dat,
    coords = coords,
    starting = list(phi = 3e-4,
                    sigma.sq = var(dat$Sale_Price)/2,
                    tau.sq   = var(dat$Sale_Price)/10),
    tuning   = list(phi = 0.05,
                    sigma.sq = 0.05,
                    tau.sq   = 0.05),
    priors   = list(
      phi.Unif    = c(1e-5, 1e-2),
      sigma.sq.IG = c(2, 2),
      tau.sq.IG   = c(2, 2)
    ),
    cov.model = "exponential",
    n.samples = 40000,
    verbose   = FALSE
  )
}


# Seeds for each chain
seeds <- c(2025, 84531, 98765, 45678)

# Run 4 chains, each with its own seed
chains <- lapply(1:4, function(i) run_chain(seeds[i]))



# --- Recover samples for each chain ---
burn <- 10000; thin <- 40
mrec_list <- lapply(chains, function(m) {
  spRecover(m, start = burn, thin = thin, verbose = FALSE)
})

# --- Extract parameter draws (sigma, tau, phi) per chain ---
theta_list <- lapply(mrec_list, function(mrec) {
  as.mcmc(mrec$p.theta.recover.samples[, c("sigma.sq","tau.sq","phi")])
})

# --- Combine into coda::mcmc.list for diagnostics ---
theta_mcmc <- mcmc.list(theta_list)

summary(theta_mcmc)
gelman.diag(theta_mcmc)   # Gelman–Rubin R-hat
#traceplot(theta_mcmc)     # Trace plots





###########################
#test code for mean
library(spBayes)

## 1) Use the SAME formula as in spLM to build X
## If your spLM call was Sale_Price ~ 1, keep ~ 1; otherwise pull from the first chain
form_sp <- tryCatch(eval(chains[[1]]$call$formula), error = function(e) ~ 1)
X <- model.matrix(form_sp, data = dat)    # n × p

## 2) For each recovered chain, make iter × n matrices for xb and w, then add
mu_all <- do.call(rbind, lapply(mrec_list, function(mrec) {
  beta_mat <- as.matrix(mrec$p.beta.recover.samples)  # iter × p
  w_raw    <- as.matrix(mrec$p.w.recover.samples)     # either iter × n  OR  n × iter

  # Ensure w has iter rows, n columns
  w_mat <- if (nrow(w_raw) == nrow(beta_mat)) w_raw else t(w_raw)

  # Compute Xβ for all draws: (iter × p) %*% (p × n) = iter × n
  xb_mat <- beta_mat %*% t(X)

  # Posterior draws of E[y | params] at observed locations
  xb_mat + w_mat   # iter × n
}))

## 3) Posterior means of Sale_Price
sale_price_mean_per_obs <- colMeans(mu_all)     # length n (one per observation)
sale_price_overall_mean <- mean(sale_price_mean_per_obs)  # single number

# Inspect
print(head(sale_price_mean_per_obs))
print(sale_price_overall_mean)
###################################
#test run


library(spBayes)
library(coda)
library(sp)

# --- Settings ---
burn <- 10000; thin <- 40

# 1) θ = (sigma.sq, tau.sq, phi)  ---------------------------------------------
theta_list <- lapply(chains, function(fit) {
  keep <- seq(burn, nrow(fit$p.theta.samples), by = thin)
  as.mcmc(fit$p.theta.samples[keep, c("sigma.sq","tau.sq","phi")])
})
theta_mcmc <- mcmc.list(theta_list)
print(summary(theta_mcmc))       # quick θ summary
print(gelman.diag(theta_mcmc))   # R-hat

# 2) Bayesian mean (posterior mean of β0)  ------------------------------------
# spRecover returns β samples (and w), not θ
mrec_list <- lapply(chains, function(fit)
  spRecover(fit, start = burn, thin = thin, verbose = FALSE)
)

# Stack β samples across chains (intercept-only model => 1 column)
beta_all <- do.call(rbind, lapply(mrec_list, function(r) as.matrix(r$p.beta.recover.samples)))
bayesian_mean_beta0 <- mean(beta_all[, 1])
bayesian_ci_beta0   <- quantile(beta_all[, 1], c(.025, .975))

# Frequentist overall mean of observed Sale_Price for comparison
overall_mean <- mean(dat$Sale_Price, na.rm = TRUE)

cat("\n--- Means ---\n")
cat("Overall sample mean       :", overall_mean, "\n")
cat("Bayesian mean (beta_0)    :", bayesian_mean_beta0, "\n")
cat("Bayesian 95% CrI (beta_0) :", bayesian_ci_beta0[1], ", ", bayesian_ci_beta0[2], "\n")
cat("Difference (Bayes - overall):", bayesian_mean_beta0 - overall_mean, "\n")

# 3) OPTIONAL: Bayesian surface mean via posterior predictive kriging ----------
# Requires your prediction grid 'grid_sp'
if (exists("grid_sp")) {
  newcoords <- coordinates(grid_sp)  # n_pred x 2

  pred_list <- lapply(chains, function(fit)
    spPredict(fit, pred.coords = newcoords, start = burn, thin = thin, verbose = FALSE)
  )

  # Combine across chains: each is n_pred x n_draws
  pred_all <- do.call(cbind, lapply(pred_list, function(p) p$p.y.predictive.samples))

  # Local Bayesian means at each grid cell
  bayes_local_means   <- rowMeans(pred_all)

  # Overall Bayesian surface mean (average of local means)
  bayesian_surface_mean <- mean(bayes_local_means)

  cat("\nBayesian surface mean (predictive over grid):", bayesian_surface_mean, "\n")
}





##############################


# 1) Prediction coordinates (match CRS/units to training coords)
newcoords <- coordinates(grid_sp)  # matrix (n_pred x 2)

# 2) Posterior predictive samples at the grid for each chain
pred_list <- lapply(chains, function(fit)
  spPredict(fit, pred.coords = newcoords, start = burn, thin = thin, verbose = FALSE)
)

# 3) Combine prediction draws across chains
# Each is an n_pred x n_draws matrix (locations x posterior draws)
pred_all <- do.call(cbind, lapply(pred_list, function(p) p$p.y.predictive.samples))

# 4) Posterior predictive mean at each grid cell (local Bayesian means)
bayes_local_means <- rowMeans(pred_all)

# 5) Overall surface mean (average of local Bayesian means)
bayesian_surface_mean <- mean(bayes_local_means)

bayesian_surface_mean


####################
#traceplots

# Uses: mrec_list from spRecover() for 4 chains
params <- c("sigma.sq","tau.sq","phi")

# 1) Extract per-chain samples (same burn/thin) → mcmc objects
theta_list <- lapply(mrec_list, function(mrec) {
  as.mcmc(mrec$p.theta.recover.samples[, params, drop = FALSE])
})

# 2) Tidy long data frame with Iteration & Chain labels
df_long <- lapply(seq_along(theta_list), function(i) {
  as.data.frame(theta_list[[i]]) |>
    mutate(Iteration = seq_len(n())) |>
    mutate(Chain = paste0("Chain ", i)) |>
    pivot_longer(all_of(params), names_to = "Parameter", values_to = "Value")
}) |>
  bind_rows()

# 3) Single figure: 3 facets (one per parameter), all 4 chains overlaid
ggplot(df_long, aes(x = Iteration, y = Value, color = Chain, group = Chain)) +
  geom_line(alpha = 0.7, linewidth = 0.3) +
  facet_wrap(~ Parameter, scales = "free_y", ncol = 3) +
  labs(title = "Trace plots for sigma.sq, tau.sq, phi (4 chains)",
       x = "Kept iteration (after burn-in/thin)", y = "Value") +
  theme_minimal() +
  theme(legend.position = "bottom")
##########################
#posterior density plots

params <- c("sigma.sq","tau.sq","phi")
theta_list <- lapply(mrec_list, function(mrec) {
  as.mcmc(mrec$p.theta.recover.samples[, params, drop = FALSE])
})

# Long data frame of all chains
df <- lapply(seq_along(theta_list), function(i) {
  as.data.frame(theta_list[[i]]) |>
    mutate(Chain = paste0("Chain ", i))
}) |> bind_rows() |>
  pivot_longer(all_of(params), names_to = "Parameter", values_to = "Value")

# (Optional) add practical range panel
df_pr <- df |>
  filter(Parameter == "phi") |>
  transmute(Chain, Parameter = "practical_range (meters)", Value = 3/Value)

df_all <- bind_rows(df, df_pr)

# Posterior summaries (pooled across chains) for CI lines
summ <- df_all |>
  group_by(Parameter) |>
  summarize(
    q025 = quantile(Value, 0.025),
    med  = median(Value),
    q975 = quantile(Value, 0.975),
    .groups = "drop"
  )

ggplot(df_all, aes(x = Value, color = Chain)) +
  geom_density(adjust = 1.0, linewidth = 0.7, alpha = 0.9) +
  facet_wrap(~ Parameter, scales = "free", ncol = 2) +
  geom_vline(data = summ, aes(xintercept = med), linetype = 2) +
  geom_vline(data = summ, aes(xintercept = q025), linetype = 3) +
  geom_vline(data = summ, aes(xintercept = q975), linetype = 3) +
  labs(title = "Posterior densities by chain",
       x = "Value", y = "Density",
       subtitle = "Dashed = median, dotted = 95% CI (pooled across chains)") +
  theme_minimal() +
  theme(legend.position = "bottom")
#############################
#AUTOCORRELATION PLOTS


library(bayesplot)

# If you already have theta_list from earlier steps:
# theta_list: list of 4 mcmc objects with columns sigma.sq, tau.sq, phi
params <- c("sigma.sq","tau.sq","phi")
it <- nrow(as.matrix(theta_list[[1]]))

# 3D draws array: iterations × chains × parameters
arr <- array(NA_real_, dim = c(it, length(theta_list), length(params)),
             dimnames = list(iteration = 1:it,
                             chain = paste0("chain", 1:length(theta_list)),
                             parameter = params))
for (i in seq_along(theta_list)) {
  arr[, i, ] <- as.matrix(theta_list[[i]])[, params, drop = FALSE]
}

# Faceted ACF lines (default lags = 50; adjust lags arg as needed)
mcmc_acf(arr, lags = 50, facet_args = list(scales = "free_y")) +
  ggplot2::labs(title = "Autocorrelation (ACF) by parameter and chain")

# (Alternative) ACF bars per parameter (chains overlaid)
mcmc_acf_bar(arr, lags = 50) +
  ggplot2::labs(title = "ACF bar plots (overlaid by chain)")

################################
#EFFECTIVE SAMPLE SIZE

ess <- effectiveSize(theta_mcmc)

# Convert to data frame
ess_df <- data.frame(Parameter = names(ess),
                     ESS = as.numeric(ess))

# Desired order: sigma first, then tau, then phi
param_order <- c("sigma.sq","tau.sq","phi")

ess_df$Parameter <- factor(ess_df$Parameter, levels = param_order)

# Threshold (example: 400)
thr <- 400

# --- Bar plot ---
ggplot(ess_df, aes(x = Parameter, y = ESS, fill = Parameter)) +
  geom_col(show.legend = FALSE, alpha = 0.8) +
  geom_hline(yintercept = thr, linetype = "dashed", color = "red") +
  geom_text(aes(label = round(ESS, 0)), vjust = -0.5, size = 3.5) +
  labs(title = "Effective Sample Size (ESS) per parameter",
       x = "Parameter", y = "ESS") +
  theme_minimal()


###############################
#GEWEKE'S TEST
# parameters in desired order
params <- c("phi","tau.sq","sigma.sq")

theta_list <- lapply(mrec_list, function(mrec)
  as.mcmc(mrec$p.theta.recover.samples[, params, drop = FALSE])
)

# Geweke diagnostics
gz_list <- lapply(theta_list, function(chain)
  geweke.diag(chain, frac1 = 0.1, frac2 = 0.5)
)

# Build tidy data frame
gz_df <- bind_rows(lapply(seq_along(gz_list), function(i) {
  tibble(Chain = paste0("Chain ", i),
         Parameter = names(gz_list[[i]]$z),
         Z = as.numeric(gz_list[[i]]$z))
}))

# Chain order: descending
gz_df <- gz_df %>%
  mutate(Chain = factor(Chain, levels = paste0("Chain ", 4:1))) %>%
  mutate(Label = paste(Parameter, Chain, sep = " · "))

# Factor levels: phi first, then tau.sq, then sigma.sq
gz_df$Label <- factor(
  gz_df$Label,
  levels = unlist(lapply(params, function(p) paste(p, paste0("Chain ", 4:1), sep = " · ")))
)

# Plot
ggplot(gz_df, aes(x = Label, y = Z, color = Chain)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_hline(yintercept = c(-1.96, 1.96), linetype = 2, color = "gray40") +  # 95%
  geom_hline(yintercept = c(-2.58, 2.58), linetype = 3, color = "gray60") +  # 99%
  geom_point(size = 2.8, alpha = 0.9) +
  coord_flip() +
  labs(title = "Geweke Z-scores",
       x = NULL, y = "Geweke Z") +
  theme_minimal() +
  theme(legend.position = "bottom")


###########
#Heidelberger and Welch

# Apply Heidelberger–Welch test per chain
hw_results <- lapply(theta_list, function(chain) heidel.diag(chain))

# Print results
hw_results

###########
#Raftery Lewis Test

# Apply Raftery–Lewis per chain
rl_results <- lapply(theta_list, function(chain) raftery.diag(chain))

# Print results
rl_results
###########
library(coda)

# theta_list = list of your 4 mcmc objects (phi, tau.sq, sigma.sq)
theta_mcmc <- mcmc.list(theta_list)

# Gelman–Rubin diagnostic
gelman_results <- gelman.diag(theta_mcmc)

# Print results
gelman_results

