# Load necessary libraries
library(sf)
library(dplyr)
library(readr)
library(ggplot2)
library(spgwr)
library(sp)
library(spdep)
library(spatialEco)
library(spatstat)

# Load shapefiles
retailers_per_muni <- st_read("~/Documents/GPEC_443/Final/ecomm.shp")
census_sf <- st_read("~/Downloads/calenviroscreen40shpf2021shp/CES4 Final Shapefile.shp")
sdi <- read_csv("~/Documents/GPEC_443/Final/SDI_censustract.csv")
muni_geoms <- st_read("~/Documents/GPEC_443/Final/sd_munis.shp")
population <- read_csv("~/Documents/GPEC_443/Final/Updated_Population_Data.csv")
by_retailer <- st_read("~/Documents/GPEC_443/Final/by_retailer.shp")

# =============================================================================
# Data Preparation, Aggregated by Muni and Aggregated by Retailer (by_retailer)
# =============================================================================

# Clean and prepare retailers_per_muni
retailers_per_muni <- retailers_per_muni %>%
  mutate(name = ifelse(name == "S.D. COUNTY", "Unincorporated Areas", name)) %>%
  mutate(name = tools::toTitleCase(tolower(name))) %>%
  st_drop_geometry()

# Filter census data for San Diego County
census_sf_sd <- census_sf %>%
  filter(County == "San Diego") %>%
  st_drop_geometry()  # Drop geometry as it's not needed

# Prepare SDI data
sdi <- sdi %>%
  rename(Tract = CENSUSTRACT_FIPS)  # Rename for joining

# Merge census data with SDI data
sdi_tracts <- census_sf_sd %>%
  left_join(sdi, by = "Tract") %>%
  dplyr::select(ApproxLoc, SDI_score) %>%
  mutate(SDI_score = as.numeric(SDI_score))

# Aggregate SDI scores by municipality
sdi_aggregated <- sdi_tracts %>%
  group_by(ApproxLoc) %>%
  summarize(average_sdi_score = mean(SDI_score, na.rm = TRUE))

# Match SDI scores to retailers_per_muni
retailers_per_muni <- retailers_per_muni %>%
  left_join(sdi_aggregated, by = c("name" = "ApproxLoc"))

# Handle unmatched rows by adding them to Unincorporated Areas
unmatched_rows <- sdi_aggregated %>%
  filter(!ApproxLoc %in% retailers_per_muni$name)

unincorporated_row <- retailers_per_muni %>%
  filter(name == "Unincorporated Areas") %>%
  mutate(average_sdi_score = mean(unmatched_rows$average_sdi_score, na.rm = TRUE))

retailers_per_muni <- retailers_per_muni %>%
  filter(name != "Unincorporated Areas") %>%
  bind_rows(unincorporated_row)

# Add population data to retailers_per_muni
retailers_per_muni <- retailers_per_muni %>%
  left_join(population, by = c("name" = "Municipality"))

# Calculate retailer rate (retailers per 100,000 population)
retailers_per_muni <- retailers_per_muni %>%
  mutate(retailers_per_100k = (retailers / Population) * 100000)

# Log-transform retailer rate
retailers_per_muni <- retailers_per_muni %>%
  mutate(log_retailers_per_100k = ifelse(
    retailers_per_100k == 0, NA, log(retailers_per_100k)
  ))

# Calculate the SDI indice (normalized SDI score)
retailers_per_muni <- retailers_per_muni %>%
  mutate(
    sdi_indice = (average_sdi_score - mean(average_sdi_score, na.rm = TRUE)) / 
      sd(average_sdi_score, na.rm = TRUE)
  )

# Aggregate geometries by municipality
aggregated_geoms <- muni_geoms %>%
  group_by(name) %>%
  summarize(geometry = st_union(geometry))

# Clean `name` column in `retailers_per_muni`
retailers_per_muni <- retailers_per_muni %>%
  mutate(name = toupper(name)) %>%
  mutate(name = ifelse(name == "UNINCORPORATED AREAS", "S.D. COUNTY", name))

# Perform the left join
muni_data <- aggregated_geoms %>%
  left_join(retailers_per_muni, by = "name")

# Filter for complete cases
analyze <- muni_data %>%
  filter(!is.na(log_retailers_per_100k) & !is.na(sdi_indice))

# Add the log_retailers column to analyze
analyze <- analyze %>%
  mutate(log_retailers = ifelse(retailers == 0, NA, log(retailers)))

# Merge on "name"
analyze_ret <- by_retailer %>%
  left_join(aggregated_geoms %>% st_drop_geometry(), by = "name")

# Add the polygon geometry from aggregated_geoms to the merged table
analyze_ret <- analyze_ret %>%
  mutate(geometry = st_geometry(aggregated_geoms)[match(name, aggregated_geoms$name)])

# Ensure the resulting table is an sf object
analyze_ret <- st_as_sf(analyze_ret)

muni_data <- muni_data %>%
  mutate(retailers_z_score = (retailers - mean(retailers, na.rm = TRUE)) / 
           sd(retailers, na.rm = TRUE))

# Ensure both are in the same CRS (Coordinate Reference System)
muni_data <- st_transform(muni_data, crs = st_crs(by_retailer))

# Perform the spatial join
final_table <- st_join(by_retailer, muni_data, left = TRUE)

# Select the desired columns
final_table <- final_table %>%
  select(lat, lon, name.x, FTR.x, average_sdi_score, Population, sdi_indice, retailers_z_score)

# saving a sf
#st_write(muni_data, "~/Documents/GPEC_443/Final/muni_sdi.shp")

# =========================
# Welch's Two Sample T-test
# =========================

t_test_z <- t.test(retailers_per_100k ~ FTR, data = muni_data)
print(t_test_z)

# =================
# Nearest Neighbors
# =================

st_crs(by_retailer)

# Extract coordinates
coords <- st_coordinates(final_table)

# Define the bounding window
window <- owin(
  xrange = range(coords[, 1]),
  yrange = range(coords[, 2])
)

# Convert to ppp object
retailer_ppp <- as.ppp(coords, W = window)

# Nearest neighbor distances
nn_dist <- nndist(retailer_ppp)
summary(nn_dist)

# Calculate observed mean nearest neighbor distance
observed_mean_distance <- mean(nn_dist)

# Calculate expected mean nearest neighbor distance for a Poisson process
area <- area.owin(retailer_ppp$window)
expected_mean_distance <- 1 / (2 * sqrt(npoints(retailer_ppp) / area))

# Compute Nearest Neighbor Index (NNI)
nni <- observed_mean_distance / expected_mean_distance

# Interpret the NNI
if (nni < 1) {
  interpretation <- "Clustered pattern"
} else if (nni == 1) {
  interpretation <- "Random distribution"
} else {
  interpretation <- "Regular pattern"
}

# Print results
cat("Observed Mean Distance:", observed_mean_distance, "\n")
cat("Expected Mean Distance:", expected_mean_distance, "\n")
cat("Nearest Neighbor Index (NNI):", nni, "\n")
cat("Interpretation:", interpretation, "\n")

# ====================
# OLS Regression, Muni 
# ====================

# Perform the OLS regression
ols_model <- lm(retailers ~ FTR + sdi_indice + Population, data = analyze)

# Print the summary of the regression model
summary(ols_model)

# Calculate residuals from ols_model2 and add them to a new data frame
residuals_df <- analyze %>%
  mutate(ols_residuals = residuals(ols_model)) %>%  # Reference the correct model
  select(name, ols_residuals) %>%  # Keep only name and residuals for joining
  st_drop_geometry() %>%           # Drop geometry to avoid conflicts
  as.data.frame()                  # Convert to a standard data frame

# Merge the residuals back into muni_data using st_join
muni_data <- muni_data %>%
  left_join(residuals_df, by = "name")

# Create the map
ols_plot <- ggplot() +
  geom_sf(
    data = muni_data %>% filter(!is.na(ols_residuals)),
    aes(fill = ols_residuals),
    color = "grey44",
    size = 0.0005
  ) +
  scale_fill_distiller(
    palette = "RdBu",
    limits = c(-1.8, 1.8),
    breaks = seq(-1.8, 1.8, by = 0.9),
    name = "OLS Residuals\n(Retailer Presence vs\nFTR, SDI, & Population)"
  ) +
  geom_sf(
    data = muni_data %>% filter(is.na(ols_residuals)),
    fill = "grey",
    color = "gray44",
    size = 0.0005
  ) +
  labs(
    title = "OLS Residuals of FTR, SDI, and Population on Brick-and-click Presence Across San Diego County",
    subtitle = "Effects of flavored tobacco restrictions, social deprivation, and population size on e-commerce tobacco retailer counts.",
    caption = "Grey: No brick-and-click tobacco retailers present.\nDate: December 2024"
  ) +
  annotation_north_arrow(
    location = "tr", 
    which_north = "true", 
    style = north_arrow_fancy_orienteering()  
  ) +
  annotation_scale(
    location = "br",  
    width_hint = 0.25 
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.12, 0.18),
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
    legend.key.size = unit(0.8, "lines"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 2)
  )

# Print and save the plot
print(ols_plot)
ggsave("~/Documents/GPEC_443/Final/ols_residuals_map_2.png", plot = ols_plot, width = 10, height = 8, dpi = 300)

# ===============
# Kernal density
# ===============

# Reproject spatial data to UTM Zone 11N
final_table <- st_transform(final_table, crs = 32611)

# Reproject muni_data to UTM Zone 11N
muni_data <- st_transform(muni_data, crs = 32611)

# Create unique points dataset by removing duplicates
final_table_unique <- final_table %>%
  distinct(geometry, .keep_all = TRUE)

# Step 1: Define the Extent from muni_data
# Get bounding box from muni_data
muni_bbox <- st_bbox(muni_data)

# Create a window from the bounding box
window <- owin(
  xrange = c(muni_bbox["xmin"], muni_bbox["xmax"]),
  yrange = c(muni_bbox["ymin"], muni_bbox["ymax"])
)

# Convert retailer points to ppp format
coords <- st_coordinates(final_table_unique)  # Use your unique points dataset
ppp_data <- ppp(x = coords[, 1], y = coords[, 2], window = window)

# Perform KDE
kde <- density(ppp_data, sigma = 1421, edge = TRUE)  # Adjust sigma as needed

# Convert KDE to raster
kde_raster <- raster(kde)

# Update CRS of the raster to match muni_data
crs(kde_raster) <- st_crs(muni_data)$proj4string

# Convert density to per 10 square kilometers
kde_raster_10km2 <- kde_raster * 1e7

# Debugging Step: Plot KDE raster before masking
kde_raster_df <- as.data.frame(kde_raster_10km2, xy = TRUE)
kde_raster_df <- kde_raster_df %>% na.omit()
colnames(kde_raster_df) <- c("x", "y", "density")

initial_plot <- ggplot() +
  geom_raster(data = kde_raster_df, aes(x = x, y = y, fill = density)) +
  geom_sf(data = muni_data, fill = NA, color = "black", size = 0.5) +
  scale_fill_viridis_c(option = "D", name = "Density (per 10 km²)", direction = -1) +
  labs(
    title = "Kernel Density of Retailers (Before Masking)",
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal()

print(initial_plot)

# Mask Raster to Municipal Boundaries
muni_sp <- as(st_geometry(muni_data), "Spatial")  # Convert sf to Spatial object
kde_raster_masked <- mask(kde_raster_10km2, muni_sp)

# Convert masked raster to a data frame for ggplot2
kde_df_masked <- as.data.frame(kde_raster_masked, xy = TRUE)
kde_df_masked <- kde_df_masked %>% na.omit()
colnames(kde_df_masked) <- c("x", "y", "density")

# Plot the KDE with Municipal Boundaries
kde_plot <- ggplot() +
  # Raster layer (KDE)
  geom_raster(data = kde_df_masked, aes(x = x, y = y, fill = density)) +
  # Municipal boundaries
  geom_sf(data = muni_data, fill = NA, color = "black", size = 0.5) +
  # Adjust color scale
  scale_fill_viridis_c(option = "D", name = "Density (per 10 km²)", direction = -1) +
  # Add labels and theme
  labs(
    title = "Kernel Density of Brick-and-click Retailers (Masked to Municipal Boundaries)",
    x = "Longitude", y = "Latitude"
  ) +
  theme_minimal()

print(kde_plot)
# Save the plot
ggsave("kde_density_plot.png", plot = kde_plot, width = 10, height = 8, dpi = 300)

# =========================================================================
# Global Moran's I, Muni, test if retailer rates exhibit spatial clustering 
# =========================================================================

# Create a new spatial weights matrix based on the filtered dataset
nb_analyze <- poly2nb(analyze, queen = TRUE)
listw_non_na <- nb2listw(nb_analyze, style = "W")

# Moran's I test
moran_test <- moran.test(analyze$retailers_per_100k, listw_non_na)

# Print the result
print(moran_test)

# Visualize spatial autocorrelation
library(spatialEco)
moran.plot(
  analyze$retailers_per_100k, 
  listw_non_na, 
  labels = analyze$name, 
  xlab = "Retailer Rates",
  ylab = "Spatially Lagged Rates",
  main = "Moran's I Plot for Retailer Rates"
)

