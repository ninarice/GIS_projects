library(sf)       # For handling spatial data
library(dplyr)    # For data manipulation
library(ggplot2)  # For visualization
library(tidyverse) # General data wrangling

# Load CSV and shapefile
tabular_data <- read.csv("~/Downloads/Assignment4Data/mex_data.csv")
shapefile <- st_read("~/Downloads/Assignment4Data/municipalities.shp")

# Rename columns in tabular_data to match shapefile
tabular_data <- tabular_data %>%
  rename(CVE_ENT = state, CVE_MUN = mun)

# Calculate homicide rate per 100,000 and log-transformed homicide rate
tabular_data <- tabular_data %>%
  mutate(
    homicide_rate = (homicides / population) * 100000,
    log_homicide_rate = ifelse(homicides == 0, NA, log(homicide_rate))  # Assign NA for zero homicides
  )

# Perform regression on log homicide rate and poverty index, excluding NAs
ols_model <- lm(log_homicide_rate ~ indice, data = tabular_data, na.action = na.exclude)

# Print regression summary
print(summary(ols_model))

# Calculate the effect of one standard deviation of the poverty index
sd_indice <- sd(tabular_data$indice, na.rm = TRUE)
effect_log_rate <- coef(ols_model)[2] * sd_indice
effect_rate <- exp(effect_log_rate)
print(paste("Effect on log homicide rate (1 SD poverty):", effect_log_rate))
print(paste("Effect on homicide rate multiplier:", effect_rate))

# Add residuals to tabular_data
tabular_data <- tabular_data %>%
  mutate(
    residuals = ifelse(is.na(log_homicide_rate), NA, residuals(ols_model))  # NA residuals for zero homicides
  )

# Ensure columns in both datasets match and add leading zeros to municipality/state codes
shapefile <- shapefile %>%
  mutate(CVE_ENT = as.character(CVE_ENT), CVE_MUN = as.character(CVE_MUN))

tabular_data <- tabular_data %>%
  mutate(
    CVE_ENT = str_pad(as.character(CVE_ENT), width = 2, pad = "0"),
    CVE_MUN = str_pad(as.character(CVE_MUN), width = 3, pad = "0")
  )

# Perform a left_join to retain all shapefile geometries
merged_data <- shapefile %>%
  left_join(tabular_data, by = c("CVE_ENT", "CVE_MUN")) %>%
  st_as_sf()

# Update data_status to classify rows
merged_data <- merged_data %>%
  mutate(
    data_status = case_when(
      is.na(residuals) ~ "No Data or Zero Homicides",  # Missing data or zero homicides
      TRUE ~ "Residuals Available"                    # Residuals calculated
    )
  )

# Check residual values
min_residual <- min(merged_data$residuals, na.rm = TRUE)
max_residual <- max(merged_data$residuals, na.rm = TRUE)
print(paste("Min residual:", min_residual))
print(paste("Max residual:", max_residual))

# Print data rows for debugging
print(paste("Shapefile rows:", nrow(shapefile)))
print(paste("Tabular data rows:", nrow(tabular_data)))
print(paste("Merged data rows:", nrow(merged_data)))

# Create the plot
ols_plot <- ggplot() +
  # Residuals layer
  geom_sf(
    data = merged_data %>% filter(data_status == "Residuals Available"),
    aes(fill = residuals),
    color = "gray90",
    size = 0.05
  ) +
  # Adjusted color scale for better differentiation
  scale_fill_distiller(
    palette = "RdBu",
    limits = c(-4, ceiling(max_residual)),
    breaks = seq(-4, ceiling(max_residual), by = 1),
    name = "OLS Residuals"
  ) +
  # No Data or Zero Homicides layer
  geom_sf(
    data = merged_data %>% filter(data_status == "No Data or Zero Homicides"),
    fill = "gray",
    color = "gray90",
    size = 0.05
  ) +
  # Titles and captions
  labs(
    title = "OLS Residuals Map for Homicide Rates in Mexican Municipalities",
    subtitle = "Residuals from the regression of log homicide rate on poverty index",
    caption = "Gray areas: No data/Zero homicides"
  ) +
  # Theme adjustments
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  )

# Save the plot
ggsave("OLS_Residuals_Map_Final100.png", plot = ols_plot, width = 10, height = 8, dpi = 300)




