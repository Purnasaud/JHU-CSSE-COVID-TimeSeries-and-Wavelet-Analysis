library(WaveletComp)
library(tidyverse)
library(lubridate)
library(zoo)

# Load the COVID dataset
covid_data <- read_csv('D:/WY-COVID Data/Gaussian Smoothing/WY Gaussian Data/WY_gaussian_smoothened_data.csv', show_col_types = FALSE)

# Output Directory
out_dir <- "D:/WY-COVID Data/Gaussian Smoothing/WY Gaussian Results"

# Extract all counties
counties <- unique(covid_data$County)

wavelet_list <- list()

# Compute wavelet for each county
for(cty in counties) {
  message("Processing county: ", cty)
  
  # Process data
  county_data <- covid_data %>%
    filter(County == cty) %>%
    select(Date, cases = gauss_smoothed_per_10k) %>%   
    filter(!is.na(cases)) %>%       
    arrange(Date) %>%
    mutate(standardized = as.numeric(scale(cases)))
  
  # Prepare for wavelet
  counties_wave <- county_data %>%
    select(Date, standardized) %>% 
    rename(date = Date, cases = standardized) %>%
    mutate(date = as.POSIXct(date, tz = "UTC")) %>%
    as.data.frame()
  
  # Compute wavelet transform
  counties_wt <- analyze.wavelet(
    my.data     = counties_wave,
    my.series   = "cases",
    loess.span  = 0.75,
    dt          = 1,
    dj          = 1/20,
    lowerPeriod = 2,
    upperPeriod = floor(nrow(counties_wave) / 3), 
    make.pval   = TRUE,
    method      = "white.noise",
    n.sim       = 1000
  )
  
  wavelet_list[[cty]] <- counties_wt
}

# Generate plots for all counties
for(cty in counties) {
  message("Plotting county: ", cty)
  counties_wt <- wavelet_list[[cty]]
  
  # Create county-specific filename
  filename <- file.path(out_dir, paste0(cty, "_wavelet_power.png"))
  
  png(
    filename = filename,
    width    = 8,
    height   = 7,
    units    = "in",
    res      = 1600
  )
  
  # Plot with each county's own max scaling
  wt.image(
    counties_wt,
    my.series = 1,
    plot.coi = TRUE,
    plot.contour = TRUE,
    col.contour = "white",
    plot.ridge = TRUE,
    col.ridge = "black",
    color.key = "quantile",
    n.levels = 100,
    color.palette = "rainbow(100, start=0, end=.7)",
    useRaster = TRUE,
    max.contour.segments = 250000,
    plot.legend = TRUE,
    legend.params = list(
      width = 1.2,
      shrink = 0.9,
      mar = 5.1,
      n.ticks = 6,
      label.digits = 1,
      label.format = "f",
      lab = NULL,
      lab.line = 2.5
    ),
    label.time.axis = TRUE,
    show.date = TRUE,
    date.format = "%Y-%m",
    timelab = "Time",
    label.period.axis = TRUE,
    periodlab = "Period (days)",
    main = paste("Morlet Wavelet Power Spectrum:", cty),
    lwd = 2,
    graphics.reset = TRUE,
    siglvl = 0.05
  )
  
  dev.off()
  message("Saved plot for ", cty, " to: ", filename)
}

message("\nProcessing complete for all counties")
