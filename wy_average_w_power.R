library(tidyverse)
library(lubridate)
library(zoo)
library(WaveletComp)

# Load the COVID data set
wy_covid <- read_csv('D:\WY-COVID Data\Gaussian Smoothing\WY Gaussian Data\WY_gaussian_smoothened_data.csv', show_col_types = FALSE)

# Output Directory
out_dir <- "D:\WY-COVID Data\Gaussian Smoothing\WY Gaussian Results\Avg Wavelet Power Plot"

# Extract all counties
counties <- unique(wy_covid$County)

# Loop through each county
for(cty in counties) {
  message("Processing county: ", cty)
  
  # Process data for the county
  county_data <- wy_covid %>%
    filter(County == cty) %>%
    mutate(Date = ymd(Date)) %>%
    arrange(Date) %>%
    select(Date, cases = gauss_smoothed_per_10k) %>%
    filter(!is.na(cases)) %>%
    mutate(standardized = as.numeric(scale(cases)))
  
  # Prepare for Wavelet Comp
  county_wave <- county_data %>%
    select(Date, standardized) %>% 
    rename(date = Date, cases = standardized) %>%
    mutate(date = as.POSIXct(date, tz = "UTC")) %>%
    as.data.frame()
  
  # Run the wavelet analysis
  county_wt <- analyze.wavelet(
    my.data     = county_wave,
    my.series   = "cases",
    loess.span  = 0.75,
    dt          = 1,
    dj          = 1/20,
    lowerPeriod = 2,
    upperPeriod = floor(nrow(county_wave) / 3),  
    make.pval   = TRUE,
    method      = "white.noise",
    n.sim       = 1000
  )
  
  # Compute the max level for each county
  maximum_power = 1.001 * max(county_wt$Power.avg, na.rm = TRUE)
  
  # Create County-Specific File Name
  filename <- file.path(out_dir, paste0(cty, "_wt_avg.png"))
  
  # Save the average-power plot
  png(
    filename = filename,
    width    = 8,
    height   = 7,
    units    = "in",
    res      = 1600
  )
  
  wt.avg(
    county_wt,
    maximum.level = maximum_power,
    siglvl = 0.05,
    sigcol = "red",
    show.legend = TRUE,
    main = paste("Average Wavelet Power Spectrum:", cty)
  )
  
  
  dev.off()
  message("Saved plot for ", cty, " to: ", filename)
}

message("\nProcessing complete for all counties")

