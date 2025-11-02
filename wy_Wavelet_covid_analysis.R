library(tidyverse)
library(lubridate)
library(zoo)           
library(WaveletComp)

# Read & prepare Wyoming daily series
wy_covid <- read_csv("D:/WY-COVID Data/Wyoming_Covid-Analysis/complete_state_wide_covid_cases.csv")

wy_daily <- wy_covid %>% 
  filter(Province_State == "Wyoming") %>% 
  mutate(Date = ymd(Date)) %>% 
  arrange(Date) %>%
  group_by(Date) %>%
  summarize(cum_cases = sum(Cases, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(Date) %>%
  mutate(
    daily_new  = cum_cases - lag(cum_cases, default = first(cum_cases)),
    daily_new  = if_else(daily_new < 0, 0, daily_new),
    daily_norm = daily_new / 5850935 * 10000,
    wy_ma3        = rollmean(daily_norm, k = 3, fill = NA, align = "right")
  )

# Build time index, detrend & standardize on the smoothed series
wy_series <- wy_daily %>%
  select(Date, cases = wy_ma3) %>%   
  filter(!is.na(cases)) %>%       
  arrange(Date) %>%
  mutate(
    t = as.numeric(Date - min(Date)) + 1
  ) %>%
  { 
    fit <- lm(cases ~ poly(t, 2), data = ., na.action = na.exclude)
    mutate(
      ., 
      detrended    = resid(fit),
      standardized = as.numeric(scale(detrended))
    )
  }

# Prepare for WaveletComp
wy_wave <- wy_series %>%
  select(Date, standardized) %>% 
  rename(
         date = Date,
         cases = standardized
         ) %>% 
  select(date, cases) %>%
  mutate(date = as.POSIXct(date, tz = "UTC")) %>%
  as.data.frame()

# Run the wavelet analysis
wy_wt <- analyze.wavelet(
  my.data     = wy_wave,
  my.series   = "cases",
  loess.span  = 0.75,
  dt          = 1,
  dj          = 1/20,
  lowerPeriod = 2,
  upperPeriod = floor(nrow(wy_wave) / 3),  
  make.pval   = TRUE,
  method      = "white.noise",
  n.sim       = 1000
)


# Compute the max level
maximum_power = 1.001*max(wy_wt$Power.avg, wy_wt$Power.avg)

# Plot and save high-res image
out_dir <- "D:/WY-COVID Data/Wyoming_Covid-Analysis/R_wavelet_plots/Wyoming"
png(
  # filename = file.path(out_dir, "Wyoming_wavelet_quantile_new_1.png"),
  filename = file.path(out_dir, "Wyoming_wt_avg_1.png"),
  width    = 8,
  height   = 7,
  units = "in",
  res = 1600
)
# wt.image(
#   wy_wt,
#   legend.params = list(lab = "wavelet power levels", label.digits = 2),
#   color.key     = "quantile",
#   periodlab     = "Period (days)",
#   show.date     = TRUE,
#   date.format   = "%Y-%m-%d",
#   date.tz       = "UTC",
#   siglvl        = 0.05,
#   n.levels      = 250,
#   main          = "Morlet Wavelet Power Spectrum: Wyoming Daily Cases"
# )

# Draw the average-power plot
wt.avg(wy_wt,
       maximum.level = maximum_power,
       siglvl = 0.05,
       sigcol = "red",
       show.legend = TRUE,
       main = "Average Wavelet Power Spectrum of Wyoming Daily COVID-19 Cases" 
       )

dev.off()
# utils::browseURL("Wyoming_wavelet_quantile_new_1.png")
utils::browseURL("Wyoming_wt_avg_1.png")

