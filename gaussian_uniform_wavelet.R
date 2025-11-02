library(WaveletComp)
library(tidyverse)
library(lubridate)
library(zoo)

input_file <- "D:/WY-COVID Data/Gaussian Smoothing/WY Gaussian Data/WY_gaussian_smoothened_data.csv"
out_dir    <- "D:/WY-COVID Data/Gaussian Smoothing/WY Gaussian Results/Uniform Wavelet Power Spectrum"

# Read data & grab every county
covid_data <- read_csv(input_file, show_col_types = FALSE)
counties   <- unique(covid_data$County)
message("Will process all counties: ", paste(counties, collapse = ", "))

# Compute wavelets & track global max
wavelet_list <- vector("list", length(counties))
names(wavelet_list) <- counties
global_max <- 0

for(cty in counties) {
  message("Computing wavelet for county: ", cty)
  cd <- covid_data %>%
    filter(County == cty) %>%
    select(Date, cases = gauss_smoothed_per_10k) %>%
    filter(!is.na(cases)) %>%
    arrange(Date) %>%
    mutate(Date = as.Date(Date))
  
  wave_df <- cd %>%
    select(Date, cases) %>%
    rename(date = Date) %>%
    mutate(date = as.POSIXct(date, tz = "UTC")) %>%
    as.data.frame()
  
  wt <- analyze.wavelet(
    my.data     = wave_df,
    my.series   = "cases",
    loess.span  = 0.75,
    dt          = 1,
    dj          = 1/20,
    lowerPeriod = 2,
    upperPeriod = 256,        
    make.pval   = TRUE,
    method      = "white.noise",
    n.sim       = 1000
  )
  
  wavelet_list[[cty]] <- list(wt = wt, dates = as.Date(wave_df$date))
  global_max <- max(global_max, max(wt$Power, na.rm = TRUE))
}

message("Global maximum wavelet power = ", round(global_max, 2))

# ggplot2 plotting function
plot_wavelet_gg <- function(wt_obj, county_name, out_dir, max_power) {
  wt     <- wt_obj$wt
  dates  <- wt_obj$dates
  period <- wt$Period
  power  <- wt$Power
  pval   <- wt$Power.pval
  
  # cap above the global max
  power[power > max_power] <- max_power
  
  # build a long data.frame
  df <- expand.grid(Period = period, Time = dates) %>%
    mutate(
      Power = as.vector(power),
      Sig   = as.vector(pval < 0.05)
    )
  
  pal <- colorRampPalette(c("blue","cyan","green","yellow","orange","red"))(100)
  
  p <- ggplot(df, aes(x = Time, y = Period, fill = Power)) +
    geom_tile() +
    geom_contour(aes(z = as.numeric(Sig)), breaks = 0.5,
                 color = "white", size = 0.3) +
    scale_x_date(
      date_breaks = "1 year",
      date_labels = "%Y",
      expand      = c(0,0)
    ) +
    scale_y_log10(
      breaks = c(2,4,8,16,32,64,128),
      labels = c(2,4,8,16,32,64,128),
      expand = c(0,0)
    ) +
    scale_fill_gradientn(
      colours = pal,
      limits  = c(0, max_power),
      name    = "Power"
    ) +
    labs(
      title = paste("Morlet Wavelet Power Spectrum —", county_name),
      x     = "Date",
      y     = "Period (days)"
    ) +
    theme_dark(base_size = 14) +
    theme(
      plot.title       = element_text(hjust = 0.5),
      panel.grid.minor = element_blank()
    )
  
  # write high-res 8×7 in PNG
  ggsave(
    filename = file.path(out_dir, paste0(county_name, "_wavelet.png")),
    plot     = p,
    width    = 8,
    height   = 7,
    units    = "in",
    dpi      = 1600
  )
}

# Render all plots via ggplot2
for(cty in counties) {
  message("Plotting county: ", cty)
  plot_wavelet_gg(wavelet_list[[cty]], cty, out_dir, global_max)
}

message("All plots saved to ", out_dir)
