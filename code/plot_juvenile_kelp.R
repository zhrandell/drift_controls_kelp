## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## visualize juvenile giant kelp from SNI subtidal timeseries
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



## load ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dat <- read.csv(file.path(data, "juvenile_kelp_SNI_subtidal_timeseries.csv"))
## END load ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## format data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## create new df with juvenile and adult macro counts 
df <- data.frame(
  MacJuv = as.numeric(dat$MacJuv),
  MacPyr = as.numeric(dat$MacPyr)
)

## drop all 0's from MacJuv, as any associated MacPyr rows
df <- df %>% filter(MacJuv != 0)

## log10 transform 
df$MacJuv_log10 <- log10(df$MacJuv)


## END data prep ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## graphing params ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
my.theme = theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line = element_line(colour = "black"),
                 axis.title=element_text(size=20),
                 axis.text=element_text(size=18),
                 plot.title = element_text(size=20))
## END graphing params ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## function to create ECDF with shaded band between two log10() x-values
plot_ecdf <- function(dat, 
                      col,
                      breaks_by = 1,
                      shade_min = log10(20),
                      shade_max = log10(80)) {
  
  vals <- dplyr::pull(dat, {{ col }})
  vals <- vals[!is.na(vals)]
  n_used <- length(vals)
  
  ecdf_fun <- ecdf(vals)
  xs <- sort(unique(c(vals, shade_min, shade_max)))
  ecdf_df  <- tibble::tibble(x = xs, y = ecdf_fun(xs))
  shade_df <- dplyr::filter(ecdf_df, x >= shade_min & x <= shade_max)
  
  xmax <- max(vals, na.rm = TRUE)
  
  y_min <- ecdf_fun(shade_min)
  y_max <- ecdf_fun(shade_max)
  perc_shaded <- 100 * (y_max - y_min)
  
  p <- ggplot() +
    geom_area(data = shade_df, aes(x = x, y = y),
              fill = "steelblue", alpha = 0.2) +
    stat_ecdf(data = tibble::tibble(v = vals), aes(x = v),
              geom = "step", linewidth = 0.75, na.rm = TRUE) +
    labs(
      title = '',
      x     = expression("Juvenile Giant Kelp per 1 m"^2),
      y     = "Cumulative proportion"
    ) +
    scale_x_continuous(
      breaks = seq(-4, xmax, by = breaks_by),
      labels = function(x){ (10^x) / 20 }
    ) +
    scale_y_continuous(limits = c(0, 1), 
                       expand = c(0, 0),
                       breaks = seq(0, 1, by = 0.1)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    my.theme +
    annotate( 'text',
              x = (shade_min + shade_max) / 2, 
              y = 0.3,
              label = 'Experimental\nrange')
  
  ## text output
  ## raw natural-scale bounds
  nat_min <- 10^shade_min
  nat_max <- 10^shade_max
  
  ## rescale to 1m^2 (divide by 20)
  nat_min_1m2 <- nat_min / 20
  nat_max_1m2 <- nat_max / 20
  
  ## log10-scale bounds rescaled to 1m^2 (divide by 20) 
  log_min_1m2 <- shade_min / 20
  log_max_1m2 <- shade_max / 20
  
  text_out <- sprintf(
    paste0(
      "Rows used: %d\n",
      "Shaded band (log10 scale) at the 20 m^2 spatial scale: [%.3f, %.3f]\n",
      "Shaded band (log10 scale) at the 1 m^2 spatial scale: [%.3f, %.3f]\n",
      "Shaded band natural scale (10^x) at the 20 m^2 spatial scale: [%.3f, %.3f]\n",
      "Shaded band natural scale (10^x) at the 1 m^2 spatial scale: [%.3f, %.3f]\n",
      "Cumulative percentage at band minimum: %.2f%%\n",
      "Cumulative percentage at band maximum: %.2f%%\n",
      "Percentage of ECDF shaded: %.2f%%\n"
    ),
    n_used,
    shade_min, shade_max,
    log_min_1m2, log_max_1m2,
    nat_min, nat_max,
    nat_min_1m2, nat_max_1m2,
    100 * y_min,
    100 * y_max,
    perc_shaded
  )
  
  cat(text_out)
  writeLines(text_out, file.path(results, "juvenile_kelp.txt"))
  
  return(p)
}
## END function to create ecdf ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## view plot and save ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p1 <- plot_ecdf(df, 
                    MacJuv_log10,
                    breaks_by = 1,
                    shade_min = log10(20),
                    shade_max = log10(80))
print(p1)


## save pdf 
ggsave(
  filename = paste0(figs, "/juvenile_kelp.pdf"),
  plot = p1,
  device = cairo_pdf,
  dpi = 1200,
  width = 10,
  height = 7,
  units = "in"
)
## END view and save ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
