## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## visualize juvenile giant kelp from SNI subtidal timeseries
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## load ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list=ls())

library(readr)
library(tidyverse)

getwd()
setwd("../")

data <- "data"
results <- "results"
figs <- "figs"

dat <- read.csv(file.path(data, "juvenile_kelp_SNI_subtidal_timeseries.csv"))
## END load ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## format data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## create new df with juvenile and adult macro counts 
df <- data.frame(
  MacJuv = as.numeric(dat$MacJuv),
  MacPyr = as.numeric(dat$MacPyr)
)


## divide by 20 to provide counts per 1m^2 
df <- df / 20 


## log10 transform 
df$MacJuv_log10 <- round(log10(df$MacJuv + 1), 3)
df$MacPyr_log10 <- round(log10(df$MacPyr + 1), 3)


## drop all 0's from MacJuv, as any associated MacPyr rows
df <- df %>% filter(MacJuv_log10 != 0)
## END data prep ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## graphing params ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
my.theme = theme(panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line = element_line(colour = "black"),
                 axis.title=element_text(size=20),
                 axis.text=element_text(size=18),
                 plot.title = element_text(size=20))


graphics.off()
windows(10,10,record=TRUE)
## END graphing params ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## function to create ecdf figure ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_ecdf <- function(data, 
                      col,
                      breaks_by = 0.2,
                      shade_to  = log10(5)) {
  
  vals <- dplyr::pull(data, {{ col }})
  vals <- vals[!is.na(vals)]
  n_used <- length(vals)
  
  ecdf_fun <- ecdf(vals)
  xs <- sort(unique(vals))
  ecdf_df <- tibble::tibble(x = xs, y = ecdf_fun(xs))
  shade_df <- dplyr::filter(ecdf_df, x <= shade_to)
  
  y_at_cut     <- ecdf_fun(shade_to)
  x_raw_at_cut <- 10^shade_to - 1
  xmax <- max(vals, na.rm = TRUE)
  
  p <- ggplot() +
    geom_area(data = shade_df, aes(x = x, y = y),
              fill = "steelblue", alpha = 0.2) +
    stat_ecdf(data = tibble::tibble(v = vals), aes(x = v),
              geom = "step", linewidth = 0.75, na.rm = TRUE) +
    geom_hline(yintercept = y_at_cut, linetype = "dashed", linewidth = 1.25) +
    labs(
      title = expression("Juvenile Giant Kelp per 1m"^2),
      x     = expression("log"[10] * "(x+1) juvenile Giant Kelp"),
      y     = "Cumulative proportion"
    ) +
    scale_x_continuous(breaks = seq(0, xmax, by = breaks_by)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) + 
    my.theme
  
  ## text output 
  text_out <- sprintf(
    "Rows used: %d\nECDF at x = %.3f: y = %.4f (%.1f%%)\nNatural-scale x prior to the log10(x+1) transformation: %.3f\n",
    n_used, shade_to, y_at_cut, 100 * y_at_cut, x_raw_at_cut
  )
  cat(text_out)
  writeLines(text_out, file.path(results, "juvenile_kelp.txt"))
  
  return(p)
}
## END function to create ecdf ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## view plot and save ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p1 <- plot_ecdf(
  df, 
  MacJuv_log10, 
  breaks_by = 0.2, 
  shade_to = log10(5)
  )

print(p1)


## save pdf 
ggsave(
  filename = file.path(figs, "juvenile_kelp.pdf"),
  plot = p1,
  device = cairo_pdf,
  dpi = 1200,
  width = 10,
  height = 10,
  units = "in"
)
## END view and save ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
