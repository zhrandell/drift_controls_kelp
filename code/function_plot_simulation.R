## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## function to plot simulated temporal dynamics code ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## set up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)
library(rlang)


setwd("../")
getwd()
results <- "results"
figs <- "figs"

# Choose the model [1] or [2]
sel.model <- c('Logistic', 'vanLeeuwen')[2]
load(paste0(results,"/ODE_toPlot_kelp_high_", sel.model,".RDA"))
load(paste0(results,"/ODE_toPlot_kelp_low_", sel.model,".RDA"))

## END set up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## graphical parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## colors
low_drift_col <- "#6183A6"
high_drift_col <- "#00688B"
low_kelp_col <- "#77896C"
high_kelp_col <- "#006400"
low_fullness_col <- "#CD9B9B"
high_fullness_col <- "#9D1309"


## graphing details hard coded
alph <- 0.5
linew <- 1
lwd_sml <- 0.75
wid1 <- 0.5
wid2 <- 1
lineCol <- "black"
ymax_loss <- 130
ymax_fill <- 45


## axis options for final, aggregated figure
both.blank <- theme(
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.title.x = element_blank(),
  axis.text.x = element_blank()
)

x.blank <- theme(axis.title.x = element_blank(), 
                 axis.text.x = element_blank())
y.blank <- theme(axis.title.y = element_blank(), 
                 axis.text.y = element_blank())


## set up custom ggplot theme
my.theme = theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.title = element_text(size = 16),
  axis.text = element_text(size = 14),
  plot.title = element_text(size = 16)
)
## END custom graphing parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## function to plot single Period for Low and High Kelp treatments ~~~~~~~~~~~~~ 
plot.dynamics <- function(high_dat, var, 
                      high_lower_ribbon, 
                      high_upper_ribbon, 
                      high_fill, ## params to plot High Kelp ribbon
                      high_init, 
                      high_period, ## params to plot High Kelp median value line across x-axis 
                      low_dat, 
                      low_lower_ribbon, 
                      low_upper_ribbon, 
                      low_fill, ## params to plot Low Kelp ribbon
                      low_init, 
                      low_period, ## params to plot Low Kelp median value line across x-axis
                      x_axis_text, 
                      y_axis_text, 
                      plot_title,
                      ymax){ ## axis label and ylim params

## High Kelp treatment  
  plot <- ggplot(data = high_dat,
                 aes(x = {{ var }})) + 
    my.theme +
    geom_ribbon(aes(
      ymin = {{ high_lower_ribbon }}, 
      ymax = {{ high_upper_ribbon }}),
      fill = high_fill, 
      alpha = alph
      ) +
    geom_line(aes(
      x = {{ high_init }}, 
      y = {{ high_period }}), 
      color = lineCol, 
      linewidth = wid2) +

## Low Kelp treatment
    geom_ribbon(data = low_dat, 
                aes(
                  ymin = {{ low_lower_ribbon }}, 
                  ymax = {{ low_upper_ribbon }}),
                fill = low_fill, 
                alpha = alph 
                ) +
    geom_line(data = low_dat, 
              aes(
                x = {{ low_init }}, 
                y = {{ low_period }}
                ), 
      color = lineCol, 
      linewidth = wid2) +

## overarching pane details
    xlab(x_axis_text) +
    ylab(y_axis_text) +
    ggtitle(plot_title) +
    ylim(0, ymax)

  return(plot)
}
## END function to plot single pane of temporal dynamics fig ~~~~~~~~~~~~~~~~~~~





## invoke function to plot Period 1 high & low drift consumed at Period 1 ~~~~~~
drift.1 <- plot.dynamics(
  
  high_dat = combined_high,
  var = S0,
  high_lower_ribbon = Smin_1, 
  high_upper_ribbon = Smax_1, 
  high_fill = high_drift_col,
  high_init = S0,
  high_period = S_P1, 
  
  low_dat = combined_low,
  low_lower_ribbon = Smin_1, 
  low_upper_ribbon = Smax_1, 
  low_fill = low_drift_col,
  low_init = S0,
  low_period = S_P1, 
  
  x_axis_text = "Initial drift",
  y_axis_text = "Drift consumed",
  plot_title =  "Period 1: t = 44hr snapshot",
  ymax = ymax_loss
)


## edit axis labels as required for final, aggregated figure
drift.1 <- drift.1 + x.blank
print(drift.1)
## END Period 1 drift consumption ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## invoke function to plot Period 2 high & low drift consumed ~~~~~~~~~~~~~~~~~~
drift.2 <- plot.dynamics(
  
  high_dat = combined_high,
  var = S0,
  high_lower_ribbon = Smin_2, 
  high_upper_ribbon = Smax_2, 
  high_fill = high_drift_col,
  high_init = S0,
  high_period = S_P2, 
  
  low_dat = combined_low,
  low_lower_ribbon = Smin_2, 
  low_upper_ribbon = Smax_2, 
  low_fill = low_drift_col,
  low_init = S0,
  low_period = S_P2, 
  
  x_axis_text = "Initial drift",
  y_axis_text = "Drift consumed",
  plot_title =  "Period 2: t = 89hr snapshot",
  ymax = ymax_loss
)


## edit axis labels as required for final, aggregated figure
drift.2 <- drift.2 + x.blank + y.blank 
print(drift.2)
## END Period 2 drift consumption ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## invoke function to plot Period 3 high & low drift consumed ~~~~~~~~~~~~~~~~~~
drift.3 <- plot.dynamics(
  
  high_dat = combined_high,
  var = S0,
  high_lower_ribbon = Smin_3, 
  high_upper_ribbon = Smax_3, 
  high_fill = high_drift_col,
  high_init = S0,
  high_period = S_P3, 
  
  low_dat = combined_low,
  low_lower_ribbon = Smin_3, 
  low_upper_ribbon = Smax_3, 
  low_fill = low_drift_col,
  low_init = S0,
  low_period = S_P3, 
  
  x_axis_text = "Initial drift",
  y_axis_text = "Drift consumed",
  plot_title =  "Period 3: t = 134hr snapshot",
  ymax = ymax_loss
)


## edit axis labels as required for final, aggregated figure
drift.3 <- drift.3 + x.blank + y.blank 
print(drift.3)
## END Period 3 drift consumption ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




## invoke function to plot Period 1 high & low kelp consumed ~~~~~~~~~~~~~~~~~~~
kelp.1 <- plot.dynamics(
  
  high_dat = combined_high,
  var = S0,
  high_lower_ribbon = Amin_1, 
  high_upper_ribbon = Amax_1, 
  high_fill = high_kelp_col,
  high_init = S0,
  high_period = A_P1, 
  
  low_dat = combined_low,
  low_lower_ribbon = Amin_1, 
  low_upper_ribbon = Amax_1, 
  low_fill = low_kelp_col,
  low_init = S0,
  low_period = A_P1, 
  
  x_axis_text = "Initial drift",
  y_axis_text = "Kelp consumed",
  plot_title =  element_blank(),
  ymax = ymax_loss
)


## edit axis labels as required for final, aggregated figure
kelp.1 <- kelp.1 + x.blank #+ y.blank 
print(kelp.1)
## END Period 1 kelp consumption ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## invoke function to plot Period 2 high & low kelp consumed ~~~~~~~~~~~~~~~~~~~
kelp.2 <- plot.dynamics(
  
  high_dat = combined_high,
  var = S0,
  high_lower_ribbon = Amin_2, 
  high_upper_ribbon = Amax_2, 
  high_fill = high_kelp_col,
  high_init = S0,
  high_period = A_P2, 
  
  low_dat = combined_low,
  low_lower_ribbon = Amin_2, 
  low_upper_ribbon = Amax_2, 
  low_fill = low_kelp_col,
  low_init = S0,
  low_period = A_P2, 
  
  x_axis_text = "Initial drift",
  y_axis_text = "Kelp consumed",
  plot_title =  element_blank(),
  ymax = ymax_loss
)


## edit axis labels as required for final, aggregated figure
kelp.2 <- kelp.2 + x.blank + y.blank 
print(kelp.2)
## END Period 2 kelp consumption ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## invoke function to plot Period 3 high & low kelp consumed ~~~~~~~~~~~~~~~~~~~
kelp.3 <- plot.dynamics(
  
  high_dat = combined_high,
  var = S0,
  high_lower_ribbon = Amin_3, 
  high_upper_ribbon = Amax_3, 
  high_fill = high_kelp_col,
  high_init = S0,
  high_period = A_P3, 
  
  low_dat = combined_low,
  low_lower_ribbon = Amin_3, 
  low_upper_ribbon = Amax_3, 
  low_fill = low_kelp_col,
  low_init = S0,
  low_period = A_P3, 
  
  x_axis_text = "Initial drift",
  y_axis_text = "Kelp consumed",
  plot_title =  element_blank(),
  ymax = ymax_loss
)


## edit axis labels as required for final, aggregated figure
kelp.3 <- kelp.3 + x.blank + y.blank 
print(kelp.3)
## END Period 3 kelp consumption ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## invoke function to plot Period 1 cumulative fullness ~~~~~~~~~~~~~~~~~~~~~~~~
fullness.1 <- plot.dynamics(
  
  high_dat = combined_high,
  var = S0,
  high_lower_ribbon = Fmin_1, 
  high_upper_ribbon = Fmax_1, 
  high_fill = high_fullness_col,
  high_init = S0,
  high_period = F_P1, 
  
  low_dat = combined_low,
  low_lower_ribbon = Fmin_1, 
  low_upper_ribbon = Fmax_1, 
  low_fill = low_fullness_col,
  low_init = S0,
  low_period = F_P1, 
  
  x_axis_text = "Initial drift",
  y_axis_text = "Cumulative stomach fullness",
  plot_title =  element_blank(),
  ymax = ymax_fill
)


## edit axis labels as required for final, aggregated figure
#fullness.1 <- fullness.1 + x.blank + y.blank 
print(fullness.1)
## END Period 1 cumulative fullness ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## invoke function to plot Period 2 cumulative fullness ~~~~~~~~~~~~~~~~~~~~~~~~
fullness.2 <- plot.dynamics(
  
  high_dat = combined_high,
  var = S0,
  high_lower_ribbon = Fmin_2, 
  high_upper_ribbon = Fmax_2, 
  high_fill = high_fullness_col,
  high_init = S0,
  high_period = F_P2, 
  
  low_dat = combined_low,
  low_lower_ribbon = Fmin_2, 
  low_upper_ribbon = Fmax_2, 
  low_fill = low_fullness_col,
  low_init = S0,
  low_period = F_P2, 
  
  x_axis_text = "Initial drift",
  y_axis_text = "Cumulative stomach fullness",
  plot_title =  element_blank(),
  ymax = ymax_fill
)


## edit axis labels as required for final, aggregated figure
fullness.2 <- fullness.2 + y.blank 
print(fullness.2)
## END Period 2 cumulative fullness ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## invoke function to plot Period 2 cumulative fullness ~~~~~~~~~~~~~~~~~~~~~~~~
fullness.3 <- plot.dynamics(
  
  high_dat = combined_high,
  var = S0,
  high_lower_ribbon = Fmin_3, 
  high_upper_ribbon = Fmax_3, 
  high_fill = high_fullness_col,
  high_init = S0,
  high_period = F_P3, 
  
  low_dat = combined_low,
  low_lower_ribbon = Fmin_3, 
  low_upper_ribbon = Fmax_3, 
  low_fill = low_fullness_col,
  low_init = S0,
  low_period = F_P3, 
  
  x_axis_text = "Initial drift",
  y_axis_text = "Cumulative stomach fullness",
  plot_title =  element_blank(),
  ymax = ymax_fill
)


## edit axis labels as required for final, aggregated figure
fullness.3 <- fullness.3 + y.blank 
print(fullness.3)
## END Period 2 cumulative fullness ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## plot all 9 panes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)
library(egg)

all.temporal.dynamics <- ggarrange(
  tag_facet(drift.1 + facet_wrap(~ "time"), tag_pool = "a"),
  tag_facet(drift.2 + facet_wrap(~ "time"), tag_pool = "b"),
  tag_facet(drift.3 + facet_wrap(~ "time"), tag_pool = "c"),
  tag_facet(kelp.1 + facet_wrap(~ "time"), tag_pool = "d"),
  tag_facet(kelp.2 + facet_wrap(~ "time"), tag_pool = "e"),
  tag_facet(kelp.3 + facet_wrap(~ "time"), tag_pool = "f"),
  tag_facet(fullness.1 + facet_wrap(~ "time"), tag_pool = "d"),
  tag_facet(fullness.2 + facet_wrap(~ "time"), tag_pool = "e"),
  tag_facet(fullness.3 + facet_wrap(~ "time"), tag_pool = "f"),
  nrow = 3,
  ncol = 3
)

windows(11, 8)
print(all.temporal.dynamics)


## save pdf 
ggplot2::ggsave(filename = paste0(figs, "/ODE_simulation_", sel.model, ".pdf"), 
                plot = all.temporal.dynamics, 
                dpi = 1200, 
                width = 11,
                height = 8, 
                units = "in")
## END plotting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
