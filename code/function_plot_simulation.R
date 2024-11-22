## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## function to plot simulated temporal dynamics code ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## set up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(tidyverse)
library(rlang)

load("ODE_toPlot_kelp_high_vanLeeuwen.RDA")
load("ODE_toPlot_kelp_low_vanLeeuwen.RDA")
## END set up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## graphical parameters ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## colors
low_drift_col <- "#6183A6"
high_drift_col <- "#00688B"
low_kelp_col <- "#77896C"
high_kelp_col <- "#006400"
low_fullness_col <- "#CD9B9B"
high_fullnes_col <- "#9D1309"


## graphing details hard coded
alph <- 0.5
linew <- 1
lwd_sml <- 0.75
wid1 <- .5
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

x.blank <- theme(axis.title.x = element_blank(), axis.text.x = element_blank())
y.blank <- theme(axis.title.y = element_blank(), axis.text.y = element_blank())


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
test.plot <- function(high_dat, var, high_lower_ribbon, high_upper_ribbon, high_fill, ## params to plot High Kelp ribbon
                      high_init, high_period, ## params to plot High Kelp median value line across x-axis 
                      low_dat, low_lower_ribbon, low_upper_ribbon, low_fill, ## params to plot Low Kelp ribbon
                      low_init, low_period, ## params to plot Low Kelp median value line across x-axis
                      x_axis_text, y_axis_text, plot_title, ymax){ ## axis label and ylim params

## High Kelp treatment  
  plot <- ggplot(data = high_dat, aes(x = {{ var }})) + my.theme +
    geom_ribbon(aes(
      ymin = {{ high_lower_ribbon }}, 
      ymax = {{ high_upper_ribbon }}),
      fill = high_fill, 
      alpha = alph) +
    geom_line(aes(
      x = {{ high_init }}, 
      y = {{ high_period }}), 
      color = lineCol, 
      linewidth = wid2) +

## Low Kelp treatment
    geom_ribbon(data = low_dat, aes(
      ymin = {{ low_lower_ribbon }}, 
      ymax = {{ low_upper_ribbon }}),
      fill = low_fill, 
      alpha = alph) +
    geom_line(data = low_dat, aes(
      x = {{ low_init }}, 
      y = {{ low_period }}), 
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





## invoke function to plot a single pane ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p1 <- test.plot(
  
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
  
  x_axis_text = "initial Drift",
  y_axis_text = "Drift consumed",
  plot_title =  "period 1: t = 44hr snapshot",
  ymax = ymax_loss
)


## plot single figure
windows(4, 4, record=T)
print(p1)


## edit axis labels as required for final, aggregated figure
p1 <- p1 + x.blank
print(p1)
## END graphing function and invocation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
