## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## function to plot simulated temporal dynamics code ~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## set up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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


## graphing details
alph <- 0.5
linew <- 1
lwd_sml <- 0.75
wid1 <- 0.5
wid2 <- 1
lineCol <- "black"

combined <- rbind(combined_high, combined_low)
ymax_loss <- max(combined[, c('Smin_1', 'Smax_1', 'Amin_1', 'Amax_1', 
                              'Smin_2', 'Smax_2', 'Amin_2', 'Amax_2',
                              'Smin_3', 'Smax_3', 'Amin_3', 'Amax_3',
                              'S_P1', 'S_P2', 'S_P3',
                              'A_P1', 'A_P2', 'A_P3')]) * 1.1
ymax_fill <- max(combined[, c('Fmin_1', 'Fmax_1',
                              'Fmin_1', 'Fmax_1',
                              'Fmin_3', 'Fmax_3',
                              'F_P1', 'F_P2', 'F_P3')]) * 1.1


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





## function to create custom legend via grob ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(grid)

# Function to return a grob for custom legend
custom_legend_grob <- function(low_col, high_col) {
  legend_df <- data.frame(
    x = 0.55,
    y = c(1.00, 0.87),
    label = c("Low", "High"),
    fill = c(low_col, high_col)
  )
  
  ggplot(legend_df, aes(x = x, y = y)) +
    geom_tile(aes(fill = fill),
              width = 0.2,
              height = 0.075,
              show.legend = FALSE) +
    geom_text(aes(label = label),
              hjust = 0,
              nudge_x = 0.25,
              size = 5) +
    scale_fill_identity() +
    theme_void() +
    coord_fixed(
      ratio = 1.5,
      xlim = c(0.5, 1.1),
      ylim = c(0.8, 1.05),  
      expand = FALSE
    )
}


## create the three legends
legend_drift <- ggplotGrob(custom_legend_grob(low_drift_col, high_drift_col))
legend_kelp <- ggplotGrob(custom_legend_grob(low_kelp_col, high_kelp_col))
legend_full <- ggplotGrob(custom_legend_grob(low_fullness_col, high_fullness_col))
## END custom legend ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





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
  plot_title =  "Period 1",
  ymax = ymax_loss
)


## edit axis labels as required for final, aggregated figure
drift.1 <- drift.1 + x.blank
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
  plot_title =  "Period 2",
  ymax = ymax_loss
)


## edit axis labels as required for final, aggregated figure
drift.2 <- drift.2 + x.blank + y.blank 
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
  plot_title =  "Period 3",
  ymax = ymax_loss
)


## edit axis labels as required for final, aggregated figure
drift.3 <- drift.3 + x.blank + y.blank +
  annotation_custom(grob = legend_drift, xmin = 200, xmax = 300, ymin = 80, ymax = 130)
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
kelp.3 <- kelp.3 + x.blank + y.blank +
  annotation_custom(grob = legend_kelp, xmin = 200, xmax = 300, ymin = 80, ymax = 130)
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
  y_axis_text = "Stomach fullness",
  plot_title =  element_blank(),
  ymax = ymax_fill
)
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
  y_axis_text = "Stomach fullness",
  plot_title =  element_blank(),
  ymax = ymax_fill
)


## edit axis labels as required for final, aggregated figure
fullness.2 <- fullness.2 + y.blank 
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
  y_axis_text = "Stomach fullness",
  plot_title =  element_blank(),
  ymax = ymax_fill
)


## edit axis labels as required for final, aggregated figure
fullness.3 <- fullness.3 + y.blank +
  annotation_custom(grob = legend_full, 
                    xmin = 200, xmax = 300, ymin = 0, ymax = 2)
## END Period 2 cumulative fullness ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## plot all 9 panes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all.temporal.dynamics <- ggarrange(
  tag_facet(drift.1 + facet_wrap(~ "time"), tag_pool = "a"),
  tag_facet(drift.2 + facet_wrap(~ "time"), tag_pool = "b"),
  tag_facet(drift.3 + facet_wrap(~ "time"), tag_pool = "c"),
  tag_facet(kelp.1 + facet_wrap(~ "time"), tag_pool = "d"),
  tag_facet(kelp.2 + facet_wrap(~ "time"), tag_pool = "e"),
  tag_facet(kelp.3 + facet_wrap(~ "time"), tag_pool = "f"),
  tag_facet(fullness.1 + facet_wrap(~ "time"), tag_pool = "g"),
  tag_facet(fullness.2 + facet_wrap(~ "time"), tag_pool = "h"),
  tag_facet(fullness.3 + facet_wrap(~ "time"), tag_pool = "i"),
  nrow = 3,
  ncol = 3
)


## view plot, if desired
#graphics.off()
#windows(11, 8)
#print(all.temporal.dynamics)


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
