## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Empirical data visualization
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## startup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## read in data
dat <- read.csv(file.path(data, "drift_kelp_loss.csv"))
## END startup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## format data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Subtract initial wet-weight of bags used in the field 
## -20g wet weight for Drift bags 
dat$Drift_Initial <- dat$Drift_Initial - 20
dat$Drift_Remaining <- dat$Drift_Remaining - 20
dat$Drift_Consumed <- dat$Drift_Initial - dat$Drift_Remaining


## -80g wet weight for kelp bags
dat$Kelp_Initial <- dat$Kelp_Initial - 80
dat$Kelp_Remaining <- dat$Kelp_Remaining - 80
dat$Kelp_Consumed <- dat$Kelp_Initial - dat$Kelp_Remaining


## total consumption  
dat$Total_Initial <- dat$Drift_Initial + dat$Kelp_Initial 
dat$Total_Consumed <- dat$Drift_Consumed + dat$Kelp_Consumed 


## preserve original treatment labels
dat$Treatment_Controls <- dat$Treatment  


## incorporate controls 
dat$Treatment[dat$Treatment == "High_Control"] <- "High"
dat$Treatment[dat$Treatment == "Low_Control"] <- "Low"


## create new column seq for sequence 
dat$seq <- ifelse(dat$Trial <= 4, "seq1", "seq2")


## assign hours
dat <- dat %>%
  mutate(hours = case_when(
    seq == "seq1" & Period == 1 ~ 48,
    seq == "seq1" & Period == 2 ~ 96,
    seq == "seq1" & Period == 3 ~ 144,
    seq == "seq2" & Period == 1 ~ 24,
    seq == "seq2" & Period == 2 ~ 48,
    TRUE ~ NA_real_
  ))


## assign group labels using sequence, treatment, and urchin count
dat <- dat %>%
  mutate(group_label = case_when(
    seq == "seq1" & Treatment == "Low"  ~ "Sequence 1, 20 urchins, Low",
    seq == "seq1" & Treatment == "High" ~ "Sequence 1, 20 urchins, High",
    seq == "seq2" & Urchins == 10 & Treatment == "Low"  ~ "Sequence 2, 10 urchins, Low",
    seq == "seq2" & Urchins == 10 & Treatment == "High" ~ "Sequence 2, 10 urchins, High",
    seq == "seq2" & Urchins == 20 & Treatment == "Low"  ~ "Sequence 2, 20 urchins, Low",
    seq == "seq2" & Urchins == 20 & Treatment == "High" ~ "Sequence 2, 20 urchins, High",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(hours), !is.na(group_label))


## define group labels for total consumption fig
dat$group_label_simple <- case_when(
  dat$seq == "seq1" ~ "Sequence 1, 48 hrs, 20 urchins",
  dat$seq == "seq2" & dat$Urchins == 20 ~ "Sequence 2, 24 hrs, 20 urchins",
  dat$seq == "seq2" & dat$Urchins == 10 ~ "Sequence 2, 24 hrs, 10 urchins"
)


## Remove rows with NA values 
dat <- dat %>%
  filter(!is.na(Drift_Initial),
         !is.na(Drift_Consumed),
         !is.na(Kelp_Initial),
         !is.na(Kelp_Consumed))


## to plot only sequence 1 data 
seq1 <- dat %>% filter(Trial %in% 1:4)


## drop un-used Treatments (urchin control, drift control)
seq1 <- seq1 %>% filter(Treatment %in% c("Low", "High")) %>% droplevels()
## END data prep ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## graphing params ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## assign colors to groups - group names must match those above 
pal_custom <- c(
  "Sequence 1, 48 hrs, 20 urchins" = "#009E73",  # Darker Orange
  "Sequence 2, 24 hrs, 20 urchins" = "#003F87",  # Sky Blue
  "Sequence 2, 24 hrs, 10 urchins" = "#D55E00"   # Bluish Green
)


## visual/text options
my.theme = theme(
  panel.grid.major = element_blank(), 
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"))


## custom text size
my.text = theme(
  axis.title = element_text(size = 15),
  axis.text = element_text(size = 14),
  plot.title = element_text(size = 16),
  strip.text = element_text(size = 14),
  legend.text = element_text(size = 15),
  legend.title = element_text(size = 15)
)


## upper right legend
upper.R.legend = theme(legend.position.inside = c(1, 1),
                       legend.justification = c(1, 1),
                       legend.background = element_rect(fill = "white", color = NA))


## legend version with more space
upper.R.legend.2 <- theme(
  legend.position = "inside",
  legend.position.inside = c(0.98, 0.98),
  legend.justification = c(1, 1),
  legend.background = element_rect(fill = "white", color = NA),
  legend.box.margin = margin(5,5,5,5)
)


## disable legend 
no.legend = theme(legend.position="none")


## add text label to facet_wrap() 
facet.label <- labeller(hours = function(x) paste0(x, " hours"))


## graphing parameters
pal_drift <- c("#00688B", "#6183A6")
pal_kelp <- c("#006400", "#77896C")
pal_seq <- c("#153d52","#96312e","#ffbf6b")
blk <- "black"
c.drift <- "#68228B"
c.kelp <- "#000000"
ptType <- 19
ptSize <- 1.25
ptCol <- "black"
ptFill <- "black"
alph <- 0.8
smoothSize <- 0.5
smoothSpan <- 1
smoothCol <- "black"
smoothCI <- 0.95
## END graphing params ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## create fig. 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## define function to create fig. 2 
fig2 <- function(data,
                 xvar,
                 yvar,
                 palette,
                 plot_title,
                 xlab_text,
                 ylab_text,
                 show_legend = TRUE) {
  
  legend_theme <- if (show_legend) {
    upper.R.legend
  } else {
    no.legend
  }
  
  ggplot(data, aes(x = {{ xvar }}, y = {{ yvar }})) + my.theme + my.text +
    geom_point(aes(color = Treatment, shape = Treatment), size = ptSize) +
    geom_smooth(
      data = data %>% filter(!Treatment_Controls %in% c("High_Control", "Low_Control")),
      aes(color = Treatment, fill = Treatment),
      linewidth = smoothSize, span = smoothSpan, level = smoothCI
    ) +
    scale_color_manual(values = palette) +
    scale_fill_manual(values = palette) +
    scale_shape_manual(values = c("Low" = 19, "High" = 17)) + 
    labs(title = plot_title) +
    ylab(ylab_text) + xlab(xlab_text) +
    facet_grid(cols = vars(hours), labeller = facet.label) +
    scale_y_continuous(limit = c(0, NA), oob = squish) +
    legend_theme
}


## create rows of fig 2 
p1 <- fig2(
  data = seq1,
  xvar = Drift_Initial,
  yvar = Drift_Consumed,
  palette = pal_drift,
  plot_title = "Drift consumed vs Drift available",
  xlab_text = "initial Drift (g)",
  ylab_text = "Drift consumed (g)",
  show_legend = TRUE
)


p2 <- fig2(
  data = seq1,
  xvar = Drift_Initial,
  yvar = Kelp_Consumed,
  palette = pal_kelp,
  plot_title = "Kelp consumed vs Drift available",
  xlab_text = "initial Drift (g)",
  ylab_text = "Kelp consumed (g)",
  show_legend = TRUE
)


p3 <- fig2(
  data = seq1,
  xvar = Kelp_Initial,
  yvar = Kelp_Consumed,
  palette = pal_kelp,
  plot_title = "Kelp consumed vs Kelp available",
  xlab_text = "initial Kelp (g)",
  ylab_text = "Kelp consumed (g)",
  show_legend = FALSE
)


## combine all rows into 1
seq1_results <- ggarrange(p1, p2, p3, nrow = 3)


## save 
ggsave(paste0(figs,"/fig2.pdf"), 
       plot = seq1_results, 
       width = 10, 
       height = 12, 
       dpi = 1200, 
       units = "in")
## END fig. 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## all consumption (drift + kelp) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all_levels <- names(pal_custom)
dat$group_label_simple <- factor(dat$group_label_simple, levels = all_levels)


smooth_data <- dat %>%
  filter(!Treatment_Controls %in% c("High_Control", "Low_Control")) %>%
  mutate(group_label_simple = factor(group_label_simple, levels = all_levels))


all_consumption <- ggplot(dat, aes(x = Total_Initial, y = Total_Consumed)) +
  my.theme +
  geom_point(aes(color = group_label_simple, fill = group_label_simple),
             pch = ptType, size = ptSize, alpha = alph) +
  geom_smooth(
    data = smooth_data,
    aes(x = Total_Initial, y = Total_Consumed,
        color = group_label_simple, fill = group_label_simple, group = group_label_simple),
    method = "loess",
    formula = y ~ x,
    span = smoothSpan,
    linewidth = smoothSize,
    level = smoothCI
  ) +
  scale_color_manual(name = NULL, values = pal_custom, drop = FALSE) +
  scale_fill_manual(name = NULL, values = pal_custom, drop = FALSE) +
  guides(fill = guide_legend(nrow = 3), color = "none") +
  labs(
    title = "Total biomass (kelp + drift) consumed across all eight experimental Trials",
    x = "Biomass available",
    y = "Biomass consumed"
  ) +
  facet_grid(cols = vars(hours), labeller = facet.label) +
  scale_y_continuous(limit = c(-1, NA), oob = squish) +
  my.text + upper.R.legend.2


## save
ggsave(paste0(figs,"/SOM_all_consumption.pdf"), 
       plot = all_consumption, 
       width = 12, 
       height = 7, 
       dpi = 1200, 
       units = "in")
## END all consumption figure ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## SOM figure with all sequences plotted; seq 2 highlighted; seq 1 grayed ~~~~~~
## create groups for text labels and color
dat <- dat %>%
  mutate(
    group_label = case_when(
      seq == "seq1" & Treatment == "Low"  ~ "Sequence 1, 20 urchins, Low",
      seq == "seq1" & Treatment == "High" ~ "Sequence 1, 20 urchins, High",
      seq == "seq2" & Urchins == 10 & Treatment == "Low"  ~ "Sequence 2, 10 urchins, Low",
      seq == "seq2" & Urchins == 10 & Treatment == "High" ~ "Sequence 2, 10 urchins, High",
      seq == "seq2" & Urchins == 20 & Treatment == "Low"  ~ "Sequence 2, 20 urchins, Low",
      seq == "seq2" & Urchins == 20 & Treatment == "High" ~ "Sequence 2, 20 urchins, High"
    )
  )


## colors for each group 
pal_custom <- c(
  "Sequence 1, 20 urchins, Low"  = "gray80",
  "Sequence 1, 20 urchins, High" = "gray20",
  "Sequence 2, 10 urchins, Low"  = "#CD3700",  # light orange
  "Sequence 2, 10 urchins, High" = "#8B0000",  # dark orange
  "Sequence 2, 20 urchins, Low"  = "#92c5de",  # light blue
  "Sequence 2, 20 urchins, High" = "#2166ac"   # dark blue
)


## function to plot all sequences
all_seq <- function(data,
                    xvar,
                    yvar,
                    palette,
                    plot_title,
                    xlab_text,
                    ylab_text,
                    show_legend = TRUE) {
  
  legend_elements <- if (show_legend) {
    list(
      theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.background = element_rect(fill = "white", color = NA),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15)
      ),
      guides(
        color = guide_legend(nrow = 2),
        fill = guide_legend(nrow = 2)
      )
    )
  } else {
    list(
      theme(legend.position = "none")
    )
  }
  
  ggplot(data, aes(x = {{ xvar }}, y = {{ yvar }})) + my.theme + my.text +
    geom_point(aes(color = group_label, fill = group_label, shape = Treatment), size = ptSize) +
    geom_smooth(
      data = data %>% filter(!Treatment_Controls %in% c("High_Control", "Low_Control")),
      aes(x = {{ xvar }}, y = {{ yvar }}, color = group_label, fill = group_label),
      linewidth = smoothSize, span = smoothSpan, level = smoothCI
    ) +
    scale_color_manual(values = palette) +
    scale_fill_manual(values = palette) +
    scale_shape_manual(values = c("Low" = 19, "High" = 17)) + 
    labs(title = plot_title, color = NULL, fill = NULL) +
    ylab(ylab_text) + xlab(xlab_text) +
    facet_grid(cols = vars(hours), labeller = facet.label) +
    scale_y_continuous(limit = c(0, NA), oob = squish) +
    legend_elements
}



## invoke function - drift as a function of drift
p4 <- all_seq(
  data = dat,
  xvar = Drift_Initial,
  yvar = Drift_Consumed,
  palette = pal_custom,
  plot_title = "Drift consumed vs available drift",
  xlab_text = "initial drift (g)",
  ylab_text = "drift consumed (g)",
  show_legend = FALSE
)

## invoke function - kelp as a function of drift 
p5 <- all_seq(
  data = dat,
  xvar = Drift_Initial,
  yvar = Kelp_Consumed,
  palette = pal_custom,
  plot_title = "Kelp consumed vs available drift",
  xlab_text = "initial drift (g)",
  ylab_text = "kelp consumed (g)",
  show_legend = FALSE
)


## invoke function - kelp as a function of kelp
p6 <- all_seq(
  data = dat,
  xvar = Kelp_Initial,
  yvar = Kelp_Consumed,
  palette = pal_custom,
  plot_title = "Kelp consumed vs available kelp",
  xlab_text = "initial kelp (g)",
  ylab_text = "kelp consumed (g)",
  show_legend = TRUE
)


## aggregate
pAll <- ggarrange(p4, p5, p6,
                  nrow = 3,
                  common.legend = TRUE,
                  legend = "bottom")


## save
ggsave(paste0(figs,"/SOM_both_sequences.pdf"), 
       plot = pAll, 
       width = 13, 
       height = 13, 
       dpi = 1200, 
       units = "in")
## END SOM figure with all seq ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
