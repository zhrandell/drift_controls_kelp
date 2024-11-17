## load simulated data to plot
load(paste0(results, "/ODE_toPlot_kelp_low_", sel.model, ".RDA"))
load(paste0(results, "/ODE_toPlot_kelp_high_", sel.model, ".RDA"))

## graphical params ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c.drift <- "#00688B"
c.kelp <- "#308014"
c.fullness <- "#CD0000"
s.drift <- "#68228B"
s.kelp <- "#000000"
s.fullness <- "#CD0000"

lowD <- "#6183A6"
highD <- "#00688B"
lowK <- "#77896C"
highK <- "#006400"
lowF <- "#CD9B9B"
highF <- "#9D1309"

alph <- 0.5
linew <- 1
lwd_sml <- 0.75
wid1 <- .5
wid2 <- 1
lineCol <- "black"
ymax_loss <- 130
ymax_fill <- 45



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
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## Plot Low and High individually in final simulated figure ~~~~~~~~~~~~~~~~~~~~
## Drift loss
p1 <- ggplot(data = combined_high, aes(S0)) +
  geom_ribbon(aes(ymin = Smin_1, ymax = Smax_1),
              fill = highD,
              alpha = alph) + my.theme +
  geom_line(aes(S0, S_P1), color = lineCol, size = wid2) +
  geom_line(aes(S0, S_P1), color = highD, size = wid1) +
  geom_ribbon(
    data = combined_low,
    aes(ymin = Smin_1, ymax = Smax_1),
    fill = lowD,
    alpha = alph
  ) + my.theme +
  geom_line(data = combined_low,
            aes(S0, S_P1),
            color = lineCol,
            size = wid2) +
  geom_line(data = combined_low,
            aes(S0, S_P1),
            color = lowD,
            size = wid1) +
  xlab("initial Drift") +
  ylab("drift consumed") +
  ggtitle("period 1: t = 44hr snapshot") +
  ylim(0, ymax_loss) +
  x.blank

p2 <- ggplot(data = combined_high, aes(S0)) +
  geom_ribbon(aes(ymin = Smin_2, ymax = Smax_2),
              fill = highD,
              alpha = alph) + my.theme +
  geom_line(aes(S0, S_P2), color = lineCol, size = wid2) +
  geom_line(aes(S0, S_P2), color = highD, size = wid1) +
  geom_ribbon(
    data = combined_low,
    aes(ymin = Smin_2, ymax = Smax_2),
    fill = lowD,
    alpha = alph
  ) + my.theme +
  geom_line(data = combined_low,
            aes(S0, S_P2),
            color = lineCol,
            size = wid2) +
  geom_line(data = combined_low,
            aes(S0, S_P2),
            color = lowD,
            size = wid1) +
  xlab("initial Drift") + 
  ggtitle("period 2: t = 89hr snapshot") +
  ylim(0, ymax_loss) +
  both.blank

p3 <- ggplot(data = combined_high, aes(S0)) +
  geom_ribbon(aes(ymin = Smin_3, ymax = Smax_3),
              fill = highD,
              alpha = alph) + my.theme +
  geom_line(aes(S0, S_P3), color = lineCol, size = wid2) +
  geom_line(aes(S0, S_P3), color = highD, size = wid1) +
  geom_ribbon(
    data = combined_low,
    aes(ymin = Smin_3, ymax = Smax_3),
    fill = lowD,
    alpha = alph
  ) + my.theme +
  geom_line(data = combined_low,
            aes(S0, S_P3),
            color = lineCol,
            size = wid2) +
  geom_line(data = combined_low,
            aes(S0, S_P3),
            color = lowD,
            size = wid1) +
  xlab("initial Drift") +
  ggtitle("period 3: t = 134hr snapshot") +
  ylim(0, ymax_loss) +
  both.blank


## Kelp loss
p4 <- ggplot(data = combined_high, aes(S0)) +
  geom_ribbon(aes(ymin = Amin_1, ymax = Amax_1),
              fill = highK,
              alpha = alph) + my.theme +
  geom_line(aes(S0, A_P1), color = lineCol, size = wid2) +
  geom_line(aes(S0, A_P1), color = highK, size = wid1) +
  geom_ribbon(
    data = combined_low,
    aes(ymin = Amin_1, ymax = Amax_1),
    fill = lowK,
    alpha = alph
  ) +
  my.theme +
  geom_line(data = combined_low,
            aes(S0, A_P1),
            color = lineCol,
            size = wid2) +
  geom_line(data = combined_low,
            aes(S0, A_P1),
            color = lowK,
            size = wid1) +
  xlab("initial Drift") +
  ylab("kelp consumed") +
  ylim(0, ymax_loss) +
  x.blank

p5 <- ggplot(data = combined_high, aes(S0)) +
  geom_ribbon(aes(ymin = Amin_2, ymax = Amax_2),
              fill = highK,
              alpha = alph) + my.theme +
  geom_line(aes(S0, A_P2), color = lineCol, size = wid2) +
  geom_line(aes(S0, A_P2), color = highK, size = wid1) +
  geom_ribbon(
    data = combined_low,
    aes(ymin = Amin_2, ymax = Amax_2),
    fill = lowK,
    alpha = alph
  ) +
  my.theme +
  geom_line(data = combined_low,
            aes(S0, A_P2),
            color = lineCol,
            size = wid2) +
  geom_line(data = combined_low,
            aes(S0, A_P2),
            color = lowK,
            size = wid1) +
  xlab("initial Drift") +
  ylab("kelp consumed") +
  ylim(0, ymax_loss) +
  both.blank

p6 <- ggplot(data = combined_high, aes(S0)) +
  geom_ribbon(aes(ymin = Amin_3, ymax = Amax_3),
              fill = highK,
              alpha = alph) + my.theme +
  geom_line(aes(S0, A_P3), color = lineCol, size = wid2) +
  geom_line(aes(S0, A_P3), color = highK, size = wid1) +
  geom_ribbon(
    data = combined_low,
    aes(ymin = Amin_3, ymax = Amax_3),
    fill = lowK,
    alpha = alph
  ) +
  my.theme +
  geom_line(data = combined_low,
            aes(S0, A_P3),
            color = lineCol,
            size = wid2) +
  geom_line(data = combined_low,
            aes(S0, A_P3),
            color = lowK,
            size = wid1) +
  xlab("initial Drift") +
  ylab("kelp consumed")  +
  ylim(0, ymax_loss) +
  both.blank


## Gut fullness
p7 <- ggplot(data = combined_high, aes(S0)) +
  geom_ribbon(aes(ymin = Fmin_1, ymax = Fmax_1),
              fill = highF,
              alpha = alph) + 
  my.theme +
  geom_line(aes(S0, F_P1), color = lineCol, size = wid2) +
  geom_line(aes(S0, F_P1), color = highF, size = wid1) +
  geom_ribbon(
    data = combined_low,
    aes(ymin = Fmin_1, ymax = Fmax_1),
    fill = lowF,
    alpha = alph
  ) +
  my.theme +
  geom_line(data = combined_low,
            aes(S0, F_P1),
            color = lineCol,
            size = wid2) +
  geom_line(data = combined_low,
            aes(S0, F_P1),
            color = lowF,
            size = wid1) +
  xlab("initial Drift") +
  ylab("cumulative stomach fullness")
  ylim(0, ymax_fill)

p8 <- ggplot(data = combined_high, aes(S0)) +
  geom_ribbon(aes(ymin = Fmin_2, ymax = Fmax_2),
              fill = highF,
              alpha = alph) + my.theme +
  geom_line(aes(S0, F_P2), color = lineCol, size = wid2) +
  geom_line(aes(S0, F_P2), color = highF, size = wid1) +
  geom_ribbon(
    data = combined_low,
    aes(ymin = Fmin_2, ymax = Fmax_2),
    fill = lowF,
    alpha = alph
  ) +
  my.theme +
  geom_line(data = combined_low,
            aes(S0, F_P2),
            color = lineCol,
            size = wid2) +
  geom_line(data = combined_low,
            aes(S0, F_P2),
            color = lowF,
            size = wid1) +
  xlab("initial Drift") +
  ylab("drift consumed") +
  ylim(0, ymax_fill) +
  y.blank

p9 <- ggplot(data = combined_high, aes(S0)) +
  geom_ribbon(aes(ymin = Fmin_3, ymax = Fmax_3),
              fill = highF,
              alpha = alph) + my.theme +
  geom_line(aes(S0, F_P3), color = lineCol, size = wid2) +
  geom_line(aes(S0, F_P3), color = highF, size = wid1) +
  geom_ribbon(
    data = combined_low,
    aes(ymin = Fmin_3, ymax = Fmax_3),
    fill = lowF,
    alpha = alph
  ) +
  my.theme +
  geom_line(data = combined_low,
            aes(S0, F_P3),
            color = lineCol,
            size = wid2) +
  geom_line(data = combined_low,
            aes(S0, F_P3),
            color = lowF,
            size = wid1) +
  xlab("initial Drift") +
  ylab("drift consumed") +
  ylim(0, ymax_fill) +
  y.blank


## plot
all9 <- ggarrange(
  tag_facet(p1 + facet_wrap(~ "time"), tag_pool = "a"),
  tag_facet(p2 + facet_wrap(~ "time"), tag_pool = "b"),
  tag_facet(p3 + facet_wrap(~ "time"), tag_pool = "c"),
  tag_facet(p4 + facet_wrap(~ "time"), tag_pool = "d"),
  tag_facet(p5 + facet_wrap(~ "time"), tag_pool = "e"),
  tag_facet(p6 + facet_wrap(~ "time"), tag_pool = "f"),
  tag_facet(p7 + facet_wrap(~ "time"), tag_pool = "g"),
  tag_facet(p8 + facet_wrap(~ "time"), tag_pool = "h"),
  tag_facet(p9 + facet_wrap(~ "time"), tag_pool = "i"),
  nrow = 3,
  ncol = 3
)

print(all9)

ggplot2::ggsave(filename = paste0(figs, "/ODE_simulation_", sel.model, ".pdf"), 
                plot = all9, 
                dpi = 1200, 
                width = 11,
                height = 8, 
                units = "in")
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~