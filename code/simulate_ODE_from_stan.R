## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Simulate ODE system with posteriors from STAN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##





## initialize ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list=ls())
  
library(deSolve)
library(reshape2)
library(egg)
library(tidyselect)
library(scales)

## set relative file paths
setwd("../")
getwd()

results <- "results"
code <- "code"
figs <- "figs"
## END initialization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




## ode() with multiple params ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## load posts_Df
setwd(results)
load("posts_new_All-vanLeeuwen.RDA")


## select every 10th row to downsample to 1000 
#posts_df <- posts_df_raw[,1:5]
#posts_df <- posts_df[seq(1, nrow(posts_df), 10), ]



## create sequences of initial conditions
len_init <- 200
S0 <- seq(1, 300, length.out = len_init)
A0 <- seq(300, 300, length.out = len_init)
H0 <- seq(0, 0, length.out = len_init)


## set time points (in hrs) for the three Periods
P1 <- 44 # 44 hrs
P2 <- 45 # 89 hrs
P3 <- 45 # 134 hrs


## time sequences to pass to ode()
t.list_P1 <- seq(1, P1, by = 1)
t.list_P2 <- seq(1, P2, by = 1)
t.list_P3 <- seq(1, P3, by = 1)


## set initial conditions 
init_P1 = c(S = S0[1], 
            A = A0[1], 
            H = H0[1])


## concactonate into list of lists 
inits_P1 <- mapply(c, S=S0, A=A0, H=H0, SIMPLIFY = F)


## list of params 
aL <- posts_df_raw$a 
vL <- posts_df_raw$v
wL <- posts_df_raw$w
qL <- posts_df_raw$q


## first set of params 
parm_list = c(a = aL[1],
              v = vL[1],
              w = wL[1],
              q = qL[1])


## specify ODE function 
resourceLoss <- function(S0, A0, H0, params) {
  with(as.list(c(S0, A0, H0)),{
 
    f_S <- S * a * (1 - (1 + exp(w + log(S / A))) / (1 + exp(log(2) + w + log(S / A)) + exp(q + 2 * log(S / A)))) 
    f_A <- A * a * ((1 + exp(w + log(S / A))) / (1 + exp(log(2) + w + log(S / A)) + exp(q + 2 * log(S / A)))) 
    Hunger = exp(- v * H)
    
    dS_dt = f_S * Hunger
    dA_dt = f_A * Hunger
    dH_dt = f_S + f_A
    
    return(list(c(dS_dt, dA_dt, dH_dt)))
  })
}


## run single ODE 
out_P1 <- ode(init_P1, 
              times=t.list_P1,
              func=resourceLoss,
              parms=parm_list)


## concactonate params into list of lists
full_parm_list <- mapply(c, a=aL, v=vL, w=wL, q=qL, SIMPLIFY = F)


## t.lists for restocking model
t.list_P1_restock <- seq(1, P1, by = 1)
t.list_P2_restock <- seq(P1+1, P1+P2, by = 1)
t.list_P3_restock <- seq(P1+P2+1, P1+P2+P3, by = 1)


## nested lapply 
outs_parms <- lapply(full_parm_list, function(x){
  lapply(inits_P1, function(y){
    
    p1 <- ode(y,
              times = t.list_P1_restock,
              func = resourceLoss,
              parms = x)
    
    p2 <- ode(c(y[1], y[2], p1[P1,4]), 
              times = t.list_P2_restock,
              func = resourceLoss,
              parms = x)
    
    p3 <- ode(c(y[1], y[2], p2[P2,4]), 
              times = t.list_P3_restock,
              func = resourceLoss,
              parms = x)
    
    return(rbind(p1,p2,p3))

  })
})

#setwd(simData)
save(outs_parms, file="ODE_vL_high-kelp.RDA")
## END ODE simulation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## calculate and plot 95% CI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd(results)
load("ODE_vL_low-kelp.RDA")
load("ODE_vL_high-kelp.RDA")


## create sequences of initial conditions
len_init <- 200
S0 <- seq(1, 20, length.out = len_init)
A0 <- seq(20, 20, length.out = len_init)
H0 <- seq(0, 0, length.out = len_init)


## set time points (in hrs) for the three Periods
P1 <- 44 # 44 hrs
P2 <- 45 # 89 hrs
P3 <- 45 # 134 hrs


## functions to extract state values from list of lists 
extract_Sloss_P1 = function(y){y[[1,1]] <- y[[1,2]] - y[[P1,2]]}
extract_Aloss_P1 = function(y){y[[1,1]] <- y[[1,3]] - y[[P1,3]]}
extract_Hfill_P1 = function(y){y[[1,1]] <- (-1*(y[[1,4]] - y[[P1,4]]))}

extract_Sloss_P2 = function(y){y[[1,1]] <- y[[1,2]] - y[[(P1+P2),2]]}
extract_Aloss_P2 = function(y){y[[1,1]] <- y[[1,3]] - y[[(P1+P2),3]]}
extract_Hfill_P2 = function(y){y[[1,1]] <- (-1*(y[[1,4]] - y[[(P1+P2),4]]))}

extract_Sloss_P3 = function(y){y[[1,1]] <- y[[1,2]] - y[[(P1+P2+P3),2]]}
extract_Aloss_P3 = function(y){y[[1,1]] <- y[[1,3]] - y[[(P1+P2+P3),3]]}
extract_Hfill_P3 = function(y){y[[1,1]] <- (-1*(y[[1,4]] - y[[(P1+P2+P3),4]]))}


## Drift loss
S_loss_P1 <- lapply(outs_parms, function(x)
  lapply(x, function(y) extract_Sloss_P1(y)))
A_loss_P1 <- lapply(outs_parms, function(x)
  lapply(x, function(y) extract_Aloss_P1(y)))
H_fill_P1 <- lapply(outs_parms, function(x)
  lapply(x, function(y) extract_Hfill_P1(y)))


S_loss_P2 <- lapply(outs_parms, function(x)
  lapply(x, function(y) extract_Sloss_P2(y)))
A_loss_P2 <- lapply(outs_parms, function(x)
  lapply(x, function(y) extract_Aloss_P2(y)))
H_fill_P2 <- lapply(outs_parms, function(x)
  lapply(x, function(y) extract_Hfill_P2(y)))


S_loss_P3 <- lapply(outs_parms, function(x)
  lapply(x, function(y) extract_Sloss_P3(y)))
A_loss_P3 <- lapply(outs_parms, function(x)
  lapply(x, function(y) extract_Aloss_P3(y)))
H_fill_P3 <- lapply(outs_parms, function(x)
  lapply(x, function(y) extract_Hfill_P3(y)))


## convert to df 
S_loss_P1 <- na.omit(as.data.frame(matrix(unlist(do.call(rbind, S_loss_P1)), ncol=len_init)))
A_loss_P1 <- na.omit(as.data.frame(matrix(unlist(do.call(rbind, A_loss_P1)), ncol=len_init)))
H_fill_P1 <- na.omit(as.data.frame(matrix(unlist(do.call(rbind, H_fill_P1)), ncol=len_init)))

S_loss_P2 <- na.omit(as.data.frame(matrix(unlist(do.call(rbind, S_loss_P2)), ncol=len_init)))
A_loss_P2 <- na.omit(as.data.frame(matrix(unlist(do.call(rbind, A_loss_P2)), ncol=len_init)))
H_fill_P2 <- na.omit(as.data.frame(matrix(unlist(do.call(rbind, H_fill_P2)), ncol=len_init)))

S_loss_P3 <- na.omit(as.data.frame(matrix(unlist(do.call(rbind, S_loss_P3)), ncol=len_init)))
A_loss_P3 <- na.omit(as.data.frame(matrix(unlist(do.call(rbind, A_loss_P3)), ncol=len_init)))
H_fill_P3 <- na.omit(as.data.frame(matrix(unlist(do.call(rbind, H_fill_P3)), ncol=len_init)))


## rarify to compensate for NA's produced during simulation 
rarify <- 990
S_loss_P1 <- S_loss_P1[sample(1:nrow(S_loss_P1), rarify), ]
A_loss_P1 <- A_loss_P1[sample(1:nrow(A_loss_P1), rarify), ]
H_fill_P1 <- H_fill_P1[sample(1:nrow(H_fill_P1), rarify), ]

S_loss_P2 <- S_loss_P2[sample(1:nrow(S_loss_P1), rarify), ]
A_loss_P2 <- A_loss_P2[sample(1:nrow(A_loss_P1), rarify), ]
H_fill_P2 <- H_fill_P2[sample(1:nrow(H_fill_P1), rarify), ]

S_loss_P3 <- S_loss_P3[sample(1:nrow(S_loss_P1), rarify), ]
A_loss_P3 <- A_loss_P3[sample(1:nrow(A_loss_P1), rarify), ]
H_fill_P3 <- H_fill_P3[sample(1:nrow(H_fill_P1), rarify), ]


## calculate median from simulations
S_P1 <- as.data.frame(matrix(apply(S_loss_P1,2,median), ncol=1))
A_P1 <- as.data.frame(matrix(apply(A_loss_P1,2,median), ncol=1))
H_P1 <- as.data.frame(matrix(apply(H_fill_P1,2,median), ncol=1))

S_P2 <- as.data.frame(matrix(apply(S_loss_P2,2,median), ncol=1))
A_P2 <- as.data.frame(matrix(apply(A_loss_P2,2,median), ncol=1))
H_P2 <- as.data.frame(matrix(apply(H_fill_P2,2,median), ncol=1))

S_P3 <- as.data.frame(matrix(apply(S_loss_P3,2,median), ncol=1))
A_P3 <- as.data.frame(matrix(apply(A_loss_P3,2,median), ncol=1))
H_P3 <- as.data.frame(matrix(apply(H_fill_P3,2,median), ncol=1))


## calculate 95 CI
S_loss_P1 <- apply(S_loss_P1, 2, quantile, c(0.025, 0.975))
A_loss_P1 <- apply(A_loss_P1, 2, quantile, c(0.025, 0.975))
H_fill_P1 <- apply(H_fill_P1, 2, quantile, c(0.025, 0.975))


S_loss_P2 <- apply(S_loss_P2, 2, quantile, c(0.025, 0.975))
A_loss_P2 <- apply(A_loss_P2, 2, quantile, c(0.025, 0.975))
H_fill_P2 <- apply(H_fill_P2, 2, quantile, c(0.025, 0.975))


S_loss_P3 <- apply(S_loss_P3, 2, quantile, c(0.025, 0.975))
A_loss_P3 <- apply(A_loss_P3, 2, quantile, c(0.025, 0.975))
H_fill_P3 <- apply(H_fill_P3, 2, quantile, c(0.025, 0.975))




## bind df 
df <- as.data.frame(t(rbind(S_loss_P1,A_loss_P1,H_fill_P1,
                            S_loss_P2,A_loss_P2,H_fill_P2,
                            S_loss_P3,A_loss_P3,H_fill_P3,
                            S0,A0,H0)))

names(df)[1]="Smin_1"
names(df)[2]="Smax_1"
names(df)[3]="Amin_1"
names(df)[4]="Amax_1"
names(df)[5]="Hmin_1"
names(df)[6]="Hmax_1"
names(df)[7]="Smin_2"
names(df)[8]="Smax_2"
names(df)[9]="Amin_2"
names(df)[10]="Amax_2"
names(df)[11]="Hmin_2"
names(df)[12]="Hmax_2"
names(df)[13]="Smin_3"
names(df)[14]="Smax_3"
names(df)[15]="Amin_3"
names(df)[16]="Amax_3"
names(df)[17]="Hmin_3"
names(df)[18]="Hmax_3"


## bind data frames with median 
means <- as.data.frame(cbind(S_P1, A_P1, H_P1,
                             S_P2, A_P2, H_P2,
                             S_P3, A_P3, H_P3))

names(means)[1]="S_P1"
names(means)[2]="A_P1"
names(means)[3]="H_P1"
names(means)[4]="S_P2"
names(means)[5]="A_P2"
names(means)[6]="H_P2"
names(means)[7]="S_P3"
names(means)[8]="A_P3"
names(means)[9]="H_P3"
## END data configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## save and load ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## use to plot one: either High or Low
#singleSim <- cbind(df, means)


## use to plot both in final figure: low kelp (A = 30) simulation
combined_low <- cbind(df, means)
save(combined_low, file="ODE_low-kelp_toPlot.RDA")


## use to plot both in final figure: high kelp (A=300) simulation
combined_high <- cbind(df, means)
save(combined_high, file="ODE_high-kelp_toPlot.RDA")


## load simulated data to plot
setwd(simData)
load("combined_low.RDA")
load("combined_high.RDA")
singleSim <- combined_high
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





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
lowH <- "#CD9B9B"
highH <- "#9D1309"

alph <- 0.5
linew <- 1
lwd_sml <- 0.75
wid1 <- .5
wid2 <- 1
lineCol <- "black"
#ymax <- 6


both.blank <- theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
                    axis.title.x = element_blank(), axis.text.x = element_blank())
x.blank <- theme(axis.title.x = element_blank(), axis.text.x = element_blank())
y.blank <- theme(axis.title.y = element_blank(), axis.text.y = element_blank())


## set up custom ggplot theme 
my.theme = theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line = element_line(colour = "black"),
                 axis.title=element_text(size=16),
                 axis.text=element_text(size=14),
                 plot.title = element_text(size=16))
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Plot individual Low or High simulation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
windows(12,5,record=T)

p1 <- ggplot(data=singleSim, aes(S0)) +
  geom_ribbon(aes(ymin=Amin_1, ymax=Amax_1), fill=s.kelp, alpha=alph) +
  geom_ribbon(aes(ymin=Hmin_1, ymax=Hmax_1), fill=s.fullness, alpha=alph) +
  geom_ribbon(aes(ymin=Smin_1, ymax=Smax_1), fill=s.drift, alpha=alph) + my.theme +
  geom_line(aes(S0, H_P1), color=lineCol, size=wid2) +
  geom_line(aes(S0, H_P1), color=s.fullness, size=wid1) +
  geom_line(aes(S0, A_P1), color=lineCol, size=wid2) +
  geom_line(aes(S0, A_P1), color=s.kelp, size=wid1) +
  geom_line(aes(S0, S_P1), color=lineCol, size=wid2) +
  geom_line(aes(S0, S_P1), color=s.drift, size=wid1) +
  xlab("initial Drift") + ylab("state value") + ggtitle("period 1: t = 44hr snapshot") #+
  #ylim(0, ymax)

p2 <- ggplot(data=singleSim, aes(S0)) +
  geom_ribbon(aes(ymin=Amin_2, ymax=Amax_2), fill=s.kelp, alpha=alph) +
  geom_ribbon(aes(ymin=Hmin_2, ymax=Hmax_2), fill=s.fullness, alpha=alph) +
  geom_ribbon(aes(ymin=Smin_2, ymax=Smax_2), fill=s.drift, alpha=alph) + my.theme +
  geom_line(aes(S0, H_P2), color=lineCol, size=wid2) +
  geom_line(aes(S0, H_P2), color=s.fullness, size=wid1) +
  geom_line(aes(S0, A_P2), color=lineCol, size=wid2) +
  geom_line(aes(S0, A_P2), color=s.kelp, size=wid1) +
  geom_line(aes(S0, S_P2), color=lineCol, size=wid2) +
  geom_line(aes(S0, S_P2), color=s.drift, size=wid1) +
  xlab("initial Drift") + ylab("state value") + ggtitle("period 2: t = 89hr snapshot") #+
  #ylim(0, ymax) + y.blank

p3 <- ggplot(data=singleSim, aes(S0)) +
  geom_ribbon(aes(ymin=Amin_3, ymax=Amax_3), fill=s.kelp, alpha=alph) +
  geom_ribbon(aes(ymin=Hmin_3, ymax=Hmax_3), fill=s.fullness, alpha=alph) +
  geom_ribbon(aes(ymin=Smin_3, ymax=Smax_3), fill=s.drift, alpha=alph) + my.theme +
  geom_line(aes(S0, H_P3), color=lineCol, size=wid2) +
  geom_line(aes(S0, H_P3), color=s.fullness, size=wid1) +
  geom_line(aes(S0, A_P3), color=lineCol, size=wid2) +
  geom_line(aes(S0, A_P3), color=s.kelp, size=wid1) +
  geom_line(aes(S0, S_P3), color=lineCol, size=wid2) +
  geom_line(aes(S0, S_P3), color=s.drift, size=wid1) +
  xlab("initial Drift") + ylab("state value") + ggtitle("period 3: t = 134hr snapshot") #+
#  ylim(0, ymax) + y.blank


## plot 
all6 <- ggarrange(tag_facet(p1 + facet_wrap(~"time"), tag_pool="a"),
                  tag_facet(p2 + facet_wrap(~"time"), tag_pool="b"),
                  tag_facet(p3 + facet_wrap(~"time"), tag_pool="c"),
                  nrow=1, ncol=3)

print(all6)
## END ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Plot Low and High individually in final simulated figure ~~~~~~~~~~~~~~~~~~~~
## Drift loss 
p1 <- ggplot(data=combined_high, aes(S0)) +
  geom_ribbon(aes(ymin=Smin_1, ymax=Smax_1), fill=highD, alpha=alph) + my.theme +
  geom_line(aes(S0, S_P1), color=lineCol, size=wid2) +
  geom_line(aes(S0, S_P1), color=highD, size=wid1) +
  geom_ribbon(data=combined_low, aes(ymin=Smin_1, ymax=Smax_1), fill=lowD, alpha=alph) + my.theme +
  geom_line(data=combined_low, aes(S0, S_P1), color=lineCol, size=wid2) +
  geom_line(data=combined_low, aes(S0, S_P1), color=lowD, size=wid1) +
  xlab("initial Drift") + ylab("drift consumed") + ggtitle("period 1: t = 44hr snapshot") +
  #ylim(0, ymax) + 
  x.blank

p2 <- ggplot(data=combined_high, aes(S0)) +
  geom_ribbon(aes(ymin=Smin_2, ymax=Smax_2), fill=highD, alpha=alph) + my.theme +
  geom_line(aes(S0, S_P2), color=lineCol, size=wid2) +
  geom_line(aes(S0, S_P2), color=highD, size=wid1) +
  geom_ribbon(data=combined_low, aes(ymin=Smin_2, ymax=Smax_2), fill=lowD, alpha=alph) + my.theme +
  geom_line(data=combined_low, aes(S0, S_P2), color=lineCol, size=wid2) +
  geom_line(data=combined_low, aes(S0, S_P2), color=lowD, size=wid1) +
  xlab("initial Drift") + ggtitle("period 2: t = 89hr snapshot") +
  #ylim(0, ymax) + 
  both.blank

p3 <- ggplot(data=combined_high, aes(S0)) +
  geom_ribbon(aes(ymin=Smin_3, ymax=Smax_3), fill=highD, alpha=alph) + my.theme +
  geom_line(aes(S0, S_P3), color=lineCol, size=wid2) +
  geom_line(aes(S0, S_P3), color=highD, size=wid1) +
  geom_ribbon(data=combined_low, aes(ymin=Smin_3, ymax=Smax_3), fill=lowD, alpha=alph) + my.theme +
  geom_line(data=combined_low, aes(S0, S_P3), color=lineCol, size=wid2) +
  geom_line(data=combined_low, aes(S0, S_P3), color=lowD, size=wid1) +
  xlab("initial Drift") + ggtitle("period 3: t = 134hr snapshot") +
  #ylim(0, ymax) + 
  both.blank


## Kelp loss
p4 <- ggplot(data=combined_high, aes(S0)) +
  geom_ribbon(aes(ymin=Amin_1, ymax=Amax_1), fill=highK, alpha=alph) + my.theme +
  geom_line(aes(S0, A_P1), color=lineCol, size=wid2) +
  geom_line(aes(S0, A_P1), color=highK, size=wid1) +
  geom_ribbon(data=combined_low, aes(ymin=Amin_1, ymax=Amax_1), fill=lowK, alpha=alph) + my.theme +
  geom_line(data=combined_low, aes(S0, A_P1), color=lineCol, size=wid2) +
  geom_line(data=combined_low, aes(S0, A_P1), color=lowK, size=wid1) +
  xlab("initial Drift") + ylab("kelp consumed") +
  #ylim(0, ymax) + 
  x.blank 

p5 <- ggplot(data=combined_high, aes(S0)) +
  geom_ribbon(aes(ymin=Amin_2, ymax=Amax_2), fill=highK, alpha=alph) + my.theme +
  geom_line(aes(S0, A_P2), color=lineCol, size=wid2) +
  geom_line(aes(S0, A_P2), color=highK, size=wid1) +
  geom_ribbon(data=combined_low, aes(ymin=Amin_2, ymax=Amax_2), fill=lowK, alpha=alph) + my.theme +
  geom_line(data=combined_low, aes(S0, A_P2), color=lineCol, size=wid2) +
  geom_line(data=combined_low, aes(S0, A_P2), color=lowK, size=wid1) +
  xlab("initial Drift") + ylab("kelp consumed") +
  #ylim(0, ymax) + 
  both.blank

p6 <- ggplot(data=combined_high, aes(S0)) +
  geom_ribbon(aes(ymin=Amin_3, ymax=Amax_3), fill=highK, alpha=alph) + my.theme +
  geom_line(aes(S0, A_P3), color=lineCol, size=wid2) +
  geom_line(aes(S0, A_P3), color=highK, size=wid1) +
  geom_ribbon(data=combined_low, aes(ymin=Amin_3, ymax=Amax_3), fill=lowK, alpha=alph) + my.theme +
  geom_line(data=combined_low, aes(S0, A_P3), color=lineCol, size=wid2) +
  geom_line(data=combined_low, aes(S0, A_P3), color=lowK, size=wid1) +
  xlab("initial Drift") + ylab("kelp consumed")  +
  #ylim(0, ymax) + 
  both.blank


## Gut fullness 
p7 <- ggplot(data=combined_high, aes(S0)) +
  geom_ribbon(aes(ymin=Hmin_1, ymax=Hmax_1), fill=highH, alpha=alph) + my.theme +
  geom_line(aes(S0, H_P1), color=lineCol, size=wid2) +
  geom_line(aes(S0, H_P1), color=highH, size=wid1) +
  geom_ribbon(data=combined_low, aes(ymin=Hmin_1, ymax=Hmax_1), fill=lowH, alpha=alph) + my.theme +
  geom_line(data=combined_low, aes(S0, H_P1), color=lineCol, size=wid2) +
  geom_line(data=combined_low, aes(S0, H_P1), color=lowH, size=wid1) +
  xlab("initial Drift") + ylab("cumulative stomach fullness") 
  #ylim(0, ymax)

p8 <- ggplot(data=combined_high, aes(S0)) +
  geom_ribbon(aes(ymin=Hmin_2, ymax=Hmax_2), fill=highH, alpha=alph) + my.theme +
  geom_line(aes(S0, H_P2), color=lineCol, size=wid2) +
  geom_line(aes(S0, H_P2), color=highH, size=wid1) +
  geom_ribbon(data=combined_low, aes(ymin=Hmin_2, ymax=Hmax_2), fill=lowH, alpha=alph) + my.theme +
  geom_line(data=combined_low, aes(S0, H_P2), color=lineCol, size=wid2) +
  geom_line(data=combined_low, aes(S0, H_P2), color=lowH, size=wid1) +
  xlab("initial Drift") + ylab("drift consumed") +
  #ylim(0, ymax) + 
  y.blank

p9 <- ggplot(data=combined_high, aes(S0)) +
  geom_ribbon(aes(ymin=Hmin_3, ymax=Hmax_3), fill=highH, alpha=alph) + my.theme +
  geom_line(aes(S0, H_P3), color=lineCol, size=wid2) +
  geom_line(aes(S0, H_P3), color=highH, size=wid1) +
  geom_ribbon(data=combined_low, aes(ymin=Hmin_3, ymax=Hmax_3), fill=lowH, alpha=alph) + my.theme +
  geom_line(data=combined_low, aes(S0, H_P3), color=lineCol, size=wid2) +
  geom_line(data=combined_low, aes(S0, H_P3), color=lowH, size=wid1) +
  xlab("initial Drift") + ylab("drift consumed") +
  #ylim(0, ymax) + 
  y.blank


## plot 
windows(12,12,record=T)

all9 <- ggarrange(tag_facet(p1 + facet_wrap(~"time"), tag_pool="a"),
                  tag_facet(p2 + facet_wrap(~"time"), tag_pool="b"),
                  tag_facet(p3 + facet_wrap(~"time"), tag_pool="c"),
                  tag_facet(p4 + facet_wrap(~"time"), tag_pool="d"),
                  tag_facet(p5 + facet_wrap(~"time"), tag_pool="e"),
                  tag_facet(p6 + facet_wrap(~"time"), tag_pool="f"),
                  tag_facet(p7 + facet_wrap(~"time"), tag_pool="d"),
                  tag_facet(p8 + facet_wrap(~"time"), tag_pool="e"),
                  tag_facet(p9 + facet_wrap(~"time"), tag_pool="f"),
                  nrow=3, ncol=3)

print(all9)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END OF SCRIPT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
