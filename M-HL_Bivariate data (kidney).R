###---------------------------------------------- Compare 3 methods in 3 cases ----------------------------------------------###
rm(list=ls())

library(frailtyHL)
library(frailtyEM)
library(joint.Cox)
library(numDeriv)
library(Copula.surv)
library(ggplot2)
library(dplyr)
setwd("C:\\directory")
source("source_code.R")

###------------------------------------- Bivariate survival data -------------------------------------###
data(kidney)
dataset <- kidney
G=length( unique(dataset$id) )

x.obs=y.obs=dx=dy=age=sex=NULL
for(i in 1:G){
	temp=dataset$id==i
	x.obs[i]=dataset[temp,"time"][1]
	y.obs[i]=dataset[temp,"time"][2]
	dx[i]=dataset[temp,"status"][1]
	dy[i]=dataset[temp,"status"][2]
	age[i]=dataset[temp,"age"][1]
	sex[i]=dataset[temp,"sex"][1]
}

P <- length(x.obs)
EVENT <-matrix(0, P, 4)

EVENT[,1] <- as.matrix(x.obs, P, 1) # t_event1
EVENT[,2] <- as.matrix( (y.obs), P, 1) # t_event2
EVENT[,3] <- as.matrix(dx, P, 1) # event1
EVENT[,4] <- as.matrix(dy, P, 1) # event2

t_event1 <- EVENT[,1]
t_event2 <- EVENT[,2]
event1 <- EVENT[,3]
event2 <- EVENT[,4]

t1max <- max(t_event1)
t1min <- min(t_event1)
t2max <- max(t_event2)
t2min <- min(t_event2)

formula <- ~ age + sex + (1| unique(dataset$id) )

X <- model.matrix(formula)
X <- data.frame(X)
X <- as.matrix(X)

Z <- cbind(age, sex)

G <- length(unique(dataset$id))
group <- unique(dataset$id)

#------------------------------------------- M-HL -------------------------------------------#
p <- dim(Z)[2]
init <- c(rep(0, 2), rep(0, 2), rep(0.0, 5), rep(0.0, 5), log(0.2), 1, 1) # beta1, beta2, spline1, spline2, eta, kappa1, kappa2

res11 <- comp_HL_gamma1(formula, EVENT, init, seed=333)
res11

#------------------------------------------- random effect figure -------------------------------------------#
v_h <- res11[[9]]
SE <- res11[[10]]

#lower (upper) bound
LB <- v_h-1.96*SE
UB <- v_h+1.96*SE 

CI <- cbind(LB, UB) #confidence interval (CI)

x <- 1:length(unique(dataset$id))

result <- cbind(x, v_h, CI)

bi_random_figure <- ggplot(result, aes(x = x, y = v_h) ) +
		geom_point(size = 1.5) +
		geom_errorbar(aes(ymin = LB, ymax = UB), width = 0.2) +
		geom_hline(yintercept = 0) +
		labs(x = "Patient number", y = "Estimated patient effects") +
		theme_minimal()+
		theme(
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			panel.border = element_rect(fill = NA, size = 1),
			panel.background = element_blank()
		)

#------------------------------------------- baseline hazard figure -------------------------------------------#
# Estimate
res_g1 <- res11$Est_g1
res_g2 <- res11$Est_g2

# Setting M-spline
time <- seq(0, t1max, 1)
Ms <- M.spline(time, 0, t1max)

# Plot
esti_g1 <- Ms%*%t(res_g1)
esti_g2 <- Ms%*%t(res_g2)

ylim <- c(0, 4)

df1 <- data.frame(Time = time, Hazard = esti_g1, Group = "First infection")
df2 <- data.frame(Time = time, Hazard = esti_g2, Group = "Second infection")
plot_data <- rbind(df1, df2)

bi_base_figure <- ggplot(plot_data, aes(x = Time, y = Hazard, color = Group)) +
	geom_line(size = 1) +
	scale_color_manual(values = c("First infection" = "black", "Second infection" = "red")) +
	scale_x_continuous(breaks = seq(0, max(plot_data$Time), by = 100), expand = c(0, 0)) + 
	scale_y_continuous(breaks = seq(0, ylim[2], by = 0.5), limits = ylim, expand = c(0, 0)) +
	labs(x = "Time (days)", y = "Baseline hazard") +
	theme_classic() +
	theme(
		legend.position = c(0.95, 0.95),
		legend.justification = c("right", "top"),
		legend.title = element_blank(),
		legend.background = element_blank()
	)
