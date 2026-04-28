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

###------------------------------------- Semi-competing risks survival data -------------------------------------###
data(dataOvarian1)
data(dataOvarian2)
dataset <- cbind(dataOvarian1[, c("t.event", "event", "group", "CXCL12")], dataOvarian2[, c("t.death", "death")])
dataset$id <- as.numeric(factor(dataset$group, levels = unique(dataset$group)))
attach(dataset)

P <- dim(dataset)[1]
EVENT <-matrix(0, P, 4)

EVENT[,1] <- as.matrix(t.event, P, 1) # t_event1
EVENT[,2] <- as.matrix(t.death, P, 1) # t_event2
EVENT[,3] <- as.matrix(event, P, 1) # event1
EVENT[,4] <- as.matrix(death, P, 1) # event2

t_event1 <- EVENT[,1]
t_event2 <- EVENT[,2]
event1 <- EVENT[,3]
event2 <- EVENT[,4]

t1max <- max(t_event1)
t1min <- min(t_event1)
t2max <- max(t_event2)
t2min <- min(t_event2)

formula <-  ~ CXCL12 + (1|id)

X <- model.matrix(formula)
X <- data.frame(X)
X <- as.matrix(X)
Z <- as.matrix(X[, c(2)])

G <- length(unique(id))
group <- id

temp1 <- dataOvarian1[, c("t.event", "event", "group", "CXCL12")]
temp2 <- dataOvarian2[, c("t.death", "death", "group", "CXCL12")]
temp2$death <- ifelse(temp2$death==1, 2, 0)

names(temp1)[names(temp1) == "t.event"] <- "time"
names(temp2)[names(temp2) == "t.death"] <- "time"
names(temp1)[names(temp1) == "event"] <- "status"
names(temp2)[names(temp2) == "death"] <- "status"

dataset <- rbind(temp1, temp2)
dataset$id <- as.numeric(factor(dataset$group, levels = unique(dataset$group)))

#----------------- M-HL -----------------#
init <- c(rep(0, 1), rep(0, 1), rep(0.0, 5), rep(0.0, 5), log(0.1), 1, 1) # beta1, beta2, spline1, spline2, eta, kappa1, kappa2

res21 <- comp_HL_gamma1(formula, EVENT, init, seed=777)
res21

#----- random effect figure -----#
v_h <- res21[[9]]
SE <- res21[[10]]

#lower (upper) bound
LB <- v_h-1.96*SE
UB <- v_h+1.96*SE 

CI <- cbind(LB, UB) #confidence interval (CI)

x <- 1:length(unique(dataset$id))

result <- cbind(x, v_h, CI)

semi_random_figure <- ggplot(result, aes(x = x, y = v_h) ) +
		geom_point(size = 1.5) +
		geom_errorbar(aes(ymin = LB, ymax = UB), width = 0.2) +
		geom_hline(yintercept = 0) +
		labs(x = "Group number", y = "Estimated group effects") +
		theme_minimal()+
		theme(
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			panel.border = element_rect(fill = NA, size = 1),
			panel.background = element_blank()
		)

#----- baseline hazard figure -----#
# Estimate
res_g1 <- res21$Est_g1
res_g2 <- res21$Est_g2

# Setting M-spline
time <- seq(0, t1max, 1)
Ms <- M.spline(time, 0, t1max)

# Plot
esti_g1 <- Ms%*%t(res_g1)
esti_g2 <- Ms%*%t(res_g2)

ylim <- c(0, 0.0016)

df1 <- data.frame(Time = time, Hazard = esti_g1, Group = "Recurrence")
df2 <- data.frame(Time = time, Hazard = esti_g2, Group = "Death")
plot_data <- rbind(df1, df2)
plot_data$Group <- factor(plot_data$Group, levels = c("Recurrence", "Death"))

semi_base_figure <- ggplot(plot_data, aes(x = Time, y = Hazard, color = Group)) +
	geom_line(size = 1) +
	scale_color_manual(values = c("Recurrence" = "black", "Death" = "red")) +
	scale_x_continuous(breaks = seq(0, max(plot_data$Time), by = 500), expand = c(0, 0)) + 
	scale_y_continuous(breaks = seq(0, ylim[2], by = 0.0002), limits = ylim, expand = c(0, 0)) +
	labs(x = "Time (days)", y = "Baseline hazard") +
	theme_classic() +
	theme(
		legend.position = c(0.95, 0.95),
		legend.justification = c("right", "top"),
		legend.title = element_blank(),
		legend.background = element_blank()
	)
