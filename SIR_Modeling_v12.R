#Adding R libraries - do install previously if needed

library("ggplot2")
library("RColorBrewer")
library("deSolve")
library("reshape2")
library("gridExtra")

# Set initial parameters for the population that are used for all plots
  # Proportion in each compartment:
  # Susceptible (S) 0.999998
  # Infected (I) 0.000002
  # Recovered (R) 0
    # For the SIR model with treatment, the starting population of I is transferred to:
    # Treated (Tr) 0.000001
    # Untreated (Ut) 0.000001
init1<- c(S = 1-1e-6, I = 2e-6, R = 0.0)
init2<- c(S = 1-1e-6, Tr = 1e-6, Ut = 1e-6, R = 0.0)
  # Time frame
times <- seq(0, 50, by = 0.01)

# Parameters for transmission rates:
  # beta: infection parameter
  # gamma: recovery parameter
  # alpha: treatment parameter
  # zeta: effect of treatment to infectivity parameter
  # eff: efficacy of treatment on recovery parameter

# Create an SIR function
sir <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    dS <- -beta * S * I
    dI <- beta * S * I - gamma * I
    dR <- gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}

# Parameters for SIR:
parameters_basic <- c(beta = 1.2,
                      gamma = 0.1)

# Solve using ode (General Solver for Ordinary Differential Equations)
out1 <- ode(y = init1, times = times, func = sir, parms = parameters_basic)
# change to data frame
out1 <- as.data.frame(out1)

head(out1)

data.sir <- melt(out1,id="time",measure=c("S","I","R"))

names(data.sir) = c("Time","Compartment","Value")

plot1 <- ggplot(data.sir)+
  geom_line(aes(x = Time, y = Value, color = Compartment), linewidth = 1.2) +
  theme_minimal() +
  xlab("Time") +
  ylab("Proportion of Population") +
  theme_classic() +
  theme(text = element_text(size = 15)) +
  ylim(0, 1) +
  scale_color_manual(values = c("#808080", "#FFAE42", "#4aa5d9")) +
  theme(plot.title = element_text(size = 12, face = "bold"))

# Create an SIR function with treatment
sir_t <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    # Ut+Tr are the same as I within SIR
    # From the whole Ut+Tr transmission rate, the effect of the Tr group gets subtracted
    dS <- -(beta * S * (Ut+Tr) - beta * (1-zeta) * S * Tr)
    dTr <- alpha* (beta * S * (Ut+Tr) - beta * (1-zeta) * S * Tr) - eff * gamma * Tr
    dUt <- (1-alpha) * (beta * S * (Ut+Tr) - beta * (1-zeta) * S * Tr) - gamma * Ut
    dR <- gamma * Ut + eff * gamma * Tr
    
    return(list(c(dS, dTr, dUt, dR)))
  })
}

# Parameters used for SIR with treatment:
#SIR with treatment, 40% get treatment, alpha = 0.4
#SIR with treatment - no effect of the drug - introduction of Tr and Ut
parameters_2 <- c(parameters_basic, alpha = 0.4, zeta = 1, eff = 1)
#SIR with treatment - high zeta  , low eff - low effect on infectivity, low efficacy on treatment
parameters_3 <- c(parameters_basic, alpha = 0.4, zeta = 0.8, eff = 1.25)
#SIR with treatment - high zeta  , high eff - low effect on infectivity, high efficacy on treatment
parameters_4 <- c(parameters_basic, alpha = 0.4, zeta = 0.8, eff = 3)
#SIR with treatment - low zeta , low eff - high effect on infectivity, low efficacy on treatment
parameters_5 <- c(parameters_basic, alpha = 0.4, zeta = 0.33, eff = 1.25)
#SIR with treatment - low zeta , high eff - high effect on infectivity, high efficacy on treatment
parameters_6 <- c(parameters_basic, alpha = 0.4, zeta = 0.33, eff = 3)

parameter_treatment <- list(parameters_2,parameters_3,parameters_4,parameters_5,parameters_6)
plot_treatment <- list()

#loop for different settings/parameters
for (i in 1:5){
  parameters <- get(paste0("parameters_",i+1))
  # Solve using ode (General Solver for Ordinary Differential Equations)
  out <- ode(y = init2, times = times, func = sir_t, parms = parameters)
  # change to data frame
  out <- as.data.frame(out)

  head(out)

  data.sir_t <- melt(out,id="time",measure=c("S","Tr","Ut","R"))

  names(data.sir_t) = c("Time","Compartment","Value")

  plot_treatment[[i]] <- ggplot(data.sir_t)+
    geom_line(aes(x = Time, y = Value, color = Compartment), linewidth = 1.2) +
    theme_minimal() +
    xlab("Time") +
    ylab("Proportion of Population") +
    theme_classic() +
    theme(text = element_text(size = 15)) +
    ylim(0, 1) +
    scale_color_manual(values = c("#808080","#FF10F0","#006942", "#4aa5d9")) +
    theme(plot.title = element_text(size = 12, face = "bold"))
}
#adding the title of the plots
plot1 <- plot1 + ggtitle("Basic SIR model")+theme(plot.title = element_text(hjust = 0.5, size = 12))
plot_treatment[[1]] <- plot_treatment[[1]]+ ggtitle("SIR model with treatment - no effect of the drug")+theme(plot.title = element_text(hjust = 0.5, size = 12))
plot_treatment[[2]] <- plot_treatment[[2]]+ ggtitle("SIR model with treatment - low effect on infectivity, low efficacy on treatment")+theme(plot.title = element_text(hjust = 0.5, size = 12))
plot_treatment[[3]] <- plot_treatment[[3]]+ ggtitle("SIR model with treatment - low effect on infectivity, high efficacy on treatment")+theme(plot.title = element_text(hjust = 0.5, size = 12))
plot_treatment[[4]] <- plot_treatment[[4]]+ ggtitle("SIR model with treatment - high effect on infectivity, low efficacy on treatment")+theme(plot.title = element_text(hjust = 0.5, size = 12))
plot_treatment[[5]] <- plot_treatment[[5]]+ ggtitle("SIR model with treatment - high effect on infectivity, high efficacy on treatment")+theme(plot.title = element_text(hjust = 0.5, size = 12))

# all plots in one grid
grid.arrange(plot1,
             plot_treatment[[1]],
             plot_treatment[[2]],
             plot_treatment[[3]],
             plot_treatment[[4]],
             plot_treatment[[5]],
             ncol=3,nrow=2)

# To access the single plots for proper human readability remove the #:
#plot1
#plot_treatment[[1]]
#plot_treatment[[2]]
#plot_treatment[[3]]
#plot_treatment[[4]]
#plot_treatment[[5]]