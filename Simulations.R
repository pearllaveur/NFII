############################################################
# ILLUSTRATION - SIMULATED DATA
# ----------------------------------------------------------
# Objective: Reproduce Figures 3 and 4
# We compare the empirical behaviour of four inequality 
# indices (Atkinson, Gini, Extended Gini, Bonferroni)
# across three distributions (Pareto, Beta, Weibull),
# calibrated to share the same theoretical inequality level
# I = 0.15.
#
# Simulation methodology:
# N = 1000 replications, sample size n = 1000.
############################################################

# Useful packages
library(ggplot2)
library(reshape2)
library(actuar)
library(pracma)
library(dplyr)
library(cubature)

n<-1000 # Sample size
N<-1000 # Number of replications
I<-0.15 # Target theoretical inequality level
set.seed(156) # Set the seed

#######################################################################
# GENERAL ESTIMATOR OF I_{kappa,G}
#######################################################################
# By selecting the parameter kappa and the weight
# function g (or its c.d.f. G), several classical indices are recovered.

Est<-function(ech, kappa, G){
  X=sort(ech) # Order statistics
  n=length(ech)
  i=1:n
  D=sum((G((n+1-i)/n)-G((n-i)/n)) * (X**kappa)) # Estimator of D_{kappa,G}
  return(1-(D**(1/kappa))/mean(ech))
}

############################################################
# GENERIC VIOLIN PLOT
############################################################
# Produces standardized violin plots to compare the empirical
# sampling distributions of the estimators across the three distributions.

plot_violin_estimator <- function(data,
                                  value_var,
                                  true_value,
                                  y_limits,
                                  x_expressions,
                                  evi_values) {
  
  data$Distribution <- factor(
    data$Distribution,
    levels = c("F1", "F2", "F3")
  )
  ymin <- y_limits[1]
  ymax <- y_limits[2]
  
  ggplot(data,
         aes(x = Distribution,
             y = .data[[value_var]],
             fill = Distribution)) +
    geom_violin(trim = TRUE,
                alpha = 0.7) +
    # Sample mean represented as black dot
    stat_summary(fun = mean,
                 geom = "point",
                 size = 2.5,
                 color = "black") +
    # True theoretical inequality level
    geom_hline(yintercept = true_value,
               linetype = "dashed",
               color = "black") +
    
    # Aesthetic framing
    geom_segment(aes(x = 0.5, xend = 0.5,
                     y = ymin, yend = ymax),
                 inherit.aes = FALSE,
                 linewidth = 0.4,
                 color = "black") +
    geom_segment(aes(x = 1.5, xend = 1.5,
                     y = ymin, yend = ymax),
                 inherit.aes = FALSE,
                 linewidth = 0.4,
                 color = "black") +
    geom_segment(aes(x = 2.5, xend = 2.5,
                     y = ymin, yend = ymax),
                 inherit.aes = FALSE,
                 linewidth = 0.4,
                 color = "black") +
    geom_segment(aes(x = 3.5, xend = 3.5,
                     y = ymin, yend = ymax),
                 inherit.aes = FALSE,
                 linewidth = 0.4,
                 color = "black") +
    geom_segment(aes(x = 0.5, xend = 3.5,
                     y = ymin, yend = ymin),
                 inherit.aes = FALSE,
                 linewidth = 0.4,
                 color = "black") +
    geom_segment(aes(x = 0.5, xend = 3.5,
                     y = ymax, yend = ymax),
                 inherit.aes = FALSE,
                 linewidth = 0.4,
                 color = "black") +
    scale_x_discrete(limits = c("F1","F2","F3"),
                     labels = x_expressions) +
    scale_fill_manual(values = c("blue", "darkgreen", "red")) +
    coord_flip(ylim = y_limits) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.y = element_text(size = 18),
      axis.text.x = element_text(size = 18)
    ) +
    
    labs(x = NULL,
         y = NULL) +
    # Display extreme-value indices
    annotate("text", x = 0.65, y = ymax-0.05, label = bquote(gamma == .(evi_values[1])), color = "blue", size =10) +
    annotate("text", x = 1.65, y = ymax-0.05, label = bquote(gamma == .(evi_values[2])), color = "darkgreen", size = 10) +
    annotate("text", x = 2.65, y = ymax-0.05, label = bquote(gamma == .(evi_values[3])), color = "red", size = 10)
}

############################################################
# 1) ATKINSON INDEX (eta = 0.9)
# Reproduce top right panel of Figure 3
############################################################

eta<-0.9

# Pareto parameter calibration (numerical root finding)
f1<-function(x){
  return(1-I-(1-x)*(1-(1-eta)*x)**(1/(eta-1)))
}
gam<-uniroot(f1,interval=c(0.05,0.95))$root 

# Beta parameter calibration
f2<-function(x){
  return(1-I-((x+1)*(x/(1-eta+x))**(1/(1-eta)))/x)
}
theta<-uniroot(f2,interval=c(0.1,10))$root

# Weibull parameter calibration
f3<-function(x){
  return(1-I-(gamma(1+(1-eta)/x)**(1/(1-eta)))/(gamma(1+1/x)))
}
b<-uniroot(f3,interval=c(0.1,10))$root

# Simulations
Pareto_samples <- as.matrix(replicate(N, rpareto1(n, shape = 1/gam, min = 1))) #F1
Beta_samples <- as.matrix(replicate(N, rbeta(n, shape1 = theta, shape2 = 1))) #F2
Weibull_samples <- as.matrix(replicate(N, rweibull(n, shape = b, scale = 1))) #F3

# Estimations
# Atkinson index corresponds to kappa = 1 - eta and weight function g(u) = 1 i.e. G(u) = u.
kappa<-1-eta
G<-function(x){
  return(x)
}
EstAt<-function(ech){
  return(Est(ech, kappa, G))
}
At_chap1<-apply(Pareto_samples,1,EstAt)
At_chap2<-apply(Beta_samples,1,EstAt)
At_chap3<-apply(Weibull_samples,1,EstAt)

# Data preparation for visualization
data_at <- data.frame(
  value = c(At_chap1, At_chap2, At_chap3),
  Distribution = rep(c("F1","F2","F3"),
                     each = length(At_chap1))
)

# Produce violin plots
p_at <- plot_violin_estimator(
  data = data_at,
  value_var = "value",
  true_value = I,
  y_limits = c(0.05, max(At_chap1)+0.01),
  x_expressions = c(
    expression(hat(I)[Atkinson*","*0.9](X[gamma])),
    expression(hat(I)[Atkinson*","*0.9](Y[theta])),
    expression(hat(I)[Atkinson*","*0.9](Z[b]))
  ),
  evi_values = c(round(gam,2),-1,0)
)

# Reproduce top right panel of Figure 3
ggsave("violin-plot-Atkinson.pdf", p_at,width=12, height=8)

# ------------------------------------------------------------------
# Asymptotic distribution of Atkinson estimator
# Reproduce Figure 4: illustration of Theorem 3 (i)
# ------------------------------------------------------------------

# Generic asymptotic variance formula
# f(param, t) = 1/E[X^t], param = distribution parameter
VarAsymp<-function(gam,k,f){
  return(f(gam,1)**2/f(gam,2) - 2*f(gam,k)*f(gam,1)/(k*f(gam,k+1)) + f(gam,k)**2/(k**2*f(gam,2*k)) - (1-1/k)**2)
}

# Under Pareto distribution
# Moment condition satisfied: gamma < 1/2
# Inverse moment function: f(gamma, t) = 1 - gamma * t
fP<-function(gam,t){
  return(1-gam*t)
}
VP_At<-VarAsymp(gam,1-eta,fP)

# Under Beta distribution
# Moment condition: kappa > -theta / 2
# Inverse moment function: f(theta, t) = 1 + t / theta
fB<-function(theta,t){
  return(1+t/theta)
}
VB_At<-VarAsymp(theta,1-eta,fB)

# Under Weibull distribution
# Moment condition: kappa > -b / 2
# Inverse moment function: f(b, t) = 1 / Gamma(1 + t / b)
fW<-function(b,t){
  return(1/gamma(1+t/b))
}
VW_At<-VarAsymp(b,1-eta,fW)

# Data frame construction of limiting Gaussian densities
x <- seq(-10, 10, length.out = 1000)
df <- bind_rows(
  data.frame(x = x, y = dnorm(x, mean = 0, sd = sqrt(VP_At)),Var="VP"),
  data.frame(x = x, y = dnorm(x, mean = 0, sd = sqrt(VB_At)),Var="VB"),
  data.frame(x = x, y = dnorm(x, mean = 0, sd = sqrt(VW_At)),Var="VW")
)

# Plot of limiting Gaussian densities
figVar<-ggplot(df, aes(x = x, y = y, color = Var)) +
  geom_line(linewidth = 0.75) +
  geom_segment(aes(x = -7, y = -0.1, xend = 7, yend = -0.1), color = "black", linewidth = 0.4) +
  geom_segment(aes(x = -7, y = 1.75, xend = 7, yend = 1.75), color = "black", linewidth = 0.4) +  
  geom_segment(aes(x = -7, y = -0.1, xend = -7, yend = 1.75), color = "black", linewidth = 0.4) +  
  geom_segment(aes(x = 7, y = -0.1, xend = 7, yend = 1.75), color = "black", linewidth = 0.4)  +
  scale_color_manual(values=c("VP"="blue","VB"="darkgreen","VW"="red")) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_minimal() +
  xlim(-7,7)+
  theme(legend.position = "none",
        axis.text.y = element_text(size = 25),
        axis.text.x = element_text(size = 25))

# Reproduce Figure 4
ggsave("A-N-Atkinson-3lois.pdf", figVar,width=12, height=8)

############################################################
# 2) GINI INDEX
# Reproduce top left panel of Figure 3
############################################################

# Pareto parameter calibration (closed-form relationship)
gam<-(2*I)/(1+I)

# Beta parameter calibration
theta<-(1-I)/(2*I)

# Weibull parameter calibration
b<- -log(2)/log(1-I)

# Simulations
Pareto_samples <- as.matrix(replicate(N, rpareto1(n, shape = 1/gam, min = 1))) #F1
Beta_samples <- as.matrix(replicate(N, rbeta(n, shape1 = theta, shape2 = 1))) #F2
Weibull_samples <- as.matrix(replicate(N, rweibull(n, shape = b, scale = 1))) #F3

# Estimations
# Gini index corresponds to kappa = 1 and weight function g(u) = 2u i.e. G(u) = u^2.
kappa<-1
G<-function(x){
  return(x^2)
}
EstGini<-function(ech){ # Coincides with the "Gini" function in the "DescTools" package
  return(Est(ech, kappa, G))
}

Gini_chap1<-apply(Pareto_samples,1,EstGini)
Gini_chap2<-apply(Beta_samples,1,EstGini)
Gini_chap3<-apply(Weibull_samples,1,EstGini)

# Data preparation for visualization
data_gini <- data.frame(
  value = c(Gini_chap1, Gini_chap2, Gini_chap3),
  Distribution = rep(c("F1","F2","F3"),
                     each = length(Gini_chap1))
)

# Produce violin plots
p_gini <- plot_violin_estimator(
  data = data_gini,
  value_var = "value",
  true_value = I,
  y_limits = c(0.05, max(At_chap1)+0.01),
  x_expressions = c(
    expression(hat(I)[Gini](X[gamma])),
    expression(hat(I)[Gini](Y[theta])),
    expression(hat(I)[Gini](Z[b]))
  ),
  evi_values=c(round(gam,2),-1,0)
)

# Reproduce top left panel of Figure 3
ggsave("violin-plot-Gini.pdf", p_gini,width=12, height=8)

############################################################
# 3) EXTENDED GINI INDEX (nu = 3)
# Reproduce bottom left panel of Figure 3
############################################################
nu<-3

# Pareto parameter calibration (numerical root finding)
f1<-function(x){
  return((nu-1)*x/(nu-x) - I)
}
gam<-uniroot(f1,interval=c(0.05,0.95))$root

# Beta parameter calibration
f2<-function(x){
  return(1-I-nu*(x+1)*beta(nu,1/x)/(x*(nu*x+1)))
}
theta<-uniroot(f2,interval=c(0.1,10))$root

# Weibull parameter calibration
f3<-function(x){
  return(1-I-nu**(-1/x))
}
b<-uniroot(f3,interval=c(0.1,10))$root 

# Simulations
Pareto_samples <- as.matrix(replicate(N, rpareto1(n, shape = 1/gam, min = 1))) #F1
Beta_samples <- as.matrix(replicate(N, rbeta(n, shape1 = theta, shape2 = 1))) #F2
Weibull_samples <- as.matrix(replicate(N, rweibull(n, shape = b, scale = 1))) #F3

# Estimations
# Extended Gini index corresponds to kappa = 1 and weight function g(u) = (nu)*u^{nu-1} i.e. G(u) = u^(nu).
kappa<-1
G<-function(x){
  return(x**(nu))
}

EstEG<-function(ech){
  return(Est(ech, kappa, G))
}
EG_chap1<-apply(Pareto_samples,1,EstEG)
EG_chap2<-apply(Beta_samples,1,EstEG)
EG_chap3<-apply(Weibull_samples,1,EstEG)

# Data preparation for visualization
data_eg <- data.frame(
  value = c(EG_chap1, EG_chap2, EG_chap3),
  Distribution = rep(c("F1","F2","F3"),
                     each = length(EG_chap1))
)

# Produce violin plots
p_eg <- plot_violin_estimator(
  data = data_eg,
  value_var = "value",
  true_value = I,
  y_limits = c(0.05, max(At_chap1)+0.01),
  x_expressions = c(
    expression(hat(I)[EG*","*3](X[gamma])),
    expression(hat(I)[EG*","*3](Y[theta])),
    expression(hat(I)[EG*","*3](Z[b]))
  ),
  evi_values=c(round(gam,2),-1,0)
)

# Reproduce bottom left panel of Figure 3
ggsave("violin-plot-EG.pdf", p_eg,width=12, height=8)

############################################################
# 4) BONFERRONI INDEX (classical one: h = 1)
# Reproduce bottom right panel of Figure 3
############################################################

# Pareto parameter calibration (numerical finding)
f <- function(y, x) {
  return(y*exp(-y)*(1-exp(-y))**(-x))
}
evaluate_equation <- function(x) {
  integral_value <- integral(function(y) f(y, x), 0, Inf) 
  result <- 1 - I-(1 - x) * integral_value
  return(abs(result))
}
g_values <- seq(0.1, 0.3, by = 0.001)
errors <- sapply(g_values, evaluate_equation)
gam <- g_values[which.min(errors)]

# Beta parameter calibration (closed-form relationship)
theta<-1/I -1

# Weibull parameter calibration (numerical finding)
f <- function(y, x) {
  return(log(1/(1-exp(-y)))*exp(-y)*y**(1/x))
}
evaluate_equation <- function(x) {
  integral_value <- integral(function(y) f(y, x), 0, Inf) 
  result <- 1 - I- integral_value/gamma(1+1/x)
  return(abs(result))
}
b_values <- seq(6, 9, by = 0.001)
errors <- sapply(b_values, evaluate_equation)
b <- b_values[which.min(errors)]

# Simulations
Pareto_samples <- as.matrix(replicate(N, rpareto1(n, shape = 1/gam, min = 1))) #F1
Beta_samples <- as.matrix(replicate(N, rbeta(n, shape1 = theta, shape2 = 1))) #F2
Weibull_samples <- as.matrix(replicate(N, rweibull(n, shape = b, scale = 1))) #F3

# Estimations
# Bonferroni index corresponds to kappa = 1 and weight function g(u) = -log(1-u) i.e. G(u) = u + (1 - u) * log(1 - u).
kappa<-1
G <- function(x){
  ifelse(x == 1, 1, x + (1 - x) * log(1 - x))
}

EstBon<-function(ech){
  return(Est(ech, kappa, G))
}
Bon_chap1<-apply(Pareto_samples,1,EstBon)
Bon_chap2<-apply(Beta_samples,1,EstBon)
Bon_chap3<-apply(Weibull_samples,1,EstBon)

# Data preparation for visualization
data_bon <- data.frame(
  value = c(Bon_chap1, Bon_chap2, Bon_chap3),
  Distribution = rep(c("F1","F2","F3"),
                     each = length(Bon_chap1))
)

# Produce violin plots
p_bon <- plot_violin_estimator(
  data = data_bon,
  value_var = "value",
  true_value = I,
  y_limits = c(0.05, max(At_chap1)+0.01),
  x_expressions = c(
    expression(hat(I)[Bonferroni*","*1](X[gamma])),
    expression(hat(I)[Bonferroni*","*1](Y[theta])),
    expression(hat(I)[Bonferroni*","*1](Z[b]))
  ),
  evi_values=c(round(gam,2),-1,0)
)

# Reproduce bottom right panel of Figure 3
ggsave("violin-plot-Bonferroni.pdf", p_bon,width=12, height=8)
