##############################################################
# DENSITIES AND LORENZ CURVES COMPARISON
# Exponential VS Pareto
# ------------------------------------------------------------
# Objective: Reproduce Figure 1
##############################################################

# Useful packages
library(ggplot2)
library(actuar)
library(ineq)

# Grid of evaluation points
x <- seq(0, 5, by = 0.01)

# Distribution parameters
rate_exp <- 1                 # Exponential parameter
alpha_par <- 2/3              # Pareto parameter (shape)
xmin_par  <- 1                # Pareto parameter (scale)

##############################################################
# DENSITIES
##############################################################

# Exponential density
dens_exp <- dexp(x, rate = rate_exp)

# Pareto density (type I)
dens_pareto <- dpareto1(x, shape = alpha_par, min = xmin_par)
dens_pareto[x < xmin_par] <- NA # Defined only for x >= xmin

# Data frame construction
df_dens <- data.frame(x = x,dens_exp = dens_exp,dens_pareto = dens_pareto)

# Plot
dens<-ggplot(df_dens, aes(x = x)) +
  geom_line(aes(y = dens_exp), color = "#30D5C8", linewidth = 1.2) +
  geom_line(aes(y = dens_pareto), color = "magenta", linewidth = 1.2) +
  # Support of the Pareto distribution
  geom_segment( x = 0, xend = xmin_par,y = 0, yend = 0,color = "magenta",linewidth = 1.2) +
  labs(x = "", y = "") +
  coord_cartesian(xlim = c(0, 5), ylim = c(0, 1)) +
  theme_minimal()

# Reproduce left panel of Figure 1
print(dens)

##############################################################
# LORENZ CURVES
##############################################################

# Grid of population shares
p <- seq(0, 1, by = 0.01)

# Exponential Lorenz curve
L_exp <- Lc.exp(p)
L_exp[length(L_exp)] <- 1  # Numerical correction at p = 1

# Pareto Lorenz curve
L_par <- 1 - (1 - p)^(1 - alpha_par) # Analytical expression

# Data frame contruction
df_lorenz <- data.frame(p = p,L_exp = L_exp,L_par = L_par)

# Perfect equality line
df_equal <- data.frame(p = c(0, 1),L = c(0, 1))

# Reference population share
p0 <- 0.9

# Interpolated Lorenz values at p0
L_exp_09 <- approx(p, L_exp, xout = p0)$y
L_par_09 <- approx(p, L_par, xout = p0)$y

# Plot
lor<-ggplot(df_lorenz, aes(x = p)) +
  # Perfect equality line
  geom_line(data = df_equal, aes(y = L),color = "#66B032",linetype = "dashed",linewidth = 0.8) +
  # Exponential Lorenz curve
  geom_line(aes(y = L_exp),color = "#30D5C8",linewidth = 1.2) +
  # Area between perfect equality line and Lorenz curve
  geom_ribbon(aes(ymin = L_exp, ymax = p),fill = "#30D5C8",alpha = 0.1) +
  # Pareto Lorenz curve
  geom_line(aes(y = L_par),color = "magenta",linewidth = 1.2) +
  # Area between perfect equality line and Lorenz curve
  geom_ribbon(aes(ymin = L_par, ymax = p),fill = "magenta",alpha = 0.1) +
  # Vertical reference segments at p0
  geom_segment(aes(x = p0, xend = p0,y = 0, yend = L_exp_09),color = "#30D5C8",linetype = "dotted",linewidth = 0.9) +
  geom_segment(aes(x = p0, xend = p0,y = 0, yend = L_par_09),color = "magenta",linetype = "dotted",linewidth = 0.9) +
  # Horizontal reference segments
  geom_segment(aes(x = 0, xend = p0,y = L_exp_09, yend = L_exp_09),color = "#30D5C8",linetype = "dotted",linewidth = 0.9) +
  geom_segment(aes(x = 0, xend = p0,y = L_par_09, yend = L_par_09),color = "magenta",linetype = "dotted",linewidth = 0.9) +
  labs(x = "", y = "") +
  coord_fixed() +
  theme_minimal()

# Reproduce right panel of Figure 1
print(lor)
