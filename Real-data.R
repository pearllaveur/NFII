################################################################################
# ILLUSTRATION - INCOME DATA (NMES 1988)
# ------------------------------------------------------------------------------
# Objective: Reproduce Figure 5 and recover numerical results in Subsection 6.2
#
# We estimate the Gini index of the data distribution,
# and the extreme-value index of the data distribution.
# We construct two distributions (densities) with 
# - same Gini index as the data set,
# - different extreme-value indices.
# We compute the grouped data Gini estimate (see (4.6)) of the data
# and compare it with the individual data Gini estimate
# using an asymptotic confidence interval
################################################################################

# Useful packages
library(evir)
library(evmix)
library(AER)
library(ReIns)
library(DescTools)
library(ggplot2)
library(actuar)

# Data preparation
data(NMES1988)
X=NMES1988$income[NMES1988$income>0]
n=length(X)

# Gini estimate
Gini_hat=Gini(X)

################################################################################
# HISTOGRAM
# Reproduce top left panel of Figure 5
################################################################################

hist<-ggplot(data.frame(X), aes(x = X,y=after_stat(density))) +
  # Aesthetic framing
  geom_segment(aes(x = -1.2, y = 0, xend = 57, yend = 0), color = "black", linewidth = 0.4) +
  geom_segment(aes(x = -1.2, y = 0.35, xend = 57, yend = 0.35), color = "black", linewidth = 0.4) +
  geom_segment(aes(x = -1.2, y = 0, xend = -1.2, yend = 0.35), color = "black", linewidth = 0.4) +
  geom_segment(aes(x = 57, y = 0, xend = 57, yend = 0.35), color = "black", linewidth = 0.4)  +
  geom_histogram(bins = 50,fill = "blue", color = "black",linewidth=0.5)+
  labs(x = "", y = "", title = "")+
  xlim(-1.2,57)+
  theme_minimal()+
  theme(axis.text = element_text(size=20))

# Reproduce top left panel of Figure 5
print(hist)

################################################################################
# Extreme-value index estimation using Hill estimator with optimal threshold
################################################################################

HPopt=Hill.kopt(X)
kopt=HPopt$kopt # Optimal number of excesses
gammaopt=HPopt$gammaopt # Estimated extreme-value index of the data distribution

################################################################################
# HILLPLOT
# Reproduce top right panel of Figure 5
################################################################################

hillplot_ggplot <- function(X,k_max){ # Plot up to kmax excesses
    HP=hillplot(X)
    df <- data.frame(k = HP$ks[2:k_max], Hill = HP$H[2:k_max],LCI = HP$cil.H[2:k_max],UCI = HP$ciu.H[2:k_max])
    ggplot(df, aes(x = k)) +
    geom_line(aes(y = Hill), color = "blue", linewidth = 1) +
    # Aesthetic framing
    geom_segment(aes(x = 0, y = 0.25, xend = 1125, yend = 0.25), color = "black", linewidth = 0.4) +  
    geom_segment(aes(x = 0, y = 0.55, xend = 1125, yend = 0.55), color = "black", linewidth = 0.4) +  
    geom_segment(aes(x = 0, y = 0.25, xend = 0, yend = 0.55), color = "black", linewidth = 0.4) +  
    geom_segment(aes(x = 1125, y = 0.25, xend = 1125, yend = 0.55), color = "black", linewidth = 0.4)  +
    annotate("point", x = kopt, y = gammaopt, shape = 17, size = 5, color = "red", fill = "red")+
    labs(x = "", y = "", title = "") +
    ylim(0.25,0.55)+
    theme_minimal()+
    theme(axis.text = element_text(size=20))
}
# Reproduce top right panel of Figure 5
print(hillplot_ggplot(X,floor(n/4))) # Plot up to n/4 excesses (25% of the data)

################################################################################
# PARETO-QUANTILE PLOT
# Verification of the extreme-value index estimate gammaopt
# Reproduce bottom left panel of Figure 5
################################################################################

Xsort=sort(X)
x_vec <- log((kopt+1) / (1:kopt))
y_vec <- log(Xsort[n - (1:kopt) + 1]) - rep(log(Xsort[n - kopt]), kopt)
line_y <- gammaopt * x_vec
df <- data.frame(x = x_vec, y = y_vec, line_y = line_y)

pqplot<-ggplot(df, aes(x = x, y = y)) +
  geom_point() +
  geom_line(aes(y = line_y),color="blue",lwd=1.1) +
  # Aesthetic framing
  geom_segment(aes(x = 0, y = 0, xend = 6, yend = 0), color = "black", linewidth = 0.4) +  
  geom_segment(aes(x = 0, y = 2.5, xend = 6, yend = 2.5), color = "black", linewidth = 0.4) +  
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = 2.5), color = "black", linewidth = 0.4) +  
  geom_segment(aes(x = 6, y = 0, xend = 6, yend = 2.5), color = "black", linewidth = 0.4)  +
  labs(x = "", y = "", title = "") +
  theme_minimal()+
  theme(axis.text = element_text(size=20))

# Reproduce bottom left panel of Figure 5
print(pqplot)

################################################################################
# CONSTRUCTION OF THE TWO DENSITIES
# Mixtures considered: p is the proportion associated with the Beta(theta,1) 
# distribution and 1-p the one associated with the Pareto(1,gamma) distribution
# Reproduce bottom right panel of Figure 5
################################################################################

# Quadratic equation in p derived from the Gini index expression of the mixture
A=function(x,y){
  return(2*x**2/((x+1)*(2*x+1)) - 2*x/(x+1) + 2/(2-y))
}
B=function(x,y){
  return(2*x/(x+1) - 4/(2-y) - x*(1-Gini_hat)/(x+1) + (1-Gini_hat)/(1-y))
}
C=function(y){
  return(2/(2-y) - (1-Gini_hat)/(1-y))
}
solve_Gini_mixt_pstar <- function(theta, gamma) {
  pstar <- c(NA, NA)
  discriminant <- B(theta, gamma)**2 - 4 * A(theta, gamma) * C(gamma)
  if (discriminant < 0) {
    return("No real solution")
  }
  p1 <- (-B(theta, gamma) - sqrt(discriminant)) / (2 * A(theta, gamma))
  p2 <- (-B(theta, gamma) + sqrt(discriminant)) / (2 * A(theta, gamma))
  pstar <- c(p1, p2)
  return(pstar[pstar > 0 & pstar < 1])
}

# Couples (theta,p) associated with gamma1 and gamma2, different from gammaopt
gamma1<-0.3
# For gamma1, thetastar1 is determined by plotting:
# p=solve_Gini_mixt_pstar(theta,gamma1), as a function of theta.
thetastar1=0.6 
pstar1=solve_Gini_mixt_pstar(thetastar1,gamma1)[1]

gamma2<-0.5
# Idem
thetastar2<-1 
pstar2<-solve_Gini_mixt_pstar(thetastar2,gamma2)[1]

# Density functions of both mixtures
f1<-function(x){ 
  return(pstar1*dbeta(x,shape1 = thetastar1, shape2 = 1)+(1-pstar1)*dpareto1(x,shape=1/gamma1,min=1,log=FALSE))
}
f2<-function(x){ 
  return(pstar2*dbeta(x,shape1 = thetastar2, shape2 = 1)+(1-pstar2)*dpareto1(x,shape=1/gamma2,min=1,log=FALSE))
}

# Plot comparison of both mixture densities
range<-seq(1, 5, length.out = 10000)
densities<-data.frame(
  f1 = f1(range),
  f2 = f2(range)
)
mixt_dens<-ggplot() +
  geom_line(data = densities, aes(x = range, y = f1),
            color = "blue", linewidth = 1) +
  geom_line(data = densities, aes(x = range, y = f2),
            color = "red", linewidth = 1) +
  theme_minimal() +
  labs(x = "", y = "", title = "") +
  ylim(0,1) +
  # Aesthetic framing
  geom_segment(aes(x = 1, y = 0, xend = 5, yend = 0), color = "black", linewidth = 0.4) +  
  geom_segment(aes(x = 1, y = 1, xend = 5, yend = 1), color = "black", linewidth = 0.4) +  
  geom_segment(aes(x = 1, y = 0, xend = 1, yend = 1), color = "black", linewidth = 0.4) +  
  geom_segment(aes(x = 5, y = 0, xend = 5, yend = 1), color = "black", linewidth = 0.4) +
  theme(axis.text = element_text(size=20))

# Reproduce bottom right panel of Figure 5
print(mixt_dens)

################################################################################
# GROUPED DATA ESTIMATION OF THE GINI INDEX OF THE INCOME DATA
# Recover numerical results in Subsection 6.2
################################################################################

# Extract histogram class boundaries
classes <- ggplot_build(hist)$data[[5]][,"xmin"] 
a <- c(0,classes[3:length(classes)]) # Start the 1st class at a_1=0
m <- length(a)-1 # Number of classes

mu_groups <- numeric(m)
p_hat <- numeric(m)

for(j in 1:m){
  mu_groups[j]=(a[j]+a[j+1])/2 # Theoretical expectation of group j distribution
  p_hat[j]=sum((X>a[j]) & (X<a[j+1]))/n # Group j frequency
}
mu_hat=sum(mu_groups*p_hat) # Vector of the j group means

weights=numeric(m)
for(j in 1:m){
  if(j<m){
    S=sum(p_hat[(j+1):m])
  }
  else{S=0}
  weights[j]=3*(a[j+1]+a[j])*S+(a[j+1]+2*a[j])*p_hat[j]
}

# Compute the grouped data Gini estimate
Gini_tilde <- 1-(sum(p_hat*weights))/(3*mu_hat)
print(Gini_tilde)

# Recall the individual data Gini estimate 
print(Gini_hat)

# Estimated asymptotic variance (4.4) of the Gini index on a generic sample
Est_VA_Gini<-function(sample){
  n=length(sample)
  muhat=mean(sample)
  X=sort(sample)
  Dhat=sum( X[n:1] * (2 * (1:n) - 1) ) / n^2
  J=function(u){
    return(1/muhat - 2*(1-u)/Dhat)
  }
  somme=0
  for(k in 1:(n-1)){
    Sk=0
    if(k>1){
      for(l in (1:(k-1))){
        Sk=Sk + J(l/n)*(l/n)*(X[l+1]-X[l])*(X[k+1]-X[k])
      }
    }
    somme=somme + (1-k/n)*J(k/n)*(2*Sk + J(k/n)*(k/n)*(X[k+1]-X[k])**2)
  }
  return(somme)
}

# Estimated asymptotic variance of the Gini index for the income data X
Vhat<-Est_VA_Gini(X)

# Asymptotic confidence interval for the Gini index (Theorem 3(ii))
IC_Gini<-function(alpha){ # alpha: level
  c((sqrt(n)*Gini_hat-qnorm(1-alpha/2)*sqrt(Vhat))/(sqrt(n)-qnorm(1-alpha/2)*sqrt(Vhat)),(sqrt(n)*Gini_hat+qnorm(1-alpha/2)*sqrt(Vhat))/(sqrt(n)+qnorm(1-alpha/2)*sqrt(Vhat)))}

# Compute the 90% level asymptotic confidence interval for the Gini index
print(IC_Gini(0.1))

# Compute the 90% level bootstrap confidence interval for the Gini index
IC_bootstrap <-Gini(X,unbiased=FALSE,conf.level=0.9)[2:3]
print(IC_bootstrap)

