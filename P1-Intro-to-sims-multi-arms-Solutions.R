#
# Solutions to 
# Practical 1: Introduction to simulation in R
# and simulations of multi-arm trials
#
# by Dr Pavel Mozgunov and Dr Dominique-Laurent Couturier
# MRC Biostatistics Unit, 22nd of October 2025
#
#
# 
##############################
# Setting up the parameters
##############################
n = 170 #sample size on each treatment
sigma = 1 # standard deviation
alpha = 0.05 # error rate

mu_0 = 0 # mean in the control group
mu_1_null = 0 # mean in the experimental group under the null
mu_1_alt = 0.3 # mean in the experimental group under the alternative

theta_null = mu_1_null-mu_0
theta_alt = mu_1_alt-mu_0


##############################
### Task 1
##############################

#simulating patient observations
x = rnorm(n,mu_1_null,sigma)
y = rnorm(n,mu_0,sigma)

#find the mean in each treatment group
xbar = mean(x)
ybar = mean(y)

#estimate the treatment effect (from the lecture)
thetahat = xbar-ybar
z = thetahat/(sigma*sqrt(2/n))

#find the critical region for the test
k = qnorm(1-alpha)

#is the test rejected
z>k

##############################
### Task 2
##############################


#exercise 1 made into a function
sim = function(n,mu_1,mu_0,sigma){
  #simulating patient observations
  x = rnorm(n,mu_1,sigma)
  y = rnorm(n,mu_0,sigma)
  
  #find the mean in each treatment group
  xbar = mean(x)
  ybar = mean(y)
  
  #estimate the treatment effect
  thetahat = xbar-ybar
  
  #output the treatment effect
  thetahat
}

#number of trials to simulate
m = 10000

#create a vector to store results
thetahat_null = rep(0,m)
thetahat_alt = rep(0,m)

#use loop to simulate multiple trials
for(i in 1:m){
  thetahat_null[i] = sim(n,mu_1=mu_1_null,mu_0,sigma)
  thetahat_alt[i] = sim(n,mu_1=mu_1_alt,mu_0,sigma)
  
}

# find corresponding z-statistic
z_null = thetahat_null/(sigma*sqrt(2/n))
z_alt = thetahat_alt/(sigma*sqrt(2/n))

#results of which hypotheses are rejected 
res_null = z_null>k
res_alt = z_alt>k

#probability of rejecting the null hypothesis
mean(res_null)
mean(res_alt)

#plotting the z-statistics to see what we have produced
hist(z_null,xlim=c(-4,6),col=rgb(1,0,0,1/4),main="Histogramm of Z-statistic")
hist(z_alt,add=T,col=rgb(0,0,1,1/4))

#adding the rejection region will be informative
abline(v=k,col="red")

##############################
### Task 3
##############################

#the standard deviation of our estimate
sigmahat_null= sqrt(mean(res_null)*(1-mean(res_null))/m)
sigmahat_null
sigmahat_alt= sqrt(mean(res_alt)*(1-mean(res_alt))/m)
sigmahat_alt

#pick a target standard deviation
#this is a reasonable pick for stability in the second decimal place
tar = 0.001

#what m should we use?
mhat_null = round((mean(res_null)*(1-mean(res_null)))/(tar^2))
mhat_null
mhat_alt = round((mean(res_alt)*(1-mean(res_alt)))/(tar^2))
mhat_alt


##############################
### Task 4
##############################

# Creat the function for the simulations [optional but convinient]
simloop = function(m,n,mu_1,mu_0,sigma)
{
  #create a vector to store results
  thetahat = rep(0,m)
  
  #use loop to simulate multiple trials
  for(i in 1:m){
    thetahat[i] = sim(n,mu_1,mu_0,sigma)
  }
  
  #output results
  thetahat
}

# Reproducing results above but just with simulation function
thetahat_null = simloop(m,n,mu_1_null,mu_0,sigma)
thetahat_alt = simloop(m,n,mu_1_alt,mu_0,sigma)

#find corresponding z-statistic
z_null = thetahat_null/(sigma*sqrt(2/n))
z_alt = thetahat_alt/(sigma*sqrt(2/n))

#results of which hypotheses are rejected, simulating every patient
res_null = z_null>k
res_alt = z_alt>k

#computation with summary statistics
thetahat_null_direct = rnorm(m,theta_null,sigma*sqrt(2/n))
thetahat_alt_direct = rnorm(m,theta_alt,sigma*sqrt(2/n))

#find corresponding z-statistic
z_null_direct = thetahat_null_direct/(sigma*sqrt(2/n))
z_alt_direct = thetahat_alt_direct/(sigma*sqrt(2/n))

#results of which hypotheses are rejected, simulating summary statistics
res_null_direct = z_null_direct>k
res_alt_direct = z_alt_direct>k


#probability of rejecting the null hypothesis with each method
#all patients
mean(res_null)
#summary statistics
mean(res_null_direct)

#all patients
mean(res_alt)
#summary statistics
mean(res_alt_direct)

#note the time difference between each simulation
#simulating all patients
system.time(simloop(m,n,mu_1_alt,mu_0,sigma))
#simulating summary statistics
system.time(rnorm(m,theta_alt,sigma*sqrt(2/n)))

# What if we want 100,000 simulations?
system.time(simloop(10^5,n,mu_1_alt,mu_0,sigma))
#simulating summary statistics
system.time(rnorm(10^5,theta_alt,sigma*sqrt(2/n)))


##############################
### Task 5
##############################

# Add mean for second experimental treatment
mu_2_null = 0
mu_2_alt = 0.3

theta_2_null = mu_2_null - mu_0
theta_2_alt = mu_2_alt - mu_0

theta_1_null = mu_1_null - mu_0
theta_1_alt = mu_1_alt - mu_0

#Option 1, simulate xbar, ybar and wbar
ybar = rnorm(m,mu_0,sigma/sqrt(n))
xbar = rnorm(m,mu_1_null,sigma/sqrt(n))
wbar = rnorm(m,mu_2_null,sigma/sqrt(n))

#computing corresponding estimates
thetahat_1_null = xbar - ybar
thetahat_2_null = wbar - ybar

#find corresponding z-statistics
z_1_null = thetahat_1_null/(sigma*sqrt(2/n))
z_2_null = thetahat_2_null/(sigma*sqrt(2/n))

#and for the Dunnett test
z_d_null = (z_1_null>=z_2_null)*z_1_null + (z_1_null<z_2_null)*z_2_null

# Displaying the distributions
plot(z_1_null,z_2_null) # positively correlated
hist(z_d_null)
abline(v=k,col="red")

# Rejection probabilities
# Treatment 1
mean(z_1_null > k)
#  Treatment 2
mean(z_2_null > k)
# Dunnett statistic
mean(z_d_null > k)

# Funding the correct critical value
# Calibration
mean(z_d_null > 1.8)
mean(z_d_null > 1.9)
mean(z_d_null > 2.0)
mean(z_d_null > 1.95)
mean(z_d_null > 1.925)
mean(z_d_null > 1.91)

# Exact solution
corr<-matrix(c(1,0.5,0.5,1),ncol=2)

require(mvtnorm)
D = qmvnorm((1-alpha),tail="lower.tail",mean=c(0,0),corr=corr)$quantile
D
mean(z_d_null > D)

##############################
### Task 6
##############################

# Variance-covariance matrix
Sigma<-(2*sigma^2/n)*corr
thetahat = rmvnorm(m,mean=c(theta_1_null,theta_2_null),Sigma)

#find the corresponding z-statistics
z = matrix(0,nrow=m,ncol=3)
z[,1] = thetahat[,1]/(sigma*sqrt(2/n))
z[,2] = thetahat[,2]/(sigma*sqrt(2/n))

#find the corresponding Dunnett test statistic
z[,3] = (z[,1]>=z[,2])*z[,1] + (z[,1]<z[,2])*z[,2]

#rejection probabilities
#local rejection
#H_01
mean(z[,1] > D)
#H_02
mean(z[,2] > D)
#H_012
mean(z[,3] > D)


##############################
### Task 7
##############################

# Under LFC
mu_2_un = 0.1
theta_2_un = mu_2_un - mu_0

thetahat = rmvnorm(m,mean=c(theta_1_alt,theta_2_un),Sigma)

# Find the corresponding z-statistics
z = matrix(0,nrow=m,ncol=3)
z[,1] = thetahat[,1]/(sigma*sqrt(2/n))
z[,2] = thetahat[,2]/(sigma*sqrt(2/n))
z[,3] = (z[,1]>=z[,2])*z[,1] + (z[,1]<z[,2])*z[,2]

#rejection probabilities
mean(z[,3] > D) # H_012 # Is it power? NO!
mean((z[,1]>z[,2]) & (z[,1]>D)) # This is power
mean((z[,2]>z[,1]) & (z[,2]>D)) # This is rejection of uninteresting arm

# Under Global Alternative

thetahat = rmvnorm(m,mean=c(theta_1_alt,theta_2_alt),Sigma)

# Find the corresponding z-statistics
z = matrix(0,nrow=m,ncol=3)
z[,1] = thetahat[,1]/(sigma*sqrt(2/n))
z[,2] = thetahat[,2]/(sigma*sqrt(2/n))
z[,3] = (z[,1]>=z[,2])*z[,1] + (z[,1]<z[,2])*z[,2]

#rejection probabilities
mean(z[,3] > D) # H_012 # Is it power? YES!


# Sample size under LFC
n2<-240
Sigma<-(2*sigma^2/n2)*corr

thetahat = rmvnorm(m,mean=c(theta_1_alt,theta_2_un),Sigma)

# Find the corresponding z-statistics under new sample size n2
z = matrix(0,nrow=m,ncol=3)
z[,1] = thetahat[,1]/(sigma*sqrt(2/n2))
z[,2] = thetahat[,2]/(sigma*sqrt(2/n2))
z[,3] = (z[,1]>=z[,2])*z[,1] + (z[,1]<z[,2])*z[,2]

mean((z[,1]>z[,2]) & (z[,1]>D)) # This is power


##############################
### Task 8
##############################
#
# Adjustment for the baseline covariate
#
######
# Point (a)
corr.base<-0.0
patient.corr<-matrix(c(1,corr.base,corr.base,1),ncol=2)
mu_base<-0

z1<-c()
z2<-c()
nsims<-10000
for(z in 1:nsims){
  x = rmvnorm(n,mean=c(mu_base,mu_1_null),patient.corr)
  y = rmvnorm(n,mean=c(mu_base,mu_0),patient.corr)
  w = rmvnorm(n,mean=c(mu_base,mu_2_null),patient.corr)
  
  x<-cbind(x,rep(1,n),rep(0,n))
  y<-cbind(y,rep(0,n),rep(0,n))
  w<-cbind(w,rep(0,n),rep(1,n))
  
  data<-rbind(x,y,w)
  
  model<-lm(data[,2] ~ data[,1] + data[,3] + data[,4])
  output<-summary(model)
  z1[z]<-output$coefficients[3,1]/output$coefficients[3,2]
  z2[z]<-output$coefficients[4,1]/output$coefficients[4,2]
}
z.d<- (z1>=z2)*z1 + (z1<z2)*z2
mean(z.d>D) 

######
# Point (b)
corr.base<-0.0
patient.corr<-matrix(c(1,corr.base,corr.base,1),ncol=2)
mu_base<-0

z1<-c()
z2<-c()
nsims<-10000
for(z in 1:nsims){
  x = rmvnorm(n,mean=c(mu_base,mu_1_alt),patient.corr)
  y = rmvnorm(n,mean=c(mu_base,mu_0),patient.corr)
  w = rmvnorm(n,mean=c(mu_base,mu_2_un),patient.corr)
  
  x<-cbind(x,rep(1,n),rep(0,n))
  y<-cbind(y,rep(0,n),rep(0,n))
  w<-cbind(w,rep(0,n),rep(1,n))
  
  data<-rbind(x,y,w)
  
  model<-lm(data[,2] ~ data[,1] + data[,3] + data[,4])
  output<-summary(model)
  z1[z]<-output$coefficients[3,1]/output$coefficients[3,2]
  z2[z]<-output$coefficients[4,1]/output$coefficients[4,2]
}
z.d<- (z1>=z2)*z1 + (z1<z2)*z2
mean(z.d>D) 

# Now with correlation
corr.base<-0.5
patient.corr<-matrix(c(1,corr.base,corr.base,1),ncol=2)
mu_base<-0

z1<-c()
z2<-c()
nsims<-10000
for(z in 1:nsims){
  x = rmvnorm(n,mean=c(mu_base,mu_1_alt),patient.corr)
  y = rmvnorm(n,mean=c(mu_base,mu_0),patient.corr)
  w = rmvnorm(n,mean=c(mu_base,mu_2_un),patient.corr)
  
  x<-cbind(x,rep(1,n),rep(0,n))
  y<-cbind(y,rep(0,n),rep(0,n))
  w<-cbind(w,rep(0,n),rep(1,n))
  
  data<-rbind(x,y,w)
  
  model<-lm(data[,2] ~ data[,1] + data[,3] + data[,4])
  output<-summary(model)
  z1[z]<-output$coefficients[3,1]/output$coefficients[3,2]
  z2[z]<-output$coefficients[4,1]/output$coefficients[4,2]
}
z.d<- (z1>=z2)*z1 + (z1<z2)*z2
mean(z.d>D) 
mean((z1>z2) & (z1>D))

######
# Point (c)
t<-2 # number of repeated observations
sigma2<-sqrt(t+corr.base*(t*t-t))/t

Sigma2<-(2*sigma2^2/n)*corr
thetahat = rmvnorm(m,mean=c(theta_1_alt,theta_2_un),Sigma)

# Find the corresponding z-statistics under new sample size n2
z = matrix(0,nrow=m,ncol=3)
z[,1] = thetahat[,1]/(sigma2*sqrt(2/n))
z[,2] = thetahat[,2]/(sigma2*sqrt(2/n))
z[,3] = (z[,1]>=z[,2])*z[,1] + (z[,1]<z[,2])*z[,2]

mean((z[,1]>z[,2]) & (z[,1]>D)) # This is power



######
# Point (*)
# Confirming the found sample size after the adjustment (and no corr)
corr.base<-0.0
patient.corr<-matrix(c(1,corr.base,corr.base,1),ncol=2)
mu_base<-0

z1<-c()
z2<-c()
nsims<-10000
for(z in 1:nsims){
  x = rmvnorm(n2,mean=c(mu_base,mu_1_alt),patient.corr)
  y = rmvnorm(n2,mean=c(mu_base,mu_0),patient.corr)
  w = rmvnorm(n2,mean=c(mu_base,mu_2_un),patient.corr)
  
  x<-cbind(x,rep(1,n2),rep(0,n2))
  y<-cbind(y,rep(0,n2),rep(0,n2))
  w<-cbind(w,rep(0,n2),rep(1,n2))
  
  data<-rbind(x,y,w)
  
  model<-lm(data[,2] ~ data[,1] + data[,3] + data[,4])
  output<-summary(model)
  z1[z]<-output$coefficients[3,1]/output$coefficients[3,2]
  z2[z]<-output$coefficients[4,1]/output$coefficients[4,2]
}
z.d<- (z1>=z2)*z1 + (z1<z2)*z2
mean(z.d>D) 
mean((z1>z2) & (z1>D))
