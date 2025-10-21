#
# Template for
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
n = xx #sample size on each treatment
sigma = xx # standard deviation
alpha = xx # error rate

mu_0 = 0 # mean in the control group
mu_1_null = xx # mean in the experimental group under the null
mu_1_alt = xx # mean in the experimental group under the alternative

theta_null = mu_1_null-mu_0
theta_alt = mu_1_alt-mu_0


##############################
### Task 1
##############################

#simulating patient observations
x = rnorm( )
y = rnorm( )

#find the mean in each treatment group
xbar = mean(x)
ybar = mean(y)

#estimate the treatment effect (from the lecture)
thetahat = xbar-ybar
z = xxx  # normalise the observed different to obtain Z-statistic
  
#find the critical region for the test
k = xxx 

#is the test rejected
z>k

##############################
### Task 2
##############################


#exercise 1 made into a function
sim = function(n,mu_1,mu_0,sigma){
  #simulating patient observations
  
  xxx
  
  #find the mean in each treatment group
  
  xxx
  
  #estimate the treatment effect
  thetahat = xbar-ybar
  
  #output the treatment effect
  thetahat
}

#number of trials to simulate
m = xxx

#create a vector to store results
thetahat_null = rep(0,m)
thetahat_alt = rep(0,m)

#use loop to simulate multiple trials
for(i in 1:m){
  
  xxx
  
}

# find corresponding z-statistic (normalising the output above)
z_null = xxx
z_alt = xxx

#results of which hypotheses are rejected (for each trial 0/1)
res_null = xxx
res_alt = xx

#probability of rejecting the null hypothesis
mean(res_null)
mean(res_alt)

#plotting the z-statistics to see what we have produced

xxx


##############################
### Task 3
##############################

xxx

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
    xxx
  }
  
  #output results
  thetahat
}

# Reproducing results above but just with simulation function
thetahat_null = xxx
thetahat_alt = xxx

#find corresponding z-statistic
z_null = xxx
z_alt = xxx

#results of which hypotheses are rejected, simulating every patient
res_null = xx
res_alt = xxx

#computation with summary statistics (means)
thetahat_null_direct = xxx
thetahat_alt_direct = xxx

#find corresponding z-statistic
z_null_direct = xxx
z_alt_direct = xxx

#results of which hypotheses are rejected, simulating summary statistics
res_null_direct = xxx
res_alt_direct = xxx


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
#simulating patients's responses
system.time(xxx)
#simulating summary statistics
system.time(xxx)
