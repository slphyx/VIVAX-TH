# slphyx@SHIFT-ENTER
library(deSolve)
library(BayesianTools) # for the MCMC modules

## data
TH.data <- data.frame(
  year = c(2012,2013,2014,2015,2016,2017,2018,2019,2020,2021),
  cases = c(18985,18456,19318,14251,14466,9612,5548,4589,3635,3084)
)

model_equations <- function(t, x, parameters) {

  with(as.list(parameters), {
    IL <- x[1]
    I0 <- x[2]
    SL <- x[3]
    S0 <- x[4]
    TL <- x[5]
    T0 <- x[6]
    h <- x[7]
    hL <- x[8]
    hh <- x[9]
    hhL <- x[10]

    dIL= (1-alpha)*(omega*lambda*(IL+I0+TL+T0)+delta)*(S0+SL) + (1-alpha)*f*SL + (omega*lambda*(IL+I0+TL+T0)+delta)*I0 - gammaL*IL - r*IL;
    dI0= -(omega*lambda*(IL+I0+TL+T0)+delta)*(I0) +  gammaL*IL - r*I0;
    dSL= -(omega*lambda*(IL+I0+TL+T0)+delta)*(SL) - f*SL +(1-beta)*sigma*TL -  gammaL*SL + r*IL+ r*TL;
    dS0= -(omega*lambda*(IL+I0+TL+T0)+delta)*(S0) +beta*sigma*TL+ sigma *T0 +  gammaL*SL + r*I0+ r*T0;
    dTL= alpha*(omega*lambda*(IL+I0+TL+T0)+delta)*(S0+SL) + alpha*f*SL - sigma*TL -r*TL -gammaL*TL+ (omega*lambda*(IL+I0+TL+T0)+delta)*T0;
    dT0= -(omega*lambda*(IL+I0+TL+T0)+delta)*(T0) - sigma*T0+gammaL*TL - r*T0;

    dh=rho*((omega*lambda*(IL+I0+TL+T0)+delta)*(S0+SL) +f*SL);
    dhL= rho*((omega*lambda*(IL+I0+TL+T0))*(S0+SL) +f*SL);
    dhh= ((omega*lambda*(IL+I0+TL+T0)+delta)*(S0+SL) +f*SL);
    dhhL= ((omega*lambda*(IL+I0+TL+T0))*(S0+SL) +f*SL);

    list(c(dIL, dI0, dSL, dS0, dTL, dT0, dh, dhL, dhh, dhhL))
  })
}

## function for optimization
# sum square error
SSQ <- function(obs, model) {
  sum((obs-model)^2)
}


control_options <- list(
  atol = 1e-12,  # Absolute tolerance, similar to TolX
  rtol = 1e-12,  # Relative tolerance, similar to TolFun
  maxsteps = 10000,  # Maximum number of steps, similar to MaxIter
  hmax = 0.1,  # Maximum step size (optional, similar to MaxStep in MATLAB)
  hini = 1e-6  # Initial step size (optional)
)


OBJ.fn <- function(pars){
  print(pars)

  times <- seq(1, 10)

  lambda <- pars[1]
  r <- 1/60
  gammaL <- 1/223
  f <- 1/72
  alpha <- 0.17
  beta <- 0.43
  sigma <- 1/15
  omega <- 0.26
  rho <- 0.5
  delta <- pars[2]
  parameters <- c(lambda = lambda,r = r,gammaL = gammaL,f = f,
                  alpha = alpha,beta=beta,sigma = sigma,omega = omega,rho =rho,delta=delta)


  IL0 <- 500
  I00 <- 500

  TL0 <- 500
  T00 <- 200
  SL0 <- ((1-beta)*sigma*TL0+r*IL0+r*TL0)/(omega*lambda*(IL0+I00+TL0+T00))
  S00 <- (beta*sigma*TL0+sigma*T00+gammaL*SL0+gammaL*I00+r*T00)/(omega*lambda*(IL0+I00+TL0+T00))
  h0 <- 19090
  hL0 <- 1000
  hh0 <- 3084
  hhL0 <-1000



  initial_state <- c(IL = IL0, I0 = I00, SL = SL0, S0 = S00, TL = TL0, T0 = T00,
                     h = h0, hL = hL0, hh = hh0, hhL = hhL0)

  result <- ode(y = initial_state, times = times, func = model_equations, parms = parameters,
                method = rkMethod("ode45"), atol = 1e-10, rtol = 1e-10, maxsteps = 10e5)

  inc <- c(result[,"h"][1],diff(result[,"h"]))

  SSQ(obs = TH.data$cases, model = inc)

}

########################
Plot.par <- function(pars){
  print(pars)

  times <- seq(1, 10)

  lambda <- pars[1]
  r <- 1/60
  gammaL <- 1/223
  f <- 1/72
  alpha <- 0.17
  beta <- 0.43
  sigma <- 1/15
  omega <- 0.26
  rho <- 0.5
  delta <- pars[2]
  parameters <- c(lambda = lambda,r = r,gammaL = gammaL,f = f,
                  alpha = alpha,beta=beta,sigma = sigma,omega = omega,rho =rho,delta=delta)


  IL0 <- 500
  I00 <- 500

  TL0 <- 500
  T00 <- 200
  SL0 <- ((1-beta)*sigma*TL0+r*IL0+r*TL0)/(omega*lambda*(IL0+I00+TL0+T00))
  S00 <- (beta*sigma*TL0+sigma*T00+gammaL*SL0+gammaL*I00+r*T00)/(omega*lambda*(IL0+I00+TL0+T00))
  h0 <- 19090
  hL0 <- 1000
  hh0 <- 3084
  hhL0 <-1000



  initial_state <- c(IL = IL0, I0 = I00, SL = SL0, S0 = S00, TL = TL0, T0 = T00,
                     h = h0, hL = hL0, hh = hh0, hhL = hhL0)

  result <- ode(y = initial_state, times = times, func = model_equations, parms = parameters,
                method = rkMethod("ode45"), atol = 1e-10, rtol = 1e-10, maxsteps = 10e5)

  inc <- c(result[,"h"][1],diff(result[,"h"]))

  ssq <- SSQ(obs = TH.data$cases, model = inc)
  plot(TH.data$year, TH.data$cases, ylim = c(0,3e4))
  lines(TH.data$year,inc)
  title(paste("SSQ:", ssq))

}

times <- seq(1, 10)

lambda <- 2.95e-6
r <- 1/60
gammaL <- 1/223
f <- 1/72
alpha <- 0.17
beta <- 0.43
sigma <- 1/15
omega <- 0.26
rho <- 0.5
delta <- 2.6e-1
parameters <- c(lambda = lambda,r = r,gammaL = gammaL,f = f,
                alpha = alpha,beta=beta,sigma = sigma,omega = omega,rho =rho,delta=delta)


IL0 <- 500
I00 <- 500

TL0 <- 500
T00 <- 200
SL0 <- ((1-beta)*sigma*TL0+r*IL0+r*TL0)/(omega*lambda*(IL0+I00+TL0+T00))
S00 <- (beta*sigma*TL0+sigma*T00+gammaL*SL0+gammaL*I00+r*T00)/(omega*lambda*(IL0+I00+TL0+T00))
h0 <- 19090
hL0 <- 1000
hh0 <- 3084
hhL0 <-1000

initial_state <- c(IL = IL0, I0 = I00, SL = SL0, S0 = S00, TL = TL0, T0 = T00,
                   h = h0, hL = hL0, hh = hh0, hhL = hhL0)


prior <- createUniformPrior(lower = c(1e-7,0),
                            upper = c(9e-6,0.5))


pois_log_likelihood <- function(pars, obs = TH.data,
                                           plot= F) {
  times <- seq(1, 10)
  lambda <- pars[1]
  r <- 1/60
  gammaL <- 1/223
  f <- 1/72
  alpha <- 0.17
  beta <- 0.43
  sigma <- 1/15
  omega <- 0.26
  rho <- 0.5
  delta <- pars[2]
  parameters <- c(lambda = lambda,r = r,gammaL = gammaL,f = f,
                  alpha = alpha,beta=beta,sigma = sigma,omega = omega,rho =rho,delta=delta)


  IL0 <- 500
  I00 <- 500

  TL0 <- 500
  T00 <- 200
  SL0 <- ((1-beta)*sigma*TL0+r*IL0+r*TL0)/(omega*lambda*(IL0+I00+TL0+T00))
  S00 <- (beta*sigma*TL0+sigma*T00+gammaL*SL0+gammaL*I00+r*T00)/(omega*lambda*(IL0+I00+TL0+T00))
  h0 <- 19090
  hL0 <- 1000
  hh0 <- 3084
  hhL0 <-1000



  initial_state <- c(IL = IL0, I0 = I00, SL = SL0, S0 = S00, TL = TL0, T0 = T00,
                     h = h0, hL = hL0, hh = hh0, hhL = hhL0)


  # run the model
  result <- ode(y = initial_state, times = times, func = model_equations, parms = parameters,
                method = rkMethod("ode45"), atol = 1e-10, rtol = 1e-10, maxsteps = 10e5)

  model.inc <- c(result[,"h"][1],diff(result[,"h"]))

  ll <- sum(dpois(x=TH.data$cases,
                  lambda=  model.inc,
                  log = T))

  if(plot == T){
    plot(TH.data$year, TH.data$cases, ylim = c(0,3e4))
    lines(TH.data$year,model.inc)
  }

  if(ll == Inf || is.na(ll) || is.nan(ll) ){
    ll = -Inf
  }

  return(ll)
}

#check LL
pois_log_likelihood(c( 2.95e-6,2.6e-1),plot=T)
pois_log_likelihood(c( 9e-6,2.6e-1),plot=T)


# run MCMC

bayesianSetup <- createBayesianSetup(pois_log_likelihood,
                                     prior = prior,
                                     names =
                                       c( "lambda",
                                          "delta"))


settings <- list(iterations =33000, nrChains = 3, message = TRUE,burnin = 3000)
mcmcout <- runMCMC(bayesianSetup = bayesianSetup, settings = settings)

summary(mcmcout)
tracePlot(mcmcout,start=2000)
plot(mcmcout)

DIC(mcmcout)

sample <- getSample(mcmcout, coda = TRUE, parametersOnly = T, 
                    start = 5000)
summary(sample)

saveRDS(mcmcout,"SS_mcmcout.rds")

