# slphyx@SHIFT-ENTER
library(deSolve)

## data
years <- 2012:2035
TH.data <- data.frame(
  year = years,
  cases = c(c(18985,18456,19318,14251,14466,9612,5548,4589,3635,3084),rep(NA,14))
)

model_equations <- function(t, x, parameters) {
  # Use 'with' to handle parameter values
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
    
    
    
    
    # Return the derivatives as a list
    list(c(dIL, dI0, dSL, dS0, dTL, dT0, dh, dhL, dhh, dhhL))
  })
}

update_parameter <- function(parameters, param_name, new_value) {
  # Check if the specified parameter exists in the vector
  if (param_name %in% names(parameters)) {
    parameters[param_name] <- new_value
  } else {
    warning(paste("Parameter", param_name, "not found in the parameters vector."))
  }
  return(parameters)
}


model_equations_SS <- function(time, x, parameters, par.name=NULL, par.val=NULL, at.time=NULL) {
  
  if(time >= at.time){
    parms <- update_parameter(parameters = parameters, param_name = par.name, new_value = par.val)
  }else{
    parms <- parameters
  }
  
  with(as.list(parms), {
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
    
    # dIL= (1-alpha)*(omega*lambda*(IL+I0+TL+T0)+delta)*(S0+SL) + (1-alpha)*f*SL + (omega*lambda*(IL+I0+TL+T0)+delta)*I0 - gammaL*IL - r*IL;
    # dI0= -(omega*lambda*(IL+I0+TL+T0)+delta)*(I0) +  gammaL*IL - r*I0;
    # dSL= -(omega*lambda*(IL+I0+TL+T0)+delta)*(SL) - f*SL +(1-beta)*sigma*TL -  gammaL*SL + r*IL+ r*TL;
    # dS0= -(omega*lambda*(IL+I0+TL+T0)+delta)*(S0) +beta*sigma*TL+ sigma *T0 +  gammaL*SL + r*I0+ r*T0;
    # dTL= alpha*(omega*lambda*(IL+I0+TL+T0)+delta)*(S0+SL) + alpha*f*SL - sigma*TL -r*TL -gammaL*TL+ (omega*lambda*(IL+I0+TL+T0)+delta)*T0;
    # dT0= -(omega*lambda*(IL+I0+TL+T0)+delta)*(T0) - sigma*T0+gammaL*TL - r*T0;
    # 
    # dh=rho*((omega*lambda*(IL+I0+TL+T0)+delta)*(S0+SL) +f*SL); 
    # dhL= rho*((omega*lambda*(IL+I0+TL+T0))*(S0+SL) +f*SL);
    # dhh= ((omega*lambda*(IL+I0+TL+T0)+delta)*(S0+SL) +f*SL);
    # dhhL= ((omega*lambda*(IL+I0+TL+T0))*(S0+SL) +f*SL);
    # 
    
    dIL= (1-alpha)*(omega*lambda*(IL+I0+TL+T0)+delta)*(S0+SL) + (1-alpha)*f*SL + (omega*lambda*(IL+I0+TL+T0)+delta)*I0 - gammaL*IL - r*IL-tau*IL-tau*IL;
    dI0= -(omega*lambda*(IL+I0+TL+T0)+delta)*(I0) +  gammaL*IL - r*I0-tau*I0;
    dTL= alpha*(omega*lambda*(IL+I0+TL+T0)+delta)*(S0+SL) + alpha*f*SL - sigma*TL -r*TL -gammaL*TL+ (omega*lambda*(IL+I0+TL+T0)+delta)*T0;
    dT0= -(omega*lambda*(IL+I0+TL+T0)+delta)*(T0) - sigma*T0+gammaL*TL - r*T0;
    dSL= -(omega*lambda*(IL+I0+TL+T0)+delta)*(SL) - f*SL +(1-beta)*sigma*TL -  gammaL*SL + r*IL+ r*TL+tau*IL;
    dS0= -(omega*lambda*(IL+I0+TL+T0)+delta)*(S0) +beta*sigma*TL+ sigma *T0 +  gammaL*SL + r*I0+ r*T0+tau*I0+tau*IL;
    dh=rho*((omega*lambda*(IL+I0+TL+T0)+delta)*(S0+SL) +f*SL); #%rho*((IL+I0+TL+T0)+delta);
    dhL=rho*((omega*lambda*(IL+I0+TL+T0))*(S0+SL) +f*SL); #% rho*((IL+I0+TL+T0));
    dhh=((omega*lambda*(IL+I0+TL+T0)+delta)*(S0+SL) +f*SL); #%((IL+I0+TL+T0)+delta);
    dhhL=((omega*lambda*(IL+I0+TL+T0))*(S0+SL) +f*SL); #%((IL+I0+TL+T0));
    
    
    
    
    # Return the derivatives as a list
    list(c(dIL, dI0, dSL, dS0, dTL, dT0, dh, dhL, dhh, dhhL))
  })
}


#############################
## for changing a parameter value
## par.name - name of the parameter , i.e., "delta"
## par.val - a vector of the values , i.e., c(0.1,0.4,0.8)
## at.time - time for changing the value. in this case will be between 1-24 for year 2012:2035
############################
Run.Scenarios <- function(par.name, par.val, at.time, xlim=c(2012,2035), ylim=c(0,3e4)){
  
  par.len <- length(par.val)
  
  times <- seq(1, 24)  # from 2012:2035
  
  lambda <- 2.651e-6 # from fitting
  r <- 1/60
  gammaL <- 1/223 
  f <- 1/72 
  alpha <- 0.17
  beta <- 0.43
  sigma <- 1/15 
  omega <- 0.26 
  rho <- 0.5
  delta <- 2.845e-1 # from fitting
  tau <- 0.001
  
  parameters <- c(lambda = lambda, r = r,gammaL = gammaL,f = f,
                  alpha = alpha,beta=beta,sigma = sigma,omega = omega,rho =rho,delta=delta, tau=tau)
  
  baseline.pars <- parameters
  
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
  
  # baseline
  result <- ode(y = initial_state, times = times, func = model_equations_SS, parms = baseline.pars, 
                method = rkMethod("ode45"), at.time = 1000)
  
  inc <- c(result[,"h"][1],diff(result[,"h"]))
  
 
  plot(years, c(TH.data$cases,rep(NULL,14)),  
       ylab = "Local reported incidence", xlab = "Year",
       xlim = xlim, ylim = ylim, main = par.name
       )
  
  
  # scenarios
  cl <- rainbow(par.len)
  
  for(i in 1:par.len){
    result.tmp <- ode(y = initial_state, times = times, func = model_equations_SS, 
                    parms = parameters, 
                    method = rkMethod("ode45"), 
                    par.name = par.name,
                    par.val = par.val[i],
                    at.time = at.time
                  )
    
    inc.tmp <- c(result.tmp[,"h"][1],diff(result.tmp[,"h"]))
    lines(years, inc.tmp, col = cl[i])
  }
  legend("topright", legend = paste(par.name, "=", sprintf("%.2f", par.val)), col = cl, lty = 1)
}


# Run.Scenarios(par.name = "alpha", par.val = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7), at.time = 11)
# Run.Scenarios(par.name = "alpha", par.val = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7), at.time = 11,
#               xlim = c(2020,2035), ylim = c(0,5000))
# 
# Run.Scenarios(par.name = "delta", par.val = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7), at.time = 11)
# Run.Scenarios(par.name = "delta", par.val = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7), at.time = 11,
#               xlim = c(2020,2035), ylim = c(0,5000))
# 
# 
# Run.Scenarios(par.name = "sigma", par.val = c(1/11, 1/7, 1/3, 1), at.time = 11,
#               xlim = c(2020,2035), ylim = c(0,5000))


### Manipulate
library(manipulate)
manipulate(
  Run.Scenarios(
    par.name = par.name,
    par.val = c(par.val1, par.val2, par.val3, par.val4),
    at.time = 11,
    xlim = c(2012,2035),
    ylim = c(0, 30000)
  ),
  par.name = picker("lambda", "r", "gammaL", "f", "alpha", "beta", 
                    "sigma", "omega", "rho", "delta", "tau", initial = "rho"),
  par.val1 = slider(min = 0, max = 1, initial = 0.2, step = 0.0001, ticks = F),
  par.val2 = slider(min = 0, max = 1, initial = 0.3, step = 0.0001, ticks = F),
  par.val3 = slider(min = 0, max = 1, initial = 0.4, step = 0.0001, ticks = F),
  par.val4 = slider(min = 0, max = 1, initial = 0.5, step = 0.0001, ticks = F)
)
