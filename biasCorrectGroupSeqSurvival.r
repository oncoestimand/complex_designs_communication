# Note: These functions are exploratory only and have not been formally validated!

simulateGroupSeqSurvivalTrialResults <- function(
  n,hazardRatio,
  informationRates, maxNumberOfEvents, allocationRatioPlanned = 1,
  boundaryEffectScaleLower= -Inf,boundaryEffectScaleUpper= Inf){
  # Simulate results of hypothetical group-sequential survival trial
  # Simulation is based on the canonical joint multivariate normal distribtion of treatment effect estimates across interims
  # 
  # Arguments:
  # - n:                      Number of simulation runs
  # - hazardRatio:            True hazard ratio
  # - informationRates:       Information rates at interim analyses (as in rpact) - vector of length K with last value =1 
  # - maxNumberOfEvents:      Planned number of events at final analysis
  # - allocationRatioPlanned: Planned allocation ratio 
  # - boundaryEffectScaleLower, boundaryEffectScaleUpper: Stopping boundaries on hazard ratio scale (MDDs)
  #---------------------------------------------------------------------------------------------------------------------------------------
  require(MASS)
  p <- length(informationRates)
  # calculate covariance matrix
  information <- informationRates*maxNumberOfEvents*allocationRatioPlanned/(1+allocationRatioPlanned)^2
  variance <- 1/information
  covMat <- outer(variance,variance,FUN=function(x,y) pmin(x,y))
  # transform boundaries to log-scale
  if (length(boundaryEffectScaleLower)==1) { 
    boundaryEffectScaleLower <- rep(boundaryEffectScaleLower,p-1)
  }
  if (length(boundaryEffectScaleUpper)==1) { 
    boundaryEffectScaleUpper <- rep(boundaryEffectScaleUpper,p-1)
  }
  boundaryEffectScaleLower <- log(boundaryEffectScaleLower)
  boundaryEffectScaleUpper <- log(boundaryEffectScaleUpper)
  # simulate group-sequential trial results
  simTrials <- mvrnorm(n,mu=rep(log(hazardRatio),p),Sigma=covMat)
  estimate <- simTrials[,p] # estimator if trial went to final analysis
  stopstage <- rep(p,n)
  for (j in ((p-1):1)){ 
    # update estimate "backwards in time" with interim results if boundary was crossed at interim
    estimate <- ifelse((simTrials[,j]<boundaryEffectScaleLower[j])|(simTrials[,j]>boundaryEffectScaleUpper[j]), simTrials[,j], estimate)
    stopstage <- ifelse((simTrials[,j]<boundaryEffectScaleLower[j])|(simTrials[,j]>boundaryEffectScaleUpper[j]), j, stopstage)
  }
  data.frame(hazardRatioEstimate=exp(estimate),logHazardRatioEstimate=estimate,stopstage=stopstage)
}


biasCorrectGroupSeqSurvival<- function(
  hazardRatioMLE, stopstage, type = c("global","stagewise"),
  n=1E6, 
  informationRates, maxNumberOfEvents, allocationRatioPlanned = 1,
  boundaryEffectScaleLower= -Inf,boundaryEffectScaleUpper= Inf,
  seed=12, 
  lower=NULL, upper=NULL){
  # Calculate global or stagewise bias-adjusted estimate for a survival trial based on simulation
  # 
  # Arguments:
  # - hazardRAtioMLE:         Standard maximum likelihood estimator of the HR at the stage where the trial stopped
  # - stopstage:              Stage at which the trial stopped in case the stagewise bias-adjusted estimat should be calculated
  # - type:                   "global" or "stagewise" adjusted estimate
  # - n:                      Number of simulation runs
  # - hazardRatio:            True hazard ratio
  # - informationRates:       Information rates at interim analyses (as in rpact) - vector of length K with last value =1 
  # - maxNumberOfEvents:      Planned number of events at final analysis
  # - allocationRatioPlanned: Planned allocation ratio 
  # - boundaryEffectScaleLower, boundaryEffectScaleUpper: Stopping boundaries on hazard ratio scale (MDDs)
  # - seed:                   Random seed to guarantee reproducibility
  # - lower, upper:           Lower of upper bound of plausible bias-adjusted estimates (only needed if default causes numerical problems)
  #---------------------------------------------------------------------------------------------------------------------------------------
  # set lower and upper bound for optimisation as estimate +/- 3*SE (if not provided)
  seStop <- sqrt((1+allocationRatioPlanned)^2/(allocationRatioPlanned*informationRates[stopstage]*maxNumberOfEvents))
  if (is.null(lower)){ lower <- exp(log(hazardRatioMLE)-3*seStop) } 
  if (is.null(upper)){ upper <- exp(log(hazardRatioMLE)+3*seStop) }
  root <- function(newHazardRatio){
    set.seed(seed)
    sim <- simulateGroupSeqSurvivalTrialResults(n,hazardRatio=newHazardRatio,
                                                informationRates, maxNumberOfEvents, allocationRatioPlanned,
                                                boundaryEffectScaleLower,boundaryEffectScaleUpper)
    if (type=="global") { 
      bias <- mean(sim$logHazardRatioEstimate)-log(newHazardRatio) 
    } else if (type=="stagewise"){
      bias <- mean(sim$logHazardRatioEstimate[sim$stopstage==stopstage])-log(newHazardRatio)  
    } 
    log(newHazardRatio)-log(hazardRatioMLE)+bias # this corresponds to formula (4.14)/(4.15) in Wassmer&Brannath (page 98) 
  }
  uniroot(root,lower=lower,upper=upper)$root
}