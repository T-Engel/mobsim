#' Calculate ENS
#' 
#' ENS is calculated from the unbiased PIE.
#' @param x community matrix or vector
#' @param ENS boolian. Should ENS be returned? If \code{ENS=FALSE}, PIE will be returned.
#' 
calc_ENS = function(x, ENS=TRUE) {

  x = drop(as.matrix(x))
  if (any(x < 0, na.rm = TRUE)) 
    stop("input data must be non-negative")
  if (length(dim(x)) > 1) {
    total = apply(x, 1, sum)
    x = sweep(x, 1, total, "/")
  } else {
    total = sum(x)
    x = x / total
  }
  x = x * x
  if (length(dim(x)) > 1) {
    H = rowSums(x, na.rm = TRUE)
  } else {
    H = sum(x, na.rm = TRUE)
  }
  # calculate PIE without replacement (for total >= 2)
  H = ifelse(total < 2, NA, (total / (total - 1) * (1 - H)))
  if (ENS) {
    # convert to effective number of species (except for PIE == 1)
    H = ifelse(H==1, NA, (1/ (1-H)))
  }     
  return(H)
}
#' Simulate lognorm community with prescribed ENS
#' 
#' Wrapper for sim_sad optimizing parameters of lnorm
#' distribution to create community with desired ENS.
#' The unbiased version of ENS is used here.
#' 
#' @param ENS_target positiv numeric. Desired ENS value of the simulated community.
#' @param S integer. Number of species to be simulated.
#' @param mean_abu mean species abundance. Will be used to calculate total number of individuals \code{n_sim}.
#' @param n_sim total number of individuals. If \code{n_sim=NULL}, it will be calculated as \code{mean_abu*S}
#' @param sims number of repetitions used during the optimisation process (using different random seeds, respectively).
sim_ENS<-function(ENS_target,S,mean_abu=10, n_sim=NULL, sims=100){
  if(is.null(n_sim)) n_sim<-S*mean_abu
  seeds=sample(1:(sims*1000), sims,replace = F)
  
  
  testsad<-function(par, S=S, seed=sample(1:1000)){
    set.seed(seed)
    com<-sim_sad(s_pool=S,n_sim=n_sim,sad_type ="lnorm",sad_coef = list(cv_abund=par[1]),fix_s_sim=T)
    return(abs(ENS_target-calc_ENS(com)))
  }
  pars<-sapply(seeds, function(seed){
    pars_opt<-optim(par=c(0),testsad,S=S,seed=seed,method = "Brent",lower = 0, upper = S,control = list(maxit=10000))
    return(c(seed=seed,par=pars_opt$par,diff=pars_opt$value))
  }
  )
  pars<-as.data.frame(t(pars))
  pars<-subset(pars,diff==min(pars$diff))
  set.seed(pars$seed)
  com<-sim_sad(s_pool=S,n_sim=n_sim,sad_type ="lnorm",sad_coef = list(cv_abund=pars$par),fix_s_sim=T)
  return(com)
}
