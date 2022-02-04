#### Instantaneous Numerical Solver ####

#' Profit Function
#'
#' The profit function passed to the \eqn{\Delta\psi} optimizer. It calculates 
#' the profit, defined as \eqn{(A - \gamma\Delta\psi^2)}
#'
#' @param par numeric, A vector of variables to be optimized, namely, \eqn{c(log(Jmax), \Delta\psi)}
#' @param psi_soil numeric, soil water potential (Mpa)
#' @param jmax numeric, Jmax (potentially acclimated to previous conditions) (umol/m2/s)
#' @param vcmax numeric, Vcmax (potentially coordinated to previously acclimated Jmax) (umol/m2/s)
#' @param par_cost numeric, List of named cost parameters: list(alpha = xx, gamma = xx)
#' @param par_photosynth numeric, Photosynthesis parameters
#' @param par_env numeric, Environmental parameters
#' @param opt_hypothesis character, Either "Lc" or "PM"
#'
#' @return Net assimilation rate after accounting for costs (profit) (umol/m2/s)
#'
#' @export
#'
# Calculate the net assimilation, i.e. assimilation - costs (umol/m2/s), given 
# - par = c(jmax, dpsi), the optimization parameters
# - Jmax
# - Vcmax
# - soil water potential (Mpa)
# - cost parameters (a, a1, b, c, d)
# - photosynthesis parameters 
# - Plant parameters
# - Env parameters
fn_profit_instantaneous = function(par, jmax, vcmax, psi_soil, par_cost, par_photosynth, par_plant, par_env, do_optim=F){
  dpsi = par
  gs = calc_gs(dpsi, psi_soil, par_plant, par_env)
  A = calc_assimilation_limiting(vcmax, jmax, gs, par_photosynth)$a
  profit = A - par_cost$gamma * dpsi^2  
  
  if (do_optim){
    return(-profit)
  } else {
    return(profit)
  }
}


optimise_shortterm <- function(fn_profit_inst, jmax, vcmax, psi_soil, par_cost, par_photosynth, par_plant, par_env, return_all = FALSE){
  
  out_optim <- optimr::optimr(
    par       = c(dpsi=1),  
    lower     = c(0.001),
    upper     = c(20),
    fn        = fn_profit_inst,
    psi_soil  = psi_soil,
    jmax      = jmax,
    vcmax     = vcmax,
    par_cost  = par_cost,
    par_photosynth = par_photosynth,
    par_plant = par_plant,
    par_env   = par_env,
    do_optim  = TRUE,
    # opt_hypothesis = opt_hypothesis,
    method    = "L-BFGS-B",
    control   = list() 
  )
  
  out_optim$value <- -out_optim$value
  
  if (return_all){
    out_optim
  } else {
    return(out_optim$par)
  }
}
