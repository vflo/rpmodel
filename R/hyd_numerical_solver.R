#### Numerical Solver ####


#' Coordinated Vcmax
#'
#' Calculates the carboxylation capacity, as coordinated to a given electron-transport limited assimilation rate.
#'
#' @param aj numeric, electron-transport limited assimilation rate (umol/m2/s) 
#' @param ci numeric, leaf-internal CO2 concentration, converted to partial pressure (Pa)
#' @param par_photosynth numeric, A list of photosynthesis parameters, which must include ...... . All concentrations must be converted to partial pressures (Pa)
#'
#' @return Coordinated Vcmax (umol/m2/s)
#'
#' @export
#'
# Calculate coordinated Vcmax (umol/m2/s), given
# - Aj (umol/m2/s)
# - ci, converted to partial pressure (Pa)
# - photosynthesis parameters (K and gamma_star), converted to partial pressures (Pa)
# - b_rd, which is part of photosynthesis parameters
calc_vcmax_coordinated_numerical = function(aj, ci, par_photosynth){
  d = par_photosynth$delta
  vcmax_coord = aj*(ci + par_photosynth$kmm)/(ci*(1-d)- (par_photosynth$gammastar+par_photosynth$kmm*d))
  return(vcmax_coord)
}



#' Profit Function
#'
#' The profit function passed to the optimizer. It calculates 
#' the profit, defined as \eqn{(A - \alpha Jmax - \gamma\Delta\psi^2)}
#'
#' @param par numeric, A vector of variables to be optimized, namely, \eqn{c(log(Jmax), \Delta\psi)}
#' @param psi_soil numeric, soil water potential (Mpa)
#' @param par_cost numeric, List of named cost parameters: list(alpha = xx, gamma = xx)
#' @param par_photosynth numeric, Photosynthesis parameters
#' @param par_env numeric, Environmental parameters
#' @param opt_hypothesis character, Either "Lc" or "PM"
#' @param gs_approximation character, Either "Ohm" or "PM". Ohm is used when gs is calculated as ohm approximation and PM is used when gs is calculated using Penman-Monteith
#'
#' @return Net assimilation rate after accounting for costs (profit) (umol/m2/s)
#'
#' @export
#'
# Calculate the net assimilation, i.e. assimilation - costs (umol/m2/s), given 
# - par = c(jmax, dpsi), the optimization parameters
# - soil water potential (Mpa)
# - cost parameters (a, a1, b, c, d)
# - photosynthesis parameters 
# - Plant parameters
# - Env parameters
# - Optimization hypothesis used (PM or LC)
fn_profit <- function(par, psi_soil, par_cost, par_photosynth, par_plant, par_env, do_optim = FALSE, opt_hypothesis, gs_approximation){
  
  jmax = exp(par[1])  # Jmax in umol/m2/s (logjmax is supplied by the optimizer)
  dpsi = par[2]       # delta Psi in MPa
  
  if (gs_approximation == "Ohm"){
    gs = calc_gs(dpsi, psi_soil, par_plant, par_env)  # gs in mol/m2/s/Mpa
    E = 1.6*gs*(par_env$vpd/par_env$patm)*1e6         # E in umol/m2/s
  } else if (gs_approximation == "PM"){
    PM_params = calc_PM_params(par_env$tc,par_env$patm, par_env$nR, par_env$LAI)
    ga = calc_ga(par_env$u, par_env$ustar, PM_params$R, par_env$tc, par_env$patm)
    gsh2o = calc_gs_PM(dpsi, psi_soil, par_plant_now, par_env_now, PM_params)
    gs = gsh2o/1.6
    E = (S*Q+dens*cp*D*ga)/(L*(S+pch*(1+ga/gsh2o)))*patm/R/(tc+273.15)*1e6 # E in umol/m2/s
  }
  
  ## light-limited assimilation
  a_j <- calc_assim_light_limited(gs, jmax, par_photosynth) # Aj in umol/m2/s
  a = a_j$a
  ci = a_j$ci
  vcmax = calc_vcmax_coordinated_numerical(a,ci, par_photosynth)
  
  costs = par_cost$alpha * jmax + 
    par_cost$gamma * dpsi^2 #((abs((-dpsi)/par_plant$psi50)))^2  
  
  benefit = 1 #(1+1/(par_photosynth$ca/40.53))/2
  
  dummy_costs = 0*exp(20*(-abs(dpsi/4)-abs(jmax/1))) # ONLY added near (0,0) for numerical stability. 
  
  if (opt_hypothesis == "PM"){
    ## Profit Maximisation
    out <- a*benefit - costs - dummy_costs
  } else if (opt_hypothesis == "LC"){
    ## Least Cost
    out <- -(costs+dummy_costs) / (a+1e-4)
  }
  
  if (do_optim){
    return(-out)
  } else {
    return(out)
  }
}

optimise_midterm_multi <- function(fn_profit, psi_soil, par_cost, par_photosynth, par_plant, par_env, return_all = FALSE, opt_hypothesis, gs_approximation){
  
  out_optim <- optimr::optimr(
    par       = c(logjmax=0, dpsi=1),  
    lower     = c(-10, .0001),
    upper     = c(10, 1e6),
    fn        = fn_profit,
    psi_soil  = psi_soil,
    par_cost  = par_cost,
    par_photosynth = par_photosynth,
    par_plant = par_plant,
    par_env   = par_env,
    do_optim  = TRUE,
    opt_hypothesis = opt_hypothesis,
    gs_approximation = gs_approximation,
    method    = "L-BFGS-B",
    control   = list( maxit = 500, maximize = TRUE, fnscale=1e4 )
  )
  
  out_optim$value <- -out_optim$value
  
  if (return_all){
    out_optim
  } else {
    return(out_optim$par)
  }
}


