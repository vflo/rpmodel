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
fn_profit_instantaneous = function(par, jmax, vcmax, psi_soil, par_cost, par_photosynth, par_plant, par_env, gs_approximation = gs_approximation, do_optim=F){
  dpsi = par
  if (gs_approximation == "Ohm"){
    gs = calc_gs(dpsi, psi_soil, par_plant, par_env)  # gs in mol/m2/s/Mpa
    E = 1.6*gs*(par_env$vpd/par_env$patm)*1e6         # E in umol/m2/s
  } else if (gs_approximation == "PM"){
    PM_params = calc_PM_params(par_env$tc,par_env$patm, par_env$nR, par_env$LAI)
    u = par_env$u
    ustar = par_env$ustar
    R = PM_params$R
    tc = par_env$tc
    patm = par_env$patm
    dens = PM_params$air_dens
    cp =PM_params$cp
    L = PM_params$L
    pch = PM_params$pch
    C = PM_params$C
    S = PM_params$S
    Q = PM_params$Q
    vpd = par_env$vpd
    D = par_env$vpd/par_env$patm
    ga = calc_ga(par_env$u, par_env$ustar, R, tc, patm)
    gsh2o = calc_gs_PM(dpsi, psi_soil, par_plant, par_env, PM_params)
    # gs = gsh2o/1.6
    # E = (S*Q+dens*cp*vpd*ga)/(L*(S+pch*(1+ga/gsh2o)))*patm/R/(tc+273.15)*1e6 # E in umol/m2/s
    foo <- data.frame(Tair = tc, pressure = patm/1000, Rn = Q, VPD = vpd/1000, Ga_h = ga,Gs_pot = gsh2o, G=0, S=0)
    E = bigleaf::potential.ET(foo,approach="Penman-Monteith",S=0,G=0)
    foo <- data.frame(Tair = tc, pressure = patm/1000, Rn = Q, VPD = vpd/1000, LE = L*E$ET_pot, Ga_h = ga, Gs_pot = gsh2o, G=0, S=0)
    gsh2o = bigleaf::surface.conductance(foo,formulation ="Penman-Monteith" ,S=0,G=0)[['Gs_mol']]
    gs = gsh2o/1.6*1e6 
    E = E[['ET_pot']]
  }
  A = calc_assimilation_limiting(vcmax, jmax, gs, par_photosynth)$a
  profit = A - par_cost$gamma * dpsi^2  
  
  if (do_optim){
    return(-profit)
  } else {
    return(profit)
  }
}


optimise_shortterm <- function(fn_profit_inst, jmax, vcmax, psi_soil, par_cost, par_photosynth, par_plant, par_env, gs_approximation = gs_approximation, return_all = FALSE){
  
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
    gs_approximation = gs_approximation,
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
