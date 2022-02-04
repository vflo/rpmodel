#' The hydraulic p-model (HyPE), numerical version
#'
#' Instantaneous photosynthesis rates for a given Vcmax and Jmax.
#'
#' @export
pmodel_hydraulics_instantaneous <- function(tc, ppfd, vpd, co2, elv, fapar, kphio, psi_soil, rdark = 0, par_plant, par_cost = NULL, vcmax, jmax, opt_hypothesis = "PM"){
  
  p = rpmodel::calc_patm(elv)
  
  par_photosynth_now <- list(
    kmm = rpmodel::calc_kmm(tc, p),  # Why does this use std. atm pressure, and not p(z)?
    gammastar = rpmodel::calc_gammastar(tc, p),
    phi0 = kphio*rpmodel::calc_ftemp_kphio(tc),
    Iabs = ppfd*fapar,
    ca = co2*p*1e-6,  # Convert to partial pressure
    patm = p,
    delta = rdark
  )
  
  par_env_now = list(
    viscosity_water = rpmodel::calc_viscosity_h2o(tc, p),  # Needs to be imported from rpmodel.R
    density_water = rpmodel::calc_density_h2o(tc, p),  # Needs to be imported from rpmodel.R
    patm = p,
    tc = tc,
    vpd = vpd
  )
  
  par_plant_now = par_plant
  
  if (!is.null(par_cost)){
    par_cost_now = par_cost
  }
  else{
    if (opt_hypothesis == "PM"){
      par_cost_now = list(
        alpha = 0.1,       # cost of Jmax
        gamma = 4          # cost of hydraulic repair
      )
    } else if (opt_hypothesis == "LC"){
      par_cost_now = list(
        alpha = .1,        # cost of Jmax
        gamma = 2          # cost of hydraulic repair
      )
    }
  }
  
  dpsi = optimise_shortterm(fn_profit_instantaneous, jmax=jmax, vcmax=vcmax, psi_soil = psi_soil, par_cost  = par_cost_now, par_photosynth = par_photosynth_now, par_plant = par_plant_now, par_env = par_env_now)
  
  gs = calc_gs(dpsi, psi_soil, par_plant_now, par_env_now)
  
  a_l = calc_assimilation_limiting(vcmax, jmax, gs, par_photosynth_now)
  
  a = a_l$a
  ci = a_l$ci
  
  profit = a - par_cost_now$gamma*dpsi^2
  
  return(list(
    jmax=jmax,
    dpsi=dpsi,
    gs=gs,
    a=a,
    ci=ci,
    ac = a_l$ac,
    aj = a_l$aj,
    chi = ci/par_photosynth_now$ca,
    vcmax=vcmax,
    profit = profit,
    chi_jmax_lim = 0
  ))
  
  
}
