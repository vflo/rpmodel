#' The hydraulic p-model (HyPE), numerical version
#'
#' Calculates the carboxylation capacity, as coordinated to a given electron-transport limited assimilation rate.
#'
#' @export
pmodel_hydraulics_numerical <- function(tc, ppfd, vpd, co2, elv, fapar, kphio, psi_soil, rdark = 0, par_plant, par_cost = NULL, opt_hypothesis = "PM"){
  
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
        gamma = 1          # cost of hydraulic repair
      )
    } else if (opt_hypothesis == "LC"){
      par_cost_now = list(
        alpha = .1,        # cost of Jmax
        gamma = 0.5          # cost of hydraulic repair
      )
    }
  }
  
  lj_dps = optimise_midterm_multi(fn_profit, psi_soil = psi_soil, par_cost  = par_cost_now, par_photosynth = par_photosynth_now, par_plant = par_plant_now, par_env = par_env_now, opt_hypothesis = opt_hypothesis)
  
  profit = fn_profit(par = lj_dps, psi_soil = psi_soil, par_cost  = par_cost_now, par_photosynth = par_photosynth_now, par_plant = par_plant_now, par_env = par_env_now, opt_hypothesis = opt_hypothesis)
  
  jmax = exp(lj_dps[1])
  dpsi = lj_dps[2]
  
  gs = calc_gs(dpsi=dpsi, psi_soil=psi_soil, par_plant = par_plant_now, par_env = par_env_now)
  
  a_j = calc_assim_light_limited(gs = gs, jmax = jmax, par_photosynth = par_photosynth_now)
  a = a_j$a
  ci = a_j$ci
  
  vcmax = calc_vcmax_coordinated_numerical(a,ci, par_photosynth_now)
  
  return(list(
    jmax=jmax,
    dpsi=dpsi,
    gs=gs,
    a=a,
    ci=ci,
    chi = ci/par_photosynth_now$ca,
    vcmax=vcmax,
    profit = profit,
    chi_jmax_lim = 0
  ))
  
  
}

