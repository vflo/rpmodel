#' The hydraulic p-model (HyPE), analytical version
#'
#' Calculates the carboxylation capacity, as coordinated to a given electron-transport limited assimilation rate.
#'
#' @export
pmodel_hydraulics_analytical <- function(tc, ppfd, vpd, co2, elv, fapar, kphio, psi_soil, rdark = 0, par_plant, par_cost){
  
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
  if (is.null(par_cost)){ 
    par_cost_now = list(alpha=0.1, gamma=4)
  }
  else{ 
    par_cost_now = par_cost
  }
  
  res = nleqslv::nleqslv(x = c(x=.8, dpsi=.5), 
                fn = derivatives, 
                method="Newton", 
                psi_soil= psi_soil, 
                par_photosynth = par_photosynth_now, 
                par_plant = par_plant_now, 
                par_env = par_env_now, 
                par_cost = par_cost_now ) #list(alpha=0.1, gamma=4))   
  
  x = res$x[1]
  dpsi = res$x[2]
  
  gs_unscaled = calc_gs(dpsi, psi_soil, par_plant_now, par_env_now) # In mol/m2/s
  
  # Scale by atmospheric pressure since ca is in Pa
  gs = gs_unscaled*1e6/par_photosynth_now$patm # Now in umol/m2/s/Pa
  
  Ajmax = calc_Aj_max(gs, x, par_photosynth_now)
  Jmax = calc_jmax_from_Ajmax(Ajmax, par_photosynth_now)
  
  ca = par_photosynth_now$ca
  Vcmax = Ajmax*(x*ca + par_photosynth_now$kmm)/(x*ca + 2*par_photosynth_now$gammastar)
  
  A = gs * ca*(1-x)
  
  return(list(
    jmax=Jmax,
    dpsi=dpsi,
    gs=gs_unscaled,
    a=A,
    ci=x*par_photosynth_now$ca,
    chi = x,
    chi_jmax_lim = chi_jmax_limited(par_photosynth_now, par_cost_now),
    vcmax=Vcmax,
    profit = A - par_cost_now$alpha*Jmax - par_cost_now$gamma*dpsi^2,
    niter = res$niter,
    nfcnt = res$nfcnt
  ))
  
  
}
