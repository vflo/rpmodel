#' The hydraulic p-model (HyPE), numerical version
#'
#' Instantaneous photosynthesis rates for a given Vcmax and Jmax.
#'
#' @export
pmodel_hydraulics_instantaneous <- function(tc, ppfd, vpd, u, ustar, nR, co2, elv, LAI, fapar, kphio, psi_soil, rdark = 0, par_plant, par_cost = NULL, vcmax, jmax, opt_hypothesis = "PM", gs_approximation = "PM"){
  
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
    vpd = vpd,
    u = u,
    ustar = ustar,
    nR = nR,
    LAI = LAI
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
  
  dpsi = optimise_shortterm(fn_profit_instantaneous, jmax=jmax, vcmax=vcmax, psi_soil = psi_soil, par_cost  = par_cost_now,
                            par_photosynth = par_photosynth_now, par_plant = par_plant_now, par_env = par_env_now, gs_approximation = gs_approximation)
  
  if (gs_approximation == "Ohm"){
    gs = calc_gs(dpsi, psi_soil, par_plant, par_env_now)  # gs in mol/m2/s/Mpa
    E = 1.6*gs*(par_env_now$vpd/par_env_now$patm)*1e6         # E in umol/m2/s
  } else if (gs_approximation == "PM"){
    PM_params = calc_PM_params(par_env_now$tc,par_env_now$patm, par_env_now$nR, par_env_now$LAI)
    u = par_env_now$u
    ustar = par_env_now$ustar
    R = PM_params$R
    tc = par_env_now$tc
    patm = par_env_now$patm
    dens = PM_params$air_dens
    cp =PM_params$cp
    L = PM_params$L
    pch = PM_params$pch
    C = PM_params$C
    S = PM_params$S
    Q = PM_params$Q
    vpd = par_env_now$vpd
    D = vpd/patm
    ga = calc_ga(u, ustar, R, tc, patm)
    gs = calc_gs_PM(dpsi, psi_soil, par_plant_now, par_env_now, PM_params)
    # E = C*(S*Q+dens*cp*vpd*ga*R*tc/patm)/(L*(S+pch*(1+ga/(1.6*gs))))*1e6 # E in umol/m2/s
    E = (S*Q+dens*cp*vpd*ga)/(L*(S+pch*(1+(ga*patm/(R*(tc+273.15)))/(1.6*gs))))*1e6 # E in umol/m2/s
  }
  
  a_l = calc_assimilation_limiting(vcmax, jmax, gs, par_photosynth_now)
  
  a = a_l$a
  ci = a_l$ci
  
  profit = a - par_cost_now$gamma*dpsi^2
  
  if (gs_approximation == "Ohm"){
    return(list(
      jmax=jmax,
      dpsi=dpsi,
      gs=gs,
      E = E,
      a=a,
      ci=ci,
      ac = a_l$ac,
      aj = a_l$aj,
      chi = ci/par_photosynth_now$ca,
      vcmax=vcmax,
      profit = profit,
      chi_jmax_lim = 0
    ))
  } else if (gs_approximation == "PM"){
    return(list(
      jmax=jmax,
      dpsi=dpsi,
      gs=gs,
      ga = ga,
      E = E,
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

  
  
}
