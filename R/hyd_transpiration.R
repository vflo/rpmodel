calc_conductivity_m = function(sapwood_perm, hv, height){
  sapwood_perm*hv/height
}

## Returns conductivity in mol/m2/s/Mpa
scale_conductivity = function(K, par_env){
  # Flow rate in m3/m2/s/Pa
  K2 = K / par_env$viscosity_water
  
  # Flow rate in mol/m2/s/Pa
  mol_h20_per_kg_h20 = 55.5
  K3 = K2 * par_env$density_water * mol_h20_per_kg_h20
  
  # Flow rate in mol/m2/s/Mpa
  K4 = K3*1e6
  
  return(K4)  
}


## Integral of the vulnerabiity curve
## integral of P(psi)dpsi  (Mpa)
integral_P_ana = function(dpsi, psi_soil, psi50, b, ...){
  ps = psi_soil/psi50
  pl = (psi_soil-dpsi)/psi50
  l2 = log(2)
  -(psi50/b)*(l2^(-1/b))*(expint::gammainc(a = 1/b, x = l2*pl^b)-expint::gammainc(a = 1/b, x = l2*ps^b))
}

integral_P_num = function(dpsi, psi_soil, psi50, b, ...){
  integrate(P, psi50=psi50, b=b, lower = psi_soil, upper = (psi_soil - dpsi), ...)$value
  # -(P(psi_soil, psi50=psi50, b=b)*dpsi) # Linearized version
}

integral_P_approx = function(dpsi, psi_soil, psi50, b, ...){
  -dpsi * P(psi_soil-dpsi/2, psi50, b)
}



integral_P = integral_P_ana

## Calculation of water parameters
calc_PM_params <- function(tc, p, nR, LAI){
  #tc air temperature ºC
  #p atmospheric pressure Pa
  #nR surface net rariation J s-1 m-2soil
  #LAI leaf area index m2leaf m-2soil
  R =  8.31446261815324 #universal gas constant m3 Pa K−1 mol−1
  dry_air_mol = 0.0289647 #kg mol-1
  R_dry_air = R/dry_air_mol # specific dry air gas constant J kg-1 K-1
  h2o_mol_mass = 0.01801528 # kg water mol-1
  R_water = R/h2o_mol_mass # specific water vapour gas constant J kg-1 K-1
  # dry air density 
  air_dens = p/(R_dry_air*(tc+273.15)) #kg m-3
  # dry air specific heat capacity at constant pressure
  cp = 1012 #J kg-1 K-1
  # Latent heat of water vaporization
  L =  (2500.8-2.36*tc+0.0016*tc^2-0.00006*tc^3)*1000 #J kg-1
  # psychrometric constant
  pch = (p*cp)/(0.622*L) #Pa K-1
  # H2O mol Kg-1
  C = 1/h2o_mol_mass #mol kg-1
  # slope of the curve relating saturation vapour pressure to temperature
  S = 4098*(611.1495*exp(17.27*tc/(tc+237.3)))/((237.3+tc)^2) #FAO Pa K-1
  Q = nR*(1-exp(-0.5*LAI))/LAI #leaf available enery J s-1 m-2leaf
  if(Q<0){Q=0}
  
  df = data.frame(R, dry_air_mol, R_dry_air, h2o_mol_mass, R_water, air_dens, cp, L, pch, C, S, Q)
  return(df)
}

## Calculation of aerodynamic conductance
calc_ga <- function(u,ustar,R,tc,p){
  # u = wind speed m s-1
  # ustar = friction velocity m s-1
  # ga = aerodynamic conductance m s-1
  if(!is.na(ustar)){
    (1/((u/ustar)+135*ustar^(-0.67))) #Thom 1972 Momentum, mass and heat exchange of vegetation
  }else{
    (u/208) #Allen et al 1998 crop evapotranspiration-guidelines for computing crop water requirements
  }
}

#' Stomatal conductance and Transpiration using PM extension
#'
#' Calculates regulated stomatal conducatnce given the leaf water potential, 
#' plant hydraulic traits, the environment.
#'
#' @param dpsi numeric, soil-to-leaf water potential difference (\eqn{\psi_s-\psi_l}), Mpa
#' @param psi_soil numeric, soil water potential, Mpa
#' @param par_plant numeric, A list of plant hydraulic parameters, which must include LAI 
#' @param par_env numeric, A list of environmental parameters which must include wind speed (u), 
#' friction velocity (ustar), net radiation (nR), viscosity_water, density_water, atmospheric pressure (patm), 
#' air temperature (tc), vapour pressure deficit (vpd)
#'
#' @return numeric, canopy stomatal conductance in molco2 m-2 s-1
#'
#'
#' @export
#'
calc_gs_PM = function(dpsi, psi_soil, par_plant, par_env, PM_params, ...){
  # Does the column u exist in par_env? If not, stop the proces and show error
  if (is.na(par_env$u)|is.null(par_env$u)) {
    stop("Wind speed (u) must be provided in par_env variable", call. = FALSE)
  }
  K = scale_conductivity(par_plant$conductivity, par_env)
  vpd = par_env$vpd
  patm = par_env$patm
  D = (vpd/patm)
  u = par_env$u
  ustar = par_env$ustar
  R = PM_params$R
  tc = par_env$tc
  dens = PM_params$air_dens
  cp =PM_params$cp
  L = PM_params$L
  pch = PM_params$pch
  C = PM_params$C
  S = PM_params$S
  Q = PM_params$Q
  ga = calc_ga(u, ustar, R, tc, patm)
  
  divid <- (S*Q+dens*cp*vpd*ga)/(L*(-K*integral_P(dpsi, psi_soil, par_plant$psi50, par_plant$b, ...))) -S -pch
  
  pch*ga/(divid)*patm/R/(tc+273.15) #Return Gs in molH20 m-2leaf  s-1
}

#' Stomatal conductance and Transpiration
#'
#' Calculates regulated stomatal conducatnce given the leaf water potential, 
#' plant hydraulic traits, and the environment.
#'
#' @param dpsi numeric, soil-to-leaf water potential difference (\eqn{\psi_s-\psi_l}), Mpa
#' @param psi_soil numeric, soil water potential, Mpa
#' @param par_plant numeric, A list of plant hydraulic parameters, which must include ...... 
#' @param par_env numeric, A list of environmental parameters
#'
#' @return numeric, stomatal conductance in mol/m2/s
#'
#'
#' @export
#'
calc_gs = function(dpsi, psi_soil, par_plant, par_env, ...){
  K = scale_conductivity(par_plant$conductivity, par_env)
  D = (par_env$vpd/par_env$patm)
  K/1.6/D * -integral_P(dpsi, psi_soil, par_plant$psi50, par_plant$b, ...)
  # # papprox = P(psi_soil-dpsi/2, par_plant$psi50, par_plant$b)
  # papprox = P(psi_soil, par_plant$psi50, par_plant$b)-Pprime(psi_soil, par_plant$psi50, par_plant$b)*dpsi/2.5
  # K/1.6/D * dpsi * papprox
  
}

## Derivative of gs wrt dpsi (mol/m2/s/Mpa)
calc_gsprime_analytical = function(dpsi, psi_soil, par_plant, par_env){
  K = scale_conductivity(par_plant$conductivity, par_env)
  D = (par_env$vpd/par_env$patm)
  K/1.6/D*P(psi_soil-dpsi, par_plant$psi50, par_plant$b)
}

## Derivative of gs wrt dpsi (mol/m2/s/Mpa)
calc_gsprime_approx = function(dpsi, psi_soil, par_plant, par_env){
  K = scale_conductivity(par_plant$conductivity, par_env)
  D = (par_env$vpd/par_env$patm)
  K/1.6/D*(P(psi_soil-dpsi/2, par_plant$psi50, par_plant$b) - Pprime(psi_soil-dpsi/2, par_plant$psi50, par_plant$b)*dpsi/2)
}


# ~~ For verification only ~~
# Numerical derivative of gs wrt dpsi (mol/m2/s/Mpa) 
calc_gsprime_numerical = function(dpsi, psi_soil, par_plant, par_env, ...){
  (calc_gs(dpsi+.01, psi_soil, par_plant, par_env, ...)-calc_gs(dpsi, psi_soil, par_plant, par_env, ...))/.01
}

calc_gsprime = calc_gsprime_analytical

#' #'@describeIn calc_gs Transpiration calculated as 1.6*gs*D (mol/m2/s)
#' calc_transpiration = function(dpsi, psi_soil, par_plant, par_env, ...){
#'   K = scale_conductivity(par_plant$conductivity, par_env)
#'   D = (par_env$vpd/par_env$patm)
#'   K * -integral_P(dpsi, psi_soil, par_plant$psi50, par_plant$b, ...)
#' }


