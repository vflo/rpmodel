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

## Calculation of water properties
calc_h2O_prop <- function(tc, p){
  R = 287.058 #Specific gas constant for dry air J kg-1 K-1
  # dry air density 
  air_dens = p/(R*(tc+273.15)) #kg m-3
  # dry air specific heat capacity at constant pressure
  cp = 1012 #J kg-1 K-1
  # Latent heat of water vaporization
  L =  2500.8-2.36*Tc+0.0016*tc^2-0.00006*tc^3 #J kg-1
  # Psychrometric constant
  psichro = (p*cp)/(0.622*L) #Pa K-1
  #H2O mol per Kg
  C = 18 #mol kg-1
}

## Calculation of aerodynamic conductance
calc_ga <- function(u,ustar){
  # u = wind speed m s-1
  # ustar = friction velocity m s-1
  if(is.na(ustar)){
    1/((u/ustar)+135*ustar^(-0.67)) #Thom 1972 Momentum, mass and heat exchange of vegetation
  }else{
    u/208
  }
}

#' Stomatal conductance and Transpiration using PM extension
#'
#' Calculates regulated stomatal conducatnce given the leaf water potential, 
#' plant hydraulic traits, the environment.
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
integral_P = integral_P_ana

#' #'@describeIn calc_gs Transpiration calculated as 1.6*gs*D (mol/m2/s)
#' calc_transpiration = function(dpsi, psi_soil, par_plant, par_env, ...){
#'   K = scale_conductivity(par_plant$conductivity, par_env)
#'   D = (par_env$vpd/par_env$patm)
#'   K * -integral_P(dpsi, psi_soil, par_plant$psi50, par_plant$b, ...)
#' }


