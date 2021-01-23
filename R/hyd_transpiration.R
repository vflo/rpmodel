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
integral_P = function(dpsi, psi_soil, psi50, b, ...){
  integrate(P, psi50=psi50, b=b, lower = psi_soil, upper = (psi_soil - dpsi), ...)$value
  # -(P(psi_soil, psi50=psi50, b=b)*dpsi) # Linearized version
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
  
}

## Derivative of gs wrt dpsi (mol/m2/s/Mpa)
calc_gsprime = function(dpsi, psi_soil, par_plant, par_env){
  K = scale_conductivity(par_plant$conductivity, par_env)
  D = (par_env$vpd/par_env$patm)
  K/1.6/D*P(psi_soil-dpsi, par_plant$psi50, par_plant$b)
}


# ~~ For verification only ~~
# Numerical derivative of gs wrt dpsi (mol/m2/s/Mpa) 
calc_gsprime_numerical = function(dpsi, psi_soil, par_plant, par_env, ...){
  (calc_gs(dpsi+.01, psi_soil, par_plant, par_env, ...)-calc_gs(dpsi, psi_soil, par_plant, par_env, ...))/.01
}


#'@describeIn calc_gs Transpiration calculated as 1.6*gs*D (mol/m2/s)
calc_transpiration = function(dpsi, psi_soil, par_plant, par_env, ...){
  K = scale_conductivity(par_plant$conductivity, par_env)
  D = (par_env$vpd/par_env$patm)
  K * -integral_P(dpsi, psi_soil, par_plant$psi50, par_plant$b, ...)
}


