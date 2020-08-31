#### Photosynthesis Rates ####

#' Carboxylation-limited assimilation rate 
#'
#' Calculates the carboxylation limited CO2 assilimation rate
#'
#' @param gs numeric, Stomatal conductance in mol/m2/s, for e.g. as calculated from \code{\link{calc_gs}} 
#' @param vcmax numeric, Carboxylation capacity (umol/m2/s)
#' @param par_photosynth numeric, A list of photosynthesis parameters, which must include ...... . All concentrations must be converted to partial pressures (Pa)
#'
#' @return Rubisco-limited assimilation rate (umol/m2/s)
#'
#' @export
#'
# Calculate Assimilation rate (umol/m2/s), given
# - gs (mol/m2/s)
# - vcmax (umol/m2/s)
# - Photosynthesis parameters with concentrations converted to partial pressures (Pa)
calc_assim_rubisco_limited <- function(gs, vcmax, par_photosynth){
  
  ca = par_photosynth$ca             # ca is in Pa
  gs = gs * 1e6/par_photosynth$patm  # convert to umol/m2/s/Pa
  
  d = par_photosynth$delta
  ## with
  A <- -1.0 * gs
  B <- gs * ca - gs * par_photosynth$kmm - vcmax*(1-d)
  C <- gs * ca * par_photosynth$kmm + vcmax * (par_photosynth$gammastar + par_photosynth$kmm*d)
  
  ci <- QUADM(A, B, C)
  # a_c <- vcmax * (ci - par$gammastar) / (ci + par$kmm)
  a_c <- gs*(ca-ci) 
  
  return(list(a=a_c, ci=ci))
}

#' Electron-transport-limited assimilation rate 
#'
#' Calculates the electron transport limited CO2 assilimation rate
#'
#' @param gs numeric, Stomatal conductance in mol/m2/s, for e.g. as calculated from \code{\link{calc_gs}} 
#' @param jmax numeric, Electron transport capacity (umol/m2/s)
#' @param par_photosynth numeric, A list of photosynthesis parameters, which must include ...... . All concentrations must be converted to partial pressures (Pa)
#'
#' @return Rubisco-limited assimilation rate (umol/m2/s)
#'
#' @export
#'
# Calculate Assimilation rate (umol/m2/s), given
# - gs (mol/m2/s)
# - jmax (umol/m2/s)
# - Photosynthesis parameters with concentrations converted to partial pressures (Pa)
# - Iabs (umol/m2/s)
calc_assim_light_limited <- function(gs, jmax, par_photosynth){
  
  ## Only light is limiting
  ## Solve Eq. system
  ## A = gs (ca- ci)
  ## A = phi0 * Iabs * jlim * (ci - gammastar)/(ci + 2*gamma_star)
  
  ## This leads to a quadratic equation:
  ## A * ci^2 + B * ci + C  = 0
  ## 0 = a + b*x + c*x^2
  
  ca = par_photosynth$ca             # ca is in Pa
  gs = gs * 1e6/par_photosynth$patm  # convert to umol/m2/s/Pa
  
  phi0iabs = par_photosynth$phi0 * par_photosynth$Iabs
  jlim = phi0iabs / sqrt(1+ (4*phi0iabs/jmax)^2)
  
  d = par_photosynth$delta 
  ## with
  A <- -1.0 * gs
  B <- gs * ca - gs * 2 * par_photosynth$gammastar - jlim*(1-d)
  C <- gs * ca * 2*par_photosynth$gammastar + jlim * (par_photosynth$gammastar + d*par_photosynth$kmm)
  
  ci <- QUADM(A, B, C)
  aj <- gs*(ca-ci)
  # vcmax_pot <- a*(ci + par$kmm)/(ci - par$gammastar)
  
  return(list(a=aj, ci=ci))
}

#' Limiting assimilation rate 
#'
#' Calculates the limiting assimilation rate, the minimum of the 
#' carboxylation limited and electron-transport-limited rates.
#'
#' @param gs numeric, Stomatal conductance in mol/m2/s, for e.g. as calculated from \code{\link{calc_gs}} 
#' @param vcmax numeric, Carboxylation capacity (umol/m2/s)
#' @param jmax numeric, Electron transport capacity (umol/m2/s)
#' @param par_photosynth numeric, A list of photosynthesis parameters, which must include ...... . All concentrations must be converted to partial pressures (Pa)
#'
#' @return Rubisco-limited assimilation rate (umol/m2/s)
#'
#' @export
#'
calc_assimilation_limiting = function(vcmax, jmax, gs, par_photosynth){
  # gs = calc_gs(dpsi, psi_soil, par_plant = par_plant, par_env = par_env)
  
  # We need not employ numerical root-finding. calculate chi independently assuming Ac and Aj, and bigger of the two will be the limiting one. Accordingly return Ac or Aj
  
  Ac = calc_assim_rubisco_limited(gs = gs, vcmax = vcmax, par_photosynth = par_photosynth)
  Aj = calc_assim_light_limited(gs = gs, jmax = jmax, par_photosynth = par_photosynth)
  
  if (Ac$ci > Aj$ci ){
    Ac
  } else {
    Aj
  }
}
