#### Analytical Solver ####

calc_Aj_max = function(gs, x, par_photosynth){
  g = par_photosynth$gammastar/par_photosynth$ca
  k = par_photosynth$kmm/par_photosynth$ca
  ca = par_photosynth$ca
  d = par_photosynth$delta
  gs*ca*(1-x)*(x+2*g)/(x*(1-d)-(g+d*k)) 
}

calc_jmax_from_Ajmax = function(ajmax, par_photosynth){
  p = par_photosynth$phi0 * par_photosynth$Iabs
  4*p/((p/ajmax)^2-1)^(1/2)
}


calc_djmax_dAjmax = function(ajmax, par_photosynth){
  p = par_photosynth$phi0 * par_photosynth$Iabs
  4*p^3/ajmax^3/((p/ajmax)^2-1)^(3/2)
}

calc_dAjmax_dchi = function(gs, x, par_photosynth){
  g = par_photosynth$gammastar/par_photosynth$ca
  k = par_photosynth$kmm/par_photosynth$ca
  ca = par_photosynth$ca
  d = par_photosynth$delta
  # gs*ca * ((d*(2*g*(k + 1) + k*(2*x - 1) + x^2) + 2*g^2 + g*(2*x - 3) - x^2)/(d*(k + x) + g - x)^2)
  gs*ca * ((d*(2*g*(k + 1) + k*(2*x - 1) + x^2) - ((x-g)^2+3*g*(1-g)))/(d*(k + x) + g - x)^2)
  #gs*ca*(3*(g-1)*g/(g-x)^2 - 1)
}

calc_dAjmax_ddpsi = function(gsprime, x, par_photosynth){
  g = par_photosynth$gammastar/par_photosynth$ca
  k = par_photosynth$kmm/par_photosynth$ca
  ca = par_photosynth$ca
  d = par_photosynth$delta
  gsprime*ca*(1-x)*(x+2*g)/(x*(1-d)-(g+d*k))
}


derivatives = function(x, psi_soil, par_photosynth, par_plant, par_env, par_cost){
  X = x[1]
  dpsi = x[2]
  
  gs = calc_gs(dpsi, psi_soil, par_plant, par_env)* 1e6/par_photosynth$patm
  gsprime = calc_gsprime(dpsi, psi_soil, par_plant, par_env)* 1e6/par_photosynth$patm
  ajmax = calc_Aj_max(gs, X, par_photosynth)
  
  ca = par_photosynth$ca
  g = par_photosynth$gammastar/par_photosynth$ca
  
  dP_dx = -gs*ca - par_cost$alpha * calc_djmax_dAjmax(ajmax, par_photosynth) * calc_dAjmax_dchi(gs, X, par_photosynth)
  
  dP_ddpsi = gsprime*ca*(1-X) - par_cost$alpha * calc_djmax_dAjmax(ajmax, par_photosynth) * calc_dAjmax_ddpsi(gsprime, X, par_photosynth) - 2*par_cost$gamma*dpsi #/par_plant$psi50^2
  # cat(c(dP_dx, dP_ddpsi), "\n")
  c(dP_dx, dP_ddpsi)
}


#' Electron-transport limited \eqn{\chi} 
#'
#' Analytically calculate \eqn{\chi} in the case of strong Jmax limitation
#'
#' @param par_photosynth numeric, A list of photosynthesis parameters, which must include ...... . All concentrations must be converted to partial pressures (Pa)
#' @param par_cost numeric, A list of cost parameters
#'
#' @return \eqn{\chi}
#'
#' @export
#'
# Analytical chi in the case of strong Jmax limitation
chi_jmax_limited = function(par_photosynth, par_cost){
  g = par_photosynth$gammastar/par_photosynth$ca
  k = par_photosynth$kmm/par_photosynth$ca
  ca = par_photosynth$ca
  b = par_photosynth$delta
  a = par_cost$alpha
  
  #(g*(1-4*a) + 2*sqrt(3)*sqrt(a*(1-4*a)*(1-g)*g))/(1-4*a)
  # (2*sqrt(-a*(4*a + b - 1)*(2*b^2*g*k + 2*b^2*g + b^2*(-k^2) - b^2*k + 2*b*g^2 - 4*b*g*k - 5*b*g + b*k - 3*g^2 + 3*g)) - 4*a*b*k - 4*a*g + b^2*(-k) - b*g + b*k + g)/((b - 1)*(4*a + b - 1))
  (2*sqrt(-a*(4*a + b - 1)*(-3*g + 2*b*g - b*k)*(-1 + b + g + b*k)) - (4*a + b - 1)*(b*k + g))/((b - 1)*(4*a + b - 1))
}

