
## Purpose

`rpmodel` provides an implementation of the P-model (Prentice et al., 2014; Wang et al., 2017; Stocker et al., 2020) for predicting acclimated photosynthetic parameters, assimilation, and dark respiration rates as a function of the environment. The main function is `rpmodel()` which returns a list of variables that are mutually consistent within the theory of the P-model (see [Usage](./articles/usage.html) ). Further functions used within `rpmodel()` are also provided through the package.

It also provides both analytical and numerical implementations of the hydraulic extension of the P-model (Joshi et al., XXX). See theory [here](https://rpubs.com/jaideep777/pmodel_hydraulics). 

## Usage

This loads the `rpmodel` package and executes the `rpmodel()` function without $J_{\text{max}}$ limitation (argument `method_jmaxlim = "none"`), and with a temperature-independent quantum yield efficiency (argument `do_ftemp_kphio = FALSE`):
```r
library(rpmodel)
out_pmodel <- rpmodel( 
  tc             = 20           # temperature, deg C
  vpd            = 1000         # Pa,
  co2            = 400          # ppm,
  elv            = 0            # m.a.s.l.,
  kphio          = 0.05         # quantum yield efficiency,
  beta           = 146,         # unit cost ratio a/b,
  fapar          = 1            # fraction  ,
  ppfd           = 300          # mol/m2/d,
  method_optci   = "prentice14",
  method_jmaxlim = "none",
  do_ftemp_kphio = FALSE 
  )
print( out_pmodel )
```

To use the hydraulic version, first make a list of the 3 required plant hydraulic traits and then run the model.


```r
plant_traits = list(
  # legacy parameters, do not change, will be deleted in future versions
  Ks0=1e-12,              # m2 
  v_huber=1e-4,           #	
  height=10,              # m
  # hydraulic traits
  conductivity_scalar=3,  # Leaf conductivity (x 10^-15 m2) 
  psi50 = -2,             # Leaf P50 (Mpa)
  b=2                     # Slope of leaf vulnerability curve 
)

out_pmodel_hydraulic <- pmodel_hydraulics_numerical( 
  tc             = 20,           # temperature, deg C
  vpd            = 1000,         # Pa,
  co2            = 400,          # ppm,
  elv            = 0,            # m.a.s.l.,
  psi_soil       = -0.5,         # soil water potential (Mpa) 
  kphio          = 0.05,         # quantum yield efficiency,
  fapar          = 1,            # fraction  ,
  ppfd           = 1000,         # umol/m2/s, # Note difference in units compared to classic model
  par_plant      = plant_traits, # plant hydraulic traits
  )
print( out_pmodel_hydraulic )
```

## Installation

### Stable release
`rpmodel` is available on CRAN [here](https://CRAN.R-project.org/package=rpmodel).

### Development release

Hydraulic pmodel development release is available at [github/jaideep777](https://github.com/jaideep777/rpmodel)

To install and load the latest version of the rpmodel package (development release, not yet on CRAN) run the following command in your R terminal: 
```r
if(!require(devtools)){install.packages(devtools)}
devtools::install_github( "jaideep777/rpmodel", ref="hydraulics", build_vignettes = TRUE )
library(rpmodel)
```

## Author and contact

Benjamin Stocker
benjamin.stocker@gmail.com

Hydraulics version:
Jaideep Joshi
jaideep777@gmail.com

## References

Stocker, B. D., Wang, H., Smith, N. G., Harrison, S. P., Keenan, T. F., Sandoval, D., Davis, T., and Prentice, I. C.: P-model v1.0: an optimality-based light use efficiency model for simulating ecosystem gross primary production, Geosci. Model Dev., 13, 1545–1581, https://doi.org/10.5194/gmd-13-1545-2020, 2020.

Wang, H., Prentice, I. C., Keenan, T. F., Davis, T. W., Wright, I. J., Cornwell, W. K.,Evans, B. J., and Peng, C.:  Towards a universal model for carbon dioxide uptake by plants, Nat Plants, 3, 734–741, 2017.

Prentice,  I. C.,  Dong,  N.,  Gleason,  S. M.,  Maire,  V.,  and Wright,  I. J.:  Balancingthe costs of carbon gain and water transport:  testing a new theoretical frameworkfor  plant  functional  ecology, Ecology  Letters,  17,  82–91, 10.1111/ele.12211, 2014.

## Acknowledgement

This project was funded by Marie Sklodowska-Curie fellowship H2020-MSCA-IF-2015, project FIBER, grant number 701329.

Hydraulic pmodel development was funded by Marie Sklodowska-Curie fellowship H2020-MSCA-IF-2019, project PlantFATE.

