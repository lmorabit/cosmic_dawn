speedoflight <<- 2.99792458e8  ## m/s
hplanck <<- 6.6260755e-27 ## erg/s
gravitational_constant <<- 6.67259e-8 ## cm^3/g/s^2
electron_charge <<- 4.8032068e-10 ## esu
mass_electron <<- 9.1093897e-28  ## grams
mass_proton <<- 1.6726231e-24  ## grams
mass_neutron <<- 1.6749286e-24  ## grams
mass_hydrogen <<- 1.6733e-24  ## grams
atomic_mass_unit <<- 1.6605402e-24  ## grams
boltzmann_constant <<- 1.380658e-16 ## erg/K
electron_volt <<- 1.6021772e-12 ## erg
thomson_xc <<- 6.6524e-25 ## cm^2
solar_mass <<- 1.99e33 ## grams
charge_electron <<- 4.8032068e-10 ## esu
vacuum_permitivity <<- 8.854187817e-12 ## farad/meter
electron_radius <<- 2.8179e-15 ## meters

## unit conversion
rad2arcsec <<- 206264.806247


## cosmological info Planck 2016
#library('astro')
H_0 <<- 67.8
omega_m <<- 0.308
omega_lam <<- 0.692

l_angdist <- function( c = speedoflight, H = H_0, M = omega_m, L = omega_lam, ... ){
	val <- angdist( ... )
	return( val )
}
