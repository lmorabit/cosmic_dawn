##################################################
##
##  A collection of functions that calculate
##  things like FWHM, FoV, resolution, etc.
##
##  based on:
##  http://www.skatelescope.org/uploaded/59513_113_Memo_Nijboer.pdf
##  with help from various websites about interferometry
##
##################################################

require('argparser')

## constants
speedoflight <- 299792458 # meters / sec
rad2deg <- 206264.8 / 60. / 60.
boltzmann <- 1.380658e-16 # erg / k
LBA_outer_diameter <- 81.34
LBA_international_diameter <- 65.0


###################################
## STATION EFFECTIVE AREA

getAeff <- function( frequency, aeff ) {
	## I would need the actual dipole locations to calculate the A_eff, which is:
	## A_eff,dipole = min{ wavelength^2./3 , pi*d^2./4. }
	## where d = the available physical area a dipole has 
	## i.e., distance to the next dipole
	## instead I'm using values from van Haarlem 2013
	freq_aeff <- c(30.,45.,60.,75.) * 1e6
	fit <- lm(aeff ~ poly(freq_aeff, 3, raw=TRUE) )
	aeff_fit <- predict(fit, data.frame(freq_aeff=frequency))
	return(aeff_fit)
}

###################################
## SEFD

getTsky <- function( wavelength, galactic_lat ){
	if (galactic_lat >= 10) Ts0 <- 60.
	Tsky <- Ts0 * wavelength^2.55
	return(Tsky)
}

getTsys <- function(frequency, Tsky){

	Freq <- c( 6.6878980892, 17.2611464968, 18.4076433121, 22.4840764331, 25.1592356688, 28.5987261146, 30.6369426752, 34.7133757962, 39.1719745223, 42.2292993631, 44.6496815287, 46.9426751592, 49.7452229299, 52.1656050955, 54.076433121, 55.6050955414, 56.6242038217, 57.898089172, 58.9171974522, 60.4458598726, 62.4840764331, 65.0318471338, 66.8152866242, 68.3439490446, 70.3821656051, 72.1656050955, 73.821656051, 75.8598726115, 78.6624203822, 80.7006369427, 82.3566878981, 84.0127388535, 85.2866242038, 86.4331210191 )
	SkyoverSys <- c( 7.9457364341, 24.0310077519, 27.519379845, 33.9147286822, 37.5968992248, 43.2170542636, 49.0310077519, 59.1085271318, 64.7286821705, 72.2868217054, 76.9379844961, 81.7829457364, 85.6589147287, 88.5658914729, 91.2790697674, 92.2480620155, 92.2480620155, 91.6666666667, 90.1162790698, 87.4031007752, 84.496124031, 81.2015503876, 75.3875968992, 69.9612403101, 62.4031007752, 56.976744186, 50.9689922481, 43.7984496124, 37.5968992248, 29.2635658915, 22.6744186047, 12.015503876, 5.4263565892, 2.9069767442 )
	freqshz <- Freq * 1e6
	interpvals <- approx(freqshz, SkyoverSys, xout=frequency, method="linear")$y / 100.
	Tsys <- Tsky / interpvals
	return(Tsys)
}

getSEFD <- function(frequency, Tsky, effective_area){
	efficiency_factor <- 1.0
	Tinstr <- getTsys(frequency, Tsky) - Tsky
	Tsys <- Tsky + Tinstr
	Ssys <- (2. * efficiency_factor * boltzmann) / effective_area * Tsys
	# convert to Jy
	Ssys <- (Ssys / 100.^2.) / 1e-23 
	return(Ssys)
}

###################################
## Sensitivity

singleStation <- function( stn_sefd, bandwidth, inttime ){
	delta_sens <- stn_sefd / sqrt(2. * bandwidth * inttime)
	return( delta_sens )
}

###################################
## Main program

getArraySens <- function( frequency, imageweight, bandwidth, inttime, ncore, nremote, nint, gal_lat, noCSCS, quiet=FALSE ){

	wavelength <- speedoflight / frequency

	## get the effective areas
	aeff_outer <- c(1559.,708.3,399.9,256.)
	aeff_international <- c(2516.,1378.,800.,512.)
	lba_outer_aeff <- getAeff( frequency, aeff_outer )
	lba_inner_aeff <- getAeff( frequency, aeff_international )
	
	## get the SEFDs
	LBA_outer_SEFD <- getSEFD( frequency, getTsky( wavelength, gal_lat ), lba_outer_aeff)
	LBA_international_SEFD <- getSEFD( frequency, getTsky( wavelength, gal_lat ), lba_inner_aeff)

	core_deltaSens <- singleStation( LBA_outer_SEFD, bandwidth, inttime )
	remote_deltaSens <- singleStation( LBA_outer_SEFD, bandwidth, inttime )
	int_deltaSens <- singleStation( LBA_international_SEFD, bandwidth, inttime )

	core_core <- (ncore * (ncore-1.)/2.) / core_deltaSens^2.
	core_remote <- ncore * nremote / (core_deltaSens*remote_deltaSens)
	core_int <- ncore * nint / (core_deltaSens*int_deltaSens)
	remote_remote <- (nremote * (nremote-1.)/2.) / remote_deltaSens^2.
	remote_int <- nremote * nint / (remote_deltaSens*int_deltaSens)
	int_int <- (nint * (nint - 1.)/2.) / int_deltaSens^2.

	if (noCSCS){
		outerArraySens <- imageweight / sqrt(2. * ( core_remote + remote_remote + core_int + remote_int + int_int ) )
	} else{
		outerArraySens <- imageweight / sqrt(2. * ( core_core + core_remote + remote_remote + core_int + remote_int + int_int ) )
	}

    if ( !quiet ) print(paste("Sensitivity is: ",format(outerArraySens[1],digits=3)," Jy",sep=""))
	return( outerArraySens[1] )

}


p <- arg_parser("Calculate array sensitivity")
p <- add_argument(p,"--weight",help="Image weighting value",default=1.0,type="double")
p <- add_argument(p,"--bandwidth",help="Bandwidth of observation",type="double")
p <- add_argument(p,"--frequency",help="Central frequency of bandwidth",type="double")
p <- add_argument(p,"--time",help="Integration time (in seconds)",type="double")
p <- add_argument(p,"--n_core",help="Number of core stations",default=24.,type="double")
p <- add_argument(p,"--n_remote",help="Number of remote stations",default=16.,type="double")
p <- add_argument(p,"--n_international",help="Number of international stations",default=9.,type="double")
p <- add_argument(p,"--latitude",help="Galactic Latitude of observation",default=30.,type="double")
p <- add_argument(p,"--noCSCS",help="Do not use CS-CS baselines in calculation (default is to use)",flag=TRUE)
argv <- parse_args(p)

mysens <- getArraySens( argv$frequency, argv$weight, argv$bandwidth, argv$time, argv$n_core, argv$n_remote, argv$n_international, argv$latitude, argv$noCSCS )

