library('geosphere')

get_baselines <- function( station_file = 'lofar_hba_baselines.list', make_plot=FALSE ){

	## read in the baselines
	A <- read.table( station_file, stringsAsFactors=FALSE, header=FALSE )

	lofar_stn <- A$V2
	## convert lat/long to degrees
	lofar_lat <- numeric()
	for ( l in A$V6 ){
        	tmp <- strsplit( l, ".", fixed=TRUE )
	        lofar_lat <- c( lofar_lat, as.numeric( tmp[[1]][1] ) + ( as.numeric( tmp[[1]][2] ) + as.numeric( paste( tmp[[1]][3], tmp[[1]][4], sep='.' ) )/60. )/60. )
	}

	lofar_long <- numeric()
	for ( l in A$V5 ){
        	tmp <- strsplit( l, ".", fixed=TRUE )
	        lofar_long <- c( lofar_long, as.numeric( tmp[[1]][1] ) + ( as.numeric( tmp[[1]][2] ) + as.numeric( paste( tmp[[1]][3], tmp[[1]][4], sep='.' ) )/60. )/60. )
	}


	baselines <- numeric()

	for ( stn in lofar_stn ){
        	p1 <- c( lofar_long[ which(lofar_stn == stn) ], lofar_lat[ which(lofar_stn == stn) ] )
	        baseline <- distCosine( p1, cbind(lofar_long, lofar_lat) )  # in meters
        	baselines <- cbind( baselines, baseline )
	}

	if ( make_plot ){

		pdf( 'baselines.pdf' )
		layout( rbind( cbind( 1, 2 ), cbind( 3, 4 ) ) )
		## dutch stations
		my_stns <- grep( 'CS', lofar_stn )
		my_bls <- baselines[ my_stns, my_stns ]
		plot_bls <- unique( my_bls )/1e3 ## convert to km
		hist( plot_bls[which(plot_bls > 0)], breaks=20, xlab="Baseline Length [km]", ylab="Number of baselines", main="Core Stations", cex=1.2 )
		box( which="plot", lwd=2 )
		## remote stations
		tmp <- grep( 'RS', lofar_stn )
		my_stns <- c( my_stns, tmp )
		my_bls <- baselines[ my_stns, my_stns ]
		plot_bls <- unique( my_bls )/1e3 ## convert to km
		hist( plot_bls[which(plot_bls > 0)], breaks=20, xlab="Baseline Length [km]", ylab="Number of baselines", main="All Dutch Stations", cex=1.2 )
		box( which="plot", lwd=2 )
		dev.off()
	}
	return( list( antennas=lofar_stn, baselines=baselines ) )

}
