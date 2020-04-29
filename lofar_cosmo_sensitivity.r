source('lofar_sensitivity.r')
source('get_baselines.r')


plot_lofar_sensitivity_vs_time <- function( n_hrs=8, sens_norm=1000, outfile='LBA_sens.pdf'  ){

    ############ PLOT OF LOFAR SENSITIVITY VS TIME
    time_sec <- seq( 10, 1e7 )
    sens <- 1./sqrt( time_sec ) * 1e6 ## convert to microjy

    ## normalize 
    time_nhrs <- n_hrs * 60 * 60
    sens_nhrs <- sens[ which( time_sec == time_nhrs ) ]
    sens_norm <- sens_norm / sens_nhrs
    sens <- sens_norm * sens

    ## plot sensitivity vs. number of observations
    the_obs <- seq( time_nhrs, 1e7, time_nhrs )
    the_sens <- sens[ which( time_sec %in% the_obs ) ]

    pdf( outfile )
    lplot( seq( 1, length(the_obs)), the_sens, x_lab=paste('Number of',n_hrs,'hour Observations'), y_lab='Sensitivity [Jy]', log='y', type='l' )
    dev.off()
}

plot_21cm_vs_redshift <- function( outfile='H21_Plot1.pdf' ){

    ############ PLOT OF 21 CM FREQ VS. REDSHIFT
    redshift <- seq( 0, 60 )
    rest_freq <- 1420.4e6 
    obs_freq <- rest_freq / ( 1 + redshift )

    pdf( outfile )
    lplot( redshift, obs_freq/1e6, x_lab='redshift', y_lab='Observed Frequency [MHz]', type='l', xlim=c(10,60), ylim=c(0,250) )
    dev.off()

    mycols <- viridis( 20 )


    pdf( gsub( 'Plot1', 'Plot2', outfile ) )
    lplot( redshift, obs_freq/1e6, x_lab='redshift', y_lab='Observed Frequency [MHz]', type='l', xlim=c(10,60), ylim=c(0,250) )
    ## get area for polygon plotting
    area <- par('usr')
    xvec <- c( 30, 30, 60, 60 )
    yvec <- c( area[3], area[4], area[4], area[3] )
    polygon( xvec, yvec, col=tcol(mycols[8],50) )
    dev.off()


    pdf( gsub( 'Plot1', 'Plot3', outfile ) )
    lplot( redshift, obs_freq/1e6, x_lab='redshift', y_lab='Observed Frequency [MHz]', type='l', xlim=c(10,60), ylim=c(0,250) )
    ## get area for polygon plotting
    area <- par('usr')
    ## plasma frequency
    xvec <- c( area[1], area[1], area[2], area[2] )
    yvec <- c( area[3], 15, 15, area[3] )
    polygon( xvec, yvec, col=tcol('black',100) )
    yvec <- c( area[3], 10, 10, area[3] )
    polygon( xvec, yvec, col=tcol('black',100) )
    text( 15.5, 3, labels='Plasma Frequency', col='white' )
    ## frequency range
    xvec <- c( 30, 30, 60, 60 )
    yvec <- c( 15, area[4], area[4], 15 )
    polygon( xvec, yvec, col=tcol(mycols[8],50) )
    dev.off()

    pdf( gsub( 'Plot1', 'Plot4', outfile ) )
    lplot( redshift, obs_freq/1e6, x_lab='redshift', y_lab='Observed Frequency [MHz]', type='l', xlim=c(10,60), ylim=c(0,250) )
    ## get area for polygon plotting
    area <- par('usr')
    ## plasma frequency
    xvec <- c( area[1], area[1], area[2], area[2] )
    yvec <- c( area[3], 15, 15, area[3] )
    polygon( xvec, yvec, col=tcol('black',100) )
    yvec <- c( area[3], 10, 10, area[3] )
    polygon( xvec, yvec, col=tcol('black',100) )
    text( 15.5, 3, labels='Plasma Frequency', col='white' )
    ## frequency range
    xvec <- c( 30, 30, 60, 60 )
    yvec <- c( 15, area[4], area[4], 15 )
    polygon( xvec, yvec, col=tcol(mycols[8],50) )
    ## lofar coverage
    xvec <- c( area[1], area[1], area[2], area[2] )
    yvec <- c( 15, 63, 63, 15 )
    polygon( xvec, yvec, col=tcol(mycols[16],50) )/data/lofar/morabito/LC6_013
    text( 14.5, 20, labels='Ideal Conditions', col=mycols[13] )
    dev.off()

    pdf( gsub( 'Plot1', 'Plot5', outfile ) )
    lplot( redshift, obs_freq/1e6, x_lab='redshift', y_lab='Observed Frequency [MHz]', type='l', xlim=c(10,60), ylim=c(0,250) )
    ## get area for polygon plotting
    area <- par('usr')
    ## plasma frequency
    xvec <- c( area[1], area[1], area[2], area[2] )
    yvec <- c( area[3], 15, 15, area[3] )
    polygon( xvec, yvec, col=tcol('black',100) )
    yvec <- c( area[3], 10, 10, area[3] )
    polygon( xvec, yvec, col=tcol('black',100) )
    text( 15.5, 3, labels='Plasma Frequency', col='white' )
    xvec <- c( 30, 30, 60, 60 )
    yvec <- c( 15, area[4], area[4], 15 )
    polygon( xvec, yvec, col=tcol(mycols[8],50) )
    xvec <- c( area[1], area[1], area[2], area[2] )
    yvec <- c( 30, 78, 78, 30 )
    polygon( xvec, yvec, col=tcol(mycols[16],50) )
    text( 15.85, 35, labels='Realistic Conditions', col=mycols[13] )
    dev.off()

}


calculate_redshift_and_size_resolution <- function( cloud_size=1., target_z=30 ){

    ## cloud size must be in Mpc
    rest_freq <- 1420.4e6 
   
    ## first calculate the angular scale in Mpc
    D_A <- angdist( z = target_z, H = H_0, M = omega_m, L=omega_lam )
    angular_size <- cloud_size / D_A * rad2arcsec ## convert to arcsec
    
    ## calculate redshift resolution
    redshift_sequence <- seq( 0, 0.2, 0.001 )
    redshift_array <- c()
    for ( ii in seq(1,length(target_z) ) ) redshift_array <- rbind( redshift_array, target_z[ii]+redshift_sequence )
    redshift_codist_vec <- cosdistCoDist( z = redshift_array )
    ## reform the array to the right shape
    redshift_codist <- array( redshift_codist_vec, dim=dim(redshift_array) )
    redshift_diff <- redshift_codist - cosdistCoDist( z = target_z )
    ## and get rid of the first column which is just the target_z
    redshift_diff <- redshift_diff[,2:length(redshift_sequence)]
    
    delta_nu <- c()
    for ( ii in seq(1,length(target_z)) ){
        best_index <- which( abs( redshift_diff[ii,] - cloud_size ) == min( abs( redshift_diff[ii,] - cloud_size ) ) )
        ## add 1 back on to account for the fact that redshift_diff is one column smaller
        delta_z <- redshift_sequence[ best_index + 1 ] 
        ## convert to frequency range
        ## in Hz
        d_nu <- rest_freq * ( delta_z ) / ( ( 1 + target_z[ii] )*( 1 + target_z[ii] ) )
        delta_nu <- c( delta_nu, d_nu )
    }

    result <- list( size_resolution=angular_size, frequency_resolution=delta_nu )
    return( result )

}


## calculate the baseline lengths of the array ;; this will set the resolution and the sensitivity
bl_info <- get_baselines()
## set the DE601-DE605 baseline to Inf
de601_ind <- which( grepl( 'DE601', bl_info$antennas ) )
de605_ind <- which( grepl( 'DE605', bl_info$antennas ) )
bl_info$baselines[de601_ind,de605_ind] <- Inf
bl_info$baselines[de605_ind,de601_ind] <- Inf

find_n_stations <- function( desired_res_arcsec, desired_frequency, my_bl_info, res_tol=0.5 ){

    ## convert baseline length to resolution
    my_bl_res <- ( speedoflight / desired_frequency ) / my_bl_info$baselines * rad2arcsec
    ## find the difference to the desired resolution
    res_diff <- my_bl_res - desired_res_arcsec

    ## check if anything is within the res_tol
    within_tol <- which( abs(res_diff) < res_tol, arr.ind=TRUE )
    if ( dim(within_tol)[1] > 0 ){
        ## there is a match to resolution!
        ## find where the resolution is *higher*
        higher_res <- which( res_diff < 0, arr.ind=TRUE )
        ## get list of telescopes in higher res 
        unnecessary_tels <- higher_res[,1]
        ant_count <- c()
        for ( ii in seq(1,length(my_bl_info$antennas) ) ) ant_count <- c( ant_count, length( which( unnecessary_tels == ii ) ) )
        valid_ant_index <- which( ant_count < length(my_bl_info$antennas)-1 )
        valid_antennas <- my_bl_info$antennas[valid_ant_index]

        ## count ncore, nremote, nint
        ncore <- length( which( grepl( 'CS', valid_antennas ) ) )
        nremote <- length( which( grepl( 'RS', valid_antennas ) ) )
        nint <- length( valid_antennas ) - ncore - nremote
    } else {
        ## there is not a match.
        ncore <- 0
        nremote <- 0
        nint <- 0
    } 

    return( list( ncore=ncore, nremote=nremote, nint=nint ) )
}

find_n_stations_ska <- function( desired_res_arcsec, desired_frequency, my_bl_info, res_tol=0.5 ){

    ## convert baseline length to resolution
    my_bl_res <- ( speedoflight / desired_frequency ) / my_bl_info$baselines * rad2arcsec
    ## find the difference to the desired resolution
    res_diff <- my_bl_res - desired_res_arcsec

    ## check if anything is within the res_tol
    within_tol <- which( abs(res_diff) < res_tol, arr.ind=TRUE )
    if ( dim(within_tol)[1] > 0 ){
        ## there is a match to resolution!
        ## find where the resolution is *higher*
        higher_res <- which( res_diff < 0, arr.ind=TRUE )
        ## get list of telescopes in higher res 
        unnecessary_tels <- higher_res[,1]
        ant_count <- c()
        for ( ii in seq(1,length(my_bl_info$antennas) ) ) ant_count <- c( ant_count, length( which( unnecessary_tels == ii ) ) )
        valid_ant_index <- which( ant_count < length(my_bl_info$antennas)-1 )
        valid_antennas <- my_bl_info$antennas[valid_ant_index]

        ntels <- length( valid_antennas )

    } else {
        ntels <- 0
    } 

    return( ntels )
}


