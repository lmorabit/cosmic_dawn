## source the files I want
source('~/software/cosmic_dawn/physical_constants.r')
library('pracma')
library('plotrix')
source('~/software/cosmic_dawn/lofar_cosmo_sensitivity.r')


library('viridis')
library('FITSio')
library('celestial')
library('fields')

## wrapper to make nice plots
lplot <- function( x, y, x_lab='', y_lab='', margins=c(5,5,2,2), leg=TRUE, xaxon=TRUE, yaxon=TRUE, axcex=1.25, dobox=TRUE, ... ){

        par( mar=margins )
        plot( x, y, ..., axes=FALSE, xlab='', ylab='' )
    if ( xaxon ){
        axis( 1, cex.axis=axcex )
        mtext( x_lab, side=1, line=3, cex=axcex )
    }
    if ( yaxon ){
            axis( 2, cex.axis=axcex )
            mtext( y_lab, side=2, line=3, cex=axcex )
    }
        if ( dobox ) box( which='plot', lwd=1.5 )
        if ( !leg ) par( mar=c(5.1,4.1,4.1,2.1) )

}

tcol <- function(color, trans = 100) {

  if (length(color) != length(trans) &
        !any(c(length(color), length(trans)) == 1))
    stop('Vector lengths not correct')
  if (length(color) == 1 & length(trans) > 1)
    color <- rep(color, length(trans))
  if (length(trans) == 1 & length(color) > 1)
    trans <- rep(trans, length(color))

  res <- paste0('#', apply(apply(rbind(col2rgb(color)), 2, function(x)
    format(as.hexmode(x), 2)), 2, paste, collapse = ''))
  res <- unlist(unname(Map(paste0, res, as.character(as.hexmode(trans)))))
  res[is.na(color)] <- NA
  return(res)
}

## define cosmology

## planck 2015
cosmo <- list( omega_m =  0.3111, omega_l = 0.6889, omega_b = 0.04897468, N_eff = 2.99, h = 0.6766, w0 = -1., wa = 0., gamma = 0.55, gamma0 = 0.55, gamma1 = 0., eta0 = 0., eta1 = 0., A_xi = 0. )

## define speedoflight in km per sec
c_kmpsec <- speedoflight * 1e-3

E_z <- function( z, c=cosmo ){
    ## dimensionless Hubble rate
    zplus1 <- z + 1.
    aa <- 1. / zplus1  
    omega_k <- 1. - c$omega_l - c$omega_m
    omega_DE <- c$omega_l * exp( 3. * c$wa * ( aa - 1. ) ) / ( aa^( 3. * ( 1. + c$w0 + c$wa ) ) )
    Ez <- sqrt( c$omega_m * zplus1^3. + omega_k * zplus1^2. + omega_DE )
    return( Ez )
}

horizon_size <- function( c=cosmo, zmax=60 ){
    zeq <- 3265. # from WMAP 9; but pretty insensitive to this value
    omega_m <- c$omega_m
    omega_l <- c$omega_l
    h <- c$h
    omega_rad <- omega_m / ( 1 + zeq )
    aa <- seq(0,1.,length.out=10000)
    z <- 1/aa - 1
    integ <- 1 / sqrt( omega_m*aa + omega_l*aa^4 + omega_rad)
    rh <- c_kmpsec / ( 100 * h ) * cumtrapz(aa, integ)
    k_horizon <- 2 * pi / rh
    z_idx <- which( z <= zmax )
    z <- z[z_idx]
    k_horizon <- k_horizon[z_idx]
    result <- list( z=z, k_horizon=k_horizon)
    return( result )
}

evolution_splines <- function( c=cosmo, zmax=60., nsamples=500 ){
    ## interpolate for functions of redshift:
    ## H(z) -- Hubble rate in km/s/Mpc
    ## r(z) -- comoving distance in Mpc
    ## D(z) -- linear growth factor
    ## f(z) -- linear growth rate

    ## redshifts
    z_samp <- seq( 0, zmax, length.out=nsamples )
    zplus1 <- 1. + z_samp
    aa <- 1. / zplus1

    ## set up parameters
    H0 <- 100. * c$h
    omega_k <- 1. - c$omega_m - c$omega_l
    omega_DE <- c$omega_l * exp( 3. * c$wa * ( aa - 1. ) ) / ( aa^( 3. * ( 1. + c$w0 + c$wa ) ) )

    ## sample hubble rate and comoving distance
    Ez <- sqrt( omega_m * zplus1^3. + omega_k * zplus1^2. + omega_DE )
    Hz <- H0 * Ez

    ##########################
    ## CHECK FROM HERE

    ## comoving distances
    r_c <- cumtrapz( z_samp, 1./Ez ) 
    if ( omega_k > 0 ){
        rz <- c_kmpsec / ( H0 * sqrt( omega_k ) ) * sinh( r_c * sqrt( omega_k ) )
    } else if ( omega_k < 0 ){
        rz <- c_kmpsec / ( H0 * sqrt( -omega_k ) ) * sin( r_c * sqrt( -omega_k ) )
    } else {

        rz <- c_kmpsec / H0 * r_c 
    }
    
    ## integrate the linear growth rate to find the linear growth factor
    ## generalised form for the growth rate
    Oma = c$omega_m * zplus1^3. / Ez^2. 
    growth_fac <- Oma^c$gamma0
    Dz <- cumtrapz(log(aa), growth_fac)
    Dz <- exp(Dz)
    
    r <- approxfun( z_samp, rz, method='linear' )    
    H <- approxfun( z_samp, Hz, method='linear' )
    D <- approxfun( z_samp, Dz, method='linear' )
    f <- approxfun( z_samp,  growth_fac, method='linear' )
    
    return( list( r=r, H=H, D=D, f=f) )

}

## get evolution splines

tmp <- evolution_splines( c=cosmo, zmax=60, nsamples = 1000 )
r <- tmp$r
H <- tmp$H
D <- tmp$D
f <- tmp$f

## set up redshift range 
z <- seq( 0, 60., length.out=1000 )
rz <- r(z)
Hz <- H(z)
wl <- ( 1. + z ) * speedoflight / 1420e6  # metres

## min and max baselines for whole array
Dmin <- 43.01293
Dmax <- 1987e3

## min and max k_perp
kperp_min <- 2 * pi * Dmin / ( rz * wl )
kperp_max <- 2 * pi * Dmax / ( rz * wl )

## min and max baselines for the core
Dmax_core <- 3704.49
kperp_min_core <- 2 * pi * Dmin / ( rz * wl )
kperp_max_core <- 2 * pi * Dmax_core / ( rz * wl )

## min and max baselines for the Dutch array
Dmin_dutch <- 2278.98
Dmax_dutch <- 120415.9
kperp_min_dutch <- 2 * pi * Dmin_dutch / ( rz * wl ) 
kperp_max_dutch <- 2 * pi * Dmax_dutch / ( rz * wl )

## for the international array
Dmin_intl <- 277852.1
kperp_min_intl <- 2 * pi * Dmin_intl / ( rz * wl )
kperp_max_intl <- 2 * pi * Dmax / ( rz * wl )

## get the horizon size
res <- horizon_size(c=cosmo, zmax=60)

mycols <- viridis( 6 )

cairo_pdf('k_perp.pdf')
## initial plot just to put something up there
lplot( res$k_horizon, res$z, type='l', col='gray48', lwd=1.8, xlim=c(1e-5,1e5), ylim=c(3,55), log='x', x_lab='k_perp', y_lab='Redshift', xaxs = "i", yaxs = "i", xaxon=FALSE) 
#plot( res$k_horizon, res$z, type='l', col='gray48', lwd=1.8, xlim=c(1e-5,1e5), ylim=c(3,60), log='xy', x_lab='k_perp', y_lab='Redshift', xaxs = "i", yaxs = "i", xaxon=FALSE) 
## plot core
polygon( c(kperp_min_core, rev(kperp_max_core)), c(z,rev(z)), col=tcol(mycols[1]), border=mycols[1] )
## plot dutch
polygon( c(kperp_min_dutch, rev(kperp_max_dutch)), c(z,rev(z)), col=tcol(mycols[3]), border=mycols[3] )
## plot intl
polygon( c(kperp_min_intl, rev(kperp_max_intl)), c(z,rev(z)), col=tcol(mycols[5]), border=mycols[5] )

## plot HBA and LBA frequency ranges in terms of redshift for 21 cm line
hba_zmin <- 1420e6 / 168e6 - 1
hba_zmax <- 1420e6 / 120e6 - 1
lba_zmin <- 1420e6 / 78e6 - 1
lba_zmax <- 1420e6 / 30e6 - 1
## shade outside these lines
polygon( c(1e-5,1e-5,1e5,1e5), c(3,hba_zmin,hba_zmin,3), col=tcol('gray48'), border='gray48', lty=2 )
polygon( c(1e-5,1e-5,1e5,1e5), c(hba_zmax,lba_zmin,lba_zmin,hba_zmax), col=tcol('gray48'), border='gray48', lty=2 )
polygon( c(1e-5,1e-5,1e5,1e5), c(lba_zmax,60,60,lba_zmax), col=tcol('gray48'), border='gray48', lty=2 )

## label HBA and LBA
text( 3e3, 32, 'LBA', font=2)
text( 3e3, 9, 'HBA', font=2)

## shade the horizon
polygon( c(1e-5,1e-5,res$k_horizon), c(3,60,res$z), col=tcol('gray48'))
text( 1e-4, 30, 'Horizon', font=2)
lines( res$k_horizon, res$z, lwd=3)

## legend and axes
legend( 'topright', c('Core','Dutch','International'), pch=22, pt.bg=tcol(mycols[c(1,3,5)]), col=mycols[c(1,3,5)], pt.cex=2 )
axis( 1, at=10^seq(-5,5,2), labels=c(expression(10^-5), expression(10^-3), expression(10^-1), expression(10^1), expression(10^3), expression(10^5)), line=0, cex=1.25 )
mtext( bquote(italic('k')~~"[Mpc"^"-1"*"]"), side=1, line=3, font=3 )
par(xpd=TRUE)
lines( c(3.4e-1,4.3e-1), c(-4,-4) )
lines( c(3.8e-1,3.8e-1), c(-4,-3.3) )
box(which='plot', lwd=2.5)
dev.off()

## now for k-parallel
## ref fig 2 of phil's paper
sb_width <- 195312.5e-6
delta_nu_MHz <- sb_width / 256
total_bw <-  sb_width * 244

k_par_min <- 2 * pi * 1420.406 *  H( res$z ) / ( delta_nu_MHz * c_kmpsec  * ( 1. + res$z )^2. )

k_par_max <- 2 * pi * 1420.406 *  H( res$z ) / ( total_bw * c_kmpsec  * ( 1. + res$z )^2. )

cairo_pdf('k_par.pdf')
lplot( k_par_max, res$z, type='l', col='gray48', lwd=2.5, xlim=c(1e-5,1e5), ylim=c(3,55), log='x', x_lab='k_perp', y_lab='Redshift', xaxs = "i", yaxs = "i", xaxon=FALSE) 
logsamp <- seq( 0, log10(15616), length.out=200)
cols <- viridis(15616, option='C')
z30_val <- c()
for ( x in logsamp ){
    mycol <- cols[10.^x]
    k_par <-  2 * pi * 1420.406 *  H( res$z ) / ( total_bw / 10.^x * c_kmpsec  * ( 1. + res$z )^2. )
    z30_val <- c( z30_val, k_par[which(abs(z-30)==min(abs(z-30)))[1]])
    lines( k_par, res$z, col=tcol(mycol), lwd=1.5 )   
}
#for ( x in 1:64 ){
#    k_par <-  2 * pi * 1420.406 *  H( res$z ) / ( sb_width / x * c_kmpsec  * ( 1. + res$z )^2. )
#    lines( k_par, res$z )   
#}
lines(  k_par_max, res$z, col='gray48', lwd=3 )
k_par_min <- 2 * pi * 1420.406 *  H( res$z ) / ( total_bw / 15616 * c_kmpsec  * ( 1. + res$z )^2. )
lines(  k_par_min, res$z, col='gray48', lwd=3 )
text(c(3e-4,8e2), c(45,45), labels=c(expression(Delta*nu*"=48 MHz"), expression(Delta*nu*"=0.76 kHz")), col='gray48', cex=1.25 )

## shade outside these lines
polygon( c(1e-5,1e-5,1e5,1e5), c(3,hba_zmin,hba_zmin,3), col=tcol('gray48'), border='gray48', lty=2 )
polygon( c(1e-5,1e-5,1e5,1e5), c(hba_zmax,lba_zmin,lba_zmin,hba_zmax), col=tcol('gray48'), border='gray48', lty=2 )
polygon( c(1e-5,1e-5,1e5,1e5), c(lba_zmax,60,60,lba_zmax), col=tcol('gray48'), border='gray48', lty=2 )

## label HBA and LBA
text( 3e3, 32, 'LBA', font=2)
text( 3e3, 9, 'HBA', font=2)


box( which='plot', lwd=2 )
axis( 1, at=10^seq(-5,5,2), labels=c(expression(10^-5), expression(10^-3), expression(10^-1), expression(10^1), expression(10^3), expression(10^5)), line=0, cex=1.25 )
mtext( bquote(italic('k')~~"[Mpc"^"-1"*"]"), side=1, line=3, font=3 )
par(xpd=TRUE)
lines( c(3.3e-1,3.3e-1), c(-4,-3.3) )
lines( c(3.6e-1,3.6e-1), c(-4,-3.3) )
dev.off()


zlims <- c( hba_zmin, hba_zmax, lba_zmin, lba_zmax)
wlz <- ( 1. + zlims ) * speedoflight / 1420e6  # metres
kperp_zmin_intl <- 2 * pi * Dmin_intl / ( r(zlims) * wlz ) 
kperp_zmax_intl <- 2 * pi * Dmax / ( r(zlims) * wlz ) 

cat( 'HBA limits are: $\\sim 10^{', log10(mean(kperp_zmin_intl[c(1,2)])), '}$ - $10^{', log10(mean(kperp_zmax_intl[c(1,2)])), '}$\n' )
cat( 'LBA limits are: $\\sim 10^{', log10(mean(kperp_zmin_intl[c(3,4)])), '}$ - $10^{', log10(mean(kperp_zmax_intl[c(3,4)])), '}$\n' )

k_par_zmin <- 2 * pi * 1420.406 *  H( zlims ) / ( delta_nu_MHz * c_kmpsec  * ( 1. + zlims )^2. )
k_par_zmax <- 2 * pi * 1420.406 *  H( zlims ) / ( total_bw * c_kmpsec  * ( 1. + zlims )^2. )

cat( 'HBA limits are: $\\sim 10^{', log10(mean(k_par_zmax[c(1,2)])), '}$ - $10^{', log10(mean(k_par_zmin[c(1,2)])), '}$\n' )
cat( 'LBA limits are: $\\sim 10^{', log10(mean(k_par_zmax[c(3,4)])), '}$ - $10^{', log10(mean(k_par_zmin[c(3,4)])), '}$\n' )

## plot a wedge



## evaluate at these redshifts
z_values <- c( 30, 35, 40, 45 )
zcols <- viridis(6,option='C')

zz_kperp_min <- c()
zz_kperp_max <- c()
zz_kpar_min <- c()
zz_kpar_max <- c()

cairo_pdf('wedge_limits.pdf')

lplot( 1, 1, xlim=c(3e0,1e3), ylim=c(2e-1,4e2), type='n', log='xy', x_lab='k_perp', y_lab='Redshift', xaxs = "i", yaxs = "i", xaxon=FALSE, yaxon=FALSE ) 
for ( zval in z_values ){
    
    ## get k-perp-min and max
    z_idx <- which( abs(z - zval) == min(abs(z - zval)) )[1]
    zcol <- zcols[which(zval == z_values)]

    perp_min <- kperp_min_intl[z_idx]
    perp_max <- kperp_max_intl[z_idx]
    ## and k-par-min and max
    par_min <- k_par_max[z_idx]
    par_max <- k_par_min[z_idx]
    
    ## theoretical wedge
    fovrad <- 1.25
    fov <- pi * fovrad^2
    fovsr <- pi * fovrad^2 * ( pi / 180. )^2
    scalefac <- sin(fovsr) * H_0 * rz * E_z(z) / ( c_kmpsec * ( 1 + z ) )
    w_kperp_min <- 1e-3 / scalefac[ which(abs(z-zval) == min(abs(z-zval)))[1]]
    w_kperp <- 10^(c(log10(w_kperp_min),5))
    w_kpar <- w_kperp * scalefac[ which(abs(z-zval) == min(abs(z-zval)))[1]]
    
    ## meets the wedge at:
    meets_wedge <- perp_max * scalefac[ which(abs(z-zval) == min(abs(z-zval)))[1]]
    cat( 'z=', zval, ' meets the wedge at k_par=', meets_wedge, '\n')
    ## find the bandwidth
    tmp <- which( abs(z30_val - meets_wedge) == min(abs(z30_val- meets_wedge)))[1]
    tmp_mhz <- total_bw / 10^logsamp[tmp]
    nchan <- ceiling(sb_width/tmp_mhz)
    cat('corresponding to ', tmp_mhz, 'MHz which is a minimum of ', nchan, ' channels\n' )
    
    zz_kperp_min <- c( zz_kperp_min, perp_min )
    zz_kperp_max <- c( zz_kperp_max, perp_max)
    
    zz_kpar_min <- c( zz_kpar_min, meets_wedge )
    zz_kpar_max <- c( zz_kpar_max, par_max )
    
    
    ## leakage
    fovrad <- 2.3
    fovsr <- pi * fovrad^2 * ( pi / 180. )^2
    scalefac <- sin(fovsr) * H_0 * rz * E_z(z) / ( c_kmpsec * ( 1 + z ) )
    wl_kperp_min <- 1e-3 / scalefac[ which(abs(z-zval) == min(abs(z-zval)))[1]]
    wl_kperp <- 10^(c(log10(wl_kperp_min),5))
    wl_kpar <- wl_kperp  * scalefac[ which(abs(z-zval) == min(abs(z-zval)))[1]]
    
    ## meets the wedge at:
    meets_wedge <- perp_max * scalefac[ which(abs(z-zval) == min(abs(z-zval)))[1]]
    cat( 'z=', zval, ' meets the leakage wedge at k_par=', meets_wedge, '\n')
    ## find the bandwidth
    tmp <- which( abs(z30_val - meets_wedge) == min(abs(z30_val- meets_wedge)))[1]
    tmp_mhz <- total_bw / 10^logsamp[tmp]
    nchan <- ceiling(sb_width/tmp_mhz)
    cat('corresponding to ', tmp_mhz, 'MHz which is a minimum of ', nchan, ' channels\n' )
    
    
    
    ## shade the allowable region
    polygon( c(perp_min,perp_min,perp_max,perp_max), c(par_min,par_max,par_max,par_min), col=tcol(zcol))
    
    ## cut off stuff outside plot
    polygon( c(w_kperp,1e5), c(w_kpar,1e-3), col=tcol('gray48') )
    lines( w_kperp, w_kpar, col=zcol, lwd=3 )
    lines( wl_kperp, wl_kpar, col=zcol, lwd=3, lty=3 )
    
}

text( 6e0, 0.125e2, 'z=30', col=tcol(zcols[1],150), font=2, cex=2)
text( 6e0, 0.25e2, 'z=35', col=tcol(zcols[2],150), font=2, cex=2)
text( 6e0, 0.5e2, 'z=40', col=tcol(zcols[3],150), font=2, cex=2)
text( 6e0, 1e2, 'z=45', col=tcol(zcols[4],150), font=2, cex=2)

text( 2e2, 8e-1, 'Foreground', font=2, cex=1.5 )
text( 2e2, 5.8e-1, 'wedge', font=2, cex=1.5 )

box( which='plot', lwd=2 )
axis( 1, at=10^seq(1,3,1), labels=c( expression(10^1), expression(10^2), expression(10^3)), line=0, cex=1.25 )
axis( 2, at=10^seq(0,2,1), labels=c( expression(10^0), expression(10^1), expression(10^2)), line=0, cex=1.25 )
mtext( bquote(italic('k')~~"[Mpc"^"-1"*"]"), side=1, line=3, font=3 )
par(xpd=TRUE)
lines( c(6.4e0,7.2e0), c(2.9e-2,2.9e-2) )
lines( c(6.75,6.75), c(2.9e-2,3.3e-2) )
mtext( bquote(italic('k')~~"[Mpc"^"-1"*"]"), side=2, line=3, font=3 )
lines( c(3.1e-2,3.6e-2), c(6.3,6.3) )
lines( c(3.1e-2,3.6e-2), c(6.55,6.55) )

dev.off()

## calculate sensitivities

## write a table
tfile <- 'table.tex'
cat( '\\begin{table}[!h]\n', file=tfile )
cat( '\\caption{$T_b$ sensitivity values for LOFAR}\n', file=tfile, append=TRUE )
cat( '\\label{t1}\n\\begin{tabular}{cccccccc}\n\\hline\n', file=tfile, append=TRUE )
cat( '$z$ & \\kperp\\ range & $\\theta$ range & \\kpar\\ range & $\\Delta\\nu$ range & $T_{b,\\textrm{min}$ & $T_{b,\\textrm{max}$ & $T_{b,\\textrm{median}$ \\\\ \n', file=tfile, append=TRUE )
cat( '& [\\permpc] & [arcsec] & [\\permpc] & [kHz] & [K] & [K] & [K] \\\\ \n \\hline \\\\[-8pt] \n', file=tfile, append=TRUE )

## reverse engineer the k-par and k-perp limits to resolution and frequency coverage for several redshift values
for ( i in seq(length(z_values)) ){

    ## sample evenly in log space
    kperp_range <- 10^seq(log10(zz_kperp_min[i]),log10(zz_kperp_max[i]),0.01)
    kpar_range <- 10^seq(log10(zz_kpar_min[i]),log10(zz_kpar_max[i]),0.01)
    
    resolution <- 2 * pi * 1.22 / kperp_range / rz[which(abs(z-z_values[i]) == min(abs(z-z_values[i])))[1]]
    resolution_asec <- resolution * 206265
    
    delta_nu <- 2 * pi * 1420e6 * Hz[which(abs(z-z_values[i]) == min(abs(z-z_values[i])))[1]] / c_kmpsec / ( 1 + z_values[i])^2. / kpar_range
    
    obs_freq <- 1420e6 / ( 1 + z_values[i])
    obs_time <- 1 ## seconds

    sens_file <- paste( 'sensitivity_matrix_linear_z',as.character(z_values[i]),'.dat', sep='' )
    if ( !file.exists(sens_file) ){
        sensitivity_matrix <- matrix( data=NA, nrow=length(resolution), ncol=length(delta_nu) )
        pb <- txtProgressBar( min=1, max=length(resolution), style=3 )
        for ( ii in seq(1,length(resolution)) ){
            for ( jj in seq(1,length(delta_nu)) ){
                set_resolution <- resolution_asec[ii]
                set_resolution_tolerance <- 0.1 * set_resolution
                stations_to_use <- find_n_stations( set_resolution, obs_freq, bl_info, res_tol=set_resolution_tolerance )
                sensitivity_matrix[ii,jj] <- getArraySens( obs_freq, 1., delta_nu[jj], obs_time, stations_to_use$ncore, stations_to_use$nremote, stations_to_use$nint, 30, FALSE, quiet=TRUE )
            }
            setTxtProgressBar( pb, ii )
        }
        close(pb)
        
        ## write it out for future use
        ## values are in JANSKY
        write.table( sensitivity_matrix, file=sens_file )
        
    } else {
        ## read in the data
        sensitivity_matrix <- as.matrix( read.table( sens_file ), as.is=TRUE )
        inf_index <- which( !is.finite( sensitivity_matrix ), arr.ind=TRUE )
        sensitivity_matrix[inf_index] <- NA
    }
    
    ## convert to temperature brightness
    tb_matrix <- 2.445395e24 * sensitivity_matrix  / ( obs_freq^2. )
    for ( ii in seq(1,length(resolution) ) ){
            tb_matrix[ii,] <- tb_matrix[ii,] / ( resolution_asec[ii]^2. )
    }
    
    time_scales <- c( 1, 800*3600, 8000*3600, 87600*3600 ) ## in seconds
    max_val <- max( tb_matrix, na.rm=TRUE )
    min_val <- min( tb_matrix / sqrt(max(time_scales)) )
    ## log min and max
    l_min <- log10( min_val )
    l_max <- log10( max_val )
    ## set the logarithmic levels 
    nlevels <- 29
    log_levels <- seq( l_min, l_max, length.out=nlevels )
    log_levels <- c( log_levels, max( log_levels ) + log_levels[2] - log_levels[1] )
    mylevels <- 10.^log_levels
    mycols <- viridis( nlevels )
    tmp_cb_labels = format( log10(mylevels), digits=1, nsmall=1 )
    cb_labels = rep('',nlevels)
    cb_labels[seq(1,29,3)] = tmp_cb_labels[seq(1,29,3)]
        
    for ( t in time_scales ){
        
        if ( t-1 == 0 ) tmp_tb <- tb_matrix else tmp_tb <- tb_matrix / sqrt( t - 1 )
        tb_breaks <- which( mylevels >= min( tmp_tb, na.rm=TRUE ) & mylevels <= max( tmp_tb, na.rm=TRUE ) )
        ## expand one each side
        tb_breaks <- seq( min(tb_breaks)-1, max(tb_breaks)+1 )
        
        if ( t/3600 < 1 ) tmp <- paste( as.character(t), 's', sep='' ) else tmp <- paste( as.character(t/3600), 'hrs', sep='' )

        outfile <- paste( 'brightness_temp_z', as.character(z_values[i]), '_t', tmp, '.pdf', sep='' )
        cairo_pdf( outfile)
        par( mar=c(5,5,5,7) )
        #image( tmp_tb, breaks=mylevels[tb_breaks], col=mycols[tb_breaks[1:(length(tb_breaks)-1)]], axes=FALSE )
        filled.contour( kperp_range, kpar_range, tmp_tb, key.title='Brightness Temperature [K]', plot.axes = { axis(1); axis(2); contour( kperp_range, kpar_range, tmp_tb, add=TRUE ) } )#, breaks=mylevels[tb_breaks], col=mycols[tb_breaks[1:(length(tb_breaks)-1)]], axes=FALSE )
        ## k-perp axis
        mtext( bquote(italic('k')~~"[Mpc"^"-1"*"]"), side=1, line=3, font=3 )
        par(xpd=TRUE)
	xrange = max(kperp_range)-min(kperp_range)
	yrange = max(kpar_range)-min(kpar_range)
	a = xrange*(0.45)+min(kperp_range)
	b = xrange*(0.462)+min(kperp_range)
	c = yrange*(-0.15) + min(kpar_range)
	d = yrange*(-0.135) + min(kpar_range)
	lines( c(a,b), c(c,c) )
	lines( c(b-a,b-a)*0.5+a, c(c,d) )
        mtext( bquote(italic('k')~~"[Mpc"^"-1"*"]"), side=2, line=3, font=3 )
	a = xrange*(-0.185) + min(kperp_range)
	b = xrange*(-0.17) + min(kperp_range)
	c = yrange*0.449 + min(kpar_range)
	d = yrange*0.449 + min(kpar_range)
        lines( c(a,b), c(c,d) )
        lines( c(a,b), c(c+yrange*0.005,c+yrange*0.005) )
        mtext( bquote(italic('z=')*.(z_values[i])), side=3, line=1, font=2, cex=1.5)
        dev.off()
        
        if ( t == max( time_scales) ){
        
            ## kperp range
            tmp1 <- paste( '$10^{', format(log10(zz_kperp_min[i]), digits=1, nsmall=1), '}$ - $10^{', format( log10(zz_kperp_max[i]), digits=1, nsmall=1 ), '}$', sep='' )
            ## res range 
            rrange <- paste( format(resolution_asec[1], digits=1, nsmall=1), ' - ', format(resolution_asec[length(resolution_asec)], digits=1, nsmall=1), sep='' )
            ## kpar range
            tmp2 <- paste( '$10^{', format(log10(zz_kpar_min[i]), digits=1, nsmall=1), '}$ - $10^{', format( log10(zz_kpar_max[i]), digits=1, nsmall=1 ), '}$', sep='' )
            ## delta_nu range
            frange <- paste( format( delta_nu[1]*1e-3, digits=1, nsmall=1 ), ' - ', format(delta_nu[length(delta_nu)]*1e-3, digits=1, nsmall=1), sep='' )
            ## t_b min
            tmin <- paste( '$10^{', format( log10(min(tmp_tb, na.rm=TRUE)), digits=1, nsmall=1), '}$', sep='' )
            ## t_b max
            tmax <- paste( '$10^{', format( log10(max(tmp_tb, na.rm=TRUE)), digits=1, nsmall=1), '}$', sep='' )
            ## t_b median
            tmedian <- paste( '$10^{', format( log10(median(tmp_tb, na.rm=TRUE)), digits=1, nsmall=1), '}$', sep='' )
            
            ss <- paste( as.character( z_values[i]), tmp1, rrange, tmp2, frange, tmin, tmax, tmedian, sep=' & ')
            ss <- paste( ss, ' \\\\ \n' , sep='')
            cat( ss, file=tfile, append=TRUE )
            
        }
    
    }

}

cat( '\\hline \n\\end{tabular}\n\\vspace*{-4pt}\n\\end{table}', file=tfile, append=TRUE )

