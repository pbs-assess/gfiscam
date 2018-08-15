 ## ------------------------------------------------------------------------- ##
## CONTROLS FOR LEADING PARAMETERS                                           ##
##  Prior descriptions:                                                      ##
##                      -0 uniform      (0,0)                                ##
##                      -1 normal       (p1=mu,p2=sig)                       ##
##                      -2 lognormal    (p1=log(mu),p2=sig)                  ##
##                      -3 beta         (p1=alpha,p2=beta)                   ##
##                      -4 gamma        (p1=alpha,p2=beta)                   ##
## ------------------------------------------------------------------------- ##
## npar
7
##____________________________________________________________________________ ##
##  ival        lb      ub      phz     prior    p1      p2      parameter name
## ____________________________________________________________________________ ##
    7.28       -5.0    15       4        0       -5.0    15        #log_ro 
    0.80        0.2     1.0     4        3       10.0    4.925373  #steepness,h
   -0.69186    -5.0    5.0      3        1       -0.7985077 0.4    #log.m
    7.09       -5.0    15       1        0       -5.0    15        #log_avgrec #log_rbar
    5.97       -5.0    15       1        0       -5.0    15        #log_recinit
    0.413297   0.001   0.999    3        3       17.08696    39.0559    #rho 
    1.22062    0.01    5.0      4        4       25.0  28.75   #kappa (precision)
## ____________________________________________________________________________ ##
## ------------------------------------------------------------------------- ##
##
## ------------------------------------------------------------------------- ##
## CONTROL PARAMETERS FOR AGE/SIZE COMPOSITION DATA FOR na_gears             ##
## ------------------------------------------------------------------------- ##
## Likelihood type for each gear:
##     -1 : multivariate logistic (dmvlogistic)
##     -2 : multinomial, sample size based on input data
##     -3 : logistic_normal, no autocorrelation, AR1, AR2.
##     -4 : logistic_normal, AR1
##     -5 : logistic_normal, AR2
## ------------------------------------------------------------------------- ##
## Number of columns == na_gears. One column for each gear with age data
   1       2    3    ## Gear Index
   1       1    1    ## Likelihood type
   0.02   0.02  0.02    ## Minimum proportion for aggregation & tail compression
   0.0   0.0  0.0    ## Small constant to add to comps & renormalize
  -1      -1   -1    ## phase for log_age_tau2 estimation.
  -3      -3   -3    ## phase for phi1 estimation: bounded (-1,1) AR1
  -2      -2   -2    ## phase for phi2 estimation: bounded (0,1)  AR2
  -2      -2   -2    ## phase for degrees of freedom for student T.
  -12345             ## int check (-12345), one value only, not one for each gear
## ------------------------------------------------------------------------- ##
## ____________________________________________________________________________ ##
## _________________________SELECTIVITY PARAMETERS_____________________________ ##
## ------------------------------------------------------------------------- ##
## SELECTIVITY PARAMETERS Columns for gear                                   ##
## NB: To mirror another gear, use (-ve) phase with the mirrored gear number.##
## OPTIONS FOR SELECTIVITY (isel_type):                                      ##
##      1) logistic selectivity parameters                                   ##
##      2) selectivity coefficients                                          ##
##      3) a constant cubic spline with age-nodes                            ##
##      4) a time varying cubic spline with age-nodes                        ##
##      5) a time varying bicubic spline with age & year nodes.              ##
##      6) fixed logistic (set isel_type=6, and estimation phase to -1)      ##
##      7) logistic function of body weight.                                 ##
##      8) logistic with weight deviations (3 parameters)  
## WHAT HAPPENED TO 9) AND 10)
##      11) logistic selectivity with 2 parameters based on mean length      ##
##      12) length-based selectivity coefficients with spline interpolation  ##
##      sig=0.05 0.10 0.15 0.20 0.30 0.40 0.50                               ##
##      wt =200. 50.0 22.2 12.5 5.56 3.12 2.00                               ##
## ------------------------------------------------------------------------- ##
1		  1		  1	   6		  6     # 1  -selectivity type ivector(isel_type) for gear
3.0		3.0		4.0	 2.055	2.055 # 2  -Age/length at 50% selectivity (logistic)
0.25	0.25	0.25 0.05   0.05  # 3  -STD at 50% selectivity (logistic)
5		  5		  5	   0		  0  	  # 4  -No. of age nodes for each gear (0=ignore)
12    3		  10	 0		  0     # 5  -No. of year nodes for 2d spline(0=ignore)
2    2     2    -4	   -5     # 6  -Phase of estimation (-1(gear#) for fixed) If neg number, it reflects a mirroring of another gear's selectivity.  
125.0 125.0	12.5  12.5	12.5  # 7  -Penalty wt for 2nd differences w=1/(2*sig^2)
50.0  50.0	200.0 200.0 200.0 # 8  -Penalty wt for dome-shaped w=1/(2*sig^2)
12.5 	12.5 	12.5 12.5 	12.5  # 9  -Penalty wt for time-varying selectivity #PCOD NOT UPDATED
 1    	1    	1    1    	1     #10  -n_sel_blocks (number of selex blocks) #PCOD NOT UPDATED
## ------------------------------------------------------------------------- ##
## Start year of each time block: 1 row for each gear
1951
1972
1972
1951 
1951
##
## ____________________________________________________________________________ ##
##                             Priors for Survey q                              ##
## ____________________________________________________________________________ ##                                                    ##
## Prior type:                                                               ##
##       0 - uninformative prior                                             ##
##       1 - normal prior density for log(q)                                 ##
##       2 - random walk in q                                                ##
## Need one column for each survey.                                          ##
## ------------------------------------------------------------------------- ##
2              # -number of surveys (nits)  
 0     1       # -prior type (see legend above)
 0     0.      # -prior log(mean) mean q =   exp(-0.6931472) =  0.5
 1     0.01    # -prior sd
## ------------------------------------------------------------------------- ##
## CONTROLS FOR FITTING TO MEAN WEIGHT DATA
## ------------------------------------------------------------------------- ##
0     # 1 = fit to annual mean weights, 0 = do not fit to annual mean weights
1	    # Number of annual mean weight series
0.2  # SD for likelihood for fitting to annual mean weight (one for each series)
## ------------------------------------------------------------------------- ##
## ------------------------------------------------------------------------- ##
## OTHER MISCELANEOUS CONTROLS                                               ##
## ------------------------------------------------------------------------- ##
  0          # 1  -verbose ADMB output (0=off, 1=on)
 1          # 2  -recruitment model (1=beverton-holt, 2=ricker)
 0.100      # 3  -std in observed catches in first phase.
 0.0707     # 4  -std in observed catches in last phase.
 0          # 5  -Assume unfished equilibrium in first year (0=FALSE, 1=TRUE)
 1.00       # 6  -Maternal effects multiplier
 0.20       # 7  -Mean fishing mortality for regularizing the estimates of Ft
 0.01       # 8  -std in mean fishing mortality in first phase
 2.00       # 9  -std in mean fishing mortality in last phase
 3          # 10 -phase for estimating m_deviations (use -1 to turn off mdevs)
 0.1        # 11 -std in deviations for natural mortality
12          # 12 -number of estimated nodes for deviations in natural mortality
 1.00          # 13 -fraction of total mortality that takes place prior to spawning
 0          # 14 -number of prospective years to start estimation from syr
 0          # 15 -switch for IFD distribution in selectivity simulations
 0          # 16 -toggle to fit to annual mean weights for commercial catch
 0          # 17 -toggle to perform "slow" fmsy test (runs model out 100 years)
            #      this produces a file called TEST_frp.rep. Note that even
            #      if you set this to 0, the routine will be run if control 13
            #      is greater than 1.
 0.00001    # 18 -precision for F for the "slow" fmsy calculations (only used if
            #      control 17 is 1). This must be a maximum of 0.0001 or the
            #      program will stop.
 2          # 19 -maximum F for the "slow" fmsy calculations (only used if
            #      control 17 is 1). If this is greater than 1, a warning will
            #      be issued that it will take a long time to run.
 1          # 20 -report B0 only in MCMC calculations if the "slow" msy
            #      routine was run. MSY-based reference points will not be
            #      output for MCMCs. This control is only used if control 13
            #      is greater than 0.
## ------------------------------------------------------------------------- ##
## MARKER FOR END OF CONTROL FILE (eofc)
## ------------------------------------------------------------------------- ##
999
