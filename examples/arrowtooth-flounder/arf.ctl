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
## ival      lb      ub     phz prior p1     p2     #parameter
3.000000  -2.00    6.00   1   0     -2. 6.00  #log_ro
0.800000   0.20    1.00   1   3     13.4 2.4  #h
-1.609400  -3.00    2.00   1   1    -1.6094 0.2   #log_m g&b 0.2
# -1.203973  -3.00    2.00   2   1    -1.2040 0.075 #log_m g&b 0.3
# -1.049822  -3.00    2.00  -2   1    -1.0498 0.075 #log_m g&b 0.35
3.000000  -2.00    6.00   1   0     -2.0 6.00  #log_rbar
3.000000  -2.00    6.00   1   0     -2. 6.00  #log_rinit
0.0588235  0.00001 0.99  -3   3     3.0000 5.00  #rho
1.470588   0.01    5.00  -2   4     1.0100 1.01  #vartheta #sig=0.2,tau=0.8
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
   1       2    4     5    ## Gear Index
   1       1    1     1    ## Likelihood type
   0.0   0.0  0.0   0.0   ## Minimum proportion for aggregation & tail compression
   0.0   0.0  0.0   0.0   ## Small constant to add to comps & renormalize
   1       2    2     2    ## phase for log_age_tau2 estimation.
   3       3    3     3    ## phase for phi1 estimation: bounded (-1,1) AR1
  -2      -2   -2    -2   ## phase for phi2 estimation: bounded (0,1)  AR2
  -2      -2   -2    -2    ## phase for degrees of freedom for student T.
  -12345                 ## int check (-12345), one value only, not one for each gear
## ------------------------------------------------------------------------- ##
##
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
##      8) logistic with weight deviations (3 parameters)                    ##
##      11) logistic selectivity with 2 parameters based on mean length      ##
##      12) length-based selectivity coefficients with spline interpolation  ##
##      sig=0.05 0.10 0.15 0.20 0.30 0.40 0.50                               ##
##      wt =200. 50.0 22.2 12.5 5.56 3.12 2.00                               ##
## ------------------------------------------------------------------------- ##
 1    1    6    1    1   # 1  -selectivity type ivector(isel_type) for gear
 2   2    9   2    2   # 2  -Age/length at 50% selectivity (logistic)
 0.5  1.  0.5  0.5  0.5 # 3  -STD at 50% selectivity (logistic)
 7    7    7    7    7   # 4  -No. of age nodes for each gear (0=ignore)
12   12   12   12   12   # 5  -No. of year nodes for 2d spline(0=ignore)
 2    2   -3   2    2   # 6  -Phase of estimation (-1 for fixed) If neg number, it reflects a mirroring of another gear's selectivity.
 2.0  2.0  2.0  2.0  2.0 # 7  -Penalty wt for 2nd differences w=1/(2*sig^2)
 2.0  2.0  2.0  2.0  2.0 # 8  -Penalty wt for dome-shaped w=1/(2*sig^2)
12.5 12.5 12.5 12.5 12.5 # 9  -Penalty wt for time-varying selectivity
 1    1    1    1    1   #10  -n_sel_blocks (number of selex blocks)
## ------------------------------------------------------------------------- ##
## Start year of each time block: 1 row for each gear
1996
1996
1996
1996
1996
##
## ------------------------------------------------------------------------- ##
## PRIORS FOR SURVEY Q                                                       ##
## Prior type:                                                               ##
##       0 - uninformative prior                                             ##
##       1 - normal prior density for log(q)                                 ##
##       2 - random walk in q                                                ##
## Need one column for each survey.                                          ##
## ------------------------------------------------------------------------- ##
4                        # -number of surveys (nits)
#0     0     0     0     # -prior type (see legend above)
 1     1     1     1     # -prior type (see legend above)
 -0.6931472     -0.6931472    -0.6931472     -0.6931472     # -prior log(mean)	mean q =   exp(-0.6931472) =  0.5
#2.0   2.0   2.0   2.0   # -prior sd
 1.0   1.0   1.0   1.0   # -prior sd
## ------------------------------------------------------------------------- ##
##
## CONTROLS FOR FITTING TO MEAN WEIGHT DATA
## ------------------------------------------------------------------------- ##
0     # 1 = fit to annual mean weights, 0 = do not fit to annual mean weights
1	    # Number of annual mean weight series
0.25  # SD for likelihood for fitting to annual mean weight (one for each series)
## ------------------------------------------------------------------------- ##
## ------------------------------------------------------------------------- ##
## OTHER MISCELANEOUS CONTROLS                                               ##
## ------------------------------------------------------------------------- ##
 0          # 1  -verbose ADMB output (0=off, 1=on)
 1          # 2  -recruitment model (1=beverton-holt, 2=ricker)
 0.2        # 3  -std in observed catches in first phase.
 0.05        # 4  -std in observed catches in last phase.
 0          # 5  -Assume unfished equilibrium in first year (0=FALSE, 1=TRUE)
 1.00       # 6  -Maternal effects multiplier
 0.10       # 7  -Mean fishing mortality for regularizing the estimates of Ft
 1.0        # 8  -std in mean fishing mortality in first phase
 2.0        # 9  -std in mean fishing mortality in last phase
-1          # 10 -phase for estimating m_deviations (use -1 to turn off mdevs)
 0.1        # 11 -std in deviations for natural mortality
12          # 12 -number of estimated nodes for deviations in natural mortality
 0.00       # 13 -fraction of total mortality that takes place prior to spawning
 0          # 14 -number of prospective years to start estimation from syr
 0          # 15 -switch for IFD distribution in selectivity simulations
 0          # 16 1= fit to annual mean weights for commercial catch; 0=do not do this
##
## ------------------------------------------------------------------------- ##
## MARKER FOR END OF CONTROL FILE (eofc)
## ------------------------------------------------------------------------- ##
999
