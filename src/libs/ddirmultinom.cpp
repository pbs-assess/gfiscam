#include "../../include/ddirmultinom.h"
#include "../../include/Logger.h"

dvariable ddirmultinom(const dvar_vector& obs,
		       const dvar_vector& p,
		       const dvariable& log_phi,
                       const dvariable& samp_size,
                       bool linear){

 static int ntimes = 1;
  if(linear){
    if(ntimes == 1)
      LOG<<"ddirmultinom: linear DM method"<<"\n";
  }else{
    if(ntimes == 1)
      LOG<<"ddirmultinom: saturating DM method"<<"\n";
  }
  ntimes++;

  // See equation 4 in Thorsen et. al. 2017
  // Model-based estimates of effective sample
  //  size in stock assessment models using the
  //  Dirichlet-multinomial distribution
  // If linear is true, use the linear method of DM
  // If linear is false, use the saturating method of DM
  RETURN_ARRAYS_INCREMENT();
  dvector obs_nums = value(obs) * value(samp_size);
  dvariable N = sum(obs_nums);
  dvariable phi = N * exp(log_phi);
  dvariable ll = gammln(N + 1.0) + // top of first term - Equation 4 or 10
                 gammln(phi) -     // top of second term
                 gammln(N + phi);  // bottom of second term
  for(int a = obs_nums.indexmin(); a <= obs_nums.indexmax(); a++){
    ll += -gammln(obs_nums(a) + 1.0) + // bottom of first term
           gammln(obs_nums(a) + phi * (p(a) + 1.0e-15)) - // top of third term
           gammln(phi * (p(a) + 1.0e-15));  // bottom of third term
  }
  RETURN_ARRAYS_DECREMENT();
  return(ll);
}
