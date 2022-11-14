#include "../../include/ddirmultinom.h"
#include "../../include/Logger.h"

dvariable ddirmultinom(const dvar_vector& obs, const dvar_vector& p, const dvariable& log_phi){
  RETURN_ARRAYS_INCREMENT();
  dvariable phi = exp(log_phi);
  dvariable N = sum(obs);
  dvariable ll = gammln(N + 1.0) +                       // top of first term (eq. 4/10 Thorson et al. 2017)
                 gammln(phi) -                           // top of second term
                 gammln(N + phi);                        // bottom of second term
  for(int a = obs.indexmin(); a <= obs.indexmax(); a++){
    ll += -gammln(N * obs(a) + 1.0) +                    // bottom of first term
           gammln(N * obs(a) + phi * (p(a))) -           // top of third term
           gammln(phi * (p(a)));                         // bottom of third term
  }
  RETURN_ARRAYS_DECREMENT();
  return(ll);
}
