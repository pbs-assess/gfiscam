#include "../../include/ddirmultinom.h"
#include "../../include/Logger.h"

dvariable ddirmultinom(const dvar_vector& obs, const dvar_vector& p, const dvariable& log_phi){
  dvariable phi = exp(log_phi);
  dvariable N = sum(obs);
  dvariable ll = gammln(N + 1.0) +
                 gammln(phi) -
                 gammln(N + phi);
  for(int a = obs.indexmin(); a <= obs.indexmax(); a++){
    ll += -gammln(obs(a) + 1.0) +
           gammln(obs(a) + phi * (p(a) + 1.0e-15)) -
           gammln(phi * (p(a) + 1.0e-15));
  }
  return(ll);
}
