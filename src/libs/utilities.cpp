#include "../../include/utilities.h"
#include "../../include/Logger.h"

#ifndef NA
#define NA -99.0
#endif

adstring stripExtension(adstring fileName){
  /*
		This function strips the file extension
		from the fileName argument and returns
		the file name without the extension.
  */
  const int length = fileName.size();
  for (int i=length; i>=0; --i){
    if (fileName(i)=='.'){
      return fileName(1,i-1);
    }
  }
  return(fileName);
}

ivector getIndex(const dvector& a, const dvector& b){
  /*
    Returns a dvector which contains the elements of dvector a where
    the elements of b are not NA
    Assumes the length of a and b are the same.
   */
  int i,j,n;
  n = 0;
  j = 1;
  for(i = a.indexmin(); i <= a.indexmax(); i++){
    if(b(i) != NA){
      n++;
    }
  }
  ivector tmp(1,n);
  for(i = a.indexmin(); i <= a.indexmax(); i++){
    if(b(i) != NA){
      tmp(j++) = a(i);
    }
  }
  return(tmp);
}

void write_proj_headers(ofstream &ofsP,
                        int syr,
                        int nyr,
                        int nfleet,
                        int n_ags,
                        int ngroup,
                        bool include_msy,
                        bool include_sbo){
  int i, ig;
  // Write the decision table headers for projection years
  ofsP<<"TAC"<<",";
  ofsP<<"B"<<nyr + 1<<",";
  ofsP<<"B"<<nyr + 2<<",";
  ofsP<<"R"<<nyr<<",";
  ofsP<<"R"<<nyr + 1<<",";
  if(include_sbo){
    ofsP<<"B0"<<",";
    ofsP<<"04B0"<<",";
    ofsP<<"03B0"<<",";
    ofsP<<"02B0"<<",";
  }
  ofsP<<"B"<<syr<<",";
  ofsP<<"B"<<nyr + 2<<"_B"<<nyr + 1<<",";
  if(include_sbo){
    ofsP<<"B"<<nyr + 2<<"_04B0"<<",";
    ofsP<<"B"<<nyr + 2<<"_02B0"<<",";
  }
  ofsP<<"B"<<nyr + 2<<"_B"<<syr<<",";
  for(ig = 1; ig <= n_ags; ig++){
    for(i = 1; i <= nfleet; i++){
      ofsP<<"fleet_"<<i<<"_sex_"<<ig<<"_F"<<nyr<<",";
      ofsP<<"fleet_"<<i<<"_sex_"<<ig<<"_F"<<nyr + 1<<",";
      ofsP<<"fleet_"<<i<<"_sex_"<<ig<<"_U"<<nyr<<",";
      ofsP<<"fleet_"<<i<<"_sex_"<<ig<<"_U"<<nyr + 1<<",";
    }
  }
  ofsP<<"UT"<<",";
  for(ig = 1; ig <= n_ags; ig++){
    // u20
    ofsP<<"U20_sex_"<<ig;
    if(ig < n_ags){
	ofsP<<",";
    }
  }
  if(include_msy){
    //MSY based ref points
    ofsP<<","<<"BMSY"<<",";
    ofsP<<"B"<<nyr + 2<<"_BMSY"<<",";
    ofsP<<"B"<<nyr + 2<<"_08BMSY"<<",";
    ofsP<<"B"<<nyr + 2<<"_04BMSY"<<",";
    for(i = 1; i <= nfleet; i++){
      ofsP<<"fleet_"<<i<<"_FMSY"<<",";
      ofsP<<"fleet_"<<i<<"_UMSY"<<",";
      for(ig = 1; ig <= n_ags; ig++){
        ofsP<<"fleet_"<<i<<"_sex_"<<ig<<"_F"<<nyr + 1<<"_FMSY"<<",";
        ofsP<<"fleet_"<<i<<"_sex_"<<ig<<"_U"<<nyr + 1<<"_UMSY";
        if(!(ig == n_ags && i == nfleet)){
          ofsP<<",";
        }
      }
    }
  }
  ofsP<<"\n";
}

void write_proj_output(ofstream &ofsP,
                       int syr,
                       int nyr,
                       int nage,
                       int nfleet,
                       int n_ags,
                       int ngroup,
                       double tac,
                       int pyr,
                       dvector p_sbt,
                       dvector p_rt,
                       d3_array p_ft,
                       d3_array p_N,
                       dvar3_array M,
                       dmatrix ma,
                       dmatrix dWt_bar,
                       dvar3_array ft,
                       double sbo,
                       dmatrix fmsy,
                       dvector bmsy,
                       bool include_msy,
                       bool include_sbo){

  int i, ig;
  double ut = tac / (tac + p_sbt(pyr));

  // Write the projection output to the file
  ofsP<<tac<<",";
  ofsP<<p_sbt(nyr + 1)<<",";
  ofsP<<p_sbt(nyr + 2)<<",";
  ofsP<<p_rt(nyr)<<",";
  ofsP<<p_rt(nyr + 1)<<",";
  if(include_sbo){
    ofsP<<sbo<<",";
    ofsP<<0.4 * sbo<<",";
    ofsP<<0.3 * sbo <<",";
    ofsP<<0.2 * sbo <<",";
  }
  ofsP<<p_sbt(syr)<<",";
  ofsP<<p_sbt(nyr + 2) / p_sbt(nyr + 1)<<",";
  if(include_sbo){
    ofsP<<p_sbt(nyr + 2) / (0.4 * sbo)<<",";
    ofsP<<p_sbt(nyr + 2) / (0.2 * sbo)<<",";
  }
  ofsP<<p_sbt(nyr + 2) / p_sbt(syr);
  for(ig = 1; ig <= n_ags; ig++){
    for(i = 1; i <= nfleet; i++){
      ofsP<<",";
      ofsP<<p_ft(ig,nyr,i)<<",";
      ofsP<<p_ft(ig,nyr + 1,i)<<",";
      ofsP<<1.0 - mfexp(-p_ft(ig,nyr,i))<<",";
      ofsP<<1.0 - mfexp(-p_ft(ig,nyr + 1,i))<<",";
    }
  }
  ofsP<<ut<<",";
  for(ig = 1; ig <= n_ags; ig++){
    // u20
    ofsP<<tac / ((p_N(ig,pyr)(3,nage) * exp(-value(M(1)(nyr,3)))) * dWt_bar(1)(3,nage))<<",";
  }
  if(include_msy){
    //MSY based ref points
    ofsP<<","<<bmsy<<",";
    ofsP<<p_sbt(nyr + 2) / bmsy<<",";
    ofsP<<p_sbt(nyr + 2) / (0.8 * bmsy)<<",";
    ofsP<<p_sbt(nyr + 2) / (0.4 * bmsy);
    // The following loop assumed only one group,
    // Which is the '1' in the fmsy(1,i) indexing
    for(i = 1; i <= nfleet; i++){
      ofsP<<",";
      ofsP<<fmsy(1, i)<<",";
      ofsP<<1. - mfexp(-fmsy(1, i));
      for(ig = 1; ig <= n_ags; ig++){
        ofsP<<",";
        ofsP<<p_ft(ig, nyr + 1,i) / fmsy(1,i)<<",";
        ofsP<<(1.0 - mfexp(-p_ft(ig, nyr + 1,i))) /
              (1.0 - mfexp(-fmsy(1, i)));
      }
    }
  }
  ofsP<<"\n";
}
