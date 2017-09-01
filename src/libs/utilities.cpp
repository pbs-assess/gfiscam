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
                        bool include_msy,
                        bool include_bo){
  // Write the decision table headers for projection years
  ofsP<<"TAC"                  <<",";
  ofsP<<"B"<<nyr+1             <<",";
  ofsP<<"B"<<nyr+2             <<",";
  if(include_bo){
    ofsP<<"B0"                   <<",";
    ofsP<<"04B0"                 <<",";
    ofsP<<"02B0"                 <<",";
  }
  ofsP<<"B"<<syr               <<",";
  ofsP<<"B"<<nyr+2<<"B"<<nyr+1 <<",";  //want probability B2016<B2015 - this will be < 1 if true
  ofsP<<"B"<<nyr+2<<"04B0"     <<",";  //want probability B2016<0.4B0 - this will be < 1 if true
  ofsP<<"B"<<nyr+2<<"02B0"     <<",";  //want probability B2016<0.4B0 - this will be < 1 if true
  ofsP<<"B"<<nyr+2<<"B"<<syr   <<",";
  ofsP<<"F"<<nyr               <<",";
  ofsP<<"F"<<nyr+1             <<",";
  ofsP<<"F"<<nyr+1<<"F"<<nyr   <<",";  //want probability F2015>F2014     - this will be > 1 if true
  ofsP<<"U"<<nyr+1             <<",";
  ofsP<<"U"<<nyr+1<<"U"<<nyr   <<",";
  if(include_msy){
    //MSY based ref points
    ofsP<<"BMSY"                 <<",";
    ofsP<<"B"<<nyr+2<<"BMSY"     <<",";  //want probability B2016<BMSY - this will be < 1 if true
    ofsP<<"B"<<nyr+2<<"08BMSY"   <<",";  //want probability B2016<0.8BMSY - this will be< 1 if true
    ofsP<<"B"<<nyr+2<<"04BMSY"   <<",";  //want probability B2016<0.4BMSY - this will be < 1 if true
    ofsP<<"FMSY"                 <<",";
    ofsP<<"F"<<nyr+1<<"FMSY"     <<",";
    ofsP<<"UMSY"                 <<",";
    ofsP<<"U"<<nyr+1<<"UMSY";            //want probability F2015>FMSY - this will be > 1 if true
  }
  ofsP<<'\n';
}

void write_proj_output(ofstream &ofsP,
                       int syr,
                       int nyr,
                       double tac,
                       int pyr,
                       dvector p_sbt,
                       dmatrix p_ft,
                       dvar_matrix ft,
                       double bo,
                       dmatrix fmsy,
                       dvector bmsy,
                       bool include_msy,
                       bool include_bo){
  // Write the projection output to the file
  ofsP<<tac                        <<",";
  ofsP<<p_sbt(pyr)                 <<",";
  ofsP<<p_sbt(pyr+1)               <<",";
  if(include_bo){
    ofsP<<bo                         <<",";
    ofsP<<0.4*bo                     <<",";
    ofsP<<0.2*bo                     <<",";
  }
  ofsP<<p_sbt(syr)                 <<",";
  ofsP<<p_sbt(pyr+1)/p_sbt(pyr)    <<",";
  ofsP<<p_sbt(pyr+1)/(0.4*bo)      <<",";
  ofsP<<p_sbt(pyr+1)/(0.2*bo)      <<",";
  ofsP<<p_sbt(pyr+1)/p_sbt(syr)    <<",";
  //ofsP<<ft(1)(1,nyr)               <<",";
  ofsP<<ft(1,nyr)                  <<",";
  ofsP<<p_ft(pyr,1)                <<",";
  ofsP<<p_ft(pyr,1)/ft(1,nyr)      <<",";
  ofsP<<(1. - mfexp(-p_ft(pyr,1))) <<",";
  ofsP<<(1. - mfexp(-p_ft(pyr,1)))/(1. - mfexp(-ft(1,nyr)));
  if(include_msy){
    //MSY based ref points
    ofsP<<","<<bmsy                <<",";
    ofsP<<p_sbt(pyr+1)/bmsy          <<",";
    ofsP<<p_sbt(pyr+1)/(0.8*bmsy)    <<",";
    ofsP<<p_sbt(pyr+1)/(0.4*bmsy)    <<",";
    ofsP<<fmsy                       <<",";
    ofsP<<p_ft(pyr,1)/fmsy           <<",";
    ofsP<<(1. - mfexp(-fmsy))        <<",";
    ofsP<<(1. - mfexp(-p_ft(pyr,1)))/(1. - mfexp(-fmsy));
  }
  ofsP<<'\n';
}
