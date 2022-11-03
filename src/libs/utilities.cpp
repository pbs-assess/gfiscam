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
			int pyr, // End year for projections
                        bool include_msy,
                        bool include_sbo){
  int i, ig, yr;
  // Write the decision table headers for projection years
  ofsP<<"TAC"<<",";
  for(yr = nyr + 1; yr <= pyr + 1; yr++){
    ofsP<<"B"<<yr<<",";
  }
  for(yr = nyr + 1; yr <= pyr; yr++){
    ofsP<<"R"<<yr<<",";
  }
  if(include_sbo){
    ofsP<<"B0"<<",";
    ofsP<<"04B0"<<",";
    ofsP<<"03B0"<<",";
    ofsP<<"02B0"<<",";
  }
  for(yr = nyr + 1; yr <= pyr; yr++){
    ofsP<<"B"<<yr + 1<<"_B"<<yr<<",";
  }
  if(include_sbo){
    for(yr = nyr; yr <= pyr; yr++){
      ofsP<<"B"<<yr + 1<<"_04B0"<<",";
      ofsP<<"B"<<yr + 1<<"_02B0"<<",";
    }
  }
  for(ig = 1; ig <= n_ags; ig++){
    for(i = 1; i <= nfleet; i++){
      for(yr = nyr + 1; yr<= pyr + 1; yr++){
	ofsP<<"F"<<yr<<"_flt"<<i<<"_sex"<<ig<<",";
	ofsP<<"U"<<yr<<"_flt"<<i<<"_sex"<<ig;
	if(include_msy || !(ig == n_ags && i == nfleet && yr == pyr)){
	  ofsP<<",";
	}
      }
    }
  }
  //ofsP<<"UT"<<",";
  if(include_msy){
    ofsP<<"BMSY"<<",";
    for(yr = nyr + 1; yr <= pyr; yr++){
      ofsP<<"B"<<yr + 1<<"_BMSY"<<",";
      ofsP<<"B"<<yr + 1<<"_08BMSY"<<",";
      ofsP<<"B"<<yr + 1<<"_04BMSY"<<",";
    }
    for(i = 1; i <= nfleet; i++){
      ofsP<<"FMSY_flt"<<i<<",";
      ofsP<<"UMSY_flt"<<i<<",";
      for(ig = 1; ig <= n_ags; ig++){
        for(yr = nyr + 1; yr <= pyr + 1; yr++){
	  ofsP<<"FMSY"<<yr<<"_flt"<<i<<"_sex"<<ig<<",";
	  ofsP<<"UMSY"<<yr<<"_flt"<<i<<"_sex"<<ig;
	  if(!(ig == n_ags && i == nfleet && yr == pyr + 1)){
	    ofsP<<",";
          }
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

  int i, ig, yr;
  //double ut = tac / (tac + p_sbt(pyr));

  // Write the projection output to the file
  ofsP<<tac<<",";
  for(yr = nyr + 1; yr <= pyr + 1; yr++){
    ofsP<<p_sbt(yr)<<",";
  }
  for(yr = nyr + 1; yr <= pyr; yr++){
    ofsP<<p_rt(yr)<<",";
  }
  if(include_sbo){
    ofsP<<sbo<<",";
    ofsP<<0.4 * sbo<<",";
    ofsP<<0.3 * sbo <<",";
    ofsP<<0.2 * sbo <<",";
  }
  for(yr = nyr + 1; yr <= pyr; yr++){
      ofsP<<p_sbt(yr + 1) / p_sbt(yr)<<",";
  }
  if(include_sbo){
    for(yr = nyr; yr <= pyr; yr++){
      ofsP<<p_sbt(yr + 1) / (0.4 * sbo)<<",";
      ofsP<<p_sbt(yr + 1) / (0.2 * sbo)<<",";
    }
  }
  for(ig = 1; ig <= n_ags; ig++){
    for(i = 1; i <= nfleet; i++){
      for(yr = nyr + 1; yr<= pyr + 1; yr++){
	//LOG<<"fleet "<<i<<", sex "<<ig<<", year "<<yr<<", F value:\n"<<p_ft(ig,yr,i)<<"\n";
	//LOG<<"fleet "<<i<<", sex "<<ig<<", year "<<yr<<", U value:\n"<<1.0 - mfexp(-p_ft(ig,yr,i))<<"\n\n";
        ofsP<<p_ft(ig,yr,i)<<",";
        ofsP<<1.0 - mfexp(-p_ft(ig,yr,i));
	if(include_msy || !(ig == n_ags && i == nfleet && yr == pyr)){
	  ofsP<<",";
	}
      }
    }
  }
  /*
  ofsP<<ut<<",";
  for(ig = 1; ig <= n_ags; ig++){
    ofsP<<tac / ((p_N(ig,pyr)(3,nage) * exp(-value(M(1)(nyr,3)))) * dWt_bar(1)(3,nage));
    if(ig != n_ags || include_msy){
      ofsP<<",";
    }
  }
  */
  if(include_msy){
    //LOG<<"fmsy in projection:\n"<<fmsy<<"\n\n";
    //LOG<<"bmsy in projection:\n"<<bmsy(1)<<"\n\n";
    //MSY based ref points
    ofsP<<bmsy(1)<<",";
    for(yr = nyr + 1; yr <= pyr; yr++){
      ofsP<<p_sbt(yr + 1) / bmsy(1)<<",";
      ofsP<<p_sbt(yr + 1) / (0.8 * bmsy(1))<<",";
      ofsP<<p_sbt(yr + 1) / (0.4 * bmsy(1))<<",";
    }
    // The following loop assumed only one group,
    // Which is the '1' in the fmsy(1,i) indexing
    for(i = 1; i <= nfleet; i++){
      ofsP<<fmsy(1, i)<<",";
      ofsP<<1.0 - mfexp(-fmsy(1, i))<<",";
      for(ig = 1; ig <= n_ags; ig++){
        for(yr = nyr + 1; yr<= pyr + 1; yr++){
          ofsP<<p_ft(ig,yr,i) / fmsy(1,i)<<",";
          ofsP<<(1.0 - mfexp(-p_ft(ig,yr,i))) / (1.0 - mfexp(-fmsy(1,i)));
	  if(!(ig == n_ags && i == nfleet && yr == pyr + 1)){
	      ofsP<<",";
          }
	}
      }
    }
  }
  ofsP<<"\n";
}
