#include "../../include/utilities.h"
#include "../../include/Logger.h"

//Extra test functions by RF to test ref points
//Called by run_FRP() in calcReferencePoints
void slow_msy(dvector& ftest,
              dvector& ye,
              dvector& be,
              double& msy,
              double& fmsy,
              double& bmsy,
              double& bo,
              int sage,
              int nage,
              int nyr,
              dvar3_array M,
              dmatrix dWt_bar,
              dmatrix ma,
              dvar_vector ro,
              dvar4_array log_sel,
              dvar_vector so,
              dvar_vector beta,
              dvector d_iscamCntrl,
              dvector pf_cntrl){
  int i;
  int j;
  int k;
  int t;
  int NF = size_count(ftest);
  int Nyr = 100; //number of years to run out the model
  ye.initialize();
  be.initialize();
  double sa;
  dvector za(sage, nage);
  za.initialize();
  dvector saf(sage, nage);
  saf.initialize();
  dvector lx(sage, nage);
  lx.initialize();
  dvector vd(sage, nage);
  dvector lw(sage, nage);
  lw.initialize();
  vd.initialize();
  dvector avg_wt(sage, nage);
  avg_wt.initialize();
  dvector avg_fec(sage, nage);
  dvector M_bar(sage, nage);

  avg_wt = dWt_bar(1);
  avg_fec = elem_prod(dWt_bar(1), ma(1));
  vd = exp(value(log_sel(1)(1)(nyr)));

  double Ro = value(ro(1));
  //double CR = value(kappa(1));
  double Mbar;
  M_bar = colsum(value(M(1).sub(pf_cntrl(3), pf_cntrl(4)))); //mean across years
  M_bar /= pf_cntrl(4) - pf_cntrl(3) + 1;
  Mbar = mean(M_bar); //mean across ages
  //RF checking that average weights are the same in the slow_MSY code as in calcReferencePoints
  // cout<<"slow_MSY: avg_wt, avg_fec, Mbar, sel, ro "<<endl;
  // cout<<avg_wt<<endl;
  // cout<< avg_fec<<endl;
  // cout<< Mbar<<endl;
  // cout<< vd<<endl;
  // cout<< Ro<<endl;

  dmatrix Nn(1, Nyr + 1, sage, nage); //Numbers at age
  dmatrix Ff(1, Nyr + 1, sage, nage); //Age-specific fishing mortality
  dmatrix Zz(1, Nyr + 1, sage, nage);
  dvector Ss(sage, nage);
  dvector stmp(sage, nage);
  dmatrix Cc(1, Nyr, sage, nage);
  dvector Ssb(1, Nyr + 1);
  dvector Bb(1, Nyr + 1);
  dvector Y(1, Nyr); //predicted catch biomass

  //dvector finaly(1,NF);
  //dvector finalb(1,NF);

  //unfished
  sa = mfexp(-Mbar);
  lx(sage) = 1.0;
  lw(sage) = 1.0;
  for(i = sage; i <= nage; i++){
    if(i >sage){
    	lx(i)=lx(i-1)*sa;
    }
    lw(i) = lx(i) * mfexp(-Mbar*d_iscamCntrl(13)); //correction for spawn timing
  }
  lx(nage) /= (1 - sa);
  lw(nage) /= (1 - sa);

  //Initialize model - same for all F scenarios
  for(j = sage; j <= nage; j++){
    Nn(1, j) = Ro * lx(j);
  }

  for(k = 1; k <= NF; k++){
    za = (Mbar + ftest(k) * vd);
    Ss = mfexp(-za);

    //  Compute spawning biomass at time of spawning.
    stmp      = mfexp(-za*d_iscamCntrl(13)); //RF this returns a vector(sage,nage) of 1s if cntrl(13) is set to zero; and a vector ~ 0.69 if cntrl(13) set to 1
    Ssb(1) = elem_prod(Nn(1), avg_fec)*stmp; //
   /*
    LOG<<"Nn "<<'\n'<<Nn(1)<<'\n';
    LOG<<"SS "<<'\n'<<Ss<<'\n';
    LOG<<"SSb "<<'\n'<<Ssb(1)<<'\n';
    LOG<<"wt "<<avg_wt<<'\n';
    */
    for(t = 1; t <= Nyr; t++){
      Nn(t + 1)(sage + 1, nage) = ++elem_prod(Nn(t)(sage,nage-1),Ss(sage, nage - 1));
      Nn(t + 1, nage) += Nn(t, nage) * Ss(nage);

      if(t ==1){
        Nn(t + 1)(sage) = Ro;
      }
      if(t >1){
        Nn(t + 1)(sage) = value(so(1)) * Ssb(t-1 ) /  (1 + value(beta(1)) * Ssb(t -1));
      }

      //this is correct ... returns bo, MSY and FMSY verified in spreadsheet calcs. FMSY ref pts do not match SM's. B0 matches SM's.
      Ssb(t + 1) = elem_prod(Nn(t + 1),avg_fec) * stmp;

      //catch
      for(j = sage; j <= nage; j++){
        Cc(t, j) = ((ftest(k) * vd(j)) / za(j)) * (1.-exp(-(za(j)))) *
          Nn(t, j) * avg_wt(j);
      }
      Y(t) = sum(Cc(t));
    }
    ye(k) = Y(Nyr);
    be(k) = Ssb(Nyr);

  }
  //ye = finaly;
  //be = finalb;
  bo = be(1); //Unfished biomass

  //get MSY and Fmsy
  msy = max(ye);
  double mtest;
  for(k = 1; k <= NF; k++){
    mtest = ye(k);
    if(mtest == msy){
      fmsy = ftest(k);
    }
    if(mtest == msy){
      bmsy = be(k);
    }
  }

  LOG<<"Ref points from running out model\n";
  LOG<<msy<<'\n';
  LOG<<fmsy<<'\n';
  LOG<<bmsy<<'\n';
  LOG<<bo<<'\n';
}


