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
              dvar_vector kappa,
              dvar4_array log_sel,
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
  double phie;
  double bhalpha; //Bev-Holt alpha parameter (called so in the rest of iscam -- it's a salmon thing)
  double bhbeta; //Bev-Holt beta parameter (called beta in the rest of iscam)
  double Mbar;
  double Kappa;
  double Ro;
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

  M_bar = colsum(value(M(1).sub(pf_cntrl(3), pf_cntrl(4)))); //mean across years
  M_bar /= pf_cntrl(4) - pf_cntrl(3) + 1;
  Mbar = mean(M_bar); //mean across ages

  Ro = value(ro(1));
  Kappa = value(kappa(1));

  //RF checking that average weights are the same in the slow_MSY code as in calcReferencePoints
  // cout<<"slow_MSY: avg_wt, avg_fec, Mbar, sel, ro "<<endl;
  cout<<"mean wt at age for slow ref pts = "<<avg_wt<<endl;
  cout<<"Mbar for slow ref pts = "<<Mbar<<endl;
  cout<< "Selectivity (fleet 1) for slow ref pts = "<<vd<<endl;
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

  //Get phie0 (unfished spawning biomass per recruit) -- Needed for stock-recruit parameters. Need to calculate it here so correct average weight is used.
  phie = lw*avg_fec;

  //Stock-recruit parameters - need to get these from inside this function (not global environment) to ensure that correct weight at age is used
 //Kappa = (4.*h)/(1.-h);
 bhalpha = Kappa/phie;
 bhbeta = (Kappa-1)/(Ro*phie);

  cout<<"bhalpha for slow ref pts = "<<bhalpha<<endl;
  cout<<"bhbeta for slow ref pts = "<<bhbeta<<endl;

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
	     if(t <sage){
		        Nn(t + 1)(sage,nage) = Ro * lx;
		        //cout<<Nn(t + 1)(sage,nage)<<endl;
            }else{
			  Nn(t + 1)(sage) = bhalpha * Ssb(t+1-sage) /  (1 + bhbeta * Ssb(t+1-sage));
			  Nn(t + 1)(sage + 1, nage) = ++elem_prod(Nn(t)(sage,nage-1),Ss(sage, nage - 1));
			  Nn(t + 1, nage) += Nn(t, nage) * Ss(nage);
          } //end if

		  //this is correct ... returns bo, MSY and FMSY verified in spreadsheet calcs. FMSY ref pts do not match SM's. B0 matches SM's.
		  Ssb(t + 1) = elem_prod(Nn(t + 1),avg_fec) * stmp;

		  //catch
		  for(j = sage; j <= nage; j++){
			Cc(t, j) = ((ftest(k) * vd(j)) / za(j)) * (1.-exp(-(za(j)))) *
			  Nn(t, j) * avg_wt(j);
		  }
		  Y(t) = sum(Cc(t));

		} //end t
			ye(k) = Y(Nyr);
			be(k) = Ssb(Nyr);

			//Testing the numbers work when sage > 1 (OK)
			//if(k == 1001) 	cout<<"ftest = "<<ftest(k)<<endl<<"Nn = "<<Nn<<endl;

  } //end k
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


